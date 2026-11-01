#pragma once
// Node-line Jacobi: the p-multigrid smoother of the cell-wise solvers on stretched elements, in
// place of point Jacobi. Elements are chained into lines along their thin axis (the reference
// axis with the shortest mean edge; two elements chain through a face that lies across the thin
// axis of both). Every free node belongs to exactly one block, a run of nodes along the thin axis
// of its element line, on the rank that counts the node. A block's matrix is the assembled
// operator restricted to it: the assembled diagonal, and off the diagonal the sum over the
// elements that hold both nodes, the line's own element and the lateral neighbours of a node
// line on a face or an edge of it (their line matrices from their stored geometry). Neighbours
// on other ranks have no geometry here and give the diagonal only, which keeps the block
// positive definite. It is banded (bandwidth P), factored once by banded Cholesky and solved by
// one thread. A scaled own-line form without the neighbours was tried and fails where lines are
// cut: rows scaled by different factors lose the near-zero row sums along a strongly coupled
// line (docs/design/cellwise_bp_fcg.md, "Node-line smoother on the GPU").
// Evidence: docs/design/cellwise_bp_fcg.md ("Smoothers for stretched elements", "Real wing
// elements and graded boxes"). Reference: scripts/bp3_smoother_study.py (NodeLines), whose
// blocks this builds in the same order; on one rank they are the same matrices.

#include "backend/distributed/unstructured/solvers/mars_cellwise_lines_host.hpp"
#include "backend/distributed/unstructured/solvers/mars_cellwise_topology.hpp"
#include "backend/distributed/unstructured/utils/mars_gpu_thrust.hpp"

#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <thrust/scan.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace mars {
namespace cellwise {

// ---- node blocks (device) --------------------------------------------------------------

namespace line_detail {

constexpr int kIndexBits = 34;   // copy index bits of the owner code (rank in the bits above)

struct ZeroValues {
    __device__ double operator()(long long) const { return 0.0; }
};

// enc = (rank << kIndexBits | copy + 1) on the counting copy of every free node, 0 elsewhere:
// summed over a node's copies it names the node's counting copy and its rank.
template <int P>
__global__ void __launch_bounds__(kThreads, kMinBlocks)
owner_code_kernel(double* __restrict__ enc, unsigned long long rank_bits, TopologyView T, ElementList set)
{
    for_each_copy_sum<P>(T, set, ZeroValues{}, [&](long long t, double, bool dirichlet, bool counted) {
        enc[t] = counted && !dirichlet ? (double)(rank_bits | (unsigned long long)(t + 1)) : 0.0;
    });
}

// own[t]: the counting copy of copy t's node when this rank counts it, else -1.
__global__ void owner_decode_kernel(const double* __restrict__ sum, long long* __restrict__ own, long long n,
                                    unsigned long long rank)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= n) return;
    const unsigned long long v = (unsigned long long)sum[t];
    own[t] = v != 0 && (v >> kIndexBits) == rank ? (long long)(v & ((1ULL << kIndexBits) - 1)) - 1 : -1;
}

// Packed upper triangle of an n x n symmetric matrix.
MARS_CW_HD inline int packed(int n, int k, int l) { return k * n - k * (k - 1) / 2 + (l - k); }

// K[packed(k, l)]: entry (k, l) of an element's matrix on its node line along axis c at
// transverse position (i, j) (on the ascending other two axes), from the element's stored
// geometry g[s Q^3 + point] (w det J J^-1 J^-T, entries 00 01 02 11 12 22) and the GLL basis B
// and its derivative DB at the Gauss points (Q x n): with V, W the values and derivatives of the
// 1D factors, the 2D sums over the transverse Gauss points are formed once per point along the
// line, then the 1D sums. line_matrices_kernel stores it for every element along its thin axis,
// M[(e n^2 + i n + j) n(n+1)/2 + packed(k, l)].
template <int P>
__device__ inline void line_matrix(const double* __restrict__ g, int c, int i, int j, const double* __restrict__ B,
                                   const double* __restrict__ DB, double (&K)[(P + 1) * (P + 2) / 2])
{
    constexpr int n = P + 1, Q = P + 2, q3 = Q * Q * Q, NP = n * (n + 1) / 2;
    int ka, kb;
    other_axes(c, ka, kb);
    auto s = [](int r, int q) {
        const int lo = r < q ? r : q, hi = r < q ? q : r;
        return lo == 0 ? hi : (lo == 1 ? 2 + hi : 5);
    };
    const double* gaa = g + s(ka, ka) * q3;
    const double* gbb = g + s(kb, kb) * q3;
    const double* gab = g + s(ka, kb) * q3;
    const double* gcc = g + s(c, c) * q3;
    const double* gac = g + s(ka, c) * q3;
    const double* gbc = g + s(kb, c) * q3;
    int stride[3];
    stride[0] = Q * Q, stride[1] = Q, stride[2] = 1;
#pragma unroll
    for (int k = 0; k < NP; ++k) K[k] = 0.0;
    for (int qc = 0; qc < Q; ++qc) {
        double vv = 0.0, ww = 0.0, vw = 0.0;   // weights of V V', W W' and V W' + W V' along the line
        for (int qa = 0; qa < Q; ++qa) {
            const double Vi = B[qa * n + i], Wi = DB[qa * n + i];
            for (int qb = 0; qb < Q; ++qb) {
                const double Vj = B[qb * n + j], Wj = DB[qb * n + j];
                const int pt = qa * stride[ka] + qb * stride[kb] + qc * stride[c];
                vv += gaa[pt] * Wi * Wi * Vj * Vj + gbb[pt] * Vi * Vi * Wj * Wj + 2.0 * gab[pt] * Wi * Vi * Wj * Vj;
                ww += gcc[pt] * Vi * Vi * Vj * Vj;
                vw += gac[pt] * Wi * Vi * Vj * Vj + gbc[pt] * Vi * Vi * Wj * Vj;
            }
        }
#pragma unroll
        for (int k = 0; k < n; ++k)
#pragma unroll
            for (int l = k; l < n; ++l) {
                const double Vk = B[qc * n + k], Wk = DB[qc * n + k], Vl = B[qc * n + l], Wl = DB[qc * n + l];
                K[packed(n, k, l)] += vv * Vk * Vl + ww * Wk * Wl + vw * (Vk * Wl + Wk * Vl);
            }
    }
}

// Entry (a, b) of an element's matrix for any two of its nodes (coordinates a, b), from its stored
// geometry g by quadrature: the couplings a lateral element adds where a node line turns a corner
// in it (build_kernel), which are not on one of its node lines.
template <int P>
__device__ inline double element_entry(const double* __restrict__ g, const int (&a)[3], const int (&b)[3],
                                       const double* __restrict__ B, const double* __restrict__ DB)
{
    constexpr int n = P + 1, Q = P + 2, q3 = Q * Q * Q;
    double sum = 0.0;
    for (int q0 = 0; q0 < Q; ++q0)
        for (int q1 = 0; q1 < Q; ++q1)
            for (int q2 = 0; q2 < Q; ++q2) {
                const int pt = (q0 * Q + q1) * Q + q2;
                const double va0 = B[q0 * n + a[0]], va1 = B[q1 * n + a[1]], va2 = B[q2 * n + a[2]];
                const double vb0 = B[q0 * n + b[0]], vb1 = B[q1 * n + b[1]], vb2 = B[q2 * n + b[2]];
                const double ga[3] = {DB[q0 * n + a[0]] * va1 * va2, va0 * DB[q1 * n + a[1]] * va2,
                                      va0 * va1 * DB[q2 * n + a[2]]};
                const double gb[3] = {DB[q0 * n + b[0]] * vb1 * vb2, vb0 * DB[q1 * n + b[1]] * vb2,
                                      vb0 * vb1 * DB[q2 * n + b[2]]};
                sum += g[pt] * ga[0] * gb[0] + g[3 * q3 + pt] * ga[1] * gb[1] + g[5 * q3 + pt] * ga[2] * gb[2] +
                       g[q3 + pt] * (ga[0] * gb[1] + ga[1] * gb[0]) + g[2 * q3 + pt] * (ga[0] * gb[2] + ga[2] * gb[0]) +
                       g[4 * q3 + pt] * (ga[1] * gb[2] + ga[2] * gb[1]);
            }
    return sum;
}

template <int P>
__global__ void line_matrices_kernel(const double* __restrict__ G, const int* __restrict__ axis,
                                     const double* __restrict__ B, const double* __restrict__ DB, long long E,
                                     double* __restrict__ M)
{
    constexpr int n = P + 1, n2 = n * n, q3 = (P + 2) * (P + 2) * (P + 2), NP = n * (n + 1) / 2;
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= E * n2) return;
    const long long e = t / n2;
    double K[NP];
    line_matrix<P>(G + e * 6 * q3, axis[e], (int)(t % n2) / n, (int)(t % n2) % n, B, DB, K);
    double* out = M + t * NP;
#pragma unroll
    for (int k = 0; k < NP; ++k) out[k] = K[k];
}

// Visits the nodes of chain `chain` (line chain / (P+1)^2, transverse position chain % (P+1)^2
// in the line's first element) in order: visit(m, s, k, i, j, copy), s the position 0 .. P in
// entry m (s = 0 is skipped after the first entry: it is the previous element's last node),
// k its coordinate along the thin axis, (i, j) the transverse position in the element's frame.
template <int P, typename Visit>
__device__ inline void walk_chain(long long chain, const int* __restrict__ start, const int* __restrict__ elem,
                                  const unsigned char* __restrict__ code, Visit&& visit)
{
    constexpr int n = P + 1, n2 = n * n, n3 = n2 * n;
    const long long l = chain / n2;
    const int i0 = (int)(chain % n2) / n, j0 = (int)(chain % n2) % n;
    for (int m = start[l]; m < start[l + 1]; ++m) {
        const int cd = code[m], a = cd & 3;
        const bool rev = (cd >> 2) & 1;
        int i = i0, j = j0, k1, k2;
        transverse(cd >> 3, P, i, j);
        other_axes(a, k1, k2);
        for (int s = m == start[l] ? 0 : 1; s <= P; ++s) {
            int x[3];
            x[a] = rev ? P - s : s;
            x[k1] = i;
            x[k2] = j;
            visit(m, s, x[a], i, j, (long long)elem[m] * n3 + (x[0] * n + x[1]) * n + x[2]);
        }
    }
}

// Every chain proposes itself for the nodes it passes; the lowest chain (longest line first,
// then transverse position) wins a node, as the reference's sequential claims.
template <int P>
__global__ void claim_kernel(const int* __restrict__ start, const int* __restrict__ elem,
                             const unsigned char* __restrict__ code, const long long* __restrict__ own,
                             unsigned long long* __restrict__ claim, long long chains)
{
    const long long chain = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (chain >= chains) return;
    walk_chain<P>(chain, start, elem, code, [&](int, int, int, int, int, long long t) {
        const long long cc = own[t];
        if (cc >= 0) atomicMin(claim + cc, (unsigned long long)chain);
    });
}

// Rows (won nodes) and blocks (runs of won nodes) of every chain.
template <int P>
__global__ void count_kernel(const int* __restrict__ start, const int* __restrict__ elem,
                             const unsigned char* __restrict__ code, const long long* __restrict__ own,
                             const unsigned long long* __restrict__ claim, long long chains, long long* __restrict__ rows,
                             long long* __restrict__ runs)
{
    const long long chain = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (chain >= chains) return;
    long long r = 0, b = 0;
    bool prev = false;
    walk_chain<P>(chain, start, elem, code, [&](int, int, int, int, int, long long t) {
        const long long cc = own[t];
        const bool won = cc >= 0 && claim[cc] == (unsigned long long)chain;
        r += won;
        b += won && !prev;
        prev = won;
    });
    rows[chain] = r;
    runs[chain] = b;
}

// Device views of ElementLines' neighbour tables and of the level's geometry and basis.
struct LineNeighbours {
    const int *face_elem, *edge_start, *edge_elem;
    const unsigned char *face_code, *edge_code;
    const double *G, *B, *DB;
};

// The chain's blocks: their rows (counting copies) and band matrices A[row * (P+1) + d] =
// A(row, row - d): the assembled diagonal, and off it the sums over the line's own element (M)
// and its lateral neighbours on this rank (line_matrix from their geometry), in a fixed order;
// then banded Cholesky in place (L(row, row - d)). At an edge of valence 3 one lateral element
// borders two consecutive line elements and the node line turns a corner in it: its couplings
// between the nodes before and after the joint are added too (element_entry). One thread per
// chain. bad is set when a pivot is not positive.
template <int P>
__global__ void build_kernel(const int* __restrict__ start, const int* __restrict__ elem,
                             const unsigned char* __restrict__ code, const long long* __restrict__ own,
                             const unsigned long long* __restrict__ claim, const double* __restrict__ M,
                             LineNeighbours nb, const double* __restrict__ diag, const long long* __restrict__ row_off,
                             const long long* __restrict__ run_off, long long chains, long long* __restrict__ row_copy,
                             long long* __restrict__ block_row, int* __restrict__ block_len, double* __restrict__ A,
                             int* __restrict__ bad)
{
    constexpr int n = P + 1, n2 = n * n, NP = n * (n + 1) / 2, q3 = (P + 2) * (P + 2) * (P + 2);
    const long long chain = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (chain >= chains) return;
    long long row = row_off[chain], b = run_off[chain] - 1, seg[P + 1], prev_seg[P + 1];
    int kseg[P + 1], prev_kseg[P + 1];
    bool prev = false;
    int cur = -1, ci = 0, cj = 0, prev_cur = -1, prev_ci = 0, prev_cj = 0;
    double K[NP];
    // visit(element, axis, c1, c2, flip) for the lateral elements of entry m's node line at
    // transverse position (i, j): the element across the face it lies on, or the other elements
    // of its edge; the line runs along their axis at coordinates (c1, c2) on their two other
    // axes, reversed when flip.
    auto laterals = [&](int m, int i, int j, auto&& visit) {
        const bool on1 = i == 0 || i == P, on2 = j == 0 || j == P;
        if (on1 && on2) {
            const int g = (i == P) + 2 * (j == P);
            for (int q = nb.edge_start[m * 4 + g]; q < nb.edge_start[m * 4 + g + 1]; ++q) {
                const int ec = nb.edge_code[q];
                visit(nb.edge_elem[q], ec & 3, (ec >> 2 & 1) * P, (ec >> 3 & 1) * P, (ec >> 4) & 1);
            }
        } else if (on1 || on2) {
            const int f = on1 ? (i == P) : 2 + (j == P), other = nb.face_elem[m * 4 + f];
            if (other >= 0) {
                const int fc = nb.face_code[m * 4 + f], an = fc & 3, kn = (fc >> 2) & 3, t = on1 ? j : i;
                int x[3], o1, o2;
                x[kn] = (fc >> 4 & 1) * P;
                x[3 - an - kn] = (fc >> 6 & 1) ? P - t : t;
                other_axes(an, o1, o2);
                visit(other, an, x[o1], x[o2], (fc >> 5) & 1);
            }
        }
    };
    // The segment's pairs within one run, off the diagonal, from Kp in coordinates along the
    // line that run the other way when flip.
    auto add = [&](const double* Kp, bool flip) {
        for (int s1 = 0; s1 < P; ++s1) {
            if (seg[s1] < 0) continue;
            for (int s2 = s1 + 1; s2 <= P && seg[s2] >= 0; ++s2) {
                int k = flip ? P - kseg[s1] : kseg[s1], l = flip ? P - kseg[s2] : kseg[s2];
                if (k > l) {
                    const int t = k;
                    k = l;
                    l = t;
                }
                A[seg[s2] * (P + 1) + (s2 - s1)] += Kp[packed(n, k, l)];
            }
        }
    };
    // Element en, lateral to the previous segment (its line along axis ap at (w1, w2), flip fp) and
    // to this one (along an at (v1, v2), flip fl): its couplings between the previous segment's
    // nodes before the joint and this segment's after it, within one run and the band.
    auto corner = [&](int en, int ap, int w1, int w2, int fp, int an, int v1, int v2, int fl) {
        const double* g = nb.G + (long long)en * 6 * q3;
        int pa, pb, ca, cb;
        other_axes(ap, pa, pb);
        other_axes(an, ca, cb);
        for (int s1 = 0; s1 < P; ++s1) {
            bool run = true;
            for (int q = s1; q <= P; ++q) run = run && prev_seg[q] >= 0;
            if (!run) continue;
            for (int s2 = 1; s2 <= P && seg[s2] >= 0; ++s2) {
                const long long d = seg[s2] - prev_seg[s1];
                if (d > P) break;
                int x[3], y[3];
                x[ap] = fp ? P - prev_kseg[s1] : prev_kseg[s1];
                x[pa] = w1;
                x[pb] = w2;
                y[an] = fl ? P - kseg[s2] : kseg[s2];
                y[ca] = v1;
                y[cb] = v2;
                A[seg[s2] * (P + 1) + d] += element_entry<P>(g, x, y, nb.B, nb.DB);
            }
        }
    };
    // Entry cur at transverse position (ci, cj): the corner couplings with the previous entry,
    // its own matrix, then its lateral elements' line matrices.
    auto flush = [&]() {
        if (prev_cur >= 0)
            laterals(cur, ci, cj, [&](int en, int an, int v1, int v2, int fl) {
                laterals(prev_cur, prev_ci, prev_cj, [&](int ep, int ap, int w1, int w2, int fp) {
                    if (ep == en) corner(en, ap, w1, w2, fp, an, v1, v2, fl);
                });
            });
        add(M + ((long long)elem[cur] * n2 + ci * n + cj) * NP, false);
        laterals(cur, ci, cj, [&](int en, int an, int v1, int v2, int fl) {
            line_matrix<P>(nb.G + (long long)en * 6 * q3, an, v1, v2, nb.B, nb.DB, K);
            add(K, fl);
        });
    };
    walk_chain<P>(chain, start, elem, code, [&](int m, int s, int k, int i, int j, long long t) {
        if (m != cur) {
            if (cur >= 0) {
                flush();
                for (int q = 0; q <= P; ++q) prev_seg[q] = seg[q], prev_kseg[q] = kseg[q];
                prev_cur = cur, prev_ci = ci, prev_cj = cj;
            }
            for (int q = 0; q <= P; ++q) seg[q] = -1;
            if (cur >= 0) {   // position 0: the previous element's last node
                seg[0] = prev ? row - 1 : -1;
                kseg[0] = (code[m] >> 2) & 1 ? P : 0;
            }
            cur = m, ci = i, cj = j;
        }
        const long long cc = own[t];
        const bool won = cc >= 0 && claim[cc] == (unsigned long long)chain;
        kseg[s] = k;
        if (won) {
            if (!prev) {
                ++b;
                block_row[b] = row;
            }
            row_copy[row] = cc;
            seg[s] = row++;
            block_len[b] = (int)(row - block_row[b]);
        }
        prev = won;
    });
    if (cur >= 0) flush();
    for (long long q = run_off[chain]; q <= b; ++q) {
        const long long r0 = block_row[q];
        const int len = block_len[q];
        for (int r = 0; r < len; ++r) {
            double* Lr = A + (r0 + r) * (P + 1);
            Lr[0] = diag[row_copy[r0 + r]];
            const int w = r < P ? r : P;
            for (int d = w; d >= 1; --d) {
                const double* Lj = A + (r0 + r - d) * (P + 1);
                double sum = Lr[d];
                for (int u = d + 1; u <= w; ++u) sum -= Lr[u] * Lj[u - d];
                Lr[d] = sum / Lj[0];
            }
            double sum = Lr[0];
            for (int d = 1; d <= w; ++d) sum -= Lr[d] * Lr[d];
            if (!(sum > 0.0)) atomicExch(bad, 1);
            Lr[0] = sqrt(sum > 0.0 ? sum : 1.0);
        }
    }
}

// Blocks in order of decreasing length, 32 to a warp: lane q % 32 of warp q / 32 holds sorted
// block q; row r of a warp's blocks is stored lane-interleaved at row slot warp_off[w] + r.
template <int P>
__global__ void interleave_kernel(const int* __restrict__ order, const long long* __restrict__ block_row,
                                  const int* __restrict__ block_len, const long long* __restrict__ warp_off,
                                  const double* __restrict__ A, const long long* __restrict__ row_copy, long long blocks,
                                  double* __restrict__ L, long long* __restrict__ copy, int* __restrict__ len)
{
    const long long q = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (q >= blocks) return;
    const int b = order[q], lane = (int)(q & 31), m = block_len[b];
    const long long base = warp_off[q >> 5], r0 = block_row[b];
    len[q] = m;
    for (int r = 0; r < m; ++r) {
        for (int d = 0; d <= P; ++d) L[((base + r) * (P + 1) + d) * 32 + lane] = A[(r0 + r) * (P + 1) + d];
        copy[(base + r) * 32 + lane] = row_copy[r0 + r];
    }
}

// s = DSS(b - y) on the counting copies of free nodes (the only copies the blocks read).
template <int P>
__global__ void __launch_bounds__(kThreads, kMinBlocks)
residual_kernel(const double* __restrict__ b, const double* __restrict__ y, const double* __restrict__ ghost,
                double* __restrict__ s, TopologyView T, ElementList set)
{
    const ResidualValues value{b, y, ghost, T.local * (P + 1) * (P + 1) * (P + 1)};
    for_each_copy_sum<P>(T, set, value, [&](long long t, double sum, bool dirichlet, bool counted) {
        if (counted && !dirichlet) s[t] = sum;
    });
}

// z = L^-T L^-1 s per block, written to the blocks' counting copies (one thread per block).
// Forward: y_r = (s_r - sum_d L(r, r-d) y_(r-d)) / L(r, r); backward from the last row, each
// solved z_r pushes L(r, r-d) z_r to the next P rows' sums. Both keep P values in registers.
template <int P>
__global__ void __launch_bounds__(kThreads)
solve_kernel(const double* __restrict__ L, const long long* __restrict__ copy, const int* __restrict__ len,
             const long long* __restrict__ warp_off, const double* __restrict__ s, double* __restrict__ z,
             long long blocks)
{
    const long long q = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (q >= blocks) return;
    const int lane = (int)(q & 31), m = len[q];
    const long long base = warp_off[q >> 5];
    double ring[P];
#pragma unroll
    for (int d = 0; d < P; ++d) ring[d] = 0.0;
    for (int r = 0; r < m; ++r) {
        const double* Lr = L + (base + r) * (P + 1) * 32 + lane;
        const long long t = copy[(base + r) * 32 + lane];
        double acc = s[t];
#pragma unroll
        for (int d = 1; d <= P; ++d) acc -= Lr[d * 32] * ring[d - 1];
        const double v = acc / Lr[0];
        z[t] = v;
#pragma unroll
        for (int d = P - 1; d > 0; --d) ring[d] = ring[d - 1];
        ring[0] = v;
    }
#pragma unroll
    for (int d = 0; d < P; ++d) ring[d] = 0.0;
    for (int r = m - 1; r >= 0; --r) {
        const double* Lr = L + (base + r) * (P + 1) * 32 + lane;
        const long long t = copy[(base + r) * 32 + lane];
        const double v = (z[t] - ring[0]) / Lr[0];
        z[t] = v;
#pragma unroll
        for (int d = 1; d <= P; ++d) ring[d - 1] = (d < P ? ring[d] : 0.0) + Lr[d * 32] * v;
    }
}

// The Chebyshev step from the block solution on the counting copies: z = DSS (every copy of a
// node gets its value; zero on Dirichlet copies), d = c1 d + c2 z (FIRST: c2 z), x = d when
// assign, else x += d.
template <int P, bool FIRST>
__global__ void __launch_bounds__(kThreads, kMinBlocks)
update_kernel(const double* __restrict__ zc, const double* __restrict__ ghost, double* __restrict__ d,
              double* __restrict__ x, double c1, double c2, bool assign, TopologyView T, ElementList set)
{
    const GhostedValues value{zc, ghost, T.local * (P + 1) * (P + 1) * (P + 1)};
    for_each_copy_sum<P>(T, set, value, [&](long long t, double sum, bool dirichlet, bool) {
        const double zt = dirichlet ? 0.0 : sum;
        const double dn = FIRST ? c2 * zt : c1 * d[t] + c2 * zt;
        d[t] = dn;
        x[t] = assign ? dn : x[t] + dn;
    });
}

inline unsigned grid_for(long long n) { return (unsigned)((n + kThreads - 1) / kThreads); }

}  // namespace line_detail

// The node-line smoother at degree P on one rank's tables. G: the level's stored geometry
// (line_matrices_kernel); B, DB: its GLL basis and derivative at the Gauss points (Q x n, on
// the device); diag: the assembled diagonal on every copy; halo: the level's exchange (null on
// one rank).
template <int P>
class NodeLines {
    static constexpr int n = P + 1, n2 = n * n, n3 = n2 * n, NP = n * (n + 1) / 2;

public:
    NodeLines(const ElementLines& lines, const UnstructuredTopology& T, UnstructuredHalo* halo, const double* G,
              const double* B, const double* DB, const double* diag, MPI_Comm comm, cudaStream_t stream)
        : T_(view(T)), halo_(halo), local_(T.local), stream_(stream)
    {
        using namespace line_detail;
        using thrust::raw_pointer_cast;
        int rank = 0;
        MARS_CELLWISE_MPI(MPI_Comm_rank(comm, &rank));
        const long long copies = local_ * n3;
        if (copies >= (1LL << kIndexBits) - 1 || rank >= (1 << (53 - kIndexBits))) {
            fprintf(stderr, "node lines: %lld copies on rank %d exceed the owner code\n", copies, rank);
            MPI_Abort(comm, 1);
        }
        element_lines_ = lines.lines();
        longest_ = 0;
        for (long long l = 0; l < element_lines_; ++l)
            longest_ = std::max(longest_, lines.start[l + 1] - lines.start[l]);
        std::vector<int> axis(local_, 0);
        for (size_t m = 0; m < lines.elem.size(); ++m) axis[lines.elem[m]] = lines.code[m] & 3;
        thrust::device_vector<int> d_start(lines.start), d_elem(lines.elem), d_axis(axis), d_face(lines.face_elem),
            d_edge_start(lines.edge_start), d_edge(lines.edge_elem);
        thrust::device_vector<unsigned char> d_code(lines.code), d_face_code(lines.face_code),
            d_edge_code(lines.edge_code);
        const int* start = raw_pointer_cast(d_start.data());
        const int* elem = raw_pointer_cast(d_elem.data());
        const unsigned char* code = raw_pointer_cast(d_code.data());

        // Owners: the counting copy of every node this rank counts.
        thrust::device_vector<double> enc(copies), sum(copies);
        if (local_ > 0) {
            static const int grid = resident_grid(owner_code_kernel<P>);
            owner_code_kernel<P><<<grid, kThreads, 0, stream_>>>(raw_pointer_cast(enc.data()),
                                                                 (unsigned long long)rank << kIndexBits, T_,
                                                                 ElementList{nullptr, local_});
            MARS_CELLWISE_CK(cudaGetLastError());
        }
        dss<P>(raw_pointer_cast(enc.data()), raw_pointer_cast(sum.data()), T, halo_, stream_);
        thrust::device_vector<long long> own(copies);
        if (copies > 0)
            owner_decode_kernel<<<grid_for(copies), kThreads, 0, stream_>>>(
                raw_pointer_cast(sum.data()), raw_pointer_cast(own.data()), copies, (unsigned long long)rank);
        MARS_CELLWISE_CK(cudaGetLastError());
        enc.clear(), enc.shrink_to_fit(), sum.clear(), sum.shrink_to_fit();

        // Every element's matrices along its thin axis.
        thrust::device_vector<double> M(local_ * n2 * NP);
        if (local_ > 0)
            line_matrices_kernel<P><<<grid_for(local_ * n2), kThreads, 0, stream_>>>(
                G, raw_pointer_cast(d_axis.data()), B, DB, local_, raw_pointer_cast(M.data()));
        MARS_CELLWISE_CK(cudaGetLastError());

        // Claims, rows and blocks per chain.
        const long long chains = element_lines_ * n2;
        thrust::device_vector<unsigned long long> claim(copies, ~0ULL);
        thrust::device_vector<long long> rows(chains), runs(chains), row_off(chains), run_off(chains);
        if (chains > 0) {
            claim_kernel<P><<<grid_for(chains), kThreads, 0, stream_>>>(start, elem, code, raw_pointer_cast(own.data()),
                                                                        raw_pointer_cast(claim.data()), chains);
            count_kernel<P><<<grid_for(chains), kThreads, 0, stream_>>>(
                start, elem, code, raw_pointer_cast(own.data()), raw_pointer_cast(claim.data()), chains,
                raw_pointer_cast(rows.data()), raw_pointer_cast(runs.data()));
        }
        MARS_CELLWISE_CK(cudaGetLastError());
        const auto policy = gpuPar.on(stream_);
        thrust::exclusive_scan(policy, rows.begin(), rows.end(), row_off.begin());
        thrust::exclusive_scan(policy, runs.begin(), runs.end(), run_off.begin());
        // Element reads below copy on the default stream, which does not wait for stream_.
        MARS_CELLWISE_CK(cudaStreamSynchronize(stream_));
        rows_ = chains > 0 ? (long long)row_off.back() + (long long)rows.back() : 0;
        blocks_ = chains > 0 ? (long long)run_off.back() + (long long)runs.back() : 0;

        thrust::device_vector<long long> row_copy(rows_), block_row(blocks_);
        thrust::device_vector<int> block_len(blocks_), bad(1, 0);
        thrust::device_vector<double> A(rows_ * (P + 1), 0.0);
        const LineNeighbours nb{raw_pointer_cast(d_face.data()),      raw_pointer_cast(d_edge_start.data()),
                                raw_pointer_cast(d_edge.data()),      raw_pointer_cast(d_face_code.data()),
                                raw_pointer_cast(d_edge_code.data()), G,
                                B,                                    DB};
        if (chains > 0)
            build_kernel<P><<<grid_for(chains), kThreads, 0, stream_>>>(
                start, elem, code, raw_pointer_cast(own.data()), raw_pointer_cast(claim.data()),
                raw_pointer_cast(M.data()), nb, diag, raw_pointer_cast(row_off.data()), raw_pointer_cast(run_off.data()),
                chains, raw_pointer_cast(row_copy.data()), raw_pointer_cast(block_row.data()),
                raw_pointer_cast(block_len.data()), raw_pointer_cast(A.data()), raw_pointer_cast(bad.data()));
        MARS_CELLWISE_CK(cudaGetLastError());
        MARS_CELLWISE_CK(cudaStreamSynchronize(stream_));
        if (bad[0]) {
            fprintf(stderr, "node lines: a block matrix is not positive definite (rank %d)\n", rank);
            MPI_Abort(comm, 1);
        }

        // Interleave the blocks by warps of similar length.
        thrust::device_vector<int> order(blocks_), keys(block_len);
        thrust::sequence(policy, order.begin(), order.end());
        thrust::stable_sort_by_key(policy, keys.begin(), keys.end(), order.begin(), thrust::greater<int>());
        MARS_CELLWISE_CK(cudaStreamSynchronize(stream_));
        const long long warps = (blocks_ + 31) / 32;
        std::vector<long long> h_off(warps + 1, 0);
        {
            std::vector<int> h_keys(keys.size());
            thrust::copy(keys.begin(), keys.end(), h_keys.begin());
            for (long long w = 0; w < warps; ++w) h_off[w + 1] = h_off[w] + h_keys[w * 32];
            mean_ = blocks_ > 0 ? (double)rows_ / (double)blocks_ : 0.0;
        }
        slots_ = h_off[warps];
        warp_off_ = h_off;
        L_.assign(slots_ * 32 * (P + 1), 0.0);
        copy_.assign(slots_ * 32, -1);
        len_.resize(blocks_);
        if (blocks_ > 0)
            interleave_kernel<P><<<grid_for(blocks_), kThreads, 0, stream_>>>(
                raw_pointer_cast(order.data()), raw_pointer_cast(block_row.data()), raw_pointer_cast(block_len.data()),
                raw_pointer_cast(warp_off_.data()), raw_pointer_cast(A.data()), raw_pointer_cast(row_copy.data()),
                blocks_, raw_pointer_cast(L_.data()), raw_pointer_cast(copy_.data()), raw_pointer_cast(len_.data()));
        MARS_CELLWISE_CK(cudaGetLastError());
        s_.assign(copies, 0.0);
        z_.assign(copies, 0.0);
        MARS_CELLWISE_CK(cudaStreamSynchronize(stream_));
    }

    // One Chebyshev step with this smoother in place of Jacobi: z = lines(DSS(b - y)) (y null:
    // b), d = c1 d + c2 z (first: c2 z), x = d when assign, else x += d. Its three passes are
    // also callable one by one (the driver times them).
    void step(const double* b, const double* y, double* d, double* x, double c1, double c2, bool first,
              bool assign) const
    {
        residual(b, y);
        solve();
        update(d, x, c1, c2, first, assign);
    }

    // s = DSS(b - y) on the counting copies.
    void residual(const double* b, const double* y) const
    {
        using namespace line_detail;
        double* s = const_cast<double*>(thrust::raw_pointer_cast(s_.data()));
        run_unstructured_pass(halo_, local_, b, y, stream_, [&](const ElementList& set, const double* ghost, int) {
            static const int grid = resident_grid(residual_kernel<P>);
            residual_kernel<P><<<grid, kThreads, 0, stream_>>>(b, y, ghost, s, T_, set);
            MARS_CELLWISE_CK(cudaGetLastError());
            return grid;
        });
    }

    // z = the block solves of s, on the counting copies.
    void solve() const
    {
        using namespace line_detail;
        if (blocks_ > 0)
            solve_kernel<P><<<grid_for(blocks_), kThreads, 0, stream_>>>(
                thrust::raw_pointer_cast(L_.data()), thrust::raw_pointer_cast(copy_.data()),
                thrust::raw_pointer_cast(len_.data()), thrust::raw_pointer_cast(warp_off_.data()),
                thrust::raw_pointer_cast(s_.data()), const_cast<double*>(thrust::raw_pointer_cast(z_.data())), blocks_);
        MARS_CELLWISE_CK(cudaGetLastError());
    }

    // z on every copy (DSS), then the Chebyshev update of d and x.
    void update(double* d, double* x, double c1, double c2, bool first, bool assign) const
    {
        using namespace line_detail;
        const double* z = thrust::raw_pointer_cast(z_.data());
        run_unstructured_pass(halo_, local_, z, nullptr, stream_, [&](const ElementList& set, const double* ghost, int) {
            if (first) {
                static const int grid = resident_grid(update_kernel<P, true>);
                update_kernel<P, true><<<grid, kThreads, 0, stream_>>>(z, ghost, d, x, c1, c2, assign, T_, set);
                MARS_CELLWISE_CK(cudaGetLastError());
                return grid;
            }
            static const int grid = resident_grid(update_kernel<P, false>);
            update_kernel<P, false><<<grid, kThreads, 0, stream_>>>(z, ghost, d, x, c1, c2, assign, T_, set);
            MARS_CELLWISE_CK(cudaGetLastError());
            return grid;
        });
    }

    long long element_lines() const { return element_lines_; }
    int longest() const { return longest_; }
    long long blocks() const { return blocks_; }
    long long rows() const { return rows_; }
    double mean_rows() const { return mean_; }
    long long slots() const { return slots_; }   // padded row slots per lane

private:
    TopologyView T_;
    UnstructuredHalo* halo_;
    long long local_;
    cudaStream_t stream_;
    long long element_lines_ = 0, blocks_ = 0, rows_ = 0, slots_ = 0;
    int longest_ = 0;
    double mean_ = 0.0;
    thrust::device_vector<double> L_, s_, z_;
    thrust::device_vector<long long> copy_, warp_off_;
    thrust::device_vector<int> len_;
};

}  // namespace cellwise
}  // namespace mars
