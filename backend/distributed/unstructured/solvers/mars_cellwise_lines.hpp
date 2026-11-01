#pragma once
// Node-line Jacobi: the p-multigrid smoother of the cell-wise solvers on stretched elements, in
// place of point Jacobi. Elements are chained into lines along their thin axis (the reference
// axis with the shortest mean edge; two elements chain through a face that lies across the thin
// axis of both). Every free node belongs to exactly one block, a run of nodes along the thin axis
// of its element line, on the rank that counts the node. A block's matrix is the sum of its own
// line's element matrices restricted to it, scaled symmetrically to the assembled diagonal, so
// the lateral neighbours enter only through the diagonal and no ghost geometry is needed. It is
// banded (bandwidth P), factored once by banded Cholesky and solved by one thread.
// Evidence and design: docs/design/cellwise_bp_fcg.md ("Smoothers for stretched elements",
// "Real wing elements and graded boxes"). Reference: scripts/bp3_smoother_study.py (NodeLines
// with own = "scaled"), whose blocks this builds in the same order.

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

// ---- element lines (host) --------------------------------------------------------------

// One rank's local elements as lines, longest first: line l is entries start[l] .. start[l+1)
// of elem and code. code bits 0-1: the element's thin axis; bit 2: the line runs from node P
// down to 0 along it; bits 3-5: the map from the line's transverse node position, given in its
// first element's frame, to this element's (bit 3 swaps the two transverse axes, bits 4 and 5
// reverse the first and the second).
struct ElementLines {
    std::vector<int> start{0}, elem;
    std::vector<unsigned char> code;
    long long lines() const { return (long long)start.size() - 1; }
};

namespace line_detail {

MARS_CW_HD inline void transverse(int code, int P, int& i, int& j)
{
    if (code & 1) {
        const int t = i;
        i = j;
        j = t;
    }
    if (code & 2) i = P - i;
    if (code & 4) j = P - j;
}

// The two reference axes other than a, ascending.
MARS_CW_HD inline void other_axes(int a, int& k1, int& k2)
{
    k1 = a == 0 ? 1 : 0;
    k2 = a == 2 ? 1 : 2;
}

// The code c with transverse(c) = transverse(outer) after transverse(inner).
inline int compose(int outer, int inner)
{
    for (int c = 0; c < 8; ++c) {
        bool same = true;
        for (int i = 0; i <= 1 && same; ++i)
            for (int j = 0; j <= 1 && same; ++j) {
                int a = i, b = j, x = i, y = j;
                transverse(inner, 1, a, b);
                transverse(outer, 1, a, b);
                transverse(c, 1, x, y);
                same = a == x && b == y;
            }
        if (same) return c;
    }
    return -1;
}

}  // namespace line_detail

// corners[(e * 8 + c) * 3 + a]: coordinate a of corner c (corner_bits order) of local element
// e; vid[e * 8 + c]: an id of that corner's vertex, equal ids for equal vertices. As
// element_lines of the reference: elements in order, a line started at every element with no
// chain neighbour on one side, closed loops opened at their lowest element, then sorted by
// length (stable).
inline ElementLines build_element_lines(const std::vector<double>& corners, const std::vector<long long>& vid,
                                        long long E)
{
    using line_detail::other_axes;
    std::vector<int> axis(E);
    for (long long e = 0; e < E; ++e) {
        double len[3];
        for (int k = 0; k < 3; ++k) {
            double sum = 0.0;
            for (int c = 0; c < 8; ++c) {
                int b[3];
                corner_bits(c, b[0], b[1], b[2]);
                if (b[k]) continue;
                b[k] = 1;
                const int d = corner_from_bits(b[0], b[1], b[2]);
                double s2 = 0.0;
                for (int a = 0; a < 3; ++a) {
                    const double t = corners[(e * 8 + d) * 3 + a] - corners[(e * 8 + c) * 3 + a];
                    s2 += t * t;
                }
                sum += std::sqrt(s2);
            }
            len[k] = sum / 4.0;
        }
        int a = 0;
        if (len[1] < len[a]) a = 1;
        if (len[2] < len[a]) a = 2;
        axis[e] = a;
    }
    // The faces across each element's thin axis, sorted by their vertices: a face listed twice
    // chains its two elements.
    struct Face {
        std::array<long long, 4> v;
        long long e;
        int side;
    };
    std::vector<Face> faces;
    faces.reserve(2 * E);
    for (long long e = 0; e < E; ++e)
        for (int side = 0; side < 2; ++side) {
            Face f{{}, e, side};
            int m = 0;
            for (int c = 0; c < 8; ++c) {
                int b[3];
                corner_bits(c, b[0], b[1], b[2]);
                if (b[axis[e]] == side) f.v[m++] = vid[e * 8 + c];
            }
            std::sort(f.v.begin(), f.v.end());
            faces.push_back(f);
        }
    std::sort(faces.begin(), faces.end(), [](const Face& x, const Face& y) {
        return x.v != y.v ? x.v < y.v : (x.e != y.e ? x.e < y.e : x.side < y.side);
    });
    std::vector<long long> nxt(2 * E, -1);   // nxt[e * 2 + side]: neighbour * 2 + its side
    for (size_t i = 0; i + 1 < faces.size(); ++i)
        if (faces[i].v == faces[i + 1].v) {
            nxt[faces[i].e * 2 + faces[i].side] = faces[i + 1].e * 2 + faces[i + 1].side;
            nxt[faces[i + 1].e * 2 + faces[i + 1].side] = faces[i].e * 2 + faces[i].side;
            ++i;
        }
    // Walk: a line entered through side `entry` (0: the axis' low face, 1: its high face).
    std::vector<char> seen(E, 0);
    std::vector<std::vector<std::pair<long long, int>>> lines;
    auto walk = [&](long long e, int entry) {
        std::vector<std::pair<long long, int>> line;
        while (!seen[e]) {
            seen[e] = 1;
            line.emplace_back(e, entry);
            const long long n = nxt[e * 2 + (1 - entry)];
            if (n < 0) break;
            e = n / 2;
            entry = (int)(n % 2);
        }
        lines.push_back(std::move(line));
    };
    for (long long e = 0; e < E; ++e)
        for (int side = 0; side < 2; ++side)
            if (!seen[e] && nxt[e * 2 + side] < 0) walk(e, side);
    for (long long e = 0; e < E; ++e)
        if (!seen[e]) walk(e, 0);
    std::stable_sort(lines.begin(), lines.end(), [](const auto& x, const auto& y) { return x.size() > y.size(); });

    ElementLines out;
    out.elem.reserve(E);
    out.code.reserve(E);
    for (const auto& line : lines) {
        int cum = 0;
        for (size_t m = 0; m < line.size(); ++m) {
            const long long e = line[m].first;
            const int entry = line[m].second, a = axis[e];
            if (m > 0) {
                // The previous element's exit face and this one's entry face are the same 4
                // vertices; their transverse corner positions give the map between the frames.
                const long long p = line[m - 1].first;
                const int pa = axis[p], pexit = 1 - line[m - 1].second;
                int pk1, pk2, k1, k2;
                other_axes(pa, pk1, pk2);
                other_axes(a, k1, k2);
                auto corner_vid = [&](long long el, int ax, int bit, int t1, int t2, int i, int j) {
                    int b[3];
                    b[ax] = bit;
                    b[t1] = i;
                    b[t2] = j;
                    return vid[el * 8 + corner_from_bits(b[0], b[1], b[2])];
                };
                int rel = -1;
                for (int c = 0; c < 8 && rel < 0; ++c) {
                    bool ok = true;
                    for (int i = 0; i <= 1 && ok; ++i)
                        for (int j = 0; j <= 1 && ok; ++j) {
                            int x = i, y = j;
                            line_detail::transverse(c, 1, x, y);
                            ok = corner_vid(p, pa, pexit, pk1, pk2, i, j) == corner_vid(e, a, entry, k1, k2, x, y);
                        }
                    if (ok) rel = c;
                }
                cum = line_detail::compose(rel, cum);
            }
            out.elem.push_back((int)e);
            out.code.push_back((unsigned char)(a | (entry << 2) | (cum << 3)));
        }
        out.start.push_back((int)out.elem.size());
    }
    return out;
}

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

// M[(e * n^2 + i n + j) * n(n+1)/2 + packed(k, l)]: entry (k, l) of element e's matrix on its
// node line along its thin axis at transverse position (i, j) (the ascending other two axes),
// from the stored geometry G[(e * 6 + s) Q^3 + point] (w det J J^-1 J^-T, entries 00 01 02 11
// 12 22) and the GLL basis B and its derivative DB at the Gauss points (Q x n): with V, W the
// values and derivatives of the 1D factors, the 2D sums over the transverse Gauss points are
// formed once per point along the line, then the 1D sums.
template <int P>
__global__ void line_matrices_kernel(const double* __restrict__ G, const int* __restrict__ axis,
                                     const double* __restrict__ B, const double* __restrict__ DB, long long E,
                                     double* __restrict__ M)
{
    constexpr int n = P + 1, n2 = n * n, Q = P + 2, q3 = Q * Q * Q, NP = n * (n + 1) / 2;
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= E * n2) return;
    const long long e = t / n2;
    const int i = (int)(t % n2) / n, j = (int)(t % n2) % n, c = axis[e];
    int ka, kb;
    other_axes(c, ka, kb);
    auto s = [](int r, int q) {
        const int lo = r < q ? r : q, hi = r < q ? q : r;
        return lo == 0 ? hi : (lo == 1 ? 2 + hi : 5);
    };
    const double* g = G + e * 6 * q3;
    const double* gaa = g + s(ka, ka) * q3;
    const double* gbb = g + s(kb, kb) * q3;
    const double* gab = g + s(ka, kb) * q3;
    const double* gcc = g + s(c, c) * q3;
    const double* gac = g + s(ka, c) * q3;
    const double* gbc = g + s(kb, c) * q3;
    int stride[3];
    stride[0] = Q * Q, stride[1] = Q, stride[2] = 1;
    double K[NP] = {};
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

// The chain's blocks: their rows (counting copies), band matrices A[row * (P+1) + d] =
// A(row, row - d) summed from the line's element matrices, scaled to the assembled diagonal
// and factored in place by banded Cholesky (L(row, row - d)). One thread per chain, so the sums
// run in a fixed order. bad is set when a pivot is not positive.
template <int P>
__global__ void build_kernel(const int* __restrict__ start, const int* __restrict__ elem,
                             const unsigned char* __restrict__ code, const long long* __restrict__ own,
                             const unsigned long long* __restrict__ claim, const double* __restrict__ M,
                             const double* __restrict__ diag, const long long* __restrict__ row_off,
                             const long long* __restrict__ run_off, long long chains, long long* __restrict__ row_copy,
                             long long* __restrict__ block_row, int* __restrict__ block_len, double* __restrict__ A,
                             double* __restrict__ scale, int* __restrict__ bad)
{
    constexpr int n = P + 1, n2 = n * n, NP = n * (n + 1) / 2;
    const long long chain = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (chain >= chains) return;
    long long row = row_off[chain], b = run_off[chain] - 1, seg[P + 1];
    int kseg[P + 1];
    bool prev = false;
    int cur = -1;
    const double* Me = nullptr;
    auto flush = [&]() {   // the element's entries for the pairs of the segment within one run
        for (int s1 = 0; s1 <= P; ++s1) {
            if (seg[s1] < 0) continue;
            for (int s2 = s1; s2 <= P && seg[s2] >= 0; ++s2) {
                const int k = kseg[s1] < kseg[s2] ? kseg[s1] : kseg[s2], l = kseg[s1] < kseg[s2] ? kseg[s2] : kseg[s1];
                A[seg[s2] * (P + 1) + (s2 - s1)] += Me[packed(n, k, l)];
            }
        }
    };
    walk_chain<P>(chain, start, elem, code, [&](int m, int s, int k, int i, int j, long long t) {
        if (m != cur) {
            if (cur >= 0) flush();
            Me = M + ((long long)elem[m] * n2 + i * n + j) * NP;
            for (int q = 0; q <= P; ++q) seg[q] = -1;
            if (cur >= 0) {   // position 0: the previous element's last node
                seg[0] = prev ? row - 1 : -1;
                kseg[0] = (code[m] >> 2) & 1 ? P : 0;
            }
            cur = m;
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
        for (int r = 0; r < len; ++r) scale[r0 + r] = sqrt(diag[row_copy[r0 + r]] / A[(r0 + r) * (P + 1)]);
        for (int r = 0; r < len; ++r)
            for (int d = 0; d <= P && d <= r; ++d) A[(r0 + r) * (P + 1) + d] *= scale[r0 + r] * scale[r0 + r - d];
        for (int r = 0; r < len; ++r) {
            double* Lr = A + (r0 + r) * (P + 1);
            const int w = r < P ? r : P;
            for (int d = w; d >= 1; --d) {
                const double* Lj = A + (r0 + r - d) * (P + 1);
                double sum = Lr[d];
                for (int t = d + 1; t <= w; ++t) sum -= Lr[t] * Lj[t - d];
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
        thrust::device_vector<int> d_start(lines.start), d_elem(lines.elem), d_axis(axis);
        thrust::device_vector<unsigned char> d_code(lines.code);
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
        thrust::device_vector<double> A(rows_ * (P + 1), 0.0), scale(rows_);
        if (chains > 0)
            build_kernel<P><<<grid_for(chains), kThreads, 0, stream_>>>(
                start, elem, code, raw_pointer_cast(own.data()), raw_pointer_cast(claim.data()),
                raw_pointer_cast(M.data()), diag, raw_pointer_cast(row_off.data()), raw_pointer_cast(run_off.data()),
                chains, raw_pointer_cast(row_copy.data()), raw_pointer_cast(block_row.data()),
                raw_pointer_cast(block_len.data()), raw_pointer_cast(A.data()), raw_pointer_cast(scale.data()),
                raw_pointer_cast(bad.data()));
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
    // b), d = c1 d + c2 z (first: c2 z), x = d when assign, else x += d.
    void step(const double* b, const double* y, double* d, double* x, double c1, double c2, bool first,
              bool assign) const
    {
        using namespace line_detail;
        double* s = const_cast<double*>(thrust::raw_pointer_cast(s_.data()));
        double* z = const_cast<double*>(thrust::raw_pointer_cast(z_.data()));
        run_unstructured_pass(halo_, local_, b, y, stream_, [&](const ElementList& set, const double* ghost, int) {
            static const int grid = resident_grid(residual_kernel<P>);
            residual_kernel<P><<<grid, kThreads, 0, stream_>>>(b, y, ghost, s, T_, set);
            MARS_CELLWISE_CK(cudaGetLastError());
            return grid;
        });
        if (blocks_ > 0)
            solve_kernel<P><<<grid_for(blocks_), kThreads, 0, stream_>>>(
                thrust::raw_pointer_cast(L_.data()), thrust::raw_pointer_cast(copy_.data()),
                thrust::raw_pointer_cast(len_.data()), thrust::raw_pointer_cast(warp_off_.data()), s, z, blocks_);
        MARS_CELLWISE_CK(cudaGetLastError());
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
