#pragma once
// Element-local (cell-wise) DSS on block-structured hex meshes: macro elements
// (blocks) of B^3 elements each, stored block after block, the elements of a block in
// lexicographic order (x slowest). Element e = blk B^3 + (lx B + ly) B + lz, its node
// (a, b, c) at e n^3 + (a n + b) n + c, with a along the block's axis 0.
//
// The scheme of M. Wichrowski, "Coalesced Matrix-Free Finite Elements in Cell-Wise
// Storage", arXiv:2607.02335 (2026), Sec. 4:
//  - Inside a block, one pass per axis sums the two copies across every interface
//    normal to that axis and writes the sum to both. Edge and vertex nodes come out
//    summed as a byproduct: each pass adds up what the previous passes summed. No
//    tables, no atomics: every copy belongs to exactly one interface of a pass.
//  - After the passes, all copies of a node inside one block hold that block's sum. The
//    skeleton (block faces, edges and corners) then adds ONE copy per block and writes
//    the total to every copy. Adding all copies would count a block's sum once per copy.
// Instead of the passes, the first stage can also run one thread per block node, which
// sums the node's copies in the block (up to 8) and writes the sum to each: every copy
// is read once and written once, and no thread depends on another.
// The skeleton tables are the unstructured ones (mars_cellwise_topology.hpp) built over
// the blocks as hexes, so the blocks may form any conforming hex mesh, in any
// orientation. One rank for now (no ghost blocks).
//
// The per-thread bodies are host-device functions, checked on the host against
// assembled sums by tests/cellwise_block_dss_host_check.cpp.

#include "backend/distributed/unstructured/solvers/mars_cellwise_hex.hpp"

#include <algorithm>

namespace mars {
namespace cellwise {

// The block skeleton: face_nbr, face_code (own frame bits 0-2, other's bits 3-5),
// edge and vertex stars of the blocks, as UnstructuredTopology stores them.
struct BlockSkeleton {
    const int *face_nbr, *edge_off, *edge_ent, *vert_off, *vert_ent;
    const unsigned char* face_code;
    long long blocks, edges, vertices;
};

template <int P, int B>
struct BlockShape {
    static constexpr int n = P + 1, n3 = n * n * n, B3 = B * B * B;
    static constexpr int Pm = B * P;   // block nodes per axis: Pm + 1
};

// The copies of block node coordinate X (0 .. Pm) along one axis: element index l and
// local coordinate a; two when X lies on an interface between elements. The first is
// the block's representative.
template <int P, int B>
MARS_CW_HD inline int axis_copies(int X, int (&l)[2], int (&a)[2])
{
    constexpr int Pm = B * P;
    if (X == Pm) {
        l[0] = B - 1;
        a[0] = P;
        return 1;
    }
    l[0] = X / P;
    a[0] = X % P;
    if (a[0] == 0 && X > 0) {
        l[1] = l[0] - 1;
        a[1] = P;
        return 2;
    }
    return 1;
}

// f(copy) for every copy of block node x of block blk, the representative first.
template <int P, int B, typename F>
MARS_CW_HD inline void for_block_copies(long long blk, const int (&x)[3], F&& f)
{
    using S = BlockShape<P, B>;
    int l0[2], a0[2], l1[2], a1[2], l2[2], a2[2];
    const int m0 = axis_copies<P, B>(x[0], l0, a0), m1 = axis_copies<P, B>(x[1], l1, a1),
              m2 = axis_copies<P, B>(x[2], l2, a2);
    for (int i = 0; i < m0; ++i)
        for (int j = 0; j < m1; ++j)
            for (int k = 0; k < m2; ++k)
                f((blk * S::B3 + ((long long)l0[i] * B + l1[j]) * B + l2[k]) * S::n3 +
                  (a0[i] * S::n + a1[j]) * S::n + a2[k]);
}

template <int P, int B>
MARS_CW_HD inline long long block_representative(long long blk, const int (&x)[3])
{
    long long r = -1;
    for_block_copies<P, B>(blk, x, [&](long long c) {
        if (r < 0) r = c;
    });
    return r;
}

// Block node coordinates of face f at (i, j), edge k at t, corner c, at block degree Pm
// (the conventions of face_node, edge_node and corner_node).
template <int Pm>
MARS_CW_HD inline void block_face_xyz(int f, int i, int j, int (&x)[3])
{
    const int axis = f / 2;
    int o0, o1;
    other_axes(axis, o0, o1);
    x[axis] = (f % 2) * Pm;
    x[o0] = i;
    x[o1] = j;
}

template <int Pm>
MARS_CW_HD inline void block_edge_xyz(int k, int t, int (&x)[3])
{
    const int axis = k / 4;
    int o0, o1;
    other_axes(axis, o0, o1);
    x[axis] = t;
    x[o0] = ((k / 2) % 2) * Pm;
    x[o1] = (k % 2) * Pm;
}

template <int Pm>
MARS_CW_HD inline void block_corner_xyz(int c, int (&x)[3])
{
    int b[3];
    corner_bits(c, b[0], b[1], b[2]);
    for (int d = 0; d < 3; ++d) x[d] = b[d] * Pm;
}

// ---- pass along z (out of place) -------------------------------------------------------
// Copy t of the whole vector: in[t] plus the other copy across its z interface, if any.
// It reads every value once and writes every value once, so it also serves as the copy
// from in to out; z faces are strided inside an element, so this pass touches every
// sector of the vector anyway.
template <int P, int B>
MARS_CW_HD inline void block_pass_z(const double* __restrict__ in, double* __restrict__ out, long long t)
{
    using S = BlockShape<P, B>;
    const int c = (int)(t % S::n), lz = (int)((t / S::n3) % B);
    double s = in[t];
    if (c == P && lz < B - 1)
        s += in[t + S::n3 - P];
    else if (c == 0 && lz > 0)
        s = in[t - S::n3 + P] + s;   // same operand order on both copies: same bits
    out[t] = s;
}

// ---- passes along x and y (in place) ---------------------------------------------------
// Interfaces normal to AXIS per block: (B - 1) B^2, each with n^2 node pairs. Thread t:
// the face node fastest (contiguous for x faces, rows of n for y faces), then the
// interfaces, the element along z next (n^3 apart).
template <int P, int B>
constexpr long long block_pairs_per_block = (long long)(B - 1) * B * B * (P + 1) * (P + 1);

template <int P, int B, int AXIS>
MARS_CW_HD inline void block_pass_xy(double* __restrict__ v, long long t)
{
    static_assert(AXIS == 0 || AXIS == 1, "z is the out-of-place pass");
    using S = BlockShape<P, B>;
    constexpr int nn = S::n * S::n;
    const int q = (int)(t % nn);
    long long r = t / nn;
    const int lz = (int)(r % B);
    r /= B;
    const int w = (int)(r % B);   // the element index along the other in-face axis
    r /= B;
    const int la = (int)(r % (B - 1));
    const long long blk = r / (B - 1);
    long long lo, hi;   // the two copies
    if constexpr (AXIS == 0) {
        const long long e = blk * S::B3 + ((long long)la * B + w) * B + lz;
        lo = e * S::n3 + P * nn + q;
        hi = (e + B * B) * S::n3 + q;
    } else {
        const long long e = blk * S::B3 + ((long long)w * B + la) * B + lz;
        const int a = q / S::n, c = q % S::n;
        lo = e * S::n3 + (a * S::n + P) * S::n + c;
        hi = (e + B) * S::n3 + a * nn + c;
    }
    const double s = v[lo] + v[hi];
    v[lo] = s;
    v[hi] = s;
}

// ---- one thread per block node (instead of the passes) --------------------------------
// Thread t: node (X, Y, Z) of block t / (Pm + 1)^3, Z fastest. Every copy belongs to
// exactly one block node, so in and out are each touched once per copy.
template <int P, int B>
constexpr long long block_nodes_per_block = (long long)(B * P + 1) * (B * P + 1) * (B * P + 1);

template <int P, int B>
MARS_CW_HD inline void block_node_sum(const double* __restrict__ in, double* __restrict__ out, long long t)
{
    constexpr int m = B * P + 1;
    const int x[3] = {(int)((t / (m * m)) % m), (int)((t / m) % m), (int)(t % m)};
    const long long blk = t / ((long long)m * m * m);
    double s = 0.0;
    for_block_copies<P, B>(blk, x, [&](long long k) { s += in[k]; });
    for_block_copies<P, B>(blk, x, [&](long long k) { out[k] = s; });
}

// ---- skeleton --------------------------------------------------------------------------
// One thread per skeleton node: interior nodes of block faces (each pair once, from its
// lower copy), then interior nodes of block edges, then block corners. These sets are
// disjoint, so no two threads write the same copy. Stars are summed in their stored
// (canonical) order, faces as lower + upper copy.
template <int P, int B>
MARS_CW_HD inline long long block_skeleton_threads(const BlockSkeleton& S)
{
    constexpr long long m = B * P - 1;
    return S.blocks * 6 * m * m + S.edges * m + S.vertices;
}

template <int P, int B>
MARS_CW_HD inline void block_skeleton_node(double* __restrict__ v, const BlockSkeleton& S, long long t)
{
    constexpr int Pm = B * P, m = Pm - 1;
    const long long faces = S.blocks * 6 * (long long)m * m;
    double s = 0.0;
    if constexpr (m > 0) {   // block faces and edges have interior nodes
    if (t < faces) {
        const int c = (int)(t / (m * m)), q = (int)(t % (m * m));
        const int nbr = S.face_nbr[c];
        if (nbr < c) return;   // the other copy visits this pair, or a boundary face
        const int code = S.face_code[c];
        const int i = 1 + q / m, j = 1 + q % m;
        int I, J, i2, j2, xa[3], xb[3];
        face_to_canonical<Pm>(code & 7, i, j, I, J);
        canonical_to_face<Pm>((code >> 3) & 7, I, J, i2, j2);
        block_face_xyz<Pm>(c % 6, i, j, xa);
        block_face_xyz<Pm>(nbr % 6, i2, j2, xb);
        s = v[block_representative<P, B>(c / 6, xa)] + v[block_representative<P, B>(nbr / 6, xb)];
        for_block_copies<P, B>(c / 6, xa, [&](long long k) { v[k] = s; });
        for_block_copies<P, B>(nbr / 6, xb, [&](long long k) { v[k] = s; });
        return;
    }
    if (t < faces + S.edges * m) {
        t -= faces;
        const long long g = t / m;
        const int T = 1 + (int)(t % m);   // canonical position on the edge
        const int begin = S.edge_off[g], end = S.edge_off[g + 1];
        if (end - begin < 2) return;
        auto at = [&](int ent, int (&x)[3]) {
            const int k = ent & 0x7fffffff;
            block_edge_xyz<Pm>(k % 12, ent < 0 ? Pm - T : T, x);
            return (long long)(k / 12);
        };
        for (int e = begin; e < end; ++e) {
            int x[3];
            const long long blk = at(S.edge_ent[e], x);
            s += v[block_representative<P, B>(blk, x)];
        }
        for (int e = begin; e < end; ++e) {
            int x[3];
            const long long blk = at(S.edge_ent[e], x);
            for_block_copies<P, B>(blk, x, [&](long long k) { v[k] = s; });
        }
        return;
    }
    }
    t -= faces + S.edges * m;
    if (t < S.vertices) {
        const int begin = S.vert_off[t], end = S.vert_off[t + 1];
        if (end - begin < 2) return;
        for (int e = begin; e < end; ++e) {
            int x[3];
            block_corner_xyz<Pm>(S.vert_ent[e] % 8, x);
            s += v[block_representative<P, B>(S.vert_ent[e] / 8, x)];
        }
        for (int e = begin; e < end; ++e) {
            int x[3];
            block_corner_xyz<Pm>(S.vert_ent[e] % 8, x);
            for_block_copies<P, B>(S.vert_ent[e] / 8, x, [&](long long k) { v[k] = s; });
        }
    }
}

#if defined(__CUDACC__) || defined(__HIPCC__)

constexpr int kBlockDssThreads = 256;

template <int P, int B>
__global__ void __launch_bounds__(kBlockDssThreads)
block_pass_z_kernel(const double* __restrict__ in, double* __restrict__ out, long long N)
{
    for (long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x; t < N; t += (long long)gridDim.x * blockDim.x)
        block_pass_z<P, B>(in, out, t);
}

template <int P, int B, int AXIS>
__global__ void __launch_bounds__(kBlockDssThreads) block_pass_xy_kernel(double* __restrict__ v, long long pairs)
{
    for (long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x; t < pairs;
         t += (long long)gridDim.x * blockDim.x)
        block_pass_xy<P, B, AXIS>(v, t);
}

template <int P, int B>
__global__ void __launch_bounds__(kBlockDssThreads)
block_node_sum_kernel(const double* __restrict__ in, double* __restrict__ out, long long nodes)
{
    for (long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x; t < nodes;
         t += (long long)gridDim.x * blockDim.x)
        block_node_sum<P, B>(in, out, t);
}

template <int P, int B>
__global__ void __launch_bounds__(kBlockDssThreads) block_skeleton_kernel(double* __restrict__ v, BlockSkeleton S,
                                                                          long long threads)
{
    for (long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x; t < threads;
         t += (long long)gridDim.x * blockDim.x)
        block_skeleton_node<P, B>(v, S, t);
}

// out = DSS(in) on S.blocks blocks of B^3 elements of degree P: the block stage (three
// passes, or one thread per block node when `nodes`), then the skeleton, on one stream.
// in and out must not alias.
template <int P, int B>
void block_dss(const double* in, double* out, const BlockSkeleton& S, bool nodes, cudaStream_t stream = 0)
{
    using Sh = BlockShape<P, B>;
    auto grid = [](long long work) {
        return (unsigned)std::min<long long>((work + kBlockDssThreads - 1) / kBlockDssThreads, 1LL << 20);
    };
    const long long N = S.blocks * Sh::B3 * Sh::n3;
    if (nodes) {
        const long long m = S.blocks * block_nodes_per_block<P, B>;
        block_node_sum_kernel<P, B><<<grid(m), kBlockDssThreads, 0, stream>>>(in, out, m);
    } else {
        block_pass_z_kernel<P, B><<<grid(N), kBlockDssThreads, 0, stream>>>(in, out, N);
        if constexpr (B > 1) {
            const long long pairs = S.blocks * block_pairs_per_block<P, B>;
            block_pass_xy_kernel<P, B, 0><<<grid(pairs), kBlockDssThreads, 0, stream>>>(out, pairs);
            block_pass_xy_kernel<P, B, 1><<<grid(pairs), kBlockDssThreads, 0, stream>>>(out, pairs);
        }
    }
    const long long sk = block_skeleton_threads<P, B>(S);
    if (sk > 0) block_skeleton_kernel<P, B><<<grid(sk), kBlockDssThreads, 0, stream>>>(out, S, sk);
}

#endif

}  // namespace cellwise
}  // namespace mars
