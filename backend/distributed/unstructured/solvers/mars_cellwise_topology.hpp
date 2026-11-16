#pragma once
// Element-local (cell-wise) DSS on unstructured hex meshes: which copies of a node
// exist, and the kernels that sum them. Design: docs/design/cellwise_unstructured_dss.md;
// executable spec: marsir-mlir/test/cellwise_unstructured_ref.py.
//
// The shared node sets of faces, edges and vertices are disjoint, so each kind is
// handled on its own: a face node has two copies related by the faces' orientations, an
// edge or vertex node one copy per element of its star (any valence). Every copy sums
// all copies of its node in the canonical order (copies sorted by their element's global
// identity, the sorted tuple of its 8 corner keys), so all copies of a node get the same
// bits, nothing is atomic, and the order does not depend on local numbering: the sums
// are the same on any number of ranks.
//
// The tables are built once on the device from each element's corner keys: global keys
// (orientation, canonical order) and dense local ids (grouping by radix sort).
//
// Credit: summing shared nodes by codimension (faces with a 3-bit orientation, edges
// with a reversal bit, vertices) is the face/line/vertex DSS of M. Wichrowski,
// "Coalesced Matrix-Free Finite Elements in Cell-Wise Storage", arXiv:2607.02335 (2026),
// Alg. 3, which applies it at the interfaces of structured macro-blocks. Here it runs at
// element granularity on any conforming hex mesh, with tables built on the GPU from
// corner SFC keys and a canonical summation order that makes the sums independent of
// the rank count.

#include "backend/distributed/unstructured/solvers/mars_cellwise_hex.hpp"
#include "backend/distributed/unstructured/solvers/mars_cellwise_krylov.hpp"
#include "backend/distributed/unstructured/solvers/mars_cellwise_layout.hpp"
#include "backend/distributed/unstructured/utils/mars_gpu_reduce.hpp"

#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/equal.h>
#include <thrust/execution_policy.h>
#include <thrust/for_each.h>
#include <thrust/gather.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/reduce.h>
#include <thrust/scan.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <thrust/transform.h>

#include <algorithm>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <type_traits>
#include <vector>

namespace mars {
namespace cellwise {

// ---- tables ----------------------------------------------------------------------------


// One rank's tables over its elements: its own (local) elements first, then the ghost
// elements of other ranks that share a vertex with them (none on one rank). Copies are
// coded as e * 6 + f (faces), e * 12 + k (edges) or e * 8 + c (corners). Ranges are in
// canonical order throughout.
struct UnstructuredTopology {
    long long elements = 0, local = 0;
    // Faces with two copies on this rank, and with one (the physical boundary on one rank).
    long long face_pairs = 0, single_faces = 0;
    // Edge stars: entries edge_ent[edge_off[g] .. edge_off[g+1]); bit 31 = reversed.
    thrust::device_vector<int> edge_off, edge_ent;
    thrust::device_vector<unsigned char> edge_dirichlet;
    // Vertex stars.
    thrust::device_vector<int> vert_off, vert_ent;
    thrust::device_vector<unsigned char> vert_dirichlet;
    // Per element: the entity bits of the entities whose counting copy (the canonically
    // first) this element holds; kInteriorBit is always set.
    thrust::device_vector<unsigned> counting;
    // Per copy, so each copy of a node can find the others: face_nbr[e * 6 + f] is the
    // other copy of the face (-1 if none), face_code its frame in bits 0-2, the other
    // copy's in bits 3-5 and the Dirichlet flag in bit 6; edge_of[e * 12 + k] is the
    // edge star (bit 31: this copy runs reversed); vert_of[e * 8 + c] the vertex star.
    thrust::device_vector<int> face_nbr, edge_of, vert_of;
    thrust::device_vector<unsigned char> face_code;

    long long edges() const { return (long long)edge_off.size() - 1; }
    long long vertices() const { return (long long)vert_off.size() - 1; }
};

namespace topo_detail {

template <typename T>
T* raw(thrust::device_vector<T>& v) { return thrust::raw_pointer_cast(v.data()); }
template <typename T>
const T* raw(const thrust::device_vector<T>& v) { return thrust::raw_pointer_cast(v.data()); }

// Lexicographic comparison of two elements' sorted corner-key tuples.
__device__ inline int compare_elements(const unsigned long long* tuples, int a, int b)
{
    for (int i = 0; i < 8; ++i) {
        const unsigned long long x = tuples[a * 8 + i], y = tuples[b * 8 + i];
        if (x != y) return x < y ? -1 : 1;
    }
    return 0;
}

__global__ void element_tuples_kernel(const unsigned long long* const* key, long long E,
                                      unsigned long long* tuples)
{
    const long long e = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= E) return;
    unsigned long long t[8];
    for (int c = 0; c < 8; ++c) t[c] = key[c][e];
    for (int i = 1; i < 8; ++i)   // insertion sort of 8 keys
        for (int j = i; j > 0 && t[j] < t[j - 1]; --j) {
            const unsigned long long s = t[j];
            t[j] = t[j - 1];
            t[j - 1] = s;
        }
    for (int c = 0; c < 8; ++c) tuples[e * 8 + c] = t[c];
}

__global__ void vertex_keys_kernel(const int* const* lid, long long E, unsigned long long* k, int* v)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= E * 8) return;
    const long long e = t / 8;
    const int c = (int)(t % 8);
    k[t] = (unsigned long long)(unsigned)lid[c][e];
    v[t] = (int)(e * 8 + c);
}

// Edge key: the two local corner ids, sorted, packed. Value bit 31: reversed (the
// corner at t = 0 has the larger GLOBAL key).
__global__ void edge_keys_kernel(const int* const* lid, const unsigned long long* const* key, long long E,
                                 unsigned long long* k, int* v)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= E * 12) return;
    const long long e = t / 12;
    const int ek = (int)(t % 12);
    const int c0 = edge_corner(ek, 0), c1 = edge_corner(ek, 1);
    const unsigned a = (unsigned)lid[c0][e], b = (unsigned)lid[c1][e];
    k[t] = a < b ? ((unsigned long long)a << 32 | b) : ((unsigned long long)b << 32 | a);
    const bool reversed = key[c0][e] > key[c1][e];
    v[t] = (int)(e * 12 + ek) | (reversed ? (int)0x80000000u : 0);
}

// Face key: the 4 local corner ids, sorted, packed into two words.
__global__ void face_keys_kernel(const int* const* lid, long long E, unsigned long long* hi,
                                 unsigned long long* lo, int* v)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= E * 6) return;
    const long long e = t / 6;
    const int f = (int)(t % 6);
    unsigned s[4];
    for (int q = 0; q < 4; ++q) s[q] = (unsigned)lid[face_corner(f, q)][e];
    for (int i = 1; i < 4; ++i)
        for (int j = i; j > 0 && s[j] < s[j - 1]; --j) {
            const unsigned w = s[j];
            s[j] = s[j - 1];
            s[j - 1] = w;
        }
    hi[t] = (unsigned long long)s[0] << 32 | s[1];
    lo[t] = (unsigned long long)s[2] << 32 | s[3];
    v[t] = (int)(e * 6 + f);
}

// Puts every star's entries in canonical order (insertion sort; stars are small).
// `per` is the number of local entities per element (12 edges, 8 corners).
__global__ void canonical_star_kernel(const int* off, int* ent, long long stars,
                                      const unsigned long long* tuples, int per)
{
    const long long g = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (g >= stars) return;
    for (int i = off[g] + 1; i < off[g + 1]; ++i)
        for (int j = i; j > off[g]; --j) {
            const int a = ent[j - 1], b = ent[j];
            const int ea = (a & 0x7fffffff) / per, eb = (b & 0x7fffffff) / per;
            if (compare_elements(tuples, ea, eb) <= 0) break;
            ent[j - 1] = b;
            ent[j] = a;
        }
}

// Face segments of size 2 become pairs, size 1 singles; pair_lo keeps the canonically
// first copy of each pair (the counting copy).
// A segment of 3 or more copies means a non-conforming or broken mesh.
__global__ void face_pairs_kernel(const int* off, const int* ent, long long faces,
                                  const unsigned long long* tuples, const unsigned long long* const* key,
                                  int* pair_slot, int* single_slot, int* pair_lo, int* single_face,
                                  int* face_nbr, unsigned char* face_code, int* bad)
{
    const long long g = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (g >= faces) return;
    const int n = off[g + 1] - off[g];
    if (n == 1) {
        const int c = ent[off[g]];
        single_face[single_slot[g]] = c;
        face_nbr[c] = -1;
        face_code[c] = 1 << 6;   // on one rank every single face is on the physical boundary
        return;
    }
    if (n != 2) {
        atomicAdd(bad, 1);
        return;
    }
    int a = ent[off[g]], b = ent[off[g] + 1];
    if (compare_elements(tuples, a / 6, b / 6) > 0) {
        const int s = a;
        a = b;
        b = s;
    }
    unsigned long long ka[4], kb[4];
    for (int s = 0; s < 4; ++s) {
        ka[s] = key[face_corner(a % 6, s)][a / 6];
        kb[s] = key[face_corner(b % 6, s)][b / 6];
    }
    const int p = pair_slot[g];
    const int ca = face_frame_code(ka), cb = face_frame_code(kb);
    pair_lo[p] = a;
    face_nbr[a] = b;
    face_nbr[b] = a;
    face_code[a] = (unsigned char)(ca | cb << 3);
    face_code[b] = (unsigned char)(cb | ca << 3);
}

// Every edge and vertex of a Dirichlet face is Dirichlet; the writes all store 1.
__global__ void dirichlet_spread_kernel(const int* single_face, long long singles,
                                        const int* edge_of, const int* vert_of,
                                        unsigned char* edge_dirichlet, unsigned char* vert_dirichlet)
{
    const long long s = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (s >= singles) return;
    const int e = single_face[s] / 6, f = single_face[s] % 6;
    for (int w = 0; w < 4; ++w) {
        edge_dirichlet[edge_of[e * 12 + face_edge(f, w)] & 0x7fffffff] = 1;
        vert_dirichlet[vert_of[e * 8 + face_corner(f, w)]] = 1;
    }
}

// Entity id of every copy, from the sorted segments: of[copy] = segment, keeping the
// copy's reversed bit (bit 31).
__global__ void segment_of_kernel(const int* off, const int* ent, long long segments, int* of)
{
    const long long g = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (g >= segments) return;
    for (int i = off[g]; i < off[g + 1]; ++i) of[ent[i] & 0x7fffffff] = (int)g | (ent[i] & (int)0x80000000u);
}

// The canonically first copy of every entity counts in dot products.
__global__ void counting_kernel(const int* pair_lo, long long pairs, const int* single_face,
                                long long singles, const int* edge_off, const int* edge_ent,
                                long long edges, const int* vert_off, const int* vert_ent,
                                long long vertices, unsigned char* flags)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    long long i = t;
    if (i < pairs) { const int c = pair_lo[i]; flags[(c / 6) * 32 + c % 6] = 1; return; }
    i -= pairs;
    if (i < singles) { const int c = single_face[i]; flags[(c / 6) * 32 + c % 6] = 1; return; }
    i -= singles;
    if (i < edges) {
        const int c = edge_ent[edge_off[i]] & 0x7fffffff;
        flags[(c / 12) * 32 + kEdgeBit + c % 12] = 1;
        return;
    }
    i -= edges;
    if (i < vertices) { const int c = vert_ent[vert_off[i]]; flags[(c / 8) * 32 + kVertexBit + c % 8] = 1; }
}

__global__ void pack_counting_kernel(const unsigned char* flags, long long E, unsigned* bits)
{
    const long long e = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= E) return;
    unsigned b = 1u << kInteriorBit;
    for (int i = 0; i < kInteriorBit; ++i) b |= (unsigned)flags[e * 32 + i] << i;
    bits[e] = b;
}

inline unsigned blocks(long long n) { return (unsigned)((n + kThreads - 1) / kThreads); }

// Sorts (key, value) pairs by key and returns CSR offsets of the equal-key segments.
template <typename Key>
void segments(thrust::device_vector<Key>& keys, thrust::device_vector<int>& vals,
              thrust::device_vector<int>& off)
{
    thrust::sort_by_key(thrust::device, keys.begin(), keys.end(), vals.begin());
    thrust::device_vector<Key> unique(keys.size());
    thrust::device_vector<int> count(keys.size());
    const auto end = thrust::reduce_by_key(thrust::device, keys.begin(), keys.end(),
                                           thrust::constant_iterator<int>(1), unique.begin(), count.begin());
    const long long n = end.second - count.begin();
    off.resize(n + 1);
    off[0] = 0;
    thrust::inclusive_scan(thrust::device, count.begin(), count.begin() + n, off.begin() + 1);
}

}  // namespace topo_detail

// Builds the tables of E elements, of which the first `local` (all when negative) are
// this rank's and the rest ghosts (see UnstructuredHalo). key[c][e]: global key of
// element e's corner c (partition-independent, e.g. its SFC key); lid[c][e]: a local
// id of the same corner (equal ids iff equal keys). Both arrays of 8 device pointers
// live on the device.
inline UnstructuredTopology build_topology(const unsigned long long* const* d_key, const int* const* d_lid,
                                           long long E, long long local = -1)
{
    using namespace topo_detail;
    if (E * 12 > 0x7fffffffLL) {   // copy codes e * 12 + k are 31-bit
        fprintf(stderr, "cell-wise topology: %lld elements on one rank, at most %lld\n", E, 0x7fffffffLL / 12);
        std::abort();
    }
    UnstructuredTopology T;
    T.elements = E;
    T.local = local < 0 ? E : local;
    thrust::device_vector<unsigned long long> tuples(E * 8);
    element_tuples_kernel<<<blocks(E), kThreads>>>(d_key, E, raw(tuples));

    // Vertices.
    {
        thrust::device_vector<unsigned long long> k(E * 8);
        thrust::device_vector<int> v(E * 8);
        vertex_keys_kernel<<<blocks(E * 8), kThreads>>>(d_lid, E, raw(k), raw(v));
        segments(k, v, T.vert_off);
        T.vert_ent = v;
        canonical_star_kernel<<<blocks(T.vertices()), kThreads>>>(raw(T.vert_off), raw(T.vert_ent),
                                                                  T.vertices(), raw(tuples), 8);
    }
    // Edges.
    {
        thrust::device_vector<unsigned long long> k(E * 12);
        thrust::device_vector<int> v(E * 12);
        edge_keys_kernel<<<blocks(E * 12), kThreads>>>(d_lid, d_key, E, raw(k), raw(v));
        segments(k, v, T.edge_off);
        T.edge_ent = v;
        canonical_star_kernel<<<blocks(T.edges()), kThreads>>>(raw(T.edge_off), raw(T.edge_ent),
                                                               T.edges(), raw(tuples), 12);
    }
    // Faces: two-word keys, sorted low word first then (stably) high word.
    thrust::device_vector<int> face_off, face_ent;
    {
        thrust::device_vector<unsigned long long> hi(E * 6), lo(E * 6);
        thrust::device_vector<int> v(E * 6);
        face_keys_kernel<<<blocks(E * 6), kThreads>>>(d_lid, E, raw(hi), raw(lo), raw(v));
        thrust::device_vector<int> perm(E * 6);
        thrust::sequence(thrust::device, perm.begin(), perm.end());
        // Sort a copy: lo and hi are gathered by perm below and must keep the original order.
        thrust::device_vector<unsigned long long> sorted(lo);
        thrust::stable_sort_by_key(thrust::device, sorted.begin(), sorted.end(), perm.begin());
        thrust::gather(thrust::device, perm.begin(), perm.end(), hi.begin(), sorted.begin());
        thrust::stable_sort_by_key(thrust::device, sorted.begin(), sorted.end(), perm.begin());
        // perm is now the face order; rebuild sorted keys and copies in that order.
        thrust::device_vector<unsigned long long> h2(E * 6), l2(E * 6);
        thrust::gather(thrust::device, perm.begin(), perm.end(), hi.begin(), h2.begin());
        thrust::gather(thrust::device, perm.begin(), perm.end(), lo.begin(), l2.begin());
        face_ent.resize(E * 6);
        thrust::gather(thrust::device, perm.begin(), perm.end(), v.begin(), face_ent.begin());
        // Segments of equal (hi, lo): a zip-free form, comparing both words.
        thrust::device_vector<int> head(E * 6);
        const unsigned long long* ph = raw(h2);
        const unsigned long long* pl = raw(l2);
        thrust::transform(thrust::device, thrust::counting_iterator<long long>(0),
                          thrust::counting_iterator<long long>(E * 6), head.begin(),
                          [ph, pl] __host__ __device__(long long i) {
                              return (int)(i == 0 || ph[i] != ph[i - 1] || pl[i] != pl[i - 1]);
                          });
        thrust::device_vector<int> seg(E * 6);
        thrust::inclusive_scan(thrust::device, head.begin(), head.end(), seg.begin());
        const long long faces = seg.back();
        face_off.resize(faces + 1);
        const int* ps = raw(seg);
        int* po = raw(face_off);
        thrust::for_each(thrust::device, thrust::counting_iterator<long long>(0),
                         thrust::counting_iterator<long long>(E * 6), [ps, po] __host__ __device__(long long i) {
                             if (i == 0 || ps[i] != ps[i - 1]) po[ps[i] - 1] = (int)i;
                         });
        face_off[faces] = (int)(E * 6);
    }
    const long long faces = (long long)face_off.size() - 1;
    // Pair and single slots by scans over the segment sizes.
    thrust::device_vector<int> is_pair(faces), is_single(faces), pair_slot(faces), single_slot(faces);
    {
        const int* po = raw(face_off);
        int* pp = raw(is_pair);
        int* ps = raw(is_single);
        thrust::for_each(thrust::device, thrust::counting_iterator<long long>(0),
                         thrust::counting_iterator<long long>(faces), [po, pp, ps] __host__ __device__(long long g) {
                             const int n = po[g + 1] - po[g];
                             pp[g] = n == 2;
                             ps[g] = n == 1;
                         });
    }
    thrust::exclusive_scan(thrust::device, is_pair.begin(), is_pair.end(), pair_slot.begin());
    thrust::exclusive_scan(thrust::device, is_single.begin(), is_single.end(), single_slot.begin());
    const long long npairs = faces ? (long long)pair_slot.back() + is_pair.back() : 0;
    const long long nsingles = faces ? (long long)single_slot.back() + is_single.back() : 0;
    T.face_pairs = npairs;
    T.single_faces = nsingles;
    thrust::device_vector<int> pair_lo(npairs), single_face(nsingles);
    T.face_nbr.resize(E * 6);
    T.face_code.resize(E * 6);
    thrust::device_vector<int> bad(1, 0);
    face_pairs_kernel<<<blocks(faces), kThreads>>>(raw(face_off), raw(face_ent), faces, raw(tuples), d_key,
                                                   raw(pair_slot), raw(single_slot), raw(pair_lo),
                                                   raw(single_face), raw(T.face_nbr), raw(T.face_code),
                                                   raw(bad));
    MARS_CELLWISE_CK(cudaGetLastError());
    if (bad[0] != 0) {
        fprintf(stderr, "cell-wise topology: %d faces shared by 3 or more elements (non-conforming mesh)\n",
                (int)bad[0]);
        std::abort();
    }

    // Dirichlet: on one rank every single face is on the physical boundary.
    T.edge_of.resize(E * 12);
    T.vert_of.resize(E * 8);
    segment_of_kernel<<<blocks(T.edges()), kThreads>>>(raw(T.edge_off), raw(T.edge_ent), T.edges(),
                                                       raw(T.edge_of));
    segment_of_kernel<<<blocks(T.vertices()), kThreads>>>(raw(T.vert_off), raw(T.vert_ent), T.vertices(),
                                                          raw(T.vert_of));
    T.edge_dirichlet.assign(T.edges(), 0);
    T.vert_dirichlet.assign(T.vertices(), 0);
    dirichlet_spread_kernel<<<blocks(nsingles), kThreads>>>(raw(single_face), nsingles, raw(T.edge_of),
                                                            raw(T.vert_of), raw(T.edge_dirichlet),
                                                            raw(T.vert_dirichlet));

    thrust::device_vector<unsigned char> flags(E * 32, 0);
    counting_kernel<<<blocks(npairs + nsingles + T.edges() + T.vertices()), kThreads>>>(
        raw(pair_lo), npairs, raw(single_face), nsingles, raw(T.edge_off), raw(T.edge_ent), T.edges(),
        raw(T.vert_off), raw(T.vert_ent), T.vertices(), raw(flags));
    T.counting.resize(E);
    pack_counting_kernel<<<blocks(E), kThreads>>>(raw(flags), E, raw(T.counting));
    MARS_CELLWISE_CK(cudaGetLastError());
    MARS_CELLWISE_CK(cudaDeviceSynchronize());
    return T;
}

// ---- the DSS -------------------------------------------------------------------------

// Device view of the tables (kernel argument).
struct TopologyView {
    long long local;
    const int *face_nbr, *edge_of, *edge_off, *edge_ent, *vert_of, *vert_off, *vert_ent;
    const unsigned char *face_code, *edge_dirichlet, *vert_dirichlet;
    const unsigned* counting;
};

inline TopologyView view(const UnstructuredTopology& T)
{
    using topo_detail::raw;
    return {T.local,            raw(T.face_nbr),  raw(T.edge_of),         raw(T.edge_off),
            raw(T.edge_ent),    raw(T.vert_of),   raw(T.vert_off),        raw(T.vert_ent),
            raw(T.face_code),   raw(T.edge_dirichlet), raw(T.vert_dirichlet), raw(T.counting)};
}

// Interior nodes of a face and an edge at degree P (36 and 6 at p = 7). The per-degree
// code below takes P as a template argument, p = 7 by default; the tables are the same at
// every degree.
template <int P>
constexpr int face_nodes = (P - 1) * (P - 1);
template <int P>
constexpr int edge_nodes = P - 1;

// ---- gather form: every copy sums its own node ------------------------------------
//
// One warp per element, elements in grid-stride order: the elements in flight form a
// contiguous window, so on an SFC-ordered mesh the neighbours' copies are still in L2
// when they are read. A warp stages its element's 512 values in shared memory with one
// coalesced load, adds the other copy of every face node (one loop), sums the edge and
// vertex stars (a second loop), then writes all 512 results with one coalesced store.
// Interior nodes need no work. Each loop runs one kind of node, so lanes do not wait on
// each other's paths, and the reference-hex index arithmetic comes from tables built
// once per thread block.
// Every copy of a node gets the same bits: a face adds its two copies (addition is
// commutative), edges and vertices sum their stars in canonical order.

template <int P>
struct HexTables {
    static constexpr int F = face_nodes<P> > 0 ? face_nodes<P> : 1;   // no zero-size arrays at p = 1
    unsigned short face_l[6 * F];    // local node of face f at position q = (i - 1) (P - 1) + j - 1
    unsigned char face_map[64 * F];  // frame codes (own | other << 3), own q -> the other's q
    unsigned short edge_l[12 * (P + 1)];   // local node of edge k at position t
    unsigned short corner_l[8];
    unsigned char bit[(P + 1) * (P + 1) * (P + 1)];   // entity bit of every local node
};

template <int P>
__device__ inline void build_hex_tables(HexTables<P>& h)
{
    constexpr int n = P + 1, F = face_nodes<P>;
    if constexpr (P > 1) {
        for (int i = threadIdx.x; i < 64 * F; i += blockDim.x) {
            const int q = i % F, code = i / F;
            int I, J, i2, j2;
            face_to_canonical<P>(code & 7, 1 + q / (P - 1), 1 + q % (P - 1), I, J);
            canonical_to_face<P>(code >> 3, I, J, i2, j2);
            h.face_map[i] = (unsigned char)((i2 - 1) * (P - 1) + j2 - 1);
        }
        for (int i = threadIdx.x; i < 6 * F; i += blockDim.x)
            h.face_l[i] = (unsigned short)face_node<P>(i / F, 1 + (i % F) / (P - 1), 1 + (i % F) % (P - 1));
    }
    for (int i = threadIdx.x; i < 12 * n; i += blockDim.x) h.edge_l[i] = (unsigned short)edge_node<P>(i / n, i % n);
    if (threadIdx.x < 8) h.corner_l[threadIdx.x] = (unsigned short)corner_node<P>(threadIdx.x);
    for (int l = threadIdx.x; l < n * n * n; l += blockDim.x)
        h.bit[l] = (unsigned char)entity_bit(local_entity<P>(l / (n * n), (l / n) % n, l % n));
    __syncthreads();
}

// One pad slot per row of P + 1 spreads the strided nodes of a face over the shared-memory banks.
template <int P>
constexpr int stage_size = (P + 1) * (P + 1) * (P + 2);
template <int P>
__device__ inline int stage_slot(int l) { return l + l / (P + 1); }

// Sum of value over star entries [begin, end) in order; entries are loaded four at a
// time so their value loads overlap.
template <typename Value, typename Index>
__device__ inline double star_sum(const Value& value, const int* __restrict__ ent, int begin, int end,
                                  const Index& index)
{
    double s = 0.0;
    for (int m = begin; m < end; m += 4) {
        long long idx[4];
#pragma unroll
        for (int j = 0; j < 4; ++j) idx[j] = m + j < end ? index(ent[m + j]) : 0;
        double v[4];
#pragma unroll
        for (int j = 0; j < 4; ++j) v[j] = m + j < end ? value(idx[j]) : 0.0;
#pragma unroll
        for (int j = 0; j < 4; ++j)
            if (m + j < end) s += v[j];
    }
    return s;
}

// The local elements a pass visits: ids[0 .. count), or 0 .. count when ids is null.
struct ElementList {
    const int* ids = nullptr;
    long long count = 0;
};

// Value i of the local copies, or of the ghost elements' copies stored after them.
struct GhostedValues {
    const double* local;
    const double* ghost;
    long long local_values;
    __device__ double operator()(long long i) const { return i < local_values ? local[i] : ghost[i - local_values]; }
};

// The residual b - y (b alone when y is null) on the local copies; the ghost copies hold
// the exchanged residual.
struct ResidualValues {
    const double* b;
    const double* y;
    const double* ghost;
    long long local_values;
    __device__ double operator()(long long i) const
    {
        return i < local_values ? b[i] - (y ? y[i] : 0.0) : ghost[i - local_values];
    }
};

// visit(t, s, dirichlet, counted) for every copy t of the elements in `set`, in
// coalesced order: s is the sum over all copies of its node, counted marks the counting
// copy. Called by every thread of the block (it synchronises the block once, to build
// the tables).
template <int P, typename Value, typename Visit>
__device__ inline void for_each_copy_sum(const TopologyView& T, const ElementList& set, const Value& value,
                                         Visit&& visit)
{
    constexpr unsigned kAll = 0xffffffffu;
    constexpr int n = P + 1, n3 = n * n * n, F = face_nodes<P>, Ed = edge_nodes<P>;
    constexpr int kFaceSlots = 6 * F, kEdgeSlots = 12 * Ed, kStarSlots = kEdgeSlots + 8;
    __shared__ HexTables<P> h;
    __shared__ double stage_all[kThreads / 32][stage_size<P>];
    build_hex_tables<P>(h);
    const int lane = threadIdx.x & 31;
    double* stage = stage_all[threadIdx.x >> 5];
    const long long warps = ((long long)gridDim.x * blockDim.x) >> 5;
    for (long long w = ((long long)blockIdx.x * blockDim.x + threadIdx.x) >> 5; w < set.count; w += warps) {
        const long long e = set.ids ? set.ids[w] : w;
        // The descriptor of the entity with bit b sits on lane b. Face lanes: d0 = other
        // copy, d1 = frame codes. Edge lanes: d0 = edge_of entry, [d1, d2) = star.
        // Vertex lanes: [d1, d2) = star.
        int d0 = 0, d1 = 0, d2 = 0, dir = 0;
        if (lane < 6) {
            d0 = T.face_nbr[e * 6 + lane];
            d1 = T.face_code[e * 6 + lane];
            dir = d1 >> 6;
        } else if (lane >= kEdgeBit && lane < kEdgeBit + 12) {
            d0 = T.edge_of[e * 12 + lane - kEdgeBit];
            const int g = d0 & 0x7fffffff;
            d1 = T.edge_off[g];
            d2 = T.edge_off[g + 1];
            dir = T.edge_dirichlet[g];
        } else if (lane >= kVertexBit && lane < kVertexBit + 8) {
            const int g = T.vert_of[e * 8 + lane - kVertexBit];
            d1 = T.vert_off[g];
            d2 = T.vert_off[g + 1];
            dir = T.vert_dirichlet[g];
        }
        const unsigned dirichlet = __ballot_sync(kAll, dir != 0);
        const unsigned counting = T.counting[e];
        const long long base = e * n3;
#pragma unroll
        for (int it = 0; it < (n3 + 31) / 32; ++it) {
            const int l = it * 32 + lane;
            if (l < n3) stage[stage_slot<P>(l)] = value(base + l);
        }
        __syncwarp();
        // Shuffles stay outside conditions: every lane must execute each of them.
        if constexpr (P > 1) {
#pragma unroll
            for (int p0 = 0; p0 < kFaceSlots; p0 += 32) {
                const int p = p0 + lane;
                const int f = p < kFaceSlots ? p / F : 0;
                const int nbr = __shfl_sync(kAll, d0, f);
                const int code = __shfl_sync(kAll, d1, f);
                if (p < kFaceSlots && nbr >= 0) {
                    const int q = h.face_map[(code & 63) * F + p - f * F];
                    stage[stage_slot<P>(h.face_l[p])] += value((long long)(nbr / 6) * n3 + h.face_l[(nbr % 6) * F + q]);
                }
            }
        }
        for (int p0 = 0; p0 < kStarSlots; p0 += 32) {
            const int p = p0 + lane;
            const int k = p < kEdgeSlots ? p / (Ed > 0 ? Ed : 1) : 0;
            const int src = p < kEdgeSlots ? kEdgeBit + k : (p < kStarSlots ? kVertexBit + p - kEdgeSlots : 0);
            const int n0 = __shfl_sync(kAll, d0, src);
            const int n1 = __shfl_sync(kAll, d1, src);
            const int n2 = __shfl_sync(kAll, d2, src);
            if (p < kEdgeSlots) {
                const int t = 1 + p - k * Ed;
                const int tc = n0 < 0 ? P - t : t;   // canonical position on the edge
                stage[stage_slot<P>(h.edge_l[k * n + t])] = star_sum(value, T.edge_ent, n1, n2, [&](int c) {
                    const int m = c & 0x7fffffff;
                    return (long long)(m / 12) * n3 + h.edge_l[(m % 12) * n + (c < 0 ? P - tc : tc)];
                });
            } else if (p < kStarSlots) {
                stage[stage_slot<P>(h.corner_l[p - kEdgeSlots])] = star_sum(
                    value, T.vert_ent, n1, n2, [&](int c) { return (long long)(c / 8) * n3 + h.corner_l[c % 8]; });
            }
        }
        __syncwarp();
#pragma unroll
        for (int it = 0; it < (n3 + 31) / 32; ++it) {
            const int l = it * 32 + lane;
            if (l >= n3) continue;
            const int bit = h.bit[l];
            visit(base + l, stage[stage_slot<P>(l)], ((dirichlet >> bit) & 1u) != 0, ((counting >> bit) & 1u) != 0);
        }
        __syncwarp();   // the next element reuses the stage
    }
}

// Both kernels are compiled for 4 thread blocks per SM (at most 64 registers). Measured
// on GH200, 64^3 rotated elements: DSS 1.21 ms (1.25 uncapped, 2.09 at 5 blocks), DSS +
// Jacobi + 2 dots 2.26 ms (3.22 uncapped, 3.15 at 5 blocks).
constexpr int kMinBlocks = 4;

template <int P = kP>
__global__ void __launch_bounds__(kThreads, kMinBlocks)
unstructured_dss_kernel(const double* __restrict__ in, const double* __restrict__ ghost, double* __restrict__ out,
                        TopologyView T, ElementList set)
{
    const GhostedValues value{in, ghost, T.local * (P + 1) * (P + 1) * (P + 1)};
    for_each_copy_sum<P>(T, set, value, [&](long long t, double s, bool, bool) { out[t] = s; });
}

// Weight of a copy in the assembled inner product: 1 on the counting copy of its node
// (the canonically first), 0 on the others, so every global node counts exactly once.
struct CountingWeight {
    const unsigned* counting;
    __device__ double operator()(const Node& nd) const
    {
        return (counting[nd.e] >> entity_bit(local_entity(nd.a, nd.b, nd.c))) & 1u ? 1.0 : 0.0;
    }
};

// z = P q: the DSS sum, zero on Dirichlet nodes, divided by the assembled diagonal
// elsewhere; with AZ / ZZ also <a, z> and <z, z>, each node counted once.
template <bool AZ, bool ZZ>
__global__ void __launch_bounds__(kThreads, kMinBlocks)
unstructured_precondition_kernel(const double* __restrict__ q, const double* __restrict__ ghost,
                                 double* __restrict__ z, const double* __restrict__ diag,
                                 const double* __restrict__ a, TopologyView T, ElementList set,
                                 double* __restrict__ partial)
{
    constexpr int NV = AZ + ZZ;
    [[maybe_unused]] double acc[NV > 0 ? NV : 1] = {};
    const GhostedValues value{q, ghost, T.local * kN3};
    for_each_copy_sum<kP>(T, set, value, [&](long long t, double s, bool dirichlet, bool counted) {
        const double v = dirichlet ? 0.0 : s / diag[t];
        z[t] = v;
        if constexpr (NV > 0) {
            const double w = counted ? 1.0 : 0.0;
            if constexpr (AZ) acc[0] += a[t] * v * w;
            if constexpr (ZZ) acc[NV - 1] += v * v * w;
        }
    });
    if constexpr (NV > 0) store_partial<NV>(acc, partial);
}

// One Chebyshev step of the p-multigrid smoother at degree P: z = Jacobi(DSS(b - y)),
// zero on Dirichlet copies, with y = A x (null while x = 0); then d = c1 d + c2 z and
// x += d (FIRST: d = c2 z, and x = d when y is null).
template <int P, bool FIRST>
__global__ void __launch_bounds__(kThreads)
unstructured_chebyshev_kernel(const double* __restrict__ b, const double* __restrict__ y,
                              const double* __restrict__ ghost, const double* __restrict__ diag,
                              double* __restrict__ d, double* __restrict__ x, double c1, double c2, TopologyView T,
                              ElementList set)
{
    const ResidualValues value{b, y, ghost, T.local * (P + 1) * (P + 1) * (P + 1)};
    for_each_copy_sum<P>(T, set, value, [&](long long t, double s, bool dirichlet, bool) {
        const double z = dirichlet ? 0.0 : s / diag[t];
        const double dn = FIRST ? c2 * z : c1 * d[t] + c2 * z;
        d[t] = dn;
        x[t] = FIRST && !y ? dn : x[t] + dn;
    });
}

// Block partials of <a, b> over the local copies at degree P, each node counted once.
template <int P>
__global__ void __launch_bounds__(kThreads)
level_wdot_kernel(const double* __restrict__ a, const double* __restrict__ b, const unsigned* __restrict__ counting,
                  long long local, double* __restrict__ partial)
{
    constexpr int n = P + 1, n3 = n * n * n;
    double acc[1] = {0.0};
    for (long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x; t < local * n3;
         t += (long long)gridDim.x * blockDim.x) {
        const int l = (int)(t % n3);
        if ((counting[t / n3] >> entity_bit(local_entity<P>(l / (n * n), (l / n) % n, l % n))) & 1u)
            acc[0] += a[t] * b[t];
    }
    store_partial<1>(acc, partial);
}

// u = 0 on every local copy of a Dirichlet node: a homogeneous Dirichlet condition on an
// element-local vector (e.g. a manufactured solution on quantized coordinates).
template <int P = kP>
__global__ void zero_dirichlet_kernel(double* __restrict__ u, TopologyView T)
{
    constexpr int n = P + 1, n3 = n * n * n;
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= T.local * n3) return;
    const long long e = t / n3;
    const int l = (int)(t % n3);
    const LocalEntity le = local_entity<P>(l / (n * n), (l / n) % n, l % n);
    bool dirichlet = false;
    if (le.kind == 1) dirichlet = (T.face_code[e * 6 + le.index] >> 6) != 0;
    else if (le.kind == 2) dirichlet = T.edge_dirichlet[T.edge_of[e * 12 + le.index] & 0x7fffffff] != 0;
    else if (le.kind == 3) dirichlet = T.vert_dirichlet[T.vert_of[e * 8 + le.index]] != 0;
    if (dirichlet) u[t] = 0.0;
}

template <int P = kP>
void zero_dirichlet(double* d_u, const TopologyView& T, cudaStream_t stream = 0)
{
    const long long n = T.local * (P + 1) * (P + 1) * (P + 1);
    if (n > 0) zero_dirichlet_kernel<P><<<topo_detail::blocks(n), kThreads, 0, stream>>>(d_u, T);
    MARS_CELLWISE_CK(cudaGetLastError());
}

template <int P = kP>
void zero_dirichlet(double* d_u, const UnstructuredTopology& T, cudaStream_t stream = 0)
{
    zero_dirichlet<P>(d_u, view(T), stream);
}

// ---- several ranks ---------------------------------------------------------------------
//
// Each rank's tables cover its own elements and, after them, ghost elements: every
// element of another rank that shares a vertex with one of its own. All copies of a
// shared node are then in the tables, so the canonical order and the counting copy come
// out the same on every rank, and the sums are bit-identical to one rank. The exchange
// fills the ghost elements' shared nodes: each rank sends its copies of the nodes it
// shares with a peer's elements and receives the peer's copies of the same nodes. Both
// sides list these nodes sorted by (peer, element id, local node), so the messages need
// no index lists; the order is checked once at setup.

namespace halo_detail {

using topo_detail::raw;
using topo_detail::blocks;

__host__ __device__ inline unsigned long long splitmix64(unsigned long long z)
{
    z += 0x9e3779b97f4a7c15ULL;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

struct KeyHash {
    __device__ unsigned long long operator()(unsigned long long id, int node) const
    {
        return splitmix64(id) ^ (unsigned long long)node;
    }
};

template <typename T>
MPI_Datatype mpi_type();
template <>
inline MPI_Datatype mpi_type<double>() { return MPI_DOUBLE; }
template <>
inline MPI_Datatype mpi_type<unsigned long long>() { return MPI_UINT64_T; }
template <>
inline MPI_Datatype mpi_type<int>() { return MPI_INT; }

// The interior nodes of a face ((p - 1)^2), an edge (p - 1) or a corner (1) at degree p,
// and the local node of the j-th of them (the node functions of mars_cellwise_hex.hpp
// with a run-time degree: the exchange is set up once per degree).
__device__ inline int entity_nodes(int kind, int p) { return kind == 1 ? (p - 1) * (p - 1) : (kind == 2 ? p - 1 : 1); }
__device__ inline int entity_node(int kind, int index, int j, int p)
{
    int x[3], o0, o1;
    if (kind == 3) {
        int b[3];
        corner_bits(index, b[0], b[1], b[2]);
        for (int k = 0; k < 3; ++k) x[k] = b[k] * p;
    } else if (kind == 1) {
        const int axis = index / 2;
        other_axes(axis, o0, o1);
        x[axis] = (index % 2) * p;
        x[o0] = 1 + j / (p - 1);
        x[o1] = 1 + j % (p - 1);
    } else {
        const int axis = index / 4;
        other_axes(axis, o0, o1);
        x[axis] = 1 + j;
        x[o0] = ((index / 2) % 2) * p;
        x[o1] = (index % 2) * p;
    }
    return (x[0] * (p + 1) + x[1]) * (p + 1) + x[2];
}

// A local element reads ghost copies iff one of its vertex stars holds a ghost element:
// every element sharing a face, edge or vertex with it shares a vertex.
__global__ void boundary_kernel(TopologyView T, int* __restrict__ flag)
{
    const long long e = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= T.local) return;
    int b = 0;
    for (int c = 0; c < 8; ++c) {
        const int g = T.vert_of[e * 8 + c];
        for (int m = T.vert_off[g]; m < T.vert_off[g + 1]; ++m) b |= T.vert_ent[m] / 8 >= T.local;
    }
    flag[e] = b;
}

// Selects the elements whose flag is (or is not) set; a functor because nvcc allows no
// device lambda inside a constructor.
struct FlagIs {
    const int* flag;
    bool set;
    __device__ bool operator()(int e) const { return (flag[e] != 0) == set; }
};

struct Entries {
    int* peer;
    unsigned long long* gid;
    int* node;
    long long* index;
};

// One item per local face, edge star and vertex star. COUNT: how many nodes the item
// sends and receives. FILL: writes them from the given offsets. A local face whose other
// copy is a ghost sends its 36 nodes to the ghost's owner and receives the ghost's. A
// star holding local and ghost copies receives every ghost copy and sends every local
// copy once to each rank that owns a ghost copy.
template <bool FILL>
__global__ void entries_kernel(TopologyView T, long long edges, long long vertices,
                               const unsigned long long* __restrict__ gid, const int* __restrict__ owner, int p,
                               long long* __restrict__ send_at, long long* __restrict__ recv_at, Entries send,
                               Entries recv)
{
    const long long n3 = (long long)(p + 1) * (p + 1) * (p + 1);
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    const long long faces = T.local * 6;
    if (t >= faces + edges + vertices) return;
    long long so = FILL ? send_at[t] : 0, ro = FILL ? recv_at[t] : 0;
    auto emit = [&](Entries& out, long long& o, int peer, long long elem, int kind, int index, long long base) {
        const int n = entity_nodes(kind, p);
        if (FILL)
            for (int j = 0; j < n; ++j) {
                const int l = entity_node(kind, index, j, p);
                out.peer[o + j] = peer;
                out.gid[o + j] = gid[elem];
                out.node[o + j] = l;
                out.index[o + j] = base + l;
            }
        o += n;
    };
    if (t < faces) {
        const long long e = t / 6;
        const int c = T.face_nbr[t];
        if (c >= 0 && c / 6 >= T.local) {
            emit(send, so, owner[c / 6], e, 1, (int)(t % 6), e * n3);
            emit(recv, ro, owner[c / 6], c / 6, 1, c % 6, (long long)(c / 6 - T.local) * n3);
        }
    } else {
        const bool edge = t < faces + edges;
        const long long g = edge ? t - faces : t - faces - edges;
        const int* off = edge ? T.edge_off : T.vert_off;
        const int* ent = edge ? T.edge_ent : T.vert_ent;
        const int per = edge ? 12 : 8, kind = edge ? 2 : 3;
        const int b = off[g], end = off[g + 1];
        bool any_local = false, any_ghost = false;
        for (int m = b; m < end; ++m) ((ent[m] & 0x7fffffff) / per < T.local ? any_local : any_ghost) = true;
        if (any_local && any_ghost)
            for (int m = b; m < end; ++m) {
                const int c = ent[m] & 0x7fffffff;
                const long long el = c / per;
                if (el >= T.local) {
                    emit(recv, ro, owner[el], el, kind, c % per, (el - T.local) * n3);
                    continue;
                }
                for (int j = b; j < end; ++j) {
                    const long long ej = (ent[j] & 0x7fffffff) / per;
                    if (ej < T.local) continue;
                    bool first = true;   // the first ghost copy owned by this peer
                    for (int i = b; i < j && first; ++i) {
                        const long long ei = (ent[i] & 0x7fffffff) / per;
                        first = !(ei >= T.local && owner[ei] == owner[ej]);
                    }
                    if (first) emit(send, so, owner[ej], el, kind, c % per, el * n3);
                }
            }
    }
    if (!FILL) {
        send_at[t] = so;
        recv_at[t] = ro;
    }
}

// Sorts the entries by (peer, element id, local node): the order both sides of every
// message agree on.
inline void sort_entries(thrust::device_vector<int>& peer, thrust::device_vector<unsigned long long>& gid,
                         thrust::device_vector<int>& node, thrust::device_vector<long long>& index)
{
    const long long n = (long long)peer.size();
    thrust::device_vector<long long> perm(n);
    thrust::sequence(thrust::device, perm.begin(), perm.end());
    thrust::device_vector<int> k32(node);
    thrust::stable_sort_by_key(thrust::device, k32.begin(), k32.end(), perm.begin());
    thrust::device_vector<unsigned long long> k64(n);
    thrust::gather(thrust::device, perm.begin(), perm.end(), gid.begin(), k64.begin());
    thrust::stable_sort_by_key(thrust::device, k64.begin(), k64.end(), perm.begin());
    thrust::gather(thrust::device, perm.begin(), perm.end(), peer.begin(), k32.begin());
    thrust::stable_sort_by_key(thrust::device, k32.begin(), k32.end(), perm.begin());
    peer.swap(k32);
    thrust::gather(thrust::device, perm.begin(), perm.end(), gid.begin(), k64.begin());
    gid.swap(k64);
    thrust::device_vector<int> sorted_node(n);
    thrust::gather(thrust::device, perm.begin(), perm.end(), node.begin(), sorted_node.begin());
    node.swap(sorted_node);
    thrust::device_vector<long long> sorted_index(n);
    thrust::gather(thrust::device, perm.begin(), perm.end(), index.begin(), sorted_index.begin());
    index.swap(sorted_index);
}

// out = b - q at the listed copies (b alone when q is null).
__global__ void pack_kernel(const double* __restrict__ b, const double* __restrict__ q,
                            const long long* __restrict__ index, long long n, double* __restrict__ out)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < n) out[t] = q ? b[index[t]] - q[index[t]] : b[index[t]];
}

__global__ void unpack_kernel(const double* __restrict__ in, const long long* __restrict__ index, long long n,
                              double* __restrict__ ghost)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < n) ghost[index[t]] = in[t];
}

// The values one item of entries_kernel sends: a star with L local copies and ghost copies
// of R ranks sends L R values per node now (every local copy to every such rank), and would
// send R with one summed value per (rank, node). An interface face sends (p - 1)^2 either way.
__global__ void send_count_kernel(TopologyView T, long long edges, long long vertices, const int* __restrict__ owner,
                                  int p, long long* __restrict__ now, long long* __restrict__ summed)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    const long long faces = T.local * 6;
    if (t >= faces + edges + vertices) return;
    if (t < faces) {
        const int c = T.face_nbr[t];
        const long long v = c >= 0 && c / 6 >= T.local ? entity_nodes(1, p) : 0;
        now[t] = v;
        summed[t] = v;
        return;
    }
    const bool edge = t < faces + edges;
    const long long g = edge ? t - faces : t - faces - edges;
    const int* off = edge ? T.edge_off : T.vert_off;
    const int* ent = edge ? T.edge_ent : T.vert_ent;
    const int per = edge ? 12 : 8;
    const int b = off[g], end = off[g + 1];
    long long L = 0, R = 0;
    for (int j = b; j < end; ++j) {
        const long long ej = (ent[j] & 0x7fffffff) / per;
        if (ej < T.local) {
            ++L;
            continue;
        }
        bool first = true;
        for (int i = b; i < j && first; ++i) {
            const long long ei = (ent[i] & 0x7fffffff) / per;
            first = !(ei >= T.local && owner[ei] == owner[ej]);
        }
        R += first;
    }
    const long long n = entity_nodes(edge ? 2 : 3, p);
    now[t] = L * R * n;
    summed[t] = L > 0 ? R * n : 0;
}

// need[g] = 1 for every ghost value g that some copy reads through src.
__global__ void need_kernel(const int* __restrict__ src, long long copies, long long local_values,
                            int* __restrict__ need)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < copies && src[t] >= local_values) atomicExch(need + src[t] - local_values, 1);
}

// A functor: nvcc allows no device lambda inside a constructor.
struct NonZero {
    __host__ __device__ bool operator()(int k) const { return k != 0; }
};

__global__ void keep_kernel(const int* __restrict__ need, const long long* __restrict__ index, long long n,
                            int* __restrict__ keep)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < n) keep[t] = need[index[t]];
}

// The entries of the counted-sum exchange (SummedHalo), one item per local face, edge star
// and vertex star. A send entry carries the sender's own copies of one node (sources, in
// canonical order); a receive entry the ghost slot its value goes to and the slots that get
// zero. Keys (element id, local node) are those of the counting element for stars and of the
// sending element for faces, so both sides list them alike.
struct SummedSend {
    int* peer;
    unsigned long long* gid;
    int* node;
    long long* begin;   // into sources
    int* len;
    long long* sources;
};
struct SummedRecv {
    int* peer;
    unsigned long long* gid;
    int* node;
    long long* target;
    long long* begin;   // into zeros
    int* len;
    long long* zeros;
};

template <bool FILL>
__global__ void summed_entries_kernel(TopologyView T, long long edges, long long vertices,
                                      const unsigned long long* __restrict__ gid, const int* __restrict__ owner,
                                      int p, int me, long long* __restrict__ send_at, long long* __restrict__ src_at,
                                      long long* __restrict__ recv_at, long long* __restrict__ zero_at,
                                      SummedSend send, SummedRecv recv)
{
    const long long n3 = (long long)(p + 1) * (p + 1) * (p + 1);
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    const long long faces = T.local * 6;
    if (t >= faces + edges + vertices) return;
    long long so = FILL ? send_at[t] : 0, sso = FILL ? src_at[t] : 0, ro = FILL ? recv_at[t] : 0,
              zo = FILL ? zero_at[t] : 0;
    auto slot = [&](long long elem, int l) {
        return elem < T.local ? elem * n3 + l : (elem - T.local) * n3 + l;
    };
    if (t < faces) {
        const long long e = t / 6;
        const int f = (int)(t % 6), c = T.face_nbr[t];
        if (c >= 0 && c / 6 >= T.local) {
            const long long ghost = c / 6;
            const int n = entity_nodes(1, p);
            const bool counted = (T.counting[e] >> f) & 1u;
            for (int j = 0; j < n; ++j) {
                if (counted) {
                    if (FILL) {
                        const int l = entity_node(1, c % 6, j, p);
                        recv.peer[ro] = owner[ghost];
                        recv.gid[ro] = gid[ghost];
                        recv.node[ro] = l;
                        recv.target[ro] = slot(ghost, l);
                        recv.begin[ro] = zo;
                        recv.len[ro] = 0;
                    }
                    ++ro;
                } else {
                    if (FILL) {
                        const int l = entity_node(1, f, j, p);
                        send.peer[so] = owner[ghost];
                        send.gid[so] = gid[e];
                        send.node[so] = l;
                        send.begin[so] = sso;
                        send.len[so] = 1;
                        send.sources[sso] = slot(e, l);
                    }
                    ++so;
                    ++sso;
                }
            }
        }
    } else {
        const bool edge = t < faces + edges;
        const long long g = edge ? t - faces : t - faces - edges;
        const int* off = edge ? T.edge_off : T.vert_off;
        const int* ent = edge ? T.edge_ent : T.vert_ent;
        const int per = edge ? 12 : 8, kind = edge ? 2 : 3;
        const int b = off[g], end = off[g + 1];
        bool any_local = false, any_ghost = false;
        for (int m = b; m < end; ++m) ((ent[m] & 0x7fffffff) / per < T.local ? any_local : any_ghost) = true;
        if (any_local && any_ghost) {
            const int first = ent[b], fc = first & 0x7fffffff;
            const long long fe = fc / per;
            const int counter = fe < T.local ? me : owner[fe];
            const int n = entity_nodes(kind, p);
            for (int j = 0; j < n; ++j) {
                // Node j in the counting element's frame, at canonical position tc; an entry
                // with the reversed bit holds it at p - 2 - tc (edge interiors, 0-based).
                const int tc = edge && first < 0 ? p - 2 - j : j;
                auto local_node = [&](int c) {
                    return entity_node(kind, (c & 0x7fffffff) % per, edge && c < 0 ? p - 2 - tc : tc, p);
                };
                const int key_node = entity_node(kind, fc % per, j, p);
                if (counter == me) {
                    for (int m = b; m < end; ++m) {
                        const long long em = (ent[m] & 0x7fffffff) / per;
                        if (em < T.local) continue;
                        bool head = true;   // the first ghost copy of its owner
                        for (int i = b; i < m && head; ++i) {
                            const long long ei = (ent[i] & 0x7fffffff) / per;
                            head = !(ei >= T.local && owner[ei] == owner[em]);
                        }
                        if (!head) continue;
                        if (FILL) {
                            recv.peer[ro] = owner[em];
                            recv.gid[ro] = gid[fe];
                            recv.node[ro] = key_node;
                            recv.target[ro] = slot(em, local_node(ent[m]));
                            recv.begin[ro] = zo;
                        }
                        int zeros = 0;
                        for (int i = m + 1; i < end; ++i) {
                            const long long ei = (ent[i] & 0x7fffffff) / per;
                            if (ei < T.local || owner[ei] != owner[em]) continue;
                            if (FILL) recv.zeros[zo + zeros] = slot(ei, local_node(ent[i]));
                            ++zeros;
                        }
                        if (FILL) recv.len[ro] = zeros;
                        ++ro;
                        zo += zeros;
                    }
                } else {
                    if (FILL) {
                        send.peer[so] = counter;
                        send.gid[so] = gid[fe];
                        send.node[so] = key_node;
                        send.begin[so] = sso;
                    }
                    int count = 0;
                    for (int m = b; m < end; ++m) {
                        const long long em = (ent[m] & 0x7fffffff) / per;
                        if (em >= T.local) continue;
                        if (FILL) send.sources[sso + count] = slot(em, local_node(ent[m]));
                        ++count;
                    }
                    if (FILL) send.len[so] = count;
                    ++so;
                    sso += count;
                }
            }
        }
    }
    if (!FILL) {
        send_at[t] = so;
        src_at[t] = sso;
        recv_at[t] = ro;
        zero_at[t] = zo;
    }
}

// out[k] = the sum of b - q (b when q is null) over entry k's sources, in their order.
__global__ void summed_pack_kernel(const double* __restrict__ b, const double* __restrict__ q,
                                   const long long* __restrict__ begin, const int* __restrict__ len,
                                   const long long* __restrict__ sources, long long n, double* __restrict__ out)
{
    const long long k = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= n) return;
    double s = 0.0;
    for (int i = 0; i < len[k]; ++i) {
        const long long c = sources[begin[k] + i];
        s += q ? b[c] - q[c] : b[c];
    }
    out[k] = s;
}

__global__ void summed_unpack_kernel(const double* __restrict__ in, const long long* __restrict__ target,
                                     const long long* __restrict__ begin, const int* __restrict__ len,
                                     const long long* __restrict__ zeros, long long n, double* __restrict__ ghost)
{
    const long long k = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= n) return;
    ghost[target[k]] = in[k];
    for (int i = 0; i < len[k]; ++i) ghost[zeros[begin[k] + i]] = 0.0;
}

}  // namespace halo_detail

// The ghost exchange of one rank's tables for values of degree p. gid: a global id of each
// of the T.elements elements (any partition-independent, unique number); owner: the rank
// that owns it.
// start() sends the local copies the peers need, the caller runs the elements that read
// no ghost copies, finish() waits and files the received copies into the ghost blocks.
// Device buffers throughout: GPU-aware MPI (on Cray MPICH, MPICH_GPU_SUPPORT_ENABLED=1).
class UnstructuredHalo {
public:
    UnstructuredHalo(const UnstructuredTopology& T, const unsigned long long* d_gid, const int* d_owner,
                     MPI_Comm comm, int p = kP)
    {
        using namespace halo_detail;
        MARS_CELLWISE_MPI(MPI_Comm_dup(comm, &comm_));
        int rank = 0, size = 1;
        MARS_CELLWISE_MPI(MPI_Comm_rank(comm_, &rank));
        MARS_CELLWISE_MPI(MPI_Comm_size(comm_, &size));
        const TopologyView v = view(T);
        const long long L = T.local;

        thrust::device_vector<int> flag(L, 0);
        if (L > 0) boundary_kernel<<<blocks(L), kThreads>>>(v, raw(flag));
        interior_ids_.resize(L);
        boundary_ids_.resize(L);
        const auto ids = thrust::counting_iterator<int>(0);
        const long long ni =
            thrust::copy_if(thrust::device, ids, ids + L, interior_ids_.begin(), FlagIs{raw(flag), false}) -
            interior_ids_.begin();
        const long long nb =
            thrust::copy_if(thrust::device, ids, ids + L, boundary_ids_.begin(), FlagIs{raw(flag), true}) -
            boundary_ids_.begin();
        interior_ = {raw(interior_ids_), ni};
        boundary_ = {raw(boundary_ids_), nb};

        // Count, scan, fill.
        const long long items = L * 6 + T.edges() + T.vertices();
        thrust::device_vector<long long> send_at(items), recv_at(items);
        Entries none{nullptr, nullptr, nullptr, nullptr};
        if (items > 0)
            entries_kernel<false><<<blocks(items), kThreads>>>(v, T.edges(), T.vertices(), d_gid, d_owner, p,
                                                               raw(send_at), raw(recv_at), none, none);
        MARS_CELLWISE_CK(cudaGetLastError());
        const long long send_total = thrust::reduce(thrust::device, send_at.begin(), send_at.end(), 0LL);
        const long long recv_total = thrust::reduce(thrust::device, recv_at.begin(), recv_at.end(), 0LL);
        long long counted = 0;
        {
            // The same send total from the star sizes, which checks the count formula that the
            // summed total comes from.
            thrust::device_vector<long long> now(items), summed(items);
            if (items > 0)
                send_count_kernel<<<blocks(items), kThreads>>>(v, T.edges(), T.vertices(), d_owner, p, raw(now),
                                                               raw(summed));
            MARS_CELLWISE_CK(cudaGetLastError());
            counted = mars::reduce(now.begin(), now.end(), 0LL, thrust::plus<long long>());
            summed_total_ = mars::reduce(summed.begin(), summed.end(), 0LL, thrust::plus<long long>());
        }
        thrust::exclusive_scan(thrust::device, send_at.begin(), send_at.end(), send_at.begin());
        thrust::exclusive_scan(thrust::device, recv_at.begin(), recv_at.end(), recv_at.begin());
        thrust::device_vector<int> speer(send_total), rpeer(recv_total), snode(send_total), rnode(recv_total);
        thrust::device_vector<unsigned long long> sgid(send_total), rgid(recv_total);
        send_index_.resize(send_total);
        recv_index_.resize(recv_total);
        if (items > 0)
            entries_kernel<true><<<blocks(items), kThreads>>>(
                v, T.edges(), T.vertices(), d_gid, d_owner, p, raw(send_at), raw(recv_at),
                Entries{raw(speer), raw(sgid), raw(snode), raw(send_index_)},
                Entries{raw(rpeer), raw(rgid), raw(rnode), raw(recv_index_)});
        MARS_CELLWISE_CK(cudaGetLastError());
        sort_entries(speer, sgid, snode, send_index_);
        sort_entries(rpeer, rgid, rnode, recv_index_);

        // Per-rank counts; the two sides of every message must agree on them.
        std::vector<long long> scount(size, 0), rcount(size, 0), expect(size, 0);
        count_by_peer(speer, scount);
        count_by_peer(rpeer, rcount);
        MARS_CELLWISE_MPI(MPI_Alltoall(scount.data(), 1, MPI_LONG_LONG, expect.data(), 1, MPI_LONG_LONG, comm_));
        int bad = 0;
        for (int p = 0; p < size; ++p) bad |= expect[p] != rcount[p];
        if (counted != send_total) bad |= 2;
        long long sofs = 0, rofs = 0;
        for (int p = 0; p < size; ++p) {
            if (scount[p] == 0 && rcount[p] == 0) continue;
            peers_.push_back({p, sofs, scount[p], rofs, rcount[p]});
            sofs += scount[p];
            rofs += rcount[p];
        }
        MARS_CELLWISE_MPI(MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_INT, MPI_BOR, comm_));
        if (bad) {
            if (rank == 0 && (bad & 1)) fprintf(stderr, "cell-wise halo: send and receive counts disagree between ranks\n");
            if (rank == 0 && (bad & 2)) fprintf(stderr, "cell-wise halo: the star sizes do not give the send count\n");
            MPI_Abort(comm_, 1);
        }
        requests_.reserve(2 * peers_.size());

        if (send_total > 0) MARS_CELLWISE_CK(cudaMalloc(&send_, send_total * sizeof(double)));
        if (recv_total > 0) MARS_CELLWISE_CK(cudaMalloc(&recv_, recv_total * sizeof(double)));
        const long long ghost_values = (T.elements - L) * (p + 1) * (p + 1) * (p + 1);
        local_values_ = L * (p + 1) * (p + 1) * (p + 1);
        ghost_values_ = ghost_values;
        if (ghost_values > 0) {
            MARS_CELLWISE_CK(cudaMalloc(&ghost_, ghost_values * sizeof(double)));
            MARS_CELLWISE_CK(cudaMemset(ghost_, 0, ghost_values * sizeof(double)));
        }
        send_total_ = send_total;
        recv_total_ = recv_total;
        check_order(sgid, snode, rgid, rnode, rank);
    }
    ~UnstructuredHalo()
    {
        cudaFree(send_);
        cudaFree(recv_);
        cudaFree(ghost_);
        MPI_Comm_free(&comm_);
    }
    UnstructuredHalo(const UnstructuredHalo&) = delete;
    UnstructuredHalo& operator=(const UnstructuredHalo&) = delete;

    bool active() const { return !peers_.empty(); }
    const ElementList& interior() const { return interior_; }
    const ElementList& boundary() const { return boundary_; }
    const double* ghost() const { return ghost_; }
    long long values_sent() const { return send_total_; }
    // What values_sent() would be with one summed value per (peer, shared node).
    long long values_sent_summed() const { return summed_total_; }
    int peers() const { return (int)peers_.size(); }

    // Overwrites every ghost value with the bytes `byte` (a test poisons the slots that a
    // narrower exchange leaves alone).
    void fill_ghost(int byte, cudaStream_t stream) const
    {
        if (ghost_values_ > 0) MARS_CELLWISE_CK(cudaMemsetAsync(ghost_, byte, ghost_values_ * sizeof(double), stream));
    }

    // Sends b - q (b when q is null) at the copies the peers need.
    void start(const double* b, const double* q, cudaStream_t stream)
    {
        using namespace halo_detail;
        if (send_total_ > 0)
            pack_kernel<<<blocks(send_total_), kThreads, 0, stream>>>(b, q, raw(send_index_), send_total_, send_);
        MARS_CELLWISE_CK(cudaGetLastError());
        MARS_CELLWISE_CK(cudaStreamSynchronize(stream));   // MPI reads send_ next
        post(send_, recv_);
    }
    void finish(cudaStream_t stream)
    {
        using namespace halo_detail;
        MARS_CELLWISE_MPI(MPI_Waitall((int)requests_.size(), requests_.data(), MPI_STATUSES_IGNORE));
        if (recv_total_ > 0)
            unpack_kernel<<<blocks(recv_total_), kThreads, 0, stream>>>(recv_, raw(recv_index_), recv_total_, ghost_);
        MARS_CELLWISE_CK(cudaGetLastError());
    }

private:
    friend class RequestedHalo;
    friend class SummedHalo;
    struct Peer {
        int rank;
        long long send_offset, send_count, recv_offset, recv_count;
    };
    static constexpr int kTag = 0x3300;   // clear of the structured halo (0x3100) and other MARS tags

    template <typename T>
    void post(const T* send, T* recv)
    {
        const MPI_Datatype type = halo_detail::mpi_type<T>();
        requests_.clear();
        for (const Peer& p : peers_)
            if (p.recv_count > 0) {
                requests_.emplace_back();
                MARS_CELLWISE_MPI(MPI_Irecv(recv + p.recv_offset, (int)p.recv_count, type, p.rank, kTag, comm_,
                                            &requests_.back()));
            }
        for (const Peer& p : peers_)
            if (p.send_count > 0) {
                requests_.emplace_back();
                MARS_CELLWISE_MPI(MPI_Isend(send + p.send_offset, (int)p.send_count, type, p.rank, kTag, comm_,
                                            &requests_.back()));
            }
    }

    static void count_by_peer(const thrust::device_vector<int>& peer, std::vector<long long>& count)
    {
        thrust::device_vector<int> rank(peer.size());
        thrust::device_vector<long long> n(peer.size());
        const auto end = thrust::reduce_by_key(thrust::device, peer.begin(), peer.end(),
                                               thrust::constant_iterator<long long>(1), rank.begin(), n.begin());
        const long long k = end.first - rank.begin();
        std::vector<int> h_rank(k);
        std::vector<long long> h_n(k);
        thrust::copy(rank.begin(), rank.begin() + k, h_rank.begin());
        thrust::copy(n.begin(), n.begin() + k, h_n.begin());
        for (long long i = 0; i < k; ++i) count[h_rank[i]] = h_n[i];
    }

    // Sends each peer the (element id, node) keys in send order; they must equal the keys
    // the peer expects, in its receive order.
    void check_order(const thrust::device_vector<unsigned long long>& sgid, const thrust::device_vector<int>& snode,
                     const thrust::device_vector<unsigned long long>& rgid, const thrust::device_vector<int>& rnode,
                     int rank)
    {
        thrust::device_vector<unsigned long long> sent(sgid.size()), expected(rgid.size());
        thrust::transform(thrust::device, sgid.begin(), sgid.end(), snode.begin(), sent.begin(), halo_detail::KeyHash{});
        thrust::transform(thrust::device, rgid.begin(), rgid.end(), rnode.begin(), expected.begin(),
                          halo_detail::KeyHash{});
        thrust::device_vector<unsigned long long> got(expected.size());
        post(thrust::raw_pointer_cast(sent.data()), thrust::raw_pointer_cast(got.data()));
        MARS_CELLWISE_MPI(MPI_Waitall((int)requests_.size(), requests_.data(), MPI_STATUSES_IGNORE));
        int bad = !thrust::equal(thrust::device, got.begin(), got.end(), expected.begin());
        MARS_CELLWISE_MPI(MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_INT, MPI_MAX, comm_));
        if (bad) {
            if (rank == 0) fprintf(stderr, "cell-wise halo: the ranks list their shared nodes in different orders\n");
            MPI_Abort(comm_, 1);
        }
    }

    MPI_Comm comm_ = MPI_COMM_NULL;
    std::vector<Peer> peers_;
    thrust::device_vector<long long> send_index_, recv_index_;
    long long send_total_ = 0, recv_total_ = 0, summed_total_ = 0, local_values_ = 0, ghost_values_ = 0;
    double *send_ = nullptr, *recv_ = nullptr, *ghost_ = nullptr;
    thrust::device_vector<int> interior_ids_, boundary_ids_;
    ElementList interior_, boundary_;
    std::vector<MPI_Request> requests_;
};

// The exchange of only the ghost copies a consumer reads through src (src[t] >= local values:
// ghost value src[t] - local values), on the peers, order and ghost slots of `full`. At setup
// each rank tells every peer which of the copies that peer sends it are read, and both sides
// keep those entries. The other ghost slots keep what an earlier exchange left there, so the
// consumer must read nothing but src. Used by the node-line update, which reads one counting
// copy per node where the full exchange sends every copy.
class RequestedHalo {
public:
    RequestedHalo(UnstructuredHalo& full, const int* d_src, long long copies, cudaStream_t stream) : full_(full)
    {
        using namespace halo_detail;
        int rank = 0;
        MARS_CELLWISE_MPI(MPI_Comm_rank(full.comm_, &rank));
        MARS_CELLWISE_CK(cudaStreamSynchronize(stream));   // src may still be being written
        thrust::device_vector<int> need(full.ghost_values_, 0), keep_send(full.send_total_),
            keep_recv(full.recv_total_);
        if (copies > 0)
            need_kernel<<<blocks(copies), kThreads>>>(d_src, copies, full.local_values_, raw(need));
        if (full.recv_total_ > 0)
            keep_kernel<<<blocks(full.recv_total_), kThreads>>>(raw(need), raw(full.recv_index_), full.recv_total_,
                                                                raw(keep_recv));
        MARS_CELLWISE_CK(cudaGetLastError());
        MARS_CELLWISE_CK(cudaDeviceSynchronize());   // MPI reads keep_recv next

        // Each peer learns which of its sends are read: the reverse direction of the messages.
        std::vector<MPI_Request> req;
        for (const auto& p : full.peers_) {
            if (p.send_count > 0) {
                req.emplace_back();
                MARS_CELLWISE_MPI(MPI_Irecv(raw(keep_send) + p.send_offset, (int)p.send_count, MPI_INT, p.rank,
                                            kMaskTag, full.comm_, &req.back()));
            }
            if (p.recv_count > 0) {
                req.emplace_back();
                MARS_CELLWISE_MPI(MPI_Isend(raw(keep_recv) + p.recv_offset, (int)p.recv_count, MPI_INT, p.rank,
                                            kMaskTag, full.comm_, &req.back()));
            }
        }
        MARS_CELLWISE_MPI(MPI_Waitall((int)req.size(), req.data(), MPI_STATUSES_IGNORE));

        send_index_.resize(full.send_total_);
        recv_index_.resize(full.recv_total_);
        send_total_ = thrust::copy_if(thrust::device, full.send_index_.begin(), full.send_index_.end(),
                                      keep_send.begin(), send_index_.begin(), NonZero{}) -
                      send_index_.begin();
        recv_total_ = thrust::copy_if(thrust::device, full.recv_index_.begin(), full.recv_index_.end(),
                                      keep_recv.begin(), recv_index_.begin(), NonZero{}) -
                      recv_index_.begin();
        send_index_.resize(send_total_);
        recv_index_.resize(recv_total_);
        long long sofs = 0, rofs = 0;
        for (const auto& p : full.peers_) {
            const long long s = mars::count(keep_send.begin() + p.send_offset,
                                            keep_send.begin() + p.send_offset + p.send_count, 1);
            const long long r = mars::count(keep_recv.begin() + p.recv_offset,
                                            keep_recv.begin() + p.recv_offset + p.recv_count, 1);
            peers_.push_back({p.rank, sofs, s, rofs, r});
            sofs += s;
            rofs += r;
        }

        // Every ghost value src reads must arrive exactly once, and each pair of ranks must
        // agree on the new message sizes.
        const long long needed = mars::count(need.begin(), need.end(), 1);
        int bad = needed != recv_total_;
        std::vector<long long> told(peers_.size());
        req.clear();
        for (size_t i = 0; i < peers_.size(); ++i) {
            req.emplace_back();
            MARS_CELLWISE_MPI(MPI_Irecv(&told[i], 1, MPI_LONG_LONG, peers_[i].rank, kMaskTag + 1, full.comm_,
                                        &req.back()));
            req.emplace_back();
            MARS_CELLWISE_MPI(MPI_Isend(&peers_[i].send_count, 1, MPI_LONG_LONG, peers_[i].rank, kMaskTag + 1,
                                        full.comm_, &req.back()));
        }
        MARS_CELLWISE_MPI(MPI_Waitall((int)req.size(), req.data(), MPI_STATUSES_IGNORE));
        for (size_t i = 0; i < peers_.size(); ++i) bad |= (told[i] != peers_[i].recv_count) << 1;
        MARS_CELLWISE_MPI(MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_INT, MPI_BOR, full.comm_));
        if (bad) {
            if (rank == 0 && (bad & 1)) fprintf(stderr, "requested halo: a ghost copy src reads is not exchanged\n");
            if (rank == 0 && (bad & 2)) fprintf(stderr, "requested halo: the ranks disagree on message sizes\n");
            MPI_Abort(full.comm_, 1);
        }
        if (send_total_ > 0) MARS_CELLWISE_CK(cudaMalloc(&send_, send_total_ * sizeof(double)));
        if (recv_total_ > 0) MARS_CELLWISE_CK(cudaMalloc(&recv_, recv_total_ * sizeof(double)));
        requests_.reserve(2 * peers_.size());
    }
    ~RequestedHalo()
    {
        cudaFree(send_);
        cudaFree(recv_);
    }
    RequestedHalo(const RequestedHalo&) = delete;
    RequestedHalo& operator=(const RequestedHalo&) = delete;

    bool active() const { return full_.active(); }
    const ElementList& interior() const { return full_.interior(); }
    const ElementList& boundary() const { return full_.boundary(); }
    const double* ghost() const { return full_.ghost_; }
    long long values_sent() const { return send_total_; }

    void fill_ghost(int byte, cudaStream_t stream) const { full_.fill_ghost(byte, stream); }

    void start(const double* b, const double* q, cudaStream_t stream)
    {
        using namespace halo_detail;
        if (send_total_ > 0)
            pack_kernel<<<blocks(send_total_), kThreads, 0, stream>>>(b, q, raw(send_index_), send_total_, send_);
        MARS_CELLWISE_CK(cudaGetLastError());
        MARS_CELLWISE_CK(cudaStreamSynchronize(stream));   // MPI reads send_ next
        requests_.clear();
        for (const auto& p : peers_)
            if (p.recv_count > 0) {
                requests_.emplace_back();
                MARS_CELLWISE_MPI(MPI_Irecv(recv_ + p.recv_offset, (int)p.recv_count, MPI_DOUBLE, p.rank, kTag,
                                            full_.comm_, &requests_.back()));
            }
        for (const auto& p : peers_)
            if (p.send_count > 0) {
                requests_.emplace_back();
                MARS_CELLWISE_MPI(MPI_Isend(send_ + p.send_offset, (int)p.send_count, MPI_DOUBLE, p.rank, kTag,
                                            full_.comm_, &requests_.back()));
            }
    }
    void finish(cudaStream_t stream)
    {
        using namespace halo_detail;
        MARS_CELLWISE_MPI(MPI_Waitall((int)requests_.size(), requests_.data(), MPI_STATUSES_IGNORE));
        if (recv_total_ > 0)
            unpack_kernel<<<blocks(recv_total_), kThreads, 0, stream>>>(recv_, raw(recv_index_), recv_total_,
                                                                        full_.ghost_);
        MARS_CELLWISE_CK(cudaGetLastError());
    }

private:
    static constexpr int kTag = 0x3301, kMaskTag = 0x3302;   // next to the full exchange's 0x3300
    UnstructuredHalo& full_;
    std::vector<UnstructuredHalo::Peer> peers_;
    thrust::device_vector<long long> send_index_, recv_index_;
    long long send_total_ = 0, recv_total_ = 0;
    double *send_ = nullptr, *recv_ = nullptr;
    std::vector<MPI_Request> requests_;
};

// The residual exchange of the node-line smoother. Only the rank that counts a shared node
// (the rank of its canonically first copy) uses that node's sum: the residual is read at
// counting copies alone. So every other rank that shares the node sends that rank one value,
// the sum of its own copies of b - q in canonical order, where the full exchange sends every
// copy to every rank holding a ghost copy. The counting rank files the value into the
// sender's first ghost copy in the node's star and zero into the sender's other ghost copies,
// so the gather adds the same terms grouped by rank: the sum changes in its rounding only. A
// face node has one copy per side, so only the non-counting side sends it, and its sum is
// unchanged. Sums at non-counting copies read stale or zero ghost values: only consumers
// that use counting copies alone may run on this exchange.
class SummedHalo {
public:
    SummedHalo(UnstructuredHalo& full, const UnstructuredTopology& T, const unsigned long long* d_gid,
               const int* d_owner, int p)
        : full_(full)
    {
        using namespace halo_detail;
        int rank = 0;
        MARS_CELLWISE_MPI(MPI_Comm_rank(full.comm_, &rank));
        const TopologyView v = view(T);
        const long long items = T.local * 6 + T.edges() + T.vertices();
        thrust::device_vector<long long> send_at(items), src_at(items), recv_at(items), zero_at(items);
        const SummedSend no_send{};
        const SummedRecv no_recv{};
        if (items > 0)
            summed_entries_kernel<false><<<blocks(items), kThreads>>>(v, T.edges(), T.vertices(), d_gid, d_owner, p,
                                                                      rank, raw(send_at), raw(src_at), raw(recv_at),
                                                                      raw(zero_at), no_send, no_recv);
        MARS_CELLWISE_CK(cudaGetLastError());
        const auto plus = thrust::plus<long long>();
        send_total_ = mars::reduce(send_at.begin(), send_at.end(), 0LL, plus);
        const long long sources = mars::reduce(src_at.begin(), src_at.end(), 0LL, plus);
        recv_total_ = mars::reduce(recv_at.begin(), recv_at.end(), 0LL, plus);
        const long long zeros = mars::reduce(zero_at.begin(), zero_at.end(), 0LL, plus);
        for (auto* at : {&send_at, &src_at, &recv_at, &zero_at})
            thrust::exclusive_scan(thrust::device, at->begin(), at->end(), at->begin());
        thrust::device_vector<int> speer(send_total_), snode(send_total_), rpeer(recv_total_), rnode(recv_total_);
        thrust::device_vector<unsigned long long> sgid(send_total_), rgid(recv_total_);
        thrust::device_vector<long long> sbegin(send_total_), rtarget(recv_total_), rbegin(recv_total_);
        thrust::device_vector<int> slen(send_total_), rlen(recv_total_);
        sources_.resize(sources);
        zeros_.resize(zeros);
        if (items > 0)
            summed_entries_kernel<true><<<blocks(items), kThreads>>>(
                v, T.edges(), T.vertices(), d_gid, d_owner, p, rank, raw(send_at), raw(src_at), raw(recv_at),
                raw(zero_at), SummedSend{raw(speer), raw(sgid), raw(snode), raw(sbegin), raw(slen), raw(sources_)},
                SummedRecv{raw(rpeer), raw(rgid), raw(rnode), raw(rtarget), raw(rbegin), raw(rlen), raw(zeros_)});
        MARS_CELLWISE_CK(cudaGetLastError());

        // Both sides in (peer, element id, node) order; the entries' other fields follow.
        auto order = [](thrust::device_vector<int>& peer, thrust::device_vector<unsigned long long>& gid,
                        thrust::device_vector<int>& node) {
            thrust::device_vector<long long> perm(peer.size());
            thrust::sequence(thrust::device, perm.begin(), perm.end());
            sort_entries(peer, gid, node, perm);
            return perm;
        };
        auto permuted = [](const thrust::device_vector<long long>& perm, auto& v) {
            std::remove_reference_t<decltype(v)> out(v.size());
            thrust::gather(thrust::device, perm.begin(), perm.end(), v.begin(), out.begin());
            v.swap(out);
        };
        const thrust::device_vector<long long> sperm = order(speer, sgid, snode);
        permuted(sperm, sbegin);
        permuted(sperm, slen);
        const thrust::device_vector<long long> rperm = order(rpeer, rgid, rnode);
        permuted(rperm, rtarget);
        permuted(rperm, rbegin);
        permuted(rperm, rlen);
        send_begin_.swap(sbegin);
        send_len_.swap(slen);
        recv_target_.swap(rtarget);
        recv_begin_.swap(rbegin);
        recv_len_.swap(rlen);

        // Counts use the existing sparse peers, including peers whose new count is zero.
        int bad = 0;
        auto counts = [&](const thrust::device_vector<int>& peer) {
            std::vector<long long> count(full.peers_.size(), 0);
            thrust::device_vector<int> ranks(peer.size());
            thrust::device_vector<long long> values(peer.size());
            const auto end = thrust::reduce_by_key(thrust::device, peer.begin(), peer.end(),
                                                   thrust::constant_iterator<long long>(1), ranks.begin(),
                                                   values.begin());
            const long long n = end.first - ranks.begin();
            std::vector<int> h_rank(n);
            std::vector<long long> h_count(n);
            thrust::copy(ranks.begin(), ranks.begin() + n, h_rank.begin());
            thrust::copy(values.begin(), values.begin() + n, h_count.begin());
            for (long long i = 0; i < n; ++i) {
                const auto it = std::lower_bound(full.peers_.begin(), full.peers_.end(), h_rank[i],
                                                 [](const auto& p, int q) { return p.rank < q; });
                if (it == full.peers_.end() || it->rank != h_rank[i] || h_count[i] > INT_MAX) bad = 1;
                else count[it - full.peers_.begin()] = h_count[i];
            }
            return count;
        };
        const auto scount = counts(speer), rcount = counts(rpeer);
        MARS_CELLWISE_MPI(MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_INT, MPI_BOR, full.comm_));
        if (bad) {
            if (rank == 0) fprintf(stderr, "summed halo: invalid peer or message exceeds the MPI count limit\n");
            MPI_Abort(full.comm_, 1);
        }
        std::vector<long long> expect(full.peers_.size());
        std::vector<MPI_Request> req;
        req.reserve(2 * full.peers_.size());
        for (size_t i = 0; i < full.peers_.size(); ++i) {
            req.emplace_back();
            MARS_CELLWISE_MPI(MPI_Irecv(&expect[i], 1, MPI_LONG_LONG, full.peers_[i].rank, kCountTag,
                                        full.comm_, &req.back()));
            req.emplace_back();
            MARS_CELLWISE_MPI(MPI_Isend(&scount[i], 1, MPI_LONG_LONG, full.peers_[i].rank, kCountTag,
                                        full.comm_, &req.back()));
        }
        MARS_CELLWISE_MPI(MPI_Waitall((int)req.size(), req.data(), MPI_STATUSES_IGNORE));
        long long sofs = 0, rofs = 0;
        for (size_t i = 0; i < full.peers_.size(); ++i) {
            bad |= expect[i] != rcount[i];
            if (scount[i] == 0 && rcount[i] == 0) continue;
            peers_.push_back({full.peers_[i].rank, sofs, scount[i], rofs, rcount[i]});
            sofs += scount[i];
            rofs += rcount[i];
        }
        MARS_CELLWISE_MPI(MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_INT, MPI_BOR, full.comm_));
        if (bad) {
            if (rank == 0) fprintf(stderr, "summed halo: send and receive counts disagree between ranks\n");
            MPI_Abort(full.comm_, 1);
        }
        requests_.reserve(2 * peers_.size());

        // The keys in send order must equal the keys the peer expects in its receive order.
        thrust::device_vector<unsigned long long> sent(send_total_), expected(recv_total_), got(recv_total_);
        thrust::transform(thrust::device, sgid.begin(), sgid.end(), snode.begin(), sent.begin(), KeyHash{});
        thrust::transform(thrust::device, rgid.begin(), rgid.end(), rnode.begin(), expected.begin(), KeyHash{});
        MARS_CELLWISE_CK(cudaDeviceSynchronize());   // MPI reads sent next
        post(raw(sent), raw(got), MPI_UINT64_T, kKeyTag);
        MARS_CELLWISE_MPI(MPI_Waitall((int)requests_.size(), requests_.data(), MPI_STATUSES_IGNORE));
        bad = !thrust::equal(thrust::device, got.begin(), got.end(), expected.begin());
        MARS_CELLWISE_MPI(MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_INT, MPI_BOR, full.comm_));
        if (bad) {
            if (rank == 0) fprintf(stderr, "summed halo: the ranks list their shared nodes in different orders\n");
            MPI_Abort(full.comm_, 1);
        }
        if (send_total_ > 0) MARS_CELLWISE_CK(cudaMalloc(&send_, send_total_ * sizeof(double)));
        if (recv_total_ > 0) MARS_CELLWISE_CK(cudaMalloc(&recv_, recv_total_ * sizeof(double)));
    }
    ~SummedHalo()
    {
        cudaFree(send_);
        cudaFree(recv_);
    }
    SummedHalo(const SummedHalo&) = delete;
    SummedHalo& operator=(const SummedHalo&) = delete;

    bool active() const { return full_.active(); }
    const ElementList& interior() const { return full_.interior(); }
    const ElementList& boundary() const { return full_.boundary(); }
    const double* ghost() const { return full_.ghost_; }
    long long values_sent() const { return send_total_; }
    void fill_ghost(int byte, cudaStream_t stream) const { full_.fill_ghost(byte, stream); }

    void start(const double* b, const double* q, cudaStream_t stream)
    {
        using namespace halo_detail;
        if (send_total_ > 0)
            summed_pack_kernel<<<blocks(send_total_), kThreads, 0, stream>>>(
                b, q, raw(send_begin_), raw(send_len_), raw(sources_), send_total_, send_);
        MARS_CELLWISE_CK(cudaGetLastError());
        MARS_CELLWISE_CK(cudaStreamSynchronize(stream));   // MPI reads send_ next
        post(send_, recv_, MPI_DOUBLE, kTag);
    }
    void finish(cudaStream_t stream)
    {
        using namespace halo_detail;
        MARS_CELLWISE_MPI(MPI_Waitall((int)requests_.size(), requests_.data(), MPI_STATUSES_IGNORE));
        if (recv_total_ > 0)
            summed_unpack_kernel<<<blocks(recv_total_), kThreads, 0, stream>>>(
                recv_, raw(recv_target_), raw(recv_begin_), raw(recv_len_), raw(zeros_), recv_total_, full_.ghost_);
        MARS_CELLWISE_CK(cudaGetLastError());
    }

private:
    static constexpr int kTag = 0x3303, kKeyTag = 0x3304, kCountTag = 0x3305;   // next to the other cell-wise exchanges
    template <typename T>
    void post(const T* send, T* recv, MPI_Datatype type, int tag)
    {
        requests_.clear();
        for (const auto& q : peers_)
            if (q.recv_count > 0) {
                requests_.emplace_back();
                MARS_CELLWISE_MPI(MPI_Irecv(recv + q.recv_offset, (int)q.recv_count, type, q.rank, tag, full_.comm_,
                                            &requests_.back()));
            }
        for (const auto& q : peers_)
            if (q.send_count > 0) {
                requests_.emplace_back();
                MARS_CELLWISE_MPI(MPI_Isend(send + q.send_offset, (int)q.send_count, type, q.rank, tag, full_.comm_,
                                            &requests_.back()));
            }
    }

    UnstructuredHalo& full_;
    std::vector<UnstructuredHalo::Peer> peers_;
    thrust::device_vector<long long> sources_, zeros_, send_begin_, recv_target_, recv_begin_;
    thrust::device_vector<int> send_len_, recv_len_;
    long long send_total_ = 0, recv_total_ = 0;
    double *send_ = nullptr, *recv_ = nullptr;
    std::vector<MPI_Request> requests_;
};

// Runs a pass over the local elements on the values b - q (q may be null): one launch on
// one rank; otherwise start the exchange, launch the elements that read no ghost copies,
// wait, launch the rest. Halo: UnstructuredHalo, RequestedHalo or SummedHalo.
// launch(set, ghost, first) returns its number of thread blocks; `first` is how many
// blocks earlier launches used (for partial sums). Returns the total.
template <typename Halo, typename Launch>
int run_unstructured_pass(Halo* halo, long long local, const double* b, const double* q,
                          cudaStream_t stream, Launch&& launch)
{
    if (!halo || !halo->active()) return launch(ElementList{nullptr, local}, nullptr, 0);
    halo->start(b, q, stream);
    const int first = launch(halo->interior(), nullptr, 0);
    halo->finish(stream);
    return first + launch(halo->boundary(), halo->ghost(), first);
}

// out = DSS(in) on every local copy of degree P. `halo` (of the same degree) is null on
// one rank.
template <int P = kP>
void dss(const double* d_in, double* d_out, const UnstructuredTopology& T, UnstructuredHalo* halo = nullptr,
         cudaStream_t stream = 0)
{
    const TopologyView v = view(T);
    run_unstructured_pass(halo, T.local, d_in, nullptr, stream, [&](const ElementList& set, const double* ghost, int) {
        static const int grid = resident_grid(unstructured_dss_kernel<P>);
        unstructured_dss_kernel<P><<<grid, kThreads, 0, stream>>>(d_in, ghost, d_out, v, set);
        MARS_CELLWISE_CK(cudaGetLastError());
        return grid;
    });
}

// The Jacobi preconditioner of the Krylov solvers (same interface as
// JacobiPreconditioner in mars_cellwise_krylov.hpp). `halo` is null on one rank.
struct UnstructuredJacobi {
    TopologyView T;
    const double* diag;
    UnstructuredHalo* halo = nullptr;

    template <bool AZ, bool ZZ>
    int operator()(const double* q, double* z, const double* a, double* partial, cudaStream_t stream) const
    {
        constexpr int NV = AZ + ZZ;
        return run_unstructured_pass(halo, T.local, q, nullptr, stream, [&](const ElementList& set, const double* ghost,
                                                                   int first) {
            double* p = partial ? partial + (long long)first * (NV > 0 ? NV : 1) : nullptr;
            static const int grid = resident_grid(unstructured_precondition_kernel<AZ, ZZ>);
            unstructured_precondition_kernel<AZ, ZZ><<<grid, kThreads, 0, stream>>>(q, ghost, z, diag, a, T, set, p);
            MARS_CELLWISE_CK(cudaGetLastError());
            return grid;
        });
    }
};

// A structured block put through the unstructured tables on several ranks: global ids
// (lattice order of the NX x NY x NZ block) and owner ranks of this rank's elements
// (its sub-block of `dec`, in (x, y, z) order) followed by its ghost elements (the
// one-element frame around the sub-block, inside the block). Returns the local count.
// The lambdas are __host__ __device__: Thrust reads their result type on the host, and
// with a __device__-only lambda it picks a kernel that was never compiled.
inline long long block_elements(const Decomposition& dec, thrust::device_vector<unsigned long long>& gid,
                                thrust::device_vector<int>& owner)
{
    const Block b = dec.local;
    const long long L = b.elements();
    const int lo[3] = {std::max(b.ox - 1, 0), std::max(b.oy - 1, 0), std::max(b.oz - 1, 0)};
    const int hi[3] = {std::min(b.ox + b.nx + 1, b.NX), std::min(b.oy + b.ny + 1, b.NY),
                       std::min(b.oz + b.nz + 1, b.NZ)};
    const int span[3] = {hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]};
    const long long V = (long long)span[0] * span[1] * span[2];
    gid.resize(V);
    const auto it = thrust::counting_iterator<long long>(0);
    thrust::transform(thrust::device, it, it + L, gid.begin(), [b] __host__ __device__(long long e) {
        const long long z = e % b.nz + b.oz, y = (e / b.nz) % b.ny + b.oy, x = e / ((long long)b.nz * b.ny) + b.ox;
        return (unsigned long long)((x * b.NY + y) * b.NZ + z);
    });
    const int l0 = lo[0], l1 = lo[1], l2 = lo[2], s1 = span[1], s2 = span[2];
    const auto frame_gid = [=] __host__ __device__(long long i) {
        const long long z = i % s2 + l2, y = (i / s2) % s1 + l1, x = i / ((long long)s2 * s1) + l0;
        return (unsigned long long)((x * b.NY + y) * b.NZ + z);
    };
    const long long G = thrust::copy_if(
        thrust::device, thrust::make_transform_iterator(it, frame_gid), thrust::make_transform_iterator(it + V, frame_gid),
        gid.begin() + L, [b] __host__ __device__(unsigned long long g) {
            const long long z = g % b.NZ, y = (g / b.NZ) % b.NY, x = g / ((unsigned long long)b.NZ * b.NY);
            return x < b.ox || x >= b.ox + b.nx || y < b.oy || y >= b.oy + b.ny || z < b.oz || z >= b.oz + b.nz;
        }) - (gid.begin() + L);
    gid.resize(L + G);
    owner.resize(L + G);
    const int P1 = dec.P[1], P2 = dec.P[2];
    thrust::transform(thrust::device, gid.begin(), gid.end(), owner.begin(), [b, P1, P2] __host__ __device__(unsigned long long g) {
        const int z = (int)(g % b.NZ), y = (int)((g / b.NZ) % b.NY), x = (int)(g / ((unsigned long long)b.NZ * b.NY));
        return ((x / b.nx) * P1 + y / b.ny) * P2 + z / b.nz;
    });
    return L;
}

}  // namespace cellwise
}  // namespace mars
