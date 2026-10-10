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
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace mars {
namespace cellwise {

// ---- tables ----------------------------------------------------------------------------

// Entity bits of an element's 32-bit masks (counting copies, Dirichlet entities): face f
// at bit f, edge k at kEdgeBit + k, corner c at kVertexBit + c, interior nodes at
// kInteriorBit. The gather kernel loads each entity's descriptor on the lane with the
// same number, so one ballot gives the Dirichlet mask.
constexpr int kEdgeBit = 8, kVertexBit = 20, kInteriorBit = 31;

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

constexpr int kFaceNodes = (kN - 2) * (kN - 2);                  // 36 per face
constexpr int kEdgeNodes = kN - 2;                               // 6 per edge

// The local entity that holds node (a, b, c): kind 0 interior, 1 face, 2 edge, 3 corner,
// and its local number (f, k or c).
struct LocalEntity {
    int kind, index;
};

__device__ inline LocalEntity local_entity(int a, int b, int c)
{
    const int x[3] = {a, b, c};
    int on = 0;
    for (int k = 0; k < 3; ++k) on += (x[k] == 0 || x[k] == kP);
    if (on == 0) return {0, 0};
    if (on == 3) return {3, corner_from_bits(x[0] == kP, x[1] == kP, x[2] == kP)};
    int axis = 0, o0, o1;
    if (on == 1) {
        while (x[axis] != 0 && x[axis] != kP) ++axis;
        return {1, axis * 2 + (x[axis] == kP)};
    }
    while (x[axis] == 0 || x[axis] == kP) ++axis;   // the edge runs along the free axis
    other_axes(axis, o0, o1);
    return {2, axis * 4 + (x[o0] == kP) * 2 + (x[o1] == kP)};
}

__device__ inline int entity_bit(const LocalEntity& le)
{
    if (le.kind == 0) return kInteriorBit;
    return (le.kind == 1 ? 0 : (le.kind == 2 ? kEdgeBit : kVertexBit)) + le.index;
}

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

struct HexTables {
    unsigned short face_l[6 * kFaceNodes];   // local node of face f at position q = (i - 1) * 6 + j - 1
    unsigned char face_map[64 * kFaceNodes];  // frame codes (own | other << 3), own q -> the other's q
    unsigned short edge_l[12 * kN];           // local node of edge k at position t
    unsigned short corner_l[8];
    unsigned char bit[kN3];                   // entity bit of every local node
};

__device__ inline void build_hex_tables(HexTables& h)
{
    for (int i = threadIdx.x; i < 64 * kFaceNodes; i += blockDim.x) {
        const int q = i % kFaceNodes, code = i / kFaceNodes;
        int I, J, i2, j2;
        face_to_canonical(code & 7, 1 + q / 6, 1 + q % 6, I, J);
        canonical_to_face(code >> 3, I, J, i2, j2);
        h.face_map[i] = (unsigned char)((i2 - 1) * 6 + j2 - 1);
    }
    for (int i = threadIdx.x; i < 6 * kFaceNodes; i += blockDim.x)
        h.face_l[i] = (unsigned short)face_node(i / kFaceNodes, 1 + (i % kFaceNodes) / 6, 1 + i % 6);
    for (int i = threadIdx.x; i < 12 * kN; i += blockDim.x) h.edge_l[i] = (unsigned short)edge_node(i / kN, i % kN);
    if (threadIdx.x < 8) h.corner_l[threadIdx.x] = (unsigned short)corner_node(threadIdx.x);
    for (int l = threadIdx.x; l < kN3; l += blockDim.x)
        h.bit[l] = (unsigned char)entity_bit(local_entity(l / kNN, (l / kN) % kN, l % kN));
    __syncthreads();
}

// One pad slot per row of 8 spreads the strided nodes of a face over the shared-memory banks.
constexpr int kStage = kN3 + kN3 / kN;
__device__ inline int stage_slot(int l) { return l + l / kN; }

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

// visit(t, s, dirichlet, counted) for every copy t of the elements in `set`, in
// coalesced order: s is the sum over all copies of its node, counted marks the counting
// copy. Called by every thread of the block (it synchronises the block once, to build
// the tables).
template <typename Value, typename Visit>
__device__ inline void for_each_copy_sum(const TopologyView& T, const ElementList& set, const Value& value,
                                         Visit&& visit)
{
    constexpr unsigned kAll = 0xffffffffu;
    constexpr int kFaceSlots = 6 * kFaceNodes, kEdgeSlots = 12 * kEdgeNodes, kStarSlots = kEdgeSlots + 8;
    __shared__ HexTables h;
    __shared__ double stage_all[kThreads / 32][kStage];
    build_hex_tables(h);
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
        const long long base = e * kN3;
#pragma unroll
        for (int it = 0; it < kN3 / 32; ++it) stage[stage_slot(it * 32 + lane)] = value(base + it * 32 + lane);
        __syncwarp();
        // Shuffles stay outside conditions: every lane must execute each of them.
#pragma unroll
        for (int p0 = 0; p0 < kFaceSlots; p0 += 32) {
            const int p = p0 + lane;
            const int f = p < kFaceSlots ? p / kFaceNodes : 0;
            const int nbr = __shfl_sync(kAll, d0, f);
            const int code = __shfl_sync(kAll, d1, f);
            if (p < kFaceSlots && nbr >= 0) {
                const int q = h.face_map[(code & 63) * kFaceNodes + p - f * kFaceNodes];
                stage[stage_slot(h.face_l[p])] += value((long long)(nbr / 6) * kN3 + h.face_l[(nbr % 6) * kFaceNodes + q]);
            }
        }
        for (int p0 = 0; p0 < kStarSlots; p0 += 32) {
            const int p = p0 + lane;
            const int k = p < kEdgeSlots ? p / kEdgeNodes : 0;
            const int src = p < kEdgeSlots ? kEdgeBit + k : (p < kStarSlots ? kVertexBit + p - kEdgeSlots : 0);
            const int n0 = __shfl_sync(kAll, d0, src);
            const int n1 = __shfl_sync(kAll, d1, src);
            const int n2 = __shfl_sync(kAll, d2, src);
            if (p < kEdgeSlots) {
                const int t = 1 + p - k * kEdgeNodes;
                const int tc = n0 < 0 ? kP - t : t;   // canonical position on the edge
                stage[stage_slot(h.edge_l[k * kN + t])] = star_sum(value, T.edge_ent, n1, n2, [&](int c) {
                    const int m = c & 0x7fffffff;
                    return (long long)(m / 12) * kN3 + h.edge_l[(m % 12) * kN + (c < 0 ? kP - tc : tc)];
                });
            } else if (p < kStarSlots) {
                stage[stage_slot(h.corner_l[p - kEdgeSlots])] = star_sum(
                    value, T.vert_ent, n1, n2, [&](int c) { return (long long)(c / 8) * kN3 + h.corner_l[c % 8]; });
            }
        }
        __syncwarp();
#pragma unroll
        for (int it = 0; it < kN3 / 32; ++it) {
            const int l = it * 32 + lane;
            const int bit = h.bit[l];
            visit(base + l, stage[stage_slot(l)], ((dirichlet >> bit) & 1u) != 0, ((counting >> bit) & 1u) != 0);
        }
        __syncwarp();   // the next element reuses the stage
    }
}

__global__ void __launch_bounds__(kThreads)
unstructured_dss_kernel(const double* __restrict__ in, const double* __restrict__ ghost, double* __restrict__ out,
                        TopologyView T, ElementList set)
{
    const GhostedValues value{in, ghost, T.local * kN3};
    for_each_copy_sum(T, set, value, [&](long long t, double s, bool, bool) { out[t] = s; });
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
__global__ void __launch_bounds__(kThreads)
unstructured_precondition_kernel(const double* __restrict__ q, const double* __restrict__ ghost,
                                 double* __restrict__ z, const double* __restrict__ diag,
                                 const double* __restrict__ a, TopologyView T, ElementList set,
                                 double* __restrict__ partial)
{
    constexpr int NV = AZ + ZZ;
    [[maybe_unused]] double acc[NV > 0 ? NV : 1] = {};
    const GhostedValues value{q, ghost, T.local * kN3};
    for_each_copy_sum(T, set, value, [&](long long t, double s, bool dirichlet, bool counted) {
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

// u = 0 on every local copy of a Dirichlet node: a homogeneous Dirichlet condition on an
// element-local vector (e.g. a manufactured solution on quantized coordinates).
__global__ void zero_dirichlet_kernel(double* __restrict__ u, TopologyView T)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= T.local * kN3) return;
    const long long e = t / kN3;
    const int l = (int)(t % kN3);
    const LocalEntity le = local_entity(l / kNN, (l / kN) % kN, l % kN);
    bool dirichlet = false;
    if (le.kind == 1) dirichlet = (T.face_code[e * 6 + le.index] >> 6) != 0;
    else if (le.kind == 2) dirichlet = T.edge_dirichlet[T.edge_of[e * 12 + le.index] & 0x7fffffff] != 0;
    else if (le.kind == 3) dirichlet = T.vert_dirichlet[T.vert_of[e * 8 + le.index]] != 0;
    if (dirichlet) u[t] = 0.0;
}

inline void zero_dirichlet(double* d_u, const UnstructuredTopology& T, cudaStream_t stream = 0)
{
    const long long n = T.local * kN3;
    if (n > 0) zero_dirichlet_kernel<<<topo_detail::blocks(n), kThreads, 0, stream>>>(d_u, view(T));
    MARS_CELLWISE_CK(cudaGetLastError());
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

// The interior nodes of a face (36), an edge (6) or a corner (1).
__device__ inline int entity_nodes(int kind) { return kind == 1 ? kFaceNodes : (kind == 2 ? kEdgeNodes : 1); }
__device__ inline int entity_node(int kind, int index, int j)
{
    if (kind == 1) return face_node(index, 1 + j / 6, 1 + j % 6);
    if (kind == 2) return edge_node(index, 1 + j);
    return corner_node(index);
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
                               const unsigned long long* __restrict__ gid, const int* __restrict__ owner,
                               long long* __restrict__ send_at, long long* __restrict__ recv_at, Entries send,
                               Entries recv)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    const long long faces = T.local * 6;
    if (t >= faces + edges + vertices) return;
    long long so = FILL ? send_at[t] : 0, ro = FILL ? recv_at[t] : 0;
    auto emit = [&](Entries& out, long long& o, int peer, long long elem, int kind, int index, long long base) {
        const int n = entity_nodes(kind);
        if (FILL)
            for (int j = 0; j < n; ++j) {
                const int l = entity_node(kind, index, j);
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
            emit(send, so, owner[c / 6], e, 1, (int)(t % 6), e * kN3);
            emit(recv, ro, owner[c / 6], c / 6, 1, c % 6, (long long)(c / 6 - T.local) * kN3);
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
                    emit(recv, ro, owner[el], el, kind, c % per, (el - T.local) * kN3);
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
                    if (first) emit(send, so, owner[ej], el, kind, c % per, el * kN3);
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

__global__ void pack_kernel(const double* __restrict__ q, const long long* __restrict__ index, long long n,
                            double* __restrict__ out)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < n) out[t] = q[index[t]];
}

__global__ void unpack_kernel(const double* __restrict__ in, const long long* __restrict__ index, long long n,
                              double* __restrict__ ghost)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < n) ghost[index[t]] = in[t];
}

}  // namespace halo_detail

// The ghost exchange of one rank's tables. gid: a global id of each of the T.elements
// elements (any partition-independent, unique number); owner: the rank that owns it.
// start() sends the local copies the peers need, the caller runs the elements that read
// no ghost copies, finish() waits and files the received copies into the ghost blocks.
// Device buffers throughout: GPU-aware MPI (on Cray MPICH, MPICH_GPU_SUPPORT_ENABLED=1).
class UnstructuredHalo {
public:
    UnstructuredHalo(const UnstructuredTopology& T, const unsigned long long* d_gid, const int* d_owner,
                     MPI_Comm comm)
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
            entries_kernel<false><<<blocks(items), kThreads>>>(v, T.edges(), T.vertices(), d_gid, d_owner,
                                                               raw(send_at), raw(recv_at), none, none);
        MARS_CELLWISE_CK(cudaGetLastError());
        const long long send_total = thrust::reduce(thrust::device, send_at.begin(), send_at.end(), 0LL);
        const long long recv_total = thrust::reduce(thrust::device, recv_at.begin(), recv_at.end(), 0LL);
        thrust::exclusive_scan(thrust::device, send_at.begin(), send_at.end(), send_at.begin());
        thrust::exclusive_scan(thrust::device, recv_at.begin(), recv_at.end(), recv_at.begin());
        thrust::device_vector<int> speer(send_total), rpeer(recv_total), snode(send_total), rnode(recv_total);
        thrust::device_vector<unsigned long long> sgid(send_total), rgid(recv_total);
        send_index_.resize(send_total);
        recv_index_.resize(recv_total);
        if (items > 0)
            entries_kernel<true><<<blocks(items), kThreads>>>(
                v, T.edges(), T.vertices(), d_gid, d_owner, raw(send_at), raw(recv_at),
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
        long long sofs = 0, rofs = 0;
        for (int p = 0; p < size; ++p) {
            if (scount[p] == 0 && rcount[p] == 0) continue;
            peers_.push_back({p, sofs, scount[p], rofs, rcount[p]});
            sofs += scount[p];
            rofs += rcount[p];
        }
        MARS_CELLWISE_MPI(MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_INT, MPI_MAX, comm_));
        if (bad) {
            if (rank == 0) fprintf(stderr, "cell-wise halo: send and receive counts disagree between ranks\n");
            MPI_Abort(comm_, 1);
        }
        requests_.reserve(2 * peers_.size());

        if (send_total > 0) MARS_CELLWISE_CK(cudaMalloc(&send_, send_total * sizeof(double)));
        if (recv_total > 0) MARS_CELLWISE_CK(cudaMalloc(&recv_, recv_total * sizeof(double)));
        const long long ghost_values = (T.elements - L) * kN3;
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
    int peers() const { return (int)peers_.size(); }

    void start(const double* q, cudaStream_t stream)
    {
        using namespace halo_detail;
        if (send_total_ > 0)
            pack_kernel<<<blocks(send_total_), kThreads, 0, stream>>>(q, raw(send_index_), send_total_, send_);
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
    long long send_total_ = 0, recv_total_ = 0;
    double *send_ = nullptr, *recv_ = nullptr, *ghost_ = nullptr;
    thrust::device_vector<int> interior_ids_, boundary_ids_;
    ElementList interior_, boundary_;
    std::vector<MPI_Request> requests_;
};

// Runs a pass over the local elements: one launch on one rank; otherwise start the
// exchange, launch the elements that read no ghost copies, wait, launch the rest.
// launch(set, ghost, first) returns its number of thread blocks; `first` is how many
// blocks earlier launches used (for partial sums). Returns the total.
template <typename Launch>
int run_unstructured_pass(UnstructuredHalo* halo, long long local, const double* q, cudaStream_t stream,
                          Launch&& launch)
{
    if (!halo || !halo->active()) return launch(ElementList{nullptr, local}, nullptr, 0);
    halo->start(q, stream);
    const int first = launch(halo->interior(), nullptr, 0);
    halo->finish(stream);
    return first + launch(halo->boundary(), halo->ghost(), first);
}

// out = DSS(in) on every local copy. `halo` is null on one rank.
inline void dss(const double* d_in, double* d_out, const UnstructuredTopology& T, UnstructuredHalo* halo = nullptr,
                cudaStream_t stream = 0)
{
    const TopologyView v = view(T);
    run_unstructured_pass(halo, T.local, d_in, stream, [&](const ElementList& set, const double* ghost, int) {
        static const int grid = resident_grid(unstructured_dss_kernel);
        unstructured_dss_kernel<<<grid, kThreads, 0, stream>>>(d_in, ghost, d_out, v, set);
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
        return run_unstructured_pass(halo, T.local, q, stream, [&](const ElementList& set, const double* ghost,
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
