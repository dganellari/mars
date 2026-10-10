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
#include <thrust/execution_policy.h>
#include <thrust/for_each.h>
#include <thrust/gather.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/reduce.h>
#include <thrust/scan.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <thrust/transform.h>

#include <cstdint>
#include <cstdio>
#include <cstdlib>

namespace mars {
namespace cellwise {

// ---- tables ----------------------------------------------------------------------------

// Entity bits of an element's 32-bit masks (counting copies, Dirichlet entities): face f
// at bit f, edge k at kEdgeBit + k, corner c at kVertexBit + c, interior nodes at
// kInteriorBit. The gather kernel loads each entity's descriptor on the lane with the
// same number, so one ballot gives the Dirichlet mask.
constexpr int kEdgeBit = 8, kVertexBit = 20, kInteriorBit = 31;

// One rank's tables over its local elements. Copies are coded as e * 6 + f (faces),
// e * 12 + k (edges) or e * 8 + c (corners). Ranges are in canonical order throughout.
struct UnstructuredTopology {
    long long elements = 0;
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

// Builds the tables of E local elements. key[c][e]: global key of element e's corner c
// (partition-independent, e.g. its SFC key); lid[c][e]: a dense local id of the same
// corner (equal ids iff equal keys). Both arrays of 8 device pointers live on the device.
inline UnstructuredTopology build_topology(const unsigned long long* const* d_key, const int* const* d_lid,
                                           long long E)
{
    using namespace topo_detail;
    if (E * 12 > 0x7fffffffLL) {   // copy codes e * 12 + k are 31-bit
        fprintf(stderr, "cell-wise topology: %lld elements on one rank, at most %lld\n", E, 0x7fffffffLL / 12);
        std::abort();
    }
    UnstructuredTopology T;
    T.elements = E;
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
    long long elements;
    const int *face_nbr, *edge_of, *edge_off, *edge_ent, *vert_of, *vert_off, *vert_ent;
    const unsigned char *face_code, *edge_dirichlet, *vert_dirichlet;
    const unsigned* counting;
};

inline TopologyView view(const UnstructuredTopology& T)
{
    using topo_detail::raw;
    return {T.elements,         raw(T.face_nbr),  raw(T.edge_of),         raw(T.edge_off),
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

// visit(t, s, dirichlet, counted) for every local copy t, in coalesced order: s is the
// sum over all copies of its node, counted marks the counting copy. Called by every
// thread of the block (it synchronises the block once, to build the tables).
template <typename Value, typename Visit>
__device__ inline void for_each_copy_sum(const TopologyView& T, const Value& value, Visit&& visit)
{
    constexpr unsigned kAll = 0xffffffffu;
    constexpr int kFaceSlots = 6 * kFaceNodes, kEdgeSlots = 12 * kEdgeNodes, kStarSlots = kEdgeSlots + 8;
    __shared__ HexTables h;
    __shared__ double stage_all[kThreads / 32][kStage];
    build_hex_tables(h);
    const int lane = threadIdx.x & 31;
    double* stage = stage_all[threadIdx.x >> 5];
    const long long warps = ((long long)gridDim.x * blockDim.x) >> 5;
    for (long long e = ((long long)blockIdx.x * blockDim.x + threadIdx.x) >> 5; e < T.elements; e += warps) {
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
unstructured_dss_kernel(const double* __restrict__ in, double* __restrict__ out, TopologyView T)
{
    const auto value = [in](long long i) { return in[i]; };
    for_each_copy_sum(T, value, [&](long long t, double s, bool, bool) { out[t] = s; });
}

inline void dss(const double* d_in, double* d_out, const UnstructuredTopology& T, cudaStream_t stream = 0)
{
    static const int grid = resident_grid(unstructured_dss_kernel);
    unstructured_dss_kernel<<<grid, kThreads, 0, stream>>>(d_in, d_out, view(T));
    MARS_CELLWISE_CK(cudaGetLastError());
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
unstructured_precondition_kernel(const double* __restrict__ q, double* __restrict__ z,
                                 const double* __restrict__ diag, const double* __restrict__ a,
                                 TopologyView T, double* __restrict__ partial)
{
    constexpr int NV = AZ + ZZ;
    [[maybe_unused]] double acc[NV > 0 ? NV : 1] = {};
    const auto value = [q](long long i) { return q[i]; };
    for_each_copy_sum(T, value, [&](long long t, double s, bool dirichlet, bool counted) {
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

// The Jacobi preconditioner of the Krylov solvers (same interface as
// JacobiPreconditioner in mars_cellwise_krylov.hpp), on one rank.
struct UnstructuredJacobi {
    TopologyView T;
    const double* diag;

    template <bool AZ, bool ZZ>
    int operator()(const double* q, double* z, const double* a, double* partial, cudaStream_t stream) const
    {
        static const int grid = resident_grid(unstructured_precondition_kernel<AZ, ZZ>);
        unstructured_precondition_kernel<AZ, ZZ><<<grid, kThreads, 0, stream>>>(q, z, diag, a, T, partial);
        MARS_CELLWISE_CK(cudaGetLastError());
        return grid;
    }
};

}  // namespace cellwise
}  // namespace mars
