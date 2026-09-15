#ifndef MARS_HO_DOF_HANDLER_TET_GPU_HPP
#define MARS_HO_DOF_HANDLER_TET_GPU_HPP

// GPU build of the tet HO DOF numbering -- device twin of HoCvfemTetDofHandler::build
// and HoTetDofHandler::build in mars_ho_dof_handler_tet.hpp.
//
// Both host builders key a node by the sorted list of (global vertex id, weight)
// pairs with zero weights dropped, and dedup through a std::map. That map is the
// only host-bound part: it caps single-GPU problem size and setup time, exactly as
// HODofHandler::build did for hex before buildGpu.
//
// Device algorithm (validated against the std::map reference before being written):
//   1. pack the key into FOUR uint64 lanes, one per barycentric slot:
//        lane = (uint32)gid << 32 | wq        with unused slots (0x7fffffff, 0)
//      so a short key pads canonically and still compares correctly.
//   2. sort by STABLE LSD RADIX over the lanes, least significant lane first.
//      Four thrust::stable_sort_by_key calls on uint64 keys stay on thrust's radix
//      path. A comparator over the whole 32-byte key would drop thrust into merge
//      sort, which is what crashed the hex path at scale.
//   3. flag unique boundaries, scan to dense ids, scatter back into elemDof.
//
// As on the hex path the device ids come out in sorted-key order and the host ids in
// first-encounter order, so the two numberings are a PERMUTATION of each other. What
// is identical is the partition: which (element, node) slots share a DOF.
//
// MEMORY: peak while building is ~48 B/node over nElem*NN nodes -- four uint64 lanes
// (32) plus perm, isNew, ids and elemDof (4 each) -- all freed before return except
// elemDof. That is ~12x elemDof itself, and it is the scale limit of this approach.

#include "mars_ho_dof_handler_tet.hpp"

#include <algorithm>
#include <cmath>

#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/gather.h>
#include <thrust/scatter.h>
#include <thrust/scan.h>
#include <thrust/transform.h>
#include <thrust/for_each.h>
#include <thrust/copy.h>
#include <thrust/sequence.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>

#include <array>
#include <cstdint>
#include <vector>

namespace mars {
namespace fem {

namespace ho_tet_gpu_detail {

// Unused slots sort last, so a 1-, 2- or 3-pair key pads to a canonical 4-slot record.
static constexpr uint32_t kUnusedGid = 0x7fffffffu;

// A weight of exactly 1 quantizes to 2^32, one bit too wide. It occurs only at an
// exact corner, where it is the sole pair; the nearest other weight a GLL grid can
// produce is O(1/n^2) away, far outside one ulp of the clamp, so clamping cannot
// merge distinct nodes.
__host__ __device__ inline uint64_t packSlot(int gid, long long wq) {
    uint64_t w = (wq < 0) ? 0ull : (uint64_t)wq;
    if (w > 0xffffffffull) w = 0xffffffffull;
    return ((uint64_t)(uint32_t)gid << 32) | w;
}

__host__ __device__ inline uint64_t unusedSlot() {
    return ((uint64_t)kUnusedGid << 32);
}

// Insertion sort of at most 4 slots -- the pairs are already nearly ordered and the
// count is fixed, so this is cheaper than anything general.
__host__ __device__ inline void sortSlots(uint64_t* s) {
    for (int i = 1; i < 4; ++i) {
        uint64_t v = s[i];
        int j = i - 1;
        while (j >= 0 && s[j] > v) { s[j + 1] = s[j]; --j; }
        s[j + 1] = v;
    }
}

// The collapsed (Duffy) barycentric weights of tensor node (ia, ib, ic), matching
// HoCvfemTetDofHandler::build exactly.
__host__ __device__ inline void duffyWeights(double a, double b, double c, double* w) {
    const double r = a * (1.0 - b) * (1.0 - c);
    const double s = b * (1.0 - c);
    const double t = c;
    w[0] = 1.0 - r - s - t; w[1] = r; w[2] = s; w[3] = t;
}

// Shared tail: sort the packed lanes, assign dense ids, scatter into elemDof.
// Returns numDof and leaves the sorted permutation in `perm` so the caller can pick
// one representative node per DOF (for dofPos).
inline long dedupAndScatter(thrust::device_vector<uint64_t> lane[4],
                            long N,
                            thrust::device_vector<int>& perm,
                            thrust::device_vector<int>& d_elemDof,
                            thrust::device_vector<int>& firstOf)
{
    perm.resize(N);
    thrust::sequence(perm.begin(), perm.end());

    thrust::device_vector<uint64_t> tmp(N);
    for (int l = 3; l >= 0; --l) {           // least significant lane first
        thrust::gather(perm.begin(), perm.end(), lane[l].begin(), tmp.begin());
        thrust::stable_sort_by_key(tmp.begin(), tmp.end(), perm.begin());
    }
    tmp.clear(); tmp.shrink_to_fit();     // 8 B/node back before the scan buffers

    // A node starts a new DOF when any lane differs from its predecessor.
    thrust::device_vector<int> isNew(N);
    const uint64_t* l0 = thrust::raw_pointer_cast(lane[0].data());
    const uint64_t* l1 = thrust::raw_pointer_cast(lane[1].data());
    const uint64_t* l2 = thrust::raw_pointer_cast(lane[2].data());
    const uint64_t* l3 = thrust::raw_pointer_cast(lane[3].data());
    const int* p = thrust::raw_pointer_cast(perm.data());
    thrust::transform(thrust::device, thrust::counting_iterator<long>(0),
                      thrust::counting_iterator<long>(N), isNew.begin(),
                      [=] __host__ __device__ (long i) -> int {
                          if (i == 0) return 1;
                          const int a = p[i], b = p[i - 1];
                          return (l0[a] != l0[b] || l1[a] != l1[b] ||
                                  l2[a] != l2[b] || l3[a] != l3[b]) ? 1 : 0;
                      });

    thrust::device_vector<int> ids(N);
    thrust::inclusive_scan(isNew.begin(), isNew.end(), ids.begin());
    thrust::transform(thrust::device, ids.begin(), ids.end(), ids.begin(),
                      [] __host__ __device__ (int v) { return v - 1; });

    const long numDof = N ? (long)ids.back() + 1 : 0;

    // elemDof[original slot] = dense id
    d_elemDof.resize(N);
    thrust::scatter(ids.begin(), ids.end(), perm.begin(), d_elemDof.begin());

    // One representative original slot per DOF, for the physical position.
    firstOf.assign(numDof, -1);
    int* fo = thrust::raw_pointer_cast(firstOf.data());
    const int* idp = thrust::raw_pointer_cast(ids.data());
    const int* isn = thrust::raw_pointer_cast(isNew.data());
    thrust::for_each(thrust::device, thrust::counting_iterator<long>(0),
                     thrust::counting_iterator<long>(N),
                     [=] __host__ __device__ (long i) {
                         if (isn[i]) fo[idp[i]] = p[i];
                     });
    return numDof;
}

}  // namespace ho_tet_gpu_detail

// ---- CVFEM collapsed GLL tensor grid (the one mars_ho_tet_perf drives) ----
inline void buildGpu(HoCvfemTetDofHandler& dof,
                     const std::vector<int>& elemCorners,
                     const std::vector<std::array<double,3>>& coords,
                     const std::vector<double>& Z)
{
    using namespace ho_tet_gpu_detail;
    const int n = (int)Z.size();
    const int NN = n * n * n;
    const long nElem = (long)elemCorners.size() / 4;
    const long N = nElem * NN;
    dof.n = n; dof.P = n - 1; dof.NN = NN;

    // Ascending-gid corner order per element: the host does this so a shared face
    // induces the same canonical Duffy grid from both sides.
    std::vector<int> h_sorted(elemCorners.size());
    for (long e = 0; e < nElem; ++e) {
        int g[4] = { elemCorners[e*4+0], elemCorners[e*4+1],
                     elemCorners[e*4+2], elemCorners[e*4+3] };
        std::sort(g, g + 4);
        for (int c = 0; c < 4; ++c) h_sorted[e*4+c] = g[c];
    }
    dof.sortedCorners = h_sorted;

    thrust::device_vector<int>    d_corners(h_sorted.begin(), h_sorted.end());
    thrust::device_vector<double> d_Z(Z.begin(), Z.end());

    thrust::device_vector<uint64_t> lane[4];
    for (int l = 0; l < 4; ++l) lane[l].resize(N);

    const int*    dc = thrust::raw_pointer_cast(d_corners.data());
    const double* dz = thrust::raw_pointer_cast(d_Z.data());
    uint64_t* L0 = thrust::raw_pointer_cast(lane[0].data());
    uint64_t* L1 = thrust::raw_pointer_cast(lane[1].data());
    uint64_t* L2 = thrust::raw_pointer_cast(lane[2].data());
    uint64_t* L3 = thrust::raw_pointer_cast(lane[3].data());
    const double Q = 4294967296.0;

    thrust::for_each(thrust::device, thrust::counting_iterator<long>(0),
                     thrust::counting_iterator<long>(N),
                     [=] __host__ __device__ (long idx) {
                         const long e  = idx / NN;
                         const int  nd = (int)(idx - e * NN);
                         const int ia = nd / (n * n), ib = (nd / n) % n, ic = nd % n;
                         double w[4];
                         duffyWeights(dz[ia], dz[ib], dz[ic], w);
                         uint64_t s[4];
                         for (int c = 0; c < 4; ++c) {
                             s[c] = (w[c] > 1e-9)
                                  ? packSlot(dc[e * 4 + c], (long long)llround(w[c] * Q))
                                  : unusedSlot();
                         }
                         sortSlots(s);
                         L0[idx] = s[0]; L1[idx] = s[1]; L2[idx] = s[2]; L3[idx] = s[3];
                     });

    thrust::device_vector<int> perm, d_elemDof, firstOf;
    dof.numDof = dedupAndScatter(lane, N, perm, d_elemDof, firstOf);
    for (int l = 0; l < 4; ++l) lane[l].clear(), lane[l].shrink_to_fit();

    dof.elemDof.resize(N);
    thrust::copy(d_elemDof.begin(), d_elemDof.end(), dof.elemDof.begin());

    // Physical positions from one representative node per DOF.
    std::vector<int> h_first(dof.numDof);
    thrust::copy(firstOf.begin(), firstOf.end(), h_first.begin());
    dof.dofPos.assign(dof.numDof, {0.0, 0.0, 0.0});
    for (long d = 0; d < dof.numDof; ++d) {
        const long idx = h_first[d];
        const long e  = idx / NN;
        const int  nd = (int)(idx - e * NN);
        const int ia = nd / (n * n), ib = (nd / n) % n, ic = nd % n;
        double w[4];
        duffyWeights(Z[ia], Z[ib], Z[ic], w);
        std::array<double,3> x = {0, 0, 0};
        for (int c = 0; c < 4; ++c)
            for (int k = 0; k < 3; ++k) x[k] += w[c] * coords[h_sorted[e*4+c]][k];
        dof.dofPos[d] = x;
    }
}

// ---- Nodal HO tet (integer barycentric indices) ----
template <typename RealType>
inline void buildGpu(HoTetDofHandler& dof,
                     const std::vector<int>& elemCorners,
                     const std::vector<std::array<double,3>>& coords,
                     const HoTetNodal<RealType>& nd)
{
    using namespace ho_tet_gpu_detail;
    const int P = nd.P, Np = nd.Np;
    const long nElem = (long)elemCorners.size() / 4;
    const long N = nElem * Np;
    dof.P = P; dof.Np = Np;

    // bary is a host table; push the three indices per node to the device.
    std::vector<int> h_bary((size_t)Np * 3);
    for (int m = 0; m < Np; ++m)
        for (int d = 0; d < 3; ++d) h_bary[(size_t)m * 3 + d] = nd.bary[m][d];

    thrust::device_vector<int> d_corners(elemCorners.begin(), elemCorners.end());
    thrust::device_vector<int> d_bary(h_bary.begin(), h_bary.end());

    thrust::device_vector<uint64_t> lane[4];
    for (int l = 0; l < 4; ++l) lane[l].resize(N);

    const int* dc = thrust::raw_pointer_cast(d_corners.data());
    const int* db = thrust::raw_pointer_cast(d_bary.data());
    uint64_t* L0 = thrust::raw_pointer_cast(lane[0].data());
    uint64_t* L1 = thrust::raw_pointer_cast(lane[1].data());
    uint64_t* L2 = thrust::raw_pointer_cast(lane[2].data());
    uint64_t* L3 = thrust::raw_pointer_cast(lane[3].data());

    thrust::for_each(thrust::device, thrust::counting_iterator<long>(0),
                     thrust::counting_iterator<long>(N),
                     [=] __host__ __device__ (long idx) {
                         const long e = idx / Np;
                         const int  m = (int)(idx - e * Np);
                         const int i = db[(long)m*3+0], j = db[(long)m*3+1], k = db[(long)m*3+2];
                         const int w[4] = { P - i - j - k, i, j, k };
                         uint64_t s[4];
                         for (int c = 0; c < 4; ++c)
                             s[c] = (w[c] > 0) ? packSlot(dc[e*4+c], (long long)w[c])
                                               : unusedSlot();
                         sortSlots(s);
                         L0[idx] = s[0]; L1[idx] = s[1]; L2[idx] = s[2]; L3[idx] = s[3];
                     });

    thrust::device_vector<int> perm, d_elemDof, firstOf;
    dof.numDof = dedupAndScatter(lane, N, perm, d_elemDof, firstOf);
    for (int l = 0; l < 4; ++l) lane[l].clear(), lane[l].shrink_to_fit();

    dof.elemDof.resize(N);
    thrust::copy(d_elemDof.begin(), d_elemDof.end(), dof.elemDof.begin());

    std::vector<int> h_first(dof.numDof);
    thrust::copy(firstOf.begin(), firstOf.end(), h_first.begin());
    dof.dofPos.assign(dof.numDof, {0.0, 0.0, 0.0});
    for (long d = 0; d < dof.numDof; ++d) {
        const long idx = h_first[d];
        const long e = idx / Np;
        const int  m = (int)(idx - e * Np);
        const int i = h_bary[(size_t)m*3+0], j = h_bary[(size_t)m*3+1], k = h_bary[(size_t)m*3+2];
        const int w[4] = { P - i - j - k, i, j, k };
        std::array<double,3> x = {0, 0, 0};
        for (int c = 0; c < 4; ++c)
            for (int q = 0; q < 3; ++q)
                x[q] += (w[c] / (double)P) * coords[elemCorners[e*4+c]][q];
        dof.dofPos[d] = x;
    }
}

}  // namespace fem
}  // namespace mars

#endif
