#pragma once
//
// Periodic BC support for unstructured AMR meshes.
//
// Strategy (user-space, no octree changes):
//   1) After mesh load (and after each AMR rebuild), scan node coordinates.
//      Nodes on the three "max" faces (x=xmax, y=ymax, z=zmax) are "slaves".
//      For each slave, find its master partner on the opposite "min" face by
//      spatial proximity in the perpendicular plane. Store partner index in
//      d_periodicPartner[slave] = master.
//   2) During DOF numbering, slave nodes do NOT get their own DOF; their
//      "DOF" is the master's DOF. The resulting linear system size = number
//      of owned interior+master DOFs (3x fewer than in a non-periodic mesh
//      along boundary faces, ~8 fewer at corners).
//   3) After per-node accumulation kernels (mass, advection, divergence,
//      gradient, etc.) the periodic-pair-sum kernel runs:
//        master_val += slave_val;
//        slave_val   = master_val;
//      so both sides of every periodic face see the same accumulated value.
//   4) Pressure null-space removal: subtract global mean (Allreduce + scale)
//      instead of pinning a single DOF.
//
// Limitations:
//   - Assumes axis-aligned box with mesh nodes coincident across periodic
//     faces (matching mesh). For TGV on a structured periodic cube this holds
//     by construction.
//   - 1D corner nodes (8 corners of cube) collapse to a single master DOF.
//   - Multi-rank: slave/master may live on different ranks. The pair-sum
//     must run *after* a reverseExchangeNodeHaloAdd (so ghost contributions
//     are summed into owners), then a forward exchangeNodeHalo redistributes
//     the merged values.
//

#include "backend/distributed/unstructured/domain.hpp"
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/binary_search.h>
#include <thrust/fill.h>
#include <thrust/transform.h>
#include <cmath>
#include <vector>

namespace mars
{
namespace fem
{

// Encodes "is this node periodic" + "which axis pair(s)" + "master partner".
// Bit layout in periodicMask:
//   bit 0: slave on x=xmax face         (partner has x=xmin)
//   bit 1: slave on y=ymax face         (partner has y=ymin)
//   bit 2: slave on z=zmax face         (partner has z=zmin)
//   bit 3: this is itself a master (some slave points to it)
// A node may have multiple bits set (e.g. edge of cube: bits 0+1).
// Masters and interior nodes have bit 3 (master) or none.
//
// d_periodicPartner[slave_node] = master_node_index (within local node array,
// including ghosts). For non-slaves: -1.

template<typename KeyType, typename RealType>
struct PeriodicMap
{
    using DevVecInt = cstone::DeviceVector<int>;
    using DevVecU8  = cstone::DeviceVector<uint8_t>;

    DevVecInt d_periodicPartner;   // [numNodes] -> master idx or -1
    DevVecU8  d_periodicMask;      // [numNodes] -> bitmask as above
    int       numSlaves     = 0;
    int       numMasters    = 0;

    // Bounds of the periodic box (recorded at build time so we can recognize
    // boundary faces after coordinate updates).
    RealType xmin = 0, xmax = 0, ymin = 0, ymax = 0, zmin = 0, zmax = 0;

    // Tolerance for "on the face" test. Set to 1e-4 * smallest_edge by default.
    RealType faceEps = RealType(1e-5);
    // Tolerance for matching partner in the perpendicular plane.
    RealType matchEps = RealType(1e-5);
};

// Kernel: classify each node and (for slaves) compute partner via spatial search.
// Two-pass approach because spatial search requires sorted master coordinates.
//
// Pass 1: build periodicMask, count slaves/masters, gather master candidates.
// Pass 2: for each slave, binary-search the master list (sorted by mapped key).
//
// The "mapped key" for partner matching collapses the slave coordinate to its
// would-be master position: if slave is on x=xmax, partner has x=xmin and the
// same (y,z). So we form a 64-bit key = quantize(y_master, z_master) for the
// x-axis pair, etc. For axis-aligned matching meshes this is exact.

template<typename RealType>
__global__ void classifyPeriodicNodesKernel(const RealType* d_x, const RealType* d_y,
                                            const RealType* d_z, size_t numNodes,
                                            RealType xmin, RealType xmax, RealType ymin,
                                            RealType ymax, RealType zmin, RealType zmax,
                                            RealType faceEps,
                                            uint8_t* d_mask)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;

    uint8_t m = 0;
    RealType x = d_x[i], y = d_y[i], z = d_z[i];
    if (fabs(x - xmax) < faceEps) m |= 0x01;
    if (fabs(y - ymax) < faceEps) m |= 0x02;
    if (fabs(z - zmax) < faceEps) m |= 0x04;
    d_mask[i] = m;
}

template<typename RealType>
__global__ void countMastersKernel(const RealType* d_x, const RealType* d_y,
                                   const RealType* d_z, size_t numNodes,
                                   RealType xmin, RealType ymin, RealType zmin,
                                   RealType faceEps, int* d_isMaster)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    RealType x = d_x[i], y = d_y[i], z = d_z[i];
    int onXmin = (fabs(x - xmin) < faceEps) ? 1 : 0;
    int onYmin = (fabs(y - ymin) < faceEps) ? 1 : 0;
    int onZmin = (fabs(z - zmin) < faceEps) ? 1 : 0;
    d_isMaster[i] = (onXmin | onYmin | onZmin);
}

// Match slaves to masters via three sorted-search passes (one per axis).
// For axis A=x: slave on x=xmax with coords (xmax, y_s, z_s) must find
// master with coords (xmin, y_s, z_s). We pack the perpendicular coords as
// a 64-bit quantized key and binary-search.
//
// Quantization: q = uint32_t((v - vmin) / (vmax - vmin) * 2^31) — gives
// ~10 decimal digits of precision, plenty for matching meshes.

template<typename RealType>
__device__ inline uint64_t quantizePair(RealType a, RealType amin, RealType amax,
                                        RealType b, RealType bmin, RealType bmax)
{
    RealType na = (a - amin) / (amax - amin);
    RealType nb = (b - bmin) / (bmax - bmin);
    uint32_t qa = uint32_t(na * RealType(2147483647.0));   // 2^31 - 1
    uint32_t qb = uint32_t(nb * RealType(2147483647.0));
    return (uint64_t(qa) << 32) | uint64_t(qb);
}

// Build (key, idx) pairs for nodes on a given "min" face.
template<typename RealType>
__global__ void packMasterKeysKernel(const RealType* d_x, const RealType* d_y,
                                     const RealType* d_z, size_t numNodes,
                                     RealType xmin, RealType xmax, RealType ymin,
                                     RealType ymax, RealType zmin, RealType zmax,
                                     RealType faceEps, int axis,
                                     uint64_t* d_keys, int* d_idx, int* d_counter)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    RealType x = d_x[i], y = d_y[i], z = d_z[i];
    bool onMin = false;
    uint64_t key = 0;
    if (axis == 0 && fabs(x - xmin) < faceEps) {
        onMin = true; key = quantizePair(y, ymin, ymax, z, zmin, zmax);
    } else if (axis == 1 && fabs(y - ymin) < faceEps) {
        onMin = true; key = quantizePair(x, xmin, xmax, z, zmin, zmax);
    } else if (axis == 2 && fabs(z - zmin) < faceEps) {
        onMin = true; key = quantizePair(x, xmin, xmax, y, ymin, ymax);
    }
    if (onMin) {
        int slot = atomicAdd(d_counter, 1);
        d_keys[slot] = key;
        d_idx[slot]  = int(i);
    }
}

// For each slave node on the "max" face for this axis, binary-search the
// master keys and record partner. Writes into d_partner only if entry == -1
// (preserves the first match — useful for corners where multiple axes apply
// but we link to the "x-master" which itself can chain to y/z masters).
template<typename RealType>
__global__ void matchSlavesKernel(const RealType* d_x, const RealType* d_y,
                                  const RealType* d_z, size_t numNodes,
                                  RealType xmin, RealType xmax, RealType ymin,
                                  RealType ymax, RealType zmin, RealType zmax,
                                  RealType faceEps, int axis,
                                  const uint64_t* d_masterKeys, const int* d_masterIdx,
                                  int numMasters, int* d_partner)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    RealType x = d_x[i], y = d_y[i], z = d_z[i];
    bool onMax = false;
    uint64_t key = 0;
    if (axis == 0 && fabs(x - xmax) < faceEps) {
        onMax = true; key = quantizePair(y, ymin, ymax, z, zmin, zmax);
    } else if (axis == 1 && fabs(y - ymax) < faceEps) {
        onMax = true; key = quantizePair(x, xmin, xmax, z, zmin, zmax);
    } else if (axis == 2 && fabs(z - zmax) < faceEps) {
        onMax = true; key = quantizePair(x, xmin, xmax, y, ymin, ymax);
    }
    if (!onMax) return;
    if (d_partner[i] != -1) return;  // already linked by an earlier axis

    // Binary search for exact key match (matching mesh assumption).
    int lo = 0, hi = numMasters;
    while (lo < hi) {
        int mid = (lo + hi) >> 1;
        if (d_masterKeys[mid] < key) lo = mid + 1; else hi = mid;
    }
    if (lo < numMasters && d_masterKeys[lo] == key) {
        d_partner[i] = d_masterIdx[lo];
    }
}

// Flatten partner chains: if partner[i]=j and partner[j]=k (chained at corner),
// rewrite partner[i]=k so a single lookup reaches the ultimate master.
// At most 3 levels of chaining (x->y->z masters at a corner).
__global__ void flattenPartnerChainKernel(int* d_partner, size_t numNodes)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    int p = d_partner[i];
    int hops = 0;
    while (p >= 0 && d_partner[p] != -1 && hops < 4) {
        p = d_partner[p];
        ++hops;
    }
    if (d_partner[i] != -1) d_partner[i] = p;
}

// After interior assembly: copy slave-row contributions into master rows.
// For node-resident scalar fields (mass, divergence, gradient component).
// Call sequence:
//   1) per-element scatter populates d_field on all nodes (slaves included)
//   2) reverseExchangeNodeHaloAdd sums ghost contributions to owners
//   3) periodicPairSumKernel: master += slave; slave = master
//   4) forward exchangeNodeHalo to redistribute
template<typename RealType>
__global__ void periodicPairSumKernel(const int* d_partner, size_t numNodes,
                                      RealType* d_field)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    int master = d_partner[i];
    if (master < 0) return;
    // Atomic add into master, then later broadcast back. Two-pass:
    // pass A: master += slave (atomic, slave zeroed)
    // pass B: slave = master (read-only)
    // For simplicity we use a single-pass scheme: in pass A only slaves write
    // (this kernel), zero out slaves; pass B reads masters.
    atomicAdd(&d_field[master], d_field[i]);
    d_field[i] = RealType(0);
}

template<typename RealType>
__global__ void periodicBroadcastKernel(const int* d_partner, size_t numNodes,
                                        RealType* d_field)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    int master = d_partner[i];
    if (master < 0) return;
    d_field[i] = d_field[master];
}

// Build the periodic map from current domain coords. Call after mesh load
// and after every AMR rebuild.
template<typename KeyType, typename RealType, typename DomainT>
void buildPeriodicMap(const DomainT& domain, PeriodicMap<KeyType, RealType>& map,
                      RealType xmin, RealType xmax, RealType ymin, RealType ymax,
                      RealType zmin, RealType zmax, RealType faceEps = 1e-5)
{
    size_t numNodes = domain.getNodeCount();
    map.xmin = xmin; map.xmax = xmax; map.ymin = ymin; map.ymax = ymax;
    map.zmin = zmin; map.zmax = zmax;
    map.faceEps = faceEps;
    map.d_periodicPartner.resize(numNodes);
    map.d_periodicMask.resize(numNodes);
    thrust::fill(thrust::device, map.d_periodicPartner.begin(),
                 map.d_periodicPartner.end(), -1);
    thrust::fill(thrust::device, map.d_periodicMask.begin(),
                 map.d_periodicMask.end(), uint8_t(0));

    const RealType* d_x = domain.getNodeX().data();
    const RealType* d_y = domain.getNodeY().data();
    const RealType* d_z = domain.getNodeZ().data();

    int block = 256, grid = int((numNodes + block - 1) / block);

    classifyPeriodicNodesKernel<<<grid, block>>>(d_x, d_y, d_z, numNodes,
                                                 xmin, xmax, ymin, ymax,
                                                 zmin, zmax, faceEps,
                                                 map.d_periodicMask.data());

    // For each axis: gather master nodes on "min" face, sort by quantized key,
    // then match every "max" face node to a master via binary search.
    for (int axis = 0; axis < 3; ++axis) {
        thrust::device_vector<uint64_t> d_keys(numNodes);
        thrust::device_vector<int>      d_idx(numNodes);
        thrust::device_vector<int>      d_counter(1, 0);
        packMasterKeysKernel<<<grid, block>>>(d_x, d_y, d_z, numNodes,
                                              xmin, xmax, ymin, ymax, zmin, zmax,
                                              faceEps, axis,
                                              thrust::raw_pointer_cast(d_keys.data()),
                                              thrust::raw_pointer_cast(d_idx.data()),
                                              thrust::raw_pointer_cast(d_counter.data()));
        int nm = 0;
        cudaMemcpy(&nm, thrust::raw_pointer_cast(d_counter.data()), sizeof(int),
                   cudaMemcpyDeviceToHost);
        if (nm == 0) continue;

        thrust::sort_by_key(d_keys.begin(), d_keys.begin() + nm, d_idx.begin());

        matchSlavesKernel<<<grid, block>>>(d_x, d_y, d_z, numNodes,
                                           xmin, xmax, ymin, ymax, zmin, zmax,
                                           faceEps, axis,
                                           thrust::raw_pointer_cast(d_keys.data()),
                                           thrust::raw_pointer_cast(d_idx.data()),
                                           nm,
                                           map.d_periodicPartner.data());
    }

    flattenPartnerChainKernel<<<grid, block>>>(map.d_periodicPartner.data(), numNodes);
    cudaDeviceSynchronize();

    // Count slaves for reporting.
    auto count_slaves = thrust::count_if(thrust::device,
                                         map.d_periodicPartner.begin(),
                                         map.d_periodicPartner.end(),
                                         [] __device__ (int v) { return v >= 0; });
    map.numSlaves = int(count_slaves);
}

// Top-level: call this after every per-node accumulation kernel that produced
// values needing periodic consistency. Handles the halo dance internally.
template<typename KeyType, typename RealType, typename DomainT>
void enforcePeriodicSum(const DomainT& domain, const PeriodicMap<KeyType, RealType>& map,
                        cstone::DeviceVector<RealType>& d_field,
                        const int* dofMapPtr = nullptr)
{
    size_t numNodes = d_field.size();
    int block = 256, grid = int((numNodes + block - 1) / block);

    // (a) sum ghost contributions into owners (cross-rank)
    domain.reverseExchangeNodeHaloAdd(d_field, dofMapPtr);

    // (b) within each rank: master += slave, slave = 0
    periodicPairSumKernel<<<grid, block>>>(map.d_periodicPartner.data(),
                                           numNodes, d_field.data());

    // (c) redistribute master values to ghosts
    domain.exchangeNodeHalo(d_field, dofMapPtr);

    // (d) broadcast master -> slave so both periodic copies hold the merged value
    periodicBroadcastKernel<<<grid, block>>>(map.d_periodicPartner.data(),
                                             numNodes, d_field.data());

    // (e) final halo sync so slave ghosts on other ranks see broadcast values
    domain.exchangeNodeHalo(d_field, dofMapPtr);
}

// Pressure null-space removal: subtract global mean.
template<typename RealType, typename DomainT>
void removeMean(const DomainT& domain, cstone::DeviceVector<RealType>& d_p,
                MPI_Comm comm)
{
    const auto& d_ownership = domain.getNodeOwnershipMap();
    size_t numNodes = d_p.size();

    // sum p over owned nodes
    RealType local_sum = thrust::transform_reduce(
        thrust::device,
        thrust::make_counting_iterator(size_t(0)),
        thrust::make_counting_iterator(numNodes),
        [d_p_ptr = d_p.data(), own_ptr = d_ownership.data()]
        __device__ (size_t i) -> RealType {
            return (own_ptr[i] == 1) ? d_p_ptr[i] : RealType(0);
        },
        RealType(0), thrust::plus<RealType>());

    long long local_n = thrust::transform_reduce(
        thrust::device,
        thrust::make_counting_iterator(size_t(0)),
        thrust::make_counting_iterator(numNodes),
        [own_ptr = d_ownership.data()]
        __device__ (size_t i) -> long long {
            return (own_ptr[i] == 1) ? 1LL : 0LL;
        },
        0LL, thrust::plus<long long>());

    RealType global_sum = 0;
    long long global_n = 0;
    MPI_Datatype mpi_real = std::is_same<RealType, double>::value ? MPI_DOUBLE : MPI_FLOAT;
    MPI_Allreduce(&local_sum, &global_sum, 1, mpi_real, MPI_SUM, comm);
    MPI_Allreduce(&local_n,   &global_n,   1, MPI_LONG_LONG, MPI_SUM, comm);

    if (global_n == 0) return;
    RealType mean = global_sum / RealType(global_n);

    thrust::transform(thrust::device, d_p.begin(), d_p.end(), d_p.begin(),
                      [mean] __device__ (RealType v) { return v - mean; });
}

}  // namespace fem
}  // namespace mars
