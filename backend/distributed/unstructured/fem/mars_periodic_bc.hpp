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
#include <thrust/for_each.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <mpi.h>
#include <unordered_map>
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

// Cross-rank periodic-pair exchange tables. For each owned periodic slave on
// this rank whose master is owned on a remote rank, we need an MPI exchange
// that ships the slave's contribution to the master's owner and then
// broadcasts the merged master value back. reverseExchangeNodeHaloAdd cannot
// bridge this because slave and master have DISTINCT SFC keys by design
// (no coord shift); they are NOT the same cstone-shared node.
template<typename KeyType, typename RealType>
struct CrossRankPeriodicMap
{
    template<typename T>
    using DevVec = cstone::DeviceVector<T>;

    // peer ranks we exchange periodic pairs with.
    std::vector<int> peers_;

    // Per-peer CSR offsets (size = peers_.size() + 1).
    std::vector<int> sendOffsets_;
    std::vector<int> recvOffsets_;

    // d_sendOwnedSlaveIds_[k]:  local OWNED slave node id (this rank) whose
    //                           master is owned by peer for bucket-of(k).
    // d_recvOwnedMasterIds_[k]: local OWNED master node id (this rank) which
    //                           is the target of a slave owned by peer for
    //                           bucket-of(k).
    DevVec<int> d_sendOwnedSlaveIds_;
    DevVec<int> d_recvOwnedMasterIds_;

    // Persistent staging buffers (mutable so const exchange methods can stage).
    mutable DevVec<RealType> sendBuf_;
    mutable DevVec<RealType> recvBuf_;

    // Epoch counter + private sub-communicator. comm_ is MPI_Comm_dup'd from
    // the user-supplied comm at first build; freed in PeriodicMap dtor.
    // Eliminates tag-collision with cstone halo (0x4d52/0x4d53 + epoch_).
    mutable int epoch_ = 0;
    MPI_Comm    comm_  = MPI_COMM_NULL;
};

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

    // Cross-rank pair tables, built at the end of buildPeriodicMap. Empty on
    // single-rank or when no cross-rank pairs exist.
    CrossRankPeriodicMap<KeyType, RealType> cross_;

    PeriodicMap() = default;
    ~PeriodicMap()
    {
        if (cross_.comm_ != MPI_COMM_NULL) { MPI_Comm_free(&cross_.comm_); }
    }
    PeriodicMap(const PeriodicMap&)            = delete;
    PeriodicMap& operator=(const PeriodicMap&) = delete;
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
// Gated intra-rank pair sum. Only accumulates when the master is locally
// OWNED. Cross-rank pairs (master is a ghost on this rank) are handled by
// crossRankPeriodicPairSum -- atomicAdd into a ghost slot here would be
// silently overwritten by the next forward exchangeNodeHalo. On single-rank
// d_ownership[master] is always 1, so the gate is a no-op (behavior matches
// the old version).
template<typename RealType>
__global__ void periodicPairSumKernel(const int*     d_partner,
                                      const uint8_t* d_ownership,
                                      size_t         numNodes,
                                      RealType*      d_field)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    int master = d_partner[i];
    if (master < 0) return;
    if (d_ownership[master] != 1) return;
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

// Per-DOF variant of periodicBroadcastKernel. Used inside the velocity CG's
// halo callback: the CG search vector p is sized numTotalDofs and indexed by
// DOF, but the periodic partner table is indexed by NODE. We remap via the
// node->dof table: copy p[dof[master_node]] -> p[dof[slave_node]] so the
// slave's DOF value tracks its master's DOF value every iteration. Without
// this, the multi-rank velocity CG references unsynced periodic-ghost DOF
// columns inside A*p and breaks down (pAp<=0 -> NaN).
template<typename RealType>
__global__ void periodicBroadcastDofKernel(const int* d_partner,
                                           const int* d_nodeToDof,
                                           size_t numNodes,
                                           RealType* d_dofField)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    int master = d_partner[i];
    if (master < 0) return;
    int sdof = d_nodeToDof[i];
    int mdof = d_nodeToDof[master];
    if (sdof < 0 || mdof < 0) return;
    if (sdof == mdof) return;   // already collapsed (same-rank pair)
    d_dofField[sdof] = d_dofField[mdof];
}

// =============================================================================
// Cross-rank periodic-pair exchange kernels (pack/unpack helpers).
// =============================================================================

template<typename RealType>
__global__ void packCrossRankSendKernel(const int* d_send_ids,
                                        const RealType* d_field,
                                        RealType* d_send_buf,
                                        int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    d_send_buf[i] = d_field[d_send_ids[i]];
}

// Per-DOF variant: pack reads d_field[nodeToDof[d_send_node_ids[i]]] so that
// per-DOF vectors (e.g. CG's Ap, sized numOwnedDofs) can be exchanged using
// the per-NODE periodic-pair table.
template<typename RealType>
__global__ void packCrossRankSendDofKernel(const int* d_send_node_ids,
                                           const int* d_nodeToDof,
                                           const RealType* d_dof_field,
                                           RealType* d_send_buf,
                                           int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    int dof = d_nodeToDof[d_send_node_ids[i]];
    d_send_buf[i] = (dof >= 0) ? d_dof_field[dof] : RealType(0);
}

template<typename RealType>
__global__ void atomicAddCrossRankRecvDofKernel(const int* d_recv_node_ids,
                                                const int* d_nodeToDof,
                                                const RealType* d_recv_buf,
                                                RealType* d_dof_field,
                                                int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    int dof = d_nodeToDof[d_recv_node_ids[i]];
    if (dof >= 0) atomicAdd(&d_dof_field[dof], d_recv_buf[i]);
}

template<typename RealType>
__global__ void overwriteCrossRankRecvDofKernel(const int* d_recv_node_ids,
                                                const int* d_nodeToDof,
                                                const RealType* d_recv_buf,
                                                RealType* d_dof_field,
                                                int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    int dof = d_nodeToDof[d_recv_node_ids[i]];
    if (dof >= 0) d_dof_field[dof] = d_recv_buf[i];
}

// atomicAdd: a single owned master can be the target of slaves from multiple
// peers (corner master with periodic incidences from 2-3 directions).
template<typename RealType>
__global__ void atomicAddCrossRankRecvKernel(const int* d_recv_ids,
                                             const RealType* d_recv_buf,
                                             RealType* d_field,
                                             int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    atomicAdd(&d_field[d_recv_ids[i]], d_recv_buf[i]);
}

template<typename RealType>
__global__ void overwriteCrossRankRecvKernel(const int* d_recv_ids,
                                             const RealType* d_recv_buf,
                                             RealType* d_field,
                                             int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    d_field[d_recv_ids[i]] = d_recv_buf[i];
}

template<typename RealType>
__global__ void zeroCrossRankSlavesKernel(const int* d_send_ids,
                                          RealType* d_field,
                                          int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    d_field[d_send_ids[i]] = RealType(0);
}

// Build the cross-rank periodic-pair table. Called from the end of
// buildPeriodicMap. Single-rank or no-cross-rank-pairs case leaves
// map.cross_ default-constructed (empty) and the per-step exchange becomes
// a no-op.
//
// Algorithm:
//   1. Walk owned slaves with REMOTE-owned master. Each is a (slave_local,
//      master_ghost) pair; we know master_ghost's SFC key (= the periodic
//      partner's image coord).
//   2. Use cstone's NodeHaloTopology to look up which peer rank delivers
//      master_ghost to us as a halo node. That peer is the master's owner.
//      Fallback: for ghosts NOT in the per-node halo (rare; periodic-only
//      seam delivered as element halo without per-node edge), broadcast the
//      key to all ranks.
//   3. MPI_Alltoallv exchanges (master_sfc_key) lists. Each receiver tries to
//      find each incoming key in its own owned-node SFC map. Hits become
//      "I own this master" entries; misses fall through.
//   4. MPI_Alltoallv exchanges ACK lists (sender-side slot indices that the
//      receiver claimed). Sender uses these to pick the local slave ids that
//      pair with each remote master.
//   5. Build final per-peer CSR send (owned slave ids) and recv (owned
//      master ids) device buffers, plus pre-sized staging.
template<typename KeyType, typename RealType, typename DomainT>
void buildCrossRankPeriodicMap(const DomainT& domain,
                               PeriodicMap<KeyType, RealType>& map,
                               MPI_Comm comm)
{
    using XR = CrossRankPeriodicMap<KeyType, RealType>;
    XR& xr = map.cross_;

    // Always rebuild from scratch (AMR may invoke us again).
    xr.peers_.clear();
    xr.sendOffsets_.assign(1, 0);
    xr.recvOffsets_.assign(1, 0);
    xr.d_sendOwnedSlaveIds_.resize(0);
    xr.d_recvOwnedMasterIds_.resize(0);

    if (domain.numRanks() <= 1) return;

    // First-build sub-communicator dup (reused across AMR rebuilds).
    if (xr.comm_ == MPI_COMM_NULL) { MPI_Comm_dup(comm, &xr.comm_); }

    const int    myRank   = domain.rank();
    const int    numRanks = domain.numRanks();
    const size_t numNodes = domain.getNodeCount();

    // Stage to host.
    std::vector<int>      h_partner(numNodes);
    std::vector<uint8_t>  h_own(numNodes);
    std::vector<KeyType>  h_sfc(numNodes);
    cudaMemcpy(h_partner.data(), map.d_periodicPartner.data(),
               numNodes * sizeof(int), cudaMemcpyDeviceToHost);
    {
        const auto& d_own = domain.getNodeOwnershipMap();
        cudaMemcpy(h_own.data(), d_own.data(),
                   numNodes * sizeof(uint8_t), cudaMemcpyDeviceToHost);
    }
    {
        const auto& d_sfc = domain.getLocalToGlobalSfcMap();
        cudaMemcpy(h_sfc.data(), d_sfc.data(),
                   numNodes * sizeof(KeyType), cudaMemcpyDeviceToHost);
    }

    // keyToOwned: SFC key -> local owned node id (for receive-side resolution).
    std::unordered_map<KeyType, int> keyToOwned;
    keyToOwned.reserve(numNodes);
    for (size_t n = 0; n < numNodes; ++n)
    {
        if (h_own[n] == 1) keyToOwned[h_sfc[n]] = int(n);
    }

    // Use cstone NodeHaloTopology to map ghost-local-id -> peer-index.
    const auto& topo = domain.getNodeHaloTopology();
    std::vector<int> h_recvNodeIds(topo.recvNodeIds_.size());
    if (!h_recvNodeIds.empty())
    {
        cudaMemcpy(h_recvNodeIds.data(), topo.recvNodeIds_.data(),
                   h_recvNodeIds.size() * sizeof(int), cudaMemcpyDeviceToHost);
    }
    std::vector<int> ghostToPeer(numNodes, -1);
    for (size_t p = 0; p < topo.peers_.size(); ++p)
    {
        int lo = topo.recvOffsets_[p];
        int hi = topo.recvOffsets_[p + 1];
        for (int j = lo; j < hi; ++j)
        {
            int g = h_recvNodeIds[j];
            if (g >= 0 && g < int(numNodes) && ghostToPeer[g] == -1)
                ghostToPeer[g] = int(p);
        }
    }

    // STEP 1: classify owned slaves.
    std::vector<std::vector<int>>     slavesByPeer(topo.peers_.size());
    std::vector<std::vector<KeyType>> masterKeysByPeer(topo.peers_.size());
    std::vector<int>     fallbackSlaves;
    std::vector<KeyType> fallbackMasterKeys;

    bool localChainTruncated = false;
    for (size_t i = 0; i < numNodes; ++i)
    {
        if (h_own[i] != 1) continue;
        int m = h_partner[i];
        if (m < 0)                 continue;  // not a slave
        if (h_own[m] == 1)         continue;  // same-rank pair (gated kernel handles)
        // Detect cross-rank corner chain (flatten stopped at a remote ghost
        // whose own partner is set on another rank).
        if (h_partner[m] != -1)    localChainTruncated = true;
        KeyType masterKey = h_sfc[m];
        int peerIdx = ghostToPeer[m];
        if (peerIdx >= 0)
        {
            slavesByPeer[peerIdx].push_back(int(i));
            masterKeysByPeer[peerIdx].push_back(masterKey);
        }
        else
        {
            fallbackSlaves.push_back(int(i));
            fallbackMasterKeys.push_back(masterKey);
        }
    }

    // Abort if any rank found a cross-rank corner chain (unsupported in v1).
    {
        int local = localChainTruncated ? 1 : 0;
        int global = 0;
        MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX, xr.comm_);
        if (global != 0)
        {
            if (myRank == 0)
            {
                std::cerr << "buildCrossRankPeriodicMap: cross-rank corner-chain"
                          << " detected (flatten stopped at a ghost whose own"
                          << " partner is set). Unsupported in this pass.\n";
            }
            MPI_Abort(xr.comm_, 73);
        }
    }

    // STEP 2: per-rank send counts. For cstone-peer ranks: use the bucket
    // count. For fallback: broadcast count to every other rank.
    std::vector<std::vector<KeyType>> perRankKeys(numRanks);
    std::vector<std::vector<int>>     perRankSlaves(numRanks);
    for (size_t p = 0; p < topo.peers_.size(); ++p)
    {
        int peer = topo.peers_[p];
        perRankKeys[peer].insert(perRankKeys[peer].end(),
                                 masterKeysByPeer[p].begin(),
                                 masterKeysByPeer[p].end());
        perRankSlaves[peer].insert(perRankSlaves[peer].end(),
                                   slavesByPeer[p].begin(),
                                   slavesByPeer[p].end());
    }
    if (!fallbackSlaves.empty())
    {
        for (int r = 0; r < numRanks; ++r)
        {
            if (r == myRank) continue;
            perRankKeys[r].insert(perRankKeys[r].end(),
                                  fallbackMasterKeys.begin(),
                                  fallbackMasterKeys.end());
            perRankSlaves[r].insert(perRankSlaves[r].end(),
                                    fallbackSlaves.begin(),
                                    fallbackSlaves.end());
        }
    }

    std::vector<int> sendCounts(numRanks, 0);
    for (int r = 0; r < numRanks; ++r) sendCounts[r] = int(perRankKeys[r].size());
    std::vector<int> recvCounts(numRanks, 0);
    MPI_Alltoall(sendCounts.data(), 1, MPI_INT,
                 recvCounts.data(), 1, MPI_INT, xr.comm_);

    std::vector<int> sendDispls(numRanks, 0), recvDispls(numRanks, 0);
    int sendTotal = 0, recvTotal = 0;
    for (int r = 0; r < numRanks; ++r)
    {
        sendDispls[r] = sendTotal; sendTotal += sendCounts[r];
        recvDispls[r] = recvTotal; recvTotal += recvCounts[r];
    }
    std::vector<KeyType> sendKeysFlat(sendTotal);
    for (int r = 0; r < numRanks; ++r)
    {
        std::copy(perRankKeys[r].begin(), perRankKeys[r].end(),
                  sendKeysFlat.begin() + sendDispls[r]);
    }
    std::vector<KeyType> recvKeysFlat(recvTotal);
    MPI_Datatype mpiKeyType = (sizeof(KeyType) == 8) ? MPI_UNSIGNED_LONG_LONG
                                                     : MPI_UNSIGNED;
    MPI_Alltoallv(sendKeysFlat.data(), sendCounts.data(), sendDispls.data(), mpiKeyType,
                  recvKeysFlat.data(), recvCounts.data(), recvDispls.data(), mpiKeyType,
                  xr.comm_);

    // STEP 3: receiver-side resolution.
    std::vector<std::vector<int>> ownedMastersFromRank(numRanks);
    std::vector<std::vector<int>> ackedSlotsFromRank(numRanks);
    for (int r = 0; r < numRanks; ++r)
    {
        int off = recvDispls[r];
        int cnt = recvCounts[r];
        ownedMastersFromRank[r].reserve(cnt);
        ackedSlotsFromRank[r].reserve(cnt);
        for (int k = 0; k < cnt; ++k)
        {
            auto it = keyToOwned.find(recvKeysFlat[off + k]);
            if (it != keyToOwned.end())
            {
                ownedMastersFromRank[r].push_back(it->second);
                ackedSlotsFromRank[r].push_back(k);
            }
        }
    }

    // STEP 4: ACK exchange (which slot indices each receiver claimed).
    std::vector<int> ackSendCounts(numRanks, 0);
    for (int r = 0; r < numRanks; ++r)
        ackSendCounts[r] = int(ackedSlotsFromRank[r].size());
    std::vector<int> ackRecvCounts(numRanks, 0);
    MPI_Alltoall(ackSendCounts.data(), 1, MPI_INT,
                 ackRecvCounts.data(), 1, MPI_INT, xr.comm_);
    std::vector<int> ackSendDispls(numRanks, 0), ackRecvDispls(numRanks, 0);
    int ackSendTotal = 0, ackRecvTotal = 0;
    for (int r = 0; r < numRanks; ++r)
    {
        ackSendDispls[r] = ackSendTotal; ackSendTotal += ackSendCounts[r];
        ackRecvDispls[r] = ackRecvTotal; ackRecvTotal += ackRecvCounts[r];
    }
    std::vector<int> ackSendFlat(ackSendTotal);
    for (int r = 0; r < numRanks; ++r)
    {
        std::copy(ackedSlotsFromRank[r].begin(), ackedSlotsFromRank[r].end(),
                  ackSendFlat.begin() + ackSendDispls[r]);
    }
    std::vector<int> ackRecvFlat(ackRecvTotal);
    MPI_Alltoallv(ackSendFlat.data(), ackSendCounts.data(), ackSendDispls.data(), MPI_INT,
                  ackRecvFlat.data(), ackRecvCounts.data(), ackRecvDispls.data(), MPI_INT,
                  xr.comm_);

    // STEP 5: build final per-peer CSR.
    std::vector<int> sendOwnedSlaveIdsHost;
    std::vector<int> recvOwnedMasterIdsHost;
    xr.peers_.clear();
    xr.sendOffsets_.assign(1, 0);
    xr.recvOffsets_.assign(1, 0);
    for (int r = 0; r < numRanks; ++r)
    {
        int nAcksFromR    = ackRecvCounts[r];
        int nMastersFromR = int(ownedMastersFromRank[r].size());
        if (nAcksFromR == 0 && nMastersFromR == 0) continue;
        xr.peers_.push_back(r);
        // Send list: slaves we own whose master r claimed.
        for (int k = 0; k < nAcksFromR; ++k)
        {
            int slotIdx = ackRecvFlat[ackRecvDispls[r] + k];
            sendOwnedSlaveIdsHost.push_back(perRankSlaves[r][slotIdx]);
        }
        xr.sendOffsets_.push_back(int(sendOwnedSlaveIdsHost.size()));
        // Recv list: local owned masters we resolved from r's keys.
        recvOwnedMasterIdsHost.insert(recvOwnedMasterIdsHost.end(),
                                       ownedMastersFromRank[r].begin(),
                                       ownedMastersFromRank[r].end());
        xr.recvOffsets_.push_back(int(recvOwnedMasterIdsHost.size()));
    }

    // Optional build-side symmetry check (gate G6).
    if (std::getenv("MARS_PERIODIC_XR_CHECK") != nullptr)
    {
        for (size_t p = 0; p < xr.peers_.size(); ++p)
        {
            int peer = xr.peers_[p];
            int sCnt = xr.sendOffsets_[p + 1] - xr.sendOffsets_[p];
            int rCntPeer = 0;
            MPI_Sendrecv(&sCnt, 1, MPI_INT, peer, 0xCAFE,
                         &rCntPeer, 1, MPI_INT, peer, 0xCAFE,
                         xr.comm_, MPI_STATUS_IGNORE);
            int expectedRecvFromPeer = xr.recvOffsets_[p + 1] - xr.recvOffsets_[p];
            if (rCntPeer != expectedRecvFromPeer)
            {
                std::cerr << "[periodic-xr rank " << myRank << "] peer=" << peer
                          << " send=" << sCnt << " peer-sends=" << rCntPeer
                          << " we-expect-recv=" << expectedRecvFromPeer << "\n";
                MPI_Abort(xr.comm_, 74);
            }
        }
    }

    // Copy CSR + pre-size staging.
    if (!sendOwnedSlaveIdsHost.empty())
    {
        xr.d_sendOwnedSlaveIds_.resize(sendOwnedSlaveIdsHost.size());
        cudaMemcpy(xr.d_sendOwnedSlaveIds_.data(), sendOwnedSlaveIdsHost.data(),
                   sendOwnedSlaveIdsHost.size() * sizeof(int),
                   cudaMemcpyHostToDevice);
    }
    if (!recvOwnedMasterIdsHost.empty())
    {
        xr.d_recvOwnedMasterIds_.resize(recvOwnedMasterIdsHost.size());
        cudaMemcpy(xr.d_recvOwnedMasterIds_.data(), recvOwnedMasterIdsHost.data(),
                   recvOwnedMasterIdsHost.size() * sizeof(int),
                   cudaMemcpyHostToDevice);
    }
    xr.sendBuf_.resize(sendOwnedSlaveIdsHost.size());
    xr.recvBuf_.resize(recvOwnedMasterIdsHost.size());
    cudaDeviceSynchronize();

    if (myRank == 0)
    {
        std::cout << "[periodic-xr] peers=" << xr.peers_.size()
                  << " send_total=" << sendOwnedSlaveIdsHost.size()
                  << " recv_total=" << recvOwnedMasterIdsHost.size() << "\n";
    }
}

// Build the periodic map from current domain coords. Call after mesh load
// and after every AMR rebuild. When comm != MPI_COMM_NULL the cross-rank
// periodic-pair tables are built too.
template<typename KeyType, typename RealType, typename DomainT>
void buildPeriodicMap(const DomainT& domain, PeriodicMap<KeyType, RealType>& map,
                      RealType xmin, RealType xmax, RealType ymin, RealType ymax,
                      RealType zmin, RealType zmax, RealType faceEps = 1e-5,
                      MPI_Comm comm = MPI_COMM_NULL)
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

    // Cross-rank pair tables: built only when a real MPI communicator was
    // passed AND there are multiple ranks. Single-rank and pre-MPI calls
    // leave map.cross_ empty.
    if (comm != MPI_COMM_NULL && domain.numRanks() > 1)
    {
        buildCrossRankPeriodicMap<KeyType, RealType>(domain, map, comm);
    }
}

// Cross-rank periodic-pair sum. Leg 1 ships d_field[owned_slave] to the
// master's owner rank and atomicAdds into d_field[owned_master]. Leg 2
// broadcasts the merged master value back to the slave's owner rank,
// overwriting d_field[owned_slave].
//
// Pre-condition: caller has already done the intra-rank gated pair sum
// (periodicPairSumKernel with ownership pointer) so same-rank pairs are
// already merged at the master slot. This function bridges the cross-rank
// gap that reverseExchangeNodeHaloAdd cannot (because slave and master
// have distinct SFC keys when periodicity is implemented via cstone Box).
//
// No-op when xr.peers_.empty() (single-rank or no cross-rank pairs).
template<typename KeyType, typename RealType>
void crossRankPeriodicPairSum(const PeriodicMap<KeyType, RealType>& map,
                              cstone::DeviceVector<RealType>& d_field)
{
    const auto& xr = map.cross_;
    if (xr.peers_.empty()) return;

    auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    const int numPeers = int(xr.peers_.size());
    int sendTotal = xr.sendOffsets_.back();
    int recvTotal = xr.recvOffsets_.back();

    // -------- Leg 1: SLAVE owner -> MASTER owner, atomicAdd --------
    if (sendTotal > 0)
    {
        int blk = 256, grd = (sendTotal + blk - 1) / blk;
        packCrossRankSendKernel<RealType><<<grd, blk>>>(
            xr.d_sendOwnedSlaveIds_.data(), d_field.data(),
            xr.sendBuf_.data(), sendTotal);
        cudaDeviceSynchronize();
    }

    {
        const int tag = 0x5041 + xr.epoch_;
        std::vector<MPI_Request> reqs; reqs.reserve(2 * numPeers);
        for (int p = 0; p < numPeers; ++p)
        {
            int peer = xr.peers_[p];
            int rcnt = xr.recvOffsets_[p + 1] - xr.recvOffsets_[p];
            if (rcnt > 0)
            {
                MPI_Request r;
                MPI_Irecv(xr.recvBuf_.data() + xr.recvOffsets_[p], rcnt,
                          mpiType, peer, tag, xr.comm_, &r);
                reqs.push_back(r);
            }
        }
        for (int p = 0; p < numPeers; ++p)
        {
            int peer = xr.peers_[p];
            int scnt = xr.sendOffsets_[p + 1] - xr.sendOffsets_[p];
            if (scnt > 0)
            {
                MPI_Request r;
                MPI_Isend(xr.sendBuf_.data() + xr.sendOffsets_[p], scnt,
                          mpiType, peer, tag, xr.comm_, &r);
                reqs.push_back(r);
            }
        }
        if (!reqs.empty()) MPI_Waitall(int(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
        ++xr.epoch_;
    }

    if (recvTotal > 0)
    {
        int blk = 256, grd = (recvTotal + blk - 1) / blk;
        atomicAddCrossRankRecvKernel<RealType><<<grd, blk>>>(
            xr.d_recvOwnedMasterIds_.data(), xr.recvBuf_.data(),
            d_field.data(), recvTotal);
        cudaDeviceSynchronize();
    }

    // Zero the owned-slave slot (mirrors intra-rank periodicPairSumKernel's
    // pass-A semantics: after accumulate, slave slot is zeroed so leg 2's
    // overwrite is clean).
    if (sendTotal > 0)
    {
        int blk = 256, grd = (sendTotal + blk - 1) / blk;
        zeroCrossRankSlavesKernel<RealType><<<grd, blk>>>(
            xr.d_sendOwnedSlaveIds_.data(), d_field.data(), sendTotal);
        cudaDeviceSynchronize();
    }

    // -------- Leg 2: MASTER owner -> SLAVE owner, overwrite --------
    if (recvTotal > 0)
    {
        int blk = 256, grd = (recvTotal + blk - 1) / blk;
        packCrossRankSendKernel<RealType><<<grd, blk>>>(
            xr.d_recvOwnedMasterIds_.data(), d_field.data(),
            xr.recvBuf_.data(), recvTotal);
        cudaDeviceSynchronize();
    }

    {
        const int tag = 0x5042 + xr.epoch_;
        std::vector<MPI_Request> reqs; reqs.reserve(2 * numPeers);
        // Roles swapped: recvBuf_ now holds master values (source); we send
        // them BACK to the slave owners. Receive into sendBuf_ on the slave
        // side because that's sized to send-totals (per-peer send offsets
        // mirror what we sent in leg 1).
        for (int p = 0; p < numPeers; ++p)
        {
            int peer = xr.peers_[p];
            int rcnt = xr.sendOffsets_[p + 1] - xr.sendOffsets_[p];
            if (rcnt > 0)
            {
                MPI_Request r;
                MPI_Irecv(xr.sendBuf_.data() + xr.sendOffsets_[p], rcnt,
                          mpiType, peer, tag, xr.comm_, &r);
                reqs.push_back(r);
            }
        }
        for (int p = 0; p < numPeers; ++p)
        {
            int peer = xr.peers_[p];
            int scnt = xr.recvOffsets_[p + 1] - xr.recvOffsets_[p];
            if (scnt > 0)
            {
                MPI_Request r;
                MPI_Isend(xr.recvBuf_.data() + xr.recvOffsets_[p], scnt,
                          mpiType, peer, tag, xr.comm_, &r);
                reqs.push_back(r);
            }
        }
        if (!reqs.empty()) MPI_Waitall(int(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
        ++xr.epoch_;
    }

    if (sendTotal > 0)
    {
        int blk = 256, grd = (sendTotal + blk - 1) / blk;
        overwriteCrossRankRecvKernel<RealType><<<grd, blk>>>(
            xr.d_sendOwnedSlaveIds_.data(), xr.sendBuf_.data(),
            d_field.data(), sendTotal);
        cudaDeviceSynchronize();
    }
}

// Per-DOF variant of crossRankPeriodicPairSum. Used by the velocity CG's
// spmvPostCallback to sum Ap across cross-rank periodic pairs after each
// SpMV: the matrix has half-strength rows on the seam (assembler couldn't
// route the periodic-image column entries cross-rank), so the resulting
// Ap slot is half what the merged equation requires. Summing slave_Ap +
// master_Ap and broadcasting back restores symmetry.
//
// d_dof_field is sized numOwnedDofs (CG's Ap layout). nodeToDof maps each
// owned periodic-pair NODE to its owned DOF index in [0, numOwnedDofs).
template<typename KeyType, typename RealType>
void crossRankPeriodicPairSumDof(const PeriodicMap<KeyType, RealType>& map,
                                 const int* d_nodeToDof,
                                 cstone::DeviceVector<RealType>& d_dof_field)
{
    const auto& xr = map.cross_;
    if (xr.peers_.empty()) return;

    auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    const int numPeers = int(xr.peers_.size());
    int sendTotal = xr.sendOffsets_.back();
    int recvTotal = xr.recvOffsets_.back();

    // ---- Leg 1: SLAVE owner -> MASTER owner (atomicAdd) ----
    if (sendTotal > 0)
    {
        int blk = 256, grd = (sendTotal + blk - 1) / blk;
        packCrossRankSendDofKernel<RealType><<<grd, blk>>>(
            xr.d_sendOwnedSlaveIds_.data(), d_nodeToDof,
            d_dof_field.data(), xr.sendBuf_.data(), sendTotal);
        cudaDeviceSynchronize();
    }

    {
        const int tag = 0x5043 + xr.epoch_;
        std::vector<MPI_Request> reqs; reqs.reserve(2 * numPeers);
        for (int p = 0; p < numPeers; ++p)
        {
            int peer = xr.peers_[p];
            int rcnt = xr.recvOffsets_[p + 1] - xr.recvOffsets_[p];
            if (rcnt > 0)
            {
                MPI_Request r;
                MPI_Irecv(xr.recvBuf_.data() + xr.recvOffsets_[p], rcnt,
                          mpiType, peer, tag, xr.comm_, &r);
                reqs.push_back(r);
            }
        }
        for (int p = 0; p < numPeers; ++p)
        {
            int peer = xr.peers_[p];
            int scnt = xr.sendOffsets_[p + 1] - xr.sendOffsets_[p];
            if (scnt > 0)
            {
                MPI_Request r;
                MPI_Isend(xr.sendBuf_.data() + xr.sendOffsets_[p], scnt,
                          mpiType, peer, tag, xr.comm_, &r);
                reqs.push_back(r);
            }
        }
        if (!reqs.empty()) MPI_Waitall(int(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
        ++xr.epoch_;
    }

    if (recvTotal > 0)
    {
        int blk = 256, grd = (recvTotal + blk - 1) / blk;
        atomicAddCrossRankRecvDofKernel<RealType><<<grd, blk>>>(
            xr.d_recvOwnedMasterIds_.data(), d_nodeToDof,
            xr.recvBuf_.data(), d_dof_field.data(), recvTotal);
        cudaDeviceSynchronize();
    }

    // ---- Leg 2: MASTER -> SLAVE (overwrite, slave inherits master's merged value) ----
    if (recvTotal > 0)
    {
        int blk = 256, grd = (recvTotal + blk - 1) / blk;
        packCrossRankSendDofKernel<RealType><<<grd, blk>>>(
            xr.d_recvOwnedMasterIds_.data(), d_nodeToDof,
            d_dof_field.data(), xr.recvBuf_.data(), recvTotal);
        cudaDeviceSynchronize();
    }

    {
        const int tag = 0x5044 + xr.epoch_;
        std::vector<MPI_Request> reqs; reqs.reserve(2 * numPeers);
        for (int p = 0; p < numPeers; ++p)
        {
            int peer = xr.peers_[p];
            int rcnt = xr.sendOffsets_[p + 1] - xr.sendOffsets_[p];
            if (rcnt > 0)
            {
                MPI_Request r;
                MPI_Irecv(xr.sendBuf_.data() + xr.sendOffsets_[p], rcnt,
                          mpiType, peer, tag, xr.comm_, &r);
                reqs.push_back(r);
            }
        }
        for (int p = 0; p < numPeers; ++p)
        {
            int peer = xr.peers_[p];
            int scnt = xr.recvOffsets_[p + 1] - xr.recvOffsets_[p];
            if (scnt > 0)
            {
                MPI_Request r;
                MPI_Isend(xr.recvBuf_.data() + xr.recvOffsets_[p], scnt,
                          mpiType, peer, tag, xr.comm_, &r);
                reqs.push_back(r);
            }
        }
        if (!reqs.empty()) MPI_Waitall(int(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
        ++xr.epoch_;
    }

    if (sendTotal > 0)
    {
        int blk = 256, grd = (sendTotal + blk - 1) / blk;
        overwriteCrossRankRecvDofKernel<RealType><<<grd, blk>>>(
            xr.d_sendOwnedSlaveIds_.data(), d_nodeToDof,
            xr.sendBuf_.data(), d_dof_field.data(), sendTotal);
        cudaDeviceSynchronize();
    }
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

    const auto& d_ownership = domain.getNodeOwnershipMap();

    // (a) sum ghost contributions into owners (cross-rank, standard cstone path)
    domain.reverseExchangeNodeHaloAdd(d_field, dofMapPtr);

    // (b) within each rank: master += slave (ownership-gated; cross-rank no-op here)
    periodicPairSumKernel<<<grid, block>>>(map.d_periodicPartner.data(),
                                           d_ownership.data(),
                                           numNodes, d_field.data());

    // (b') cross-rank pair sum (no-op when map.cross_.peers_.empty()).
    crossRankPeriodicPairSum<KeyType, RealType>(map, d_field);

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
