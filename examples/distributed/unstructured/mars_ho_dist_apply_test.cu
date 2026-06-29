// Multi-rank gate for the HIGH-ORDER matrix-free apply (Phase B step 3).
// Builds the distributed HO DOF (buildDistributed + resolveHoDofOwnership + HoHalo),
// then runs ONE distributed matvec with the same single-rank kernel and checks the
// cross-rank assembly via the proven invariant A*1 == 0:
//   u = 1 on OWNED DOF, 0 on ghost  ->  forward()  (ghosts <- owner's 1)
//   apply owned elements (elemDof from buildDistributed is owned-only)
//   reverseAdd()  (ghost contributions summed into owners)
//   max|y| over OWNED DOF must be ~0  (the unconstrained HO Laplacian sums every
//   row to zero; a wrong cross-rank assembly leaves O(1) residual at interfaces).
// Run: srun -N1 -n4 ./mars_ho_dist_apply_test --ncells=16 --p=2
//
// OPTIONAL GPU NUMBERING (flag-gated; host path is the default and is UNCHANGED):
//   MARS_HO_GPU_NUMBERING=1 ./mars_ho_dist_apply_test --ncells=16 --p=2
//   ./mars_ho_dist_apply_test --gpu-numbering --ncells=16 --p=2
// The GPU path (buildDistributedGpu) does the SAME numbering with thrust. The
// local DOF ids may differ from host by a permutation, so correctness is checked
// by permutation-INVARIANT quantities, not elemDof element-wise.
//
// SELF-CHECK (A/B host vs GPU on the same config; prints the invariants):
//   ./mars_ho_dist_apply_test --self-check --ncells=16 --p=2
// Asserts: numDof/nEdge/nFace equal, and the MULTISET of DofKeys equal (sort both,
// compare). Then runs the A.1 apply gate with the GPU numbering -> end-to-end check.

#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/utils/mars_generate_cube.hpp"
#include "backend/distributed/unstructured/fem/mars_ho_dof_handler.hpp"
#include "backend/distributed/unstructured/fem/mars_ho_dof_handler_gpu.hpp"
#include "backend/distributed/unstructured/fem/mars_ho_halo.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_basis.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_apply.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_matfree.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_matfree_shfl.hpp"

#include <cuda_runtime.h>
#include <mpi.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/count.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/functional.h>
#include <thrust/transform_reduce.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <map>
#include <string>
#include <vector>
#include <array>

using namespace mars;
using namespace mars::fem;

using KeyType = uint64_t; using RealType = double; using ElemTag = HexTag;
using Domain  = ElementDomain<ElemTag, RealType, KeyType, cstone::GpuTag>;

// Timed-loop length, overridable by --iters=N (default 50). Higher counts damp
// the per-matvec timer noise on the throughput sweep. Warm-ups stay fixed at 5.
static int g_iters = 50;

// One row of the p=1..8 throughput sweep, filled by runDistApply and printed as a
// single clean table by main after the loop. apply-only (comm-excluded) GDOF/s is
// the throughput number; gflops = flopPerDof * gdofs. flopPerDof is the verified
// store-d_G count 3*P*(P+1)^2*(12*(P+1)+7)/P^3 (matches the kernel loops exactly).
struct SweepRow {
    int    p          = 0;
    long   gOwnedDof  = 0;     // global unique owned DOF (the sweep's DOF total)
    double dofPerGpu  = 0.0;   // gOwnedDof / numRanks
    double gdofsApply = 0.0;   // apply-only throughput (GDOF/s)
    double flopPerDof = 0.0;
    double gflops     = 0.0;   // FP64 GFLOP/s = flopPerDof * gdofsApply
    bool   pass       = false;
};
static std::vector<SweepRow> g_sweepRows;

// One row of the --esweep elements-per-block table: throughput vs E at a fixed p, so
// the user can read off the occupancy sweet spot. smemKB is the store-d_G dynamic smem
// per block at that (p,E). a1 is the per-E A.1 gate value (must PASS for every row).
struct ESweepRow {
    int    p       = 0;
    int    eblock  = 0;
    double dofPerGpu = 0.0;
    double gdofsApply = 0.0;
    double smemKB  = 0.0;
    double a1      = 0.0;
    bool   pass    = false;
};
static std::vector<ESweepRow> g_esweepRows;
// The bounded candidate E set the apply is instantiated for (mirrors the header switch).
static const int kEblockCandidates[] = {2, 4, 8, 16, 32};

// ---- Piece 1: build the per-corner ownership/shared/gid tables ON DEVICE ----
// These replace the host loops in extractDistDof that D2H'd the ownership map + SFC map
// + recv/send node-id lists and rebuilt cornerOwner/sharedCorner/cornerGid on the CPU.
// Every input here is already device-resident in the cstone domain, so the result stays
// on device and feeds the device-input numbering with NO round-trip.

// cornerGid[i] = (long) SFC key[i]  (the device SFC map is the global id).
__global__ static void ho_build_corner_gid(const KeyType* __restrict__ d_sfc,
                                           long* __restrict__ d_gid, long n)
{
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) d_gid[i] = (long)d_sfc[i];
}

// cornerOwner[i] = owned(i) ? myRank : -1  (ghosts get their true owner below).
__global__ static void ho_init_corner_owner(const uint8_t* __restrict__ d_own,
                                            int* __restrict__ d_cowner, int myRank, long n)
{
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) d_cowner[i] = d_own[i] ? myRank : -1;
}

// Ghost ownership + shared flag from the node-halo RECV list: each recv entry nd is a
// ghost owned by peers_[p], where p is the CSR segment holding the entry. recvOffsets +
// peers are tiny host->device arrays; the per-entry peer index is an upper_bound over
// the offsets. Also marks recv nodes shared. (atomic not needed: each nd appears once.)
__global__ static void ho_fill_ghost_owner_shared(const int* __restrict__ d_recv,
                                                  const int* __restrict__ d_recvOff,
                                                  const int* __restrict__ d_peers,
                                                  int nPeers, long nrecv, long nodeCount,
                                                  int* __restrict__ d_cowner,
                                                  uint8_t* __restrict__ d_shared)
{
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nrecv) return;
    // segment p such that d_recvOff[p] <= i < d_recvOff[p+1]
    int lo = 0, hi = nPeers;
    while (lo < hi) { int mid = (lo + hi) >> 1; if (d_recvOff[mid + 1] <= i) lo = mid + 1; else hi = mid; }
    int nd = d_recv[i];
    if (nd >= 0 && nd < (int)nodeCount) { d_cowner[nd] = d_peers[lo]; d_shared[nd] = 1; }
}

// Mark SEND nodes shared too (owned nodes we ship to a peer). One entry per send slot.
__global__ static void ho_mark_send_shared(const int* __restrict__ d_send, long nsend,
                                           long nodeCount, uint8_t* __restrict__ d_shared)
{
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nsend) return;
    int nd = d_send[i];
    if (nd >= 0 && nd < (int)nodeCount) d_shared[nd] = 1;
}

// A.1-gate device helpers (zero-host scale path). The owned set is captured ONCE as a
// packed 1-bit-per-DOF mask (numDof/8 bytes, ~78 MB at 625M DOF -- vs the 2.5 GB int
// dofOwner), built from the resolved device dofOwner BEFORE ownDev is freed, and survives
// through the apply. The A.1 max + owned-count are then device reductions returning only a
// scalar to the host (for the print + MPI_Allreduce). No per-DOF D2H.

// bit d of d_mask set iff d_owner[d] == myRank. Packed 32 bits per word -> the mask is
// (numDof+31)/32 words = ceil(numDof/8) bytes. Word/bit indexing (d>>5, d&31) is used
// IDENTICALLY by the read below, so the layout is endian-independent.
__global__ static void ho_pack_owned_mask(const int* __restrict__ d_owner, long n,
                                          int myRank, unsigned* __restrict__ d_mask)
{
    long d = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (d >= n) return;
    if (d_owner[d] == myRank) atomicOr(d_mask + (d >> 5), 1u << (d & 31));
}

// device max over owned DOF of |d_y[d]|, selected by the packed owned mask. transform over
// [0,n) -> (owned ? |y| : 0), reduced with max. Returns the same scalar the host loop did.
struct OwnedAbsMax {
    const double*   y;
    const unsigned* mask;
    __host__ __device__ double operator()(long d) const {
        const bool owned = (mask[d >> 5] >> (d & 31)) & 1u;
        return owned ? fabs(y[d]) : 0.0;
    }
};

// Distributed DOF inputs pulled from the domain (host side).
struct DistDof {
    size_t nodeCount = 0, nOwnedElem = 0;
    std::vector<std::array<int,8>> elemCorners;   // OWNED elements only -> numbering + apply
    std::vector<long>    cornerGid;
    std::vector<int>     cornerOwner;
    std::vector<uint8_t> sharedCorner;
    std::vector<int>     elemOwner;               // owned: myRank (numbering provisional)
    std::vector<int>     peers;

    // ZERO-MPI ownership (MARS_HO_LOCAL_OWNERSHIP): cstone halo elements + their true
    // owners, kept SEPARATE so the numbering/apply stay owned-only. Empty otherwise.
    bool                           localOwn = false;
    std::vector<std::array<int,8>> haloElemCorners;
    std::vector<int>               haloElemOwner;

    // Device-resident node-halo CSR + cornerGid for the TARGETED P2P ownership resolve
    // (MARS_HO_GPU_OWNERSHIP_P2P). Pulled from the live topology in extractDistDof so the
    // resolve never re-touches the domain. Offsets stay host (CSR), node-id lists + gid
    // stay on device. Empty/null on single rank or when targeting is off.
    thrust::device_vector<long> d_cornerGid;
    thrust::device_vector<int>  d_cornerOwner;            // per-corner owning rank (targeting)
    thrust::device_vector<int>  d_recvNodeIds, d_sendNodeIds;
    std::vector<int>            recvOffsets, sendOffsets;

    // Device handles into the LIVE domain (valid for the run; the domain outlives the
    // apply). The metric corner gather reads the 8 connectivity columns + node coords
    // straight off these -- no coord D2H, no host corner pack. startE selects the owned
    // element slice [startE, startE+nOwnedElem) inside the full-range connectivity.
    size_t           startE = 0;
    const KeyType*   d_conn[8] = {nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr};
    const RealType*  d_nodeX = nullptr, *d_nodeY = nullptr, *d_nodeZ = nullptr;

    // DEVICE-INPUT numbering (Piece 1): when set, the numbering reads cornerGid/owner/
    // shared/elemOwner straight off these device arrays (built by device kernels) + the
    // device connectivity -- no host corner tables, no H2D. Empty on the host paths.
    bool                           deviceInputs = false;
    thrust::device_vector<uint8_t> d_sharedCorner;
    thrust::device_vector<int>     d_elemOwner;
};

// wantDeviceInputs: drive the GPU-native default path. When true the per-corner tables
// (gid/owner/shared) + elemOwner are built with DEVICE kernels straight off the domain's
// device ownership map, SFC map and node-halo recv/send lists, and the host corner tables
// + their D2H/host loops are SKIPPED entirely. The numbering then runs the device-input
// wrapper (no H2D). When false (MARS_HO_HOST / self-check / single rank / local-ownership)
// the host tables are built exactly as before for the host numbering/resolve fallbacks.
static DistDof extractDistDof(Domain& domain, int rank, int numRanks, bool wantDeviceInputs)
{
    DistDof D;
    MPI_Barrier(MPI_COMM_WORLD);
    double thalo0 = MPI_Wtime();
    domain.getNodeOwnershipMap();                 // trigger halo build -> post-sync count
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) { printf("[build]   node-halo build %.1fs (in extractDistDof)\n", MPI_Wtime()-thalo0); fflush(stdout); }
    const size_t nodeCount = domain.getNodeCount();
    const size_t startE = domain.startIndex(), endE = domain.endIndex();
    D.nodeCount = nodeCount; D.nOwnedElem = endE - startE;
    D.deviceInputs = wantDeviceInputs;

    const auto& conn = domain.getElementToNodeConnectivity();
    auto C = [&](int c)->const KeyType* {
        return thrust::raw_pointer_cast(
            (c==0?std::get<0>(conn):c==1?std::get<1>(conn):c==2?std::get<2>(conn):c==3?std::get<3>(conn):
             c==4?std::get<4>(conn):c==5?std::get<5>(conn):c==6?std::get<6>(conn):std::get<7>(conn)).data()); };
    // Device connectivity columns + owned-element start are kept for BOTH paths: the
    // metric corner gather (always device) and the device-input numbering pack.
    D.startE = startE;
    for (int c = 0; c < 8; ++c) D.d_conn[c] = C(c);

    if (wantDeviceInputs) {
        // ---- Piece 1: build cornerGid/cornerOwner/sharedCorner + elemOwner ON DEVICE ----
        const KeyType* d_sfc = thrust::raw_pointer_cast(domain.getLocalToGlobalSfcMap().data());
        const uint8_t* d_own = thrust::raw_pointer_cast(domain.getNodeOwnershipMap().data());
        const int blk = 256;
        auto grid = [&](long m){ return (unsigned)((m + blk - 1) / blk); };

        D.d_cornerGid.resize(nodeCount);
        ho_build_corner_gid<<<grid((long)nodeCount), blk>>>(
            d_sfc, thrust::raw_pointer_cast(D.d_cornerGid.data()), (long)nodeCount);

        D.d_cornerOwner.resize(nodeCount);
        ho_init_corner_owner<<<grid((long)nodeCount), blk>>>(
            d_own, thrust::raw_pointer_cast(D.d_cornerOwner.data()), rank, (long)nodeCount);

        D.d_sharedCorner.assign(nodeCount, 0);

        const auto& topo = domain.getNodeHaloTopology();
        D.peers = topo.peers_;
        D.recvOffsets = topo.recvOffsets_;
        D.sendOffsets = topo.sendOffsets_;
        size_t nrecv = topo.recvOffsets_.empty() ? 0 : (size_t)topo.recvOffsets_.back();
        size_t nsend = topo.sendOffsets_.empty() ? 0 : (size_t)topo.sendOffsets_.back();

        // device-to-device copies of the node-halo lists (also reused by the P2P resolve)
        D.d_recvNodeIds.resize(nrecv > 0 ? nrecv : 1);
        D.d_sendNodeIds.resize(nsend > 0 ? nsend : 1);
        if (nrecv) cudaMemcpy(thrust::raw_pointer_cast(D.d_recvNodeIds.data()),
                              thrust::raw_pointer_cast(topo.recvNodeIds_.data()),
                              nrecv * sizeof(int), cudaMemcpyDeviceToDevice);
        if (nsend) cudaMemcpy(thrust::raw_pointer_cast(D.d_sendNodeIds.data()),
                              thrust::raw_pointer_cast(topo.sendNodeIds_.data()),
                              nsend * sizeof(int), cudaMemcpyDeviceToDevice);

        // small CSR offsets + peer list -> device (MPI metadata; tiny, not avoidable)
        thrust::device_vector<int> d_recvOff(D.recvOffsets.begin(), D.recvOffsets.end());
        thrust::device_vector<int> d_peers(D.peers.begin(), D.peers.end());
        if (nrecv)
            ho_fill_ghost_owner_shared<<<grid((long)nrecv), blk>>>(
                thrust::raw_pointer_cast(D.d_recvNodeIds.data()),
                thrust::raw_pointer_cast(d_recvOff.data()),
                thrust::raw_pointer_cast(d_peers.data()), (int)D.peers.size(),
                (long)nrecv, (long)nodeCount,
                thrust::raw_pointer_cast(D.d_cornerOwner.data()),
                thrust::raw_pointer_cast(D.d_sharedCorner.data()));
        if (nsend)
            ho_mark_send_shared<<<grid((long)nsend), blk>>>(
                thrust::raw_pointer_cast(D.d_sendNodeIds.data()), (long)nsend,
                (long)nodeCount, thrust::raw_pointer_cast(D.d_sharedCorner.data()));

        D.d_elemOwner.assign(D.nOwnedElem, rank);   // owned elements: provisional myRank
        cudaDeviceSynchronize();

        // Node coords stay device-resident for the metric gather (no D2H).
        double tcc0 = MPI_Wtime();
        domain.cacheNodeCoordinates();
        if (rank == 0) { printf("[build]   cacheNodeCoords %.1fs (in extractDistDof)\n", MPI_Wtime()-tcc0); fflush(stdout); }
        D.d_nodeX = thrust::raw_pointer_cast(domain.getNodeX().data());
        D.d_nodeY = thrust::raw_pointer_cast(domain.getNodeY().data());
        D.d_nodeZ = thrust::raw_pointer_cast(domain.getNodeZ().data());
        return D;
    }

    // ---- HOST-INPUT path (fallback): build the host corner tables as before ----
    std::vector<KeyType> h_sfc(nodeCount);
    cudaMemcpy(h_sfc.data(), thrust::raw_pointer_cast(domain.getLocalToGlobalSfcMap().data()),
               nodeCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
    std::vector<uint8_t> h_own(nodeCount);
    cudaMemcpy(h_own.data(), thrust::raw_pointer_cast(domain.getNodeOwnershipMap().data()),
               nodeCount * sizeof(uint8_t), cudaMemcpyDeviceToHost);

    D.cornerOwner.assign(nodeCount, -1);
    for (size_t i = 0; i < nodeCount; ++i) if (h_own[i]) D.cornerOwner[i] = rank;
    D.sharedCorner.assign(nodeCount, 0);
    if (numRanks > 1) {
        const auto& topo = domain.getNodeHaloTopology();
        D.peers = topo.peers_;
        size_t nrecv = topo.recvOffsets_.empty() ? 0 : (size_t)topo.recvOffsets_.back();
        std::vector<int> h_recv(nrecv);
        if (nrecv) cudaMemcpy(h_recv.data(), thrust::raw_pointer_cast(topo.recvNodeIds_.data()),
                              nrecv * sizeof(int), cudaMemcpyDeviceToHost);
        size_t nsend = topo.sendOffsets_.empty() ? 0 : (size_t)topo.sendOffsets_.back();
        std::vector<int> h_send(nsend);
        if (nsend) cudaMemcpy(h_send.data(), thrust::raw_pointer_cast(topo.sendNodeIds_.data()),
                              nsend * sizeof(int), cudaMemcpyDeviceToHost);
        for (size_t p = 0; p < topo.peers_.size(); ++p)
            for (int i = topo.recvOffsets_[p]; i < topo.recvOffsets_[p+1]; ++i) {
                int nd = h_recv[i];
                if (nd >= 0 && nd < (int)nodeCount) D.cornerOwner[nd] = topo.peers_[p];
            }
        for (int i = 0; i < (int)nrecv; ++i) if (h_recv[i]>=0 && h_recv[i]<(int)nodeCount) D.sharedCorner[h_recv[i]] = 1;
        for (int i = 0; i < (int)nsend; ++i) if (h_send[i]>=0 && h_send[i]<(int)nodeCount) D.sharedCorner[h_send[i]] = 1;

        // Keep the node-halo CSR DEVICE-resident for the targeted P2P ownership resolve:
        // the resolve maps a shared key's corners -> per-corner peer sets straight from
        // these lists. Offsets are host (CSR); node-id arrays are device-to-device copies.
        D.recvOffsets = topo.recvOffsets_;
        D.sendOffsets = topo.sendOffsets_;
        D.d_recvNodeIds.resize(nrecv > 0 ? nrecv : 1);
        D.d_sendNodeIds.resize(nsend > 0 ? nsend : 1);
        // cstone DeviceVector exposes raw ptr + size (no thrust iterators) -> device-to-
        // device copy via raw pointers, matching the domain.hpp exchange pattern.
        if (nrecv) cudaMemcpy(thrust::raw_pointer_cast(D.d_recvNodeIds.data()),
                              thrust::raw_pointer_cast(topo.recvNodeIds_.data()),
                              nrecv * sizeof(int), cudaMemcpyDeviceToDevice);
        if (nsend) cudaMemcpy(thrust::raw_pointer_cast(D.d_sendNodeIds.data()),
                              thrust::raw_pointer_cast(topo.sendNodeIds_.data()),
                              nsend * sizeof(int), cudaMemcpyDeviceToDevice);
    }

    D.cornerGid.assign(nodeCount, 0);
    for (size_t i = 0; i < nodeCount; ++i) D.cornerGid[i] = (long)h_sfc[i];
    if (numRanks > 1) {
        D.d_cornerGid = thrust::device_vector<long>(D.cornerGid.begin(), D.cornerGid.end());
        // cornerOwner is finalized above (owned -> myRank, ghosts <- node-halo recv list);
        // the targeted resolve needs it on device to find an OWNED corner per shared key.
        D.d_cornerOwner = thrust::device_vector<int>(D.cornerOwner.begin(), D.cornerOwner.end());
    }

    // OWNED elements [startE,endE) -> the numbering + apply.
    D.elemCorners.resize(D.nOwnedElem);
    for (int c = 0; c < 8; ++c) {
        std::vector<KeyType> col(D.nOwnedElem);
        cudaMemcpy(col.data(), C(c) + startE, D.nOwnedElem * sizeof(KeyType), cudaMemcpyDeviceToHost);
        for (size_t e = 0; e < D.nOwnedElem; ++e) D.elemCorners[e][c] = (int)col[e];
    }
    D.elemOwner.assign(D.nOwnedElem, rank);

    // ZERO-MPI HO ownership (MARS_HO_LOCAL_OWNERSHIP): also pull cstone's HALO elements
    // (the element indices OUTSIDE [startE,endE) within [0,getElementCount())). Their
    // corners share the SAME local-node id space as owned corners (createElementToNodeLocalIdMap
    // fills all 8 columns over the full element range through the one shared SFC->localId
    // map), so the connectivity columns C(c) are valid over the full range and the halo
    // corner ids index straight into cornerGid/cornerOwner. The halo is a LOCAL view of
    // our peers' boundary elements -> we resolve shared edge/face owners with NO key
    // exchange. Default path leaves haloElem* empty (owned-only behaviour preserved).
    D.localOwn = (numRanks > 1) && std::getenv("MARS_HO_LOCAL_OWNERSHIP") != nullptr;
    if (D.localOwn) {
        const auto ranges = domain.haloElementRanges();   // {[0,startE),[endE,elemCount)}
        size_t nHalo = 0;
        for (auto& r : ranges) nHalo += r.second - r.first;
        D.haloElemCorners.resize(nHalo);
        D.haloElemOwner.resize(nHalo);
        for (int c = 0; c < 8; ++c) {
            const KeyType* col = C(c);
            size_t hi = 0;
            for (auto& r : ranges) {
                std::vector<KeyType> buf(r.second - r.first);
                cudaMemcpy(buf.data(), col + r.first, buf.size() * sizeof(KeyType), cudaMemcpyDeviceToHost);
                for (size_t i = 0; i < buf.size(); ++i) D.haloElemCorners[hi + i][c] = (int)buf[i];
                hi += buf.size();
            }
        }
        // Halo element owner = owner of its LOWEST-SFC corner. cstone places each element
        // by the argmin corner SFC key, so element-owning-rank == that corner's P1 owner.
        // The min-key corner is the argmin over all 8, NOT local slot 0 -- using slot 0
        // blindly mislabels elements. cornerOwner is the already-resolved local table
        // (ghosts filled from the node-halo recv list) -> pure local lookup, no MPI.
        for (size_t e = 0; e < nHalo; ++e) {
            int minC = D.haloElemCorners[e][0];
            for (int c = 1; c < 8; ++c)
                if (D.cornerGid[D.haloElemCorners[e][c]] < D.cornerGid[minC]) minC = D.haloElemCorners[e][c];
            D.haloElemOwner[e] = D.cornerOwner[minC];
        }
    }

    // Node coords stay DEVICE-resident: the metric corner gather reads them straight off
    // the domain. Keep the device pointers; no D2H. (cacheNodeCoordinates makes the
    // getNodeX/Y/Z() arrays valid + stable for the run.)
    domain.cacheNodeCoordinates();
    D.d_nodeX = thrust::raw_pointer_cast(domain.getNodeX().data());
    D.d_nodeY = thrust::raw_pointer_cast(domain.getNodeY().data());
    D.d_nodeZ = thrust::raw_pointer_cast(domain.getNodeZ().data());
    return D;
}

// Permutation-independent A/B comparison of the GPU vs host numbering. The local
// DOF ids legitimately differ (host = std::map insertion order, GPU = sorted-key
// order), so we compare the INVARIANTS that cross-rank correctness depends on:
// numDof/nEdge/nFace (the unique counts) and the MULTISET of DofKeys. Returns true
// if everything matches on this rank.
static bool selfCheckNumbering(const HODofHandler& host, const HODofHandler& gpu, int rank)
{
    bool ok = true;
    if (host.numDof != gpu.numDof || host.nEdge != gpu.nEdge || host.nFace != gpu.nFace) ok = false;
    if (rank == 0)
        printf("[self-check] numDof host=%ld gpu=%ld | nEdge host=%ld gpu=%ld | nFace host=%ld gpu=%ld  [%s]\n",
               host.numDof, gpu.numDof, host.nEdge, gpu.nEdge, host.nFace, gpu.nFace,
               (host.numDof==gpu.numDof && host.nEdge==gpu.nEdge && host.nFace==gpu.nFace) ? "MATCH" : "MISMATCH");

    if (host.dofKey.size() != gpu.dofKey.size()) {
        if (rank == 0) printf("[self-check] dofKey size differs host=%zu gpu=%zu\n",
                              host.dofKey.size(), gpu.dofKey.size());
        return false;
    }

    auto packed = [](const HODofHandler::DofKey& k) {
        return std::array<long,6>{ (long)k.kind, k.g0, k.g1, k.g2, k.g3, (long)k.pos }; };
    std::vector<std::array<long,6>> hk(host.dofKey.size()), gk(gpu.dofKey.size());
    for (size_t i = 0; i < host.dofKey.size(); ++i) hk[i] = packed(host.dofKey[i]);
    for (size_t i = 0; i < gpu.dofKey.size();  ++i) gk[i] = packed(gpu.dofKey[i]);
    std::sort(hk.begin(), hk.end());
    std::sort(gk.begin(), gk.end());

    long mismatch = 0; long firstBad = -1;
    for (size_t i = 0; i < hk.size(); ++i)
        if (hk[i] != gk[i]) { ++mismatch; if (firstBad < 0) firstBad = (long)i; }
    if (mismatch) ok = false;
    if (rank == 0)
        printf("[self-check] DofKey multiset: %ld / %zu mismatched (firstBad slot=%ld)  [%s]\n",
               mismatch, hk.size(), firstBad, mismatch == 0 ? "MATCH" : "MISMATCH");

    // dofShared / dofBoundary are per-DOF flags keyed on local ids -> not directly
    // comparable. But their TOTALS are permutation-invariant -> compare counts.
    long hShared=0,gShared=0,hBnd=0,gBnd=0;
    for (auto v : host.dofShared)   hShared += v;
    for (auto v : gpu.dofShared)    gShared += v;
    for (auto v : host.dofBoundary) hBnd += v;
    for (auto v : gpu.dofBoundary)  gBnd += v;
    if (hShared != gShared || hBnd != gBnd) ok = false;
    if (rank == 0)
        printf("[self-check] shared count host=%ld gpu=%ld | boundary count host=%ld gpu=%ld  [%s]\n",
               hShared, gShared, hBnd, gBnd, (hShared==gShared && hBnd==gBnd) ? "MATCH" : "MISMATCH");

    // Owner histogram (rank -> #DOF owned) is permutation-invariant too.
    {
        std::map<int,long> ho, go;
        for (int o : host.dofOwner) ho[o]++;
        for (int o : gpu.dofOwner)  go[o]++;
        bool ownMatch = (ho == go);
        if (!ownMatch) ok = false;
        if (rank == 0) printf("[self-check] owner histogram  [%s]\n", ownMatch ? "MATCH" : "MISMATCH");
    }
    return ok;
}

enum class Numbering { Host, Gpu, SelfCheck };

template<int P>
static void runDistApply(DistDof& D, int rank, int numRanks, Numbering mode)
{
    const int n = P + 1, N3 = n*n*n;

    HODofHandler dof;       // the handler actually used downstream (host or GPU built)
    HODofHandler dofHost;   // only populated in SelfCheck, for the A/B compare

    // GPU ownership resolve is now the DEFAULT: device-resident targeted P2P + GPUDirect,
    // verified bit-identical to the all-peer reference on regular AND irregular partitions,
    // ~0.7s at 625M DOF/rank (vs 19s host). Needs GPU numbering (the device columns) + >1 rank.
    // Debug overrides: MARS_HO_HOST forces the host resolve; MARS_HO_GPU_OWNERSHIP forces the
    // GPU all-peer ground truth; MARS_HO_VERIFY_OWNERSHIP runs targeted+all-peer and aborts on
    // mismatch. (MARS_HO_GPU_OWNERSHIP_P2P stays accepted but is redundant -- P2P is default.)
    const bool wantHostOwn = std::getenv("MARS_HO_HOST") != nullptr;
    const bool wantGpuOwn  = std::getenv("MARS_HO_GPU_OWNERSHIP") != nullptr;
    const bool wantVerify  = std::getenv("MARS_HO_VERIFY_OWNERSHIP") != nullptr;
    // targeted P2P = default whenever GPU numbering is in play and no override is set
    const bool wantP2P     = (mode == Numbering::Gpu) && !wantHostOwn && !wantGpuOwn && !wantVerify;
    const bool needOwnDev  = (mode == Numbering::Gpu) && !wantHostOwn;
    HoOwnershipDeviceData ownDev;
    bool haveOwnDev = false;

    if (mode == Numbering::Host) {
        if (rank == 0) { printf("[build] numbering = HOST (--host-numbering / MARS_HO_HOST_NUMBERING)\n"); fflush(stdout); }
        double tb0 = MPI_Wtime();
        dof.buildDistributed(D.elemCorners, (long)D.nodeCount, P, D.cornerGid, D.cornerOwner, D.elemOwner, rank, D.sharedCorner);
        double tb1 = MPI_Wtime();
        if (rank == 0) { printf("[build] buildDistributed (host) %.3fs (numDof=%ld)\n", tb1-tb0, dof.numDof); fflush(stdout); }
    } else if (mode == Numbering::Gpu) {
        cudaDeviceSynchronize();
        double tb0 = MPI_Wtime();
        haveOwnDev = needOwnDev && numRanks > 1;
        if (D.deviceInputs) {
            // DEVICE-INPUT numbering: corner tables + connectivity already device-resident
            // (Piece 1) -> zero H2D. Same core -> bit-identical numbering / A.1.
            if (rank == 0) { printf("[build] numbering = GPU (device inputs, no H2D)\n"); fflush(stdout); }
            buildDistributedGpuDevice<KeyType>(dof, D.d_conn, (long)D.startE, (long)D.nOwnedElem,
                (long)D.nodeCount, P,
                thrust::raw_pointer_cast(D.d_cornerGid.data()),
                thrust::raw_pointer_cast(D.d_cornerOwner.data()),
                thrust::raw_pointer_cast(D.d_sharedCorner.data()),
                thrust::raw_pointer_cast(D.d_elemOwner.data()),
                rank, haveOwnDev ? &ownDev : nullptr);
        } else {
            if (rank == 0) { printf("[build] numbering = GPU (host inputs)\n"); fflush(stdout); }
            buildDistributedGpu(dof, D.elemCorners, (long)D.nodeCount, P, D.cornerGid, D.cornerOwner, D.elemOwner, rank, D.sharedCorner,
                                haveOwnDev ? &ownDev : nullptr);
        }
        cudaDeviceSynchronize();
        double tb1 = MPI_Wtime();
        if (rank == 0) { printf("[build] buildDistributedGpu %.3fs (numDof=%ld)\n", tb1-tb0, dof.numDof); fflush(stdout); }
    } else { // SelfCheck: build BOTH, time each, then compare invariants. Use GPU for downstream.
        if (rank == 0) { printf("[build] numbering = SELF-CHECK (host vs GPU A/B)\n"); fflush(stdout); }
        double th0 = MPI_Wtime();
        dofHost.buildDistributed(D.elemCorners, (long)D.nodeCount, P, D.cornerGid, D.cornerOwner, D.elemOwner, rank, D.sharedCorner);
        double th1 = MPI_Wtime();
        cudaDeviceSynchronize();
        double tg0 = MPI_Wtime();
        buildDistributedGpu(dof, D.elemCorners, (long)D.nodeCount, P, D.cornerGid, D.cornerOwner, D.elemOwner, rank, D.sharedCorner);
        cudaDeviceSynchronize();
        double tg1 = MPI_Wtime();
        double thMax=0, tgMax=0, dh=th1-th0, dg=tg1-tg0;
        MPI_Allreduce(&dh, &thMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&dg, &tgMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (rank == 0)
            printf("[self-check] build time: host %.3fs  gpu %.3fs  speedup %.2fx\n",
                   thMax, tgMax, tgMax > 0 ? thMax/tgMax : 0.0);
        bool ok = selfCheckNumbering(dofHost, dof, rank);
        int allOk = ok ? 1 : 0, gAllOk = 1;
        MPI_Allreduce(&allOk, &gAllOk, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        if (rank == 0)
            printf("[self-check] numbering invariants (all ranks): %s -- now running A.1 gate with GPU numbering\n",
                   gAllOk ? "PASS" : "FAIL");
    }

    // Auto-skip the host A.1 cross-check at scale (it costs minutes past ~50M DOF/rank, and
    // needs the host dofOwner). Known here (depends only on numDof + env) so the resolved
    // dofOwner D2H below is gated by it -- at scale the device A.1 is the only gate, and it
    // reads the resolved owner straight from device (ownDev.dofOwner / the owned bitmask).
    const bool skipHostGate = std::getenv("MARS_HO_SKIP_HOST_GATE") != nullptr || dof.numDof > 50000000;

    double tb1 = MPI_Wtime();
    if (D.localOwn) {
        // ZERO-MPI resolve: reuse cstone's per-element halo. Feed owned + halo elements
        // (owned first, halo appended) with their owners -- each rank reads its peers'
        // boundary-element contributions LOCALLY, so no key exchange. Result is bit-
        // identical to the host resolveHoDofOwnership (min-rank-among-holders).
        std::vector<std::array<int,8>> allCorners = D.elemCorners;
        std::vector<int>               allOwner   = D.elemOwner;
        allCorners.insert(allCorners.end(), D.haloElemCorners.begin(), D.haloElemCorners.end());
        allOwner.insert(allOwner.end(),     D.haloElemOwner.begin(),   D.haloElemOwner.end());
        resolveHoDofOwnershipLocal(allCorners, D.cornerGid, allOwner, P,
                                   dof.dofShared, dof.dofKey, dof.dofOwner, rank);
        if (rank == 0) {
            printf("[build] resolveOwnership = LOCAL zero-MPI (MARS_HO_LOCAL_OWNERSHIP), halo elems=%zu\n",
                   D.haloElemCorners.size()); fflush(stdout);
        }
    } else if (haveOwnDev) {
        // device columns are resident in ownDev; atomicMin onto its d_dofOwner, then copy
        // ONLY dofOwner back (the halo build still runs on host). Three device routes:
        //   VERIFY -> targeted + all-peer, assert bit-identical (MPI_Abort on mismatch)
        //   P2P    -> targeted only (drops the all-peer fan-out volume)
        //   else   -> all-peer ground truth (MARS_HO_GPU_OWNERSHIP)
        int* d_own = thrust::raw_pointer_cast(ownDev.dofOwner.data());
        const long* d_kind = thrust::raw_pointer_cast(ownDev.dofKind.data());
        const long* d_g0 = thrust::raw_pointer_cast(ownDev.dofG0.data());
        const long* d_g1 = thrust::raw_pointer_cast(ownDev.dofG1.data());
        const long* d_g2 = thrust::raw_pointer_cast(ownDev.dofG2.data());
        const long* d_g3 = thrust::raw_pointer_cast(ownDev.dofG3.data());
        const int*  d_pos = thrust::raw_pointer_cast(ownDev.dofPos.data());
        const int*  d_sh = thrust::raw_pointer_cast(ownDev.dofShared.data());
        if (wantVerify) {
            resolveHoDofOwnershipGpuVerify(ownDev.numDof, d_kind, d_g0, d_g1, d_g2, d_g3, d_pos, d_sh, d_own,
                rank, D.peers, thrust::raw_pointer_cast(D.d_cornerGid.data()),
                thrust::raw_pointer_cast(D.d_cornerOwner.data()), (int)D.nodeCount,
                thrust::raw_pointer_cast(D.d_recvNodeIds.data()), D.recvOffsets,
                thrust::raw_pointer_cast(D.d_sendNodeIds.data()), D.sendOffsets);
            if (rank == 0) { printf("[build] resolveOwnership = GPU VERIFY (targeted vs all-peer)\n"); fflush(stdout); }
        } else if (wantP2P) {
            resolveHoDofOwnershipGpuTargeted(ownDev.numDof, d_kind, d_g0, d_g1, d_g2, d_g3, d_pos, d_sh, d_own,
                rank, D.peers, thrust::raw_pointer_cast(D.d_cornerGid.data()),
                thrust::raw_pointer_cast(D.d_cornerOwner.data()), (int)D.nodeCount,
                thrust::raw_pointer_cast(D.d_recvNodeIds.data()), D.recvOffsets,
                thrust::raw_pointer_cast(D.d_sendNodeIds.data()), D.sendOffsets);
            if (rank == 0) { printf("[build] resolveOwnership = GPU TARGETED P2P (default)\n"); fflush(stdout); }
        } else {
            resolveHoDofOwnershipGpu(ownDev.numDof, d_kind, d_g0, d_g1, d_g2, d_g3, d_pos, d_sh, d_own, rank, D.peers);
            if (rank == 0) { printf("[build] resolveOwnership = GPU all-peer (MARS_HO_GPU_OWNERSHIP)\n"); fflush(stdout); }
        }
        // Resolved owner D2H ONLY for the small-scale host A.1 cross-check (it loops host
        // dofOwner). At scale (skipHostGate) the host vector is left stale -- the device A.1
        // reads the resolved owner from ownDev.dofOwner / the owned bitmask, so the ~2.5 GB
        // int D2H is skipped. ownDev.dofOwner stays the source of truth either way.
        if (!skipHostGate)
            thrust::copy(ownDev.dofOwner.begin(), ownDev.dofOwner.end(), dof.dofOwner.begin());
        // ownDev (resolved dofOwner + the key columns) is KEPT RESIDENT through the device
        // halo build below -- buildDevice keys the boundary DOF + matches request keys from
        // these same columns with no re-upload. The D.d_* corner/halo columns the resolve
        // used are done now; free them here (they are NOT needed by the halo build). The
        // ownDev free is deferred until AFTER buildDevice, before the apply allocates d_G.
        D.d_cornerGid   = thrust::device_vector<long>();
        D.d_cornerOwner = thrust::device_vector<int>();
        D.d_recvNodeIds = thrust::device_vector<int>();
        D.d_sendNodeIds = thrust::device_vector<int>();
    } else {
        resolveHoDofOwnership(dof.dofShared, dof.dofKey, dof.dofOwner, rank, D.peers);
    }
    double tb2 = MPI_Wtime();
    if (rank == 0) { printf("[build] resolveOwnership %.3fs\n", tb2-tb1); fflush(stdout); }

    // HALO build. DEVICE buildDevice is the default whenever the GPU numbering kept the
    // device columns (haveOwnDev) and MARS_HO_HOST is not set -> the per-peer send/recv
    // lists are constructed ON DEVICE (receiver-driven, GPUDirect), bit-identical to host
    // build(), and stay device-resident so forwardDevice/reverseAddDevice need no
    // uploadDevice. MARS_HO_HOST (wantHostOwn) keeps the host build() (the fallback).
    HoHalo<RealType> halo;
    const bool haloOnDevice = haveOwnDev && !wantHostOwn;
    if (haloOnDevice) {
        halo.buildDevice(ownDev.numDof,
                         thrust::raw_pointer_cast(ownDev.dofKind.data()),
                         thrust::raw_pointer_cast(ownDev.dofG0.data()),
                         thrust::raw_pointer_cast(ownDev.dofG1.data()),
                         thrust::raw_pointer_cast(ownDev.dofG2.data()),
                         thrust::raw_pointer_cast(ownDev.dofG3.data()),
                         thrust::raw_pointer_cast(ownDev.dofPos.data()),
                         thrust::raw_pointer_cast(ownDev.dofOwner.data()),
                         thrust::raw_pointer_cast(ownDev.dofBoundary.data()),
                         rank, D.peers);
        if (rank == 0) { printf("[build] HoHalo = GPU buildDevice (device-resident lists)\n"); fflush(stdout); }
    } else {
        halo.build((int)dof.numDof, dof.dofOwner, dof.dofKey, dof.dofBoundary, rank, D.peers);
    }
    // elemDof stays DEVICE-resident: move it OUT of ownDev before the free below so the
    // apply reads it directly (no host download in buildDistributedGpu, no H2D here). The
    // apply needs a d_elemDof of exactly this size anyway, so this is not extra peak HBM.
    thrust::device_vector<int> d_elemDofDev;
    if (haveOwnDev) d_elemDofDev = std::move(ownDev.elemDof);

    // Owned set captured ON DEVICE before the ownDev free, as a packed 1-bit-per-DOF mask
    // (ceil(numDof/32) words = ceil(numDof/8) bytes, ~78 MB at 625M DOF -- NOT the 2.5 GB
    // int dofOwner). Tiny vs d_G, so it survives the ownDev free WITHOUT touching the
    // ownDev-before-d_G ceiling. The device A.1 readback selects owned DOF from this mask.
    // locOwned (this rank's owned-DOF count) is a thrust::count over the same device owner,
    // also before the free, feeding gOwned's MPI_Allreduce. On the host fallback (no ownDev)
    // both are derived from the host dofOwner instead (small-scale only).
    thrust::device_vector<unsigned> d_ownedMask;
    long locOwned = 0;
    if (haveOwnDev) {
        const long nd = ownDev.numDof;
        const int* d_own = thrust::raw_pointer_cast(ownDev.dofOwner.data());
        d_ownedMask.assign((nd + 31) / 32, 0u);
        const int blk = 256;
        ho_pack_owned_mask<<<(unsigned)((nd + blk - 1) / blk), blk>>>(
            d_own, nd, rank, thrust::raw_pointer_cast(d_ownedMask.data()));
        locOwned = (long)thrust::count(ownDev.dofOwner.begin(), ownDev.dofOwner.end(), rank);
        cudaDeviceSynchronize();
    }

    // Free the resolve's device columns now -- AFTER the halo build, BEFORE the apply
    // allocates d_G. ownDev is ~53 B/DOF (~33 GB at 625M/rank); leaving it resident blows
    // the per-GPU ceiling (apply OOM). The halo's own lists (d_sendDof_/d_recvDof_) are
    // separate device allocations and survive this free.
    ownDev = HoOwnershipDeviceData{};
    double tb3 = MPI_Wtime();
    if (rank == 0) { printf("[build] HoHalo %.3fs\n", tb3-tb2); fflush(stdout); }

    const size_t nEl  = D.nOwnedElem;
    const long   nDof = dof.numDof;

    auto op = buildHoCvfemOperators(P);
    ho_cvfem_upload_operators(P, op.Btil.data(), op.Dtil.data(),
                              op.D.data(), op.W.data(), op.xi.data(), op.zeta.data());

    // RECOMPUTE path (MARS_HO_RECOMPUTE / --recompute): skip the d_G array AND the
    // metric precompute -- the apply recomputes the per-SCS-point metric inline from the
    // 8 corner coords each matvec. d_corners is kept either way (store-d_G needs it to
    // FILL d_G; recompute reads it every matvec), so the gather is unconditional. The
    // store-d_G path is UNCHANGED and stays the default.
    const bool wantRecompute = std::getenv("MARS_HO_RECOMPUTE") != nullptr;
    // OPTIONAL fp32 metric storage (MARS_HO_FP32_METRIC): store d_G as float instead
    // of double -- d_G is ~75% of the apply's DRAM stream, so halving it is the top
    // throughput lever on the memory-bound apply (~+12-15% apply GDOF/s). Compute stays
    // fp64 (the float metric is promoted to double before the FMAs). Only the STORE-d_G
    // path has a d_G to narrow; recompute forms the metric inline in double and ignores
    // this flag. Default unset = fp64, bit-identical to before.
    const bool wantFp32Metric = !wantRecompute && std::getenv("MARS_HO_FP32_METRIC") != nullptr;
    // MARS_HO_SHFL: route the apply through the register+warp-shuffle kernel (mars_cvfem_ho_matfree_shfl.hpp).
    // With MARS_HO_RECOMPUTE it picks the MF-shuffle variant (same sweep, metric regenerated inline from
    // d_corners, no d_G); alone it is the store-d_G shuffle. fp64 only -- ignored under fp32.
    const bool wantShfl = !wantFp32Metric && std::getenv("MARS_HO_SHFL") != nullptr;
    // MARS_HO_NEK: register-line (hoisted normal column) variant of the shuffle, for an A/B vs SHFL.
    const bool wantNek = !wantRecompute && !wantFp32Metric && std::getenv("MARS_HO_NEK") != nullptr;
    if (rank == 0) {
        printf("[apply] operator metric path = %s  [metric precision = %s]\n",
               wantRecompute ? "RECOMPUTE (no d_G, inline metric)" : "STORE-d_G (default)",
               wantRecompute ? "fp64 (inline)" : (wantFp32Metric ? "fp32 (d_G=float)" : "fp64 (d_G=double)"));
        fflush(stdout);
    }
    if (rank == 0 && wantShfl) {
        printf("[apply]   MARS_HO_SHFL -> register+warp-shuffle apply%s (all p=1..8: p<=4 single-warp, p>=5 padded multi-warp)\n",
               wantRecompute ? " [MF: metric recomputed inline, no d_G]" : "");
        fflush(stdout);
    }
    if (rank == 0 && wantNek) {
        printf("[apply]   MARS_HO_NEK -> register-line (hoisted normal column), FULL sweep p=1..8 (p<=4 single-warp, p>=5 multi-warp)\n");
        fflush(stdout);
    }

    // elemDof: device-resident (GPU default path) -> use it straight, no alloc/H2D. Only
    // the host fallback (host numbering / self-check) allocates a d_elemDof and uploads it.
    const bool elemDofResident = (d_elemDofDev.size() == (size_t)nEl * N3);
    int* d_elemDof = elemDofResident ? thrust::raw_pointer_cast(d_elemDofDev.data()) : nullptr;
    double *d_corners=nullptr, *d_G=nullptr, *d_u=nullptr, *d_y=nullptr;
    // fp32 metric path stores d_G as float (half the bytes). Only one of d_G / d_G32 is
    // ever allocated -- wantFp32Metric selects which, and the alloc size below picks
    // sizeof(float) vs sizeof(double) to match.
    float *d_G32=nullptr;
    const size_t gLen = nEl * (size_t)(3 * P * n * n) * 3;
    const size_t gElemBytes = wantFp32Metric ? sizeof(float) : sizeof(double);
    // Checked allocation: report OOM as the per-GPU ceiling, not a null-deref crash.
    // (-n1 ceiling probe, or balanced runs where all ranks OOM together.) d_G is NOT
    // allocated on the recompute path -- that is the whole point (it frees the biggest
    // allocation, so the per-GPU DOF ceiling rises).
    cudaError_t me = cudaSuccess;
    MPI_Barrier(MPI_COMM_WORLD); double t_a0 = MPI_Wtime();
    if (!elemDofResident) me = cudaMalloc(&d_elemDof, sizeof(int) * nEl * N3);
    if (me==cudaSuccess) me = cudaMalloc(&d_corners, sizeof(double) * nEl * 24);
    double t_dG0 = MPI_Wtime();
    // d_G alloc size follows the metric precision: float (fp32) or double (default). The
    // recompute path allocates neither -- that is the whole point (frees the biggest array).
    if (me==cudaSuccess && !wantRecompute) {
        if (wantFp32Metric) me = cudaMalloc(&d_G32, gElemBytes * gLen);
        else                me = cudaMalloc(&d_G,   gElemBytes * gLen);
    }
    double t_dG1 = MPI_Wtime();
    if (me==cudaSuccess) me = cudaMalloc(&d_u,       sizeof(double) * nDof);
    if (me==cudaSuccess) me = cudaMalloc(&d_y,       sizeof(double) * nDof);
    if (rank == 0 && me==cudaSuccess) {
        printf("[apply]   alloc %.1fs total  (d_G %.0f GB malloc = %.1fs)\n",
               MPI_Wtime()-t_a0, gElemBytes*(double)gLen/1e9, t_dG1-t_dG0); fflush(stdout); }
    if (me != cudaSuccess) {
        const double gBytes = wantRecompute ? 0.0 : gElemBytes*(double)gLen;
        double needGB = (gBytes + sizeof(double)*2.0*nDof
                         + sizeof(int)*(double)nEl*N3 + sizeof(double)*(double)nEl*24) / 1e9;
        // OOM on ANY rank must kill the whole job, not return: the MPI_Barrier below would
        // otherwise deadlock forever on the rank that bailed out here (this is what hung the
        // 448-node trillion run for hours -- an imbalanced rank OOM'd on d_G, returned, and
        // rank 0 sat at the barrier). Print from the FAILING rank so the offenders are visible.
        fprintf(stderr, "[apply] rank %d DEVICE OOM (needs ~%.1f GB, d_G=%.1f GB, nDof=%ld nEl=%zu): %s"
                        "  -- per-GPU ceiling exceeded\n",
                rank, needGB, gBytes/1e9, nDof, nEl, cudaGetErrorString(me));
        fflush(stderr);
        cudaGetLastError();
        if (!elemDofResident) cudaFree(d_elemDof);
        cudaFree(d_corners); cudaFree(d_G); cudaFree(d_G32); cudaFree(d_u); cudaFree(d_y);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (!elemDofResident)
        cudaMemcpy(d_elemDof, dof.elemDof.data(), sizeof(int) * nEl * N3, cudaMemcpyHostToDevice);
    // Gather the metric corner coords ON DEVICE from the domain's connectivity + coords
    // (no host corner pack, no coord D2H, no H2D). The device handles are set for every
    // rank in extractDistDof straight off the live domain. Both paths need d_corners.
    MPI_Barrier(MPI_COMM_WORLD); double t_g0 = MPI_Wtime();
    ho_gather_corner_coords_launch<double, KeyType>(
        D.d_conn[0], D.d_conn[1], D.d_conn[2], D.d_conn[3],
        D.d_conn[4], D.d_conn[5], D.d_conn[6], D.d_conn[7],
        D.d_nodeX, D.d_nodeY, D.d_nodeZ, d_corners, D.startE, nEl);
    cudaDeviceSynchronize();
    if (rank == 0) { printf("[apply]   corner gather %.1fs\n", MPI_Wtime()-t_g0); fflush(stdout); }
    // Store-d_G ONLY: fill d_G once from the corners (the apply then reads it). The
    // recompute path skips this -- it forms the metric inside the apply each matvec.
    double t_p0 = MPI_Wtime();
    if (!wantRecompute) {
        // Write the metric at the selected precision into the matching buffer. The float
        // launcher narrows the store; geometry math inside the kernel stays double.
        if (wantFp32Metric) ho_cvfem_metric_perpoint_launch<double, P, float >(d_corners, d_G32, nEl);
        else                ho_cvfem_metric_perpoint_launch<double, P, double>(d_corners, d_G,   nEl);
    }
    cudaDeviceSynchronize();
    if (rank == 0) { printf("[apply]   metric precompute %.1fs\n", MPI_Wtime()-t_p0); fflush(stdout); }

    // One apply dispatch for every call site below (gates + timed loops): store-d_G reads
    // d_G, recompute reads d_corners and forms the metric inline. Same signature shape, so
    // the gate/timing code is identical on both paths -- the timed loop therefore measures
    // whichever path this run selected (the head-to-head is two runs at the same p).
    // d_elemList/count select a SUBSET of elements (interior or boundary) for the overlap
    // matvec; nullptr/0 = the full owned range (the original, bit-identical dispatch).
    // RUNTIME elements-per-block override (MARS_HO_EBLOCK in {2,4,8,16,32}). UNSET (or any
    // value outside the set) -> the per-P default, bit-identical to before. Only the
    // store-d_G fp64 default path is overridable (the E-sweep is store-d_G/Laplacian-only,
    // matching the apply the ncu occupancy study profiled); recompute/fp32 keep the default
    // E. The device's true optin smem cap feeds the feasibility guard so an over-smem E is
    // rejected (default launch) rather than launched. dev/cap queried once per order.
    int g_eblock = 0;
    if (const char* e = std::getenv("MARS_HO_EBLOCK")) g_eblock = std::atoi(e);
    int dev = 0; cudaGetDevice(&dev);
    int smemOptinB = 0;
    cudaDeviceGetAttribute(&smemOptinB, cudaDevAttrMaxSharedMemoryPerBlockOptin, dev);
    const size_t smemCap = (size_t)smemOptinB;
    const bool eblockOverridable = !wantRecompute && !wantFp32Metric;
    if (rank == 0 && g_eblock != 0 && eblockOverridable) {
        const bool feas = ho_cvfem_eblock_feasible<P>(g_eblock, smemCap);
        printf("[apply]   MARS_HO_EBLOCK=%d (store-d_G) -> %s (smem=%.1f KB/block, cap=%.0f KB)\n",
               g_eblock, feas ? "applied" : "INFEASIBLE -> falling back to default E",
               ho_cvfem_smem_bytes<P>(g_eblock)/1024.0, smemCap/1024.0);
        fflush(stdout);
    }
    auto applyOp = [&](double* u, double* y, const int* d_elemList = nullptr, size_t count = 0) {
        if (wantRecompute && wantShfl)              // MF-shuffle: inline metric on the shuffle kernel (no d_G)
            ho_cvfem_apply_launch_shfl<double, P, double, true>(u, y, d_elemDof, nullptr, nEl, 0, d_elemList, count, d_corners);
        else if (wantRecompute)
            ho_cvfem_apply_recompute_launch<double, P>(u, y, d_elemDof, d_corners, nEl, 0, d_elemList, count);
        else if (wantNek)
            ho_cvfem_apply_launch_nek<double, P, double>(u, y, d_elemDof, d_G, nEl, 0, d_elemList, count);
        else if (wantShfl)
            ho_cvfem_apply_launch_shfl<double, P, double>(u, y, d_elemDof, d_G, nEl, 0, d_elemList, count);
        else if (wantFp32Metric)
            ho_cvfem_apply_launch<double, P, float >(u, y, d_elemDof, d_G32, nEl, 0, d_elemList, count);
        else if (g_eblock != 0)
            // launch_E silently falls back to the per-P default for out-of-set/infeasible E,
            // so the env-unset and infeasible cases stay bit-identical to the default apply.
            ho_cvfem_apply_launch_E<double, P, double>(u, y, d_elemDof, d_G, nEl, g_eblock, smemCap, 0, d_elemList, count);
        else
            ho_cvfem_apply_launch<double, P, double>(u, y, d_elemDof, d_G,   nEl, 0, d_elemList, count);
    };
    // skipHostGate was computed above (numbering-time) so it could gate the resolved-owner
    // D2H. The host A.1 cross-check + the small-scale host paths below still consult it here.

    // gOwned (global owned-DOF count) feeds the throughput print below -> always compute.
    // On the device path locOwned was already a thrust::count over the resolved device owner
    // (before the ownDev free); here only the host fallback (no ownDev) loops host dofOwner.
    // The owned bitmask is built here for the host fallback too (small-scale only) so the
    // device A.1 readback has a uniform mask -> no per-DOF host loop on either path.
    if (!haveOwnDev) {
        for (long d = 0; d < nDof; ++d) if (dof.dofOwner[d] == rank) ++locOwned;
        std::vector<unsigned> h_mask((nDof + 31) / 32, 0u);
        for (long d = 0; d < nDof; ++d)
            if (dof.dofOwner[d] == rank) h_mask[d >> 5] |= (1u << (d & 31));
        d_ownedMask = h_mask;   // H2D (small-scale fallback only)
    }
    long gOwned = 0;
    MPI_Allreduce(&locOwned, &gOwned, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    const unsigned* d_mask = thrust::raw_pointer_cast(d_ownedMask.data());

    // Host-halo A.1: a host-side cross-check of the device path. SLOW at scale (host halo
    // + 2x nDof H2D/D2H + host reductions over ~600M DOF/rank -> minutes), so skip it on
    // big runs via MARS_HO_SKIP_HOST_GATE. The device A.1 below is the real gate + numbers.
    if (!skipHostGate) {
        if (rank == 0) { printf("[gate] metric ready, entering host A.1\n"); fflush(stdout); }
        std::vector<double> u(nDof, 0.0);
        for (long d = 0; d < nDof; ++d) if (dof.dofOwner[d] == rank) u[d] = 1.0;
        halo.forward(u);
        cudaMemcpy(d_u, u.data(), sizeof(double) * nDof, cudaMemcpyHostToDevice);
        cudaMemset(d_y, 0, sizeof(double) * nDof);
        applyOp(d_u, d_y);
        cudaError_t err = cudaDeviceSynchronize();
        std::vector<double> y(nDof);
        cudaMemcpy(y.data(), d_y, sizeof(double) * nDof, cudaMemcpyDeviceToHost);
        halo.reverseAdd(y);
        if (rank == 0) { printf("[gate] host apply+halo done, entering Allreduce (waits on slowest rank)\n"); fflush(stdout); }
        double locMax = 0;
        for (long d = 0; d < nDof; ++d)
            if (dof.dofOwner[d] == rank) locMax = std::max(locMax, std::abs(y[d]));
        double gMax = 0;
        MPI_Allreduce(&locMax, &gMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        const char* tag = (mode==Numbering::Host) ? "host-numbering"
                        : (mode==Numbering::Gpu)  ? "gpu-numbering" : "gpu-numbering(self-check)";
        if (rank == 0)
            printf("p=%d  owned DOF=%ld  max|A.1| over owned (host halo, %s) = %.3e   [%s]%s\n",
                   P, gOwned, tag, gMax, gMax < 1e-8 ? "PASS" : "FAIL",
                   err != cudaSuccess ? "  (CUDA error!)" : "");
    } else if (rank == 0) {
        printf("[gate] metric ready; host A.1 SKIPPED (MARS_HO_SKIP_HOST_GATE) -> device A.1 is the gate\n"); fflush(stdout);
    }

    // ---- device-halo matvec: same A.1, but forward/reverseAdd on the GPU (scaling path) ----
    halo.uploadDevice();

    // OVERLAP (MARS_HO_OVERLAP=1): hide the forward halo behind the interior apply. Classify
    // each owned element INTERIOR (no ghost DOF) or BOUNDARY (>=1 ghost DOF), then the matvec
    // runs: forwardStart (non-blocking) -> apply INTERIOR -> forwardFinish (wait) -> apply
    // BOUNDARY -> reverseAdd. The split is EXACT: an element is boundary iff any of its DOFs
    // is in the halo recv list (the ghosts the forward fills). Default OFF -> the blocking
    // path stays the reference. reverseAdd stays blocking for now (note: it is also
    // overlappable -- split it the same way and apply boundary elems whose contributions go
    // only to OWNED DOF first; deferred to keep this change focused on the forward).
    const bool wantOverlap = std::getenv("MARS_HO_OVERLAP") != nullptr;
    // MARS_HO_STAGE_TIMING: per-stage cudaEvent breakdown of the overlapped matvec
    // (interior-GPU vs exposed forward-halo wait vs reverseAdd) so we can see WHERE the
    // full-vs-apply gap is -- the forward overlap, or the un-overlapped reverseAdd round-trip.
    const bool wantStageTiming = std::getenv("MARS_HO_STAGE_TIMING") != nullptr;
    thrust::device_vector<int> d_interiorElems, d_boundaryElems;
    int* d_interior = nullptr; int* d_boundary = nullptr;
    size_t nInterior = 0, nBoundary = 0;
    if (wantOverlap) {
        // ghost flag over all DOF: 1 where a DOF is a halo ghost (in recvDof_).
        thrust::device_vector<uint8_t> d_ghostFlag(nDof, (uint8_t)0);
        if (halo.nRecv_ > 0) {
            const int blk = 256, g = (halo.nRecv_ + blk - 1) / blk;
            ho_mark_ghost_dofs_kernel<<<g, blk>>>(halo.d_recvDof_, halo.nRecv_,
                                                  thrust::raw_pointer_cast(d_ghostFlag.data()));
        }
        // per-element boundary flag, then compact element ids into two lists.
        thrust::device_vector<uint8_t> d_isBoundary(nEl > 0 ? nEl : 1, (uint8_t)0);
        if (nEl > 0) {
            const int blk = 256; const size_t g = (nEl + blk - 1) / blk;
            ho_classify_boundary_elems_kernel<P><<<(unsigned)g, blk>>>(
                d_elemDof, thrust::raw_pointer_cast(d_ghostFlag.data()), nEl,
                thrust::raw_pointer_cast(d_isBoundary.data()));
        }
        cudaDeviceSynchronize();
        nBoundary = (size_t)thrust::count(d_isBoundary.begin(), d_isBoundary.begin() + nEl, (uint8_t)1);
        nInterior = nEl - nBoundary;
        d_interiorElems.resize(nInterior > 0 ? nInterior : 1);
        d_boundaryElems.resize(nBoundary > 0 ? nBoundary : 1);
        auto first = thrust::counting_iterator<int>(0);
        auto last  = thrust::counting_iterator<int>((int)nEl);
        const uint8_t* bflag = thrust::raw_pointer_cast(d_isBoundary.data());
        thrust::copy_if(thrust::device, first, last, d_interiorElems.begin(),
                        [bflag] __device__ (int e) { return bflag[e] == 0; });
        thrust::copy_if(thrust::device, first, last, d_boundaryElems.begin(),
                        [bflag] __device__ (int e) { return bflag[e] != 0; });
        cudaDeviceSynchronize();
        d_interior = thrust::raw_pointer_cast(d_interiorElems.data());
        d_boundary = thrust::raw_pointer_cast(d_boundaryElems.data());
        long locSplit[2] = {(long)nInterior, (long)nBoundary}, gSplit[2] = {0,0};
        MPI_Allreduce(locSplit, gSplit, 2, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
        if (rank == 0) {
            printf("[matvec] OVERLAP enabled: interior=%ld boundary=%ld elems (global; boundary %.1f%%)\n",
                   gSplit[0], gSplit[1], 100.0 * (double)gSplit[1] / (double)(gSplit[0] + gSplit[1]));
            fflush(stdout);
        }
    }

    // One full distributed matvec. OVERLAP path hides forward behind the interior apply;
    // the blocking path is the original forward -> apply-all -> reverseAdd. d_y zeroed by
    // the caller. reverseAdd stays blocking on both.
    auto applyMatvec = [&](double* u, double* y) {
        if (wantOverlap) {
            // Overlap BOTH halos: split the interior (owned-only) elements in two and hide the
            // forward behind the first half, the reverse behind the second. Forward brings ghost
            // u (boundary needs it); reverse sends ghost y contributions (only boundary writes
            // them) -- so the order is fwdStart -> interiorA -> fwdFinish -> boundary -> revStart
            // -> interiorB -> revFinish. The trillion run showed the reverse, not the forward,
            // is the exposed comm, so hiding it is the win.
            const size_t nIntA = nInterior / 2, nIntB = nInterior - nIntA;
            halo.forwardDeviceStart(u);
            if (nIntA) applyOp(u, y, d_interior, nIntA);              // hide forward
            halo.forwardDeviceFinish(u);
            if (nBoundary) applyOp(u, y, d_boundary, nBoundary);      // needs ghost u
            halo.reverseAddDeviceStart(y);
            if (nIntB) applyOp(u, y, d_interior + nIntA, nIntB);      // hide reverse
            halo.reverseAddDeviceFinish(y);
        } else {
            halo.forwardDevice(u);
            applyOp(u, y);
            halo.reverseAddDevice(y);
        }
    };

    // u-init ON DEVICE: fill d_u = 1.0 over ALL numDof (owned + ghost), no host loop, no H2D.
    // This is BIT-IDENTICAL to the old 'owned=1, ghost=0, then forwardDevice': forwardDevice
    // overwrites every ghost with its owner's value (here 1.0), so after the forward both
    // schemes give u==1.0 on every DOF going into the apply. Filling all=1 just skips the
    // (now redundant) ghost-zeroing -- the forward is the single source of ghost values.
    thrust::fill(thrust::device_pointer_cast(d_u),
                 thrust::device_pointer_cast(d_u) + nDof, 1.0);
    // Single instrumented matvec. On the blocking path the three stages are timed
    // separately (forward / apply / reverseAdd); on the overlap path the forward time is
    // hidden inside the interior apply, so only the combined forward+apply stage + reverseAdd
    // are meaningful -- print the merged interior+halo+boundary as one [matvec] line there.
    cudaMemset(d_y, 0, sizeof(double) * nDof);
    if (!wantOverlap) {
        MPI_Barrier(MPI_COMM_WORLD); double t_fwd0 = MPI_Wtime();
        halo.forwardDevice(d_u);
        cudaDeviceSynchronize();
        if (rank==0) { printf("[matvec] forward halo    %.2fs\n", MPI_Wtime()-t_fwd0); fflush(stdout); }
        double t_ap0 = MPI_Wtime();
        applyOp(d_u, d_y);
        cudaDeviceSynchronize();
        if (rank==0) { printf("[matvec] apply           %.2fs\n", MPI_Wtime()-t_ap0); fflush(stdout); }
        double t_rev0 = MPI_Wtime();
        halo.reverseAddDevice(d_y);
        cudaDeviceSynchronize();
        if (rank==0) { printf("[matvec] reverseAdd halo %.2fs\n", MPI_Wtime()-t_rev0); fflush(stdout); }
    } else {
        // (c)-schedule overlap: split the interior in two -- interiorA hides the forward halo,
        // interiorB hides the reverse halo. Events bracket each interior half; the host wall-times
        // each Finish (its MPI_Waitall). exposed = max(0, wait - interiorHalf): 0 => fully hidden.
        const size_t nIntA = nInterior / 2, nIntB = nInterior - nIntA;
        cudaEvent_t ea0, ea1, eb0, eb1;
        cudaEventCreate(&ea0); cudaEventCreate(&ea1); cudaEventCreate(&eb0); cudaEventCreate(&eb1);
        MPI_Barrier(MPI_COMM_WORLD); double t_ov0 = MPI_Wtime();
        halo.forwardDeviceStart(d_u);
        cudaEventRecord(ea0); if (nIntA) applyOp(d_u, d_y, d_interior, nIntA); cudaEventRecord(ea1);
        double tf0 = MPI_Wtime(); halo.forwardDeviceFinish(d_u); double twaitFwd = MPI_Wtime() - tf0;
        if (nBoundary) applyOp(d_u, d_y, d_boundary, nBoundary);
        halo.reverseAddDeviceStart(d_y);
        cudaEventRecord(eb0); if (nIntB) applyOp(d_u, d_y, d_interior + nIntA, nIntB); cudaEventRecord(eb1);
        double tr0 = MPI_Wtime(); halo.reverseAddDeviceFinish(d_y); double twaitRev = MPI_Wtime() - tr0;
        cudaDeviceSynchronize();
        if (rank == 0) {
            printf("[matvec] overlapped fwd+rev %.3fms\n", (MPI_Wtime()-t_ov0)*1e3);
            if (wantStageTiming) {
                float msIa = 0, msIb = 0;
                cudaEventElapsedTime(&msIa, ea0, ea1);
                cudaEventElapsedTime(&msIb, eb0, eb1);
                double expFwd = (twaitFwd*1e3 > msIa) ? (twaitFwd*1e3 - msIa) : 0.0;
                double expRev = (twaitRev*1e3 > msIb) ? (twaitRev*1e3 - msIb) : 0.0;
                // "hidden" if the exposed residual is within 2% of the interior half it overlaps:
                // the host Waitall can sit ~interiorB long while the MPI itself finishes early, so
                // a tiny residual still means the halo is effectively buried.
                const char* hidFwd = (expFwd < 0.02 * msIa) ? "YES" : "NO";
                const char* hidRev = (expRev < 0.02 * msIb) ? "YES" : "NO";
                printf("[stage] interiorA=%.3fms fwdWait=%.3fms (exposed %.3f hidden=%s) | interiorB=%.3fms revWait=%.3fms (exposed %.3f hidden=%s)\n",
                       msIa, twaitFwd*1e3, expFwd, hidFwd,
                       msIb, twaitRev*1e3, expRev, hidRev);
            }
            fflush(stdout);
        }
        cudaEventDestroy(ea0); cudaEventDestroy(ea1); cudaEventDestroy(eb0); cudaEventDestroy(eb1);
    }
    // readback ON DEVICE: max over d of (owned ? |d_y[d]| : 0), selected by the packed owned
    // mask. Same value the host loop produced, but only the scalar crosses to host (no nDof
    // D2H). transform_reduce over [0,nDof) with OwnedAbsMax.
    double locMaxD = thrust::transform_reduce(
        thrust::device,
        thrust::counting_iterator<long>(0), thrust::counting_iterator<long>(nDof),
        OwnedAbsMax{d_y, d_mask}, 0.0, thrust::maximum<double>());
    double gMaxD = 0; MPI_Allreduce(&locMaxD, &gMaxD, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    // timing: full device matvec (forward+apply+reverseAdd) vs apply-only -> MDOF/s + comm%
    // The "full" matvec uses applyMatvec: with MARS_HO_OVERLAP it is the overlapped path
    // (forward hidden behind the interior apply), so the full-vs-apply gap (comm%) shrinks
    // toward 0 -- that gap IS the orange-toward-blue weak-scaling improvement.
    const int warm = 5, iters = g_iters;
    for (int it = 0; it < warm; ++it) {
        cudaMemset(d_y, 0, sizeof(double)*nDof);
        applyMatvec(d_u, d_y);
    }
    cudaDeviceSynchronize(); MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();
    for (int it = 0; it < iters; ++it) {
        cudaMemset(d_y, 0, sizeof(double)*nDof);
        applyMatvec(d_u, d_y);
    }
    cudaDeviceSynchronize(); double tFull = (MPI_Wtime() - t0) / iters;
    MPI_Barrier(MPI_COMM_WORLD); double t1 = MPI_Wtime();
    for (int it = 0; it < iters; ++it) { cudaMemset(d_y,0,sizeof(double)*nDof); applyOp(d_u,d_y); }
    cudaDeviceSynchronize(); double tApply = (MPI_Wtime() - t1) / iters;

    // per-stage attribution (MARS_HO_STAGE_TIMING): time interior-only apply (the overlap split)
    // and reverseAdd-only, averaged, so the full-vs-apply gap decomposes into forward vs reverse.
    double tInterior = 0, tRev = 0;
    if (wantStageTiming) {
        if (wantOverlap && nInterior) {
            MPI_Barrier(MPI_COMM_WORLD); double t2 = MPI_Wtime();
            for (int it = 0; it < iters; ++it) { cudaMemset(d_y,0,sizeof(double)*nDof); applyOp(d_u,d_y,d_interior,nInterior); }
            cudaDeviceSynchronize(); tInterior = (MPI_Wtime() - t2) / iters;
        }
        MPI_Barrier(MPI_COMM_WORLD); double t3 = MPI_Wtime();
        for (int it = 0; it < iters; ++it) { halo.reverseAddDevice(d_y); }
        cudaDeviceSynchronize(); tRev = (MPI_Wtime() - t3) / iters;
    }

    double maxFull = 0, maxApply = 0;
    MPI_Allreduce(&tFull,  &maxFull,  1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&tApply, &maxApply, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    double mdofsFull  = (double)gOwned / maxFull  / 1e6;
    double mdofsApply = (double)gOwned / maxApply / 1e6;
    double commFrac   = (maxFull - maxApply) / maxFull * 100.0;

    if (wantStageTiming) {
        double maxInterior = 0, maxRev = 0;
        MPI_Allreduce(&tInterior, &maxInterior, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&tRev,      &maxRev,      1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (rank == 0) {
            // fwdGap = full - applyAll - reverseAdd: the forward-halo residual, isolated from the
            // (unoverlapped) reverseAdd cost. If reverseAdd dominates, the gap is the SECOND halo
            // round-trip -- not the forward overlap, which the per-stage [stage] line confirms.
            double fwdGap = (maxFull - maxApply - maxRev) * 1e3;
            printf("[bd] full=%.3fms applyAll=%.3fms interior=%.3fms reverseAdd=%.3fms fwdGap=%.3fms commFrac=%.1f%%\n",
                   maxFull*1e3, maxApply*1e3, maxInterior*1e3, maxRev*1e3, fwdGap, commFrac);
            fflush(stdout);
        }
    }

    // FP64 work intensity of the store-d_G apply, in FLOP per globally-unique DOF.
    // Verified against the kernel loop counts: per element the 4-step sweep does
    // 3*P*(P+1)^2*(12*(P+1)+7) flops; a structured tensor-product hex contributes
    // P^3 unique DOF/element, so dividing by P^3 gives the per-DOF intensity that
    // pairs with gOwned (~ncells^3 * P^3). GFLOP/s = FLOP/DOF * apply-only GDOF/s.
    const double Pd = (double)P, Nd = (double)(P + 1);
    const double flopPerDof = 3.0 * Pd * Nd * Nd * (12.0 * Nd + 7.0) / (Pd * Pd * Pd);
    const double gdofsApply = mdofsApply / 1e3;            // MDOF/s -> GDOF/s
    const double gflops     = flopPerDof * gdofsApply;     // FP64 GFLOP/s

    if (rank == 0) {
        // A.1 (A*const=0) is fp32-SAFE by construction: the metric multiplies grad(const)=0,
        // so the stored-precision choice cannot move it -> the gate stays at <1e-8 for both
        // fp64 and fp32 d_G. (There is no host-vs-device operator-parity gate in this driver;
        // it runs only A.1 gates, so no tolerance loosening is needed for the fp32 path.)
        const char* pathTag = wantRecompute ? "RECOMPUTE (no d_G)"
                            : (wantFp32Metric ? "STORE-d_G(fp32)" : "STORE-d_G");
        printf("p=%d  device A.1 = %.3e [%s]   full matvec %.0f MDOF/s | apply-only %.0f MDOF/s | comm %.1f%%  [path=%s]\n",
               P, gMaxD, gMaxD < 1e-8 ? "PASS" : "FAIL", mdofsFull, mdofsApply, commFrac, pathTag);
        printf("       THROUGHPUT(p=%d): apply-only %.3f GDOF/s | %.1f FLOP/DOF | %.1f GFLOP/s (FP64) | %.3g DOF/GPU  [path=%s]\n",
               P, gdofsApply, flopPerDof, gflops, (double)gOwned / (double)numRanks, pathTag);

        // Record the row for the end-of-sweep table (printed once by main).
        SweepRow row;
        row.p = P; row.gOwnedDof = gOwned; row.dofPerGpu = (double)gOwned / (double)numRanks;
        row.gdofsApply = gdofsApply; row.flopPerDof = flopPerDof; row.gflops = gflops;
        row.pass = (gMaxD < 1e-8);
        g_sweepRows.push_back(row);
        // measured wall-clock (slowest rank, averaged over `iters` real matvecs) -> the headline number
        printf("       measured: %.2f ms/matvec (slowest of %d ranks, %d iters) | sustained %.3f TDOF/s | %.3g global DOF  [path=%s]\n",
               maxFull*1e3, numRanks, iters, mdofsFull/1e6, (double)gOwned, pathTag);

        // ---- STORAGE report (the recompute payoff). Per-DOF = total operator state /
        // this rank's owned DOF. d_corners (24 doubles/elem) is the recompute state; d_G
        // (3*P*(P+1)^2 vec3 doubles/elem) is the store-d_G state. nEl/nDofRank are this
        // rank's local counts -- the ratio is intrinsic (doesn't depend on the rank split).
        const double nDofRank = (double)nDof;
        const double dgBytes   = gElemBytes * (double)gLen;               // d_G array bytes (fp32 or fp64)
        const double cornBytes = sizeof(double) * (double)nEl * 24.0;     // d_corners bytes
        const double dgPerDof   = dgBytes   / nDofRank;
        const double cornPerDof = cornBytes / nDofRank;
        // Assembled 7-nnz graph CVFEM reference (p=1): ~84 B/DOF (7 doubles vals + 7 ints
        // cols + a row ptr amortized). Stated as the assembled baseline for the comparison;
        // only meaningful at p=1 (the 7-nnz stencil is the low-order graph operator).
        const double assembled7nnzPerDof = 84.0;
        printf("       STORAGE (p=%d):  store-d_G = %.1f B/DOF  |  recompute(corners) = %.1f B/DOF"
               "  |  ratio store/recompute = %.1fx  |  assembled 7-nnz ref = %.0f B/DOF\n",
               P, dgPerDof, cornPerDof, cornPerDof > 0 ? dgPerDof/cornPerDof : 0.0, assembled7nnzPerDof);

        // ---- per-GPU DOF CEILING. Operator state bytes/DOF on each path (NOT counting the
        // u/y vectors, which both paths pay): store-d_G pays d_G + d_corners + elemDof;
        // recompute pays d_corners + elemDof only (no d_G). Plus the apply's own u+y = 16
        // B/DOF (2 doubles), common to both. We report the implied max DOF/GPU as
        // hbm / total-bytes-per-DOF. HBM is env-overridable (MARS_HO_HBM_GB, default 80).
        const double elemDofPerDof = sizeof(int) * (double)nEl * (double)N3 / nDofRank;
        const double uyPerDof      = sizeof(double) * 2.0;                // d_u + d_y
        double hbmGB = 80.0;
        if (const char* e = std::getenv("MARS_HO_HBM_GB")) hbmGB = std::atof(e);
        const double storeTotPerDof = dgPerDof + cornPerDof + elemDofPerDof + uyPerDof;
        const double recompTotPerDof =          cornPerDof + elemDofPerDof + uyPerDof;
        const double storeCeil  = hbmGB * 1e9 / storeTotPerDof;
        const double recompCeil = hbmGB * 1e9 / recompTotPerDof;
        printf("       CEILING (HBM=%.0f GB):  store-d_G %.1f B/DOF -> ~%.2g DOF/GPU  |  recompute %.1f B/DOF -> ~%.2g DOF/GPU"
               "  |  recompute lifts ceiling ~%.1fx (frees d_G)\n",
               hbmGB, storeTotPerDof, storeCeil, recompTotPerDof, recompCeil,
               recompTotPerDof > 0 ? storeTotPerDof/recompTotPerDof : 0.0);
        printf("       NOTE: the p=1 throughput WIN of recompute vs store-d_G and vs the assembled 7-nnz operator is TO BE MEASURED\n"
               "             on Alps (run store-d_G and --recompute at the same p, compare ms/matvec) -- not asserted here.\n");
    }

    // ---- E-SWEEP (--esweep / MARS_HO_ESWEEP): for THIS order p, time the store-d_G apply
    // across the FEASIBLE elements-per-block set {2,4,8,16,32} so the occupancy sweet spot is
    // readable from one job. Store-d_G / Laplacian only (the apply the ncu occupancy study
    // profiled). d_u is already 1.0 on every DOF (the A.1 input), so each E reuses it: zero
    // d_y -> apply-only(E) -> reverseAdd -> owned-max = the per-E A.1 (must PASS). Skips any
    // (p,E) over the smem carveout / 1024-thread cap and records nothing for it (never
    // launched). The default-E timing above is untouched; this is purely additive.
    if (std::getenv("MARS_HO_ESWEEP") != nullptr) {
        if (!eblockOverridable) {
            if (rank == 0) printf("[esweep] SKIPPED: E-sweep is store-d_G fp64 only (drop --recompute / MARS_HO_FP32_METRIC)\n");
        } else {
            for (int eb : kEblockCandidates) {
                if (!ho_cvfem_eblock_feasible<P>(eb, smemCap)) {
                    if (rank == 0)
                        printf("[esweep] p=%d E=%d INFEASIBLE (smem %.1f KB > cap %.0f KB, or block %d>1024) -- skipped\n",
                               P, eb, ho_cvfem_smem_bytes<P>(eb)/1024.0, smemCap/1024.0,
                               eb * /*tpe*/ (int)(HoCvfemLaunchDefault<P>::Block / HoCvfemLaunchDefault<P>::Elems));
                    continue;
                }
                // A.1 for this E: apply-only on d_u(=1) + reverseAdd, then owned max.
                cudaMemset(d_y, 0, sizeof(double) * nDof);
                ho_cvfem_apply_launch_E<double, P, double>(d_u, d_y, d_elemDof, d_G, nEl, eb, smemCap);
                halo.reverseAddDevice(d_y);
                cudaDeviceSynchronize();
                double locA1 = thrust::transform_reduce(
                    thrust::device,
                    thrust::counting_iterator<long>(0), thrust::counting_iterator<long>(nDof),
                    OwnedAbsMax{d_y, d_mask}, 0.0, thrust::maximum<double>());
                double gA1 = 0; MPI_Allreduce(&locA1, &gA1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

                // apply-only throughput at this E (same warm/iters protocol as the main timing).
                for (int it = 0; it < 5; ++it) {
                    cudaMemset(d_y, 0, sizeof(double) * nDof);
                    ho_cvfem_apply_launch_E<double, P, double>(d_u, d_y, d_elemDof, d_G, nEl, eb, smemCap);
                }
                cudaDeviceSynchronize(); MPI_Barrier(MPI_COMM_WORLD);
                double te0 = MPI_Wtime();
                for (int it = 0; it < g_iters; ++it) {
                    cudaMemset(d_y, 0, sizeof(double) * nDof);
                    ho_cvfem_apply_launch_E<double, P, double>(d_u, d_y, d_elemDof, d_G, nEl, eb, smemCap);
                }
                cudaDeviceSynchronize();
                double tE = (MPI_Wtime() - te0) / g_iters, maxTE = 0;
                MPI_Allreduce(&tE, &maxTE, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                double gdofsE = (double)gOwned / maxTE / 1e9;

                if (rank == 0) {
                    ESweepRow r;
                    r.p = P; r.eblock = eb; r.dofPerGpu = (double)gOwned / (double)numRanks;
                    r.gdofsApply = gdofsE; r.smemKB = ho_cvfem_smem_bytes<P>(eb) / 1024.0;
                    r.a1 = gA1; r.pass = (gA1 < 1e-8);
                    g_esweepRows.push_back(r);
                    printf("[esweep] p=%d E=%2d  %.4g DOF/GPU  %.3f GDOF/s  %.1f KB/block  A.1=%.3e [%s]\n",
                           P, eb, r.dofPerGpu, gdofsE, r.smemKB, gA1, r.pass ? "PASS" : "FAIL");
                    fflush(stdout);
                }
            }
        }
    }

    halo.freeDevice();
    // d_elemDof is freed only if it was cudaMalloc'd here (host path). On the device path
    // it aliases d_elemDofDev, which frees itself when this scope ends.
    if (!elemDofResident) cudaFree(d_elemDof);
    cudaFree(d_corners); cudaFree(d_G); cudaFree(d_G32); cudaFree(d_u); cudaFree(d_y);
}

// Build the cube mesh + cstone domain + distributed DOF for ONE order, then run the
// timed apply. Factored out of main so the throughput sweep can call it once per p with
// a per-p ncells (each order needs a fresh mesh: holding DOF/GPU roughly constant means
// ncells must FALL as p rises, since DOF/elem grows as p^3). Single-order runs call it
// once. Everything (mesh-gen through the cudaFree in runDistApply) is rebuilt per call,
// so there is no cross-order state leak.
static void runOrder(int P, size_t ncells, int rank, int numRanks,
                     Numbering mode, bool cliIrregular)
{
    MPI_Barrier(MPI_COMM_WORLD);
    double tgen0 = MPI_Wtime();
    auto [gn, ge, gx, gy, gz, lconn] =
        generateCubeElementPartition<RealType, KeyType>(ncells, rank, numRanks);
    (void)gn; (void)ge;
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) { printf("[build] meshGen %.1fs\n", MPI_Wtime()-tgen0); fflush(stdout); }
    // --irregular: perturb interior node coords -> distorted hexes + an irregular SFC
    // partition. This is the "genuine unstructured" stress for the halo's geometric 2h
    // search and the ownership's corner multiplicities, vs the clean cube. The offset is
    // a deterministic hash of the integer grid index, so a node shared across ranks gets
    // the SAME offset everywhere (no mesh tearing). Global-boundary nodes stay fixed so
    // the cube stays a cube; amplitude 0.25h keeps every hex valid (no inverted detJ).
    if (cliIrregular) {
        const double h = 1.0 / (double)ncells, amp = 0.25 * h;
        auto idx = [&](RealType c){ return (long)((double)c / h + 0.5); };   // coords in [0,1]
        auto noise = [](long ix, long iy, long iz, int ax) {
            uint64_t k = (uint64_t)(ix*73856093L) ^ (uint64_t)(iy*19349663L)
                       ^ (uint64_t)(iz*83492791L) ^ (uint64_t)((unsigned)ax*2654435761u);
            k = (k ^ (k>>33)) * 0xff51afd7ed558ccdULL; k ^= k>>33;
            return ((double)(k & 0xffffffULL) / (double)0xffffffULL) * 2.0 - 1.0; // [-1,1]
        };
        for (size_t n = 0; n < gx.size(); ++n) {
            long ix = idx(gx[n]), iy = idx(gy[n]), iz = idx(gz[n]);
            if (ix<=0||ix>=(long)ncells || iy<=0||iy>=(long)ncells || iz<=0||iz>=(long)ncells) continue;
            gx[n] += (RealType)(amp * noise(ix,iy,iz,0));
            gy[n] += (RealType)(amp * noise(ix,iy,iz,1));
            gz[n] += (RealType)(amp * noise(ix,iy,iz,2));
        }
        if (rank == 0) { printf("[build] --irregular: jittered interior nodes <=%.3g (0.25h) -> distorted hexes / irregular partition\n", amp); fflush(stdout); }
    }
    typename Domain::HostCoordsTuple h_coords{std::move(gx), std::move(gy), std::move(gz)};
    typename Domain::HostConnectivityTuple h_conn{
        std::move(lconn[0]), std::move(lconn[1]), std::move(lconn[2]), std::move(lconn[3]),
        std::move(lconn[4]), std::move(lconn[5]), std::move(lconn[6]), std::move(lconn[7])};
    // Global decomposition bucketSize. cstone's global octree is replicated on
    // every rank and its per-leaf count Allreduce (MPI_UNSIGNED) must stay under
    // 2 GiB, or Cray MPICH's chunked recursive-doubling collective overflows its
    // 32-bit byte count and truncates. Leaves ~ 3 * ncells^3 / bucket (octree
    // over-refine ~3x), each leaf 4 bytes -> raise the bucket for huge meshes so
    // leaves*4 stays well under 2^31. bucketSizeFocus stays 8 (local resolution
    // unchanged); only the decomposition tree coarsens. Override via env.
    int gbucket = 64;
    if (const char* e = std::getenv("MARS_GLOBAL_BUCKETSIZE")) {
        gbucket = std::atoi(e);
    } else {
        const double need = 3.0 * (double)ncells * (double)ncells * (double)ncells / 4.0e8;
        while (gbucket < need && gbucket < 8192) gbucket *= 2;
    }
    if (rank == 0) { printf("[build] global bucketSize=%d (msg-safe replicated tree)\n", gbucket); fflush(stdout); }
    MPI_Barrier(MPI_COMM_WORLD);
    double tdom0 = MPI_Wtime();
    Domain domain(h_coords, h_conn, rank, numRanks, gbucket, false, 8u);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) { printf("[build] cstone domain.sync %.1fs\n", MPI_Wtime()-tdom0); fflush(stdout); }

    // Device-input setup = the GPU-native default: GPU numbering, >1 rank, and neither the
    // host fallback (MARS_HO_HOST) nor the zero-MPI local-ownership resolve (which both
    // consume the HOST corner tables). The host tables are skipped entirely in that case.
    const bool wantDeviceInputs = (mode == Numbering::Gpu) && (numRanks > 1)
                                  && std::getenv("MARS_HO_HOST") == nullptr
                                  && std::getenv("MARS_HO_LOCAL_OWNERSHIP") == nullptr;
    double te0 = MPI_Wtime();
    DistDof D = extractDistDof(domain, rank, numRanks, wantDeviceInputs);
    double te1 = MPI_Wtime();
    if (rank == 0) { printf("[build] extractDistDof %.1fs (nodeCount=%zu nOwnedElem=%zu)\n",
                            te1 - te0, D.nodeCount, D.nOwnedElem); fflush(stdout); }

    if (rank == 0)
        printf("\n== HO distributed matrix-free apply gate (A.1 = HO-Laplacian.const = 0) ==  ncells=%zu ranks=%d  mode=%s\n",
               (size_t)ncells, numRanks,
               mode==Numbering::Host ? "host" : mode==Numbering::Gpu ? "gpu" : "self-check");
    // p=1..8: the matfree kernel infra already supports p=8 (HO_CVFEM_MAX_P=8, the
    // __constant__ tables and HoCvfemLaunchDefault<8> are all sized for it; p=8 smem is
    // ~12.7 KB/block, far under the Hopper carveout). Only the dispatch was capped at 4.
    switch (P) {
        case 1: runDistApply<1>(D, rank, numRanks, mode); break;
        case 2: runDistApply<2>(D, rank, numRanks, mode); break;
        case 3: runDistApply<3>(D, rank, numRanks, mode); break;
        case 4: runDistApply<4>(D, rank, numRanks, mode); break;
        case 5: runDistApply<5>(D, rank, numRanks, mode); break;
        case 6: runDistApply<6>(D, rank, numRanks, mode); break;
        case 7: runDistApply<7>(D, rank, numRanks, mode); break;
        case 8: runDistApply<8>(D, rank, numRanks, mode); break;
        default: if (rank == 0) printf("unsupported p=%d (build adds 1..8)\n", P);
    }
}

// Hold DOF/GPU roughly constant across orders so the throughput-vs-p comparison is fair.
// Global owned DOF ~ (ncells^3) * p^3 (each hex contributes p^3 unique DOF), so the per-p
// element count must scale as 1/p^3 to keep ncells^3 * p^3 fixed. Given a p=1 baseline
// ncells n1, the per-p edge count is n_p = round(n1 / p), which keeps n_p^3 * p^3 = n1^3
// to within rounding -> DOF/GPU within a few percent across the sweep. We clamp n_p so it
// stays a positive multiple-of-nothing integer >= 1.
static size_t ncellsForOrder(size_t n1, int p)
{
    if (p <= 1) return n1;
    size_t n = (size_t)llround((double)n1 / (double)p);
    return n < 1 ? 1 : n;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, numRanks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
    int devCount = 0; cudaGetDeviceCount(&devCount); if (devCount > 0) cudaSetDevice(rank % devCount);

    size_t ncells = 16; int P = 2;
    bool cliGpu = false, cliHost = false, cliSelfCheck = false, cliIrregular = false, cliSweep = false;
    for (int i = 1; i < argc; ++i) { std::string a = argv[i];
        if (a.rfind("--ncells=",0)==0) ncells = std::stoull(a.substr(9));
        else if (a.rfind("--p=",0)==0) P = std::stoi(a.substr(4));
        else if (a.rfind("--iters=",0)==0) g_iters = std::max(1, std::stoi(a.substr(8)));
        else if (a == "--gpu-numbering")  cliGpu = true;   // accepted; GPU is the default now
        else if (a == "--host-numbering") cliHost = true;
        else if (a == "--self-check")     cliSelfCheck = true;
        else if (a == "--irregular")      cliIrregular = true;
        // --sweep: run p=1..8 in ONE job, holding DOF/GPU ~constant (ncells per p, see
        // ncellsForOrder), and print a single throughput table at the end. --ncells sets
        // the p=1 baseline edge count; each higher order uses round(ncells/p).
        else if (a == "--sweep")          cliSweep = true;
        // --esweep: for each order p (the same p-loop as --sweep), additionally time the
        // store-d_G apply across the FEASIBLE elements-per-block set {2,4,8,16,32} and print
        // a throughput-vs-E table -- the occupancy sweet spot per order in one job. Sets the
        // MARS_HO_ESWEEP env so runDistApply (a template, no flag threaded) sees it uniformly,
        // and implies --sweep so it walks p=1..8.
        else if (a == "--esweep")         { cliSweep = true; setenv("MARS_HO_ESWEEP", "1", 1); }
        // Operator path: --MF = recompute, no d_G (matrix-free); --PA = store-d_G (partial
        // assembly, the default). Set the env so runDistApply (a template, not threaded a
        // flag) sees it uniformly with MARS_HO_RECOMPUTE. setenv overwrites.
        else if (a == "--recompute" || a == "--MF") setenv("MARS_HO_RECOMPUTE", "1", 1);
        else if (a == "--PA")             { /* explicit store-d_G; this is the default */ } }

    // Default = GPU numbering (radix-uint64 dedup): bit-identical to host (A.1 matches),
    // ~10x faster build, lower host-memory peak, and the SAME per-GPU ceiling -- the
    // ceiling is the apply's d_G, and the numbering's device scratch is freed before d_G
    // is allocated. Opt OUT to host (the self-check oracle / CPU-only debugging) with
    // --host-numbering or MARS_HO_HOST_NUMBERING. --self-check runs both and A/Bs them.
    const char* envHost = std::getenv("MARS_HO_HOST_NUMBERING");
    bool wantHost = cliHost || (envHost && std::atoi(envHost) != 0);
    (void)cliGpu;
    Numbering mode = cliSelfCheck ? Numbering::SelfCheck
                   : wantHost      ? Numbering::Host
                                   : Numbering::Gpu;

    if (cliSweep) {
        if (rank == 0) {
            const bool sweepFp32 = std::getenv("MARS_HO_RECOMPUTE") == nullptr
                                 && std::getenv("MARS_HO_FP32_METRIC") != nullptr;
            printf("\n== HO matrix-free LAPLACIAN throughput sweep (store-d_G / --PA, metric %s) ==  "
                   "ranks=%d  p=1..8  baseline ncells=%zu  iters=%d\n",
                   sweepFp32 ? "fp32" : "fp64", numRanks, ncells, g_iters);
            printf("   holding DOF/GPU ~constant: ncells(p) = round(%zu / p)\n", ncells); fflush(stdout);
        }
        for (int p = 1; p <= 8; ++p) {
            size_t np = ncellsForOrder(ncells, p);
            if (rank == 0) { printf("\n--- sweep p=%d  ncells=%zu ---\n", p, np); fflush(stdout); }
            runOrder(p, np, rank, numRanks, mode, cliIrregular);
        }
        if (rank == 0) {
            printf("\n================ HO matrix-free LAPLACIAN throughput sweep (store-d_G, apply-only) ================\n");
            printf(" %2s | %12s | %14s | %12s | %10s | %8s\n",
                   "p", "DOF/GPU", "apply GDOF/s", "GFLOP/s(FP64)", "FLOP/DOF", "A.1");
            printf("----+--------------+----------------+--------------+------------+----------\n");
            for (const auto& r : g_sweepRows)
                printf(" %2d | %12.4g | %14.3f | %12.1f | %10.1f | %8s\n",
                       r.p, r.dofPerGpu, r.gdofsApply, r.gflops, r.flopPerDof, r.pass ? "PASS" : "FAIL");
            printf("===================================================================================================\n");
            fflush(stdout);

            // E-sweep table: throughput vs elements-per-block per order, so the occupancy
            // sweet spot is one read-off. Only printed when --esweep collected rows.
            if (!g_esweepRows.empty()) {
                printf("\n================ HO LAPLACIAN elements-per-block (E) sweep (store-d_G, apply-only) ================\n");
                printf(" %2s | %3s | %12s | %14s | %14s | %10s | %6s\n",
                       "p", "E", "DOF/GPU", "apply GDOF/s", "smem KB/block", "A.1", "gate");
                printf("----+-----+--------------+----------------+----------------+------------+--------\n");
                for (const auto& r : g_esweepRows)
                    printf(" %2d | %3d | %12.4g | %14.3f | %14.1f | %10.3e | %6s\n",
                           r.p, r.eblock, r.dofPerGpu, r.gdofsApply, r.smemKB, r.a1, r.pass ? "PASS" : "FAIL");
                printf("==================================================================================================\n");
                fflush(stdout);
            }
        }
        MPI_Finalize();
        return 0;
    }

    runOrder(P, ncells, rank, numRanks, mode, cliIrregular);
    MPI_Finalize();
    return 0;
}
