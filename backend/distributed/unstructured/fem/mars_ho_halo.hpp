#pragma once

// Multi-rank halo for HIGH-ORDER DOF (corners + shared edge/face GLL DOF).
//
// Same construction as the P1 Option-A node halo (domain.cu buildFromCstoneHalos):
// RECEIVER-DRIVEN. Each rank requests its ghost DOF (dofOwner != me) from the
// owner by the canonical DofKey from HODofHandler::buildDistributed; the owner
// maps the received key -> its local DOF -> send list. So A.send[B] == B.recv[A]
// by construction (no truncation, no ownership over-claim). Interior DOF are
// element-local (never shared) -> never enter the halo.
//
// Peers: reuse the P1 halo peer list (domain.getNodeHaloTopology().peers_) -- a
// high-order DOF is shared only with mesh neighbours, which are exactly the P1
// halo peers.
//
// v1 exchange operates on a HOST DOF vector (for the multi-rank correctness gates:
// forward-a-constant, A*1, A*linear). The scaling path swaps forward/reverseAdd
// for a device gather/scatter + CUDA-aware MPI, identical to exchangeNodeHalo.

#include <mpi.h>
#include <vector>
#include <array>
#include <map>
#include <algorithm>
#include <utility>
#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include "backend/distributed/unstructured/fem/mars_ho_dof_handler.hpp"
#include "backend/distributed/unstructured/fem/mars_hex_face_orientations.hpp"

#ifdef __CUDACC__
#include <cstdint>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/copy.h>
#include <thrust/sort.h>
#include <thrust/gather.h>
#include <thrust/sequence.h>
#include <thrust/scan.h>
#include <thrust/transform.h>
#include <thrust/count.h>
#include <thrust/inner_product.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/execution_policy.h>
#endif

namespace mars {
namespace fem {

// Resolve ownership of shared HO edge/face DOF by MIN-RANK-AMONG-HOLDERS.
// buildDistributed sets these to provisional owner = myRank and flags them in
// dofShared. Each rank broadcasts its shared-DOF keys to all peers; for each local
// shared DOF a peer also holds, owner = min(current, peer rank). The owner is thus
// always a rank that CONTAINS the DOF -> no orphans, and all holders agree (the
// holders of an edge/face are mutual neighbours). Corners (cstone P1 ownership) and
// interiors (element-local) are already correct and not in dofShared -> untouched.
inline void resolveHoDofOwnership(const std::vector<uint8_t>&               dofShared,
                                  const std::vector<HODofHandler::DofKey>&  dofKey,
                                  std::vector<int>&                         dofOwner,
                                  int                                       myRank,
                                  const std::vector<int>&                   peers)
{
    (void)myRank;
    const int np = (int)peers.size();
    auto pack = [](const HODofHandler::DofKey& k) {
        return std::array<long,6>{ (long)k.kind, k.g0, k.g1, k.g2, k.g3, (long)k.pos }; };

    // sorted (key -> local dof), NOT std::map: a single O(n log n) sort + cache-friendly
    // binary search replaces the millions of tree inserts/lookups that dominated this
    // stage at scale. Keys are unique per DOF, so lower_bound is an exact match.
    std::vector<std::pair<std::array<long,6>, int>> myShared;
    std::vector<std::array<long,6>>   myKeys;
    for (int d = 0; d < (int)dofShared.size(); ++d)
        if (dofShared[d]) { auto k = pack(dofKey[d]); myShared.push_back({k, d}); myKeys.push_back(k); }
    std::sort(myShared.begin(), myShared.end(),
              [](const auto& a, const auto& b){ return a.first < b.first; });
    int myCnt = (int)myKeys.size();

    std::vector<int> rc(np, 0);
    {
        std::vector<MPI_Request> rq; rq.reserve(2 * np);
        for (int i=0;i<np;++i){ MPI_Request r; MPI_Irecv(&rc[i],1,MPI_INT,peers[i],0x4850,MPI_COMM_WORLD,&r); rq.push_back(r); }
        for (int i=0;i<np;++i){ MPI_Request r; MPI_Isend(&myCnt,1,MPI_INT,peers[i],0x4850,MPI_COMM_WORLD,&r); rq.push_back(r); }
        if (!rq.empty()) MPI_Waitall((int)rq.size(), rq.data(), MPI_STATUSES_IGNORE);
    }

    std::vector<std::vector<std::array<long,6>>> peerKeys(np);
    {
        std::vector<MPI_Request> rq; rq.reserve(2 * np);
        for (int i=0;i<np;++i){ peerKeys[i].resize(rc[i]); if (rc[i]) { MPI_Request r; MPI_Irecv(peerKeys[i].data(), rc[i]*6, MPI_LONG, peers[i], 0x4851, MPI_COMM_WORLD, &r); rq.push_back(r);} }
        for (int i=0;i<np;++i){ if (myCnt) { MPI_Request r; MPI_Isend(myKeys.data(), myCnt*6, MPI_LONG, peers[i], 0x4851, MPI_COMM_WORLD, &r); rq.push_back(r);} }
        if (!rq.empty()) MPI_Waitall((int)rq.size(), rq.data(), MPI_STATUSES_IGNORE);
    }

    for (int i=0;i<np;++i) {
        int pr = peers[i];
        for (auto& k : peerKeys[i]) {
            auto it = std::lower_bound(myShared.begin(), myShared.end(), k,
                        [](const auto& a, const std::array<long,6>& key){ return a.first < key; });
            if (it != myShared.end() && it->first == k && pr < dofOwner[it->second]) dofOwner[it->second] = pr;
        }
    }
}

// ZERO-MPI resolveHoDofOwnership. Same result as resolveHoDofOwnership above
// (MIN-RANK-AMONG-HOLDERS) but with NO key exchange: it reuses cstone's per-element
// halo, which already places a copy of every rank-boundary element on each rank that
// touches it. Each rank therefore SEES the elements its peers own (as halo elements),
// so it can read those peers' contributions to a shared edge/face LOCALLY instead of
// asking for them over MPI.
//
// HOW A HALO ELEMENT'S EDGE/FACE MAPS TO A LOCAL SHARED DOF:
//   A halo element's 8 corner LOCAL ids live in the SAME local-node id space as owned
//   elements (createElementToNodeLocalIdMap fills all 8 columns over the full element
//   range through the one shared SFC->localId map). So for each of its 12 edges / 6
//   faces we recompute the canonical DofKey (kind, SORTED global corner ids, pos)
//   with the EXACT convention buildDistributed uses -- the global-id-only key is
//   identical on every rank that holds the entity. We binary-search that key in the
//   local sorted shared-key list (the same `myShared` lookup the MPI path runs on
//   peer keys). A hit is a shared DOF this rank also holds; we then take
//   min(dofOwner[hit], elemOwner[haloElem]). A miss means the halo edge/face is not
//   one of my shared DOF -> skip.
//
// WHY THIS IS EXACT: the only ranks that hold a given shared edge/face are the ranks
// owning an element that contains it; cstone delivers exactly those elements (owned +
// halo) to me, and the owner of each is in elemOwner. So min over {my owned elements'
// myRank} U {halo elements' true owners} == min over ALL holder ranks == the same
// value the MPI broadcast computes. Owned edges/faces start at provisional myRank (set
// by buildDistributed); halo elements pull the owner DOWN to the true minimum.
//
// PRECONDITION: elemCorners/elemOwner MUST include the halo elements with their true
// owners (elemOwner[e] = cornerOwner of e's lowest-SFC corner, the cstone placement
// key). extractDistDof fills these under MARS_HO_LOCAL_OWNERSHIP. Owned elements keep
// elemOwner = myRank.
inline void resolveHoDofOwnershipLocal(const std::vector<std::array<int,8>>& elemCorners,
                                       const std::vector<long>&               cornerGid,
                                       const std::vector<int>&                elemOwner,
                                       int                                    order,
                                       const std::vector<uint8_t>&            dofShared,
                                       const std::vector<HODofHandler::DofKey>& dofKey,
                                       std::vector<int>&                      dofOwner,
                                       int                                    myRank)
{
    (void)myRank;
    const int P   = order;
    const int pm1 = P - 1;

    // Same hex topology tables / corner convention as buildDistributed -- the DofKey
    // must be byte-identical, so the edge/face traversal must match exactly.
    static const int EDGES[12][2] = {
        {0,1},{1,2},{2,3},{3,0}, {4,5},{5,6},{6,7},{7,4}, {0,4},{1,5},{2,6},{3,7} };
    static const int FACES[6][4] = {
        {0,1,2,3}, {4,5,6,7}, {0,1,5,4}, {3,2,6,7}, {1,2,6,5}, {0,3,7,4} };

    auto pack = [](const HODofHandler::DofKey& k) {
        return std::array<long,6>{ (long)k.kind, k.g0, k.g1, k.g2, k.g3, (long)k.pos }; };

    // sorted (key -> local dof) over MY shared edge/face DOF -- the lookup target, same
    // structure as the MPI path's myShared. Keys are unique per DOF -> exact match.
    std::vector<std::pair<std::array<long,6>, int>> myShared;
    for (int d = 0; d < (int)dofShared.size(); ++d)
        if (dofShared[d]) myShared.push_back({pack(dofKey[d]), d});
    std::sort(myShared.begin(), myShared.end(),
              [](const auto& a, const auto& b){ return a.first < b.first; });

    // One pass over EVERY element (owned + halo). For each edge/face DOF the element
    // defines, recompute the canonical key and atomicMin its holder's rank onto the
    // matching local shared DOF. Owned elements contribute myRank (already the
    // provisional), halo elements contribute their true owner.
    const long nElem = (long)elemCorners.size();
    for (long e = 0; e < nElem; ++e) {
        const int  owner = elemOwner[e];
        const auto& c    = elemCorners[e];

        // EDGES: key = {1, lowGid, highGid, -1, -1, posCanon}. The (p-1) interior DOF
        // of one physical edge differ only in pos, so we resolve all (p-1) of them.
        for (int ed = 0; ed < 12; ++ed) {
            int cA = c[EDGES[ed][0]], cB = c[EDGES[ed][1]];
            long gA = cornerGid[cA], gB = cornerGid[cB];
            long klo = (gA<=gB)?gA:gB, khi = (gA<=gB)?gB:gA;
            for (int t = 1; t <= pm1; ++t) {
                int posCanon = (gA<=gB) ? (t-1) : (pm1-t);
                std::array<long,6> key{1, klo, khi, -1, -1, (long)posCanon};
                auto it = std::lower_bound(myShared.begin(), myShared.end(), key,
                            [](const auto& a, const std::array<long,6>& kk){ return a.first < kk; });
                if (it != myShared.end() && it->first == key && owner < dofOwner[it->second])
                    dofOwner[it->second] = owner;
            }
        }

        // FACES: key = {2, sorted 4 global corner ids, posCanon}. Sweep the (p-1)^2
        // interior (v1,v2) cells; hex_face_canonical_pos gives the orientation-invariant
        // pos exactly as buildDistributed (same global-id frame on every holder).
        for (int f = 0; f < 6; ++f) {
            const int* fc = FACES[f];
            int lcc[4] = { c[fc[0]], c[fc[1]], c[fc[2]], c[fc[3]] };
            long s[4]  = { cornerGid[lcc[0]], cornerGid[lcc[1]], cornerGid[lcc[2]], cornerGid[lcc[3]] };
            std::sort(s, s+4);
            for (int v1 = 1; v1 <= pm1; ++v1)
            for (int v2 = 1; v2 <= pm1; ++v2) {
                int posCanon = hex_face_canonical_pos(P, v1, v2, lcc[0], lcc[1], lcc[2], lcc[3]);
                std::array<long,6> key{2, s[0], s[1], s[2], s[3], (long)posCanon};
                auto it = std::lower_bound(myShared.begin(), myShared.end(), key,
                            [](const auto& a, const std::array<long,6>& kk){ return a.first < kk; });
                if (it != myShared.end() && it->first == key && owner < dofOwner[it->second])
                    dofOwner[it->second] = owner;
            }
        }
    }
}

#ifdef __CUDACC__

// ALL-GPU resolveHoDofOwnership: device-resident, point-to-point, CUDA-aware MPI.
// Bit-identical result to the host resolveHoDofOwnership above (MIN-RANK-AMONG-HOLDERS)
// but with NO host round-trip on the per-DOF arrays. The inputs are the device columns
// buildDistributedGpu already produces (d_dofKind/G0..G3/pos/shared/owner) -- the GPU
// path can hand them straight here instead of copying them down and re-sorting on host.
//
// DESIGN = TARGETED P2P (not zero-MPI). Zero-MPI is NOT possible with the current
// numbering: a rank numbers only its OWNED elements, so it cannot know which OTHER ranks
// hold a given shared edge/face without an exchange (ghost-element owners are not plumbed
// in; elemOwner == myRank everywhere). So a peer key exchange is unavoidable -- the same
// conclusion both infra analyses reached.
//
// TARGETING = broadcast each shared key to ALL P1 peers (exactly what the host path does).
// dofShared is already gated on every defining corner being P1-shared, so the shared set
// is the rank SURFACE (tiny vs numDof) -- broadcasting it to the P1 peer list is cheap and,
// critically, CANNOT miss a holder: any rank that holds a shared edge/face shares all its
// corners and is therefore a P1 peer that receives the key. (Per-corner intersection would
// be fewer messages but needs a global->local corner inverse map; a wrong inverse drops a
// holder -> wrong owner -> orphan. We keep the host's safe all-peer fan-out.)
//
// SPEED: the dominant cost is the MATCH (per received key, a search over the local shared
// keys; the all-peer path receives ~np*nShared keys). Two changes make it DRAM-friendly:
//   1. The 6 key fields (kind,g0..g3,pos = 48 bytes) are stored INTERLEAVED as one array
//      of stride 6: key i at [6*i .. 6*i+5]. One probe then reads ~one cache line instead
//      of touching 6 separate column arrays (6 cache lines, uncoalesced).
//   2. We need EQUALITY matching, not lexicographic order. So instead of a 6-pass
//      lexicographic radix sort, we sort by a SINGLE 64-bit hash of the key (ONE radix
//      pass). The match binary-searches the hash and scans the equal-hash run, FULL-
//      comparing all 6 fields -- a hash collision can never produce a false match, so the
//      result stays bit-identical to the lexicographic path.
//
// MPI buffers are device pointers (GPUDirect), same mechanics as HoHalo::forwardDevice.
// Tags 0x4852 (counts) / 0x4853 (key columns); the host resolve uses 0x4850/0x4851, so a
// mixed run never crosses wires.

// 64-bit hash of one 6-field key. Order-preserving is NOT needed (we sort by it only to
// group equal keys for binary search); we only need EQUAL keys to hash EQUAL (trivially
// true) and a low collision rate so the equal-hash scan in the match stays short. A
// collision is still correct -- the match full-compares all 6 fields on a hash hit.
#define HO_KEY_HASH_MIX(h, v) ((h) ^= (v) + 0x9e3779b97f4a7c15ull + ((h) << 6) + ((h) >> 2))
__device__ __forceinline__ uint64_t hoKeyHash(uint64_t k0, uint64_t k1, uint64_t k2,
                                              uint64_t k3, uint64_t k4, uint64_t k5)
{
    uint64_t h = 0x100000001b3ull;
    HO_KEY_HASH_MIX(h, k0); HO_KEY_HASH_MIX(h, k1); HO_KEY_HASH_MIX(h, k2);
    HO_KEY_HASH_MIX(h, k3); HO_KEY_HASH_MIX(h, k4); HO_KEY_HASH_MIX(h, k5);
    return h;
}
#undef HO_KEY_HASH_MIX

// Per-stage timing, gated by env MARS_HO_RESOLVE_TIMING (quiet in production). Each stage
// cudaDeviceSynchronizes before stop so the time is real device work, not just launch.
struct HoResolveTimer {
    bool   on;
    int    rank;
    double t0;
    HoResolveTimer(int r) : rank(r) { on = std::getenv("MARS_HO_RESOLVE_TIMING") != nullptr; t0 = MPI_Wtime(); }
    void lap(const char* stage) {
        if (!on) { t0 = MPI_Wtime(); return; }
        cudaDeviceSynchronize();
        double t1 = MPI_Wtime();
        if (rank == 0) { printf("[ho-resolve] %-22s %8.3f ms\n", stage, (t1 - t0) * 1e3); fflush(stdout); }
        t0 = t1;
    }
};

// Pack one shared DOF's key (kind,g0..g3,pos) INTERLEAVED at [6*s .. 6*s+5] + its hash and
// local dof id at slot s. Interleaved so the later binary search reads 48 contiguous bytes
// per probe (one cache line) instead of 6 strided column reads.
__global__ void hoResolvePackKernel(const long* __restrict__ kind, const long* __restrict__ g0,
                                    const long* __restrict__ g1,  const long* __restrict__ g2,
                                    const long* __restrict__ g3,  const int*  __restrict__ pos,
                                    const int*  __restrict__ sharedDof, const long* __restrict__ slot,
                                    long numDof,
                                    uint64_t* __restrict__ keyIL,      // interleaved, stride 6
                                    uint64_t* __restrict__ keyHash,
                                    int* __restrict__ localDof)
{
    long d = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (d >= numDof || !sharedDof[d]) return;
    long s = slot[d];                          // exclusive-scan rank among shared DOF
    uint64_t k0 = (uint64_t)kind[d], k1 = (uint64_t)g0[d], k2 = (uint64_t)g1[d];
    uint64_t k3 = (uint64_t)g2[d],   k4 = (uint64_t)g3[d], k5 = (uint64_t)pos[d];
    uint64_t* o = keyIL + 6L * s;
    o[0] = k0; o[1] = k1; o[2] = k2; o[3] = k3; o[4] = k4; o[5] = k5;
    keyHash[s]  = hoKeyHash(k0, k1, k2, k3, k4, k5);
    localDof[s] = (int)d;
}

// For each received key, binary-search the local SORTED shared keys BY HASH; on an exact
// 6-field match, atomicMin the sender rank onto that DOF's owner. Both the local and the
// received keys are INTERLEAVED (stride 6), so each probe reads one contiguous 48-byte key.
// The sorted hash array (sHash) is the search key; on a hash hit we scan the equal-hash run
// and full-compare -- a collision never produces a false match (result bit-identical to a
// lexicographic search).
__global__ void hoResolveMatchKernel(const uint64_t* __restrict__ recvIL,   // interleaved, stride 6
                                     int nRecv, int senderRank,
                                     const uint64_t* __restrict__ sortedIL,  // interleaved, stride 6
                                     const uint64_t* __restrict__ sHash,
                                     const int* __restrict__ sortedLocalDof, int nShared,
                                     int* __restrict__ dofOwner)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nRecv) return;
    const uint64_t* r = recvIL + 6L * i;
    uint64_t k0 = r[0], k1 = r[1], k2 = r[2], k3 = r[3], k4 = r[4], k5 = r[5];
    uint64_t h = hoKeyHash(k0, k1, k2, k3, k4, k5);
    // lower_bound on the hash over the hash-sorted shared keys.
    int lo = 0, hi = nShared;
    while (lo < hi) { int m = (lo + hi) >> 1; if (sHash[m] < h) lo = m + 1; else hi = m; }
    for (int j = lo; j < nShared && sHash[j] == h; ++j) {
        const uint64_t* s = sortedIL + 6L * j;
        if (s[0]==k0 && s[1]==k1 && s[2]==k2 && s[3]==k3 && s[4]==k4 && s[5]==k5) {
            atomicMin(&dofOwner[sortedLocalDof[j]], senderRank);
            return;                                   // keys are unique per DOF
        }
    }
}

// Gather the interleaved keys for the picked slots into a contiguous interleaved send slab
// (same stride-6 layout the match expects). One thread per picked key copies its 6 words.
__global__ void hoKeyGatherInterleaveKernel(const uint64_t* __restrict__ sortedIL,
                                            const int* __restrict__ pickIdx, int nPick,
                                            uint64_t* __restrict__ out)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nPick) return;
    const uint64_t* s = sortedIL + 6L * pickIdx[i];
    uint64_t* o = out + 6L * i;
    o[0]=s[0]; o[1]=s[1]; o[2]=s[2]; o[3]=s[3]; o[4]=s[4]; o[5]=s[5];
}

// Build the interleaved sorted local keys + hash + localDof from the per-DOF device columns.
// Sorts by a single 64-bit hash (ONE radix pass) -- equality matching only needs equal keys
// grouped, not lexicographic order. Returns nShared; fills d_sortedIL (6*nShared), d_sHash,
// d_sLocal. d_slot must hold the exclusive-scan slot of each shared DOF.
inline long hoBuildSortedSharedKeys(long numDof, const long* d_dofKind, const long* d_dofG0,
                                    const long* d_dofG1, const long* d_dofG2, const long* d_dofG3,
                                    const int* d_dofPos, const int* d_dofShared,
                                    const thrust::device_vector<long>& d_slot, long nShared,
                                    thrust::device_vector<uint64_t>& d_sortedIL,
                                    thrust::device_vector<uint64_t>& d_sHash,
                                    thrust::device_vector<int>&      d_sLocal)
{
    d_sortedIL.resize(6L * (nShared > 0 ? nShared : 1));
    d_sHash.resize(nShared > 0 ? nShared : 1);
    d_sLocal.resize(nShared > 0 ? nShared : 1);
    if (nShared <= 0) return nShared;

    // pack interleaved (unsorted) + hash + local id
    thrust::device_vector<uint64_t> packIL(6L * nShared);
    thrust::device_vector<uint64_t> packHash(nShared);
    thrust::device_vector<int>      packLocal(nShared);
    {
        int b = 256; long g = (numDof + b - 1) / b;
        hoResolvePackKernel<<<g, b>>>(d_dofKind, d_dofG0, d_dofG1, d_dofG2, d_dofG3, d_dofPos,
                                      d_dofShared, thrust::raw_pointer_cast(d_slot.data()), numDof,
                                      thrust::raw_pointer_cast(packIL.data()),
                                      thrust::raw_pointer_cast(packHash.data()),
                                      thrust::raw_pointer_cast(packLocal.data()));
    }

    // sort a permutation by hash (one radix pass), then gather the interleaved keys +
    // hash + local id by it. We sort the index, not the 48-byte keys, so the heavy payload
    // moves exactly once.
    thrust::device_vector<int> idx(nShared);
    thrust::sequence(idx.begin(), idx.end());
    thrust::stable_sort_by_key(packHash.begin(), packHash.end(), idx.begin());  // packHash now sorted
    d_sHash.swap(packHash);
    thrust::gather(idx.begin(), idx.end(), packLocal.begin(), d_sLocal.begin());
    // gather the interleaved keys: each output key i = packIL[6*idx[i] .. +5]
    {
        int b = 256, g = (int)((nShared + b - 1) / b);
        hoKeyGatherInterleaveKernel<<<g, b>>>(thrust::raw_pointer_cast(packIL.data()),
                                              thrust::raw_pointer_cast(idx.data()), (int)nShared,
                                              thrust::raw_pointer_cast(d_sortedIL.data()));
    }
    return nShared;
}

// Resolve on device. d_dofOwner is mutated in place (atomicMin). All other inputs are the
// device columns from buildDistributedGpu. Returns with the device columns intact.
inline void resolveHoDofOwnershipGpu(long numDof,
                                     const long* d_dofKind, const long* d_dofG0, const long* d_dofG1,
                                     const long* d_dofG2,   const long* d_dofG3, const int* d_dofPos,
                                     const int*  d_dofShared, int* d_dofOwner,
                                     int myRank, const std::vector<int>& peers)
{
    (void)myRank;
    const int np = (int)peers.size();
    HoResolveTimer timer(myRank);

    // count shared DOF and assign each a contiguous slot. d_dofShared is a 0/1 lane
    // (atomicOr'd with 1 in buildDistributedGpu), so an exclusive scan gives each shared
    // DOF its slot and the last (slot + flag) is the shared count -- no host reduce.
    thrust::device_vector<long> d_slot(numDof);
    {
        thrust::device_ptr<const int> sh(d_dofShared);
        thrust::exclusive_scan(sh, sh + numDof, d_slot.begin(), (long)0);
    }
    long nShared = 0;
    {
        int lastShared = 0; long lastSlot = 0;
        thrust::copy(thrust::device_pointer_cast(d_dofShared) + (numDof - 1),
                     thrust::device_pointer_cast(d_dofShared) + numDof, &lastShared);
        thrust::copy(d_slot.begin() + (numDof - 1), d_slot.begin() + numDof, &lastSlot);
        nShared = lastSlot + (lastShared ? 1 : 0);
    }
    timer.lap("scan-slots");

    // pack + hash-sort the shared keys into ONE interleaved array (stride 6) -- one radix
    // pass on a 64-bit hash, not the old 6-pass lexicographic sort.
    thrust::device_vector<uint64_t> sortedIL, sHash;
    thrust::device_vector<int>      sLocal;
    hoBuildSortedSharedKeys(numDof, d_dofKind, d_dofG0, d_dofG1, d_dofG2, d_dofG3, d_dofPos,
                            d_dofShared, d_slot, nShared, sortedIL, sHash, sLocal);
    d_slot = thrust::device_vector<long>();   // numDof-sized scratch no longer needed
    timer.lap("pack+hash-sort");

    const int myCnt = (int)nShared;

    // exchange shared-key counts with every P1 peer (host ints; 4 bytes/peer, trivial)
    std::vector<int> rc(np, 0);
    {
        std::vector<MPI_Request> rq; rq.reserve(2 * np);
        for (int i = 0; i < np; ++i) { MPI_Request r; MPI_Irecv(&rc[i], 1, MPI_INT, peers[i], 0x4852, MPI_COMM_WORLD, &r); rq.push_back(r); }
        for (int i = 0; i < np; ++i) { MPI_Request r; MPI_Isend(&myCnt, 1, MPI_INT, peers[i], 0x4852, MPI_COMM_WORLD, &r); rq.push_back(r); }
        if (!rq.empty()) MPI_Waitall((int)rq.size(), rq.data(), MPI_STATUSES_IGNORE);
    }
    timer.lap("count-exchange");

    // GPUDirect exchange of the interleaved keys. One contiguous device recv slab per peer,
    // 6 uint64/key. Send the SAME local keys to every peer (broadcast); each peer searches
    // its own local keys and atomicMins our rank where it matches. The sorted local keys are
    // already interleaved, so they ARE the send buffer (no per-column copy needed).
    std::vector<int> roff(np + 1, 0);
    for (int i = 0; i < np; ++i) roff[i + 1] = roff[i] + rc[i];
    const int totalRecv = roff[np];

    thrust::device_vector<uint64_t> recv(6 * (size_t)(totalRecv > 0 ? totalRecv : 1));
    uint64_t* d_recv = thrust::raw_pointer_cast(recv.data());

    cudaDeviceSynchronize();   // sortedIL is built by an async kernel -> finish before GPUDirect reads it
    {
        std::vector<MPI_Request> rq; rq.reserve(2 * np);
        for (int i = 0; i < np; ++i)
            if (rc[i] > 0) { MPI_Request r; MPI_Irecv(d_recv + 6L * roff[i], rc[i] * 6, MPI_UNSIGNED_LONG_LONG, peers[i], 0x4853, MPI_COMM_WORLD, &r); rq.push_back(r); }
        for (int i = 0; i < np; ++i)
            if (myCnt > 0) { MPI_Request r; MPI_Isend(thrust::raw_pointer_cast(sortedIL.data()), myCnt * 6, MPI_UNSIGNED_LONG_LONG, peers[i], 0x4853, MPI_COMM_WORLD, &r); rq.push_back(r); }
        if (!rq.empty()) MPI_Waitall((int)rq.size(), rq.data(), MPI_STATUSES_IGNORE);
    }
    timer.lap("key-exchange");

    // match each peer's received keys against our local sorted keys -> atomicMin owner.
    if (nShared > 0) {
        const uint64_t* sIL  = thrust::raw_pointer_cast(sortedIL.data());
        const uint64_t* sH   = thrust::raw_pointer_cast(sHash.data());
        const int*      sLoc = thrust::raw_pointer_cast(sLocal.data());
        for (int i = 0; i < np; ++i) {
            if (rc[i] == 0) continue;
            const uint64_t* base = d_recv + 6L * roff[i];
            int nr = rc[i];
            int b = 256, g = (nr + b - 1) / b;
            hoResolveMatchKernel<<<g, b>>>(base, nr, peers[i], sIL, sH, sLoc, (int)nShared, d_dofOwner);
        }
        cudaDeviceSynchronize();
    }
    timer.lap("match");
}

// ============================================================================
// TARGETED P2P resolveHoDofOwnership (env MARS_HO_GPU_OWNERSHIP_P2P).
//
// Same result as the all-peer resolveHoDofOwnershipGpu above (MIN-RANK-AMONG-HOLDERS),
// but it DROPS the exchange volume: the all-peer path broadcasts every shared key to
// EVERY P1 peer (O(shared x peers)); this path sends each shared key only to the peers
// that can possibly hold it.
//
// WHY TARGETING IS EXACT (cannot drop a holder). Target each shared key by the SEND list
// of a corner THIS RANK OWNS, NOT by recv+send masks. The proof:
//   Let H be any holder of a shared edge/face DOF. H owns an element that contains the
//   entity, so H shares ALL of the entity's corners with us. Pick a corner C_own of the
//   entity that THIS rank owns. Then H has C_own as a GHOST it receives FROM this rank
//   -> H is in this rank's node-halo SEND list for C_own. So the SEND-list peers of
//   C_own are a SUPERSET of every holder of the DOF. Sending the key to them reaches
//   every holder. (The earlier recv+send-mask INTERSECTION was BROKEN: for a GHOST
//   corner C the recv mask only carries the single owner q, never the OTHER ranks that
//   also ghost C, so a 3+-rank-corner co-ghost holder was dropped from the intersection
//   -> silent orphan. The owned-corner send list has no such gap.)
//
// WHY MIN-RANK IS GLOBALLY CONSISTENT (symmetry). atomicMin is order-free, so the result
// is correct iff every holder's rank reaches every co-holder:
//   - A holder that OWNS a corner of the DOF sends the key via that corner's send list,
//     which (by the proof above) reaches ALL co-holders.
//   - A holder that owns NO corner of the DOF falls back to the all-peer broadcast, which
//     trivially reaches all co-holders.
// Either way every holder's rank reaches every co-holder -> the atomicMin over received
// keys is bit-identical to the all-peer result.
//
// GLOBAL->LOCAL CORNER INVERSE: the key carries g0..g3 = SORTED GLOBAL corner ids
// (cornerGid is local->global). To find a corner's local id (to read its owner + send
// mask) we build the inverse ONLY over node-halo-shared corners (non-empty peer
// participation) -- the only corners a shared edge/face key can reference. The inverse is
// a sorted (gid -> localCorner) table; each key's g* is binary-searched in it.
//
// ALL-PEERS FALLBACK (mandatory, never drop): a key gets the all-peer mask when it has NO
// resolvable owned corner -- either no corner of the entity is owned by this rank, or a
// global->local corner lookup misses (corner not in our shared-corner table, or >64
// peers so a bit can't be set). The fallback can only ADD messages, never drop a holder.
//
// SAFETY GUARD: with MARS_HO_VERIFY_OWNERSHIP the caller runs BOTH this targeted path and
// the all-peer resolveHoDofOwnershipGpu on a clone of d_dofOwner and asserts bit-identical
// (MPI_Abort on mismatch) -- the proof the targeting never drops a holder.
//
// Peer-set representation: a 64-bit mask, bit p = "peers[p] ghosts this corner from me".
// P1 halo peer counts are small (mesh face/edge neighbours), so 64 covers real runs; >64
// peers falls back to all-peer per the rule above. Tags 0x4854 (counts) / 0x4855 (key
// columns) -- distinct from the all-peer 0x4852/0x4853 and the host 0x4850/0x4851.
// ============================================================================

// Build, per LOCAL corner, the 64-bit "SEND-to-peer" mask from the node-halo SEND CSR.
// A corner appearing in peer p's SEND node list is ghosted BY peer p FROM this rank
// (peer p holds it as a ghost) -> bit p. For a corner THIS rank owns, these bits are
// exactly the peers that ghost it = a SUPERSET of every holder of any entity built on
// that corner. (Only the SEND list -- NOT recv -- so we never confuse "peer q owns this
// ghost corner of mine" with "peer p also ghosts it".) Bits beyond 63 are not
// representable; keys whose owned corner needs such a peer fall back to all-peer.
__global__ void hoCornerSendMaskKernel(const int* __restrict__ sendNodeIds, int n,
                                       int peerBit, uint64_t* __restrict__ cornerSendMask)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    int lc = sendNodeIds[i];
    if (lc >= 0) atomicOr((unsigned long long*)&cornerSendMask[lc], (unsigned long long)(1ull << peerBit));
}

// Mark (0/1) every local corner appearing in a node list, regardless of peer. Used to
// flag corners that participate in the node halo (recv OR send) so the global->local
// inverse covers every corner a shared key can reference.
__global__ void hoMarkFlagKernel(const int* __restrict__ nodeIds, int n, int* __restrict__ flag)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    int lc = nodeIds[i];
    if (lc >= 0) flag[lc] = 1;
}

// Build the sorted global->local inverse over SHARED corners only. Emit (gid, localId)
// for every corner that participates in the node halo (recv OR send); the host sorts by
// gid for binary search. The participation flag is supplied as a 0/1 lane so the same
// emit serves both the recv+send union (for resolvability) and is consistent with slots.
__global__ void hoSharedCornerEmitKernel(const long* __restrict__ cornerGid,
                                         const int* __restrict__ cornerFlag,
                                         int nCornerLocal,
                                         const long* __restrict__ slot,
                                         long* __restrict__ outGid, int* __restrict__ outLocal)
{
    int lc = blockIdx.x * blockDim.x + threadIdx.x;
    if (lc >= nCornerLocal || cornerFlag[lc] == 0) return;
    long s = slot[lc];
    outGid[s]   = cornerGid[lc];
    outLocal[s] = lc;
}

// For each shared key, compute its target peer mask = the SEND mask of a corner THIS rank
// OWNS. The key's g0..g3 are the entity's sorted GLOBAL corner ids; we binary-search each
// in the sorted inverse, read its local id, and take the FIRST corner with
// cornerOwner[lc]==myRank -> its cornerSendMask is the target (a superset of all holders).
// If NO owned corner is found (none owned, or a corner lookup misses), the key falls back
// to ALL peers -- safe, can only add messages. keyIL = the interleaved sorted keys (stride
// 6): word 0 = kind, words 1..4 = g0..g3, word 5 = pos.
__global__ void hoKeyOwnedCornerTargetKernel(const uint64_t* __restrict__ keyIL,
                                             int nShared,
                                             const long* __restrict__ invGid, const int* __restrict__ invLocal,
                                             int nInv,
                                             const int* __restrict__ cornerOwner,
                                             const uint64_t* __restrict__ cornerSendMask,
                                             int myRank, uint64_t allPeersMask,
                                             uint64_t* __restrict__ keyMask)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nShared) return;

    const uint64_t* k = keyIL + 6L * i;
    // kind 1 = edge (2 corners g0,g1), kind 2 = face (4 corners g0..g3). Corner DOF
    // (kind 0) and interiors (kind 3) are never in the shared set; guard kind anyway.
    long kind = (long)k[0];
    int  nc   = (kind == 2) ? 4 : 2;
    long g[4] = { (long)k[1], (long)k[2], (long)k[3], (long)k[4] };

    uint64_t mask = allPeersMask;   // fallback: no owned corner found -> all peers
    for (int c = 0; c < nc; ++c) {
        long gid = g[c];
        int lo = 0, hi = nInv, found = -1;
        while (lo < hi) {
            int m = (lo + hi) >> 1;
            if (invGid[m] < gid) lo = m + 1;
            else if (invGid[m] > gid) hi = m;
            else { found = m; break; }
        }
        if (found < 0) continue;                       // corner not resolvable -> try next
        int lc = invLocal[found];
        if (cornerOwner[lc] == myRank) {               // an OWNED corner of this entity
            mask = cornerSendMask[lc];                 // its ghosting peers = holder superset
            break;
        }
    }
    keyMask[i] = mask;
}

// SINGLE-PASS per-peer bucketing. The old path ran copy_if over all nShared keys ONCE PER
// PEER (np passes -> O(np*nShared)). Instead each key emits one (peerBit, keyIndex) entry
// per set bit of its target mask; sorting those entries by peerBit groups every peer's keys
// contiguously, so the whole bucketing is O(nShared + totalOut). These two kernels feed that
// expand-then-radix-sort.

// Per key, count set bits of its target mask restricted to the representable peers (< 64
// AND < np). The exclusive scan of these counts gives each key its write offset into the
// expanded (peerBit, keyIndex) arrays, and the total is the expanded length.
__global__ void hoKeyMaskPopcountKernel(const uint64_t* __restrict__ keyMask, int nShared,
                                        int np, int* __restrict__ outCount)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nShared) return;
    uint64_t m = keyMask[i];
    if (np < 64) m &= ((np <= 0) ? 0ull : ((1ull << np) - 1ull));   // ignore non-peer bits
    outCount[i] = __popcll((unsigned long long)m);
}

// Per key, write one (peerBit, keyIndex) entry per set bit of its (representable) mask,
// starting at its scan offset. After this the entries cover every (peer, key) pair exactly
// once; a radix sort by peerBit then buckets them.
__global__ void hoKeyExpandKernel(const uint64_t* __restrict__ keyMask, int nShared,
                                  int np, const int* __restrict__ scanOff,
                                  int* __restrict__ outPeer, int* __restrict__ outKeyIdx)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nShared) return;
    uint64_t m = keyMask[i];
    if (np < 64) m &= ((np <= 0) ? 0ull : ((1ull << np) - 1ull));
    int o = scanOff[i];
    while (m) {
        int p = __ffsll((unsigned long long)m) - 1;   // lowest set bit
        outPeer[o]   = p;
        outKeyIdx[o] = i;
        ++o;
        m &= m - 1;                                    // clear lowest set bit
    }
}

// TARGETED resolve. Identical inputs/outputs to resolveHoDofOwnershipGpu, plus the
// node-halo CSR + per-corner owner so it can target by an OWNED corner's send list.
// d_cornerGid = device local->global corner id table; d_cornerOwner = device per-corner
// owning rank (myRank for owned corners, the owning peer for ghost corners). nodeHalo* =
// the topology's recv/send node lists + offsets (host offsets, device node-id arrays),
// peer order MATCHING `peers`.
inline void resolveHoDofOwnershipGpuTargeted(
        long numDof,
        const long* d_dofKind, const long* d_dofG0, const long* d_dofG1,
        const long* d_dofG2,   const long* d_dofG3, const int* d_dofPos,
        const int*  d_dofShared, int* d_dofOwner,
        int myRank, const std::vector<int>& peers,
        const long* d_cornerGid, const int* d_cornerOwner, int nCornerLocal,
        const int* d_recvNodeIds, const std::vector<int>& recvOffsets,
        const int* d_sendNodeIds, const std::vector<int>& sendOffsets)
{
    const int np = (int)peers.size();
    HoResolveTimer timer(myRank);

    // ---- count shared DOF + assign slots (identical to the all-peer path) ----
    thrust::device_vector<long> d_slot(numDof);
    {
        thrust::device_ptr<const int> sh(d_dofShared);
        thrust::exclusive_scan(sh, sh + numDof, d_slot.begin(), (long)0);
    }
    long nShared = 0;
    {
        int lastShared = 0; long lastSlot = 0;
        thrust::copy(thrust::device_pointer_cast(d_dofShared) + (numDof - 1),
                     thrust::device_pointer_cast(d_dofShared) + numDof, &lastShared);
        thrust::copy(d_slot.begin() + (numDof - 1), d_slot.begin() + numDof, &lastSlot);
        nShared = lastSlot + (lastShared ? 1 : 0);
    }
    timer.lap("scan-slots");

    // ---- pack + hash-sort the shared keys into ONE interleaved array (same as all-peer) ----
    thrust::device_vector<uint64_t> sortedIL, sHash;
    thrust::device_vector<int>      sLocal;
    hoBuildSortedSharedKeys(numDof, d_dofKind, d_dofG0, d_dofG1, d_dofG2, d_dofG3, d_dofPos,
                            d_dofShared, d_slot, nShared, sortedIL, sHash, sLocal);
    d_slot = thrust::device_vector<long>();   // numDof-sized scratch no longer needed
    timer.lap("pack+hash-sort");

    // ---- per-corner SEND mask (bit p = peers[p] ghosts this corner FROM me) ----
    // ONLY the send list: for a corner this rank OWNS, these bits are exactly the peers
    // that ghost it = a SUPERSET of every holder of any entity built on that corner.
    // (Using recv too -- the old broken path -- mixes "peer q owns this ghost corner of
    // mine" into the mask and at a 3+-rank corner drops the co-ghost holders.) Peers
    // beyond bit 63 cannot be set -> keys whose owned corner needs them fall back below.
    thrust::device_vector<uint64_t> d_cornerSendMask(nCornerLocal, 0ull);
    {
        uint64_t* p_mask = thrust::raw_pointer_cast(d_cornerSendMask.data());
        for (int p = 0; p < np && p < 64; ++p) {
            int sc = sendOffsets[p + 1] - sendOffsets[p];
            if (sc > 0) { int b = 256, g = (sc + b - 1) / b;
                hoCornerSendMaskKernel<<<g, b>>>(d_sendNodeIds + sendOffsets[p], sc, p, p_mask); }
        }
        cudaDeviceSynchronize();
    }

    // allPeersMask = every representable peer bit set (the safe fallback target). With
    // np>64 peers we cannot mask them all; set all 64 bits AND route every key with no
    // resolvable owned corner to all peers below (the >=64 send-slab branch).
    const uint64_t allPeersMask = (np >= 64) ? ~0ull
                                             : ((np == 0) ? 0ull : ((1ull << np) - 1ull));

    // ---- corner participation flag (recv OR send) for the global->local inverse ----
    // The inverse must resolve ANY corner a shared key references so we can read its owner
    // + send mask, including ghost corners (recv-only) that anchor an entity whose OTHER
    // corner this rank owns. So flag corners touched by recv OR send, not by send alone.
    thrust::device_vector<int> d_cornerFlag(nCornerLocal, 0);
    {
        int* p_flag = thrust::raw_pointer_cast(d_cornerFlag.data());
        auto flagFrom = [&](const int* nodeIds, const std::vector<int>& off) {
            for (int p = 0; p < np; ++p) {
                int c = off[p + 1] - off[p];
                if (c > 0) { int b = 256, g = (c + b - 1) / b;
                    hoMarkFlagKernel<<<g, b>>>(nodeIds + off[p], c, p_flag); }
            }
        };
        flagFrom(d_recvNodeIds, recvOffsets);
        flagFrom(d_sendNodeIds, sendOffsets);
        cudaDeviceSynchronize();
    }

    // ---- build the sorted global->local inverse over PARTICIPATING corners ----
    thrust::device_vector<long> d_invSlot(nCornerLocal);
    thrust::exclusive_scan(d_cornerFlag.begin(), d_cornerFlag.end(), d_invSlot.begin(), (long)0);
    long nInv = 0;
    {
        int lastFlag = 0; long lastSlot = 0;
        thrust::copy(d_cornerFlag.begin() + (nCornerLocal - 1), d_cornerFlag.end(), &lastFlag);
        thrust::copy(d_invSlot.begin() + (nCornerLocal - 1), d_invSlot.begin() + nCornerLocal, &lastSlot);
        nInv = lastSlot + lastFlag;
    }
    thrust::device_vector<long> d_invGid(nInv > 0 ? nInv : 1);
    thrust::device_vector<int>  d_invLocal(nInv > 0 ? nInv : 1);
    if (nInv > 0) {
        int b = 256, g = (nCornerLocal + b - 1) / b;
        hoSharedCornerEmitKernel<<<g, b>>>(d_cornerGid, thrust::raw_pointer_cast(d_cornerFlag.data()),
                                           nCornerLocal, thrust::raw_pointer_cast(d_invSlot.data()),
                                           thrust::raw_pointer_cast(d_invGid.data()),
                                           thrust::raw_pointer_cast(d_invLocal.data()));
        // sort the inverse by gid so the per-key search is a binary search
        thrust::stable_sort_by_key(d_invGid.begin(), d_invGid.end(), d_invLocal.begin());
        cudaDeviceSynchronize();
    }
    timer.lap("corner-masks+inverse");

    // ---- per-key target mask = the SEND mask of an OWNED corner of the key's entity ----
    thrust::device_vector<uint64_t> d_keyMask(nShared > 0 ? nShared : 1, 0ull);
    if (nShared > 0) {
        int b = 256, g = (int)((nShared + b - 1) / b);
        hoKeyOwnedCornerTargetKernel<<<g, b>>>(
            thrust::raw_pointer_cast(sortedIL.data()), (int)nShared,
            thrust::raw_pointer_cast(d_invGid.data()), thrust::raw_pointer_cast(d_invLocal.data()),
            (int)nInv, d_cornerOwner, thrust::raw_pointer_cast(d_cornerSendMask.data()),
            myRank, allPeersMask, thrust::raw_pointer_cast(d_keyMask.data()));
        cudaDeviceSynchronize();
    }
    timer.lap("per-key-target");

    // ---- SINGLE-PASS per-peer bucketing (was np passes of copy_if) ----
    // Expand each key to one (peerBit, keyIndex) entry per set mask bit, radix-sort the
    // entries by peerBit (groups every peer's keys contiguously), then gather the
    // interleaved keys in that order ONCE. Per-peer slabs are slices of the one big slab.
    // Peers >= 64 (only when np > 64) are not in the mask -> they get ALL shared keys, the
    // safe all-peer fallback, appended after the bucketed peers.
    std::vector<int> sendCnt(np, 0);
    const int nBucketPeers = (np < 64) ? np : 64;     // peers representable in the 64-bit mask
    thrust::device_vector<uint64_t> bucketSlab;        // [6 * totalOut], peer-major
    std::vector<int> bucketOff(nBucketPeers + 1, 0);   // word-key offsets per bucketed peer
    if (nShared > 0 && nBucketPeers > 0) {
        // per-key set-bit count -> exclusive scan = write offsets; total = expanded length
        thrust::device_vector<int> cnt(nShared);
        {
            int b = 256, g = (int)((nShared + b - 1) / b);
            hoKeyMaskPopcountKernel<<<g, b>>>(thrust::raw_pointer_cast(d_keyMask.data()),
                                              (int)nShared, np, thrust::raw_pointer_cast(cnt.data()));
        }
        thrust::device_vector<int> off(nShared);
        thrust::exclusive_scan(cnt.begin(), cnt.end(), off.begin(), 0);
        int lastCnt = 0, lastOff = 0;
        thrust::copy(cnt.begin() + (nShared - 1), cnt.end(), &lastCnt);
        thrust::copy(off.begin() + (nShared - 1), off.begin() + nShared, &lastOff);
        const int totalOut = lastOff + lastCnt;

        if (totalOut > 0) {
            thrust::device_vector<int> expPeer(totalOut), expKey(totalOut);
            {
                int b = 256, g = (int)((nShared + b - 1) / b);
                hoKeyExpandKernel<<<g, b>>>(thrust::raw_pointer_cast(d_keyMask.data()), (int)nShared,
                                            np, thrust::raw_pointer_cast(off.data()),
                                            thrust::raw_pointer_cast(expPeer.data()),
                                            thrust::raw_pointer_cast(expKey.data()));
            }
            // group by peer: one radix sort of expKey keyed on expPeer
            thrust::stable_sort_by_key(expPeer.begin(), expPeer.end(), expKey.begin());
            // gather the interleaved keys in the sorted (peer-major) order -> one big slab
            bucketSlab.resize(6L * totalOut);
            {
                int b = 256, g = (totalOut + b - 1) / b;
                hoKeyGatherInterleaveKernel<<<g, b>>>(thrust::raw_pointer_cast(sortedIL.data()),
                                                      thrust::raw_pointer_cast(expKey.data()), totalOut,
                                                      thrust::raw_pointer_cast(bucketSlab.data()));
            }
            // per-peer counts: reduce_by_key over the sorted peer ids -> (peer, count) segs
            thrust::device_vector<int> peerId(nBucketPeers), peerCnt(nBucketPeers);
            auto endp = thrust::reduce_by_key(expPeer.begin(), expPeer.end(),
                                              thrust::make_constant_iterator(1),
                                              peerId.begin(), peerCnt.begin());
            int nSeg = (int)(endp.first - peerId.begin());
            std::vector<int> hId(nSeg), hCnt(nSeg);
            thrust::copy(peerId.begin(), peerId.begin() + nSeg, hId.begin());
            thrust::copy(peerCnt.begin(), peerCnt.begin() + nSeg, hCnt.begin());
            for (int s = 0; s < nSeg; ++s) sendCnt[hId[s]] = hCnt[s];
        }
        cudaDeviceSynchronize();
    }
    for (int p = 0; p < nBucketPeers; ++p) bucketOff[p + 1] = bucketOff[p] + sendCnt[p];

    // Peers >= 64 (high peer-count fallback) get the full interleaved key array directly.
    for (int p = 64; p < np; ++p) sendCnt[p] = (int)nShared;
    timer.lap("per-peer-compaction");

    // ---- exchange per-peer counts (asymmetric now: each peer gets a different count) ----
    std::vector<int> rc(np, 0);
    {
        std::vector<MPI_Request> rq; rq.reserve(2 * np);
        for (int i = 0; i < np; ++i) { MPI_Request r; MPI_Irecv(&rc[i], 1, MPI_INT, peers[i], 0x4854, MPI_COMM_WORLD, &r); rq.push_back(r); }
        for (int i = 0; i < np; ++i) { MPI_Request r; MPI_Isend(&sendCnt[i], 1, MPI_INT, peers[i], 0x4854, MPI_COMM_WORLD, &r); rq.push_back(r); }
        if (!rq.empty()) MPI_Waitall((int)rq.size(), rq.data(), MPI_STATUSES_IGNORE);
    }
    timer.lap("count-exchange");

    // ---- GPUDirect exchange of the targeted interleaved keys ----
    // A bucketed peer (p < 64) sends its slice of bucketSlab; a fallback peer (p >= 64)
    // sends the full sorted interleaved key array. Both are stride-6 interleaved.
    std::vector<int> roff(np + 1, 0);
    for (int i = 0; i < np; ++i) roff[i + 1] = roff[i] + rc[i];
    const int totalRecv = roff[np];
    thrust::device_vector<uint64_t> recv(6 * (size_t)(totalRecv > 0 ? totalRecv : 1));
    uint64_t* d_recv = thrust::raw_pointer_cast(recv.data());

    cudaDeviceSynchronize();   // bucketSlab/sortedIL built by async kernels -> finish before GPUDirect reads them
    {
        std::vector<MPI_Request> rq; rq.reserve(2 * np);
        for (int i = 0; i < np; ++i)
            if (rc[i] > 0) { MPI_Request r; MPI_Irecv(d_recv + 6L * roff[i], rc[i] * 6, MPI_UNSIGNED_LONG_LONG, peers[i], 0x4855, MPI_COMM_WORLD, &r); rq.push_back(r); }
        for (int i = 0; i < np; ++i) {
            if (sendCnt[i] <= 0) continue;
            const uint64_t* src = (i < nBucketPeers)
                ? thrust::raw_pointer_cast(bucketSlab.data()) + 6L * bucketOff[i]
                : thrust::raw_pointer_cast(sortedIL.data());     // p >= 64 fallback: full array
            MPI_Request r; MPI_Isend((void*)src, sendCnt[i] * 6, MPI_UNSIGNED_LONG_LONG, peers[i], 0x4855, MPI_COMM_WORLD, &r); rq.push_back(r);
        }
        if (!rq.empty()) MPI_Waitall((int)rq.size(), rq.data(), MPI_STATUSES_IGNORE);
    }
    timer.lap("key-exchange");

    // ---- match each peer's received keys -> atomicMin owner (same kernel as all-peer) ----
    if (nShared > 0) {
        const uint64_t* sIL  = thrust::raw_pointer_cast(sortedIL.data());
        const uint64_t* sH   = thrust::raw_pointer_cast(sHash.data());
        const int*      sLoc = thrust::raw_pointer_cast(sLocal.data());
        for (int i = 0; i < np; ++i) {
            if (rc[i] == 0) continue;
            const uint64_t* base = d_recv + 6L * roff[i];
            int nr = rc[i];
            int b = 256, g = (nr + b - 1) / b;
            hoResolveMatchKernel<<<g, b>>>(base, nr, peers[i], sIL, sH, sLoc, (int)nShared, d_dofOwner);
        }
        cudaDeviceSynchronize();
    }
    timer.lap("match");
}

// VERIFIER driver (env MARS_HO_VERIFY_OWNERSHIP). Runs BOTH the targeted resolve and the
// all-peer ground truth on independent clones of d_dofOwner and asserts they are bit-
// identical, MPI_Abort on mismatch. This is the proof the targeting never drops a holder.
// Returns the verified targeted result in d_dofOwner. The all-peer path is the reference.
inline void resolveHoDofOwnershipGpuVerify(
        long numDof,
        const long* d_dofKind, const long* d_dofG0, const long* d_dofG1,
        const long* d_dofG2,   const long* d_dofG3, const int* d_dofPos,
        const int*  d_dofShared, int* d_dofOwner,
        int myRank, const std::vector<int>& peers,
        const long* d_cornerGid, const int* d_cornerOwner, int nCornerLocal,
        const int* d_recvNodeIds, const std::vector<int>& recvOffsets,
        const int* d_sendNodeIds, const std::vector<int>& sendOffsets)
{
    // clone the provisional owner for the reference (all-peer) run
    thrust::device_vector<int> d_ref(numDof);
    thrust::copy(thrust::device_pointer_cast(d_dofOwner),
                 thrust::device_pointer_cast(d_dofOwner) + numDof, d_ref.begin());

    // targeted into d_dofOwner (in place), all-peer into the clone
    resolveHoDofOwnershipGpuTargeted(numDof, d_dofKind, d_dofG0, d_dofG1, d_dofG2, d_dofG3,
                                     d_dofPos, d_dofShared, d_dofOwner, myRank, peers,
                                     d_cornerGid, d_cornerOwner, nCornerLocal,
                                     d_recvNodeIds, recvOffsets, d_sendNodeIds, sendOffsets);
    resolveHoDofOwnershipGpu(numDof, d_dofKind, d_dofG0, d_dofG1, d_dofG2, d_dofG3,
                             d_dofPos, d_dofShared, thrust::raw_pointer_cast(d_ref.data()),
                             myRank, peers);

    long localMismatch = thrust::inner_product(
        thrust::device_pointer_cast(d_dofOwner), thrust::device_pointer_cast(d_dofOwner) + numDof,
        d_ref.begin(), (long)0, thrust::plus<long>(),
        [] __device__ (int a, int b) { return a != b ? 1 : 0; });

    long globalMismatch = 0;
    MPI_Allreduce(&localMismatch, &globalMismatch, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    if (globalMismatch != 0) {
        if (myRank == 0)
            fprintf(stderr, "[MARS_HO_VERIFY_OWNERSHIP] FAIL: targeted vs all-peer dofOwner "
                            "mismatch (%ld DOF globally) -- targeting dropped a holder. Aborting.\n",
                    globalMismatch);
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 73);
    }
    if (myRank == 0) {
        printf("[MARS_HO_VERIFY_OWNERSHIP] PASS: targeted == all-peer (bit-identical dofOwner)\n");
        fflush(stdout);
    }
}

// ============================================================================
// DEVICE HoHalo::build (HoHalo::buildDevice). Same receiver-driven contract as the
// host build(): ghosts (dofOwner != me) grouped by owner -> recv lists + request keys;
// GPUDirect exchange of counts+keys; each peer matches my requests against ITS local
// boundary-DOF key table and replies with the matched local DOF in MY request order.
// The result (sendDof_/recvDof_ + offsets) is bit-identical to host build() and stays
// DEVICE-resident so forwardDevice/reverseAddDevice consume it with no re-upload.
//
// Reuses the resolve's exact patterns: interleaved 6-word keys (one cache line/probe),
// a 64-bit hash sort (ONE radix pass) over the local boundary keys, and a binary-search-
// with-full-compare match. The only structural difference vs the resolve is WHAT the
// match emits: not an atomicMin onto an owner, but the matched LOCAL DOF id written at
// the request's slot -> the send list, in the requester's recv order.
// ============================================================================

// Pack one BOUNDARY DOF's key (kind,g0..g3,pos) INTERLEAVED at [6*s..6*s+5] + its hash and
// local dof id at slot s. dofBoundary is the host build()'s gate (only boundary DOF are
// keyed), so this builds the SAME (key->local) table the host sorts, just on device.
__global__ void hoHaloPackBoundaryKernel(const long* __restrict__ kind, const long* __restrict__ g0,
                                         const long* __restrict__ g1,  const long* __restrict__ g2,
                                         const long* __restrict__ g3,  const int*  __restrict__ pos,
                                         const uint8_t* __restrict__ dofBoundary, const long* __restrict__ slot,
                                         long numDof,
                                         uint64_t* __restrict__ keyIL,      // interleaved, stride 6
                                         uint64_t* __restrict__ keyHash,
                                         int* __restrict__ localDof)
{
    long d = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (d >= numDof || !dofBoundary[d]) return;
    long s = slot[d];                          // exclusive-scan rank among boundary DOF
    uint64_t k0 = (uint64_t)kind[d], k1 = (uint64_t)g0[d], k2 = (uint64_t)g1[d];
    uint64_t k3 = (uint64_t)g2[d],   k4 = (uint64_t)g3[d], k5 = (uint64_t)pos[d];
    uint64_t* o = keyIL + 6L * s;
    o[0] = k0; o[1] = k1; o[2] = k2; o[3] = k3; o[4] = k4; o[5] = k5;
    keyHash[s]  = hoKeyHash(k0, k1, k2, k3, k4, k5);
    localDof[s] = (int)d;
}

// Per boundary DOF, decide if it is a GHOST to be requested from a peer (dofOwner != me,
// owner >= 0, owner is in the peer list) and, if so, emit (peerSlot=ownerPeerIndex) else -1.
// peerOf maps a rank -> its index in the peer list (or -1); built host-side, tiny. This is
// the device twin of the host loop that fills recvLocal[peerIdx]/reqKeys[peerIdx].
__global__ void hoHaloGhostPeerKernel(const long* __restrict__ slot, const uint8_t* __restrict__ dofBoundary,
                                      const int* __restrict__ dofOwner, long numDof, int myRank,
                                      const int* __restrict__ rankToPeer, int rankToPeerN,
                                      int* __restrict__ ghostPeer)   // [nBoundary], -1 if not a ghost
{
    long d = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (d >= numDof || !dofBoundary[d]) return;
    long s = slot[d];
    int o = dofOwner[d];
    int pe = -1;
    if (o != myRank && o >= 0 && o < rankToPeerN) pe = rankToPeer[o];
    ghostPeer[s] = pe;
}

// Pack the interleaved 6-word key for a list of LOCAL DOF ids, straight from the per-DOF
// columns. Used to assemble the request-key send buffer in recv order (out[6*i..] = key of
// localDof[i]) -- the device twin of host build() copying reqKeys in recvLocal order.
__global__ void hoHaloKeyFromLocalKernel(const int* __restrict__ localDof, int n,
                                         const long* __restrict__ kind, const long* __restrict__ g0,
                                         const long* __restrict__ g1,  const long* __restrict__ g2,
                                         const long* __restrict__ g3,  const int*  __restrict__ pos,
                                         uint64_t* __restrict__ out)   // interleaved, stride 6
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    long d = (long)localDof[i];
    uint64_t* o = out + 6L * i;
    o[0] = (uint64_t)kind[d]; o[1] = (uint64_t)g0[d]; o[2] = (uint64_t)g1[d];
    o[3] = (uint64_t)g2[d];   o[4] = (uint64_t)g3[d]; o[5] = (uint64_t)pos[d];
}

// Match received request keys against MY sorted boundary keys; write the matched LOCAL DOF
// at the request's slot (so the reply is in the requester's recv order). A miss writes -1
// (host build() simply skips a miss; the compaction below drops -1s). Hash binary-search +
// full 6-field compare -> a hash collision never produces a wrong match (bit-identical to a
// lexicographic search, same as the resolve's match).
__global__ void hoHaloMatchSendKernel(const uint64_t* __restrict__ recvIL,   // interleaved, stride 6
                                      int nRecv,
                                      const uint64_t* __restrict__ sortedIL,  // interleaved, stride 6
                                      const uint64_t* __restrict__ sHash,
                                      const int* __restrict__ sortedLocalDof, int nBoundary,
                                      int* __restrict__ outLocalDof)          // [nRecv], -1 on miss
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nRecv) return;
    const uint64_t* r = recvIL + 6L * i;
    uint64_t k0 = r[0], k1 = r[1], k2 = r[2], k3 = r[3], k4 = r[4], k5 = r[5];
    uint64_t h = hoKeyHash(k0, k1, k2, k3, k4, k5);
    int lo = 0, hi = nBoundary;
    while (lo < hi) { int m = (lo + hi) >> 1; if (sHash[m] < h) lo = m + 1; else hi = m; }
    int hit = -1;
    for (int j = lo; j < nBoundary && sHash[j] == h; ++j) {
        const uint64_t* s = sortedIL + 6L * j;
        if (s[0]==k0 && s[1]==k1 && s[2]==k2 && s[3]==k3 && s[4]==k4 && s[5]==k5) {
            hit = sortedLocalDof[j]; break;                  // keys unique per DOF
        }
    }
    outLocalDof[i] = hit;
}

template<typename RealType>
__global__ void hoHaloGatherKernel(const RealType* __restrict__ vec, const int* __restrict__ idx,
                                   RealType* __restrict__ buf, int n)
{ int i = blockIdx.x*blockDim.x + threadIdx.x; if (i < n) buf[i] = vec[idx[i]]; }

template<typename RealType>
__global__ void hoHaloScatterKernel(RealType* __restrict__ vec, const int* __restrict__ idx,
                                    const RealType* __restrict__ buf, int n)
{ int i = blockIdx.x*blockDim.x + threadIdx.x; if (i < n) vec[idx[i]] = buf[i]; }

template<typename RealType>
__global__ void hoHaloScatterAddKernel(RealType* __restrict__ vec, const int* __restrict__ idx,
                                       const RealType* __restrict__ buf, int n)
{ int i = blockIdx.x*blockDim.x + threadIdx.x; if (i < n) atomicAdd(&vec[idx[i]], buf[i]); }
#endif

template<typename RealType>
class HoHalo {
public:
    using DofKey = HODofHandler::DofKey;

    std::vector<int> peers_;
    std::vector<int> sendOffsets_{0}, recvOffsets_{0};  // CSR per peer
    std::vector<int> sendDof_, recvDof_;                // local DOF indices

    static std::array<long,6> packKey(const DofKey& k)
    { return { (long)k.kind, k.g0, k.g1, k.g2, k.g3, (long)k.pos }; }

    void build(int numDof,
               const std::vector<int>&     dofOwner,
               const std::vector<DofKey>&  dofKey,
               const std::vector<uint8_t>& dofBoundary,
               int myRank,
               const std::vector<int>&     candidatePeers)
    {
        const int np = (int)candidatePeers.size();
        std::map<int,int> peerIdx;
        for (int i = 0; i < np; ++i) peerIdx[candidatePeers[i]] = i;

        // Only BOUNDARY DOF are ever exchanged, so key only those. Keying all numDof
        // builds a numDof-sized std::map (~78 GB at 650M DOF) that OOMs the host at scale.
        // sorted (key -> local), NOT std::map: cache-friendly binary search replaces the
        // millions of tree ops that made this the dominant setup cost at scale.
        std::vector<std::pair<std::array<long,6>, int>> keyToLocal;
        for (int d = 0; d < numDof; ++d) if (dofBoundary[d]) keyToLocal.push_back({packKey(dofKey[d]), d});
        std::sort(keyToLocal.begin(), keyToLocal.end(),
                  [](const auto& a, const auto& b){ return a.first < b.first; });

        // Ghost DOF (owner != me) grouped by owner -> my recv lists + the keys I request.
        // Ghosts are always boundary, so the same gate keeps this O(surface), not O(numDof).
        std::vector<std::vector<int>>            recvLocal(np);
        std::vector<std::vector<std::array<long,6>>> reqKeys(np);
        for (int d = 0; d < numDof; ++d) {
            if (!dofBoundary[d]) continue;
            int o = dofOwner[d];
            if (o == myRank || o < 0) continue;          // owned or interior
            auto it = peerIdx.find(o);
            if (it == peerIdx.end()) continue;           // owner not a neighbour (should not happen)
            recvLocal[it->second].push_back(d);
            reqKeys[it->second].push_back(packKey(dofKey[d]));
        }

        // Exchange request counts (symmetric peers -> no deadlock).
        std::vector<int> reqSendCnt(np), reqRecvCnt(np, 0);
        for (int i = 0; i < np; ++i) reqSendCnt[i] = (int)reqKeys[i].size();
        exchangeCounts(candidatePeers, reqSendCnt, reqRecvCnt, 0x4849);

        // Exchange request keys (6 longs each). Keys received from peer i = the DOF
        // peer i wants from me, in peer i's recv-slot order.
        std::vector<std::vector<std::array<long,6>>> gotKeys(np);
        {
            std::vector<MPI_Request> rq; rq.reserve(2 * np);
            for (int i = 0; i < np; ++i) {
                gotKeys[i].resize(reqRecvCnt[i]);
                if (reqRecvCnt[i] > 0) {
                    MPI_Request r;
                    MPI_Irecv(gotKeys[i].data(), reqRecvCnt[i] * 6, MPI_LONG,
                              candidatePeers[i], 0x484a, MPI_COMM_WORLD, &r);
                    rq.push_back(r);
                }
            }
            for (int i = 0; i < np; ++i) {
                if (reqSendCnt[i] > 0) {
                    MPI_Request r;
                    MPI_Isend(reqKeys[i].data(), reqSendCnt[i] * 6, MPI_LONG,
                              candidatePeers[i], 0x484a, MPI_COMM_WORLD, &r);
                    rq.push_back(r);
                }
            }
            if (!rq.empty()) MPI_Waitall((int)rq.size(), rq.data(), MPI_STATUSES_IGNORE);
        }

        // My send list to peer i = its requested keys mapped to my local DOF, IN THE
        // RECEIVED ORDER (matches peer i's recv slot order) -> forward/reverse align.
        std::vector<std::vector<int>> sendLocal(np);
        for (int i = 0; i < np; ++i) {
            sendLocal[i].reserve(gotKeys[i].size());
            for (auto& k : gotKeys[i]) {
                auto it = std::lower_bound(keyToLocal.begin(), keyToLocal.end(), k,
                            [](const auto& a, const std::array<long,6>& key){ return a.first < key; });
                if (it != keyToLocal.end() && it->first == k) sendLocal[i].push_back(it->second);
            }
        }

        // Compact into CSR over peers with any traffic.
        peers_.clear(); sendOffsets_.assign(1, 0); recvOffsets_.assign(1, 0);
        sendDof_.clear(); recvDof_.clear();
        for (int i = 0; i < np; ++i) {
            if (sendLocal[i].empty() && recvLocal[i].empty()) continue;
            peers_.push_back(candidatePeers[i]);
            for (int v : sendLocal[i]) sendDof_.push_back(v);
            for (int v : recvLocal[i]) recvDof_.push_back(v);
            sendOffsets_.push_back((int)sendDof_.size());
            recvOffsets_.push_back((int)recvDof_.size());
        }
    }

    // FORWARD: owner -> ghost. Ghost DOF slots are OVERWRITTEN with the owner's value.
    void forward(std::vector<RealType>& vec) const
    {
        const int ns = sendOffsets_.back(), nr = recvOffsets_.back();
        std::vector<RealType> sbuf(ns), rbuf(nr);
        for (int i = 0; i < ns; ++i) sbuf[i] = vec[sendDof_[i]];
        exchangeVals(sbuf, rbuf, sendOffsets_, recvOffsets_, 0x484b);
        for (int i = 0; i < nr; ++i) vec[recvDof_[i]] = rbuf[i];
    }

    // REVERSE-ADD: ghost -> owner. Ghost contributions are SUMMED into owner slots.
    void reverseAdd(std::vector<RealType>& vec) const
    {
        const int ns = recvOffsets_.back(), nr = sendOffsets_.back();
        std::vector<RealType> sbuf(ns), rbuf(nr);
        for (int i = 0; i < ns; ++i) sbuf[i] = vec[recvDof_[i]];
        exchangeVals(sbuf, rbuf, recvOffsets_, sendOffsets_, 0x484c);
        for (int i = 0; i < nr; ++i) vec[sendDof_[i]] += rbuf[i];
    }

    // ---- device path (for scaling): gather/scatter + CUDA-aware MPI ----
    int*      d_sendDof_ = nullptr;
    int*      d_recvDof_ = nullptr;
    RealType* d_sendBuf_ = nullptr;
    RealType* d_recvBuf_ = nullptr;
    int       nSend_ = 0, nRecv_ = 0;

    // In-flight requests for the SPLIT (non-blocking) forward used by the overlap matvec.
    // forwardDeviceStart fills these; forwardDeviceFinish waits on them + scatters. Kept
    // here (not local) so the interior apply can run between the two calls.
    mutable std::vector<MPI_Request> fwdReq_;
    mutable std::vector<MPI_Request> revReq_;

#ifdef __CUDACC__
    void uploadDevice()
    {
        // buildDevice already leaves the lists device-resident -> nothing to upload.
        // (Re-uploading would alloc from the host mirror and leak the device lists.)
        if (d_sendDof_ != nullptr) return;
        nSend_ = sendOffsets_.back();
        nRecv_ = recvOffsets_.back();
        cudaMalloc(&d_sendDof_, sizeof(int)      * (nSend_ > 0 ? nSend_ : 1));
        cudaMalloc(&d_recvDof_, sizeof(int)      * (nRecv_ > 0 ? nRecv_ : 1));
        cudaMalloc(&d_sendBuf_, sizeof(RealType) * (nSend_ > 0 ? nSend_ : 1));
        cudaMalloc(&d_recvBuf_, sizeof(RealType) * (nRecv_ > 0 ? nRecv_ : 1));
        if (nSend_) cudaMemcpy(d_sendDof_, sendDof_.data(), sizeof(int)*nSend_, cudaMemcpyHostToDevice);
        if (nRecv_) cudaMemcpy(d_recvDof_, recvDof_.data(), sizeof(int)*nRecv_, cudaMemcpyHostToDevice);
    }
    void freeDevice()
    {
        cudaFree(d_sendDof_); cudaFree(d_recvDof_); cudaFree(d_sendBuf_); cudaFree(d_recvBuf_);
        d_sendDof_ = d_recvDof_ = nullptr; d_sendBuf_ = d_recvBuf_ = nullptr;
    }

    // DEVICE build of the per-peer send/recv lists, bit-identical to host build() but with
    // NO host round-trip on the per-DOF arrays. Inputs are the device columns the GPU
    // numbering already produced (d_dofKind/G0..G3/pos + the RESOLVED d_dofOwner + the
    // boundary mask d_dofBoundary). On return d_sendDof_/d_recvDof_ + nSend_/nRecv_ + the
    // GPUDirect buffers are READY -- forwardDevice/reverseAddDevice run with no uploadDevice.
    //
    // BIT-IDENTICAL TO HOST build() (same DOF in send/recv, same order, A.send[B]==B.recv[A]):
    //  * recv list per peer = boundary GHOST DOF (dofOwner != me) of that peer in ASCENDING
    //    local-id order. Host: iterates d in ascending order, appends to recvLocal[peer].
    //    Device: emit (peer, localDof) for each ghost, then a STABLE sort by peer -> the
    //    secondary (local-id) order is preserved == host order.
    //  * request keys to a peer = the keys of MY recv DOF in MY recv order (gathered from
    //    the recv list). Host sends reqKeys[i] in the same order as recvLocal[i].
    //  * send list per peer = peer's received request keys mapped to MY local DOF IN THE
    //    RECEIVED ORDER. The match writes the matched local DOF at the request's slot, so
    //    received order is preserved exactly (== the peer's recv order -> forward/reverse
    //    align). A miss writes -1 and is dropped, exactly like host's skipped lower_bound.
    //  * peers_/offsets = only peers with traffic, in candidatePeers order. Host and device
    //    both walk candidatePeers in order and skip empty peers.
    // The match is hash-binary-search + full 6-field compare -> never a false match; the
    // only freedom vs host (std::lower_bound) is collision-handling, which the full compare
    // removes. Keys are unique per DOF, so the matched local DOF is unique too.
    void buildDevice(long numDof,
                     const long* d_dofKind, const long* d_dofG0, const long* d_dofG1,
                     const long* d_dofG2,   const long* d_dofG3, const int* d_dofPos,
                     const int*  d_dofOwner, const uint8_t* d_dofBoundary,
                     int myRank, const std::vector<int>& candidatePeers)
    {
        const int np = (int)candidatePeers.size();
        const bool tOn = std::getenv("MARS_HO_HALO_TIMING") != nullptr;
        double t0 = MPI_Wtime();
        auto lap = [&](const char* stage) {
            if (!tOn) { t0 = MPI_Wtime(); return; }
            cudaDeviceSynchronize();
            double t1 = MPI_Wtime();
            if (myRank == 0) { printf("[ho-halo] %-22s %8.3f ms\n", stage, (t1 - t0) * 1e3); fflush(stdout); }
            t0 = t1;
        };

        // ---- count boundary DOF + assign each a contiguous slot (exclusive scan) ----
        thrust::device_vector<long> d_slot(numDof > 0 ? numDof : 1);
        if (numDof > 0) {
            thrust::device_ptr<const uint8_t> bd(d_dofBoundary);
            // exclusive_scan over a uint8 lane with a long accumulator -> slot per boundary DOF.
            thrust::exclusive_scan(bd, bd + numDof, d_slot.begin(), (long)0);
        }
        long nBoundary = 0;
        if (numDof > 0) {
            uint8_t lastB = 0; long lastSlot = 0;
            thrust::copy(thrust::device_pointer_cast(d_dofBoundary) + (numDof - 1),
                         thrust::device_pointer_cast(d_dofBoundary) + numDof, &lastB);
            thrust::copy(d_slot.begin() + (numDof - 1), d_slot.begin() + numDof, &lastSlot);
            nBoundary = lastSlot + (lastB ? 1 : 0);
        }
        lap("scan-slots");

        // ---- pack + hash-sort MY boundary keys into one interleaved (stride-6) array ----
        // This is the device twin of host build()'s sorted (key->local) table.
        thrust::device_vector<uint64_t> packIL(6L * (nBoundary > 0 ? nBoundary : 1));
        thrust::device_vector<uint64_t> packHash(nBoundary > 0 ? nBoundary : 1);
        thrust::device_vector<int>      packLocal(nBoundary > 0 ? nBoundary : 1);
        // ghostPeer[s] = peer index of boundary DOF s if it is a ghost we request, else -1.
        thrust::device_vector<int>      ghostPeer(nBoundary > 0 ? nBoundary : 1, -1);
        if (nBoundary > 0) {
            int b = 256; long g = (numDof + b - 1) / b;
            hoHaloPackBoundaryKernel<<<g, b>>>(d_dofKind, d_dofG0, d_dofG1, d_dofG2, d_dofG3, d_dofPos,
                                               d_dofBoundary, thrust::raw_pointer_cast(d_slot.data()), numDof,
                                               thrust::raw_pointer_cast(packIL.data()),
                                               thrust::raw_pointer_cast(packHash.data()),
                                               thrust::raw_pointer_cast(packLocal.data()));
            // rank -> peer index (tiny host table); ghosts read it to find their owner's peer.
            int maxRank = myRank;
            for (int p : candidatePeers) maxRank = std::max(maxRank, p);
            std::vector<int> h_rankToPeer(maxRank + 1, -1);
            for (int i = 0; i < np; ++i) h_rankToPeer[candidatePeers[i]] = i;
            thrust::device_vector<int> d_rankToPeer(h_rankToPeer.begin(), h_rankToPeer.end());
            hoHaloGhostPeerKernel<<<g, b>>>(thrust::raw_pointer_cast(d_slot.data()), d_dofBoundary,
                                            d_dofOwner, numDof, myRank,
                                            thrust::raw_pointer_cast(d_rankToPeer.data()), (int)d_rankToPeer.size(),
                                            thrust::raw_pointer_cast(ghostPeer.data()));
        }
        d_slot = thrust::device_vector<long>();   // numDof-sized scratch no longer needed

        // sort a permutation by hash (one radix pass), then gather keys+hash+local by it.
        thrust::device_vector<uint64_t> sortedIL(6L * (nBoundary > 0 ? nBoundary : 1));
        thrust::device_vector<uint64_t> sHash(nBoundary > 0 ? nBoundary : 1);
        thrust::device_vector<int>      sLocal(nBoundary > 0 ? nBoundary : 1);
        if (nBoundary > 0) {
            thrust::device_vector<int> idx(nBoundary);
            thrust::sequence(idx.begin(), idx.end());
            thrust::stable_sort_by_key(packHash.begin(), packHash.end(), idx.begin());
            sHash.swap(packHash);
            thrust::gather(idx.begin(), idx.end(), packLocal.begin(), sLocal.begin());
            int b = 256, g = (int)((nBoundary + b - 1) / b);
            hoKeyGatherInterleaveKernel<<<g, b>>>(thrust::raw_pointer_cast(packIL.data()),
                                                  thrust::raw_pointer_cast(idx.data()), (int)nBoundary,
                                                  thrust::raw_pointer_cast(sortedIL.data()));
        }
        lap("pack+hash-sort");

        // ---- build MY recv lists: ghosts grouped by peer, ascending local-id within a peer ----
        // Stable sort (peerIdx, localDof) by peerIdx -> per-peer recv DOF in ascending local
        // id (the host order). The request KEYS sent to a peer are the keys of these recv DOF
        // in this order, so the peer's reply (matched send list) lands in this same order.
        std::vector<int> recvCnt(np, 0);
        thrust::device_vector<int>      d_recvDofSorted(nBoundary > 0 ? nBoundary : 1);
        thrust::device_vector<uint64_t> reqKeysIL(6L * (nBoundary > 0 ? nBoundary : 1));
        int nGhost = 0;
        if (nBoundary > 0) {
            // compact (peer, localDof) for ghosts (peer >= 0). packLocal is the boundary
            // local id in slot order; ghostPeer is the peer index (or -1).
            nGhost = (int)thrust::count_if(ghostPeer.begin(), ghostPeer.end(),
                                           [] __device__ (int pe) { return pe >= 0; });
            thrust::device_vector<int> gPeer(nGhost > 0 ? nGhost : 1), gLocal(nGhost > 0 ? nGhost : 1);
            if (nGhost > 0) {
                thrust::copy_if(ghostPeer.begin(), ghostPeer.end(), gPeer.begin(),
                                [] __device__ (int pe) { return pe >= 0; });
                thrust::copy_if(packLocal.begin(), packLocal.end(), ghostPeer.begin(), gLocal.begin(),
                                [] __device__ (int pe) { return pe >= 0; });
                // STABLE sort by peer -> ascending-local-id within each peer is preserved
                // (packLocal is already ascending in slot order == boundary scan order).
                thrust::stable_sort_by_key(gPeer.begin(), gPeer.end(), gLocal.begin());
                thrust::copy(gLocal.begin(), gLocal.end(), d_recvDofSorted.begin());
                // per-peer recv counts via reduce_by_key over the sorted peer ids
                thrust::device_vector<int> segPeer(np), segCnt(np);
                auto e = thrust::reduce_by_key(gPeer.begin(), gPeer.end(),
                                               thrust::make_constant_iterator(1),
                                               segPeer.begin(), segCnt.begin());
                int nSeg = (int)(e.first - segPeer.begin());
                std::vector<int> hSeg(nSeg), hCnt(nSeg);
                thrust::copy(segPeer.begin(), segPeer.begin() + nSeg, hSeg.begin());
                thrust::copy(segCnt.begin(),  segCnt.begin()  + nSeg, hCnt.begin());
                for (int s = 0; s < nSeg; ++s) recvCnt[hSeg[s]] = hCnt[s];
                // gather the request keys (interleaved) for the recv DOF, in recv order. We
                // re-key from the per-DOF columns: gLocal holds the LOCAL DOF id (== d), so
                // its key is the same one packed above. Reuse the sorted table via a search
                // would cost a probe; instead pack directly from the columns by local id.
                int b = 256, g = (nGhost + b - 1) / b;
                hoHaloKeyFromLocalKernel<<<g, b>>>(thrust::raw_pointer_cast(gLocal.data()), nGhost,
                                                   d_dofKind, d_dofG0, d_dofG1, d_dofG2, d_dofG3, d_dofPos,
                                                   thrust::raw_pointer_cast(reqKeysIL.data()));
            }
        }
        // free the boundary-key build scratch we no longer need (keep sorted* for the match)
        packIL = thrust::device_vector<uint64_t>();
        packLocal = thrust::device_vector<int>();
        ghostPeer = thrust::device_vector<int>();
        lap("recv-lists+req-keys");

        // ---- exchange request counts (symmetric peers -> no deadlock) ----
        std::vector<int> reqRecvCnt(np, 0);
        {
            std::vector<MPI_Request> rq; rq.reserve(2 * np);
            for (int i = 0; i < np; ++i) { MPI_Request r; MPI_Irecv(&reqRecvCnt[i], 1, MPI_INT, candidatePeers[i], 0x4862, MPI_COMM_WORLD, &r); rq.push_back(r); }
            for (int i = 0; i < np; ++i) { MPI_Request r; MPI_Isend(&recvCnt[i],    1, MPI_INT, candidatePeers[i], 0x4862, MPI_COMM_WORLD, &r); rq.push_back(r); }
            if (!rq.empty()) MPI_Waitall((int)rq.size(), rq.data(), MPI_STATUSES_IGNORE);
        }
        lap("count-exchange");

        // ---- GPUDirect exchange of request keys (6 uint64 each) ----
        // I SEND my request keys (reqKeysIL, per-peer slices in recv order) and RECEIVE each
        // peer's requests into one device slab; I then match the received keys against MY
        // boundary table and reply with the matched local DOF.
        std::vector<int> reqSendOff(np + 1, 0);
        for (int i = 0; i < np; ++i) reqSendOff[i + 1] = reqSendOff[i] + recvCnt[i];   // my req key send layout
        std::vector<int> gotOff(np + 1, 0);
        for (int i = 0; i < np; ++i) gotOff[i + 1] = gotOff[i] + reqRecvCnt[i];
        const int totalGot = gotOff[np];
        thrust::device_vector<uint64_t> gotIL(6L * (totalGot > 0 ? totalGot : 1));
        uint64_t* d_got = thrust::raw_pointer_cast(gotIL.data());
        const uint64_t* d_req = thrust::raw_pointer_cast(reqKeysIL.data());

        cudaDeviceSynchronize();   // reqKeysIL built by async kernels -> finish before GPUDirect reads it
        {
            std::vector<MPI_Request> rq; rq.reserve(2 * np);
            for (int i = 0; i < np; ++i)
                if (reqRecvCnt[i] > 0) { MPI_Request r; MPI_Irecv(d_got + 6L * gotOff[i], reqRecvCnt[i] * 6, MPI_UNSIGNED_LONG_LONG, candidatePeers[i], 0x4863, MPI_COMM_WORLD, &r); rq.push_back(r); }
            for (int i = 0; i < np; ++i)
                if (recvCnt[i] > 0)    { MPI_Request r; MPI_Isend((void*)(d_req + 6L * reqSendOff[i]), recvCnt[i] * 6, MPI_UNSIGNED_LONG_LONG, candidatePeers[i], 0x4863, MPI_COMM_WORLD, &r); rq.push_back(r); }
            if (!rq.empty()) MPI_Waitall((int)rq.size(), rq.data(), MPI_STATUSES_IGNORE);
        }
        lap("key-exchange");

        // ---- match received requests -> matched local DOF in received order (= send list) ----
        thrust::device_vector<int> d_sendMatch(totalGot > 0 ? totalGot : 1, -1);
        if (totalGot > 0 && nBoundary > 0) {
            int b = 256, g = (totalGot + b - 1) / b;
            hoHaloMatchSendKernel<<<g, b>>>(d_got, totalGot,
                                            thrust::raw_pointer_cast(sortedIL.data()),
                                            thrust::raw_pointer_cast(sHash.data()),
                                            thrust::raw_pointer_cast(sLocal.data()), (int)nBoundary,
                                            thrust::raw_pointer_cast(d_sendMatch.data()));
            cudaDeviceSynchronize();
        }
        lap("match");

        // ---- compact into device CSR over peers with traffic (host offsets, device DOF) ----
        // sendDof for peer i = d_sendMatch[gotOff[i]..) with -1 misses dropped. recvDof for
        // peer i = d_recvDofSorted[recv slice]. Build the host offsets first (skip empty
        // peers, candidatePeers order), then pack the two device DOF arrays contiguously.
        peers_.clear(); sendOffsets_.assign(1, 0); recvOffsets_.assign(1, 0);

        // per-peer matched-send counts (drop -1 misses). Host build() drops a missed
        // lower_bound the same way; a miss is a key the peer requested but I don't hold,
        // which for a correct receiver-driven contract should not happen, but we mirror
        // host's tolerant skip rather than assume.
        std::vector<int> sendCnt(np, 0);
        std::vector<int> recvStart(np, 0);   // start of peer i's recv slice in d_recvDofSorted
        {
            int acc = 0;
            for (int i = 0; i < np; ++i) { recvStart[i] = acc; acc += recvCnt[i]; }
        }
        // count matches per peer on host (totalGot is the rank surface, tiny -> cheap copy).
        std::vector<int> h_match(totalGot > 0 ? totalGot : 1, -1);
        if (totalGot > 0) thrust::copy(d_sendMatch.begin(), d_sendMatch.begin() + totalGot, h_match.begin());
        for (int i = 0; i < np; ++i) {
            int c = 0;
            for (int s = gotOff[i]; s < gotOff[i + 1]; ++s) if (h_match[s] >= 0) ++c;
            sendCnt[i] = c;
        }

        // first pass: peers + offsets (only peers with any traffic)
        std::vector<int> keptPeer;                 // index into candidatePeers
        for (int i = 0; i < np; ++i) {
            if (sendCnt[i] == 0 && recvCnt[i] == 0) continue;
            peers_.push_back(candidatePeers[i]);
            keptPeer.push_back(i);
            sendOffsets_.push_back(sendOffsets_.back() + sendCnt[i]);
            recvOffsets_.push_back(recvOffsets_.back() + recvCnt[i]);
        }
        nSend_ = sendOffsets_.back();
        nRecv_ = recvOffsets_.back();
        lap("compact-offsets");

        // ---- pack the device DOF arrays into the final d_sendDof_/d_recvDof_ ----
        // Drop the -1 send misses with a per-peer copy_if into the peer's send slice; copy
        // the recv slice straight (no misses there). The arrays land in peers_ order.
        cudaMalloc(&d_sendDof_, sizeof(int) * (nSend_ > 0 ? nSend_ : 1));
        cudaMalloc(&d_recvDof_, sizeof(int) * (nRecv_ > 0 ? nRecv_ : 1));
        cudaMalloc(&d_sendBuf_, sizeof(RealType) * (nSend_ > 0 ? nSend_ : 1));
        cudaMalloc(&d_recvBuf_, sizeof(RealType) * (nRecv_ > 0 ? nRecv_ : 1));
        {
            thrust::device_ptr<int> sd(d_sendDof_), rd(d_recvDof_);
            int sOff = 0, rOff = 0;
            for (size_t kp = 0; kp < keptPeer.size(); ++kp) {
                int i = keptPeer[kp];
                // send: drop -1 misses from this peer's received-match slice
                if (gotOff[i + 1] > gotOff[i]) {
                    auto e = thrust::copy_if(d_sendMatch.begin() + gotOff[i], d_sendMatch.begin() + gotOff[i + 1],
                                             sd + sOff, [] __device__ (int v) { return v >= 0; });
                    sOff += (int)(e - (sd + sOff));
                }
                // recv: contiguous slice, no misses
                if (recvCnt[i] > 0) {
                    thrust::copy(d_recvDofSorted.begin() + recvStart[i],
                                 d_recvDofSorted.begin() + recvStart[i] + recvCnt[i], rd + rOff);
                    rOff += recvCnt[i];
                }
            }
        }
        cudaDeviceSynchronize();
        lap("pack-device-csr");

        // Mirror the final lists to the host sendDof_/recvDof_ so the HOST forward()/
        // reverseAdd() (the small-run A.1 oracle gate) still works. This is the rank
        // SURFACE only (nSend_+nRecv_, not numDof) -> a tiny D2H, not a per-DOF round-trip.
        // The device path proper never reads these.
        sendDof_.resize(nSend_);
        recvDof_.resize(nRecv_);
        if (nSend_) cudaMemcpy(sendDof_.data(), d_sendDof_, sizeof(int) * nSend_, cudaMemcpyDeviceToHost);
        if (nRecv_) cudaMemcpy(recvDof_.data(), d_recvDof_, sizeof(int) * nRecv_, cudaMemcpyDeviceToHost);
        lap("host-mirror");
    }

    // FORWARD on device: owner -> ghost (overwrite; recv DOF are unique). GPUDirect MPI.
    void forwardDevice(RealType* d_vec, cudaStream_t stream = 0) const
    {
        auto mpiT = (sizeof(RealType) == 8) ? MPI_DOUBLE : MPI_FLOAT;
        if (nSend_) { int b=256, g=(nSend_+b-1)/b; hoHaloGatherKernel<RealType><<<g,b,0,stream>>>(d_vec, d_sendDof_, d_sendBuf_, nSend_); }
        cudaStreamSynchronize(stream);
        std::vector<MPI_Request> rq; rq.reserve(2 * peers_.size());
        for (size_t p=0;p<peers_.size();++p){ int c=recvOffsets_[p+1]-recvOffsets_[p]; if(c){ MPI_Request r; MPI_Irecv(d_recvBuf_+recvOffsets_[p], c, mpiT, peers_[p], 0x4860, MPI_COMM_WORLD, &r); rq.push_back(r);} }
        for (size_t p=0;p<peers_.size();++p){ int c=sendOffsets_[p+1]-sendOffsets_[p]; if(c){ MPI_Request r; MPI_Isend(d_sendBuf_+sendOffsets_[p], c, mpiT, peers_[p], 0x4860, MPI_COMM_WORLD, &r); rq.push_back(r);} }
        if(!rq.empty()) MPI_Waitall((int)rq.size(), rq.data(), MPI_STATUSES_IGNORE);
        if (nRecv_) { int b=256, g=(nRecv_+b-1)/b; hoHaloScatterKernel<RealType><<<g,b,0,stream>>>(d_vec, d_recvDof_, d_recvBuf_, nRecv_); }
    }

    // SPLIT forward for the OVERLAP matvec: START posts the gather + the GPUDirect
    // Irecv/Isend but does NOT wait, so the caller can run interior compute while the
    // halo is in flight. This is exactly forwardDevice's first half: the gather kernel
    // is queued on `stream`, then cudaStreamSynchronize guarantees d_sendBuf_ is filled
    // before MPI reads it (GPUDirect reads device memory directly), then the non-blocking
    // sends/recvs are posted into fwdReq_. FINISH must be called before the recv DOF
    // (ghost slots) are read. Same buffers/tags/order as forwardDevice -> bit-identical.
    void forwardDeviceStart(RealType* d_vec, cudaStream_t stream = 0) const
    {
        auto mpiT = (sizeof(RealType) == 8) ? MPI_DOUBLE : MPI_FLOAT;
        if (nSend_) { int b=256, g=(nSend_+b-1)/b; hoHaloGatherKernel<RealType><<<g,b,0,stream>>>(d_vec, d_sendDof_, d_sendBuf_, nSend_); }
        cudaStreamSynchronize(stream);
        fwdReq_.clear(); fwdReq_.reserve(2 * peers_.size());
        for (size_t p=0;p<peers_.size();++p){ int c=recvOffsets_[p+1]-recvOffsets_[p]; if(c){ MPI_Request r; MPI_Irecv(d_recvBuf_+recvOffsets_[p], c, mpiT, peers_[p], 0x4860, MPI_COMM_WORLD, &r); fwdReq_.push_back(r);} }
        for (size_t p=0;p<peers_.size();++p){ int c=sendOffsets_[p+1]-sendOffsets_[p]; if(c){ MPI_Request r; MPI_Isend(d_sendBuf_+sendOffsets_[p], c, mpiT, peers_[p], 0x4860, MPI_COMM_WORLD, &r); fwdReq_.push_back(r);} }
    }

    // FINISH the split forward: wait on the in-flight requests, then scatter the received
    // ghost values into d_vec. After this returns, every ghost (recv DOF) slot of d_vec
    // holds its owner's value -- identical post-state to the blocking forwardDevice.
    void forwardDeviceFinish(RealType* d_vec, cudaStream_t stream = 0) const
    {
        if (!fwdReq_.empty()) MPI_Waitall((int)fwdReq_.size(), fwdReq_.data(), MPI_STATUSES_IGNORE);
        fwdReq_.clear();
        if (nRecv_) { int b=256, g=(nRecv_+b-1)/b; hoHaloScatterKernel<RealType><<<g,b,0,stream>>>(d_vec, d_recvDof_, d_recvBuf_, nRecv_); }
    }

    // REVERSE-ADD on device: ghost -> owner (atomic; an owned DOF may go to many peers).
    void reverseAddDevice(RealType* d_vec, cudaStream_t stream = 0) const
    {
        auto mpiT = (sizeof(RealType) == 8) ? MPI_DOUBLE : MPI_FLOAT;
        if (nRecv_) { int b=256, g=(nRecv_+b-1)/b; hoHaloGatherKernel<RealType><<<g,b,0,stream>>>(d_vec, d_recvDof_, d_recvBuf_, nRecv_); }
        cudaStreamSynchronize(stream);
        std::vector<MPI_Request> rq; rq.reserve(2 * peers_.size());
        for (size_t p=0;p<peers_.size();++p){ int c=sendOffsets_[p+1]-sendOffsets_[p]; if(c){ MPI_Request r; MPI_Irecv(d_sendBuf_+sendOffsets_[p], c, mpiT, peers_[p], 0x4861, MPI_COMM_WORLD, &r); rq.push_back(r);} }
        for (size_t p=0;p<peers_.size();++p){ int c=recvOffsets_[p+1]-recvOffsets_[p]; if(c){ MPI_Request r; MPI_Isend(d_recvBuf_+recvOffsets_[p], c, mpiT, peers_[p], 0x4861, MPI_COMM_WORLD, &r); rq.push_back(r);} }
        if(!rq.empty()) MPI_Waitall((int)rq.size(), rq.data(), MPI_STATUSES_IGNORE);
        if (nSend_) { int b=256, g=(nSend_+b-1)/b; hoHaloScatterAddKernel<RealType><<<g,b,0,stream>>>(d_vec, d_sendDof_, d_sendBuf_, nSend_); }
    }

    // SPLIT reverse-add for the OVERLAP matvec: START gathers the ghost contributions and posts
    // the GPUDirect Irecv/Isend but does NOT wait, so the caller can run independent owned-only
    // compute while the reverse halo is in flight. The gathered ghost slots are written ONLY by
    // boundary elements, so START must follow the boundary apply. Same buffers/tags/order as
    // reverseAddDevice -> bit-identical result. FINISH must run before d_vec's owned DOF are read.
    void reverseAddDeviceStart(RealType* d_vec, cudaStream_t stream = 0) const
    {
        auto mpiT = (sizeof(RealType) == 8) ? MPI_DOUBLE : MPI_FLOAT;
        if (nRecv_) { int b=256, g=(nRecv_+b-1)/b; hoHaloGatherKernel<RealType><<<g,b,0,stream>>>(d_vec, d_recvDof_, d_recvBuf_, nRecv_); }
        cudaStreamSynchronize(stream);
        revReq_.clear(); revReq_.reserve(2 * peers_.size());
        for (size_t p=0;p<peers_.size();++p){ int c=sendOffsets_[p+1]-sendOffsets_[p]; if(c){ MPI_Request r; MPI_Irecv(d_sendBuf_+sendOffsets_[p], c, mpiT, peers_[p], 0x4861, MPI_COMM_WORLD, &r); revReq_.push_back(r);} }
        for (size_t p=0;p<peers_.size();++p){ int c=recvOffsets_[p+1]-recvOffsets_[p]; if(c){ MPI_Request r; MPI_Isend(d_recvBuf_+recvOffsets_[p], c, mpiT, peers_[p], 0x4861, MPI_COMM_WORLD, &r); revReq_.push_back(r);} }
    }

    // FINISH the split reverse-add: wait on the in-flight requests, then atomically add the
    // received contributions into d_vec's owned (send) DOF. Identical post-state to reverseAddDevice.
    void reverseAddDeviceFinish(RealType* d_vec, cudaStream_t stream = 0) const
    {
        if (!revReq_.empty()) MPI_Waitall((int)revReq_.size(), revReq_.data(), MPI_STATUSES_IGNORE);
        revReq_.clear();
        if (nSend_) { int b=256, g=(nSend_+b-1)/b; hoHaloScatterAddKernel<RealType><<<g,b,0,stream>>>(d_vec, d_sendDof_, d_sendBuf_, nSend_); }
    }
#endif

private:
    static void exchangeCounts(const std::vector<int>& peers,
                               const std::vector<int>& sc, std::vector<int>& rc, int tag)
    {
        const int np = (int)peers.size();
        std::vector<MPI_Request> rq; rq.reserve(2 * np);
        for (int i = 0; i < np; ++i) { MPI_Request r; MPI_Irecv(&rc[i],1,MPI_INT,peers[i],tag,MPI_COMM_WORLD,&r); rq.push_back(r); }
        for (int i = 0; i < np; ++i) { MPI_Request r; MPI_Isend(&sc[i],1,MPI_INT,peers[i],tag,MPI_COMM_WORLD,&r); rq.push_back(r); }
        if (!rq.empty()) MPI_Waitall((int)rq.size(), rq.data(), MPI_STATUSES_IGNORE);
    }

    void exchangeVals(const std::vector<RealType>& sbuf, std::vector<RealType>& rbuf,
                      const std::vector<int>& sOff, const std::vector<int>& rOff, int tag) const
    {
        auto mpiT = (sizeof(RealType) == 8) ? MPI_DOUBLE : MPI_FLOAT;
        std::vector<MPI_Request> rq; rq.reserve(2 * peers_.size());
        for (size_t p = 0; p < peers_.size(); ++p) {
            int rcnt = rOff[p+1] - rOff[p];
            if (rcnt > 0) { MPI_Request r; MPI_Irecv((void*)(rbuf.data()+rOff[p]), rcnt, mpiT, peers_[p], tag, MPI_COMM_WORLD, &r); rq.push_back(r); }
        }
        for (size_t p = 0; p < peers_.size(); ++p) {
            int scnt = sOff[p+1] - sOff[p];
            if (scnt > 0) { MPI_Request r; MPI_Isend((void*)(sbuf.data()+sOff[p]), scnt, mpiT, peers_[p], tag, MPI_COMM_WORLD, &r); rq.push_back(r); }
        }
        if (!rq.empty()) MPI_Waitall((int)rq.size(), rq.data(), MPI_STATUSES_IGNORE);
    }
};

} // namespace fem
} // namespace mars
