#ifndef MARS_ACM_DIST_HPP
#define MARS_ACM_DIST_HPP

// Distributed ACM support: globally-consistent aggregation + per-level halo exchange + the replicated
// global coarsest. WHY: a rank-local hierarchy has no global coarse correction -- information crosses one
// halo layer per Krylov iteration, which stalls on elongated geometry (measured: 1-rank 57 iters, 2-rank
// 1980-cap at the same Stokes system). The fix is seam-consistent coarsening so ONE global hierarchy exists:
//   - aggregate OWNED nodes only; ghost nodes take their OWNER's aggregate (global id) via halo exchange;
//   - the coarse level then has the same shape as level 0: owned rows + a ghost ring (identity rows) + a
//     halo topology of its own -> the build recurses with the same machinery;
//   - the coarsest is REPLICATED: Allgatherv the (small) coarse triplets, every rank factors the same dense
//     LU, per-V-cycle Allgatherv of the owned RHS slices -> everyone holds the full coarse solution.
// Single-rank: none of this is instantiated (all callers gate on numRanks > 1).

#include <mpi.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/binary_search.h>
#include <thrust/copy.h>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <type_traits>

namespace mars {

// Per-level halo: exchange ghost-slot values from their owner ranks. Index lists are LOCAL node ids at
// that level (a "node" = an aggregate on coarse levels). Vectors at a level are ncomp-interleaved
// (ncomp=4 for the coupled [u,v,w,p] blocks, ncomp=1 for node scalars like aggregate ids).
struct AcmLevelHalo {
    std::vector<int> peers;                    // peer ranks, ascending
    std::vector<int> sendOff, recvOff;         // CSR offsets per peer into the id lists
    thrust::device_vector<int> sendIds;        // local OWNED ids to pack, grouped by peer
    thrust::device_vector<int> recvIds;        // local GHOST ids to unpack into, grouped by peer
    mutable thrust::device_vector<double> sendBuf, recvBuf;   // pack/unpack staging (grown on demand)
    bool empty() const { return peers.empty(); }
};

template<typename RealType>
__global__ void acmHaloPackKernel(const RealType* v, const int* ids, RealType* buf, int n, int ncomp)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x; if (i >= n) return;
    for (int c = 0; c < ncomp; ++c) buf[i * ncomp + c] = v[ids[i] * ncomp + c];
}
template<typename RealType>
__global__ void acmHaloUnpackKernel(RealType* v, const int* ids, const RealType* buf, int n, int ncomp)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x; if (i >= n) return;
    for (int c = 0; c < ncomp; ++c) v[ids[i] * ncomp + c] = buf[i * ncomp + c];
}

// owner -> ghost refresh of an ncomp-interleaved level vector (CUDA-aware MPI, device buffers).
template<typename RealType>
inline void acmHaloExchange(const AcmLevelHalo& H, RealType* d_v, int ncomp, int tag = 0x414d)
{
    if (H.empty()) return;
    static_assert(std::is_same_v<RealType, double>, "AcmLevelHalo staging buffers are double");
    const int nSend = H.sendOff.back(), nRecv = H.recvOff.back();
    auto& sbuf = H.sendBuf; auto& rbuf = H.recvBuf;      // per-halo staging (owned by the level, not static)
    if ((int)sbuf.size() < nSend * ncomp) sbuf.resize((size_t)nSend * ncomp);
    if ((int)rbuf.size() < nRecv * ncomp) rbuf.resize((size_t)nRecv * ncomp);
    const int blk = 256;
    if (nSend > 0)
        acmHaloPackKernel<RealType><<<(nSend + blk - 1) / blk, blk>>>(
            d_v, thrust::raw_pointer_cast(H.sendIds.data()), thrust::raw_pointer_cast(sbuf.data()), nSend, ncomp);
    cudaDeviceSynchronize();
    const MPI_Datatype mt = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    std::vector<MPI_Request> reqs; reqs.reserve(2 * H.peers.size());
    for (size_t p = 0; p < H.peers.size(); ++p) {
        int rc = H.recvOff[p + 1] - H.recvOff[p];
        if (rc > 0) { MPI_Request r;
            MPI_Irecv(thrust::raw_pointer_cast(rbuf.data()) + (size_t)H.recvOff[p] * ncomp, rc * ncomp, mt,
                      H.peers[p], tag, MPI_COMM_WORLD, &r); reqs.push_back(r); }
    }
    for (size_t p = 0; p < H.peers.size(); ++p) {
        int sc = H.sendOff[p + 1] - H.sendOff[p];
        if (sc > 0) { MPI_Request r;
            MPI_Isend(thrust::raw_pointer_cast(sbuf.data()) + (size_t)H.sendOff[p] * ncomp, sc * ncomp, mt,
                      H.peers[p], tag, MPI_COMM_WORLD, &r); reqs.push_back(r); }
    }
    if (!reqs.empty()) MPI_Waitall((int)reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);
    if (nRecv > 0)
        acmHaloUnpackKernel<RealType><<<(nRecv + blk - 1) / blk, blk>>>(
            d_v, thrust::raw_pointer_cast(H.recvIds.data()), thrust::raw_pointer_cast(rbuf.data()), nRecv, ncomp);
    cudaDeviceSynchronize();
}

// Build the NEXT level's halo from the exchanged global aggregate ids.
//   gAggGhost: host, one entry per LOCAL GHOST node at THIS level = the owner's GLOBAL aggregate id.
//   aggOffsets: host, size numRanks+1 = Exscan'd owned-aggregate counts (global id ranges per rank).
// Outputs: ghostAggGid (sorted unique global ids = the coarse ghost ring, its LOCAL index = position + nCoarseOwned)
// and the halo (recv: those ghosts from their owners; send: what other ranks request of me, via Alltoallv).
inline void acmBuildCoarseHalo(const std::vector<long long>& gAggGhost,
                               const std::vector<long long>& aggOffsets,
                               int rank, int numRanks,
                               std::vector<long long>& ghostAggGid, AcmLevelHalo& H,
                               int nCoarseOwned)
{
    // unique sorted ghost aggregate ids (host; ghost counts are small)
    ghostAggGid.assign(gAggGhost.begin(), gAggGhost.end());
    std::sort(ghostAggGid.begin(), ghostAggGid.end());
    ghostAggGid.erase(std::unique(ghostAggGid.begin(), ghostAggGid.end()), ghostAggGid.end());
    // a ghost slot the level halo did not cover would arrive as -1 -> ownerOf(-1) corrupts the heap.
    // Loud abort beats silent UB (the domain halo coverage is an invariant, not a given).
    if (!ghostAggGid.empty() && (ghostAggGid.front() < 0 || ghostAggGid.back() >= aggOffsets.back())) {
        std::fprintf(stderr, "[acm-dist][r%d] ghost aggregate id out of range -> level halo did not cover a ghost slot\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 7);
    }

    auto ownerOf = [&](long long gid) {
        return (int)(std::upper_bound(aggOffsets.begin(), aggOffsets.end(), gid) - aggOffsets.begin()) - 1;
    };
    // recv side: ghost aggregates grouped by owner (ghostAggGid is sorted -> already grouped by owner)
    std::vector<int> recvCnt(numRanks, 0);
    for (long long g : ghostAggGid) ++recvCnt[ownerOf(g)];
    // tell every rank which of ITS aggregates I need (global ids); receive the symmetric requests
    std::vector<int> sendCnt(numRanks, 0);
    MPI_Alltoall(recvCnt.data(), 1, MPI_INT, sendCnt.data(), 1, MPI_INT, MPI_COMM_WORLD);
    std::vector<int> rdis(numRanks + 1, 0), sdis(numRanks + 1, 0);
    for (int r = 0; r < numRanks; ++r) { rdis[r + 1] = rdis[r] + recvCnt[r]; sdis[r + 1] = sdis[r] + sendCnt[r]; }
    std::vector<long long> reqOut(ghostAggGid);                     // sorted = grouped by owner rank
    std::vector<long long> reqIn(sdis[numRanks]);
    MPI_Alltoallv(reqOut.data(), recvCnt.data(), rdis.data(), MPI_LONG_LONG,
                  reqIn.data(), sendCnt.data(), sdis.data(), MPI_LONG_LONG, MPI_COMM_WORLD);

    // assemble the halo lists (skip empty peers)
    H = AcmLevelHalo{};
    std::vector<int> h_send, h_recv;
    for (int r = 0; r < numRanks; ++r) {
        if (r == rank || (recvCnt[r] == 0 && sendCnt[r] == 0)) continue;
        H.peers.push_back(r);
        H.sendOff.push_back((int)h_send.size());
        for (int k = sdis[r]; k < sdis[r + 1]; ++k)
            h_send.push_back((int)(reqIn[k] - aggOffsets[rank]));   // my OWNED aggregate, local id
        H.recvOff.push_back((int)h_recv.size());
        for (int k = rdis[r]; k < rdis[r + 1]; ++k) {
            // position of this gid in ghostAggGid -> local ghost index nCoarseOwned + pos
            int pos = (int)(std::lower_bound(ghostAggGid.begin(), ghostAggGid.end(), reqOut[k]) - ghostAggGid.begin());
            h_recv.push_back(nCoarseOwned + pos);
        }
    }
    H.sendOff.push_back((int)h_send.size());
    H.recvOff.push_back((int)h_recv.size());
    H.sendIds.assign(h_send.begin(), h_send.end());
    H.recvIds.assign(h_recv.begin(), h_recv.end());
}

// Replicated global coarsest: every rank holds the full (small) coarse matrix. Triplets are gathered in
// GLOBAL DOF ids (4*globalAgg+comp); the caller factors the dense matrix once (existing getrf path).
struct AcmGlobalCoarse {
    int nGlobalDof = 0;                        // 4 * total owned aggregates across ranks
    long long myDofBegin = 0;                  // first global DOF owned by this rank
    int myDofCount = 0;
    std::vector<int> dofCounts, dofDispls;     // per-rank owned-DOF counts/displacements (for Allgatherv)
    thrust::device_vector<double> dense;       // column-major nGlobalDof^2 (assembled once per Picard)
};

// Gather owned coarse triplets (host) into the replicated dense matrix. rows/cols are GLOBAL DOF ids.
inline void acmGatherGlobalCoarse(const std::vector<long long>& gRow, const std::vector<long long>& gCol,
                                  const std::vector<double>& val, AcmGlobalCoarse& G, int numRanks)
{
    int nLoc = (int)val.size();
    std::vector<int> cnt(numRanks), dis(numRanks + 1, 0);
    MPI_Allgather(&nLoc, 1, MPI_INT, cnt.data(), 1, MPI_INT, MPI_COMM_WORLD);
    for (int r = 0; r < numRanks; ++r) dis[r + 1] = dis[r] + cnt[r];
    std::vector<long long> aRow(dis[numRanks]), aCol(dis[numRanks]);
    std::vector<double> aVal(dis[numRanks]);
    MPI_Allgatherv(gRow.data(), nLoc, MPI_LONG_LONG, aRow.data(), cnt.data(), dis.data(), MPI_LONG_LONG, MPI_COMM_WORLD);
    MPI_Allgatherv(gCol.data(), nLoc, MPI_LONG_LONG, aCol.data(), cnt.data(), dis.data(), MPI_LONG_LONG, MPI_COMM_WORLD);
    MPI_Allgatherv(val.data(), nLoc, MPI_DOUBLE, aVal.data(), cnt.data(), dis.data(), MPI_DOUBLE, MPI_COMM_WORLD);
    const int n = G.nGlobalDof;
    std::vector<double> h_dense((size_t)n * n, 0.0);
    for (size_t k = 0; k < aVal.size(); ++k)
        h_dense[(size_t)aCol[k] * n + aRow[k]] += aVal[k];   // col-major, duplicates (seam pairs) SUM
    G.dense.assign(h_dense.begin(), h_dense.end());
}

}  // namespace mars

#endif  // MARS_ACM_DIST_HPP
