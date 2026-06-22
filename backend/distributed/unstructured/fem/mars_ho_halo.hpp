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
#include "backend/distributed/unstructured/fem/mars_ho_dof_handler.hpp"

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

    std::map<std::array<long,6>, int> myShared;
    std::vector<std::array<long,6>>   myKeys;
    for (int d = 0; d < (int)dofShared.size(); ++d)
        if (dofShared[d]) { auto k = pack(dofKey[d]); myShared.emplace(k, d); myKeys.push_back(k); }
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
            auto it = myShared.find(k);
            if (it != myShared.end() && pr < dofOwner[it->second]) dofOwner[it->second] = pr;
        }
    }
}

#ifdef __CUDACC__
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
        std::map<std::array<long,6>, int> keyToLocal;
        for (int d = 0; d < numDof; ++d) if (dofBoundary[d]) keyToLocal.emplace(packKey(dofKey[d]), d);

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
                auto it = keyToLocal.find(k);
                if (it != keyToLocal.end()) sendLocal[i].push_back(it->second);
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

#ifdef __CUDACC__
    void uploadDevice()
    {
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
