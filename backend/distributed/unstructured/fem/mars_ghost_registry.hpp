#pragma once

// Named, per-interface ghost registry: the generic form of the receiver-driven halo build that
// mars_ho_halo.hpp's HoHalo already proves at scale, lifted off HODofHandler::DofKey so ANY entity
// set with a canonical global key can have its own exchange pattern.
//
// Why this exists: MARS today has two hand-written exchange paths (NodeHaloTopology/exchangeNodeHalo
// for P1 nodes, HoHalo for high-order DOF). Both work, but each new ghost set means writing another
// class. STK instead lets callers do create_ghosting(name) + change_ghosting(entities) and then
// communicate_field_data over exactly that set -- the gap flagged in internal-notes/stk-mars.md.
// This is that primitive.
//
// Receiver-driven, exactly like HoHalo: each rank asks the owner for its ghost entries BY KEY, and
// the owner's send list is built from the received keys IN RECEIVED ORDER. So send[A->B] and
// recv[B<-A] agree slot-for-slot by construction -- no count negotiation, no ownership over-claim,
// hence none of the MPI_ERR_TRUNCATE class of failures that a sender-driven build invites.
//
// Deliberately depends on <mpi.h> + std only, NOT domain.hpp: that keeps it host-unit-testable
// without a CUDA toolchain (see the multi-rank gates in the test), and lets the device path be
// layered on later the same way HoHalo did (host vectors first, gather/scatter kernels after).

#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

// nvcc/hipcc define these; USE_CUDA does NOT distinguish them, it is set for every translation
// unit in a CUDA build including plain .cpp ones (the MPI gate is exactly that). Guarding kernels
// on USE_CUDA hands __global__ and <<<>>> to the host compiler -- that broke the build once.
#if defined(__CUDACC__) || defined(__HIPCC__)
#define MARS_GR_CUDA 1
#include <thrust/device_vector.h>
#endif

namespace mars
{
namespace fem
{

// KeyT is the canonical global identity of an entity: it must be trivially copyable, orderable, and
// identical on every rank that sees the entity (an SFC key, or a packed tuple like HoHalo's
// std::array<long,6>). Local indices are NOT usable across ranks -- that is the whole point of
// keying the request.
#ifdef MARS_GR_CUDA
namespace detail
{

// idx < 0 marks a slot the owner could not resolve; it contributes nothing in either direction,
// matching the host loops.
template<typename T>
__global__ void ghostGatherKernel(const T* v, const int* idx, int n, T* out)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    out[i] = (idx[i] >= 0) ? v[idx[i]] : T{};
}

template<typename T>
__global__ void ghostScatterKernel(const T* in, const int* idx, int n, T* v)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    if (idx[i] >= 0) v[idx[i]] = in[i];
}

// atomicAdd, not `+=`: sendIdx_ is indexed per (peer, entity), so an owned entity ghosted by
// several peers appears several times and their contributions collide. The host loop gets away
// with a plain += only because it is serial.
template<typename T>
__global__ void ghostScatterAddKernel(const T* in, const int* idx, int n, T* v)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    if (idx[i] >= 0) atomicAdd(&v[idx[i]], in[i]);
}

}  // namespace detail
#endif  // MARS_GR_CUDA

template<typename KeyT>
class GhostRegistry
{
public:
    explicit GhostRegistry(std::string name = {})
        : name_(std::move(name))
    {
    }

    ~GhostRegistry() { freeComm(); }

    // Move-only: the duplicated communicator is an owned resource, and copying it would double-free
    // (or, worse, leave two registries exchanging on the same comm).
    GhostRegistry(const GhostRegistry&)            = delete;
    GhostRegistry& operator=(const GhostRegistry&) = delete;

    GhostRegistry(GhostRegistry&& o) noexcept { moveFrom(std::move(o)); }
    GhostRegistry& operator=(GhostRegistry&& o) noexcept
    {
        if (this != &o)
        {
            freeComm();
            moveFrom(std::move(o));
        }
        return *this;
    }

    const std::string& name() const { return name_; }

    // Build the exchange pattern.
    //   owner[i]      : rank owning entity i; < 0 means "not shared" (interior). owner[i] == myRank
    //                   means this rank owns it and may have to SEND it.
    //   key[i]        : canonical global key of entity i.
    //   participates  : 1 if entity i may take part in this ghosting. Only these are keyed and
    //                   scanned -- keeping the build O(interface), not O(numEntities), which is what
    //                   makes it affordable at scale (HoHalo hit a ~78 GB host map by keying all).
    //   candidatePeers: ranks to consider. A superset is fine (peers with no traffic are dropped);
    //                   reuse the node halo's peer list, since a shared entity's owner is a
    //                   mesh neighbour.
    //
    // comm is DUPLICATED, so this registry's messages cannot collide with any other registry's or
    // with the caller's own traffic -- that isolation is what allows several named ghostings to be
    // live at once, which fixed-tag HoHalo cannot do.
    void build(size_t numEntities,
               const std::vector<int>& owner,
               const std::vector<KeyT>& key,
               const std::vector<uint8_t>& participates,
               int myRank,
               const std::vector<int>& candidatePeers,
               MPI_Comm comm = MPI_COMM_WORLD)
    {
        freeComm();
        MPI_Comm_dup(comm, &comm_);
        myRank_ = myRank;

        const int np = static_cast<int>(candidatePeers.size());

        // Sorted (key -> local id) over participating entities only. Sorted vector + binary search,
        // not std::map: HoHalo found the tree ops dominated setup cost at scale.
        std::vector<std::pair<KeyT, int>> keyToLocal;
        for (size_t e = 0; e < numEntities; ++e)
            if (participates[e]) keyToLocal.push_back({key[e], static_cast<int>(e)});
        std::sort(keyToLocal.begin(), keyToLocal.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });

        // Peers are walked in sorted order so both sides post their messages in the same peer
        // order; the caller's list may be in any order.
        std::vector<int> peers = candidatePeers;
        std::sort(peers.begin(), peers.end());
        auto slot = [&](int rank) -> int {
            auto it = std::lower_bound(peers.begin(), peers.end(), rank);
            if (it == peers.end() || *it != rank) return -1;
            return static_cast<int>(it - peers.begin());
        };

        // Ghosts (owned elsewhere) grouped by owner -> my recv lists + the keys I request.
        std::vector<std::vector<int>>  recvLocal(np);
        std::vector<std::vector<KeyT>> reqKeys(np);

        for (size_t e = 0; e < numEntities; ++e)
        {
            if (!participates[e]) continue;
            const int o = owner[e];
            if (o < 0 || o == myRank) continue;   // interior, or mine to send rather than receive
            const int s = slot(o);
            if (s < 0) continue;                  // owner not among the candidate peers
            recvLocal[s].push_back(static_cast<int>(e));
            reqKeys[s].push_back(key[e]);
        }

        // 1) request counts, 2) request keys. Symmetric over candidatePeers -> no deadlock.
        std::vector<int> reqSendCnt(np), reqRecvCnt(np, 0);
        for (int i = 0; i < np; ++i) reqSendCnt[i] = static_cast<int>(reqKeys[i].size());
        exchangeCounts(peers, reqSendCnt, reqRecvCnt);

        std::vector<std::vector<KeyT>> gotKeys(np);
        {
            std::vector<MPI_Request> rq;
            rq.reserve(2 * np);
            for (int i = 0; i < np; ++i)
            {
                gotKeys[i].resize(reqRecvCnt[i]);
                if (reqRecvCnt[i] > 0)
                {
                    MPI_Request r;
                    MPI_Irecv(gotKeys[i].data(), static_cast<int>(reqRecvCnt[i] * sizeof(KeyT)),
                              MPI_BYTE, peers[i], kTagKeys, comm_, &r);
                    rq.push_back(r);
                }
            }
            for (int i = 0; i < np; ++i)
            {
                if (reqSendCnt[i] > 0)
                {
                    MPI_Request r;
                    MPI_Isend(reqKeys[i].data(), static_cast<int>(reqSendCnt[i] * sizeof(KeyT)),
                              MPI_BYTE, peers[i], kTagKeys, comm_, &r);
                    rq.push_back(r);
                }
            }
            if (!rq.empty()) MPI_Waitall(static_cast<int>(rq.size()), rq.data(), MPI_STATUSES_IGNORE);
        }

        // My send list to peer i = its requested keys resolved to my local ids, IN RECEIVED ORDER,
        // so my send slot j lines up with peer i's recv slot j.
        std::vector<std::vector<int>> sendLocal(np);
        for (int i = 0; i < np; ++i)
        {
            sendLocal[i].reserve(gotKeys[i].size());
            for (const KeyT& k : gotKeys[i])
            {
                auto it = std::lower_bound(keyToLocal.begin(), keyToLocal.end(), k,
                                           [](const auto& a, const KeyT& kk) { return a.first < kk; });
                // A key we cannot resolve would silently shorten the send list and break the
                // slot-for-slot pairing, so push a sentinel and count it instead.
                if (it != keyToLocal.end() && it->first == k) sendLocal[i].push_back(it->second);
                else { sendLocal[i].push_back(-1); ++unresolvedIncoming_; }
            }
        }

        // Compact to CSR over peers that actually carry traffic.
        peers_.clear();
        sendOffsets_.assign(1, 0);
        recvOffsets_.assign(1, 0);
        sendIdx_.clear();
        recvIdx_.clear();
        ghostOwner_.clear();
        for (int i = 0; i < np; ++i)
        {
            if (sendLocal[i].empty() && recvLocal[i].empty()) continue;
            peers_.push_back(peers[i]);
            for (int v : sendLocal[i]) sendIdx_.push_back(v);
            for (int v : recvLocal[i])
            {
                recvIdx_.push_back(v);
                ghostOwner_.push_back({v, peers[i]});   // O(surface), never O(numEntities)
            }
            sendOffsets_.push_back(static_cast<int>(sendIdx_.size()));
            recvOffsets_.push_back(static_cast<int>(recvIdx_.size()));
        }
        std::sort(ghostOwner_.begin(), ghostOwner_.end());

        // A rank only learns about an unresolvable key when it is the one being ASKED, so a
        // requester would otherwise sit on a silently-zeroed ghost slot forever. One int of
        // all-reduce at build time (build is already collective) gives every rank the same verdict.
        long localUnresolved = static_cast<long>(unresolvedIncoming_);
        MPI_Allreduce(&localUnresolved, &globalUnresolved_, 1, MPI_LONG, MPI_SUM, comm_);
    }

    // ---- STK-shaped queries -------------------------------------------------------------------

    const std::vector<int>& peers() const { return peers_; }

    // Local ids this rank SENDS to / RECEIVES from peers_[p]  (STK: Ghosting::send_list/receive_list).
    std::pair<const int*, const int*> sendList(int p) const
    {
        return {sendIdx_.data() + sendOffsets_[p], sendIdx_.data() + sendOffsets_[p + 1]};
    }
    std::pair<const int*, const int*> recvList(int p) const
    {
        return {recvIdx_.data() + recvOffsets_[p], recvIdx_.data() + recvOffsets_[p + 1]};
    }

    int sendCount(int p) const { return sendOffsets_[p + 1] - sendOffsets_[p]; }
    int recvCount(int p) const { return recvOffsets_[p + 1] - recvOffsets_[p]; }
    int totalSend() const { return sendOffsets_.back(); }
    int totalRecv() const { return recvOffsets_.back(); }

    // STK: parallel_owner_rank(entity). Ghost entries carry their owner; anything else is ours.
    // Stored per-ghost (O(interface)) rather than as a per-entity array, which would reintroduce
    // the O(numEntities) host allocation HoHalo explicitly avoids.
    int ownerRank(int localId) const
    {
        auto it = std::lower_bound(ghostOwner_.begin(), ghostOwner_.end(),
                                   std::pair<int, int>{localId, -1});
        if (it != ghostOwner_.end() && it->first == localId) return it->second;
        return myRank_;
    }

    bool isGhost(int localId) const { return ownerRank(localId) != myRank_; }

    // Keys THIS rank was asked for but could not resolve -- i.e. a peer believes it ghosts an
    // entity that this rank does not mark as participating. Counted on the OWNER side (the rank
    // being asked), so a pure requester sees 0 here even when its own requests went unanswered:
    // use isConsistent()/globalUnresolvedRequests() for a verdict every rank can act on. The
    // matching send slots hold -1 and are skipped by forward/reverseAdd.
    size_t unresolvedIncomingRequests() const { return unresolvedIncoming_; }

    // Summed over all ranks in the communicator, so every rank agrees. Non-zero means the two
    // sides of some interface disagree about identity or participation -- `participates` must mark
    // a shared entity on the OWNER as well as on every rank that ghosts it.
    long globalUnresolvedRequests() const { return globalUnresolved_; }
    bool isConsistent() const { return globalUnresolved_ == 0; }

    // ---- exchange -----------------------------------------------------------------------------

    // FORWARD: owner -> ghost, ghost slots OVERWRITTEN with the owner's value (STK:
    // communicate_field_data / copy_owned_to_shared over this ghosting).
    //
    // The `Host` suffix is the exception, not the rule: MARS is device-by-default, so the
    // unsuffixed forward()/reverseAdd() below are the GPU ones. These stay as the reference the
    // MPI gates run against without a GPU.
    template<typename T>
    void forwardHost(std::vector<T>& v) const
    {
        const int ns = totalSend(), nr = totalRecv();
        std::vector<T> sbuf(ns), rbuf(nr);
        for (int i = 0; i < ns; ++i) sbuf[i] = (sendIdx_[i] >= 0) ? v[sendIdx_[i]] : T{};
        exchangeValsHost(sbuf, rbuf, sendOffsets_, recvOffsets_, kTagFwd);
        for (int i = 0; i < nr; ++i) v[recvIdx_[i]] = rbuf[i];
    }

    // REVERSE-ADD: ghost -> owner, ghost contributions SUMMED into the owner's slot (the assembly
    // direction).
    template<typename T>
    void reverseAddHost(std::vector<T>& v) const
    {
        const int ns = totalRecv(), nr = totalSend();
        std::vector<T> sbuf(ns), rbuf(nr);
        for (int i = 0; i < ns; ++i) sbuf[i] = v[recvIdx_[i]];
        exchangeValsHost(sbuf, rbuf, recvOffsets_, sendOffsets_, kTagRev);
        for (int i = 0; i < nr; ++i)
            if (sendIdx_[i] >= 0) v[sendIdx_[i]] += rbuf[i];
    }

#ifdef MARS_GR_CUDA
    // ---- device exchange (the default path) ---------------------------------------------------

    // Same semantics as forwardHost/reverseAddHost, on a device pointer, with MPI moving device
    // buffers directly. Index lists are uploaded once and reused -- the build is host-side setup,
    // but nothing per-step crosses the bus.
    template<typename T>
    void forward(T* d_v) const
    {
        ensureIndicesUploaded();
        const int ns = totalSend(), nr = totalRecv();
        thrust::device_vector<T> d_s(ns), d_r(nr);
        launchGather(d_v, thrust::raw_pointer_cast(d_sendIdx_.data()), ns,
                     thrust::raw_pointer_cast(d_s.data()));
        exchangeVals(thrust::raw_pointer_cast(d_s.data()),
                           thrust::raw_pointer_cast(d_r.data()), sendOffsets_, recvOffsets_,
                           kTagFwd);
        launchScatter(thrust::raw_pointer_cast(d_r.data()),
                      thrust::raw_pointer_cast(d_recvIdx_.data()), nr, d_v, /*add=*/false);
    }

    template<typename T>
    void reverseAdd(T* d_v) const
    {
        ensureIndicesUploaded();
        const int ns = totalRecv(), nr = totalSend();
        thrust::device_vector<T> d_s(ns), d_r(nr);
        launchGather(d_v, thrust::raw_pointer_cast(d_recvIdx_.data()), ns,
                     thrust::raw_pointer_cast(d_s.data()));
        exchangeVals(thrust::raw_pointer_cast(d_s.data()),
                           thrust::raw_pointer_cast(d_r.data()), recvOffsets_, sendOffsets_,
                           kTagRev);
        launchScatter(thrust::raw_pointer_cast(d_r.data()),
                      thrust::raw_pointer_cast(d_sendIdx_.data()), nr, d_v, /*add=*/true);
    }
#endif  // MARS_GR_CUDA

private:
    // Tags are still distinct per message KIND, but collision across registries is prevented by the
    // duplicated communicator, not by these values.
    static constexpr int kTagCounts = 1;
    static constexpr int kTagKeys   = 2;
    static constexpr int kTagFwd    = 3;
    static constexpr int kTagRev    = 4;

    // Freeing a communicator after MPI_Finalize is illegal and aborts the job. A registry held as a
    // member of a long-lived object (ElementDomain, a solver) is destroyed during scope teardown,
    // which can run after finalize, so check rather than trust call order.
    void freeComm()
    {
        if (comm_ == MPI_COMM_NULL) return;
        int finalized = 0;
        MPI_Finalized(&finalized);
        if (!finalized) MPI_Comm_free(&comm_);
        comm_ = MPI_COMM_NULL;
    }

    void moveFrom(GhostRegistry&& o) noexcept
    {
        name_        = std::move(o.name_);
        peers_       = std::move(o.peers_);
        sendOffsets_ = std::move(o.sendOffsets_);
        recvOffsets_ = std::move(o.recvOffsets_);
        sendIdx_     = std::move(o.sendIdx_);
        recvIdx_     = std::move(o.recvIdx_);
        ghostOwner_  = std::move(o.ghostOwner_);
        myRank_      = o.myRank_;
        unresolvedIncoming_ = o.unresolvedIncoming_;
        globalUnresolved_   = o.globalUnresolved_;
        comm_        = o.comm_;
        o.comm_      = MPI_COMM_NULL;
#ifdef MARS_GR_CUDA
        // Carry the uploaded index lists across the move, or the moved-to registry would keep
        // idxUploaded_ = false and silently re-upload -- correct, but a wasted transfer per move.
        d_sendIdx_     = std::move(o.d_sendIdx_);
        d_recvIdx_     = std::move(o.d_recvIdx_);
        idxUploaded_   = o.idxUploaded_;
        o.idxUploaded_ = false;
#endif
    }

    void exchangeCounts(const std::vector<int>& peers, const std::vector<int>& send,
                        std::vector<int>& recv) const
    {
        const int np = static_cast<int>(peers.size());
        std::vector<MPI_Request> rq;
        rq.reserve(2 * np);
        for (int i = 0; i < np; ++i)
        {
            MPI_Request r;
            MPI_Irecv(&recv[i], 1, MPI_INT, peers[i], kTagCounts, comm_, &r);
            rq.push_back(r);
        }
        for (int i = 0; i < np; ++i)
        {
            MPI_Request r;
            MPI_Isend(&send[i], 1, MPI_INT, peers[i], kTagCounts, comm_, &r);
            rq.push_back(r);
        }
        if (!rq.empty()) MPI_Waitall(static_cast<int>(rq.size()), rq.data(), MPI_STATUSES_IGNORE);
    }

    template<typename T>
    void exchangeValsHost(const std::vector<T>& sbuf, std::vector<T>& rbuf,
                      const std::vector<int>& sOff, const std::vector<int>& rOff, int tag) const
    {
        const int np = static_cast<int>(peers_.size());
        std::vector<MPI_Request> rq;
        rq.reserve(2 * np);
        for (int i = 0; i < np; ++i)
        {
            const int n = rOff[i + 1] - rOff[i];
            if (n <= 0) continue;
            MPI_Request r;
            MPI_Irecv(rbuf.data() + rOff[i], static_cast<int>(n * sizeof(T)), MPI_BYTE, peers_[i],
                      tag, comm_, &r);
            rq.push_back(r);
        }
        for (int i = 0; i < np; ++i)
        {
            const int n = sOff[i + 1] - sOff[i];
            if (n <= 0) continue;
            MPI_Request r;
            MPI_Isend(sbuf.data() + sOff[i], static_cast<int>(n * sizeof(T)), MPI_BYTE, peers_[i],
                      tag, comm_, &r);
            rq.push_back(r);
        }
        if (!rq.empty()) MPI_Waitall(static_cast<int>(rq.size()), rq.data(), MPI_STATUSES_IGNORE);
    }

    std::string name_;
#ifdef MARS_GR_CUDA
    static constexpr int kBlock = 256;

    void ensureIndicesUploaded() const
    {
        if (idxUploaded_) return;
        d_sendIdx_ = sendIdx_;
        d_recvIdx_ = recvIdx_;
        idxUploaded_ = true;
    }

    template<typename T>
    static void launchGather(const T* v, const int* idx, int n, T* out)
    {
        if (n == 0) return;
        detail::ghostGatherKernel<T><<<(n + kBlock - 1) / kBlock, kBlock>>>(v, idx, n, out);
        cudaDeviceSynchronize();
    }

    template<typename T>
    static void launchScatter(const T* in, const int* idx, int n, T* v, bool add)
    {
        if (n == 0) return;
        const int nb = (n + kBlock - 1) / kBlock;
        if (add) detail::ghostScatterAddKernel<T><<<nb, kBlock>>>(in, idx, n, v);
        else     detail::ghostScatterKernel<T><<<nb, kBlock>>>(in, idx, n, v);
        cudaDeviceSynchronize();
    }

    // Mirrors exchangeVals, but the buffers are device pointers handed straight to MPI.
    template<typename T>
    void exchangeVals(const T* sbuf, T* rbuf, const std::vector<int>& sOff,
                            const std::vector<int>& rOff, int tag) const
    {
        const int np = static_cast<int>(peers_.size());
        std::vector<MPI_Request> rq;
        rq.reserve(2 * np);
        for (int i = 0; i < np; ++i)
        {
            const int n = rOff[i + 1] - rOff[i];
            if (n <= 0) continue;
            MPI_Request r;
            MPI_Irecv(rbuf + rOff[i], n * int(sizeof(T)), MPI_BYTE, peers_[i], tag, comm_, &r);
            rq.push_back(r);
        }
        for (int i = 0; i < np; ++i)
        {
            const int n = sOff[i + 1] - sOff[i];
            if (n <= 0) continue;
            MPI_Request r;
            MPI_Isend(const_cast<T*>(sbuf) + sOff[i], n * int(sizeof(T)), MPI_BYTE, peers_[i], tag,
                      comm_, &r);
            rq.push_back(r);
        }
        if (!rq.empty()) MPI_Waitall(int(rq.size()), rq.data(), MPI_STATUSES_IGNORE);
    }

    mutable thrust::device_vector<int> d_sendIdx_, d_recvIdx_;
    mutable bool idxUploaded_ = false;
#endif  // MARS_GR_CUDA

    std::vector<int> peers_;
    std::vector<int> sendOffsets_{0}, recvOffsets_{0};   // CSR over peers_
    std::vector<int> sendIdx_, recvIdx_;                 // local entity ids
    std::vector<std::pair<int, int>> ghostOwner_;        // sorted (local id, owner rank), ghosts only
    int      myRank_     = 0;
    size_t   unresolvedIncoming_ = 0;
    long     globalUnresolved_   = 0;
    MPI_Comm comm_       = MPI_COMM_NULL;
};

} // namespace fem
} // namespace mars
