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
#include <cstring>
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

// One thread per DESTINATION entity, summing the slots that target it -- no atomics.
//
// sendIdx_ is indexed per (peer, entity), so an owned entity ghosted by several peers is written
// by several slots. The obvious answer is atomicAdd, and contention is low enough (1-3 peers per
// entity) that it would not be slow. It is avoided for DETERMINISM: atomicAdd on doubles sums in
// whatever order the hardware happens to finish in, so reverseAdd would not be bitwise
// reproducible run to run -- and the device-vs-host gate compares with EXPECT_DOUBLE_EQ.
//
// The collision pattern comes from the ghosting topology, which is fixed at build time, so it is
// inverted ONCE into (dstIdx, csrOff, csrSlot) and costs nothing per call.
template<typename T>
__global__ void ghostGatherAddKernel(const T* in, const int* dstIdx, const int* csrOff,
                                     const int* csrSlot, int nDst, T* v)
{
    int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= nDst) return;
    T sum{};
    for (int k = csrOff[e]; k < csrOff[e + 1]; ++k) sum += in[csrSlot[k]];
    v[dstIdx[e]] += sum;
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
        keyToLocal_ = keyToLocal;

        // Peers are walked in sorted order so both sides post their messages in the same peer
        // order; the caller's list may be in any order.
        std::vector<int> peers = candidatePeers;
        std::sort(peers.begin(), peers.end());
        peers.erase(std::unique(peers.begin(), peers.end()), peers.end());
        peers.erase(std::remove(peers.begin(), peers.end(), myRank), peers.end());
        // Retained: peers_ below is COMPACTED to peers that actually have traffic, so an update
        // that gives a silent peer its first ghost would find no channel. keyToLocal_ is retained
        // so the owner can resolve later additions without a second full scan.
        candidatePeers_ = peers;
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

        // Order each peer's request list BY KEY, not by local index.
        //
        // This is the invariant the whole incremental path rests on. With local-index order, slot j
        // means "the j-th entity I happen to store", so inserting one entity renumbers every later
        // slot and both sides would have to renegotiate positions. Sorted by key, slot j means "the
        // j-th smallest key I share with this peer" -- the owner resolves the keys in the order it
        // receives them, so both sides land in the same order without any position ever going on
        // the wire. STK's comm list is kept sorted by EntityKey for exactly this reason.
        for (int i = 0; i < np; ++i)
        {
            std::vector<int> ord(reqKeys[i].size());
            for (size_t j = 0; j < ord.size(); ++j) ord[j] = static_cast<int>(j);
            std::sort(ord.begin(), ord.end(),
                      [&](int a, int b) { return reqKeys[i][a] < reqKeys[i][b]; });
            std::vector<KeyT> k2(ord.size());
            std::vector<int>  l2(ord.size());
            for (size_t j = 0; j < ord.size(); ++j)
            {
                k2[j] = reqKeys[i][ord[j]];
                l2[j] = recvLocal[i][ord[j]];
            }
            reqKeys[i].swap(k2);
            recvLocal[i].swap(l2);
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
        // Per-peer lists are the master state from here on; the CSR is a derived view, rebuilt by
        // flatten(). Editing vectors-of-vectors is what makes an incremental update cheap.
        sendIdxP_ = sendLocal;
        recvIdxP_ = recvLocal;
        recvKeysP_ = reqKeys;
        sendKeysP_.assign(np, {});
        for (int i = 0; i < np; ++i) sendKeysP_[i] = gotKeys[i];
        flatten();

        // A rank only learns about an unresolvable key when it is the one being ASKED, so a
        // requester would otherwise sit on a silently-zeroed ghost slot forever. One int of
        // all-reduce at build time (build is already collective) gives every rank the same verdict.
        long localUnresolved = static_cast<long>(unresolvedIncoming_);
        MPI_Allreduce(&localUnresolved, &globalUnresolved_, 1, MPI_LONG, MPI_SUM, comm_);
    }

    // ---- incremental update (STK: change_ghosting(ghosting, add, remove)) ---------------------
    //
    // Adds and removes ghosts WITHOUT rebuilding the match. `add` names entities this rank now
    // wants to receive; `removeKeys` names ghosts it no longer wants.
    //
    // Removing is NOT local bookkeeping. If this rank drops a ghost and the owner keeps sending it,
    // the per-peer send and recv counts disagree and the next exchange misaligns -- silent
    // corruption, not a clean failure. So the delta goes on the wire and the owner edits its send
    // list in place. STK does the same (BulkData::internal_change_ghosting -> comm_sync_send_recv).
    //
    // ONE round, not two: the payload is self-describing and MPI_Probe gives the size, so unlike
    // build() there is no separate count exchange. Point-to-point over candidatePeers_ throughout;
    // no collective, nothing O(P).
    struct Delta
    {
        KeyT key;
        int  localId;
        int  owner;
    };

    void update(const std::vector<Delta>& add, const std::vector<KeyT>& removeKeys)
    {
        const int np = static_cast<int>(candidatePeers_.size());
        if (np == 0) return;
        auto slot = [&](int rank) -> int {
            auto it = std::lower_bound(candidatePeers_.begin(), candidatePeers_.end(), rank);
            if (it == candidatePeers_.end() || *it != rank) return -1;
            return static_cast<int>(it - candidatePeers_.begin());
        };

        // Group my changes by the peer that owns them. Removes are located by key in the recv
        // lists, so the caller does not have to remember which peer a ghost came from.
        std::vector<std::vector<KeyT>> addK(np), remK(np);
        std::vector<std::vector<int>>  addL(np);
        for (const Delta& d : add)
        {
            const int p = slot(d.owner);
            if (p < 0) continue;                       // owner outside the candidate set
            addK[p].push_back(d.key);
            addL[p].push_back(d.localId);
        }
        for (const KeyT& k : removeKeys)
            for (int p = 0; p < np; ++p)
            {
                auto it = std::lower_bound(recvKeysP_[p].begin(), recvKeysP_[p].end(), k);
                if (it != recvKeysP_[p].end() && *it == k) { remK[p].push_back(k); break; }
            }

        // One self-describing message per peer: [nAdd][nRemove] then the keys.
        std::vector<std::vector<char>> out(np);
        for (int p = 0; p < np; ++p)
        {
            const int na = static_cast<int>(addK[p].size()), nr = static_cast<int>(remK[p].size());
            out[p].resize(2 * sizeof(int) + (na + nr) * sizeof(KeyT));
            char* w = out[p].data();
            std::memcpy(w, &na, sizeof(int));                     w += sizeof(int);
            std::memcpy(w, &nr, sizeof(int));                     w += sizeof(int);
            if (na) std::memcpy(w, addK[p].data(), na * sizeof(KeyT)); w += na * sizeof(KeyT);
            if (nr) std::memcpy(w, remK[p].data(), nr * sizeof(KeyT));
        }
        std::vector<MPI_Request> rq;
        rq.reserve(np);
        for (int p = 0; p < np; ++p)
        {
            rq.emplace_back();
            MPI_Isend(out[p].data(), static_cast<int>(out[p].size()), MPI_BYTE, candidatePeers_[p],
                      kTagDelta, comm_, &rq.back());
        }

        // Probe for the size instead of exchanging counts first -- this is the round build() pays
        // that update() does not.
        for (int p = 0; p < np; ++p)
        {
            MPI_Status st;
            MPI_Probe(candidatePeers_[p], kTagDelta, comm_, &st);
            int nbytes = 0;
            MPI_Get_count(&st, MPI_BYTE, &nbytes);
            std::vector<char> in(nbytes);
            MPI_Recv(in.data(), nbytes, MPI_BYTE, candidatePeers_[p], kTagDelta, comm_,
                     MPI_STATUS_IGNORE);
            int na = 0, nr = 0;
            const char* r = in.data();
            std::memcpy(&na, r, sizeof(int)); r += sizeof(int);
            std::memcpy(&nr, r, sizeof(int)); r += sizeof(int);
            std::vector<KeyT> gotAdd(na), gotRem(nr);
            if (na) { std::memcpy(gotAdd.data(), r, na * sizeof(KeyT)); r += na * sizeof(KeyT); }
            if (nr) std::memcpy(gotRem.data(), r, nr * sizeof(KeyT));

            // I am the OWNER for these: stop sending the removed keys, start sending the added.
            std::sort(gotRem.begin(), gotRem.end());
            eraseKeys(sendKeysP_[p], sendIdxP_[p], gotRem);
            std::sort(gotAdd.begin(), gotAdd.end());
            std::vector<int> resolved(gotAdd.size(), -1);
            for (size_t j = 0; j < gotAdd.size(); ++j)
            {
                auto it = std::lower_bound(keyToLocal_.begin(), keyToLocal_.end(), gotAdd[j],
                                           [](const auto& a, const KeyT& k) { return a.first < k; });
                if (it != keyToLocal_.end() && it->first == gotAdd[j]) resolved[j] = it->second;
                else ++unresolvedIncoming_;   // same contract as build(): -1 slots contribute nothing
            }
            insertKeys(sendKeysP_[p], sendIdxP_[p], gotAdd, resolved);
        }
        if (!rq.empty()) MPI_Waitall(static_cast<int>(rq.size()), rq.data(), MPI_STATUSES_IGNORE);

        // My own recv side needs no message -- I already know what I asked for, and applying the
        // identical sorted edit keeps me aligned with the owner by construction.
        for (int p = 0; p < np; ++p)
        {
            std::sort(remK[p].begin(), remK[p].end());
            eraseKeys(recvKeysP_[p], recvIdxP_[p], remK[p]);
            std::vector<int> ord(addK[p].size());
            for (size_t j = 0; j < ord.size(); ++j) ord[j] = static_cast<int>(j);
            std::sort(ord.begin(), ord.end(),
                      [&](int a, int b) { return addK[p][a] < addK[p][b]; });
            std::vector<KeyT> ak(ord.size());
            std::vector<int>  al(ord.size());
            for (size_t j = 0; j < ord.size(); ++j) { ak[j] = addK[p][ord[j]]; al[j] = addL[p][ord[j]]; }
            insertKeys(recvKeysP_[p], recvIdxP_[p], ak, al);
        }
        flatten();
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
                      thrust::raw_pointer_cast(d_recvIdx_.data()), nr, d_v);
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
        launchGatherAdd(thrust::raw_pointer_cast(d_r.data()), d_v);
    }
#endif  // MARS_GR_CUDA

private:
    // Tags are still distinct per message KIND, but collision across registries is prevented by the
    // duplicated communicator, not by these values.
    static constexpr int kTagCounts = 1;
    static constexpr int kTagKeys   = 2;
    static constexpr int kTagFwd    = 3;
    static constexpr int kTagRev    = 4;
    static constexpr int kTagDelta  = 5;

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

    // Both lists stay sorted by key and parallel, so an edit is a linear merge on each -- no
    // search per element, and no positions ever exchanged.
    static void eraseKeys(std::vector<KeyT>& keys, std::vector<int>& idx,
                          const std::vector<KeyT>& drop)
    {
        if (drop.empty() || keys.empty()) return;
        std::vector<KeyT> k2;
        std::vector<int>  i2;
        k2.reserve(keys.size());
        i2.reserve(idx.size());
        size_t d = 0;
        for (size_t j = 0; j < keys.size(); ++j)
        {
            while (d < drop.size() && drop[d] < keys[j]) ++d;
            if (d < drop.size() && drop[d] == keys[j]) continue;
            k2.push_back(keys[j]);
            i2.push_back(idx[j]);
        }
        keys.swap(k2);
        idx.swap(i2);
    }

    // Duplicates are dropped: re-adding a ghost that is already there must be a no-op, or the two
    // sides would disagree on the slot count.
    static void insertKeys(std::vector<KeyT>& keys, std::vector<int>& idx,
                           const std::vector<KeyT>& addK, const std::vector<int>& addI)
    {
        if (addK.empty()) return;
        std::vector<KeyT> k2;
        std::vector<int>  i2;
        k2.reserve(keys.size() + addK.size());
        i2.reserve(idx.size() + addI.size());
        size_t a = 0, b = 0;
        while (a < keys.size() || b < addK.size())
        {
            if (b >= addK.size() || (a < keys.size() && keys[a] < addK[b]))
            { k2.push_back(keys[a]); i2.push_back(idx[a]); ++a; }
            else if (a < keys.size() && keys[a] == addK[b])
            { k2.push_back(keys[a]); i2.push_back(idx[a]); ++a; ++b; }
            else
            { k2.push_back(addK[b]); i2.push_back(addI[b]); ++b; }
        }
        keys.swap(k2);
        idx.swap(i2);
    }

    // Rebuild the CSR view from the per-peer master state. peers_ keeps only peers with traffic,
    // which is why candidatePeers_ is retained separately.
    void flatten()
    {
        peers_.clear();
        sendIdx_.clear();
        recvIdx_.clear();
        sendOffsets_.assign(1, 0);
        recvOffsets_.assign(1, 0);
        ghostOwner_.clear();
        for (size_t i = 0; i < candidatePeers_.size(); ++i)
        {
            if (sendIdxP_[i].empty() && recvIdxP_[i].empty()) continue;
            peers_.push_back(candidatePeers_[i]);
            for (int v : sendIdxP_[i]) sendIdx_.push_back(v);
            for (int v : recvIdxP_[i])
            {
                recvIdx_.push_back(v);
                ghostOwner_.push_back({v, candidatePeers_[i]});   // O(surface)
            }
            sendOffsets_.push_back(static_cast<int>(sendIdx_.size()));
            recvOffsets_.push_back(static_cast<int>(recvIdx_.size()));
        }
        std::sort(ghostOwner_.begin(), ghostOwner_.end());
#ifdef MARS_GR_CUDA
        // The device index lists and the reverseAdd inverse map are now stale.
        idxUploaded_ = false;
#endif
    }

    void moveFrom(GhostRegistry&& o) noexcept
    {
        name_        = std::move(o.name_);
        peers_       = std::move(o.peers_);
        sendOffsets_ = std::move(o.sendOffsets_);
        recvOffsets_ = std::move(o.recvOffsets_);
        sendIdx_     = std::move(o.sendIdx_);
        recvIdx_     = std::move(o.recvIdx_);
        candidatePeers_ = std::move(o.candidatePeers_);
        sendIdxP_    = std::move(o.sendIdxP_);
        recvIdxP_    = std::move(o.recvIdxP_);
        sendKeysP_   = std::move(o.sendKeysP_);
        recvKeysP_   = std::move(o.recvKeysP_);
        keyToLocal_  = std::move(o.keyToLocal_);
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
        d_dstIdx_      = std::move(o.d_dstIdx_);
        d_csrOff_      = std::move(o.d_csrOff_);
        d_csrSlot_     = std::move(o.d_csrSlot_);
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

        // Invert sendIdx_ once: for every distinct owned entity, the list of send slots that
        // contribute to it. Sorted by entity, so the summation order is fixed and the result is
        // reproducible. Slots holding -1 (the owner could not resolve them) are dropped here
        // rather than tested in the kernel.
        std::vector<std::pair<int, int>> pairs;   // (local entity id, slot)
        pairs.reserve(sendIdx_.size());
        for (size_t i = 0; i < sendIdx_.size(); ++i)
            if (sendIdx_[i] >= 0) pairs.emplace_back(sendIdx_[i], int(i));
        std::sort(pairs.begin(), pairs.end());

        std::vector<int> dst, off{0}, slot;
        slot.reserve(pairs.size());
        for (size_t i = 0; i < pairs.size();)
        {
            const int e = pairs[i].first;
            dst.push_back(e);
            while (i < pairs.size() && pairs[i].first == e) { slot.push_back(pairs[i].second); ++i; }
            off.push_back(int(slot.size()));
        }
        d_dstIdx_  = dst;
        d_csrOff_  = off;
        d_csrSlot_ = slot;
        idxUploaded_ = true;
    }

    template<typename T>
    static void launchGather(const T* v, const int* idx, int n, T* out)
    {
        if (n == 0) return;
        detail::ghostGatherKernel<T><<<(n + kBlock - 1) / kBlock, kBlock>>>(v, idx, n, out);
        // REQUIRED, not leftover: MPI reads `out` next and is not stream-ordered against this
        // kernel. The scatter/gather-add on the other side need no sync -- whatever touches the
        // field afterwards is stream-ordered behind them.
        cudaDeviceSynchronize();
    }

    template<typename T>
    static void launchScatter(const T* in, const int* idx, int n, T* v)
    {
        if (n == 0) return;
        detail::ghostScatterKernel<T><<<(n + kBlock - 1) / kBlock, kBlock>>>(in, idx, n, v);
    }

    template<typename T>
    void launchGatherAdd(const T* in, T* v) const
    {
        const int nDst = int(d_dstIdx_.size());
        if (nDst == 0) return;
        detail::ghostGatherAddKernel<T><<<(nDst + kBlock - 1) / kBlock, kBlock>>>(
            in, thrust::raw_pointer_cast(d_dstIdx_.data()),
            thrust::raw_pointer_cast(d_csrOff_.data()),
            thrust::raw_pointer_cast(d_csrSlot_.data()), nDst, v);
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
    // Inverse of sendIdx_: entity -> contributing slots, so reverseAdd needs no atomics.
    mutable thrust::device_vector<int> d_dstIdx_, d_csrOff_, d_csrSlot_;
    mutable bool idxUploaded_ = false;
#endif  // MARS_GR_CUDA

    std::vector<int> peers_;
    std::vector<int> sendOffsets_{0}, recvOffsets_{0};   // CSR over peers_
    std::vector<int> sendIdx_, recvIdx_;                 // local entity ids
    std::vector<std::pair<int, int>> ghostOwner_;        // sorted (local id, owner rank), ghosts only
    // Master state for incremental edits: per peer, key-sorted and parallel to the index lists.
    std::vector<int>                   candidatePeers_;  // FULL list; peers_ is compacted
    std::vector<std::vector<int>>      sendIdxP_, recvIdxP_;
    std::vector<std::vector<KeyT>>     sendKeysP_, recvKeysP_;
    std::vector<std::pair<KeyT, int>>  keyToLocal_;      // sorted; resolves later additions
    int      myRank_     = 0;
    size_t   unresolvedIncoming_ = 0;
    long     globalUnresolved_   = 0;
    MPI_Comm comm_       = MPI_COMM_NULL;
};

} // namespace fem
} // namespace mars
