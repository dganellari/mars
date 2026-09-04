#ifndef MARS_COARSE_SEARCH_HPP
#define MARS_COARSE_SEARCH_HPP

// Parallel AABB-overlap search -- the MARS equivalent of stk::search::coarse_search.
//
// The geometric pairing step non-conformal interfaces and overset need: given DOMAIN boxes (the
// queries, e.g. slave-side face AABBs) and RANGE boxes (the targets, e.g. master-side face AABBs),
// both distributed, return every overlapping (domain, range) pair with the owning rank of each side.
//
//   coarseSearch      DEVICE. The default, because MARS is device-by-default. Boxes, routed
//                     queries and pairs all stay in device memory; MPI moves device pointers.
//                     The only host traffic is the peer envelopes and the message counts, which
//                     have to be on the host because MPI counts are.
//   coarseSearchHost  the reference. It exists so the communication plan can be pinned down with
//                     mpirun on a laptop -- no GPU, no mesh file -- and so the device version has
//                     something to be compared against element-for-element. Same sorted result.
//
// NO O(P) METADATA. Every exchange is point-to-point over `candidatePeers`, the same contract
// GhostRegistry::build uses: a superset is fine, peers with no traffic cost two tiny messages and
// drop out. An earlier draft used MPI_Allgather for the envelopes and MPI_Alltoall for the counts;
// that is O(P) per rank in three places and does not survive the rank counts this is aimed at.
// Peers come from cstone's own discovery (cstone/traversal/peers.hpp) or from the domain's existing
// halo peer list -- the ranks that could hold a partner are the ranks you already talk to.
//
// Routing is one envelope per PEER, not a distributed tree. Interface and overset boxes live on a
// SURFACE, so an envelope prunes most of the traffic and the fan-out stays small. A tree buys
// nothing until the query set stops being a surface.
//
// STK PARITY, checked against Trilinos develop (stk_search/CoarseSearch.hpp, BoundingBox.hpp,
// CommonSearchUtil.hpp) rather than from memory:
//   - The overlap predicate is IDENTICAL. STK's box-box disjoint test is
//         amax[d] < bmin[d] || bmax[d] < amin[d]
//     which is what `Aabb::overlaps` computes, and it is CLOSED: at amax == bmin the strict `<`
//     is false, so touching boxes DO intersect. That is the behaviour a conformal seam needs.
//   - `enforceSymmetry` mirrors STK's `enforceSearchResultSymmetry`, and like STK it defaults to
//     TRUE: the pair is delivered to the range owner as well as the domain owner. The DG-IP
//     interface assembly needs exactly that -- master rows carry opposing-element columns and
//     vice versa, so both sides have to know about the pair.
//   - We always sort and dedup. STK makes that optional (`sortSearchResults`, default false).
//     Ours is the stricter guarantee, deliberately: it makes the result independent of message
//     arrival order, which is what lets the device path be diffed against the host one.

#include <mpi.h>

// USE_CUDA says the PROJECT is built with CUDA; it is set for every translation unit, including
// plain .cpp ones the host compiler handles -- and the MPI gate is exactly that. Guarding the
// kernels on it hands __global__ and <<<>>> to g++. __CUDACC__/__HIPCC__ are what actually mean
// "a device compiler is reading this".
#if defined(__CUDACC__) || defined(__HIPCC__)
#define MARS_CS_CUDA 1
#endif

#ifdef MARS_CS_CUDA
#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/scan.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#define MARS_CS_FN __host__ __device__
#else
#define MARS_CS_FN
#endif

#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

namespace mars
{
namespace fem
{

template<typename T>
struct Aabb
{
    T lo[3];
    T hi[3];

    MARS_CS_FN static Aabb empty()
    {
        Aabb b;
        for (int d = 0; d < 3; ++d)
        {
            b.lo[d] = std::numeric_limits<T>::max();
            b.hi[d] = std::numeric_limits<T>::lowest();
        }
        return b;
    }

    MARS_CS_FN void extend(const Aabb& o)
    {
        for (int d = 0; d < 3; ++d)
        {
            lo[d] = o.lo[d] < lo[d] ? o.lo[d] : lo[d];
            hi[d] = o.hi[d] > hi[d] ? o.hi[d] : hi[d];
        }
    }

    // Closed test: boxes that only touch count as overlapping. Interface faces on a conformal seam
    // meet exactly at their shared edge, and dropping those would lose real partners.
    MARS_CS_FN bool overlaps(const Aabb& o) const
    {
        for (int d = 0; d < 3; ++d)
            if (hi[d] < o.lo[d] || o.hi[d] < lo[d]) return false;
        return true;
    }
};

// STK calls this IdentProc: which entity, and which rank owns it. The rank matters because the
// partner usually lives somewhere else, and whoever consumes the pair has to know where to go for
// the data -- that is the handoff to the ghosting.
struct IdentProc
{
    uint64_t id;
    int      proc;

    MARS_CS_FN bool operator==(const IdentProc& o) const { return id == o.id && proc == o.proc; }
    MARS_CS_FN bool operator!=(const IdentProc& o) const { return !(*this == o); }
    MARS_CS_FN bool operator<(const IdentProc& o) const
    {
        return proc != o.proc ? proc < o.proc : id < o.id;
    }
};

struct SearchPair
{
    IdentProc domain;
    IdentProc range;

    MARS_CS_FN bool operator<(const SearchPair& o) const
    {
        return domain != o.domain ? domain < o.domain : range < o.range;
    }
    MARS_CS_FN bool operator==(const SearchPair& o) const
    {
        return domain == o.domain && range == o.range;
    }
};

// One query box travelling to a peer that might hold a partner: the box, plus enough identity to
// send the answer home.
template<typename T>
struct QueryBox
{
    Aabb<T>  box;
    uint64_t id;
    int      proc;
};

namespace detail
{

// Peer list with self removed and duplicates dropped. Self needs no MPI -- a rank's own range
// boxes are already local, and posting a message to yourself just to test them is waste.
inline std::vector<int> normalizePeers(const std::vector<int>& candidatePeers, int myRank)
{
    std::vector<int> p;
    p.reserve(candidatePeers.size());
    for (int r : candidatePeers)
        if (r != myRank) p.push_back(r);
    std::sort(p.begin(), p.end());
    p.erase(std::unique(p.begin(), p.end()), p.end());
    return p;
}

// Swap one AABB with every peer. Symmetric over the peer list -- all receives posted before any
// send -- so it cannot deadlock however the peers are ordered.
template<typename T>
inline std::vector<Aabb<T>> exchangeEnvelopes(const Aabb<T>& mine, const std::vector<int>& peers,
                                              MPI_Comm comm)
{
    const int np = int(peers.size());
    std::vector<Aabb<T>> theirs(np);
    std::vector<MPI_Request> rq;
    rq.reserve(2 * np);
    const int kTag = 0x5E00;
    for (int i = 0; i < np; ++i)
    {
        rq.emplace_back();
        MPI_Irecv(&theirs[i], sizeof(Aabb<T>), MPI_BYTE, peers[i], kTag, comm, &rq.back());
    }
    for (int i = 0; i < np; ++i)
    {
        rq.emplace_back();
        MPI_Isend(const_cast<Aabb<T>*>(&mine), sizeof(Aabb<T>), MPI_BYTE, peers[i], kTag, comm,
                  &rq.back());
    }
    if (!rq.empty()) MPI_Waitall(int(rq.size()), rq.data(), MPI_STATUSES_IGNORE);
    return theirs;
}

// Counts, peer-to-peer. This is what replaces MPI_Alltoall: 2*|peers| tiny messages instead of an
// O(P) collective every rank in the job has to enter.
//
// Split into post/wait so the caller can put real work in the shadow of the round trip. `sendCnt`
// and `recvCnt` must stay alive until waitCounts returns -- the sends read from one and the
// receives write into the other.
inline void postCounts(const std::vector<int>& sendCnt, std::vector<int>& recvCnt,
                       const std::vector<int>& peers, MPI_Comm comm, int tag,
                       std::vector<MPI_Request>& rq)
{
    const int np = int(peers.size());
    recvCnt.assign(np, 0);
    rq.clear();
    rq.reserve(2 * np);
    for (int i = 0; i < np; ++i)
    {
        rq.emplace_back();
        MPI_Irecv(&recvCnt[i], 1, MPI_INT, peers[i], tag, comm, &rq.back());
    }
    for (int i = 0; i < np; ++i)
    {
        rq.emplace_back();
        MPI_Isend(const_cast<int*>(&sendCnt[i]), 1, MPI_INT, peers[i], tag, comm, &rq.back());
    }
}

inline void waitCounts(std::vector<MPI_Request>& rq)
{
    if (!rq.empty()) MPI_Waitall(int(rq.size()), rq.data(), MPI_STATUSES_IGNORE);
    rq.clear();
}

inline std::vector<int> exchangeCounts(const std::vector<int>& sendCnt,
                                       const std::vector<int>& peers, MPI_Comm comm, int tag)
{
    const int np = int(peers.size());
    std::vector<int> recvCnt(np, 0);
    std::vector<MPI_Request> rq;
    rq.reserve(2 * np);
    for (int i = 0; i < np; ++i)
    {
        rq.emplace_back();
        MPI_Irecv(&recvCnt[i], 1, MPI_INT, peers[i], tag, comm, &rq.back());
    }
    for (int i = 0; i < np; ++i)
    {
        rq.emplace_back();
        MPI_Isend(const_cast<int*>(&sendCnt[i]), 1, MPI_INT, peers[i], tag, comm, &rq.back());
    }
    if (!rq.empty()) MPI_Waitall(int(rq.size()), rq.data(), MPI_STATUSES_IGNORE);
    return recvCnt;
}

}  // namespace detail

// ---------------------------------------------------------------------------
// Host reference.
// ---------------------------------------------------------------------------
//
// `candidatePeers` is a superset of the ranks that might hold a partner; a peer with no traffic
// costs two tiny messages. On return `out` holds the pairs whose DOMAIN box this rank owns, sorted
// and deduplicated, so the answer does not depend on message arrival order or on rank count.
//
// Collective over the peer graph: every rank must call it, including ranks contributing nothing.
template<typename T>
void coarseSearchHost(const std::vector<Aabb<T>>& domainBoxes,
                      const std::vector<uint64_t>& domainIds,
                      const std::vector<Aabb<T>>& rangeBoxes,
                      const std::vector<uint64_t>& rangeIds,
                      const std::vector<int>& candidatePeers,
                      MPI_Comm comm,
                      std::vector<SearchPair>& out,
                      bool enforceSymmetry = true)
{
    int myRank = 0;
    MPI_Comm_rank(comm, &myRank);
    out.clear();

    const std::vector<int> peers = detail::normalizePeers(candidatePeers, myRank);
    const int np = int(peers.size());

    // Partners on my own rank need no messages at all.
    for (size_t i = 0; i < domainBoxes.size(); ++i)
        for (size_t j = 0; j < rangeBoxes.size(); ++j)
            if (domainBoxes[i].overlaps(rangeBoxes[j]))
                out.push_back(
                    SearchPair{IdentProc{domainIds[i], myRank}, IdentProc{rangeIds[j], myRank}});

    Aabb<T> myEnv = Aabb<T>::empty();
    for (const auto& b : rangeBoxes) myEnv.extend(b);
    const std::vector<Aabb<T>> peerEnv = detail::exchangeEnvelopes(myEnv, peers, comm);

    std::vector<std::vector<QueryBox<T>>> toSend(np);
    std::vector<int> sendCnt(np, 0);
    for (int i = 0; i < np; ++i)
    {
        for (size_t b = 0; b < domainBoxes.size(); ++b)
            if (domainBoxes[b].overlaps(peerEnv[i]))
                toSend[i].push_back(QueryBox<T>{domainBoxes[b], domainIds[b], myRank});
        sendCnt[i] = int(toSend[i].size());
    }
    const std::vector<int> recvCnt = detail::exchangeCounts(sendCnt, peers, comm, 0x5E01);

    std::vector<std::vector<QueryBox<T>>> recvBuf(np);
    std::vector<MPI_Request> rq;
    rq.reserve(2 * np);
    const int kQueryTag = 0x5EA1;
    for (int i = 0; i < np; ++i)
        if (recvCnt[i] > 0)
        {
            recvBuf[i].resize(recvCnt[i]);
            rq.emplace_back();
            MPI_Irecv(recvBuf[i].data(), recvCnt[i] * int(sizeof(QueryBox<T>)), MPI_BYTE, peers[i],
                      kQueryTag, comm, &rq.back());
        }
    for (int i = 0; i < np; ++i)
        if (sendCnt[i] > 0)
        {
            rq.emplace_back();
            MPI_Isend(toSend[i].data(), sendCnt[i] * int(sizeof(QueryBox<T>)), MPI_BYTE, peers[i],
                      kQueryTag, comm, &rq.back());
        }
    if (!rq.empty()) MPI_Waitall(int(rq.size()), rq.data(), MPI_STATUSES_IGNORE);

    // Test what arrived against my range boxes; the pair belongs to the sender.
    std::vector<std::vector<SearchPair>> reply(np);
    std::vector<int> repSendCnt(np, 0);
    for (int i = 0; i < np; ++i)
    {
        for (const auto& q : recvBuf[i])
            for (size_t j = 0; j < rangeBoxes.size(); ++j)
                if (q.box.overlaps(rangeBoxes[j]))
                    reply[i].push_back(
                        SearchPair{IdentProc{q.id, q.proc}, IdentProc{rangeIds[j], myRank}});
        repSendCnt[i] = int(reply[i].size());
        // STK's enforceSearchResultSymmetry: the RANGE owner keeps the pair too, not just the
        // domain owner. I generated these, so my copy is simply the one I already have.
        if (enforceSymmetry) out.insert(out.end(), reply[i].begin(), reply[i].end());
    }
    const std::vector<int> repRecvCnt = detail::exchangeCounts(repSendCnt, peers, comm, 0x5E02);

    std::vector<std::vector<SearchPair>> repRecv(np);
    std::vector<MPI_Request> rq2;
    rq2.reserve(2 * np);
    const int kReplyTag = 0x5EA2;
    for (int i = 0; i < np; ++i)
        if (repRecvCnt[i] > 0)
        {
            repRecv[i].resize(repRecvCnt[i]);
            rq2.emplace_back();
            MPI_Irecv(repRecv[i].data(), repRecvCnt[i] * int(sizeof(SearchPair)), MPI_BYTE,
                      peers[i], kReplyTag, comm, &rq2.back());
        }
    for (int i = 0; i < np; ++i)
        if (repSendCnt[i] > 0)
        {
            rq2.emplace_back();
            MPI_Isend(reply[i].data(), repSendCnt[i] * int(sizeof(SearchPair)), MPI_BYTE, peers[i],
                      kReplyTag, comm, &rq2.back());
        }
    if (!rq2.empty()) MPI_Waitall(int(rq2.size()), rq2.data(), MPI_STATUSES_IGNORE);

    for (int i = 0; i < np; ++i) out.insert(out.end(), repRecv[i].begin(), repRecv[i].end());
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
}

// ---------------------------------------------------------------------------
// Device path -- the default. Only visible to a device compiler; see the MARS_CS_CUDA note above.
//
// Written batched from the start, not "correct now, tune later":
//   - ONE kernel builds the whole (peers x boxes) routing table; no per-peer launch.
//   - ONE exclusive scan over that table. Because it is peer-major, the scan doubles as the
//     compaction offset AND gives every peer's count as a difference of two boundary values.
//   - ONE D2H for all np+1 boundaries, instead of two 4-byte reads per peer.
//   - Queries are received into ONE contiguous buffer, so the local test is a single
//     count/scan/fill over everything rather than one per source rank.
//   - The self-test runs BETWEEN posting the query messages and waiting on them, so the local
//     work overlaps the transfer instead of following it.
//   - Scratch lives in the object. A sliding interface re-searches every step and must not
//     allocate per call.
// ---------------------------------------------------------------------------
#ifdef MARS_CS_CUDA

namespace detail
{

// A named functor, NOT a lambda. nvcc refuses an extended __device__ lambda whose return type has
// to be queried from host code, which is what thrust::transform_reduce does to its operators.
template<typename T>
struct AabbUnion
{
    MARS_CS_FN Aabb<T> operator()(const Aabb<T>& a, const Aabb<T>& b) const
    {
        Aabb<T> r = a;
        r.extend(b);
        return r;
    }
};

// Whole routing table in one launch: flags[r*n + i] is "my box i must go to peer r".
template<typename T>
__global__ void flagAllPeersKernel(const Aabb<T>* boxes, size_t nBoxes, const Aabb<T>* env, int np,
                                   uint8_t* flags)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nBoxes) return;
    const Aabb<T> b = boxes[i];
    for (int r = 0; r < np; ++r) flags[size_t(r) * nBoxes + i] = b.overlaps(env[r]);
}

// The scan is over the whole peer-major table, so scan[r*n+i] is already the GLOBAL slot: each
// peer's queries land contiguously with no per-peer base to add.
template<typename T>
__global__ void gatherAllPeersKernel(const Aabb<T>* boxes, const uint64_t* ids, size_t nBoxes,
                                     const uint8_t* flags, const int* scan, int np, int myRank,
                                     QueryBox<T>* out)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nBoxes) return;
    QueryBox<T> q;
    q.box  = boxes[i];
    q.id   = ids[i];
    q.proc = myRank;
    for (int r = 0; r < np; ++r)
    {
        const size_t o = size_t(r) * nBoxes + i;
        if (flags[o]) out[scan[o]] = q;
    }
}

// np+1 segment boundaries of a scan, so all the counts come home in ONE copy.
__global__ void segmentBoundsKernel(const int* scan, size_t stride, int np, int* bounds)
{
    int r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r > np) return;
    bounds[r] = scan[size_t(r) * stride];
}

// Same idea, but the segments are the per-peer query ranges rather than a fixed stride.
__global__ void gatherAtKernel(const int* scan, const int* at, int n, int* outv)
{
    int r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r >= n) return;
    outv[r] = scan[at[r]];
}

template<typename T>
__global__ void countPairsKernel(const QueryBox<T>* q, size_t nQ, const Aabb<T>* range,
                                 size_t nRange, int* counts)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nQ) return;
    const Aabb<T> b = q[i].box;
    int c = 0;
    for (size_t j = 0; j < nRange; ++j)
        if (b.overlaps(range[j])) ++c;
    counts[i] = c;
}

template<typename T>
__global__ void fillPairsKernel(const QueryBox<T>* q, size_t nQ, const Aabb<T>* range,
                                const uint64_t* rangeIds, size_t nRange, const int* offsets,
                                int myRank, SearchPair* out)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nQ) return;
    const Aabb<T> b = q[i].box;
    const IdentProc dom{q[i].id, q[i].proc};
    int w = offsets[i];
    for (size_t j = 0; j < nRange; ++j)
        if (b.overlaps(range[j]))
        {
            out[w].domain = dom;
            out[w].range  = IdentProc{rangeIds[j], myRank};
            ++w;
        }
}

// My own boxes as queries, without materialising a flag array of all ones.
template<typename T>
__global__ void selfQueriesKernel(const Aabb<T>* boxes, const uint64_t* ids, size_t nBoxes,
                                  int myRank, QueryBox<T>* out)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nBoxes) return;
    out[i].box  = boxes[i];
    out[i].id   = ids[i];
    out[i].proc = myRank;
}

}  // namespace detail

// Reusable searcher. Hold ONE of these per interface and call search() every timestep -- the
// device scratch persists, so a sliding re-search costs no allocation.
template<typename T>
class CoarseSearch
{
public:
    // Same contract and same sorted result as coarseSearchHost. All inputs are device pointers;
    // the result lands in a device vector and stays there.
    //
    // Collective over the peer graph: every rank must call it, including ranks contributing
    // nothing.
    void search(const Aabb<T>* d_domainBoxes, const uint64_t* d_domainIds, size_t nDomain,
                const Aabb<T>* d_rangeBoxes, const uint64_t* d_rangeIds, size_t nRange,
                const std::vector<int>& candidatePeers, MPI_Comm comm,
                thrust::device_vector<SearchPair>& d_out, bool enforceSymmetry = true)
    {
        int myRank = 0;
        MPI_Comm_rank(comm, &myRank);
        d_out.clear();

        const std::vector<int> peers = detail::normalizePeers(candidatePeers, myRank);
        const int np = int(peers.size());

        // Envelope first: the peers cannot be asked anything until they know my extent.
        Aabb<T> myEnv = (nRange == 0)
                            ? Aabb<T>::empty()
                            : thrust::reduce(thrust::device, d_rangeBoxes, d_rangeBoxes + nRange,
                                             Aabb<T>::empty(), detail::AabbUnion<T>{});
        const std::vector<Aabb<T>> peerEnv = detail::exchangeEnvelopes(myEnv, peers, comm);

        // ---- route: one flag kernel, one scan, one D2H -------------------------------------
        h_sendCnt_.assign(np, 0);
        int totalSend = 0;
        if (np > 0 && nDomain > 0)
        {
            d_env_.assign(peerEnv.begin(), peerEnv.end());
            const size_t tab = size_t(np) * nDomain;
            d_flags_.resize(tab + 1);
            thrust::fill(thrust::device, d_flags_.end() - 1, d_flags_.end(), uint8_t(0));
            launch(nDomain, [&](int nb, int blk) {
                detail::flagAllPeersKernel<T><<<nb, blk>>>(
                    d_domainBoxes, nDomain, thrust::raw_pointer_cast(d_env_.data()), np,
                    thrust::raw_pointer_cast(d_flags_.data()));
            });
            d_scan_.resize(tab + 1);
            thrust::exclusive_scan(thrust::device, d_flags_.begin(), d_flags_.end(),
                                   d_scan_.begin(), 0);
            d_bounds_.resize(np + 1);
            detail::segmentBoundsKernel<<<1, np + 1>>>(thrust::raw_pointer_cast(d_scan_.data()),
                                                       nDomain, np,
                                                       thrust::raw_pointer_cast(d_bounds_.data()));
            h_bounds_.resize(np + 1);
            cudaMemcpy(h_bounds_.data(), thrust::raw_pointer_cast(d_bounds_.data()),
                       (np + 1) * sizeof(int), cudaMemcpyDeviceToHost);   // the ONE count copy
            for (int i = 0; i < np; ++i) h_sendCnt_[i] = h_bounds_[i + 1] - h_bounds_[i];
            totalSend = h_bounds_[np];
            d_sendAll_.resize(totalSend);
            if (totalSend > 0)
                launch(nDomain, [&](int nb, int blk) {
                    detail::gatherAllPeersKernel<T><<<nb, blk>>>(
                        d_domainBoxes, d_domainIds, nDomain,
                        thrust::raw_pointer_cast(d_flags_.data()),
                        thrust::raw_pointer_cast(d_scan_.data()), np, myRank,
                        thrust::raw_pointer_cast(d_sendAll_.data()));
                });
        }
        else { d_sendAll_.clear(); h_bounds_.assign(np + 1, 0); }

        // Post the count round trip and put the self-test in its shadow: my own boxes against my
        // own range boxes needs no message from anyone, so it is free latency-hiding.
        std::vector<MPI_Request> cntRq;
        detail::postCounts(h_sendCnt_, h_recvCnt_, peers, comm, 0x5E01, cntRq);

        d_selfQ_.resize(nDomain);
        if (nDomain > 0)
            launch(nDomain, [&](int nb, int blk) {
                detail::selfQueriesKernel<T><<<nb, blk>>>(
                    d_domainBoxes, d_domainIds, nDomain, myRank,
                    thrust::raw_pointer_cast(d_selfQ_.data()));
            });
        pairInto(thrust::raw_pointer_cast(d_selfQ_.data()), nDomain, d_rangeBoxes, d_rangeIds,
                 nRange, myRank, d_out, nullptr);

        detail::waitCounts(cntRq);
        h_recvOff_.assign(np + 1, 0);
        for (int i = 0; i < np; ++i) h_recvOff_[i + 1] = h_recvOff_[i] + h_recvCnt_[i];
        const int totalRecv = h_recvOff_[np];
        d_recvAll_.resize(totalRecv);

        // ---- post, then do the self-test WHILE the queries are in flight -------------------
        std::vector<MPI_Request> rq;
        rq.reserve(2 * np);
        const int kQueryTag = 0x5EA1;
        for (int i = 0; i < np; ++i)
            if (h_recvCnt_[i] > 0)
            {
                rq.emplace_back();
                MPI_Irecv(thrust::raw_pointer_cast(d_recvAll_.data()) + h_recvOff_[i],
                          h_recvCnt_[i] * int(sizeof(QueryBox<T>)), MPI_BYTE, peers[i], kQueryTag,
                          comm, &rq.back());
            }
        for (int i = 0; i < np; ++i)
            if (h_sendCnt_[i] > 0)
            {
                rq.emplace_back();
                MPI_Isend(thrust::raw_pointer_cast(d_sendAll_.data()) + h_bounds_[i],
                          h_sendCnt_[i] * int(sizeof(QueryBox<T>)), MPI_BYTE, peers[i], kQueryTag,
                          comm, &rq.back());
            }

        if (!rq.empty()) MPI_Waitall(int(rq.size()), rq.data(), MPI_STATUSES_IGNORE);

        // ---- one pass over ALL incoming queries, per-peer counts from the same scan ---------
        h_repSendCnt_.assign(np, 0);
        d_pairsAll_.clear();
        if (totalRecv > 0 && nRange > 0)
        {
            d_at_.assign(h_recvOff_.begin(), h_recvOff_.end());
            pairInto(thrust::raw_pointer_cast(d_recvAll_.data()), totalRecv, d_rangeBoxes,
                     d_rangeIds, nRange, myRank, d_pairsAll_, &h_repSendCnt_, np);
        }
        h_repOff_.assign(np + 1, 0);
        for (int i = 0; i < np; ++i) h_repOff_[i + 1] = h_repOff_[i] + h_repSendCnt_[i];

        // Symmetry: I am the range owner for every pair I just made, so I keep them too.
        if (enforceSymmetry && !d_pairsAll_.empty())
        {
            const size_t at = d_out.size();
            d_out.resize(at + d_pairsAll_.size());
            thrust::copy(d_pairsAll_.begin(), d_pairsAll_.end(), d_out.begin() + at);
        }

        const std::vector<int> repRecvCnt = detail::exchangeCounts(h_repSendCnt_, peers, comm,
                                                                   0x5E02);
        std::vector<int> repRecvOff(np + 1, 0);
        for (int i = 0; i < np; ++i) repRecvOff[i + 1] = repRecvOff[i] + repRecvCnt[i];
        const int totalRepRecv = repRecvOff[np];

        const size_t base = d_out.size();
        d_out.resize(base + totalRepRecv);
        std::vector<MPI_Request> rq2;
        rq2.reserve(2 * np);
        const int kReplyTag = 0x5EA2;
        for (int i = 0; i < np; ++i)
            if (repRecvCnt[i] > 0)
            {
                rq2.emplace_back();
                MPI_Irecv(thrust::raw_pointer_cast(d_out.data()) + base + repRecvOff[i],
                          repRecvCnt[i] * int(sizeof(SearchPair)), MPI_BYTE, peers[i], kReplyTag,
                          comm, &rq2.back());
            }
        for (int i = 0; i < np; ++i)
            if (h_repSendCnt_[i] > 0)
            {
                rq2.emplace_back();
                MPI_Isend(thrust::raw_pointer_cast(d_pairsAll_.data()) + h_repOff_[i],
                          h_repSendCnt_[i] * int(sizeof(SearchPair)), MPI_BYTE, peers[i],
                          kReplyTag, comm, &rq2.back());
            }
        if (!rq2.empty()) MPI_Waitall(int(rq2.size()), rq2.data(), MPI_STATUSES_IGNORE);

        thrust::sort(thrust::device, d_out.begin(), d_out.end());
        d_out.erase(thrust::unique(thrust::device, d_out.begin(), d_out.end()), d_out.end());
    }

private:
    static constexpr int kBlock = 256;

    template<typename F>
    static void launch(size_t n, F&& f)
    {
        if (n == 0) return;
        f(int((n + kBlock - 1) / kBlock), kBlock);
    }

    // Count/scan/fill over one query array. When `perSeg` is given, the per-segment pair counts
    // come out of the SAME scan -- no second pass, no extra copy per peer.
    void pairInto(const QueryBox<T>* d_q, size_t nQ, const Aabb<T>* d_range,
                  const uint64_t* d_rangeIds, size_t nRange, int myRank,
                  thrust::device_vector<SearchPair>& d_dst, std::vector<int>* perSeg,
                  int nSeg = 0)
    {
        if (nQ == 0 || nRange == 0) return;
        d_cnt_.resize(nQ + 1);
        thrust::fill(thrust::device, d_cnt_.end() - 1, d_cnt_.end(), 0);
        launch(nQ, [&](int nb, int blk) {
            detail::countPairsKernel<T><<<nb, blk>>>(d_q, nQ, d_range, nRange,
                                                     thrust::raw_pointer_cast(d_cnt_.data()));
        });
        thrust::exclusive_scan(thrust::device, d_cnt_.begin(), d_cnt_.end(), d_cnt_.begin(), 0);

        int total = 0;
        if (perSeg)
        {
            d_segv_.resize(nSeg + 1);
            detail::gatherAtKernel<<<1, nSeg + 1>>>(thrust::raw_pointer_cast(d_cnt_.data()),
                                                    thrust::raw_pointer_cast(d_at_.data()),
                                                    nSeg + 1,
                                                    thrust::raw_pointer_cast(d_segv_.data()));
            std::vector<int> h(nSeg + 1);
            cudaMemcpy(h.data(), thrust::raw_pointer_cast(d_segv_.data()), (nSeg + 1) * sizeof(int),
                       cudaMemcpyDeviceToHost);   // one copy for every segment count
            for (int i = 0; i < nSeg; ++i) (*perSeg)[i] = h[i + 1] - h[i];
            total = h[nSeg];
        }
        else
        {
            cudaMemcpy(&total, thrust::raw_pointer_cast(d_cnt_.data()) + nQ, sizeof(int),
                       cudaMemcpyDeviceToHost);
        }
        if (total == 0) return;

        const size_t at = d_dst.size();
        d_dst.resize(at + total);
        launch(nQ, [&](int nb, int blk) {
            detail::fillPairsKernel<T><<<nb, blk>>>(
                d_q, nQ, d_range, d_rangeIds, nRange, thrust::raw_pointer_cast(d_cnt_.data()),
                myRank, thrust::raw_pointer_cast(d_dst.data()) + at);
        });
    }

    // Persistent scratch -- the point of the class.
    thrust::device_vector<Aabb<T>>     d_env_;
    thrust::device_vector<uint8_t>     d_flags_;
    thrust::device_vector<int>         d_scan_, d_bounds_, d_cnt_, d_at_, d_segv_;
    thrust::device_vector<QueryBox<T>> d_sendAll_, d_recvAll_, d_selfQ_;
    thrust::device_vector<SearchPair>  d_pairsAll_;
    std::vector<int> h_sendCnt_, h_recvCnt_, h_recvOff_, h_bounds_, h_repSendCnt_, h_repOff_;
};

// Convenience for one-shot callers. A per-timestep consumer should hold a CoarseSearch<T> instead,
// so the scratch survives between calls.
template<typename T>
void coarseSearch(const Aabb<T>* d_domainBoxes, const uint64_t* d_domainIds, size_t nDomain,
                  const Aabb<T>* d_rangeBoxes, const uint64_t* d_rangeIds, size_t nRange,
                  const std::vector<int>& candidatePeers, MPI_Comm comm,
                  thrust::device_vector<SearchPair>& d_out, bool enforceSymmetry = true)
{
    CoarseSearch<T> cs;
    cs.search(d_domainBoxes, d_domainIds, nDomain, d_rangeBoxes, d_rangeIds, nRange, candidatePeers,
              comm, d_out, enforceSymmetry);
}

#endif  // MARS_CS_CUDA

}  // namespace fem
}  // namespace mars

#endif  // MARS_COARSE_SEARCH_HPP
