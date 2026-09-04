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
#include <thrust/iterator/counting_iterator.h>
#include <thrust/scan.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <thrust/transform_reduce.h>
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
                      std::vector<SearchPair>& out)
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
// ---------------------------------------------------------------------------
#ifdef MARS_CS_CUDA

namespace detail
{

template<typename T>
__global__ void flagOverlapKernel(const Aabb<T>* boxes, size_t nBoxes, Aabb<T> env, uint8_t* flags)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nBoxes) return;
    flags[i] = boxes[i].overlaps(env);
}

template<typename T>
__global__ void gatherQueriesKernel(const Aabb<T>* boxes, const uint64_t* ids, size_t nBoxes,
                                    const uint8_t* flags, const int* scan, int myRank,
                                    QueryBox<T>* outQ)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nBoxes || !flags[i]) return;
    QueryBox<T> q;
    q.box         = boxes[i];
    q.id          = ids[i];
    q.proc        = myRank;
    outQ[scan[i]] = q;
}

// Two-pass local test. Counting first lets pass 2 write without atomics, so the pair order depends
// only on the input order and not on scheduling -- which is what lets the device result be diffed
// against the host one rather than merely compared as a set.
template<typename T>
__global__ void countPairsKernel(const QueryBox<T>* q, size_t nQ, const Aabb<T>* range,
                                 size_t nRange, int* counts)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nQ) return;
    int c = 0;
    for (size_t j = 0; j < nRange; ++j)
        if (q[i].box.overlaps(range[j])) ++c;
    counts[i] = c;
}

template<typename T>
__global__ void fillPairsKernel(const QueryBox<T>* q, size_t nQ, const Aabb<T>* range,
                                const uint64_t* rangeIds, size_t nRange, const int* offsets,
                                int myRank, SearchPair* out)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nQ) return;
    int w = offsets[i];
    for (size_t j = 0; j < nRange; ++j)
        if (q[i].box.overlaps(range[j]))
        {
            out[w].domain = IdentProc{q[i].id, q[i].proc};
            out[w].range  = IdentProc{rangeIds[j], myRank};
            ++w;
        }
}

// One device pass: every query in `d_q` against every range box, compacted. Used for both the
// local self-test and each peer's incoming queries.
template<typename T>
inline void pairQueriesOnDevice(const QueryBox<T>* d_q, size_t nQ, const Aabb<T>* d_range,
                                const uint64_t* d_rangeIds, size_t nRange, int myRank,
                                thrust::device_vector<SearchPair>& d_pairs)
{
    d_pairs.clear();
    if (nQ == 0 || nRange == 0) return;
    const int kBlock = 256;
    const int nb     = int((nQ + kBlock - 1) / kBlock);
    thrust::device_vector<int> d_cnt(nQ + 1, 0);
    countPairsKernel<T>
        <<<nb, kBlock>>>(d_q, nQ, d_range, nRange, thrust::raw_pointer_cast(d_cnt.data()));
    cudaDeviceSynchronize();
    thrust::exclusive_scan(thrust::device, d_cnt.begin(), d_cnt.end(), d_cnt.begin(), 0);
    const int total = d_cnt[nQ];
    if (total == 0) return;
    d_pairs.resize(total);
    fillPairsKernel<T><<<nb, kBlock>>>(d_q, nQ, d_range, d_rangeIds, nRange,
                                       thrust::raw_pointer_cast(d_cnt.data()), myRank,
                                       thrust::raw_pointer_cast(d_pairs.data()));
    cudaDeviceSynchronize();
}

}  // namespace detail

// Device coarse search. All inputs are device pointers; the result lands in a device vector and
// stays there. Same contract as coarseSearchHost, same sorted result.
//
// Collective over the peer graph: every rank must call it, including ranks contributing nothing.
template<typename T>
void coarseSearch(const Aabb<T>* d_domainBoxes, const uint64_t* d_domainIds, size_t nDomain,
                  const Aabb<T>* d_rangeBoxes, const uint64_t* d_rangeIds, size_t nRange,
                  const std::vector<int>& candidatePeers, MPI_Comm comm,
                  thrust::device_vector<SearchPair>& d_out)
{
    int myRank = 0;
    MPI_Comm_rank(comm, &myRank);
    d_out.clear();
    const int kBlock = 256;

    const std::vector<int> peers = detail::normalizePeers(candidatePeers, myRank);
    const int np = int(peers.size());

    // Self-test first: my own domain boxes become QueryBoxes pointing at me, no MPI involved.
    thrust::device_vector<QueryBox<T>> d_selfQ(nDomain);
    if (nDomain > 0)
    {
        thrust::device_vector<uint8_t> d_all(nDomain, 1);
        thrust::device_vector<int> d_iota(nDomain);
        thrust::sequence(thrust::device, d_iota.begin(), d_iota.end(), 0);
        int nb = int((nDomain + kBlock - 1) / kBlock);
        detail::gatherQueriesKernel<T><<<nb, kBlock>>>(
            d_domainBoxes, d_domainIds, nDomain, thrust::raw_pointer_cast(d_all.data()),
            thrust::raw_pointer_cast(d_iota.data()), myRank,
            thrust::raw_pointer_cast(d_selfQ.data()));
        cudaDeviceSynchronize();
    }
    detail::pairQueriesOnDevice<T>(thrust::raw_pointer_cast(d_selfQ.data()), nDomain, d_rangeBoxes,
                                   d_rangeIds, nRange, myRank, d_out);

    // My range envelope, reduced on device. One struct per peer is the only geometry that ever
    // touches the host.
    Aabb<T> myEnv = thrust::transform_reduce(
        thrust::device, thrust::counting_iterator<size_t>(0),
        thrust::counting_iterator<size_t>(nRange),
        [d_rangeBoxes] __device__(size_t i) { return d_rangeBoxes[i]; }, Aabb<T>::empty(),
        [] __host__ __device__(const Aabb<T>& a, const Aabb<T>& b) {
            Aabb<T> r = a;
            r.extend(b);
            return r;
        });
    const std::vector<Aabb<T>> peerEnv = detail::exchangeEnvelopes(myEnv, peers, comm);

    // Route, on device, one peer at a time.
    std::vector<thrust::device_vector<QueryBox<T>>> d_send(np);
    std::vector<int> sendCnt(np, 0);
    thrust::device_vector<uint8_t> d_flags(nDomain);
    thrust::device_vector<int> d_scan(nDomain);
    for (int i = 0; i < np; ++i)
    {
        if (nDomain == 0) continue;
        int nb = int((nDomain + kBlock - 1) / kBlock);
        detail::flagOverlapKernel<T><<<nb, kBlock>>>(d_domainBoxes, nDomain, peerEnv[i],
                                                     thrust::raw_pointer_cast(d_flags.data()));
        cudaDeviceSynchronize();
        thrust::exclusive_scan(thrust::device, d_flags.begin(), d_flags.end(), d_scan.begin(), 0);
        // Exclusive scan, so the total is the last offset plus the last flag. Two 4-byte reads per
        // peer; unavoidable, because the MPI count has to be on the host.
        int last     = d_scan[nDomain - 1];
        int lastFlag = int(d_flags[nDomain - 1]);
        sendCnt[i]   = last + lastFlag;
        if (sendCnt[i] == 0) continue;
        d_send[i].resize(sendCnt[i]);
        detail::gatherQueriesKernel<T><<<nb, kBlock>>>(
            d_domainBoxes, d_domainIds, nDomain, thrust::raw_pointer_cast(d_flags.data()),
            thrust::raw_pointer_cast(d_scan.data()), myRank,
            thrust::raw_pointer_cast(d_send[i].data()));
        cudaDeviceSynchronize();
    }
    const std::vector<int> recvCnt = detail::exchangeCounts(sendCnt, peers, comm, 0x5E01);

    std::vector<thrust::device_vector<QueryBox<T>>> d_recv(np);
    std::vector<MPI_Request> rq;
    rq.reserve(2 * np);
    const int kQueryTag = 0x5EA1;
    for (int i = 0; i < np; ++i)
        if (recvCnt[i] > 0)
        {
            d_recv[i].resize(recvCnt[i]);
            rq.emplace_back();
            MPI_Irecv(thrust::raw_pointer_cast(d_recv[i].data()),
                      recvCnt[i] * int(sizeof(QueryBox<T>)), MPI_BYTE, peers[i], kQueryTag, comm,
                      &rq.back());
        }
    for (int i = 0; i < np; ++i)
        if (sendCnt[i] > 0)
        {
            rq.emplace_back();
            MPI_Isend(thrust::raw_pointer_cast(d_send[i].data()),
                      sendCnt[i] * int(sizeof(QueryBox<T>)), MPI_BYTE, peers[i], kQueryTag, comm,
                      &rq.back());
        }
    if (!rq.empty()) MPI_Waitall(int(rq.size()), rq.data(), MPI_STATUSES_IGNORE);

    std::vector<thrust::device_vector<SearchPair>> d_reply(np);
    std::vector<int> repSendCnt(np, 0);
    for (int i = 0; i < np; ++i)
    {
        detail::pairQueriesOnDevice<T>(thrust::raw_pointer_cast(d_recv[i].data()), d_recv[i].size(),
                                       d_rangeBoxes, d_rangeIds, nRange, myRank, d_reply[i]);
        repSendCnt[i] = int(d_reply[i].size());
    }
    const std::vector<int> repRecvCnt = detail::exchangeCounts(repSendCnt, peers, comm, 0x5E02);

    std::vector<thrust::device_vector<SearchPair>> d_repRecv(np);
    std::vector<MPI_Request> rq2;
    rq2.reserve(2 * np);
    const int kReplyTag = 0x5EA2;
    for (int i = 0; i < np; ++i)
        if (repRecvCnt[i] > 0)
        {
            d_repRecv[i].resize(repRecvCnt[i]);
            rq2.emplace_back();
            MPI_Irecv(thrust::raw_pointer_cast(d_repRecv[i].data()),
                      repRecvCnt[i] * int(sizeof(SearchPair)), MPI_BYTE, peers[i], kReplyTag, comm,
                      &rq2.back());
        }
    for (int i = 0; i < np; ++i)
        if (repSendCnt[i] > 0)
        {
            rq2.emplace_back();
            MPI_Isend(thrust::raw_pointer_cast(d_reply[i].data()),
                      repSendCnt[i] * int(sizeof(SearchPair)), MPI_BYTE, peers[i], kReplyTag, comm,
                      &rq2.back());
        }
    if (!rq2.empty()) MPI_Waitall(int(rq2.size()), rq2.data(), MPI_STATUSES_IGNORE);

    size_t off   = d_out.size();
    size_t total = off;
    for (int i = 0; i < np; ++i) total += d_repRecv[i].size();
    d_out.resize(total);
    for (int i = 0; i < np; ++i)
    {
        if (d_repRecv[i].empty()) continue;
        thrust::copy(d_repRecv[i].begin(), d_repRecv[i].end(), d_out.begin() + off);
        off += d_repRecv[i].size();
    }
    thrust::sort(thrust::device, d_out.begin(), d_out.end());
    d_out.erase(thrust::unique(thrust::device, d_out.begin(), d_out.end()), d_out.end());
}

#endif  // MARS_CS_CUDA

}  // namespace fem
}  // namespace mars

#endif  // MARS_COARSE_SEARCH_HPP
