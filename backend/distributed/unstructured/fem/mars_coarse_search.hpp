#ifndef MARS_COARSE_SEARCH_HPP
#define MARS_COARSE_SEARCH_HPP

// Parallel AABB-overlap search -- the MARS equivalent of stk::search::coarse_search.
//
// This is the geometric pairing step that non-conformal interfaces and overset need: given
// DOMAIN boxes (the queries, e.g. slave-side face AABBs) and RANGE boxes (the targets, e.g.
// master-side face AABBs), both distributed, return every overlapping (domain, range) pair
// with the owning rank of each side.
//
// TWO implementations, same semantics and same result:
//   coarseSearch        host + MPI. The oracle. It exists so the algorithm can be pinned down
//                       with mpirun on a laptop -- no GPU, no mesh file -- and so the device
//                       version has something to be compared against element-for-element.
//   coarseSearchDevice  the production path (USE_CUDA). Boxes, routed queries and pairs all
//                       stay in device memory; MPI moves device pointers directly. The only
//                       host traffic is the per-rank envelope and the MPI message counts,
//                       which have to be on the host because MPI counts are.
//
// Why envelopes and not a distributed tree: at the coarse level STK is doing the same thing --
// decide which ranks could possibly hold a partner, ship the query there, test locally. One
// AABB per rank is O(P) metadata, and for interface/overset work the boxes live on a surface,
// so the envelopes are small and the fan-out is a handful of ranks. A distributed octree buys
// nothing until the query set stops being a surface.

#include <mpi.h>

#ifdef USE_CUDA
#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/scan.h>
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
            lo[d] = std::min(lo[d], o.lo[d]);
            hi[d] = std::max(hi[d], o.hi[d]);
        }
    }

    // Closed test: boxes that only touch count as overlapping. Interface faces on a conformal
    // seam meet exactly at their shared edge, and dropping those would lose real partners.
    MARS_CS_FN bool overlaps(const Aabb& o) const
    {
        for (int d = 0; d < 3; ++d)
            if (hi[d] < o.lo[d] || o.hi[d] < lo[d]) return false;
        return true;
    }
};

// STK calls this IdentProc: which entity, and which rank owns it. The rank matters because the
// partner usually lives somewhere else, and whoever consumes the pair has to know where to go
// for the data -- that is the handoff to the ghosting.
struct IdentProc
{
    uint64_t id;
    int      proc;

    MARS_CS_FN bool operator==(const IdentProc& o) const { return id == o.id && proc == o.proc; }
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
    MARS_CS_FN bool operator==(const SearchPair& o) const { return domain == o.domain && range == o.range; }
};

MARS_CS_FN inline bool operator!=(const IdentProc& a, const IdentProc& b) { return !(a == b); }

namespace detail
{

// One AABB per rank over that rank's range boxes. A rank holding none reports the empty box,
// which overlaps nothing -- so it is skipped without needing a separate "has data" flag.
template<typename T>
inline std::vector<Aabb<T>> gatherRangeEnvelopes(const std::vector<Aabb<T>>& range, MPI_Comm comm)
{
    int nRanks = 0;
    MPI_Comm_size(comm, &nRanks);

    Aabb<T> local = Aabb<T>::empty();
    for (const auto& b : range) local.extend(b);

    std::vector<Aabb<T>> all(nRanks);
    MPI_Allgather(&local, sizeof(Aabb<T>), MPI_BYTE, all.data(), sizeof(Aabb<T>), MPI_BYTE, comm);
    return all;
}

}  // namespace detail

// One query box travelling to a rank that might hold a partner: the box, plus enough identity
// to send the answer home.
template<typename T>
struct QueryBox
{
    Aabb<T>  box;
    uint64_t id;
    int      proc;
};

// Every rank supplies its own domain and range boxes with matching id arrays. On return, `out`
// holds the pairs whose DOMAIN box this rank owns, sorted and deduplicated -- so the result is
// deterministic and independent of rank count, which is what makes it testable.
//
// Collective: every rank must call it, including ranks with nothing to contribute.
template<typename T>
void coarseSearch(const std::vector<Aabb<T>>& domainBoxes,
                  const std::vector<uint64_t>& domainIds,
                  const std::vector<Aabb<T>>& rangeBoxes,
                  const std::vector<uint64_t>& rangeIds,
                  MPI_Comm comm,
                  std::vector<SearchPair>& out)
{
    int myRank = 0, nRanks = 0;
    MPI_Comm_rank(comm, &myRank);
    MPI_Comm_size(comm, &nRanks);
    out.clear();

    const std::vector<Aabb<T>> envelopes = detail::gatherRangeEnvelopes(rangeBoxes, comm);

    // Which of my domain boxes goes to which rank. A box overlapping several envelopes is sent
    // to each of them; the dedup at the end removes any pair found twice.
    std::vector<std::vector<QueryBox<T>>> toSend(nRanks);
    for (size_t i = 0; i < domainBoxes.size(); ++i)
        for (int r = 0; r < nRanks; ++r)
            if (domainBoxes[i].overlaps(envelopes[r]))
                toSend[r].push_back(QueryBox<T>{domainBoxes[i], domainIds[i], myRank});

    std::vector<int> sendCnt(nRanks, 0), recvCnt(nRanks, 0);
    for (int r = 0; r < nRanks; ++r) sendCnt[r] = static_cast<int>(toSend[r].size());
    MPI_Alltoall(sendCnt.data(), 1, MPI_INT, recvCnt.data(), 1, MPI_INT, comm);

    // Ship the queries. Non-blocking both ways; a rank with nothing to say still posts nothing
    // and still reaches the Waitall, so no one is left hanging.
    std::vector<std::vector<QueryBox<T>>> recvBuf(nRanks);
    std::vector<MPI_Request> reqs;
    reqs.reserve(2 * nRanks);
    const int kQueryTag = 0x5EA1;
    for (int r = 0; r < nRanks; ++r)
    {
        if (recvCnt[r] > 0)
        {
            recvBuf[r].resize(recvCnt[r]);
            reqs.emplace_back();
            MPI_Irecv(recvBuf[r].data(), recvCnt[r] * int(sizeof(QueryBox<T>)), MPI_BYTE, r,
                      kQueryTag, comm, &reqs.back());
        }
    }
    for (int r = 0; r < nRanks; ++r)
    {
        if (sendCnt[r] > 0)
        {
            reqs.emplace_back();
            MPI_Isend(toSend[r].data(), sendCnt[r] * int(sizeof(QueryBox<T>)), MPI_BYTE, r,
                      kQueryTag, comm, &reqs.back());
        }
    }
    if (!reqs.empty()) MPI_Waitall(int(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);

    // Local test: every query that landed here against every range box I own. This brute-force
    // inner loop is the piece a device tree replaces later; the communication above does not
    // change when it does.
    std::vector<std::vector<SearchPair>> reply(nRanks);
    for (int r = 0; r < nRanks; ++r)
        for (const auto& q : recvBuf[r])
            for (size_t j = 0; j < rangeBoxes.size(); ++j)
                if (q.box.overlaps(rangeBoxes[j]))
                    reply[q.proc].push_back(
                        SearchPair{IdentProc{q.id, q.proc}, IdentProc{rangeIds[j], myRank}});

    // Send each pair back to the rank that owns its domain box, so a caller only ever sees
    // pairs it asked for.
    std::vector<int> repSendCnt(nRanks, 0), repRecvCnt(nRanks, 0);
    for (int r = 0; r < nRanks; ++r) repSendCnt[r] = static_cast<int>(reply[r].size());
    MPI_Alltoall(repSendCnt.data(), 1, MPI_INT, repRecvCnt.data(), 1, MPI_INT, comm);

    std::vector<std::vector<SearchPair>> repRecv(nRanks);
    std::vector<MPI_Request> reqs2;
    reqs2.reserve(2 * nRanks);
    const int kReplyTag = 0x5EA2;
    for (int r = 0; r < nRanks; ++r)
    {
        if (repRecvCnt[r] > 0)
        {
            repRecv[r].resize(repRecvCnt[r]);
            reqs2.emplace_back();
            MPI_Irecv(repRecv[r].data(), repRecvCnt[r] * int(sizeof(SearchPair)), MPI_BYTE, r,
                      kReplyTag, comm, &reqs2.back());
        }
    }
    for (int r = 0; r < nRanks; ++r)
    {
        if (repSendCnt[r] > 0)
        {
            reqs2.emplace_back();
            MPI_Isend(reply[r].data(), repSendCnt[r] * int(sizeof(SearchPair)), MPI_BYTE, r,
                      kReplyTag, comm, &reqs2.back());
        }
    }
    if (!reqs2.empty()) MPI_Waitall(int(reqs2.size()), reqs2.data(), MPI_STATUSES_IGNORE);

    for (int r = 0; r < nRanks; ++r) out.insert(out.end(), repRecv[r].begin(), repRecv[r].end());
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
}

// ---------------------------------------------------------------------------
// Device path. Same semantics as coarseSearch above, same result -- that host
// version is the oracle the GPU one is tested against, which is the only reason
// it exists. Everything here stays in device memory: the boxes, the routed
// queries, the pairs. MPI moves device pointers directly (GPUDirect), so a pair
// never touches the host between the two ends.
// ---------------------------------------------------------------------------
#ifdef USE_CUDA

namespace detail
{

// Per-(box, rank) overlap flags against the gathered envelopes: flag[r*n + i] tells whether my
// domain box i has to travel to rank r. One kernel for the whole P x n table, so the routing
// decision never leaves the device.
template<typename T>
__global__ void flagQueryTargetsKernel(const Aabb<T>* boxes, size_t nBoxes,
                                       const Aabb<T>* envelopes, int nRanks, uint8_t* flags)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nBoxes) return;
    for (int r = 0; r < nRanks; ++r) flags[size_t(r) * nBoxes + i] = boxes[i].overlaps(envelopes[r]);
}

template<typename T>
__global__ void gatherQueriesKernel(const Aabb<T>* boxes, const uint64_t* ids, size_t nBoxes,
                                    const uint8_t* flags, const int* scan, int myRank,
                                    QueryBox<T>* outQ)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nBoxes || !flags[i]) return;
    QueryBox<T> q;
    q.box  = boxes[i];
    q.id   = ids[i];
    q.proc = myRank;
    outQ[scan[i]] = q;
}

// Two-pass local test. Pass 1 counts so pass 2 can write without atomics into a compacted
// array -- the pair order then depends only on the input order, not on scheduling, which is
// what lets the device result be compared against the host one.
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

}  // namespace detail

// Device coarse search. All inputs are device pointers; the result lands in a device vector and
// stays there. Collective -- every rank calls it, including ranks contributing nothing.
template<typename T>
void coarseSearchDevice(const Aabb<T>* d_domainBoxes, const uint64_t* d_domainIds, size_t nDomain,
                        const Aabb<T>* d_rangeBoxes, const uint64_t* d_rangeIds, size_t nRange,
                        MPI_Comm comm, thrust::device_vector<SearchPair>& d_out)
{
    int myRank = 0, nRanks = 0;
    MPI_Comm_rank(comm, &myRank);
    MPI_Comm_size(comm, &nRanks);
    d_out.clear();
    const int kBlock = 256;

    // Envelope of my range boxes, reduced on device. Only this one small struct per rank goes
    // through the host, because MPI_Allgather of P structs is not worth a device path.
    Aabb<T> localEnv = thrust::transform_reduce(
        thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(nRange),
        [d_rangeBoxes] __device__(size_t i) { return d_rangeBoxes[i]; }, Aabb<T>::empty(),
        [] __host__ __device__(const Aabb<T>& a, const Aabb<T>& b) {
            Aabb<T> r = a;
            for (int d = 0; d < 3; ++d)
            {
                r.lo[d] = a.lo[d] < b.lo[d] ? a.lo[d] : b.lo[d];
                r.hi[d] = a.hi[d] > b.hi[d] ? a.hi[d] : b.hi[d];
            }
            return r;
        });

    std::vector<Aabb<T>> hEnv(nRanks);
    MPI_Allgather(&localEnv, sizeof(Aabb<T>), MPI_BYTE, hEnv.data(), sizeof(Aabb<T>), MPI_BYTE, comm);
    thrust::device_vector<Aabb<T>> d_env(hEnv.begin(), hEnv.end());

    // Route: flag, scan, gather -- all on device, one send buffer per destination rank.
    thrust::device_vector<uint8_t> d_flags(size_t(nRanks) * nDomain, 0);
    if (nDomain > 0)
    {
        int nb = int((nDomain + kBlock - 1) / kBlock);
        detail::flagQueryTargetsKernel<T><<<nb, kBlock>>>(
            d_domainBoxes, nDomain, thrust::raw_pointer_cast(d_env.data()), nRanks,
            thrust::raw_pointer_cast(d_flags.data()));
        cudaDeviceSynchronize();
    }

    std::vector<thrust::device_vector<QueryBox<T>>> d_send(nRanks);
    std::vector<int> sendCnt(nRanks, 0), recvCnt(nRanks, 0);
    thrust::device_vector<int> d_scan(nDomain);
    for (int r = 0; r < nRanks; ++r)
    {
        if (nDomain == 0) continue;
        const uint8_t* fr = thrust::raw_pointer_cast(d_flags.data()) + size_t(r) * nDomain;
        thrust::exclusive_scan(thrust::device, fr, fr + nDomain, d_scan.begin(), 0);
        int last = 0, lastFlag = 0;
        cudaMemcpy(&last, thrust::raw_pointer_cast(d_scan.data()) + nDomain - 1, sizeof(int),
                   cudaMemcpyDeviceToHost);
        // The scan is exclusive, so the total needs the final flag added back. Two 4-byte reads
        // per rank -- the only D2H in the routing, and unavoidable: MPI counts live on the host.
        uint8_t lf = 0;
        cudaMemcpy(&lf, fr + nDomain - 1, sizeof(uint8_t), cudaMemcpyDeviceToHost);
        lastFlag   = lf;
        sendCnt[r] = last + lastFlag;
        if (sendCnt[r] == 0) continue;
        d_send[r].resize(sendCnt[r]);
        int nb = int((nDomain + kBlock - 1) / kBlock);
        detail::gatherQueriesKernel<T><<<nb, kBlock>>>(
            d_domainBoxes, d_domainIds, nDomain, fr, thrust::raw_pointer_cast(d_scan.data()),
            myRank, thrust::raw_pointer_cast(d_send[r].data()));
    }
    cudaDeviceSynchronize();
    MPI_Alltoall(sendCnt.data(), 1, MPI_INT, recvCnt.data(), 1, MPI_INT, comm);

    // Queries move device-to-device.
    std::vector<thrust::device_vector<QueryBox<T>>> d_recv(nRanks);
    std::vector<MPI_Request> reqs;
    reqs.reserve(2 * nRanks);
    const int kQueryTag = 0x5EA1;
    for (int r = 0; r < nRanks; ++r)
        if (recvCnt[r] > 0)
        {
            d_recv[r].resize(recvCnt[r]);
            reqs.emplace_back();
            MPI_Irecv(thrust::raw_pointer_cast(d_recv[r].data()),
                      recvCnt[r] * int(sizeof(QueryBox<T>)), MPI_BYTE, r, kQueryTag, comm,
                      &reqs.back());
        }
    for (int r = 0; r < nRanks; ++r)
        if (sendCnt[r] > 0)
        {
            reqs.emplace_back();
            MPI_Isend(thrust::raw_pointer_cast(d_send[r].data()),
                      sendCnt[r] * int(sizeof(QueryBox<T>)), MPI_BYTE, r, kQueryTag, comm,
                      &reqs.back());
        }
    if (!reqs.empty()) MPI_Waitall(int(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);

    // Local test per source rank, then the pairs go home.
    std::vector<thrust::device_vector<SearchPair>> d_reply(nRanks);
    std::vector<int> repSendCnt(nRanks, 0), repRecvCnt(nRanks, 0);
    for (int r = 0; r < nRanks; ++r)
    {
        const size_t nQ = d_recv[r].size();
        if (nQ == 0 || nRange == 0) continue;
        thrust::device_vector<int> d_cnt(nQ + 1, 0);
        int nb = int((nQ + kBlock - 1) / kBlock);
        detail::countPairsKernel<T><<<nb, kBlock>>>(thrust::raw_pointer_cast(d_recv[r].data()), nQ,
                                                    d_rangeBoxes, nRange,
                                                    thrust::raw_pointer_cast(d_cnt.data()));
        cudaDeviceSynchronize();
        thrust::exclusive_scan(thrust::device, d_cnt.begin(), d_cnt.end(), d_cnt.begin(), 0);
        int total = d_cnt[nQ];
        if (total == 0) continue;
        d_reply[r].resize(total);
        detail::fillPairsKernel<T><<<nb, kBlock>>>(thrust::raw_pointer_cast(d_recv[r].data()), nQ,
                                                   d_rangeBoxes, d_rangeIds, nRange,
                                                   thrust::raw_pointer_cast(d_cnt.data()), myRank,
                                                   thrust::raw_pointer_cast(d_reply[r].data()));
        repSendCnt[r] = total;
    }
    cudaDeviceSynchronize();
    MPI_Alltoall(repSendCnt.data(), 1, MPI_INT, repRecvCnt.data(), 1, MPI_INT, comm);

    std::vector<thrust::device_vector<SearchPair>> d_repRecv(nRanks);
    std::vector<MPI_Request> reqs2;
    reqs2.reserve(2 * nRanks);
    const int kReplyTag = 0x5EA2;
    for (int r = 0; r < nRanks; ++r)
        if (repRecvCnt[r] > 0)
        {
            d_repRecv[r].resize(repRecvCnt[r]);
            reqs2.emplace_back();
            MPI_Irecv(thrust::raw_pointer_cast(d_repRecv[r].data()),
                      repRecvCnt[r] * int(sizeof(SearchPair)), MPI_BYTE, r, kReplyTag, comm,
                      &reqs2.back());
        }
    for (int r = 0; r < nRanks; ++r)
        if (repSendCnt[r] > 0)
        {
            reqs2.emplace_back();
            MPI_Isend(thrust::raw_pointer_cast(d_reply[r].data()),
                      repSendCnt[r] * int(sizeof(SearchPair)), MPI_BYTE, r, kReplyTag, comm,
                      &reqs2.back());
        }
    if (!reqs2.empty()) MPI_Waitall(int(reqs2.size()), reqs2.data(), MPI_STATUSES_IGNORE);

    size_t total = 0;
    for (int r = 0; r < nRanks; ++r) total += d_repRecv[r].size();
    d_out.resize(total);
    size_t off = 0;
    for (int r = 0; r < nRanks; ++r)
    {
        if (d_repRecv[r].empty()) continue;
        thrust::copy(d_repRecv[r].begin(), d_repRecv[r].end(), d_out.begin() + off);
        off += d_repRecv[r].size();
    }
    // Sorted + deduped on device, so the device result is comparable to the host one
    // element-for-element rather than only as a set.
    thrust::sort(thrust::device, d_out.begin(), d_out.end());
    d_out.erase(thrust::unique(thrust::device, d_out.begin(), d_out.end()), d_out.end());
}

#endif  // USE_CUDA

}  // namespace fem
}  // namespace mars

#endif  // MARS_COARSE_SEARCH_HPP
