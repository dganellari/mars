#pragma once

#if defined(__CUDACC__) || defined(__HIPCC__)
#include "../utils/mars_gpu_runtime.hpp"
#endif
#include "mars_ghost_registry.hpp"
#include "mars_nonconformal_scalar.hpp"
#include <climits>
#include <type_traits>

namespace mars::fem
{
// This registry imports whole donor elements, including the off-face corner used by the gradient.
// Its metadata setup uses the existing host registry API; geometry stays on the device.
// A scalar payload reuses the persistent exchange for solution ghosts.
template<class RealType, class Payload = InterfaceTetFace<RealType>>
class NonconformalDonorExchange
{
    using Face = Payload;
    static_assert(std::is_trivially_copyable_v<Face>);

public:
    explicit NonconformalDonorExchange(MPI_Comm comm) { MPI_Comm_dup(comm, &comm_); }
    NonconformalDonorExchange(const NonconformalDonorExchange&)            = delete;
    NonconformalDonorExchange& operator=(const NonconformalDonorExchange&) = delete;
    ~NonconformalDonorExchange()
    {
        int finalized = 0;
        MPI_Finalized(&finalized);
        if (!finalized) MPI_Comm_free(&comm_);
    }

    // All ranks participate, including ranks with no interface. Peer lists must be reciprocal.
    bool buildFromHost(const std::vector<uint64_t>& ids, const std::vector<int>& owners, std::vector<int> peers)
    {
        ready_ = false;
        int rank, size;
        MPI_Comm_rank(comm_, &rank);
        MPI_Comm_size(comm_, &size);
        bool valid = ids.size() == owners.size() && ids.size() <= size_t(INT_MAX / sizeof(Face));
        std::sort(peers.begin(), peers.end());
        peers.erase(std::unique(peers.begin(), peers.end()), peers.end());
        peers.erase(std::remove(peers.begin(), peers.end(), rank), peers.end());
        for (int peer : peers)
            valid = valid && peer >= 0 && peer < size;
        auto sorted = ids;
        std::sort(sorted.begin(), sorted.end());
        valid = valid && std::adjacent_find(sorted.begin(), sorted.end()) == sorted.end();
        for (int owner : owners)
            valid = valid && owner >= 0 && owner < size &&
                    (owner == rank || std::binary_search(peers.begin(), peers.end(), owner));
        if (!allValid(valid)) return false;

        // MPI discovers incoming neighbours without a rank-sized all-gather. Check before the
        // registry posts receives: a one-sided peer list would otherwise leave it waiting forever.
        MPI_Comm graph;
        int degree = int(peers.size());
        MPI_Dist_graph_create(comm_, 1, &rank, &degree, peers.empty() ? &rank : peers.data(), MPI_UNWEIGHTED,
                              MPI_INFO_NULL, 0, &graph);
        int incoming, outgoing, weighted;
        MPI_Dist_graph_neighbors_count(graph, &incoming, &outgoing, &weighted);
        std::vector<int> sources(incoming), destinations(outgoing);
        MPI_Dist_graph_neighbors(graph, incoming, sources.empty() ? &rank : sources.data(), MPI_UNWEIGHTED, outgoing,
                                 destinations.empty() ? &rank : destinations.data(), MPI_UNWEIGHTED);
        MPI_Comm_free(&graph);
        std::sort(sources.begin(), sources.end());
        if (!allValid(sources == peers)) return false;

        // A fresh registry also clears unresolved requests from a previous rejected plan.
        registry_ = GhostRegistry<uint64_t>{"nonconformal-donors"};
        registry_.build(ids.size(), owners, ids, std::vector<uint8_t>(ids.size(), 1), rank, peers, comm_);
        valid        = registry_.isConsistent();
        sendOffsets_ = {0};
        recvOffsets_ = {0};
        std::vector<int> send, recv;
        for (size_t p = 0; p < registry_.peers().size(); ++p)
        {
            if (registry_.sendCount(int(p)))
            {
                auto list = registry_.sendList(int(p));
                for (auto i = list.first; i != list.second; ++i)
                {
                    valid = valid && *i >= 0 && size_t(*i) < owners.size();
                    if (*i >= 0 && size_t(*i) < owners.size()) valid = valid && owners[*i] == rank;
                    send.push_back(*i);
                }
            }
            if (registry_.recvCount(int(p)))
            {
                auto list = registry_.recvList(int(p));
                recv.insert(recv.end(), list.first, list.second);
            }
            valid =
                valid && send.size() <= size_t(INT_MAX / sizeof(Face)) && recv.size() <= size_t(INT_MAX / sizeof(Face));
            if (!valid) break;
            sendOffsets_.push_back(int(send.size()));
            recvOffsets_.push_back(int(recv.size()));
        }
        if (!allValid(valid)) return false;
        count_ = ids.size();
#if defined(__CUDACC__) || defined(__HIPCC__)
        try
        {
            d_send_ = send;
            d_recv_ = recv;
            d_out_.resize(send.size());
            d_in_.resize(recv.size());
        }
        catch (...)
        {
            valid = false;
        }
        if (!allValid(valid)) return false;
#endif
        ready_ = true;
        return true;
    }

    bool forwardHost(std::vector<Face>& faces) const
    {
        if (!allValid(ready_ && faces.size() == count_)) return false;
        registry_.forwardHost(faces);
        return true;
    }

#if defined(__CUDACC__) || defined(__HIPCC__)
    // Collective completion boundary: MPI cannot consume a pending GPU gather. Only this stream
    // is synchronized, and buffers are retained between refreshes. Requires GPU-aware MPI.
    bool forward(Face* faces, size_t count, cudaStream_t stream = 0)
    {
        if (!allValid(ready_ && count == count_ && (faces || count == 0))) return false;
        if (!d_send_.empty())
            detail::ghostGatherKernel<<<(d_send_.size() + 255) / 256, 256, 0, stream>>>(
                faces, thrust::raw_pointer_cast(d_send_.data()), int(d_send_.size()),
                thrust::raw_pointer_cast(d_out_.data()));
        auto launch     = cudaGetLastError();
        auto completion = cudaStreamSynchronize(stream);
        if (!allValid(launch == cudaSuccess && completion == cudaSuccess)) return false;
        std::vector<MPI_Request> requests;
        requests.reserve(2 * registry_.peers().size());
        for (size_t p = 0; p < registry_.peers().size(); ++p)
            if (int n = recvOffsets_[p + 1] - recvOffsets_[p])
            {
                requests.emplace_back();
                MPI_Irecv(thrust::raw_pointer_cast(d_in_.data()) + recvOffsets_[p], n * int(sizeof(Face)), MPI_BYTE,
                          registry_.peers()[p], 1, comm_, &requests.back());
            }
        for (size_t p = 0; p < registry_.peers().size(); ++p)
            if (int n = sendOffsets_[p + 1] - sendOffsets_[p])
            {
                requests.emplace_back();
                MPI_Isend(thrust::raw_pointer_cast(d_out_.data()) + sendOffsets_[p], n * int(sizeof(Face)), MPI_BYTE,
                          registry_.peers()[p], 1, comm_, &requests.back());
            }
        if (!requests.empty()) MPI_Waitall(int(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
        if (!d_recv_.empty())
            detail::ghostScatterKernel<<<(d_recv_.size() + 255) / 256, 256, 0, stream>>>(
                thrust::raw_pointer_cast(d_in_.data()), thrust::raw_pointer_cast(d_recv_.data()), int(d_recv_.size()),
                faces);
        launch     = cudaGetLastError();
        completion = cudaStreamSynchronize(stream);
        return allValid(launch == cudaSuccess && completion == cudaSuccess);
    }
#endif
private:
    bool allValid(bool valid) const
    {
        int local = valid, global = 0;
        MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, comm_);
        return global != 0;
    }
    MPI_Comm comm_;
    GhostRegistry<uint64_t> registry_{"nonconformal-donors"};
    size_t count_ = 0;
    bool ready_   = false;
    std::vector<int> sendOffsets_, recvOffsets_;
#if defined(__CUDACC__) || defined(__HIPCC__)
    thrust::device_vector<int> d_send_, d_recv_;
    thrust::device_vector<Face> d_out_, d_in_;
#endif
};
} // namespace mars::fem
