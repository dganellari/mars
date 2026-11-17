#pragma once

#include <climits>
#include <limits>
#include <stdexcept>

#include "mars_amr_lepp_execution.hpp"
#include "backend/distributed/unstructured/mars_sfc_ownership.hpp"

#if defined(__CUDACC__) || defined(__HIPCC__)
#define MARS_AMR_PLAN_HD __host__ __device__
#else
#define MARS_AMR_PLAN_HD
#endif

namespace mars::amr
{
template<class KeyType, class Owner>
struct AmrNodeOwner
{
    Owner owner;
    int ranks;
    MARS_AMR_PLAN_HD int operator()(KeyType key) const { return ranks == 1 ? 0 : owner(key); }
};

namespace amr_detail
{
struct Route
{
    int rank;
    uint64_t index;
};
struct Peer
{
    int rank;
    uint64_t begin, count;
};

// Named operations avoid nvcc's constructor-lambda and private-type restrictions.
struct BuildRoutes
{
    Route* routes;
    uint64_t* flags;
    const int* destinations;
    size_t count;
    int ranks;
    MARS_AMR_PLAN_HD void operator()(size_t i) const
    {
        flags[i] = i < count && (destinations[i] < 0 || destinations[i] >= ranks);
        if (i < count) routes[i] = {destinations[i], i};
    }
};
struct RouteLess
{
    int rank, ranks;
    MARS_AMR_PLAN_HD bool operator()(const Route& a, const Route& b) const
    {
        int ar = a.rank == rank ? ranks : a.rank;
        int br = b.rank == rank ? ranks : b.rank;
        return ar < br || (ar == br && a.index < b.index);
    }
};
struct MarkPeerStarts
{
    const Route* routes;
    uint64_t* flags;
    size_t count;
    MARS_AMR_PLAN_HD void operator()(size_t i) const
    {
        flags[i] = i < count && (i == 0 || routes[i - 1].rank != routes[i].rank);
    }
};
struct BuildPeers
{
    const Route* routes;
    const uint64_t* flags;
    const uint64_t* offsets;
    uint64_t* order;
    Peer* peers;
    size_t count;
    MARS_AMR_PLAN_HD void operator()(size_t i) const
    {
        order[i] = routes[i].index;
        if (flags[i])
        {
            size_t end = i + 1;
            while (end < count && routes[end].rank == routes[i].rank)
                ++end;
            peers[offsets[i]] = {routes[i].rank, i, end - i};
        }
    }
};
} // namespace amr_detail

inline void require_collectively(bool valid, const char* message, MPI_Comm comm)
{
    int invalid = !valid;
    MPI_Allreduce(MPI_IN_PLACE, &invalid, 1, MPI_INT, MPI_MAX, comm);
    if (invalid) throw std::runtime_error(message);
}

class AmrCommunicator
{
public:
    explicit AmrCommunicator(MPI_Comm parent) { MPI_Comm_dup(parent, &comm_); }
    ~AmrCommunicator() { MPI_Comm_free(&comm_); }
    AmrCommunicator(const AmrCommunicator&)            = delete;
    AmrCommunicator& operator=(const AmrCommunicator&) = delete;
    operator MPI_Comm() const { return comm_; }

private:
    MPI_Comm comm_ = MPI_COMM_NULL;
};

// Discover incoming peers once. Repeated forward/reverse traffic uses the same
// permutation, counts and buffers; only peer metadata is copied to the host.
template<class Execution, class Packet>
class AmrSparsePlan
{
    static_assert(std::is_trivially_copyable_v<Packet>);
    template<class T>
    using Vector = typename Execution::template Vector<T>;
    using Route  = amr_detail::Route;
    using Peer   = amr_detail::Peer;

public:
    AmrSparsePlan(const Vector<int>& destinations, MPI_Comm comm, bool host_staged)
        : comm_(comm)
        , host_staged_(host_staged)
        , order_(destinations.size())
        , output_(destinations.size())
    {
        int rank, ranks;
        MPI_Comm_rank(comm_, &rank);
        MPI_Comm_size(comm_, &ranks);
        size_t n = destinations.size();
        Vector<Route> routes(n);
        Vector<uint64_t> flags(n + 1), offsets(n + 1);
        auto route  = Execution::data(routes);
        auto dest   = Execution::data(destinations);
        auto flag   = Execution::data(flags);
        auto offset = Execution::data(offsets);
        Execution::for_each(n + 1, amr_detail::BuildRoutes{route, flag, dest, n, ranks});
        require_collectively(Execution::sum(flag, n) == 0, "Invalid AMR peer rank", comm_);
        Execution::sort(route, n, amr_detail::RouteLess{rank, ranks});
        Execution::for_each(n + 1, amr_detail::MarkPeerStarts{route, flag, n});
        Execution::scan(flag, n + 1, offset);
        size_t peers = Execution::read(offset + n);
        Vector<Peer> peer_records(peers);
        auto peer  = Execution::data(peer_records);
        auto order = Execution::data(order_);
        Execution::for_each(n, amr_detail::BuildPeers{route, flag, offset, order, peer, n});
        send_peers_.resize(peers);
        Execution::download(peer, peers, send_peers_.data());
        // Keep local records on device, outside both MPI and staging buffers.
        remote_send_size_ = n;
        if (!send_peers_.empty() && send_peers_.back().rank == rank)
        {
            local_size_ = send_peers_.back().count;
            remote_send_size_ -= local_size_;
            send_peers_.pop_back();
        }
        peers = send_peers_.size();
        std::vector<int> ranks_out(peers), ones(peers, 1);
        std::vector<uint64_t> counts(peers), incoming;
        bool valid = true;
        for (size_t p = 0; p < peers; ++p)
        {
            ranks_out[p] = send_peers_[p].rank;
            counts[p]    = send_peers_[p].count;
            valid &= counts[p] <= size_t(INT_MAX) / sizeof(Packet);
        }
        require_collectively(valid, "AMR peer payload exceeds MPI count capacity", comm_);
        auto messages = mars::sparseExchange<uint64_t>(ranks_out, ones, counts.data(), incoming, 0, comm_);
        size_t total  = 0;
        for (const auto& message : messages)
        {
            if (message.count != 1)
            {
                valid = false;
                continue;
            }
            uint64_t count = incoming[message.offset];
            valid &= count <= size_t(INT_MAX) / sizeof(Packet) && count <= (SIZE_MAX / sizeof(Packet)) - total;
            recv_peers_.push_back({message.source, total, count});
            if (valid) total += count;
        }
        require_collectively(valid, "Invalid AMR incoming peer count", comm_);
        require_collectively(send_peers_.size() + recv_peers_.size() <= size_t(INT_MAX), "Too many AMR peer requests",
                             comm_);
        require_collectively(local_size_ <= SIZE_MAX / sizeof(Packet) &&
                                 total <= SIZE_MAX / sizeof(Packet) - local_size_,
                             "AMR local and incoming payload exceeds index capacity", comm_);
        remote_receive_size_ = total;
        send_.resize(n);
        received_.resize(total + local_size_);
        if (host_staged_)
        {
            h_send_.resize(remote_send_size_);
            h_received_.resize(total);
        }
        requests_.reserve(send_peers_.size() + recv_peers_.size());
    }

    const Vector<Packet>& forward(const Vector<Packet>& values)
    {
        require_collectively(values.size() == order_.size(), "AMR forward plan size changed", comm_);
        auto source = Execution::data(values);
        auto target = Execution::data(send_);
        auto order  = Execution::data(order_);
        Execution::for_each(order_.size(), [=] MARS_AMR_PLAN_HD(size_t i) { target[i] = source[order[i]]; });
        communicate(false, target, Execution::data(received_));
        return received_;
    }

    const Vector<Packet>& reverse(const Vector<Packet>& values)
    {
        require_collectively(values.size() == received_.size(), "AMR reverse plan size changed", comm_);
        communicate(true, Execution::data(values), Execution::data(send_));
        auto source = Execution::data(send_);
        auto target = Execution::data(output_);
        auto order  = Execution::data(order_);
        Execution::for_each(order_.size(), [=] MARS_AMR_PLAN_HD(size_t i) { target[order[i]] = source[i]; });
        return output_;
    }

    size_t incoming_size() const { return received_.size(); }
    size_t outgoing_peers() const { return send_peers_.size(); }
    size_t outgoing_records() const { return remote_send_size_; }
    size_t local_records() const { return local_size_; }

private:
    void communicate(bool reverse, const Packet* source, Packet* target)
    {
        const auto& sends       = reverse ? recv_peers_ : send_peers_;
        const auto& receives    = reverse ? send_peers_ : recv_peers_;
        auto& h_source          = reverse ? h_received_ : h_send_;
        auto& h_target          = reverse ? h_send_ : h_received_;
        Packet* original_target = target;
        if (local_size_)
            Execution::copy(source + (reverse ? remote_receive_size_ : remote_send_size_), local_size_,
                            target + (reverse ? remote_send_size_ : remote_receive_size_));
        Execution::synchronize();
        if (host_staged_)
        {
            Execution::download(source, h_source.size(), h_source.data());
            source = h_source.data();
            target = h_target.data();
        }
        requests_.clear();
        for (const Peer& peer : receives)
        {
            MPI_Request request;
            MPI_Irecv(target + peer.begin, int(peer.count * sizeof(Packet)), MPI_BYTE, peer.rank, 1, comm_, &request);
            requests_.push_back(request);
        }
        for (const Peer& peer : sends)
        {
            MPI_Request request;
            MPI_Isend(source + peer.begin, int(peer.count * sizeof(Packet)), MPI_BYTE, peer.rank, 1, comm_, &request);
            requests_.push_back(request);
        }
        MPI_Waitall(int(requests_.size()), requests_.data(), MPI_STATUSES_IGNORE);
        if (host_staged_) Execution::upload(target, h_target.size(), original_target);
    }

    MPI_Comm comm_;
    bool host_staged_;
    size_t remote_send_size_ = 0, remote_receive_size_ = 0, local_size_ = 0;
    Vector<uint64_t> order_;
    Vector<Packet> send_, received_, output_;
    std::vector<Peer> send_peers_, recv_peers_;
    std::vector<MPI_Request> requests_;
    std::vector<Packet> h_send_, h_received_;
};
} // namespace mars::amr

#undef MARS_AMR_PLAN_HD
