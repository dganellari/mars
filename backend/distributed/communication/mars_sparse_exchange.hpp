#ifndef MARS_SPARSE_EXCHANGE_HPP
#define MARS_SPARSE_EXCHANGE_HPP

#include "mars_mpi.hpp"
#include <memory>

namespace mars::mpi {

template<class Values, class Buffer, class Indices>
void pack_indexed(const char* label, Values values, Buffer buffer, Indices indices) {
    using Policy = Kokkos::RangePolicy<Kokkos::IndexType<Integer>>;
    Kokkos::parallel_for(label, Policy(0, indices.extent(0)), KOKKOS_LAMBDA(Integer i) {
        buffer(i) = values(indices(i));
    });
}

template<class Values, class Buffer, class Indices>
void unpack_indexed(const char* label, Values values, Buffer buffer, Indices indices) {
    using Policy = Kokkos::RangePolicy<Kokkos::IndexType<Integer>>;
    Kokkos::parallel_for(label, Policy(0, indices.extent(0)), KOKKOS_LAMBDA(Integer i) {
        values(indices(i)) = buffer(i);
    });
}

// Construct collectively once; exchanges on the resulting channel are peer-only.
class ExchangeChannel {
public:
    explicit ExchangeChannel(MPI_Comm parent, bool isolate = true) : owned_(isolate) {
        if (parent == MPI_COMM_NULL) throw std::invalid_argument("Sparse exchange requires an MPI context");
        if (isolate) {
            MPI_OR_THROW(MPI_Comm_dup, parent, &comm_);
            MPI_OR_THROW(MPI_Comm_set_errhandler, comm_, MPI_ERRORS_RETURN);
        } else comm_ = parent;
    }
    ~ExchangeChannel() { if (owned_ && comm_ != MPI_COMM_NULL) MPI_Comm_free(&comm_); }
    ExchangeChannel(const ExchangeChannel&) = delete;
    ExchangeChannel& operator=(const ExchangeChannel&) = delete;
    MPI_Comm get() const { return comm_; }
private:
    MPI_Comm comm_ = MPI_COMM_NULL;
    bool owned_;
};

struct PeerSegment {
    int rank;
    int count;
    std::size_t offset;
};

inline std::vector<PeerSegment> sparse_segments(const Integer* scan, int ranks, Integer block,
                                               MPI_Comm comm) {
    if (block <= 0 || scan[0] != 0) abort_invalid_count(comm, "invalid exchange layout");
    std::vector<PeerSegment> peers;
    for (int p = 0; p < ranks; ++p) {
        if (scan[p] < 0 || scan[p + 1] < scan[p] ||
            static_cast<std::uint64_t>(scan[p + 1]) > std::numeric_limits<std::size_t>::max() / block)
            abort_invalid_count(comm, "invalid exchange displacement");
        const auto count = static_cast<std::uint64_t>(scan[p + 1] - scan[p]) * block;
        if (count) peers.push_back({p, checked_mpi_count(count, comm),
                                   static_cast<std::size_t>(scan[p]) * static_cast<std::size_t>(block)});
    }
    return peers;
}

// Owns both MPI buffers until wait(), including the host staging allocations.
// A plan may have one exchange in flight and must die before MPI/Kokkos finalize.
template<class T>
class SparseExchange {
public:
    using View = ViewVectorType<T>;
    SparseExchange(MPI_Comm comm, int tag, std::vector<PeerSegment> send,
                   std::vector<PeerSegment> recv, std::size_t send_size, std::size_t recv_size)
        : comm_(comm), tag_(tag), send_(std::move(send)), recv_(std::move(recv)),
          send_buffer_("halo_send", send_size), recv_buffer_("halo_recv", recv_size),
          requests_(send_.size() + recv_.size(), MPI_REQUEST_NULL) {
        static_assert(mpi_traits<T>::is_mpi_native_type(), "Sparse exchange needs a native MPI scalar");
        validate(send_, send_size);
        validate(recv_, recv_size);
        // The portable tag range also works on graph communicators without inherited attributes.
        if (tag < 0 || tag > 32767) abort_invalid_count(comm_, "invalid exchange tag");
#ifdef MARS_NO_RDMA
        h_send_ = typename View::HostMirror("halo_host_send", send_size);
        h_recv_ = typename View::HostMirror("halo_host_recv", recv_size);
#endif
    }
    ~SparseExchange() { wait(); }
    SparseExchange(const SparseExchange&) = delete;
    SparseExchange& operator=(const SparseExchange&) = delete;

    const View& send_buffer() const { return send_buffer_; }
    const View& recv_buffer() const { return recv_buffer_; }
    std::size_t peer_count() const { return send_.size() + recv_.size(); }
    bool active() const { return active_; }

    void start(bool reverse = false) {
        if (active_) throw std::logic_error("Exchange already in flight");
        reverse_ = reverse;
        Kokkos::fence("halo producers before MPI");
        T* source = reverse ? recv_buffer_.data() : send_buffer_.data();
        T* target = reverse ? send_buffer_.data() : recv_buffer_.data();
#ifdef MARS_NO_RDMA
        if (reverse) Kokkos::deep_copy(h_recv_, recv_buffer_);
        else Kokkos::deep_copy(h_send_, send_buffer_);
        source = reverse ? h_recv_.data() : h_send_.data();
        target = reverse ? h_send_.data() : h_recv_.data();
#endif
        const auto& incoming = reverse ? send_ : recv_;
        const auto& outgoing = reverse ? recv_ : send_;
        active_ = true;
        std::size_t i = 0;
        for (const auto& peer : incoming)
            check(MPI_Irecv(target + peer.offset, peer.count, mpi_traits<T>::mpi_type(), peer.rank,
                            tag_, comm_, &requests_[i++]));
        for (const auto& peer : outgoing)
            check(MPI_Isend(source + peer.offset, peer.count, mpi_traits<T>::mpi_type(), peer.rank,
                            tag_, comm_, &requests_[i++]));
    }

    void wait() noexcept {
        if (!active_) return;
        check(MPI_Waitall(static_cast<int>(requests_.size()), requests_.data(), MPI_STATUSES_IGNORE));
#ifdef MARS_NO_RDMA
        if (reverse_) Kokkos::deep_copy(send_buffer_, h_send_);
        else Kokkos::deep_copy(recv_buffer_, h_recv_);
#endif
        active_ = false;
    }

private:
    void validate(const std::vector<PeerSegment>& peers, std::size_t extent) {
        int ranks = size(comm_);
        int previous = -1;
        std::size_t end = 0;
        for (const auto& peer : peers) {
            if (peer.rank <= previous || peer.rank >= ranks || peer.count <= 0 || peer.offset < end ||
                peer.offset > extent || static_cast<std::size_t>(peer.count) > extent - peer.offset)
                abort_invalid_count(comm_, "invalid sparse peer segment");
            previous = peer.rank;
            end = peer.offset + peer.count;
        }
        checked_mpi_count(requests_.size(), comm_);
    }
    void check(int error) const noexcept {
        // A rank-local exception after posting requests could strand the other peers.
        if (error != MPI_SUCCESS) { MPI_Abort(comm_, error); std::abort(); }
    }
    MPI_Comm comm_;
    int tag_;
    std::vector<PeerSegment> send_, recv_;
    View send_buffer_, recv_buffer_;
#ifdef MARS_NO_RDMA
    typename View::HostMirror h_send_, h_recv_;
#endif
    std::vector<MPI_Request> requests_;
    bool active_ = false;
    bool reverse_ = false;
};

} // namespace mars::mpi
#endif
