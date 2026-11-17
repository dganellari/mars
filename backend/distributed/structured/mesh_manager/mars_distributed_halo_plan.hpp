#ifndef MARS_DISTRIBUTED_HALO_PLAN_HPP
#define MARS_DISTRIBUTED_HALO_PLAN_HPP

#include "mars_sparse_exchange.hpp"
#include <array>
#include <tuple>

namespace mars {

namespace detail {
template<class DeviceView, bool Ghost>
struct HaloIndices {
    DeviceView view;
    ViewVectorType<Integer> indices;
    MARS_INLINE_FUNCTION void operator()(Integer i) const {
        if constexpr (Ghost) indices(i) = view.get_ghost_local_index(i);
        else indices(i) = view.get_boundary_local_index(i);
    }
};
} // namespace detail

// The layout is a snapshot. Rebuild the plan when the mesh, numbering or block changes.
// Construction is collective; start/wait only communicate with the listed peers.
template<class... T>
class HaloExchangePlan {
public:
    using UserTuple = ViewsTuple<T...>;
    template<class Handler>
    explicit HaloExchangePlan(const Handler& handler, bool isolate = true)
        : channel_(handler.get_context()->distributed->mpi_comm(), isolate),
          send_(mpi::sparse_segments(handler.get_view_scan_send_mirror().data(),
                                    num_ranks(handler.get_context()), handler.get_block(), channel_.get())),
          recv_(mpi::sparse_segments(handler.get_view_scan_recv_mirror().data(),
                                    num_ranks(handler.get_context()), handler.get_block(), channel_.get())),
          boundary_indices_("halo_boundary_indices", handler.get_boundary_dof_size()),
          ghost_indices_("halo_ghost_indices", handler.get_ghost_dof_size()),
          dof_size_(handler.get_dof_size()) {
        const auto view = handler.device_view();
        const auto boundary = boundary_indices_, ghosts = ghost_indices_;
        using Policy = Kokkos::RangePolicy<Kokkos::IndexType<Integer>>;
        Kokkos::parallel_for("halo_boundary_indices", Policy(0, boundary.extent(0)),
                            detail::HaloIndices<decltype(view), false>{view, boundary});
        Kokkos::parallel_for("halo_ghost_indices", Policy(0, ghosts.extent(0)),
                            detail::HaloIndices<decltype(view), true>{view, ghosts});
    }
    ~HaloExchangePlan() { wait(); Kokkos::fence("halo plan lifetime"); }
    HaloExchangePlan(const HaloExchangePlan&) = delete;
    HaloExchangePlan& operator=(const HaloExchangePlan&) = delete;

    template<Integer... fields>
    void start_gather(const UserTuple& data) { start<false, fields...>(data); }
    template<Integer... fields>
    void start_scatter(const UserTuple& data) { start<true, fields...>(data); }
    bool active() const { return active_; }
    std::size_t buffer_allocations() const { return allocations_; }
    std::size_t peer_count() const { return send_.size() + recv_.size(); }

    // Gather overwrites ghosts. Scatter returns contributions in boundary order.
    // Returned buffers are borrowed until the next exchange of the same field.
    UserTuple wait() {
        UserTuple result;
        if (!active_) return result;
        each([&](auto index) {
            constexpr auto i = decltype(index)::value;
            if (!selected_[i]) return;
            auto& exchange = *std::get<i>(exchanges_);
            exchange.wait();
            if (reverse_) std::get<i>(result) = exchange.send_buffer();
            else {
                const auto values = std::get<i>(data_);
                const auto buffer = exchange.recv_buffer();
                const auto indices = ghost_indices_;
                mpi::unpack_indexed("halo_unpack", values, buffer, indices);
            }
        });
        active_ = false;
        data_ = {};
        return result;
    }

private:
    template<class F, std::size_t... I>
    static void each_impl(F&& f, std::index_sequence<I...>) { (f(std::integral_constant<std::size_t, I>{}), ...); }
    template<class F>
    static void each(F&& f) { each_impl(std::forward<F>(f), std::index_sequence_for<T...>{}); }

    template<bool Reverse, Integer... fields>
    void start(const UserTuple& data) {
        static_assert(((fields >= 0 && fields < sizeof...(T)) && ...), "Invalid halo field index");
        if (active_) throw std::logic_error("Halo exchange already in flight");
        selected_.fill(false);
        if constexpr (sizeof...(fields) == 0) selected_.fill(true);
        else ((selected_[fields] = true), ...);
        // Allocate and validate all selected fields before posting any requests.
        each([&](auto index) {
            constexpr auto i = decltype(index)::value;
            if (!selected_[i]) return;
            if (std::get<i>(data).extent(0) != dof_size_)
                mpi::abort_invalid_count(channel_.get(), "halo field does not match plan layout");
            auto& exchange = std::get<i>(exchanges_);
            if (!exchange) {
                using Scalar = typename std::tuple_element<i, std::tuple<T...>>::type;
                exchange = std::make_unique<mpi::SparseExchange<Scalar>>(
                    channel_.get(), i, send_, recv_, boundary_indices_.extent(0), ghost_indices_.extent(0));
                ++allocations_;
            }
            const auto buffer = Reverse ? exchange->recv_buffer() : exchange->send_buffer();
            const auto indices = Reverse ? ghost_indices_ : boundary_indices_;
            const auto values = std::get<i>(data);
            mpi::pack_indexed("halo_pack", values, buffer, indices);
        });
        data_ = data;
        reverse_ = Reverse;
        active_ = true;
        each([&](auto index) {
            constexpr auto i = decltype(index)::value;
            if (selected_[i]) std::get<i>(exchanges_)->start(Reverse);
        });
    }

    mpi::ExchangeChannel channel_;
    std::vector<mpi::PeerSegment> send_, recv_;
    ViewVectorType<Integer> boundary_indices_, ghost_indices_;
    std::size_t dof_size_;
    std::tuple<std::unique_ptr<mpi::SparseExchange<T>>...> exchanges_;
    std::array<bool, sizeof...(T)> selected_{};
    UserTuple data_;
    std::size_t allocations_ = 0;
    bool active_ = false, reverse_ = false;
};
} // namespace mars
#endif
