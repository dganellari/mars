#pragma once

#include <mars.hpp>
#include <type_traits>

namespace mars::test {

template<class T>
concept HasContext = requires(const T& value) { value.get_context(); };

inline void check_device_snapshots_and_plans() {
    using Mesh = DistributedQuad4Mesh<MortonKey<Unsigned>>;
    using Handler = DofHandler<Mesh, 2, 2>;
    using Data = DM<Handler, double>;
    static_assert(!std::is_polymorphic_v<Mesh::DeviceView>);
    static_assert(!std::is_polymorphic_v<Handler::DeviceView>);
    static_assert(!std::is_polymorphic_v<Data::DeviceView>);
    static_assert(!std::is_polymorphic_v<decltype(std::declval<Mesh::DeviceView>().elem(0))>);
    static_assert(!HasContext<Mesh::DeviceView> && !HasContext<Handler::DeviceView> && !HasContext<Data::DeviceView>);
    static_assert(sizeof(Mesh::DeviceView) < sizeof(Mesh));
    static_assert(sizeof(Handler::DeviceView) < sizeof(Handler));
    static_assert(std::is_copy_constructible_v<Data>);

    const auto retained = [] {
        auto context = make_context(proc_allocation{}, MPI_COMM_WORLD);
        Mesh mesh(context);
        generate_distributed_cube(mesh, 5, 3, 0);
        Handler handler(mesh);
        handler.enumerate_dofs();
        Data first(handler), second(handler);
        const auto view = handler.device_view();
        const auto a = first.get_dof_data<0>(), b = second.get_dof_data<0>();
        Kokkos::parallel_for("initialize_distinct_channels", view.get_dof_size(), KOKKOS_LAMBDA(Integer i) {
            a(i) = view.is_owned(i) ? view.local_to_global(i) + 0.25 : -1;
            b(i) = view.is_owned(i) ? view.local_to_global(i) + 0.75 : -1;
        });
        // Different posting order must not cross-match the two plans' tags.
        if (rank(context) % 2) { first.start_gather<>(); second.start_gather<>(); }
        else { second.start_gather<>(); first.start_gather<>(); }
        second.wait();
        first.wait();
        int errors = 0;
        Kokkos::parallel_reduce("isolated_halo_channels", view.get_dof_size(),
            KOKKOS_LAMBDA(Integer i, int& n) {
                n += a(i) != view.local_to_global(i) + 0.25;
                n += b(i) != view.local_to_global(i) + 0.75;
            }, errors);
        EXPECT_EQ(errors, 0);
        const auto snapshot = first.device_view();
        // Replacing host mesh storage cannot change an already constructed DOF layout.
        generate_distributed_cube(mesh, 7, 2, 0);
        EXPECT_EQ(snapshot.get_dof_handler().get_XMax(), 10);
        EXPECT_EQ(snapshot.get_dof_handler().get_YMax(), 6);
        // Destruction must complete an outstanding operation before freeing its buffers/channel.
        first.start_gather<>();
        return snapshot;
    }();
    int errors = 0;
    const auto view = retained.get_dof_handler();
    Kokkos::parallel_reduce("retained_field_snapshot", view.get_dof_size(),
        KOKKOS_LAMBDA(Integer i, int& n) {
            n += retained.get_dof_data<0>(i) != view.local_to_global(i) + 0.25;
            double point[3] = {};
            view.get_local_dof_coordinates(i, point);
            n += point[0] < 0 || point[0] > 1 || point[1] < 0 || point[1] > 1;
        }, errors);
    EXPECT_EQ(errors, 0);
}

inline void check_manager_copies() {
    using Mesh = DistributedQuad4Mesh<MortonKey<Unsigned>>;
    using Handler = DofHandler<Mesh, 2, 2>;
    using Data = DM<Handler, double>;
    using Cells = UserData<Mesh, double, int>;
    static_assert(std::is_copy_constructible_v<Data> && std::is_copy_constructible_v<Cells>);
    auto context = make_context(proc_allocation{}, MPI_COMM_WORLD);
    Mesh mesh(context);
    generate_distributed_cube(mesh, 5, 3, 0);
    Handler handler(mesh);
    handler.enumerate_dofs();
    auto original = std::make_unique<Data>(handler);
    auto alias = *original;
    EXPECT_EQ(&alias.exchange_plan(), &original->exchange_plan());
    if (rank(context) == 0) {
        const auto rank_local_copy = alias;
        EXPECT_EQ(rank_local_copy.get_dof_data<0>().data(), alias.get_dof_data<0>().data());
    }
    Kokkos::parallel_for("legacy_manager_capture", handler.get_dof_size(), KOKKOS_LAMBDA(Integer i) {
        const auto& dofs = alias.get_dof_handler();
        alias.get_dof_data<0>(i) = dofs.is_owned(i) ? dofs.local_to_global(i) + 0.25 : -1;
    });
    original->start_gather<>();
    EXPECT_THROW(alias.start_gather<>(), std::logic_error);
    original.reset();
    EXPECT_TRUE(alias.exchange_plan().active());
    alias.wait();
    const auto allocations = alias.exchange_plan().buffer_allocations();
    alias.gather_ghost_data<>();
    EXPECT_EQ(alias.exchange_plan().buffer_allocations(), allocations);
    int errors = 0;
    Kokkos::parallel_reduce("copied_manager_values", handler.get_dof_size(),
        KOKKOS_LAMBDA(Integer i, int& n) {
            n += alias.get_data<0>(i) != alias.get_dof_handler().local_to_global(i) + 0.25;
        }, errors);
    EXPECT_EQ(errors, 0);

    auto cells = std::make_unique<Cells>(mesh);
    auto before_exchange = *cells;
    const auto keys = mesh.get_view_sfc(), ghost_keys = mesh.get_view_ghost();
    Kokkos::parallel_for("legacy_cell_capture", mesh.get_chunk_size(), KOKKOS_LAMBDA(Integer i) {
        before_exchange.set_elem_data<0>(i, keys(i) + 0.5);
        before_exchange.set_elem_data<1>(i, keys(i) + 7);
    });
    cells->exchange_ghost_data(context);
    auto after_exchange = *cells;
    before_exchange.exchange_ghost_data(context);
    EXPECT_EQ(before_exchange.get_ghost_data<0>().data(), after_exchange.get_ghost_data<0>().data());
    cells.reset();
    after_exchange.exchange_ghost_data(context);
    errors = 0;
    Kokkos::parallel_reduce("copied_cell_values", mesh.get_ghost_size(),
        KOKKOS_LAMBDA(Integer i, int& n) {
            n += before_exchange.get_ghost_elem_data<0>(i) != ghost_keys(i) + 0.5;
            n += after_exchange.get_ghost_elem_data<1>(i) != ghost_keys(i) + 7;
        }, errors);
    EXPECT_EQ(errors, 0);
}

inline void check_sparse_subset_plan() {
    int rank = mpi::rank(MPI_COMM_WORLD), ranks = mpi::size(MPI_COMM_WORLD);
    if (ranks >= 2 && rank < 2) {
        const int peer = 1 - rank;
        mpi::SparseExchange<int> exchange(MPI_COMM_WORLD, 2718, {{peer, 1, 0}}, {{peer, 1, 0}}, 1, 1);
        const auto send = exchange.send_buffer(), recv = exchange.recv_buffer();
        for (int step = 0; step < 3; ++step) {
            Kokkos::deep_copy(send, rank + 10 * step);
            exchange.start();
            EXPECT_TRUE(exchange.active());
            EXPECT_THROW(exchange.start(), std::logic_error);
            exchange.wait();
            const auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, recv);
            EXPECT_EQ(host(0), peer + 10 * step);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

} // namespace mars::test
