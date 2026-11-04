#include <mars.hpp>
#include "../mars_structured_runtime.hpp"

#include <iostream>

#if !defined(MARS_ENABLE_KOKKOS) || !defined(MARS_ENABLE_MPI)
#error "Install MARS with the structured Kokkos backend and MPI to run this consumer."
#endif

int main(int argc, char** argv) {
    mars::test::StructuredRuntime env(argc, argv);
    int rank = 0;
    int ranks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    auto context = mars::make_context(mars::proc_allocation{}, MPI_COMM_WORLD);
    mars::DistributedQuad4Mesh<mars::MortonKey<mars::Unsigned>> mesh(context);
    mars::generate_distributed_cube(mesh, 4, 3, 0);

    using DofHandler = mars::DofHandler<decltype(mesh), 1>;
    DofHandler dofs(mesh);
    dofs.enumerate_dofs();
    mars::DM<DofHandler, double> data(dofs);
    const auto view = dofs.device_view();
    auto field = data.template get_dof_data<0>();
    Kokkos::parallel_for("initialize_owned_field", dofs.get_dof_size(), KOKKOS_LAMBDA(const mars::Integer i) {
        double point[3] = {};
        view.get_local_dof_coordinates(i, point);
        field(i) = view.is_owned(i) ? 2.0 + 3.0 * point[0] - point[1] : -100.0;
    });
    data.template gather_ghost_data<>();

    int errors = 0;
    Kokkos::parallel_reduce("check_gathered_field", dofs.get_dof_size(),
        KOKKOS_LAMBDA(const mars::Integer i, int& local_errors) {
            double point[3] = {};
            view.get_local_dof_coordinates(i, point);
            const double expected = 2.0 + 3.0 * point[0] - point[1];
            if (Kokkos::abs(field(i) - expected) > 1.e-12) ++local_errors;
        }, errors);

    long long owned = dofs.get_owned_dof_size();
    long long ghosts = dofs.get_dof_size() - owned;
    long long global_owned = 0;
    long long global_ghosts = 0;
    MPI_Allreduce(&owned, &global_owned, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&ghosts, &global_ghosts, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
    if (dofs.get_global_dof_size() != 20 || global_owned != 20 || (ranks > 1 && global_ghosts == 0)) ++errors;
    int global_errors = 0;
    MPI_Allreduce(&errors, &global_errors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << "Structured installed consumer: " << Kokkos::DefaultExecutionSpace::name()
                  << ", " << ranks << " MPI ranks, " << global_errors << " errors\n";
    }
    return global_errors == 0 ? 0 : 1;
}
