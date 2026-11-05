#include "../mars_structured_runtime.hpp"
#include <cmath>
#include <iomanip>

namespace mars::test {
template<class DofView, class Field> struct FillBenchmarkField {
    DofView dofs;
    Field field;
    int iteration;
    MARS_INLINE_FUNCTION void operator()(Integer i) const {
        double point[3] = {};
        dofs.get_local_dof_coordinates(i, point);
        field(i) = dofs.is_owned(i) ? iteration + 2.0 + 3.0 * point[0] - point[1] + point[2] : -100.0;
    }
};
template<class DofView, class Field> struct CheckBenchmarkField {
    DofView dofs;
    Field field;
    int iteration;
    MARS_INLINE_FUNCTION void operator()(Integer i, int& errors) const {
        double point[3] = {};
        dofs.get_local_dof_coordinates(i, point);
        const double expected = iteration + 2.0 + 3.0 * point[0] - point[1] + point[2];
        errors += Kokkos::abs(field(i) - expected) > 1.e-12;
    }
};

template<class MeshType>
void measure(const char* name, int nx, int ny, int nz, bool report) {
    auto context = make_context(proc_allocation{}, MPI_COMM_WORLD);
    MeshType mesh(context);
    Kokkos::fence();
    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();
    generate_distributed_cube(mesh, nx, ny, nz);
    Kokkos::fence();
    double mesh_time = MPI_Wtime() - start;
    DofHandler<MeshType, 1> dofs(mesh);
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    dofs.enumerate_dofs();
    Kokkos::fence();
    double dof_time = MPI_Wtime() - start;
    DM<decltype(dofs), double> data(dofs);
    auto view = dofs.device_view();
    auto field = data.template get_dof_data<0>();
    const auto fill = [&](int iteration) {
        Kokkos::parallel_for("benchmark_fill", dofs.get_dof_size(),
                            FillBenchmarkField<decltype(view), decltype(field)>{view, field, iteration});
        Kokkos::fence();
    };
    long long cells = mesh.get_view_sfc().extent(0), owned = dofs.get_owned_dof_size();
    long long global_cells = 0, global_owned = 0;
    MPI_Allreduce(&cells, &global_cells, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&owned, &global_owned, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
    int errors = global_cells != 1LL * nx * ny * (nz ? nz : 1) ||
                 global_owned != 1LL * (nx + 1) * (ny + 1) * (nz ? nz + 1 : 1);
    fill(0);
    data.template gather_ghost_data<>();
    Kokkos::fence();
    double halo_time = 0;
    constexpr int iterations = 20;
    for (int iteration = 1; iteration <= iterations; ++iteration) {
        fill(iteration);
        MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();
        data.template gather_ghost_data<>();
        Kokkos::fence();
        halo_time += MPI_Wtime() - start;
        int mismatches = 0;
        Kokkos::parallel_reduce("benchmark_check", dofs.get_dof_size(),
            CheckBenchmarkField<decltype(view), decltype(field)>{view, field, iteration}, mismatches);
        errors += mismatches;
    }
    int global_errors = 0, rank = 0, ranks = 0;
    MPI_Allreduce(&errors, &global_errors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    if (global_errors) MPI_Abort(MPI_COMM_WORLD, 3);
    double local[] = {mesh_time, dof_time, halo_time / iterations}, maximum[3] = {};
    MPI_Reduce(local, maximum, 3, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0 && report) {
        std::cout << std::setprecision(17) << "MARS_STRUCTURED_PERF {\"case\":\"" << name
                  << "\",\"ranks\":" << ranks << ",\"nx\":" << nx << ",\"ny\":" << ny << ",\"nz\":" << nz
                  << ",\"cells\":" << global_cells << ",\"owned_dofs\":" << global_owned
                  << ",\"iterations\":" << iterations << ",\"errors\":" << global_errors
                  << ",\"mesh_seconds\":" << maximum[0] << ",\"dof_seconds\":" << maximum[1]
                  << ",\"halo_seconds\":" << maximum[2] << "}" << std::endl;
    }
}
} // namespace mars::test

int main(int argc, char** argv) {
    mars::test::StructuredRuntime runtime(argc, argv);
    using namespace mars;
    for (bool report : {false, true}) {
        test::measure<DistributedQuad4Mesh<MortonKey<Unsigned>>>("morton-quad-rect", 128, 96, 0, report);
        test::measure<DistributedQuad4Mesh<HilbertKey<Unsigned>>>("hilbert-quad-rect", 128, 96, 0, report);
        test::measure<DistributedQuad4Mesh<HilbertKey<Unsigned>>>("hilbert-quad-full", 128, 128, 0, report);
        test::measure<DistributedHex8Mesh<MortonKey<Unsigned>>>("morton-hex-rect", 32, 24, 16, report);
        test::measure<DistributedHex8Mesh<HilbertKey<Unsigned>>>("hilbert-hex-rect", 32, 24, 16, report);
        test::measure<DistributedHex8Mesh<HilbertKey<Unsigned>>>("hilbert-hex-full", 32, 32, 32, report);
    }
    return 0;
}
