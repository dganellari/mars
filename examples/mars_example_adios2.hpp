#include "mars_context.hpp"
#include "mars_globals.hpp"
// #include <bits/c++config.h>
#include <adios2.h>
#include <exception>
#include <iostream>
#include <tuple>
#include <type_traits>
#include <utility>

#ifdef WITH_MPI

#include "mars_mpi_guard.hpp"

#ifdef WITH_KOKKOS
#include <KokkosBlas1_sum.hpp>
#include "mars_boundary_conditions.hpp"
#include "mars_distributed_data_management.hpp"
#include "mars_distributed_mesh_generation.hpp"
#include "mars_dm_interpolate.hpp"
#include "mars_laplace_ex.hpp"
#include "mars_poisson_operator.hpp"
#include "mars_precon_conjugate_grad.hpp"
#include "mars_quad4.hpp"
#include "mars_serial_mesh_type.hpp"
#endif  // WITH_KOKKOS
#endif

namespace mars {

    // using VectorReal = mars::ViewVectorType<Real>;

    // // Definitions...
    // // Write method in mars_model_test.hpp
    // bool write(VectorReal &x) {
    //     std::cout << "Writing results to disk..." << std::endl;

    //     Integer n_nodes = mesh.n_nodes();

    //     VectorReal::HostMirror x_host("x_host", n_nodes);
    //     VectorReal::HostMirror rhs_host("rhs_host", n_nodes);
    //     Kokkos::deep_copy(x_host, x);
    //     Kokkos::deep_copy(rhs_host, rhs);

    //     SMesh serial_mesh;
    //     convert_parallel_mesh_to_serial(serial_mesh, mesh);

    //     VTUMeshWriter<SMesh> w;

    //     if (!w.write("solution.vtu", serial_mesh, x_host)) {
    //         return false;
    //     }

    //     if (!w.write("rhs.vtu", serial_mesh, rhs_host)) {
    //         return false;
    //     }

    //     return true;
    // }
}  // namespace mars