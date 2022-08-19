#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include "mars_config.hpp"
#include "mars_err.hpp"

#include "mars.hpp"
#include "mars_env.hpp"

#ifdef WITH_PAR_MOONOLITH
#include <mpi.h>

#include "mars_moonolith_test.hpp"
#endif  // WITH_PAR_MOONOLITH

#ifdef WITH_MPI
#ifdef WITH_KOKKOS_KERNELS
#include "mars_advection.hpp"
#include "mars_constant_viscosity_stokes.hpp"
#include "mars_poisson.hpp"
#include "mars_test_kokkos.hpp"
#include "mars_test_mpi.hpp"
#include "mars_variable_viscosity_stokes.hpp"
#endif  // WITH_KokkosKernels
#endif  // WITH_MPI
#include <chrono>
using namespace std::chrono;

// void test_mars_mesh_generation(const int x) {
//     using namespace mars;

//     high_resolution_clock::time_point t1 = high_resolution_clock::now();

//     Mesh1 mesh;
//     generate_line(mesh, x);

//     high_resolution_clock::time_point t2 = high_resolution_clock::now();
//     auto duration = duration_cast<seconds>(t2 - t1).count();

//     std::cout << "Generation took: " << duration << " seconds." << std::endl;

//     std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;
//     std::cout << "n_nodes: " << mesh.n_nodes() << std::endl;

//     if (x < 100) {
//         VTKMeshWriter<Mesh1> w;
//         w.write("build_line" + std::to_string(x) + ".vtu", mesh);
//     }

//     return mesh;
// }

void test_mars_mesh_generation(const int x) {
    Kokkos::Timer timer;

    mars::ParallelMesh3 pMesh;
    mars::generate_cube(pMesh, x, x, x);

    // Kokkos::fence();

    double time = timer.seconds();

    std::cout << "Generation 3D kokkos took: " << time << " seconds." << std::endl;
}

int main(int argc, char *argv[]) {
    mars::Env env(argc, argv);
    test_mars_mesh_generation(100);
    return env.exit_code();
}
