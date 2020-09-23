#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <err.h>
#include <iostream>
#include <numeric>

#include "mars_benchmark.hpp"
#include "mars_bisection.hpp"
#include "mars_lagrange_element.hpp"
#include "mars_lepp_benchmark.hpp"
#include "mars_mesh.hpp"
#include "mars_mesh_partition.hpp"
#include "mars_partitioned_bisection.hpp"
#include "mars_prelepp_benchmark.hpp"
#include "mars_quality.hpp"
#include "mars_simplex.hpp"
#include "mars_utils.hpp"
#include "mars_vtk_writer.hpp"

#include "mars_longest_edge.hpp"
#include "mars_memory.hpp"
#include "mars_mesh_reader.hpp"
#include "mars_mesh_writer.hpp"
#include "mars_oldest_edge.hpp"
#include "mars_ranked_edge.hpp"
#include "mars_test.hpp"

#include "mars_mesh_generation.hpp"

#ifdef WITH_KOKKOS
#include "mars_lepp_benchmark_kokkos.hpp"
#include "mars_test_kokkos.hpp"
#endif // WITH_KOKKOS
int main(int argc, char *argv[]) {
  using namespace mars;
  Kokkos::initialize(argc, argv);

#ifdef MARS_USE_CUDA
  cudaDeviceSetLimit(cudaLimitStackSize,
                     32768); // set stack to 32KB only for cuda since it is
                             // not yet supported in kokkos.
#endif

  {
    using PMesh = ParallelMesh2;
    using SMesh = Mesh2;

    Integer nx = 10, ny = 10, nz = 0;
    PMesh mesh;
    generate_cube(mesh, nx, ny, nz);

    ParallelBisection<PMesh> bisection(&mesh);

    ViewVectorType<Integer> marked("marked", 1);

    bisection.refine(marked);

    SMesh serial_mesh;
    convert_parallel_mesh_to_serial(serial_mesh, mesh);

    std::cout << "n_active_elements: " << serial_mesh.n_active_elements()
              << std::endl;
    std::cout << "n_nodes: " << serial_mesh.n_nodes() << std::endl;

    VTKMeshWriter<SMesh> w;
    w.write("mesh.vtu", serial_mesh);
  }

  Kokkos::finalize();
}