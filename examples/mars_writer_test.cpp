#include "mars_pvtu_writer.hpp"

#include <algorithm>
#include <cassert>
#include <cstdlib>
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
#include "mars_mesh_reader.hpp"
#include "mars_mesh_writer.hpp"
#include "mars_oldest_edge.hpp"
#include "mars_ranked_edge.hpp"
#include "mars_test.hpp"
#include <err.h>

#include "mars_mesh_generation.hpp"

#ifdef WITH_KOKKOS
#include "mars_lepp_benchmark_kokkos.hpp"
#include "mars_test_kokkos.hpp"
#endif // WITH_KOKKOS

#ifdef WITH_PAR_MOONOLITH
#include "mars_moonolith_test.hpp"
#include <mpi.h>
#endif // WITH_PAR_MOONOLITH

#ifdef WITH_MPI
#include "mars_advection.hpp"
#include "mars_distributed_mesh_generation.hpp"
#include "mars_par_bisection.hpp"
#include "mars_par_mesh.hpp"
#include "mars_poisson.hpp"
#include "test.hpp"
#include "mars_test_mpi.hpp"
#endif // WITH_MPI
#include <chrono>

using namespace std::chrono;


int main(int argc, char *argv[]) {
  using namespace mars;

#ifdef WITH_PAR_MOONOLITH
  MPI_Init(&argc, &argv);
#endif


  int level = 1;
  int refine_level = 1;
  std::string filename = "../data/write/tetrakis.MFEM";
  if (argc > 1) {
    char *end_ptr = argv[1];
    level = strtol(argv[1], &end_ptr, 10);
    if (*end_ptr != '\0' || end_ptr == argv[1])
      warnx("'%s' could not be (completely) converted to long", argv[1]);

    if (argc == 3) {
      char *end_ptr = argv[2];
      refine_level = strtol(argv[2], &end_ptr, 10);
      if (*end_ptr != '\0' || end_ptr == argv[1])
        warnx("'%s' could not be (completely) converted to long", argv[1]);
    }

  } else
    std::cout
        << "No level of refinement was specified. Setting the default to 1!"
        << std::endl;

 
#ifdef WITH_KOKKOS
  Kokkos::initialize(argc, argv);
  {

#ifdef MARS_USE_CUDA
    cudaDeviceSetLimit(cudaLimitStackSize,
                       32768); // set stack to 32KB only for cuda since it is
                               // not yet supported in kokkos.
#endif
    mesh_test(argc, argv, level);

  }

  Kokkos::finalize();

#endif

#ifdef WITH_PAR_MOONOLITH
  run_mars_moonolith_test();
#endif // WITH_PAR_MOONOLITH


#ifdef WITH_PAR_MOONOLITH
  return MPI_Finalize();
#else
  return 0;
#endif
}
