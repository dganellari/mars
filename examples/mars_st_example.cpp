#include "Kokkos_ArithTraits.hpp"
#include "Kokkos_Bitset.hpp"
#include "Kokkos_Parallel.hpp"
#include "Kokkos_Parallel_Reduce.hpp"
#include "mars_context.hpp"
#include "mars_globals.hpp"
#include <exception>
#include <iostream>
#include <tuple>
#include <type_traits>
#include <utility>

#include "mars_mpi_guard.hpp"

int main(int argc, char *argv[]) {
  using namespace mars;
  Kokkos::initialize(argc, argv);

  {

#ifdef MARS_USE_CUDA
    cudaDeviceSetLimit(cudaLimitStackSize,
                       32768); // set stack to 32KB only for cuda since it is
                               // not yet supported in kokkos.
#endif
  }

  Kokkos::finalize();
}