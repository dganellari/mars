#include "mars_instance.hpp"

#include "mars_base.hpp"
#include "mars_config.hpp"

#ifdef MARS_ENABLE_MPI
#include <mpi.h>
#endif  // MARS_ENABLE_MPI

#ifdef MARS_ENABLE_KOKKOS
#include <Kokkos_Core.hpp>
#endif  // MARS_ENABLE_KOKKOS

namespace mars {
    void MARS::init(int argc, char* argv[]) {
#ifdef MARS_ENABLE_MPI
        MPI_Init(&argc, &argv);
#else  // MARS_ENABLE_MPI
#ifdef MARS_ENABLE_PAR_MOONOLITH
        MPI_Init(&argc, &argv);
#endif  // MARS_ENABLE_PAR_MOONOLITH
#endif  // MARS_ENABLE_MPI

#ifdef MARS_ENABLE_KOKKOS
        Kokkos::initialize(argc, argv);
#ifdef MARS_ENABLE_CUDA
        cudaDeviceSetLimit(cudaLimitStackSize,
                           32768);  // set stack to 32KB only for cuda since it is not yet supported in kokkos.
#endif                              // MARS_ENABLE_CUDA
#endif                              // MARS_ENABLE_KOKKOS
    }

    int MARS::finalize() {
#ifdef MARS_ENABLE_KOKKOS
        Kokkos::finalize();
#endif  // MARS_ENABLE_KOKKOS
#ifdef MARS_ENABLE_MPI
        return MPI_Finalize();
#else
#ifdef MARS_ENABLE_PAR_MOONOLITH
        return MPI_Finalize();
#else
        return 0;
#endif  // MARS_ENABLE_PAR_MOONOLITH
#endif  // MARS_ENABLE_MPI
    }

}  // namespace mars