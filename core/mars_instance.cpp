#include "mars_instance.hpp"

#include "mars_base.hpp"
#include "mars_config.hpp"

#ifdef WITH_MPI
#include <mpi.h>
#endif  // WITH_MPI

#ifdef WITH_KOKKOS
#include <Kokkos_Core.hpp>
#endif  // WITH_KOKKOS

namespace mars {
    void MARS::init(int argc, char* argv[]) {
#ifdef WITH_MPI
        MPI_Init(&argc, &argv);
#else  // WITH_MPI
#ifdef WITH_PAR_MOONOLITH
        MPI_Init(&argc, &argv);
#endif  // WITH_PAR_MOONOLITH
#endif  // WITH_MPI

#ifdef WITH_KOKKOS
        Kokkos::initialize(argc, argv);
#ifdef MARS_USE_CUDA
        cudaDeviceSetLimit(cudaLimitStackSize,
                           32768);  // set stack to 32KB only for cuda since it is not yet supported in kokkos.
#endif                              // MARS_USE_CUDA
#endif                              // WITH_KOKKOS
    }

    int MARS::finalize() {
#ifdef WITH_KOKKOS
        Kokkos::finalize();
#endif  // WITH_KOKKOS
#ifdef WITH_MPI
        return MPI_Finalize();
#else
#ifdef WITH_PAR_MOONOLITH
        return MPI_Finalize();
#else
        return 0;
#endif  // WITH_PAR_MOONOLITH
#endif  // WITH_MPI
    }

}  // namespace mars