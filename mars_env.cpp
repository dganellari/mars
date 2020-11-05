#include "mars_env.hpp"

#ifdef WITH_MPI
#include "mars_mpi_guard.hpp"
#endif  // WITH_MPI

#ifdef WITH_KOKKOS
#include <Kokkos_Core.hpp>
#endif  // WITH_KOKKOS

namespace mars {

    class Env::Impl {
    public:
        Impl(int argc, char *argv[])
            : error_code(0)
#ifdef WITH_MPI
              ,
              guard(argc, argv, false)
#endif  // WITH_MPI
        {
#ifdef WITH_KOKKOS
            Kokkos::initialize(argc, argv);
#endif  // WITH_KOKKOS

#ifdef MARS_USE_CUDA
            cudaDeviceSetLimit(cudaLimitStackSize,
                               32768);  // set stack to 32KB only for cuda since it is
                                        // not yet supported in kokkos.
#endif
        }

        ~Impl() {
#ifdef WITH_KOKKOS
            Kokkos::finalize();
#endif  // WITH_KOKKOS
        }

        int error_code;
#ifdef WITH_MPI
        marsenv::mpi_guard guard;
#endif  // WITH_MPI
    };

    Env::Env(int argc, char *argv[]) : impl_(std::make_unique<Impl>(argc, argv)) {}

    Env::~Env() {}

    Integer Env::exit_code() { return impl_->error_code; }

}  // namespace mars