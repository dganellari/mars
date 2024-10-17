#include "mars_env.hpp"

#ifdef MARS_ENABLE_KOKKOS
#include <Kokkos_Core.hpp>
#endif  // MARS_ENABLE_KOKKOS

#include <memory>

namespace mars {

    class Env::Impl {
    public:
#ifdef MARS_ENABLE_MPI
        static void silence_warning(const MPI_Comm &) {}

        Impl(int argc, char *argv[], MPI_Comm comm) : error_code(0) {
            silence_warning(comm);
#ifdef MARS_ENABLE_KOKKOS
            Kokkos::initialize(argc, argv);
#endif  // MARS_ENABLE_KOKKOS
        }
#endif  // MARS_ENABLE_MPI

        Impl(int argc, char *argv[])
            : error_code(0)
#ifdef MARS_ENABLE_MPI
              ,
              guard(std::make_unique<marsenv::mpi_guard>(argc, argv, false))
#endif  // MARS_ENABLE_MPI
        {
#ifdef MARS_ENABLE_KOKKOS
            Kokkos::initialize(argc, argv);
#endif  // MARS_ENABLE_KOKKOS

#ifdef MARS_ENABLE_CUDA
            cudaDeviceSetLimit(cudaLimitStackSize,
                               32768);  // set stack to 32KB only for cuda since it is
                                        // not yet supported in kokkos.
#endif
        }

        ~Impl() {
#ifdef MARS_ENABLE_KOKKOS
            Kokkos::finalize();
#endif  // MARS_ENABLE_KOKKOS
        }

        int error_code;
#ifdef MARS_ENABLE_MPI
        std::unique_ptr<marsenv::mpi_guard> guard;
#endif  // MARS_ENABLE_MPI
    };

    Env::Env(int argc, char *argv[]) : impl_(std::make_unique<Impl>(argc, argv)) {}
#ifdef MARS_ENABLE_MPI
    Env::Env(int argc, char *argv[], MPI_Comm comm) : impl_(std::make_unique<Impl>(argc, argv, comm)) {}
#endif

    Env::~Env() {}

    Integer Env::exit_code() { return impl_->error_code; }

}  // namespace mars
