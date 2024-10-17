#ifndef MARS_ENV_HPP
#define MARS_ENV_HPP

#include <memory>
#include "mars_base.hpp"

#ifdef MARS_ENABLE_MPI
#include "mars_mpi_guard.hpp"
#include "mars_context.hpp"
#include "mars_execution_context.hpp"
#endif  // MARS_ENABLE_MPI

namespace mars {
    class Env {
    public:
        class Impl;

        Env(int argc, char *argv[]);
#ifdef MARS_ENABLE_MPI
        Env(int argc, char *argv[], MPI_Comm comm);
#endif
        ~Env();
        Integer exit_code();

    private:
        std::unique_ptr<Impl> impl_;
    };
}  // namespace mars

#endif  // MARS_ENV_HPP
