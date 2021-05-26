#ifndef MARS_ENV_HPP
#define MARS_ENV_HPP

#include <memory>
#include "mars_base.hpp"

#ifdef WITH_MPI
#include <mpi.h>
#endif  // WITH_MPI

namespace mars {
    class Env {
    public:
        class Impl;

        Env(int argc, char *argv[]);
#ifdef WITH_MPI
        Env(int argc, char *argv[], MPI_Comm comm);
#endif
        ~Env();
        Integer exit_code();

    private:
        std::unique_ptr<Impl> impl_;
    };
}  // namespace mars

#endif  // MARS_ENV_HPP
