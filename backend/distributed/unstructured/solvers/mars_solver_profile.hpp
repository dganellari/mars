#pragma once

#include <array>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cuda_runtime.h>
#include <mpi.h>

namespace mars::fem {

class SolverProfile
{
public:
    enum Phase { Context, Assembly, Prepare, Setup, Solve, Finish, Preconditioner, Correction, Count };

    void initialize(MPI_Comm comm)
    {
        if (initialized_) return;
        comm_ = comm;
        int rank = 0, flag = 0;
        MPI_Comm_rank(comm_, &rank);
        if (rank == 0)
        {
            const char* value = std::getenv("MARS_OUTLET_PROFILE");
            flag = value && value[0] == '1' && value[1] == '\0';
        }
        // A rank-local environment flag must not select collective participation.
        MPI_Bcast(&flag, 1, MPI_INT, 0, comm_);
        enabled_ = flag != 0;
        initialized_ = true;
    }

    bool enabled() const { return enabled_; }

    void begin_step()
    {
        if (!enabled_) return;
        ++step_;
        milliseconds_.fill(0);
        calls_.fill(0);
    }

    double stamp() const
    {
        if (!enabled_) return 0;
        const auto error = cudaDeviceSynchronize();
        if (error != cudaSuccess)
        {
            std::cerr << "ERROR: outlet profiling: " << cudaGetErrorString(error) << '\n';
            MPI_Abort(comm_, 1);
            std::abort();
        }
        return MPI_Wtime();
    }

    double lap(Phase phase, double start)
    {
        if (!enabled_) return 0;
        const double end = stamp();
        milliseconds_[phase] += 1000 * (end - start);
        calls_[phase] += 1;
        return end;
    }

    class Scope
    {
    public:
        Scope(SolverProfile& profile, Phase phase, bool report = false)
            : profile_(profile), phase_(phase), start_(profile.stamp()), report_(report) {}
        ~Scope()
        {
            profile_.lap(phase_, start_);
            if (report_) profile_.report();
        }
        Scope(const Scope&) = delete;
        Scope& operator=(const Scope&) = delete;
    private:
        SolverProfile& profile_;
        Phase phase_;
        double start_;
        bool report_;
    };

    void report() const
    {
        if (!enabled_) return;
        std::array<double, Count> maximum{}, minimum_calls{}, maximum_calls{};
        MPI_Reduce(milliseconds_.data(), maximum.data(), Count, MPI_DOUBLE, MPI_MAX, 0, comm_);
        MPI_Reduce(calls_.data(), minimum_calls.data(), Count, MPI_DOUBLE, MPI_MIN, 0, comm_);
        MPI_Reduce(calls_.data(), maximum_calls.data(), Count, MPI_DOUBLE, MPI_MAX, 0, comm_);
        int rank = 0;
        MPI_Comm_rank(comm_, &rank);
        if (rank != 0) return;
        constexpr const char* names[Count] = {
            "context", "assembly", "hypre_prepare", "hypre_setup", "hypre_solve",
            "hypre_finish", "preconditioner", "correction"};
        std::ostringstream line;
        line << "[outlet-profile] step=" << step_ << " aggregation=rank_max";
        for (int phase = 0; phase < Count; ++phase)
            line << ' ' << names[phase] << "_ms=" << std::fixed << std::setprecision(3) << maximum[phase]
                 << ' ' << names[phase] << "_calls_min=" << std::setprecision(0) << minimum_calls[phase]
                 << ' ' << names[phase] << "_calls_max=" << maximum_calls[phase];
        std::cout << line.str() << '\n';
    }

private:
    bool initialized_ = false;
    bool enabled_ = false;
    int step_ = 0;
    MPI_Comm comm_ = MPI_COMM_WORLD;
    std::array<double, Count> milliseconds_{};
    std::array<double, Count> calls_{};
};

} // namespace mars::fem
