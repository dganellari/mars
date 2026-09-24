// GPT/Codex, 2026-09-12: real MPI with a test-only CUDA synchronization stub.
#include "backend/distributed/unstructured/solvers/mars_solver_profile.hpp"
#include <string>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, ranks = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    auto require = [&](bool ok) {
        if (!ok) MPI_Abort(MPI_COMM_WORLD, 2);
    };
    using Profile = mars::fem::SolverProfile;
    // Deliberately disagree across ranks: only root's setting may govern reductions.
    setenv("MARS_OUTLET_PROFILE", rank == 0 ? "0" : "1", 1);
    Profile disabled;
    disabled.initialize(MPI_COMM_WORLD);
    disabled.initialize(MPI_COMM_WORLD);
    require(!disabled.enabled());
    std::ostringstream output;
    auto* original = std::cout.rdbuf(output.rdbuf());
    disabled.begin_step();
    { Profile::Scope scope(disabled, Profile::Preconditioner, true); }
    disabled.lap(Profile::Assembly, disabled.stamp());
    disabled.report();
    require(output.str().empty() && outlet_profile_test_syncs == 0);

    setenv("MARS_OUTLET_PROFILE", rank == 0 ? "1" : "0", 1);
    Profile enabled;
    enabled.initialize(MPI_COMM_WORLD);
    require(enabled.enabled());
    enabled.begin_step();
    {
        Profile::Scope step(enabled, Profile::Correction, true);
        for (int call = 0; call <= rank; ++call)
        {
            Profile::Scope scope(enabled, Profile::Preconditioner);
            const double start = enabled.stamp();
            enabled.lap(Profile::Solve, start);
        }
    }
    require(outlet_profile_test_syncs > 0);
    if (rank == 0)
    {
        require(output.str().find("step=1 aggregation=rank_max") != std::string::npos);
        require(output.str().find("correction_calls_min=1 correction_calls_max=1") != std::string::npos);
        require(output.str().find("hypre_solve_calls_min=1 hypre_solve_calls_max=" + std::to_string(ranks)) != std::string::npos);
        require(output.str().find("preconditioner_calls_min=1 preconditioner_calls_max=" + std::to_string(ranks)) != std::string::npos);
    }
    else require(output.str().empty());
    output.str("");
    enabled.begin_step();
    enabled.report();
    if (rank == 0)
    {
        require(output.str().find("step=2 aggregation=rank_max") != std::string::npos);
        require(output.str().find("hypre_solve_calls_min=0 hypre_solve_calls_max=0") != std::string::npos);
    }
    std::cout.rdbuf(original);
    if (rank == 0) std::cout << "PASS: outlet profiler, real MPI ranks=" << ranks
                             << ", CUDA synchronization stubbed\n";
    MPI_Finalize();
}
