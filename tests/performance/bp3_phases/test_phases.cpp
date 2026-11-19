#include "phases.hpp"
#include <fstream>
#include <iostream>

namespace
{
int waits = 0;
void sync() { ++waits; }
void require(bool ok, const char* message)
{
    if (!ok)
    {
        std::cerr << message << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}
} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    require(argc == 3, "usage: phase_test mode receipt-prefix");
    const std::string mode = argv[1], prefix = argv[2];
    unsetenv("BP3_PHASE_PREFIX");
    unsetenv("BP3_PHASE_DETAIL");
    if (mode != "off") setenv("BP3_PHASE_PREFIX", prefix.c_str(), 1);
    if (mode == "detail") setenv("BP3_PHASE_DETAIL", "1", 1);
    bp3_phase::begin(sync);
    if (mode == "unbalanced")
    {
        bp3_phase::Scope coarse(bp3_phase::Coarse);
        bp3_phase::end();
        return 1;
    }
    {
        bp3_phase::Scope cycle(bp3_phase::Cycle, true, true);
        {
            bp3_phase::Scope smoother(bp3_phase::Smoother, true, true);
            bp3_phase::Scope op(bp3_phase::Operator, true, true);
        }
        bp3_phase::entry_wait();
        {
            bp3_phase::Scope coarse(bp3_phase::Coarse, false, true);
            bp3_phase::Scope hypre(bp3_phase::Hypre, false, true);
        }
    }
    bp3_phase::end();
    bp3_phase::write();
    using namespace bp3_phase;
    require(waits == (mode == "off" ? 0 : mode == "detail" ? 13 : 7), "Wrong completion boundaries");
    require(depth == 0 && !active, "Capture did not finish");
    if (mode == "off")
    {
        require(stats[Solve].calls == 0, "Disabled capture ran");
        std::ifstream receipt(prefix + "-rank0.json");
        require(!receipt, "Disabled capture wrote a file");
    }
    else
    {
        require(stats[Solve].calls == 1 && stats[Hypre].calls == 1, "Missing calls");
        require(stats[Operator].calls == (mode == "detail" ? 1 : 0), "Wrong detail mode");
        long long exclusive = 0;
        for (auto& s : stats)
        {
            require(s.exclusive >= 0 && s.exclusive <= s.inclusive, "Invalid nested accounting");
            exclusive += s.exclusive;
        }
        require(exclusive == stats[Solve].inclusive, "Exclusive intervals do not partition the solve");
        require(stats[Coarse].inclusive >= stats[Hypre].inclusive, "Hypre escaped coarse interval");
    }
    MPI_Finalize();
    std::cout << "PASS phase boundaries " << mode << '\n';
}
