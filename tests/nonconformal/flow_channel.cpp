#include "examples/distributed/unstructured/mars_nonconformal_flow_runtime.hpp"
#include "channel_mesh.hpp"

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankId);
    MPI_Comm_size(MPI_COMM_WORLD, &rankCount);
    try
    {
        auto m = mesh(2, 1.);
        for (auto& x : m.x)
            x.x += 1;
        m.boundary.resize(m.x.size());
        for (int e = 0; e < int(m.elements.size()); ++e)
            for (int opp = 0; opp < 4; ++opp)
            {
                bool planes[6]{true, true, true, true, true, true};
                for (int k = 0; k < 4; ++k)
                    if (k != opp)
                    {
                        const auto& x = m.x[m.cells[4 * e + k]];
                        planes[0]     = planes[0] && x.x == 0;
                        planes[1]     = planes[1] && x.x == 2;
                        planes[2]     = planes[2] && x.y == 0;
                        planes[3]     = planes[3] && x.y == 1;
                        planes[4]     = planes[4] && x.z == 0;
                        planes[5]     = planes[5] && x.z == 1;
                    }
                for (int j = 0; j < 6; ++j)
                    if (planes[j]) m.outer.push_back({e, opp, j < 2 ? j : 2});
            }
        ncflow::Channel flow(m, true, .1, .1);
        check(flow.solve(m, 800, 1e-7), "analytic nonconformal plug channel");
        if (rankId == 0) std::puts("PASS: analytic stationary flow with pressure, coupled momentum and mass balance");
    }
    catch (const std::exception& e)
    {
        std::fprintf(stderr, "%s\n", e.what());
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
}
