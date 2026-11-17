#include "geometry_migration_checks.hpp"
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    try
    {
        for (bool staged : {false, true})
            amr_test::check_geometry_migration<mars::amr::LeppHostExecution>(MPI_COMM_WORLD, staged);
    }
    catch (const std::exception& e)
    {
        std::fprintf(stderr, "%s\n", e.what());
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
}
