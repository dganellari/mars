#include "field_migration_checks.hpp"
#include <cstdio>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    try
    {
        for (bool staged : {false, true})
            amr_test::check_field_sources<mars::amr::LeppHostExecution>(MPI_COMM_WORLD, staged);
        if (rank == 0) std::puts("Field migration: authority, empty ranks, reused fields and rejection controls PASS");
    }
    catch (const std::exception& e)
    {
        std::fprintf(stderr, "%s\n", e.what());
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
}
