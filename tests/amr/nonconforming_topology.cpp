#include "nonconforming_topology_checks.hpp"
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    try
    {
        for (bool staged : {false, true})
            amr_test::check_nonconforming_topology<mars::amr::LeppHostExecution>(MPI_COMM_WORLD, staged);
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
}
