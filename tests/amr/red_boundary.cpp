#include "red_boundary_checks.hpp"
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    try
    {
        for (bool staged : {false, true})
            amr_test::check_red_boundary<mars::amr::LeppHostExecution>(MPI_COMM_WORLD, staged,
                                                                       argc == 2 ? argv[1] : "");
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
}
