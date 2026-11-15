#include "red_boundary_checks.hpp"
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, devices = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    try
    {
        cudaCheckCall(cudaGetDeviceCount(&devices));
        mars::amr::require_collectively(devices > 0, "No GPU for red boundary checks", MPI_COMM_WORLD);
        cudaCheckCall(cudaSetDevice(rank % devices));
        for (bool staged : {false, true})
            amr_test::check_red_boundary<mars::amr::LeppDeviceExecution>(MPI_COMM_WORLD, staged);
        cudaCheckCall(cudaDeviceSynchronize());
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
}
