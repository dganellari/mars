#include "geometry_migration_checks.hpp"
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, devices = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    try
    {
        int local_rank;
        MPI_Comm local;
        MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &local);
        MPI_Comm_rank(local, &local_rank);
        cudaCheckCall(cudaGetDeviceCount(&devices));
        mars::amr::require_collectively(devices > 0 && (devices == 1 || local_rank < devices),
                                       "No GPU for geometry migration checks", MPI_COMM_WORLD);
        const int device = devices == 1 ? 0 : local_rank;
        cudaCheckCall(cudaSetDevice(device));
        MPI_Comm_free(&local);
        char host[MPI_MAX_PROCESSOR_NAME], pci[32];
        int length;
        MPI_Get_processor_name(host, &length);
#if defined(__HIPCC__)
        cudaCheckCall(hipDeviceGetPCIBusId(pci, sizeof(pci), device));
#else
        cudaCheckCall(cudaDeviceGetPCIBusId(pci, sizeof(pci), device));
#endif
        std::printf("[amr-geometry-device] rank=%d host=%s pci=%s\n", rank, host, pci);
        std::fflush(stdout);
        for (bool staged : {false, true})
            amr_test::check_geometry_migration<mars::amr::LeppDeviceExecution>(MPI_COMM_WORLD, staged);
        cudaCheckCall(cudaDeviceSynchronize());
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
}
