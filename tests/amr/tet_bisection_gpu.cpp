#include "tet_checks.hpp"
#include <cstdlib>
#include <string>

int main(int argc, char** argv)
{
    using namespace amr_test;
    MPI_Init(&argc, &argv);
    int rank = 0, ranks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    try
    {
        select_device(rank);
        bool full_halo = false, staged = false;
        std::string exodus, public_mesh;
        for (int i = 1; i < argc; ++i)
        {
            if (std::string(argv[i]) == "--full-halo")
                full_halo = true;
            else if (std::string(argv[i]) == "--staged")
                staged = true;
            else if (std::string(argv[i]) == "--public-mesh" && i + 1 < argc)
                public_mesh = argv[++i];
            else if (std::string(argv[i]) == "--exodus" && i + 1 < argc)
                exodus = argv[++i];
            else
                require(false, "Unknown test argument");
        }
        setenv("MARS_HALO_EXCHANGE", full_halo ? "full" : "needed", 1);
        if (rank == 0) std::cout << "Node halo mode: " << (full_halo ? "full" : "default one-layer") << "\n";
        require(exodus.empty() || public_mesh.empty(), "Select one imported mesh fixture");
        if (!public_mesh.empty())
            check_public_manager(public_mesh, rank, ranks, staged);
        else if (exodus.empty())
            run_gpu_checks(rank, ranks, staged);
        else
            check_exodus_manager(exodus, rank, ranks, staged);
        if (rank == 0)
            std::cout << "Conforming tetrahedral AMR " << (exodus.empty() ? "" : "Exodus ") << "GPU/MPI gate PASS on "
                      << ranks << " ranks\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << "rank " << rank << ": " << error.what() << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    MPI_Finalize();
    return 0;
}
