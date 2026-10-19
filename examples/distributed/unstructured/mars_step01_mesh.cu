// Step 1: the mesh. Read a hexahedral mesh onto the GPUs and see how MARS splits it between ranks.
// Walkthrough: docs/step-01-mesh.md.
//
//   python3 scripts/generate_hex_cube.py -n 16 --output cube16
//   mpirun -np 4 ./mars_step01_mesh --mesh=cube16

#include "backend/distributed/unstructured/domain.hpp"
#include <thrust/reduce.h>
#include <mpi.h>
#include <cstdio>
#include <string>
#include <vector>

using namespace mars;

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    abortAllRanksOnUncaughtException();

    int rank, numRanks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    // One GPU per rank, chosen before anything touches the device.
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount > 0) cudaSetDevice(rank % deviceCount);

    std::string meshDir;
    for (int i = 1; i < argc; ++i)
        if (std::string(argv[i]).rfind("--mesh=", 0) == 0) meshDir = std::string(argv[i]).substr(7);
    if (meshDir.empty())
    {
        if (rank == 0) std::printf("usage: %s --mesh=DIR\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    {
        // Hexahedra, double precision, 64-bit space-filling-curve keys, stored on the GPU. The constructor
        // reads the mesh, sorts the elements along the space-filling curve and gives each rank one piece
        // of the curve, plus a halo: the neighbouring elements owned by other ranks.
        using Domain = ElementDomain<HexTag, double, uint64_t, cstone::execution::Gpu>;
        Domain domain(meshDir, rank, numRanks);

        // Every element a rank holds is a column entry of the connectivity: its own elements are the
        // range [startIndex, endIndex), the halo elements are before and after it.
        const auto& conn   = domain.getElementToNodeConnectivity();
        long held          = std::get<0>(conn).size();
        long owned         = domain.localElementCount();
        long halo          = held - owned;

        // Nodes: every node of every element held, numbered locally from 0. A node shared with another
        // rank is owned by exactly one of them; the ownership map says which (1 = mine, 0 = a ghost).
        long nodes         = domain.getNodeCount();
        const auto& owner  = domain.getNodeOwnershipMap();
        long ownedNodes    = thrust::reduce(thrust::device, owner.data(), owner.data() + nodes, 0L);

        // Counts per rank, gathered so rank 0 prints them in order.
        long mine[4] = {owned, halo, nodes, ownedNodes};
        std::vector<long> all(4 * numRanks);
        MPI_Gather(mine, 4, MPI_LONG, all.data(), 4, MPI_LONG, 0, MPI_COMM_WORLD);

        // Owned elements and owned nodes add up to the whole mesh: nothing is counted twice.
        long global[2] = {owned, ownedNodes};
        MPI_Allreduce(MPI_IN_PLACE, global, 2, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

        if (rank == 0)
        {
            std::printf("rank | elements | halo | nodes | owned nodes\n");
            for (int r = 0; r < numRanks; ++r)
                std::printf("%4d | %8ld | %4ld | %5ld | %11ld\n", r, all[4 * r], all[4 * r + 1], all[4 * r + 2],
                            all[4 * r + 3]);
            std::printf("mesh: %ld elements, %ld nodes\n", global[0], global[1]);
        }
    }

    MPI_Finalize();
    return 0;
}
