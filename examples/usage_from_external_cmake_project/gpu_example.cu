// A GPU program in its own CMake project, built against an installed MARS (find_package(Mars)): the
// appendix of the step-by-step tutorial. It generates the unit cube on every rank, lets cornerstone
// distribute it, and sums the volumes of the rank's own hexahedra on the GPU. The total must be 1.
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/utils/mars_generate_cube.hpp"
#include "backend/distributed/unstructured/utils/mars_gpu_reduce.hpp"
#include <thrust/functional.h>
#include <thrust/iterator/counting_iterator.h>
#include <mpi.h>
#include <cmath>
#include <cstdio>

using Domain = mars::ElementDomain<mars::HexTag, double, uint64_t, cstone::execution::Gpu>;

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, numRanks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount > 0) cudaSetDevice(rank % deviceCount);

    double volume = 0;
    {
        // 16^3 hexahedra of the unit cube; each rank generates its part, cornerstone redistributes them
        const size_t n = 16;
        [[maybe_unused]] auto [nodes, elements, x, y, z, conn] = mars::generateBoxElementPartition<double, uint64_t>(
            {n, n, n}, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 0.0, rank, numRanks);
        Domain::HostCoordsTuple coords{std::move(x), std::move(y), std::move(z)};
        Domain::HostConnectivityTuple connectivity{std::move(conn[0]), std::move(conn[1]), std::move(conn[2]),
                                                   std::move(conn[3]), std::move(conn[4]), std::move(conn[5]),
                                                   std::move(conn[6]), std::move(conn[7])};
        Domain domain(coords, connectivity, rank, numRanks, 64, true);

        // the elements are axis-aligned boxes: the volume is the product of the diagonal's extents
        const auto& c   = domain.getElementToNodeConnectivity();
        const uint64_t *c0 = std::get<0>(c).data(), *c6 = std::get<6>(c).data();
        const double *px = domain.getNodeX().data(), *py = domain.getNodeY().data(), *pz = domain.getNodeZ().data();
        const double local = mars::transformReduce(
            thrust::counting_iterator<size_t>(domain.startIndex()), thrust::counting_iterator<size_t>(domain.endIndex()),
            [=] __host__ __device__(size_t e) {
                return fabs(px[c6[e]] - px[c0[e]]) * fabs(py[c6[e]] - py[c0[e]]) * fabs(pz[c6[e]] - pz[c0[e]]);
            },
            0.0, thrust::plus<double>());
        MPI_Allreduce(&local, &volume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }

    const bool ok = std::fabs(volume - 1.0) < 1e-12;
    if (rank == 0) std::printf("MARS GPU example: volume = %.15f (%s)\n", volume, ok ? "ok" : "WRONG");
    MPI_Finalize();
    return ok ? 0 : 1;
}
