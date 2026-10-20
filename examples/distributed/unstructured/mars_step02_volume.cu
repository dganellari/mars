// Step 2: a first loop over the elements on the GPU. Compute every element's volume from its corner
// coordinates, then the volume of the whole mesh and its first moment, integral of x over the domain.
// Walkthrough: docs/step-02-volume.md.
//
//   python3 scripts/generate_hex_cube.py -n 16 --output cube16
//   mpirun -np 4 ./mars_step02_volume --mesh=cube16

#include "backend/distributed/unstructured/domain.hpp"
#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <mpi.h>
#include <cmath>
#include <cstdio>
#include <string>

using namespace mars;

// Volume of the tetrahedron (a, b, c, d) from the triple product.
__host__ __device__ inline double tetVolume(const double* a, const double* b, const double* c, const double* d)
{
    double u[3], v[3], w[3];
    for (int k = 0; k < 3; ++k)
    {
        u[k] = b[k] - a[k];
        v[k] = c[k] - a[k];
        w[k] = d[k] - a[k];
    }
    double det = u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0]) + u[2] * (v[0] * w[1] - v[1] * w[0]);
    return fabs(det) / 6.0;
}

// Volume (moment = false) or integral of x (moment = true) of element e. The hexahedron splits into six
// tetrahedra around its diagonal from corner 0 to corner 6: exact for elements with flat faces. The struct
// holds raw device pointers, because it is copied to the GPU and called there.
struct ElementIntegral
{
    const uint64_t* corner[8]; // connectivity columns: local node id of each corner
    const double *x, *y, *z;   // node coordinates
    bool moment;

    __host__ __device__ double operator()(size_t e) const
    {
        double p[8][3];
        for (int k = 0; k < 8; ++k)
        {
            const uint64_t n = corner[k][e];
            p[k][0] = x[n];
            p[k][1] = y[n];
            p[k][2] = z[n];
        }
        const int tets[6][2] = {{1, 2}, {2, 3}, {3, 7}, {7, 4}, {4, 5}, {5, 1}};
        double sum = 0.0;
        for (int t = 0; t < 6; ++t)
        {
            const double* b = p[tets[t][0]];
            const double* c = p[tets[t][1]];
            double v = tetVolume(p[0], b, c, p[6]);
            // a tetrahedron's centroid integrates any linear function exactly
            double cx = 0.25 * (p[0][0] + b[0] + c[0] + p[6][0]);
            sum += moment ? v * cx : v;
        }
        return sum;
    }
};

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    abortAllRanksOnUncaughtException();

    int rank, numRanks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

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
        // true: keep the node coordinates exactly as read. By default MARS keeps only each node's
        // space-filling-curve key and decodes coordinates from it, rounded to the curve's grid.
        using Domain = ElementDomain<HexTag, double, uint64_t, cstone::execution::Gpu>;
        Domain domain(meshDir, rank, numRanks, true);

        // The 8 corner columns of the connectivity (local node ids) and the node coordinates, on the GPU.
        const auto& conn = domain.getElementToNodeConnectivity();
        ElementIntegral volume{{std::get<0>(conn).data(), std::get<1>(conn).data(), std::get<2>(conn).data(),
                                std::get<3>(conn).data(), std::get<4>(conn).data(), std::get<5>(conn).data(),
                                std::get<6>(conn).data(), std::get<7>(conn).data()},
                               domain.getNodeX().data(), domain.getNodeY().data(), domain.getNodeZ().data(),
                               false};
        ElementIntegral momentX = volume;
        momentX.moment          = true;

        // Only the rank's own elements, [startIndex, endIndex): the halo elements belong to other ranks,
        // and adding them here would count them twice in the global sum.
        thrust::counting_iterator<size_t> first(domain.startIndex()), last(domain.endIndex());
        double local[2];
        local[0] = thrust::transform_reduce(thrust::device, first, last, volume, 0.0, thrust::plus<double>());
        local[1] = thrust::transform_reduce(thrust::device, first, last, momentX, 0.0, thrust::plus<double>());

        double global[2];
        MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        if (rank == 0) std::printf("volume = %.15f\nintegral of x = %.15f\n", global[0], global[1]);
    }

    MPI_Finalize();
    return 0;
}
