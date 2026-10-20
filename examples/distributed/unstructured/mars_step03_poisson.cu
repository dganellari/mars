// Step 3: solve a PDE. Poisson's equation -Δu = f on the unit cube with u = 0 on the boundary, linear
// finite elements on tetrahedra, conjugate gradients on the GPU, on one or many ranks. The exact solution
// is known, u = sin(πx) sin(πy) sin(πz), so the program prints the error: refining the mesh shows the
// order of convergence. Walkthrough: docs/step-03-poisson.md.
//
//   python3 scripts/generate_tet_cube.py --nx 16 --ny 16 --nz 16 --output tet16
//   mpirun -np 4 ./mars_step03_poisson --mesh=tet16

#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_fem.hpp"
#include "backend/distributed/unstructured/solvers/mars_cg_solver.hpp"
#include <thrust/execution_policy.h>
#include <thrust/fill.h>
#include <thrust/for_each.h>
#include <thrust/functional.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/reduce.h>
#include <thrust/transform_reduce.h>
#include <mpi.h>
#include <cmath>
#include <cstdio>
#include <string>

using namespace mars;
using namespace mars::fem;

__host__ __device__ inline double exactSolution(double x, double y, double z)
{
    return sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z);
}

// f = -Δu for the exact solution above. The assembler calls it at its quadrature points.
struct Source
{
    __host__ __device__ double operator()(double x, double y, double z) const
    {
        return 3.0 * M_PI * M_PI * exactSolution(x, y, z);
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
        if (rank == 0) std::printf("usage: %s --mesh=DIR   (a tetrahedral mesh of the unit cube)\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    bool converged = false;
    {
        // 1. The mesh, as in step 1, now of tetrahedra.
        using Domain = ElementDomain<TetTag, double, uint64_t, cstone::execution::Gpu>;
        Domain domain(meshDir, rank, numRanks);
        const size_t nodeCount = domain.getNodeCount();

        // 2. Degrees of freedom: one per node. Each rank numbers its owned nodes 0 .. numOwned-1 and its
        //    ghost nodes after them, so the vectors below hold the owned values first.
        cstone::DeviceVector<int> nodeToDof(nodeCount);
        const int numOwned = buildDofMapping<uint64_t>(domain.getNodeOwnershipMap().data(), nodeToDof.data(),
                                                          nodeCount);
        long numGlobal = numOwned;
        MPI_Allreduce(MPI_IN_PLACE, &numGlobal, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

        // 3. The stiffness matrix A (rows of the owned nodes) and the load vector b = ∫ f φ_i.
        SparseMatrix<int, double, cstone::execution::Gpu> A;
        TetStiffnessAssembler<double, uint64_t>().assemble(domain, A, nodeToDof.data(), numOwned);

        cstone::DeviceVector<double> b(numOwned);
        TetMassAssembler<double, uint64_t>().assembleRHS(domain, b, Source{}, nodeToDof.data(), numOwned);

        // 4. u = 0 on the boundary: a boundary row becomes the identity row with zero right-hand side.
        //    The boundary is the six faces of the unit cube.
        {
            const double *x = domain.getNodeX().data(), *y = domain.getNodeY().data(), *z = domain.getNodeZ().data();
            const int* dof    = nodeToDof.data();
            const int* rowPtr = A.rowOffsetsPtr();
            const int* col    = A.colIndicesPtr();
            double* val       = A.valuesPtr();
            double* rhs       = b.data();
            thrust::for_each(thrust::device, thrust::counting_iterator<size_t>(0),
                             thrust::counting_iterator<size_t>(nodeCount), [=] __device__(size_t n) {
                                 const int i = dof[n];
                                 if (i >= numOwned) return; // a ghost: its owner sets its row
                                 const double eps = 1e-9;
                                 const bool onBoundary = x[n] < eps || x[n] > 1 - eps || y[n] < eps ||
                                                         y[n] > 1 - eps || z[n] < eps || z[n] > 1 - eps;
                                 if (!onBoundary) return;
                                 for (int j = rowPtr[i]; j < rowPtr[i + 1]; ++j)
                                     val[j] = col[j] == i ? 1.0 : 0.0;
                                 rhs[i] = 0.0;
                             });
        }

        // 5. Conjugate gradients. On several ranks the matrix-vector product needs the ghost values of
        //    the search direction, which the domain's node halo exchange brings from their owners.
        ConjugateGradientSolver<double, int, cstone::execution::Gpu> cg(2000, 1e-10);
        cg.setVerbose(false);
        cg.setOwnedSize(numOwned);
        if (numRanks > 1)
        {
            const int* dofMap = nodeToDof.data();
            cg.setHaloExchangeCallback(
                [&domain, dofMap](cstone::DeviceVector<double>& p) { domain.exchangeNodeHalo(p, dofMap); });
        }
        cstone::DeviceVector<double> u(nodeCount);
        thrust::fill(thrust::device, u.data(), u.data() + nodeCount, 0.0);
        converged = cg.solve(A, b, u);

        // 6. The error at the owned nodes, against the exact solution.
        const double *x = domain.getNodeX().data(), *y = domain.getNodeY().data(), *z = domain.getNodeZ().data();
        const int* dof   = nodeToDof.data();
        const double* uh = u.data();
        double maxError  = thrust::transform_reduce(
            thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(nodeCount),
            [=] __host__ __device__(size_t n) {
                const int i = dof[n];
                return i < numOwned ? fabs(uh[i] - exactSolution(x[n], y[n], z[n])) : 0.0;
            },
            0.0, thrust::maximum<double>());
        MPI_Allreduce(MPI_IN_PLACE, &maxError, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        if (rank == 0)
            std::printf("unknowns = %ld\nCG iterations = %d (%s)\nmax error = %.6e\n", numGlobal, cg.getIterations(),
                        converged ? "converged" : "NOT converged", maxError);
    }

    MPI_Finalize();
    return converged ? 0 : 1;
}
