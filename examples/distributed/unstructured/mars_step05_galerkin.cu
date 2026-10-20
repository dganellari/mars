// Step 5: your own finite-element operator. The Poisson problem of steps 3 and 4, -Δu = f on the unit cube
// with u = sin(πx) sin(πy) sin(πz), but the element matrix and vector are written here, in a dozen lines,
// and MARS assembles them on the GPU: linear Galerkin elements on hexahedra (Q1), or on tetrahedra (P1) with
// --elements=tet, where the result must match step 3. Walkthrough: docs/step-05-galerkin.md.
//
//   python3 scripts/generate_hex_cube.py -n 16 --output cube16
//   mpirun -np 4 ./mars_step05_galerkin --mesh=cube16

#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_fem.hpp"
#include "backend/distributed/unstructured/fem/mars_element_assembly.hpp"
#include "backend/distributed/unstructured/solvers/mars_cg_solver.hpp"
#include <thrust/execution_policy.h>
#include <thrust/fill.h>
#include <thrust/for_each.h>
#include <thrust/functional.h>
#include <thrust/iterator/counting_iterator.h>
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

// The element matrix K_ij = ∫ ∇φ_i · ∇φ_j and vector F_i = ∫ f φ_i of -Δu = f, by quadrature. ElementValues
// gives, at each quadrature point, the basis values φ_i, their gradients and the weight times |det J|.
template<typename ElementTag>
struct PoissonElement
{
    static constexpr int N = ElementTag::NodesPerElement;

    __host__ __device__ void operator()(const double (&X)[N][3], double (&K)[N][N], double (&F)[N]) const
    {
        ElementValues<ElementTag, double> v;
        for (int q = 0; q < v.numQuadPoints; ++q)
        {
            v.reinit(X, q);
            const double f = 3.0 * M_PI * M_PI * exactSolution(v.x[0], v.x[1], v.x[2]);
            for (int i = 0; i < N; ++i)
            {
                F[i] += f * v.phi[i] * v.JxW;
                for (int j = 0; j < N; ++j)
                    K[i][j] += (v.grad[i][0] * v.grad[j][0] + v.grad[i][1] * v.grad[j][1] + v.grad[i][2] * v.grad[j][2]) * v.JxW;
            }
        }
    }
};

// Steps 3 and 4 with the matrix and right-hand side from assembleElements; true if CG converged.
template<typename ElementTag>
bool solve(const std::string& meshDir, int rank, int numRanks)
{
    using Domain = ElementDomain<ElementTag, double, uint64_t, cstone::execution::Gpu>;
    Domain domain(meshDir, rank, numRanks);
    const size_t nodeCount = domain.getNodeCount();

    cstone::DeviceVector<int> nodeToDof(nodeCount);
    const int numOwned = buildDofMapping<uint64_t>(domain.getNodeOwnershipMap().data(), nodeToDof.data(), nodeCount);
    long numGlobal = numOwned;
    MPI_Allreduce(MPI_IN_PLACE, &numGlobal, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

    // The one new line: assemble A and b from the element kernel above.
    SparseMatrix<int, double, cstone::execution::Gpu> A;
    cstone::DeviceVector<double> b;
    assembleElements(domain, A, b, nodeToDof.data(), numOwned, PoissonElement<ElementTag>{});

    const double *x = domain.getNodeX().data(), *y = domain.getNodeY().data(), *z = domain.getNodeZ().data();
    const int* dof = nodeToDof.data();
    {
        const int* rowPtr = A.rowOffsetsPtr();
        const int* col    = A.colIndicesPtr();
        double* val       = A.valuesPtr();
        double* rhs       = b.data();
        thrust::for_each(thrust::device, thrust::counting_iterator<size_t>(0),
                         thrust::counting_iterator<size_t>(nodeCount), [=] __device__(size_t n) {
                             const int i = dof[n];
                             if (i >= numOwned) return;
                             const double eps = 1e-9;
                             const bool onBoundary = x[n] < eps || x[n] > 1 - eps || y[n] < eps || y[n] > 1 - eps ||
                                                     z[n] < eps || z[n] > 1 - eps;
                             if (!onBoundary) return;
                             for (int j = rowPtr[i]; j < rowPtr[i + 1]; ++j)
                                 val[j] = col[j] == i ? 1.0 : 0.0;
                             rhs[i] = 0.0;
                         });
    }

    ConjugateGradientSolver<double, int, cstone::execution::Gpu> cg(2000, 1e-10);
    cg.setVerbose(false);
    cg.setOwnedSize(numOwned);
    if (numRanks > 1)
        cg.setHaloExchangeCallback(
            [&domain, dof](cstone::DeviceVector<double>& p) { domain.exchangeNodeHalo(p, dof); });
    cstone::DeviceVector<double> u(nodeCount);
    thrust::fill(thrust::device, u.data(), u.data() + nodeCount, 0.0);
    const bool converged = cg.solve(A, b, u);

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
        std::printf("elements = %s\nunknowns = %ld\nCG iterations = %d (%s)\nmax error = %.6e\n", ElementTag::Name,
                    numGlobal, cg.getIterations(), converged ? "converged" : "NOT converged", maxError);
    return converged;
}

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

    std::string meshDir, elements = "hex";
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        if (a.rfind("--mesh=", 0) == 0) meshDir = a.substr(7);
        else if (a.rfind("--elements=", 0) == 0) elements = a.substr(11);
    }
    if (meshDir.empty() || (elements != "hex" && elements != "tet"))
    {
        if (rank == 0) std::printf("usage: %s --mesh=DIR [--elements=hex|tet]   (a mesh of the unit cube)\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    const bool converged =
        elements == "tet" ? solve<TetTag>(meshDir, rank, numRanks) : solve<HexTag>(meshDir, rank, numRanks);

    MPI_Finalize();
    return converged ? 0 : 1;
}
