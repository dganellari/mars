// Step 6: the heat equation, a first time-dependent problem. ∂u/∂t = Δu on the unit cube, u = 0 on the
// boundary, starting from u0 = sin(πx) sin(πy) sin(πz), whose exact solution decays as u = e^(-3π²t) u0.
// Implicit Euler in time, (M + Δt K) u^(n+1) = M u^n: both matrices come from element kernels as in
// step 5, are assembled once, and serve every time step. Walkthrough: docs/step-06-heat.md.
//
//   python3 scripts/generate_hex_cube.py -n 16 --output cube16
//   mpirun -np 4 ./mars_step06_heat --mesh=cube16 --dt=0.01

#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_fem.hpp"
#include "backend/distributed/unstructured/fem/mars_element_assembly.hpp"
#include "backend/distributed/unstructured/solvers/mars_cg_solver.hpp"
#include "backend/distributed/unstructured/solvers/mars_spmv_plan.hpp"
#include <thrust/execution_policy.h>
#include <thrust/for_each.h>
#include <thrust/functional.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/transform_reduce.h>
#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

using namespace mars;
using namespace mars::fem;

__host__ __device__ inline double initialValue(double x, double y, double z)
{
    return sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z);
}

// The element matrix of M + Δt K: mass plus Δt times stiffness.
template<typename ElementTag>
struct HeatElement
{
    static constexpr int N = ElementTag::NodesPerElement;
    double dt;

    __host__ __device__ void operator()(const double (&X)[N][3], double (&K)[N][N], double (&)[N]) const
    {
        ElementValues<ElementTag, double> v;
        for (int q = 0; q < v.numQuadPoints; ++q)
        {
            v.reinit(X, q);
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j)
                    K[i][j] += (v.phi[i] * v.phi[j] + dt * (v.grad[i][0] * v.grad[j][0] + v.grad[i][1] * v.grad[j][1] +
                                                            v.grad[i][2] * v.grad[j][2])) *
                               v.JxW;
        }
    }
};

// The element mass matrix M_ij = ∫ φ_i φ_j, for the right-hand side M u^n.
template<typename ElementTag>
struct MassElement
{
    static constexpr int N = ElementTag::NodesPerElement;

    __host__ __device__ void operator()(const double (&X)[N][3], double (&K)[N][N], double (&)[N]) const
    {
        ElementValues<ElementTag, double> v;
        for (int q = 0; q < v.numQuadPoints; ++q)
        {
            v.reinit(X, q);
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j)
                    K[i][j] += v.phi[i] * v.phi[j] * v.JxW;
        }
    }
};

// true if CG converged at every step
template<typename ElementTag>
bool solve(const std::string& meshDir, double dt, double endTime, int rank, int numRanks)
{
    using Domain = ElementDomain<ElementTag, double, uint64_t, cstone::execution::Gpu>;
    using Matrix = SparseMatrix<int, double, cstone::execution::Gpu>;
    using Vector = cstone::DeviceVector<double>;
    Domain domain(meshDir, rank, numRanks);
    const size_t nodeCount = domain.getNodeCount();

    cstone::DeviceVector<int> nodeToDof(nodeCount);
    const int numOwned = buildDofMapping<uint64_t>(domain.getNodeOwnershipMap().data(), nodeToDof.data(), nodeCount);
    long numGlobal = numOwned;
    MPI_Allreduce(MPI_IN_PLACE, &numGlobal, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

    // the fewest steps of at most dt that end exactly at endTime (dt shrinks to endTime / steps)
    const int steps = std::max(1, int(std::ceil(endTime / dt - 1e-9)));
    dt              = endTime / steps;

    // 1. The two matrices, assembled once. A = M + Δt K is the matrix of every step; M gives its
    //    right-hand side. The element vectors are not used here.
    Matrix A, M;
    Vector unused;
    assembleElements(domain, A, unused, nodeToDof.data(), numOwned, HeatElement<ElementTag>{dt});
    assembleElements(domain, M, unused, nodeToDof.data(), numOwned, MassElement<ElementTag>{});

    const double *x = domain.getNodeX().data(), *y = domain.getNodeY().data(), *z = domain.getNodeZ().data();
    const int* dof = nodeToDof.data();

    // 2. u = 0 on the boundary: identity rows in A (the right-hand side is zeroed at every step).
    cstone::DeviceVector<uint8_t> boundary(numOwned, 0);
    {
        const int* rowPtr = A.rowOffsetsPtr();
        const int* col    = A.colIndicesPtr();
        double* val       = A.valuesPtr();
        uint8_t* onB      = boundary.data();
        thrust::for_each(thrust::device, thrust::counting_iterator<size_t>(0),
                         thrust::counting_iterator<size_t>(nodeCount), [=] __device__(size_t n) {
                             const int i = dof[n];
                             if (i >= numOwned) return;
                             const double eps = 1e-9;
                             if (!(x[n] < eps || x[n] > 1 - eps || y[n] < eps || y[n] > 1 - eps || z[n] < eps ||
                                   z[n] > 1 - eps))
                                 return;
                             onB[i] = 1;
                             for (int j = rowPtr[i]; j < rowPtr[i + 1]; ++j)
                                 val[j] = col[j] == i ? 1.0 : 0.0;
                         });
    }

    // 3. The initial condition at every local node, ghosts included, indexed by DOF.
    Vector u(nodeCount), b(numOwned);
    {
        double* uh = u.data();
        thrust::for_each(thrust::device, thrust::counting_iterator<size_t>(0),
                         thrust::counting_iterator<size_t>(nodeCount),
                         [=] __device__(size_t n) { uh[dof[n]] = initialValue(x[n], y[n], z[n]); });
    }

    // 4. Time steps. M u^n reads the ghost values of u, so they are refreshed from their owners first; the
    //    product is set up once (SpmvPlan) and applied at every step.
    ConjugateGradientSolver<double, int, cstone::execution::Gpu> cg(2000, 1e-10);
    cg.setVerbose(false);
    cg.setOwnedSize(numOwned);
    if (numRanks > 1)
        cg.setHaloExchangeCallback([&domain, dof](Vector& p) { domain.exchangeNodeHalo(p, dof); });
    cusparseHandle_t handle;
    cusparseCreate(&handle);
    bool converged = true;
    int iterations = 0;
    {
        SpmvPlan<double, Matrix, Vector> massTimes(handle, M, u, b);
        for (int step = 0; step < steps; ++step)
        {
            if (numRanks > 1) domain.exchangeNodeHalo(u, dof);
            massTimes.apply(u, b);
            double* rhs        = b.data();
            const uint8_t* onB = boundary.data();
            thrust::for_each(thrust::device, thrust::counting_iterator<int>(0), thrust::counting_iterator<int>(numOwned),
                             [=] __device__(int i) {
                                 if (onB[i]) rhs[i] = 0.0;
                             });
            converged = cg.solve(A, b, u) && converged;
            iterations += cg.getIterations();
        }
    }
    cusparseDestroy(handle);

    // 5. The error at the owned nodes against the exact solution at the end time.
    const double decay = std::exp(-3.0 * M_PI * M_PI * steps * dt);
    const double* uh   = u.data();
    double maxError    = thrust::transform_reduce(
        thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(nodeCount),
        [=] __host__ __device__(size_t n) {
            const int i = dof[n];
            return i < numOwned ? fabs(uh[i] - decay * initialValue(x[n], y[n], z[n])) : 0.0;
        },
        0.0, thrust::maximum<double>());
    MPI_Allreduce(MPI_IN_PLACE, &maxError, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if (rank == 0)
        std::printf("elements = %s\nunknowns = %ld\ntime steps = %d (dt = %g, t = %g)\nCG iterations = %d (%s)\n"
                    "max error = %.6e\n",
                    ElementTag::Name, numGlobal, steps, dt, steps * dt, iterations,
                    converged ? "converged" : "NOT converged", maxError);
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
    double dt = 0.01, endTime = 0.1;
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        if (a.rfind("--mesh=", 0) == 0) meshDir = a.substr(7);
        else if (a.rfind("--elements=", 0) == 0) elements = a.substr(11);
        else if (a.rfind("--dt=", 0) == 0) dt = std::stod(a.substr(5));
        else if (a.rfind("--time=", 0) == 0) endTime = std::stod(a.substr(7));
    }
    // finite, positive, and at most a million steps, before the step count becomes an int
    const bool timesOk = std::isfinite(dt) && std::isfinite(endTime) && dt > 0 && endTime > 0 &&
                         std::ceil(endTime / dt - 1e-9) <= 1e6;
    if (meshDir.empty() || (elements != "hex" && elements != "tet") || !timesOk)
    {
        if (rank == 0)
            std::printf("usage: %s --mesh=DIR [--elements=hex|tet] [--dt=0.01] [--time=0.1]   (a mesh of the unit cube)\n",
                        argv[0]);
        MPI_Finalize();
        return 1;
    }

    const bool converged = elements == "tet" ? solve<TetTag>(meshDir, dt, endTime, rank, numRanks)
                                             : solve<HexTag>(meshDir, dt, endTime, rank, numRanks);

    MPI_Finalize();
    return converged ? 0 : 1;
}
