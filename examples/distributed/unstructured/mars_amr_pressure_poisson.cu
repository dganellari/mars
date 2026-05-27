// AMR pressure-Poisson driver. Solves
//
//     -Laplacian(p) = -div(v . grad v)
//
// on the input mesh, where v is a prescribed analytic velocity field
//
//     v(x, y, z) = (Uinf, 0, 0) + gamma * (-(y-yc), (x-xc), 0) / ((x-xc)^2 + (y-yc)^2 + eps)
//
// i.e. uniform free-stream plus a 2D potential vortex aligned with a vertical
// line through (xc, yc). Center configurable via --vortex-center=xc,yc,zc
// (zc is currently unused — vortex is z-translation-invariant).
// For the NASA wing-tip mesh this approximates the kinematic effect of a
// trailing tip vortex; for a unit cube it's a smooth test problem.
//
// Same CVFEM operator + AMR loop as mars_amr_cvfem_graph; only the RHS
// integration kernel differs. Output: VTU/PVTU/PVD frames of the pressure
// field at each AMR level, ParaView-loadable as a movie.

#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_fem.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_assembler.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_utils.hpp"
#include "backend/distributed/unstructured/fem/mars_sparsity_builder.hpp"
#include "backend/distributed/unstructured/fem/mars_perf_counters.hpp"
#include "backend/distributed/unstructured/fem/mars_sparse_matrix.hpp"
#include "backend/distributed/unstructured/solvers/mars_cg_solver.hpp"
#ifdef MARS_ENABLE_HYPRE
#include "backend/distributed/unstructured/solvers/mars_hypre_pcg_solver.hpp"
#endif
#include "backend/distributed/unstructured/amr/mars_amr.hpp"
#include "backend/distributed/unstructured/utils/mars_vtu_parallel_writer.hpp"

#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/inner_product.h>
#include <thrust/transform.h>
#include <thrust/system/cuda/execution_policy.h>

#include <memory>
#include <mpi.h>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <limits>

using namespace mars;
using namespace mars::fem;
using namespace mars::amr;

// PhaseTimer: cudaDeviceSynchronize + MPI_Barrier on entry/lap so wall-clock
// reflects all-rank work, not just rank-0. Stores millisecond timings keyed
// by phase name; printed at end of each level.
struct PhaseTimer
{
    using clk = std::chrono::high_resolution_clock;
    clk::time_point t0;
    std::vector<std::pair<std::string, float>> phases;
    bool sync;

    explicit PhaseTimer(bool sync_ = true) : sync(sync_) { reset(); }

    void reset()
    {
        if (sync) { cudaDeviceSynchronize(); MPI_Barrier(MPI_COMM_WORLD); }
        t0 = clk::now();
    }

    void lap(const std::string& name)
    {
        if (sync) { cudaDeviceSynchronize(); MPI_Barrier(MPI_COMM_WORLD); }
        auto t1 = clk::now();
        phases.emplace_back(name, std::chrono::duration<float, std::milli>(t1 - t0).count());
        t0 = t1;
    }

    void report(int rank, const std::string& header)
    {
        if (rank != 0) return;
        std::cout << "  [" << header << "] phase breakdown:\n";
        float total = 0;
        for (auto& [n, ms] : phases) total += ms;
        for (auto& [n, ms] : phases)
        {
            std::cout << "    " << std::left << std::setw(28) << n
                      << std::right << std::fixed << std::setprecision(2)
                      << std::setw(10) << ms << " ms ("
                      << std::setw(5) << std::setprecision(1)
                      << (100.0f * ms / std::max(total, 1e-3f)) << "%)\n";
        }
        std::cout << "    " << std::left << std::setw(28) << "TOTAL"
                  << std::right << std::fixed << std::setprecision(2)
                  << std::setw(10) << total << " ms\n";
    }
};

// GPU kernel: apply Dirichlet BC on boundary nodes (zero on all faces of bounding box)
template<typename RealType>
__global__ void applyBCKernel(const RealType* nodeX,
                              const RealType* nodeY,
                              const RealType* nodeZ,
                              const uint8_t* ownership,
                              const int* nodeToDof,
                              uint8_t* isBoundaryDof,
                              size_t numNodes,
                              RealType xmin,
                              RealType xmax,
                              RealType ymin,
                              RealType ymax,
                              RealType zmin,
                              RealType zmax,
                              RealType eps)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;

    if (ownership[i] == 1)
    {
        int dof = nodeToDof[i];
        if (dof >= 0)
        {
            RealType x = nodeX[i], y = nodeY[i], z = nodeZ[i];
            bool onBdry = (fabs(x - xmin) < eps || fabs(x - xmax) < eps || fabs(y - ymin) < eps ||
                           fabs(y - ymax) < eps || fabs(z - zmin) < eps || fabs(z - zmax) < eps);
            isBoundaryDof[dof] = onBdry ? 1 : 0;
        }
    }
}

// Prescribed analytic velocity field used as the source for the
// pressure-Poisson RHS. Free-stream along x plus a 2D potential vortex
// centered on the vertical line through (xc, yc). ∇·v = 0 by construction.
template<typename RealType>
__device__ inline void prescribedVelocity(RealType x, RealType y, RealType /*z*/,
                                          RealType xc, RealType yc,
                                          RealType Uinf, RealType gamma, RealType eps,
                                          RealType& vx, RealType& vy, RealType& vz)
{
    RealType dx = x - xc;
    RealType dy = y - yc;
    RealType r2 = dx * dx + dy * dy + eps;
    vx = Uinf + gamma * (-dy) / r2;
    vy = gamma * dx / r2;
    vz = RealType(0);
}

// Analytic f = -div(v . grad v) at the element centroid for the prescribed
// field above. Uses div(v . grad v) = (∂_i v_j)(∂_j v_i) since ∇·v = 0.
// For our 2D field (v_z = 0, no z-dependence):
//     f = -[ (∂_x v_x)^2 + (∂_y v_y)^2 + 2 (∂_x v_y)(∂_y v_x) ]
// Derivatives use dx = x-xc, dy = y-yc; vortex contribution dominates
// near r = sqrt(dx^2 + dy^2) = 0; eps regularizes the singularity.
template<typename RealType>
__device__ inline RealType pressurePoissonRhs(RealType x, RealType y, RealType /*z*/,
                                              RealType xc, RealType yc,
                                              RealType Uinf, RealType gamma, RealType eps)
{
    RealType dx = x - xc;
    RealType dy = y - yc;
    RealType r2 = dx * dx + dy * dy + eps;
    RealType r4 = r2 * r2;
    RealType dvx_dx = gamma * RealType(2) * dx * dy / r4;
    RealType dvy_dy = -dvx_dx;
    RealType dvx_dy = gamma * (dy * dy - dx * dx - eps) / r4;
    RealType dvy_dx = gamma * (dy * dy - dx * dx + eps) / r4;

    return -(dvx_dx * dvx_dx + dvy_dy * dvy_dy + RealType(2) * dvy_dx * dvx_dy);
}

// Pressure-Poisson RHS integration. Evaluates f at the element centroid
// (one-point quadrature, exact for constant f), multiplies by V_elem,
// scatters V_elem * f / 8 to each owned-corner DOF. Same lumped pattern
// as integrateConstantSourceKernel; only the source value is per-element.
template<typename KeyType, typename RealType>
__global__ void integratePressurePoissonRhsKernel(const KeyType* conn0,
                                                  const KeyType* conn1,
                                                  const KeyType* conn2,
                                                  const KeyType* conn3,
                                                  const KeyType* conn4,
                                                  const KeyType* conn5,
                                                  const KeyType* conn6,
                                                  const KeyType* conn7,
                                                  const RealType* nodeX,
                                                  const RealType* nodeY,
                                                  const RealType* nodeZ,
                                                  const int* nodeToDof,
                                                  const uint8_t* ownership,
                                                  RealType* rhs,
                                                  size_t numElements,
                                                  RealType xcVortex,
                                                  RealType ycVortex,
                                                  RealType Uinf,
                                                  RealType gamma,
                                                  RealType eps_vortex)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements) return;

    KeyType n[8] = {conn0[e], conn1[e], conn2[e], conn3[e],
                    conn4[e], conn5[e], conn6[e], conn7[e]};

    RealType x[8], y[8], z[8];
    for (int i = 0; i < 8; ++i)
    {
        x[i] = nodeX[n[i]];
        y[i] = nodeY[n[i]];
        z[i] = nodeZ[n[i]];
    }
    // Element centroid + axis-aligned-box volume estimate (consistent with
    // integrateConstantSourceKernel above).
    RealType xc = 0, yc = 0, zc = 0;
    for (int i = 0; i < 8; ++i) { xc += x[i]; yc += y[i]; zc += z[i]; }
    xc *= RealType(0.125); yc *= RealType(0.125); zc *= RealType(0.125);

    RealType xmin = x[0], xmax = x[0], ymin = y[0], ymax = y[0], zmin = z[0], zmax = z[0];
    for (int i = 1; i < 8; ++i)
    {
        if (x[i] < xmin) xmin = x[i]; if (x[i] > xmax) xmax = x[i];
        if (y[i] < ymin) ymin = y[i]; if (y[i] > ymax) ymax = y[i];
        if (z[i] < zmin) zmin = z[i]; if (z[i] > zmax) zmax = z[i];
    }
    RealType V = (xmax - xmin) * (ymax - ymin) * (zmax - zmin);

    RealType f = pressurePoissonRhs<RealType>(xc, yc, zc, xcVortex, ycVortex, Uinf, gamma, eps_vortex);
    RealType contrib = f * V * RealType(0.125);

    for (int i = 0; i < 8; ++i)
    {
        if (ownership[n[i]] == 1)
        {
            int dof = nodeToDof[n[i]];
            if (dof >= 0)
                atomicAdd(&rhs[dof], contrib);
        }
    }
}

// GPU kernel: zero boundary rows, set diagonal=1, rhs=0
template<typename RealType>
__global__ void enforceBCKernel(const uint8_t* isBoundaryDof,
                                const int* rowPtr,
                                const int* colInd,
                                const int* diagPtr,
                                RealType* values,
                                RealType* rhs,
                                int numDofs)
{
    int dof = blockIdx.x * blockDim.x + threadIdx.x;
    if (dof >= numDofs || !isBoundaryDof[dof]) return;

    // Zero entire row
    for (int j = rowPtr[dof]; j < rowPtr[dof + 1]; ++j)
        values[j] = RealType(0);

    // Set diagonal to 1
    int dp = diagPtr[dof];
    if (dp >= 0) values[dp] = RealType(1);

    // Zero RHS
    rhs[dof] = RealType(0);
}

// Linear solver selector. CG is the default; Hypre PCG+BoomerAMG is opt-in
// through --solver=hypre and only available when MARS is built with
// MARS_ENABLE_HYPRE=ON. The CG path is unchanged from the original driver.
enum class SolverKind { CG, Hypre };

// GPU-native Poisson solve on an AMR-managed domain. Pipeline:
//   1) GPU DOF mapping (buildDofMappingGpu)
//   2) GPU sparsity (CvfemSparsityBuilder::buildGraphSparsity, sized for owned DOFs)
//   3) GPU full 8x8 element assembly (assembleFull)
//   4) GPU bounding box reduction + MPI_Allreduce for global box
//   5) GPU BC marking + GPU BC enforcement (applyBCKernel + enforceBCKernel)
//   6) Linear solve (CG or Hypre PCG+BoomerAMG)
//   7) GPU DOF -> per-node scatter
// No host downloads of node/element data.
template<typename KeyType, typename RealType>
void solvePressurePoisson(ElementDomain<HexTag, RealType, KeyType, cstone::GpuTag>& domain,
                          RealType xcVortex,
                          RealType ycVortex,
                          RealType Uinf,
                          RealType gamma,
                          RealType eps_vortex,
                          CvfemKernelVariant kernelVariant,
                          int blockSize,
                          int maxIter,
                          RealType tolerance,
                          SolverKind solverKind,
                          int rank,
                          cstone::DeviceVector<RealType>& d_nodeSolution,
                          cstone::DeviceVector<RealType>& d_errorPerElement)
{
    PhaseTimer pt;

    size_t nodeCount    = domain.getNodeCount();
    size_t elementCount = domain.getElementCount();

    // domain getters trigger lazy halo + topology + adjacency + coord caching.
    // Time them as one bucket: "lazy domain prep". Halo+topology was already
    // done at amr.initialize on level 0; this captures only the post-AMR-rebuild
    // re-trigger + coord cache.
    const auto& d_nodeOwnership = domain.getNodeOwnershipMap();
    const auto& d_conn          = domain.getElementToNodeConnectivity();
    domain.cacheNodeCoordinates();
    const auto& d_x = domain.getNodeX();
    const auto& d_y = domain.getNodeY();
    const auto& d_z = domain.getNodeZ();
    pt.lap("lazy domain prep");

    // 1) GPU DOF mapping: owned nodes -> [0, numOwnedDofs), ghost -> [numOwnedDofs, nodeCount)
    cstone::DeviceVector<int> d_node_to_dof(nodeCount);
    int numOwnedDofs   = buildDofMappingGpu<KeyType>(d_nodeOwnership.data(), d_node_to_dof.data(), nodeCount);
    int numTotalDofs   = static_cast<int>(nodeCount);
    pt.lap("DOF mapping");

    // 2) GPU sparsity (FULL 27-NNZ pattern matching assembleFull's 8x8 stencil).
    //    Using graph (7-NNZ) pattern would silently drop entries in addValue,
    //    producing a corrupt matrix.
    cstone::DeviceVector<int> d_rowPtr(numTotalDofs + 1);
    cstone::DeviceVector<int> d_colInd;
    cstone::DeviceVector<int> d_diagPtr(numTotalDofs);

    int nnz = CvfemSparsityBuilder<KeyType>::buildFullSparsity(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(), std::get<2>(d_conn).data(),
        std::get<3>(d_conn).data(), std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
        std::get<6>(d_conn).data(), std::get<7>(d_conn).data(), elementCount, d_node_to_dof.data(), numTotalDofs,
        d_rowPtr.data(), nullptr, nullptr, 0);

    d_colInd.resize(nnz);
    CvfemSparsityBuilder<KeyType>::buildFullSparsity(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(), std::get<2>(d_conn).data(),
        std::get<3>(d_conn).data(), std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
        std::get<6>(d_conn).data(), std::get<7>(d_conn).data(), elementCount, d_node_to_dof.data(), numTotalDofs,
        d_rowPtr.data(), d_colInd.data(), d_diagPtr.data(), 0);

    cstone::DeviceVector<RealType> d_values(nnz, RealType(0));
    cstone::DeviceVector<RealType> d_rhs(numTotalDofs, RealType(0));
    pt.lap("sparsity build");

    // CSR wrapper for the assembler
    using MatrixType = CSRMatrix<RealType>;
    MatrixType* d_matrix;
    cudaMalloc(&d_matrix, sizeof(MatrixType));
    MatrixType h_matrix{d_rowPtr.data(), d_colInd.data(), d_values.data(), d_diagPtr.data(), numTotalDofs, nnz};
    cudaMemcpy(d_matrix, &h_matrix, sizeof(MatrixType), cudaMemcpyHostToDevice);

    // 3) GPU field setup + assembly. For Poisson: gamma=1, phi=0, beta=0
    cstone::DeviceVector<RealType> d_gamma(nodeCount, RealType(1));
    cstone::DeviceVector<RealType> d_phi(nodeCount, RealType(0));
    cstone::DeviceVector<RealType> d_beta(nodeCount, RealType(0));
    cstone::DeviceVector<RealType> d_grad_phi_x(nodeCount, RealType(0));
    cstone::DeviceVector<RealType> d_grad_phi_y(nodeCount, RealType(0));
    cstone::DeviceVector<RealType> d_grad_phi_z(nodeCount, RealType(0));

    cstone::DeviceVector<RealType> d_areaVec_x(elementCount * 12);
    cstone::DeviceVector<RealType> d_areaVec_y(elementCount * 12);
    cstone::DeviceVector<RealType> d_areaVec_z(elementCount * 12);
    precomputeAreaVectorsGpu<KeyType, RealType>(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(), std::get<2>(d_conn).data(),
        std::get<3>(d_conn).data(), std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
        std::get<6>(d_conn).data(), std::get<7>(d_conn).data(), elementCount, d_x.data(), d_y.data(), d_z.data(),
        d_areaVec_x.data(), d_areaVec_y.data(), d_areaVec_z.data());

    cstone::DeviceVector<RealType> d_mdot(elementCount * 12, RealType(0));

    typename CvfemHexAssembler<KeyType, RealType>::Config config;
    config.blockSize = blockSize;
    config.variant   = kernelVariant;

    CvfemHexAssembler<KeyType, RealType>::assembleFull(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(), std::get<2>(d_conn).data(),
        std::get<3>(d_conn).data(), std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
        std::get<6>(d_conn).data(), std::get<7>(d_conn).data(), elementCount, d_x.data(), d_y.data(), d_z.data(),
        d_gamma.data(), d_phi.data(), d_beta.data(), d_grad_phi_x.data(), d_grad_phi_y.data(), d_grad_phi_z.data(),
        d_mdot.data(), d_areaVec_x.data(), d_areaVec_y.data(), d_areaVec_z.data(),
        d_node_to_dof.data(), d_nodeOwnership.data(), d_matrix, d_rhs.data(), config);
    cudaDeviceSynchronize();
    pt.lap("assembly (LHS+RHS)");

    // Pressure-Poisson RHS: integrate f = -div(v . grad v) at element centroids
    // (one-point quadrature, exact for piecewise-constant f), scatter
    // V_elem * f / 8 to each owned corner. Lumped weighting matches
    // the boundary CG solve's expected residual form.
    {
        int eBlocks = (elementCount + blockSize - 1) / blockSize;
        integratePressurePoissonRhsKernel<KeyType, RealType><<<eBlocks, blockSize>>>(
            std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
            std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
            std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
            std::get<6>(d_conn).data(), std::get<7>(d_conn).data(),
            d_x.data(), d_y.data(), d_z.data(),
            d_node_to_dof.data(), d_nodeOwnership.data(),
            d_rhs.data(), elementCount, xcVortex, ycVortex, Uinf, gamma, eps_vortex);
        cudaDeviceSynchronize();
    }
    pt.lap("pressure-Poisson RHS integration");

    // 4) GPU bounding box (per-rank) + MPI_Allreduce for global box
    auto xb = thrust::device_pointer_cast(d_x.data());
    auto yb = thrust::device_pointer_cast(d_y.data());
    auto zb = thrust::device_pointer_cast(d_z.data());
    RealType lxmin = thrust::reduce(thrust::device, xb, xb + nodeCount, RealType(1e30),  thrust::minimum<RealType>());
    RealType lxmax = thrust::reduce(thrust::device, xb, xb + nodeCount, RealType(-1e30), thrust::maximum<RealType>());
    RealType lymin = thrust::reduce(thrust::device, yb, yb + nodeCount, RealType(1e30),  thrust::minimum<RealType>());
    RealType lymax = thrust::reduce(thrust::device, yb, yb + nodeCount, RealType(-1e30), thrust::maximum<RealType>());
    RealType lzmin = thrust::reduce(thrust::device, zb, zb + nodeCount, RealType(1e30),  thrust::minimum<RealType>());
    RealType lzmax = thrust::reduce(thrust::device, zb, zb + nodeCount, RealType(-1e30), thrust::maximum<RealType>());

    auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    RealType xmin, xmax, ymin, ymax, zmin, zmax;
    MPI_Allreduce(&lxmin, &xmin, 1, mpiType, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&lxmax, &xmax, 1, mpiType, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&lymin, &ymin, 1, mpiType, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&lymax, &ymax, 1, mpiType, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&lzmin, &zmin, 1, mpiType, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&lzmax, &zmax, 1, mpiType, MPI_MAX, MPI_COMM_WORLD);
    RealType eps = RealType(1e-10) * std::max({xmax - xmin, ymax - ymin, zmax - zmin});

    // 5) GPU BC: mark boundary DOFs, then enforce row=0/diag=1/rhs=0
    cstone::DeviceVector<uint8_t> d_isBdryDof(numOwnedDofs, 0);
    int nodeBlocks = (nodeCount + blockSize - 1) / blockSize;
    applyBCKernel<<<nodeBlocks, blockSize>>>(d_x.data(), d_y.data(), d_z.data(), d_nodeOwnership.data(),
                                              d_node_to_dof.data(), d_isBdryDof.data(), nodeCount,
                                              xmin, xmax, ymin, ymax, zmin, zmax, eps);

    int dofBlocks = (numOwnedDofs + blockSize - 1) / blockSize;
    enforceBCKernel<<<dofBlocks, blockSize>>>(d_isBdryDof.data(), d_rowPtr.data(), d_colInd.data(), d_diagPtr.data(),
                                               d_values.data(), d_rhs.data(), numOwnedDofs);
    cudaDeviceSynchronize();
    pt.lap("BC enforcement");

    // 6) CG solve. SparseMatrix is m x n: m=numOwnedDofs rows (one per owned DOF),
    // n=numTotalDofs columns (entries reference both owned and ghost DOFs).
    using Matrix = SparseMatrix<int, RealType, cstone::GpuTag>;
    Matrix A;
    A.allocate(numOwnedDofs, numTotalDofs, nnz);
    cudaMemcpy(A.rowOffsetsPtr(),  d_rowPtr.data(),  (numOwnedDofs + 1) * sizeof(int),     cudaMemcpyDeviceToDevice);
    cudaMemcpy(A.colIndicesPtr(),  d_colInd.data(),  nnz * sizeof(int),                    cudaMemcpyDeviceToDevice);
    cudaMemcpy(A.valuesPtr(),      d_values.data(),  nnz * sizeof(RealType),               cudaMemcpyDeviceToDevice);

    // CG solve. b has size m (owned RHS), x is auto-resized to n by solver.
    cstone::DeviceVector<RealType> b(numOwnedDofs), x(numTotalDofs);
    cudaMemcpy(b.data(), d_rhs.data(), numOwnedDofs * sizeof(RealType), cudaMemcpyDeviceToDevice);

    // Initial guess for x: if d_nodeSolution arrives non-empty (from the AMR
    // transfer step at the previous level), scatter its node-keyed values into
    // x's DOF layout. Otherwise zero. Boundary DOFs get overwritten by BC
    // enforcement before solve.
    cudaMemset(x.data(), 0, numTotalDofs * sizeof(RealType));
    if (d_nodeSolution.size() == nodeCount)
    {
        thrust::for_each(thrust::device,
                          thrust::counting_iterator<size_t>(0),
                          thrust::counting_iterator<size_t>(nodeCount),
                          [nodeToDof = d_node_to_dof.data(),
                           nodeSol   = d_nodeSolution.data(),
                           xOut      = x.data()] __device__(size_t i) {
                              int dof = nodeToDof[i];
                              if (dof >= 0) xOut[dof] = nodeSol[i];
                          });
    }

    int nRanks = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &nRanks);

    if (solverKind == SolverKind::CG)
    {
        ConjugateGradientSolver<RealType, int, cstone::GpuTag> solver(maxIter, tolerance);
        solver.setVerbose(false);
        solver.setOwnedSize(numOwnedDofs); // enables MPI_Allreduce in dot products

        // Distributed CG halo exchange: delegate to ElementDomain's GPU-native
        // node-level halo built on top of cstone's per-element halo.
        if (nRanks > 1)
        {
            const int* dofMapPtr = d_node_to_dof.data();
            solver.setHaloExchangeCallback(
                [&domain, dofMapPtr](cstone::DeviceVector<RealType>& p)
                {
                    domain.exchangeNodeHalo(p, dofMapPtr);
                });
        }

        pt.lap("CG setup (Matrix alloc + halo cb)");

        solver.solve(A, b, x);
        cudaDeviceSynchronize();
        pt.lap("CG solve");
    }
    else // SolverKind::Hypre
    {
#ifdef MARS_ENABLE_HYPRE
        // Build a contiguous global owned-DOF numbering for the interior system.
        // Each rank owns [globalRowStart, globalRowEnd); ghost local columns
        // [numOwnedDofs, numTotalDofs) are remapped to their owner-rank global IDs
        // via a halo exchange of the per-node global DOF index.
        using IdxType = int64_t;
        IdxType numOwnedLocal = static_cast<IdxType>(numOwnedDofs);
        IdxType globalRowStart = 0;
        MPI_Exscan(&numOwnedLocal, &globalRowStart, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
        if (rank == 0) globalRowStart = 0;
        IdxType globalRowEnd     = globalRowStart + numOwnedLocal;
        IdxType numInteriorGlobal = 0;
        MPI_Allreduce(&numOwnedLocal, &numInteriorGlobal, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);

        // Per-node global DOF ID: owned nodes get their global id, ghosts start
        // at -1 and are filled by the halo exchange from the owner rank. Use
        // RealType to ride the existing exchangeNodeHalo path; double exactly
        // represents integers up to 2^53 (>> any Gordon-Bell target).
        cstone::DeviceVector<RealType> d_nodeGlobalDof(nodeCount, RealType(-1));
        {
            const int* nodeToDofPtr = d_node_to_dof.data();
            const uint8_t* ownPtr   = d_nodeOwnership.data();
            RealType rs             = static_cast<RealType>(globalRowStart);
            thrust::for_each(thrust::device,
                              thrust::counting_iterator<size_t>(0),
                              thrust::counting_iterator<size_t>(nodeCount),
                              [nodeToDofPtr, ownPtr, rs,
                               out = d_nodeGlobalDof.data()] __device__(size_t i)
                              {
                                  if (ownPtr[i] == 1)
                                  {
                                      int dof = nodeToDofPtr[i];
                                      if (dof >= 0) out[i] = rs + static_cast<RealType>(dof);
                                  }
                              });
            domain.exchangeNodeHalo(d_nodeGlobalDof);
        }

        // Local-to-global DOF map indexed by local DOF id in [0, numTotalDofs).
        // Walk the per-node global DOF on host (one bulk download, no hot path).
        std::vector<RealType> h_nodeGlobalDof(nodeCount);
        thrust::copy(thrust::device_pointer_cast(d_nodeGlobalDof.data()),
                     thrust::device_pointer_cast(d_nodeGlobalDof.data() + nodeCount),
                     h_nodeGlobalDof.begin());
        std::vector<int> h_nodeToDof(nodeCount);
        thrust::copy(thrust::device_pointer_cast(d_node_to_dof.data()),
                     thrust::device_pointer_cast(d_node_to_dof.data() + nodeCount),
                     h_nodeToDof.begin());

        std::vector<IdxType> localToGlobalDof(numTotalDofs, IdxType(-1));
        for (size_t i = 0; i < nodeCount; ++i)
        {
            int dof = h_nodeToDof[i];
            if (dof >= 0)
            {
                localToGlobalDof[dof] = static_cast<IdxType>(h_nodeGlobalDof[i]);
            }
        }

        pt.lap("Hypre setup (global DOF map + halo)");

        // The wrapper expects A of size m=numOwnedDofs rows, n=numTotalDofs cols
        // (already what we built). It internally re-uses comm_=MPI_COMM_WORLD.
        // SparseMatrix typedefs in the wrapper are <int, RealType, GpuTag>; pass
        // matching index type. We already use int rowPtr/colInd; cast IdxType
        // localToGlobalDof to the wrapper's KeyType template argument.
        mars::fem::HyprePCGSolver<RealType, int, cstone::GpuTag>
            hypreSolver(MPI_COMM_WORLD, maxIter, tolerance,
                        mars::fem::HyprePCGSolver<RealType, int, cstone::GpuTag>::BOOMERAMG);
        hypreSolver.setVerbose(rank == 0);

        bool converged = hypreSolver.template solve<IdxType>(
            A, b, x,
            static_cast<int>(globalRowStart), static_cast<int>(globalRowEnd),
            0, static_cast<int>(numInteriorGlobal),
            localToGlobalDof);
        (void)converged;
        cudaDeviceSynchronize();
        pt.lap("Hypre PCG+BoomerAMG solve");
#else
        if (rank == 0)
            std::cerr << "ERROR: --solver=hypre requested but MARS was built without MARS_ENABLE_HYPRE\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    }

    // DOF -> per-node solution
    d_nodeSolution.resize(nodeCount);
    cudaMemset(d_nodeSolution.data(), 0, nodeCount * sizeof(RealType));
    thrust::for_each(thrust::device,
                      thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(nodeCount),
                      [nodeToDof = d_node_to_dof.data(), sol = x.data(),
                       nodeSol = d_nodeSolution.data()] __device__(size_t i)
                      {
                          int dof = nodeToDof[i];
                          if (dof >= 0) nodeSol[i] = sol[dof];
                      });

    // Halo-exchange d_nodeSolution so the error indicator (per-element gradient)
    // sees correct ghost-node values for boundary elements.
    domain.exchangeNodeHalo(d_nodeSolution);

    // Error indicator (per-element)
    d_errorPerElement = HexErrorIndicator<KeyType, RealType>::computeError(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(), std::get<2>(d_conn).data(),
        std::get<3>(d_conn).data(), std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
        std::get<6>(d_conn).data(), std::get<7>(d_conn).data(), d_nodeSolution.data(), d_x.data(), d_y.data(),
        d_z.data(), d_node_to_dof.data(), elementCount, blockSize);
    pt.lap("error indicator");

    cudaFree(d_matrix);

    // Collective: must be called from all ranks. Sum only OWNED elements
    // (range [startIdx, endIdx)) so halo elements aren't double-counted
    // across ranks - they appear in multiple ranks' element arrays.
    auto& cstoneDom = domain.getDomain();
    RealType errNorm = ErrorIndicator<KeyType, RealType>::globalErrorNormOwned(
        d_errorPerElement, cstoneDom.startIndex(), cstoneDom.endIndex());

    // Global ||u||_2: each rank sums squares over its OWNED DOFs (no double-count
    // since owned DOFs are unique across ranks), then MPI_Allreduce(SUM).
    auto x_ptr            = thrust::device_pointer_cast(x.data());
    RealType localSumSq   = thrust::inner_product(thrust::device, x_ptr, x_ptr + numOwnedDofs, x_ptr, RealType(0));
    RealType globalSumSq  = 0;
    auto mpiTypeR         = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    MPI_Allreduce(&localSumSq, &globalSumSq, 1, mpiTypeR, MPI_SUM, MPI_COMM_WORLD);
    RealType solNorm = std::sqrt(globalSumSq);

    // Global owned-DOF count (sum of per-rank numOwnedDofs)
    int globalDofs = 0;
    MPI_Allreduce(&numOwnedDofs, &globalDofs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0)
    {
        std::cout << "    Solve: " << globalDofs << " global DOFs (rank0=" << numOwnedDofs
                  << "), " << nnz << " NNZ, ||u||=" << std::scientific
                  << solNorm << ", ||err||=" << errNorm << "\n"
                  << std::defaultfloat;
    }

    // Phase breakdown for this level's solve (everything inside solvePressurePoisson)
    pt.report(rank, "solve phases");
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, numRanks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount > 0) cudaSetDevice(rank % deviceCount);

    // Parse arguments
    std::string meshFile;
    CvfemKernelVariant kernelVariant = CvfemKernelVariant::Tensor;
    int blockSize                    = 256;
    int bucketSize                   = 64;
    double Uinf                      = 1.0;     // free-stream speed
    double gammaVortex               = 0.5;     // 2D potential vortex strength
    double epsVortex                 = 1e-6;    // vortex-core regularization
    // Vortex line is at (xcVortex, ycVortex, *). Defaults place it at the
    // cube16 mesh center so the singular spike is in the interior, not at a
    // corner where the colormap gets dominated by one outlier cell. Override
    // with --vortex-center=x,y,z for wing-tip (e.g. 0.9,0.44,0).
    double xcVortex                  = 0.5;
    double ycVortex                  = 0.5;
    double zcVortex                  = 0.0;     // currently unused
    int maxIter                      = 1000;
    double tolerance                 = 1e-10;
    int amrLevels                    = 3;
    double refineFraction            = 0.3;
    double coarsenFraction           = 0.03;
    std::string vtuPrefix;  // empty = no output
    int gammaFrames                  = 0;       // 0 = no sweep (single AMR-level movie only)
    double gammaSweepMin             = 0.5;     // start of gamma sweep
    double gammaSweepMax             = 10.0;    // end of gamma sweep
    SolverKind solverKind            = SolverKind::CG;   // default: in-tree CG (unchanged)

    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg.find("--mesh=") == 0)
            meshFile = arg.substr(7);
        else if (arg.find("--kernel=") == 0)
        {
            std::string v = arg.substr(9);
            if (v == "tensor") kernelVariant = CvfemKernelVariant::Tensor;
            else if (v == "original") kernelVariant = CvfemKernelVariant::Original;
            else if (v == "optimized") kernelVariant = CvfemKernelVariant::Optimized;
            else if (v == "shmem") kernelVariant = CvfemKernelVariant::Shmem;
            else if (v == "tensor_colored") kernelVariant = CvfemKernelVariant::TensorColored;
            else if (v == "tensor_aos") kernelVariant = CvfemKernelVariant::TensorAoS;
            else if (v == "wmma_tensor") kernelVariant = CvfemKernelVariant::WmmaTensor;
        }
        else if (arg.find("--block-size=") == 0)
            blockSize = std::stoi(arg.substr(13));
        else if (arg.find("--bucket-size=") == 0)
            bucketSize = std::stoi(arg.substr(14));
        else if (arg.find("--uinf=") == 0)
            Uinf = std::stod(arg.substr(7));
        else if (arg.find("--gamma=") == 0)
            gammaVortex = std::stod(arg.substr(8));
        else if (arg.find("--eps-vortex=") == 0)
            epsVortex = std::stod(arg.substr(13));
        else if (arg.find("--vortex-center=") == 0)
        {
            std::string v = arg.substr(16);
            size_t c1 = v.find(','), c2 = (c1 == std::string::npos ? std::string::npos : v.find(',', c1 + 1));
            if (c1 == std::string::npos || c2 == std::string::npos)
            {
                if (rank == 0)
                    std::cerr << "Error: --vortex-center=x,y,z requires 3 comma-separated values\n";
                MPI_Finalize();
                return 1;
            }
            xcVortex = std::stod(v.substr(0, c1));
            ycVortex = std::stod(v.substr(c1 + 1, c2 - c1 - 1));
            zcVortex = std::stod(v.substr(c2 + 1));
        }
        else if (arg.find("--max-iter=") == 0)
            maxIter = std::stoi(arg.substr(11));
        else if (arg.find("--tol=") == 0)
            tolerance = std::stod(arg.substr(6));
        else if (arg.find("--amr-levels=") == 0)
            amrLevels = std::stoi(arg.substr(13));
        else if (arg.find("--refine-frac=") == 0)
            refineFraction = std::stod(arg.substr(14));
        else if (arg.find("--coarsen-frac=") == 0)
            coarsenFraction = std::stod(arg.substr(15));
        else if (arg.find("--vtu-output=") == 0)
            vtuPrefix = arg.substr(13);
        else if (arg.find("--gamma-frames=") == 0)
            gammaFrames = std::stoi(arg.substr(15));
        else if (arg.find("--gamma-min=") == 0)
            gammaSweepMin = std::stod(arg.substr(12));
        else if (arg.find("--gamma-max=") == 0)
            gammaSweepMax = std::stod(arg.substr(12));
        else if (arg.find("--solver=") == 0)
        {
            std::string v = arg.substr(9);
            if (v == "cg") solverKind = SolverKind::CG;
            else if (v == "hypre") solverKind = SolverKind::Hypre;
            else
            {
                if (rank == 0) std::cerr << "Error: --solver must be cg or hypre, got '" << v << "'\n";
                MPI_Finalize();
                return 1;
            }
        }
        else if (arg[0] != '-' && meshFile.empty())
            meshFile = arg;
    }

    if (meshFile.empty())
    {
        if (rank == 0)
        {
            std::cout << "Usage: " << argv[0] << " [options]\n\n";
            std::cout << "Options:\n";
            std::cout << "  --mesh=FILE          Mesh file (.mesh or .exo) [REQUIRED]\n";
            std::cout << "  --kernel=VARIANT     tensor, original, optimized, shmem, wmma_tensor\n";
            std::cout << "  --uinf=VALUE         Free-stream speed Uinf (default: 1.0)\n";
            std::cout << "  --gamma=VALUE        Tip-vortex strength gamma (default: 0.5)\n";
            std::cout << "  --eps-vortex=VALUE   Vortex-core regularization eps (default: 1e-6)\n";
            std::cout << "  --vortex-center=X,Y,Z  Vortex line center (default: 0.5,0.5,0 for cube)\n";
            std::cout << "  --max-iter=N         CG max iterations (default: 1000)\n";
            std::cout << "  --tol=VALUE          CG tolerance (default: 1e-10)\n";
            std::cout << "  --amr-levels=N       Max AMR levels (default: 3)\n";
            std::cout << "  --refine-frac=VALUE  Refine fraction (default: 0.3)\n";
            std::cout << "  --coarsen-frac=VALUE Coarsen fraction (default: 0.03)\n";
            std::cout << "  --bucket-size=N      Cornerstone bucket size (default: 64)\n";
            std::cout << "  --block-size=N       CUDA block size (default: 256)\n";
            std::cout << "  --vtu-output=PREFIX  Write per-AMR-level VTU/PVTU/PVD frames (default: off)\n";
            std::cout << "  --gamma-frames=N     After AMR, sweep gamma in N frames for time-animation\n";
            std::cout << "  --gamma-min=VALUE    Sweep start (default: 0.5)\n";
            std::cout << "  --gamma-max=VALUE    Sweep end (default: 10.0)\n";
            std::cout << "  --solver=KIND        cg | hypre (default: cg; hypre requires MARS_ENABLE_HYPRE=ON)\n";
        }
        MPI_Finalize();
        return 1;
    }

    using KeyType  = uint64_t;
    using RealType = double;

    if (rank == 0)
    {
        std::cout << "\n========================================\n";
        std::cout << "MARS AMR Pressure-Poisson Solver\n";
        std::cout << "========================================\n";
        std::cout << "Problem: -Laplacian(p) = -div(v . grad v), p = 0 on boundary\n";
        std::cout << "Velocity: v = (Uinf, 0, 0) + gamma*(-(y-yc), (x-xc), 0) / ((x-xc)^2 + (y-yc)^2 + eps)\n";
        std::cout << "  Uinf   = " << Uinf << "\n";
        std::cout << "  gamma  = " << gammaVortex << "\n";
        std::cout << "  eps    = " << epsVortex << "\n";
        std::cout << "  center = (" << xcVortex << ", " << ycVortex << ", " << zcVortex << ")\n";
        std::cout << "Mesh: " << meshFile << "\n";
        std::cout << "Kernel: " << CvfemHexAssembler<KeyType, RealType>::variantName(kernelVariant) << "\n";
        std::cout << "Solver: " << (solverKind == SolverKind::CG ? "cg" : "hypre (PCG+BoomerAMG)") << "\n";
        std::cout << "MPI ranks: " << numRanks << "\n";
        std::cout << "AMR levels: " << amrLevels << "\n";
        std::cout << "Bucket size: " << bucketSize << "\n";
        std::cout << "========================================\n\n";
    }

    // Configure and initialize AMR
    AmrManager<HexTag, KeyType, RealType>::Config amrConfig;
    amrConfig.maxLevels       = amrLevels;
    amrConfig.refineFraction  = refineFraction;
    amrConfig.coarsenFraction = coarsenFraction;
    amrConfig.blockSize       = blockSize;
    amrConfig.bucketSize      = bucketSize;

    AmrManager<HexTag, KeyType, RealType> amr(amrConfig);
    amr.initialize(meshFile, rank, numRanks);
    auto initT = amr.initTimings();

    if (rank == 0)
    {
        std::cout << "Initial mesh:\n";
        std::cout << "  Elements: " << amr.domain().getElementCount() << "\n";
        std::cout << "  Nodes:    " << amr.domain().getNodeCount() << "\n";
        std::cout << "  Init breakdown (ms):\n";
        std::cout << "    domain (file read + sync):    " << std::fixed << initT.domainSyncTimeMs << "\n";
        std::cout << "    halo + node topology:         " << initT.haloTopoTimeMs   << "\n";
        std::cout << "    adjacency CSR (e2n + n2e):    " << initT.adjacencyTimeMs  << "\n";
        std::cout << "    coord cache (SFC decode):     " << initT.coordCacheTimeMs << "\n";
        std::cout << "    AMR octree state:             " << initT.octreeTimeMs     << "\n";
        std::cout << "    TOTAL:                        " << initT.totalMs          << "\n\n";
    }

    // AMR loop
    cstone::DeviceVector<RealType> d_nodeSolution;
    cstone::DeviceVector<RealType> d_errorPerElement;
    AmrStats stats;
    stats.elementsRefined = 1;

    // Optional VTU output: one PVTU per AMR level + a PVD time-collection that
    // ParaView opens as an animation.
    std::unique_ptr<fem::VTUParallelWriter<KeyType, RealType>> vtuWriter;
    if (!vtuPrefix.empty())
    {
        vtuWriter = std::make_unique<fem::VTUParallelWriter<KeyType, RealType>>(vtuPrefix);
        if (rank == 0)
            std::cout << "VTU output enabled: " << vtuPrefix << "_step*.pvtu, "
                      << vtuPrefix << ".pvd\n\n";
    }

    auto totalStart = std::chrono::high_resolution_clock::now();

    for (int level = 0; level <= amrLevels; ++level)
    {
        auto levelStart = std::chrono::high_resolution_clock::now();

        if (rank == 0)
        {
            std::cout << "--- AMR Level " << level << " ---\n";
            std::cout << "  Elements: " << amr.domain().getElementCount()
                      << " (local: " << amr.domain().localElementCount() << ")\n";
            std::cout << "  Nodes:    " << amr.domain().getNodeCount() << "\n";
        }

        // Solve on current mesh (graph-based assembly)
        std::cerr << "MAIN r" << rank << " before solve L" << level << std::endl; std::cerr.flush();
        solvePressurePoisson<KeyType, RealType>(amr.domain(), xcVortex, ycVortex, Uinf,
                                                 gammaVortex, epsVortex,
                                                 kernelVariant, blockSize, maxIter, tolerance,
                                                 solverKind,
                                                 rank, d_nodeSolution, d_errorPerElement);
        cudaDeviceSynchronize();
        MPI_Barrier(MPI_COMM_WORLD);
        std::cerr << "MAIN r" << rank << " solve L" << level << " done (post-barrier)" << std::endl; std::cerr.flush();

        // VTU frame after each solved level. d_nodeSolution at this point
        // already had exchangeNodeHalo() called inside solvePressurePoisson, so
        // ghost slots hold correct owner-rank values.
        if (vtuWriter)
        {
            vtuWriter->writeFrame(level, static_cast<double>(level), amr.domain(), d_nodeSolution, "u");
            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == 0)
                std::cout << "  VTU: wrote " << vtuPrefix << "_step"
                          << std::setw(4) << std::setfill('0') << level << ".pvtu (+ rank pieces)\n";
        }

        auto levelEnd  = std::chrono::high_resolution_clock::now();
        float levelMs = std::chrono::duration<float, std::milli>(levelEnd - levelStart).count();
        if (rank == 0) std::cout << "    Level time: " << std::fixed << levelMs << " ms\n";

        if (level >= amrLevels) break;

        // Adapt mesh (mark + refine + rebuild + transfer). adaptMesh fills
        // stats.totalTimeMs which is the total wall-clock for the whole AMR
        // step (mark, refine, cstone-rebuild, transfer, halo exchange).
        cstone::DeviceVector<RealType> d_newSolution;
        auto amrStart = std::chrono::high_resolution_clock::now();
        stats = amr.adaptMesh(d_errorPerElement.data(), d_nodeSolution.data(), d_newSolution);
        cudaDeviceSynchronize();
        MPI_Barrier(MPI_COMM_WORLD);
        auto amrEnd = std::chrono::high_resolution_clock::now();
        float amrMs = std::chrono::duration<float, std::milli>(amrEnd - amrStart).count();
        AmrManager<HexTag, KeyType, RealType>::printStats(stats, rank);
        if (rank == 0) std::cout << "    AMR wall: " << std::fixed << amrMs << " ms\n";

        d_nodeSolution = std::move(d_newSolution);

        if (!amr.shouldContinue(stats))
        {
            if (rank == 0) std::cout << "\n  AMR converged.\n";
            break;
        }

        if (rank == 0) std::cout << "\n";
    }

    // Gamma sweep: re-solve on the final AMR mesh with varying vortex strength.
    // Gives a real time-animation movie (pressure builds as vortex spins up) on
    // top of the AMR snapshots. Step numbers continue from the AMR frames so
    // the PVD plays the whole thing as one continuous animation.
    if (gammaFrames > 0 && vtuWriter)
    {
        int baseStep = amr.currentLevel() + 1;
        if (rank == 0)
        {
            std::cout << "\n========================================\n";
            std::cout << "Gamma sweep: " << gammaFrames << " frames, gamma "
                      << gammaSweepMin << " -> " << gammaSweepMax << "\n";
            std::cout << "========================================\n";
        }

        for (int g = 0; g < gammaFrames; ++g)
        {
            double t       = (gammaFrames == 1) ? 0.0 : double(g) / double(gammaFrames - 1);
            double gammaG  = gammaSweepMin + t * (gammaSweepMax - gammaSweepMin);
            int step       = baseStep + g;

            if (rank == 0)
                std::cout << "--- Sweep frame " << g << "/" << gammaFrames
                          << " (gamma = " << gammaG << ") ---\n";

            solvePressurePoisson<KeyType, RealType>(amr.domain(), xcVortex, ycVortex, Uinf,
                                                     RealType(gammaG), epsVortex,
                                                     kernelVariant, blockSize, maxIter, tolerance,
                                                     solverKind,
                                                     rank, d_nodeSolution, d_errorPerElement);
            cudaDeviceSynchronize();
            MPI_Barrier(MPI_COMM_WORLD);

            // Synthetic time coordinate: AMR frames are 0..baseStep-1 in
            // integer time, sweep frames continue from baseStep. ParaView
            // animates strictly on this axis.
            vtuWriter->writeFrame(step, double(step), amr.domain(), d_nodeSolution, "u");
            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == 0)
                std::cout << "  VTU: wrote " << vtuPrefix << "_step"
                          << std::setw(4) << std::setfill('0') << step << ".pvtu\n";
        }
    }

    auto totalEnd = std::chrono::high_resolution_clock::now();
    float totalMs = std::chrono::duration<float, std::milli>(totalEnd - totalStart).count();

    // Final statistics
    if (rank == 0)
    {
        size_t nNodes = d_nodeSolution.size();
        auto sol_ptr  = thrust::device_pointer_cast(d_nodeSolution.data());
        RealType solMin = thrust::reduce(thrust::device, sol_ptr, sol_ptr + nNodes,
                                          std::numeric_limits<RealType>::max(), thrust::minimum<RealType>());
        RealType solMax = thrust::reduce(thrust::device, sol_ptr, sol_ptr + nNodes,
                                          std::numeric_limits<RealType>::lowest(), thrust::maximum<RealType>());
        RealType solNorm = std::sqrt(
            thrust::inner_product(thrust::device, sol_ptr, sol_ptr + nNodes, sol_ptr, RealType(0)));

        std::cout << "\n========================================\n";
        std::cout << "Final Solution\n";
        std::cout << "========================================\n";
        std::cout << "  Min:     " << std::scientific << solMin << "\n";
        std::cout << "  Max:     " << solMax << "\n";
        std::cout << "  L2 norm: " << solNorm << "\n";
        std::cout << "  AMR levels: " << amr.currentLevel() << "\n";
        std::cout << "  Total time: " << std::fixed << totalMs << " ms\n";
        std::cout << "========================================\n";
    }

    MPI_Finalize();
    return 0;
}
