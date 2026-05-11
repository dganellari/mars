// AMR CVFEM Graph Assembly: adaptive mesh refinement with graph-based CVFEM.
// Based on mars_cvfem_graph.cu (the primary CVFEM driver) + AMR loop.
// All operations GPU-native, using cstone octree for domain decomposition.

#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_fem.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_assembler.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_utils.hpp"
#include "backend/distributed/unstructured/fem/mars_sparsity_builder.hpp"
#include "backend/distributed/unstructured/fem/mars_perf_counters.hpp"
#include "backend/distributed/unstructured/fem/mars_sparse_matrix.hpp"
#include "backend/distributed/unstructured/solvers/mars_cg_solver.hpp"
#include "backend/distributed/unstructured/amr/mars_amr.hpp"

#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/inner_product.h>
#include <thrust/transform.h>
#include <thrust/system/cuda/execution_policy.h>
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

// Integrate constant source term against trilinear hex8 shape functions.
// Each element with volume V contributes sourceTerm * V / 8 to each of its
// 8 corner-node owned DOFs. Adds with sign matching the system Aφ = b
// produced by assembleFull (residual-form): b_i = sourceTerm * V_i.
template<typename KeyType, typename RealType>
__global__ void integrateConstantSourceKernel(const KeyType* conn0,
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
                                              RealType sourceTerm)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements) return;

    KeyType n[8] = {conn0[e], conn1[e], conn2[e], conn3[e],
                    conn4[e], conn5[e], conn6[e], conn7[e]};

    // Hex volume via box approximation: |x_max-x_min| * |y_max-y_min| * |z_max-z_min|
    // For axis-aligned regular hex this is exact; for distorted hex it's a fast estimate
    RealType x[8], y[8], z[8];
    for (int i = 0; i < 8; ++i)
    {
        x[i] = nodeX[n[i]];
        y[i] = nodeY[n[i]];
        z[i] = nodeZ[n[i]];
    }
    RealType xmin = x[0], xmax = x[0], ymin = y[0], ymax = y[0], zmin = z[0], zmax = z[0];
    for (int i = 1; i < 8; ++i)
    {
        if (x[i] < xmin) xmin = x[i]; if (x[i] > xmax) xmax = x[i];
        if (y[i] < ymin) ymin = y[i]; if (y[i] > ymax) ymax = y[i];
        if (z[i] < zmin) zmin = z[i]; if (z[i] > zmax) zmax = z[i];
    }
    RealType V = (xmax - xmin) * (ymax - ymin) * (zmax - zmin);
    RealType contrib = sourceTerm * V * RealType(0.125); // / 8

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

// GPU-native Poisson solve on an AMR-managed domain. Pipeline:
//   1) GPU DOF mapping (buildDofMappingGpu)
//   2) GPU sparsity (CvfemSparsityBuilder::buildGraphSparsity, sized for owned DOFs)
//   3) GPU full 8x8 element assembly (assembleFull)
//   4) GPU bounding box reduction + MPI_Allreduce for global box
//   5) GPU BC marking + GPU BC enforcement (applyBCKernel + enforceBCKernel)
//   6) CG solve
//   7) GPU DOF -> per-node scatter
// No host downloads of node/element data.
template<typename KeyType, typename RealType>
void solvePoissonGraph(ElementDomain<HexTag, RealType, KeyType, cstone::GpuTag>& domain,
                       RealType sourceTerm,
                       CvfemKernelVariant kernelVariant,
                       int blockSize,
                       int maxIter,
                       RealType tolerance,
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

    // Source term: integrate sourceTerm against trilinear shape functions per element.
    // Each owned node's RHS gets contribution sourceTerm * V_elem / 8 from each
    // adjacent element. Lumped-mass weighting; trapezoidal rule over the hex.
    {
        int eBlocks = (elementCount + blockSize - 1) / blockSize;
        integrateConstantSourceKernel<KeyType, RealType><<<eBlocks, blockSize>>>(
            std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
            std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
            std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
            std::get<6>(d_conn).data(), std::get<7>(d_conn).data(),
            d_x.data(), d_y.data(), d_z.data(),
            d_node_to_dof.data(), d_nodeOwnership.data(),
            d_rhs.data(), elementCount, sourceTerm);
        cudaDeviceSynchronize();
    }
    pt.lap("source integration");

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

    ConjugateGradientSolver<RealType, int, cstone::GpuTag> solver(maxIter, tolerance);
    solver.setVerbose(false);
    solver.setOwnedSize(numOwnedDofs); // enables MPI_Allreduce in dot products

    // Distributed CG halo exchange: delegate to ElementDomain's GPU-native
    // node-level halo built on top of cstone's per-element halo.
    int nRanks = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &nRanks);
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

    // Phase breakdown for this level's solve (everything inside solvePoissonGraph)
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
    double sourceTerm                = 1.0;
    int maxIter                      = 1000;
    double tolerance                 = 1e-10;
    int amrLevels                    = 3;
    double refineFraction            = 0.3;
    double coarsenFraction           = 0.03;

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
        else if (arg.find("--source=") == 0)
            sourceTerm = std::stod(arg.substr(9));
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
            std::cout << "  --source=VALUE       Source term f (default: 1.0)\n";
            std::cout << "  --max-iter=N         CG max iterations (default: 1000)\n";
            std::cout << "  --tol=VALUE          CG tolerance (default: 1e-10)\n";
            std::cout << "  --amr-levels=N       Max AMR levels (default: 3)\n";
            std::cout << "  --refine-frac=VALUE  Refine fraction (default: 0.3)\n";
            std::cout << "  --coarsen-frac=VALUE Coarsen fraction (default: 0.03)\n";
            std::cout << "  --bucket-size=N      Cornerstone bucket size (default: 64)\n";
            std::cout << "  --block-size=N       CUDA block size (default: 256)\n";
        }
        MPI_Finalize();
        return 1;
    }

    using KeyType  = uint64_t;
    using RealType = double;

    if (rank == 0)
    {
        std::cout << "\n========================================\n";
        std::cout << "MARS AMR CVFEM Graph Solver\n";
        std::cout << "========================================\n";
        std::cout << "Problem: -du = " << sourceTerm << ", u = 0 on boundary\n";
        std::cout << "Mesh: " << meshFile << "\n";
        std::cout << "Kernel: " << CvfemHexAssembler<KeyType, RealType>::variantName(kernelVariant) << "\n";
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
        solvePoissonGraph<KeyType, RealType>(amr.domain(), sourceTerm, kernelVariant, blockSize, maxIter, tolerance,
                                              rank, d_nodeSolution, d_errorPerElement);
        cudaDeviceSynchronize();
        MPI_Barrier(MPI_COMM_WORLD);
        std::cerr << "MAIN r" << rank << " solve L" << level << " done (post-barrier)" << std::endl; std::cerr.flush();

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
