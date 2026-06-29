// Tet AMR driver: parallel of mars_amr_cvfem_graph.cu, but for linear
// tetrahedra. Same flow (init AMR, per-level Poisson solve with full 4x4
// element stencil, mark + refine + transfer + rebuild loop). Hex-specific
// kernel variants are dropped; tet has only the graph kernel for now.
#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_fem.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_tet_assembler.hpp"
#include "backend/distributed/unstructured/fem/mars_perf_counters.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_utils.hpp"
#include "backend/distributed/unstructured/fem/mars_sparse_matrix.hpp"
#include "backend/distributed/unstructured/fem/mars_sparsity_builder.hpp"
#include "backend/distributed/unstructured/solvers/mars_cg_solver.hpp"
#include "backend/distributed/unstructured/amr/mars_amr.hpp"

#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/inner_product.h>
#include <thrust/system/cuda/execution_policy.h>
#include <mpi.h>
#include <chrono>
#include <iomanip>
#include <limits>

using namespace mars;
using namespace mars::fem;
using namespace mars::amr;

// PhaseTimer: cudaDeviceSynchronize + MPI_Barrier on each lap so wall time
// reflects all-rank work, not just rank-0.
struct PhaseTimer
{
    using clk = std::chrono::high_resolution_clock;
    clk::time_point t0;
    std::vector<std::pair<std::string, float>> phases;

    PhaseTimer() { reset(); }
    void reset()
    {
        cudaDeviceSynchronize();
        MPI_Barrier(MPI_COMM_WORLD);
        t0 = clk::now();
    }
    void lap(const std::string& name)
    {
        cudaDeviceSynchronize();
        MPI_Barrier(MPI_COMM_WORLD);
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

// Mark owned DOFs whose node sits on the global bounding box boundary.
template<typename RealType>
__global__ void applyBCKernel(const RealType* nodeX,
                              const RealType* nodeY,
                              const RealType* nodeZ,
                              const uint8_t* ownership,
                              const int* nodeToDof,
                              uint8_t* isBoundaryDof,
                              size_t numNodes,
                              RealType xmin, RealType xmax,
                              RealType ymin, RealType ymax,
                              RealType zmin, RealType zmax,
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

// Integrate constant source term against linear tet shape functions.
// Tet volume V = |det J| / 6 where J is the affine map matrix; each of the
// 4 corner nodes gets V/4 (lumped) contribution.
template<typename KeyType, typename RealType>
__global__ void integrateConstantSourceTetKernel(const KeyType* conn0,
                                                  const KeyType* conn1,
                                                  const KeyType* conn2,
                                                  const KeyType* conn3,
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

    KeyType n[4] = {conn0[e], conn1[e], conn2[e], conn3[e]};
    RealType px[4], py[4], pz[4];
    for (int i = 0; i < 4; ++i)
    {
        px[i] = nodeX[n[i]];
        py[i] = nodeY[n[i]];
        pz[i] = nodeZ[n[i]];
    }

    RealType a0 = px[1] - px[0], a1 = py[1] - py[0], a2 = pz[1] - pz[0];
    RealType b0 = px[2] - px[0], b1 = py[2] - py[0], b2 = pz[2] - pz[0];
    RealType c0 = px[3] - px[0], c1 = py[3] - py[0], c2 = pz[3] - pz[0];
    RealType det = a0 * (b1 * c2 - b2 * c1)
                 - a1 * (b0 * c2 - b2 * c0)
                 + a2 * (b0 * c1 - b1 * c0);
    RealType V       = fabs(det) / RealType(6);
    RealType contrib = sourceTerm * V * RealType(0.25); // / 4

    for (int i = 0; i < 4; ++i)
    {
        if (ownership[n[i]] == 1)
        {
            int dof = nodeToDof[n[i]];
            if (dof >= 0) atomicAdd(&rhs[dof], contrib);
        }
    }
}

template<typename RealType>
__global__ void enforceBCKernel(const uint8_t* isBoundaryDof,
                                const int* rowPtr,
                                const int* /*colInd*/,
                                const int* diagPtr,
                                RealType* values,
                                RealType* rhs,
                                int numDofs)
{
    int dof = blockIdx.x * blockDim.x + threadIdx.x;
    if (dof >= numDofs || !isBoundaryDof[dof]) return;

    for (int j = rowPtr[dof]; j < rowPtr[dof + 1]; ++j)
        values[j] = RealType(0);

    int dp = diagPtr[dof];
    if (dp >= 0) values[dp] = RealType(1);

    rhs[dof] = RealType(0);
}

template<typename KeyType, typename RealType>
void solvePoissonGraphTet(ElementDomain<TetTag, RealType, KeyType, cstone::GpuTag>& domain,
                          RealType sourceTerm,
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

    const auto& d_nodeOwnership = domain.getNodeOwnershipMap();
    const auto& d_conn          = domain.getElementToNodeConnectivity();
    domain.cacheNodeCoordinates();
    const auto& d_x = domain.getNodeX();
    const auto& d_y = domain.getNodeY();
    const auto& d_z = domain.getNodeZ();
    pt.lap("lazy domain prep");

    cstone::DeviceVector<int> d_node_to_dof(nodeCount);
    int numOwnedDofs = buildDofMappingGpu<KeyType>(d_nodeOwnership.data(), d_node_to_dof.data(), nodeCount);
    int numTotalDofs = static_cast<int>(nodeCount);
    pt.lap("DOF mapping");

    // FULL 4x4 sparsity for tet (16-NNZ/element pattern that matches
    // CvfemTetAssembler::assembleFull). Using graph sparsity would silently
    // drop entries via addValue and corrupt the matrix.
    cstone::DeviceVector<int> d_rowPtr(numTotalDofs + 1);
    cstone::DeviceVector<int> d_colInd;
    cstone::DeviceVector<int> d_diagPtr(numTotalDofs);

    int nnz = CvfemTetSparsityBuilder<KeyType>::buildFullSparsity(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
        std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
        elementCount, d_node_to_dof.data(), numTotalDofs,
        d_rowPtr.data(), nullptr, nullptr, 0);

    d_colInd.resize(nnz);
    CvfemTetSparsityBuilder<KeyType>::buildFullSparsity(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
        std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
        elementCount, d_node_to_dof.data(), numTotalDofs,
        d_rowPtr.data(), d_colInd.data(), d_diagPtr.data(), 0);

    cstone::DeviceVector<RealType> d_values(nnz, RealType(0));
    cstone::DeviceVector<RealType> d_rhs(numTotalDofs, RealType(0));
    pt.lap("sparsity build");

    using MatrixType = CSRMatrix<RealType>;
    MatrixType* d_matrix;
    cudaMalloc(&d_matrix, sizeof(MatrixType));
    MatrixType h_matrix{d_rowPtr.data(), d_colInd.data(), d_values.data(), d_diagPtr.data(), numTotalDofs, nnz};
    cudaMemcpy(d_matrix, &h_matrix, sizeof(MatrixType), cudaMemcpyHostToDevice);

    // Poisson fields: gamma=1, phi=0, beta=0, mdot=0.
    cstone::DeviceVector<RealType> d_gamma(nodeCount, RealType(1));
    cstone::DeviceVector<RealType> d_phi(nodeCount, RealType(0));
    cstone::DeviceVector<RealType> d_beta(nodeCount, RealType(0));
    cstone::DeviceVector<RealType> d_grad_phi_x(nodeCount, RealType(0));
    cstone::DeviceVector<RealType> d_grad_phi_y(nodeCount, RealType(0));
    cstone::DeviceVector<RealType> d_grad_phi_z(nodeCount, RealType(0));
    cstone::DeviceVector<RealType> d_mdot(elementCount * 6, RealType(0));

    typename CvfemTetAssembler<KeyType, RealType>::Config config;
    config.blockSize = blockSize;

    CvfemTetAssembler<KeyType, RealType>::assembleFull(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
        std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
        elementCount,
        d_x.data(), d_y.data(), d_z.data(),
        d_gamma.data(), d_phi.data(), d_beta.data(),
        d_grad_phi_x.data(), d_grad_phi_y.data(), d_grad_phi_z.data(),
        d_mdot.data(),
        d_node_to_dof.data(), d_nodeOwnership.data(),
        d_matrix, d_rhs.data(), config);
    cudaDeviceSynchronize();
    pt.lap("assembly (LHS+RHS)");

    {
        int eBlocks = (elementCount + blockSize - 1) / blockSize;
        integrateConstantSourceTetKernel<KeyType, RealType><<<eBlocks, blockSize>>>(
            std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
            std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
            d_x.data(), d_y.data(), d_z.data(),
            d_node_to_dof.data(), d_nodeOwnership.data(),
            d_rhs.data(), elementCount, sourceTerm);
        cudaDeviceSynchronize();
    }
    pt.lap("source integration");

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

    using Matrix = SparseMatrix<int, RealType, cstone::GpuTag>;
    Matrix A;
    A.allocate(numOwnedDofs, numTotalDofs, nnz);
    cudaMemcpy(A.rowOffsetsPtr(),  d_rowPtr.data(),  (numOwnedDofs + 1) * sizeof(int),     cudaMemcpyDeviceToDevice);
    cudaMemcpy(A.colIndicesPtr(),  d_colInd.data(),  nnz * sizeof(int),                    cudaMemcpyDeviceToDevice);
    cudaMemcpy(A.valuesPtr(),      d_values.data(),  nnz * sizeof(RealType),               cudaMemcpyDeviceToDevice);

    cstone::DeviceVector<RealType> b(numOwnedDofs), x(numTotalDofs);
    cudaMemcpy(b.data(), d_rhs.data(), numOwnedDofs * sizeof(RealType), cudaMemcpyDeviceToDevice);

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
    solver.setOwnedSize(numOwnedDofs);

    int nRanks = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &nRanks);
    if (nRanks > 1)
    {
        const int* dofMapPtr = d_node_to_dof.data();
        solver.setHaloExchangeCallback(
            [&domain, dofMapPtr](cstone::DeviceVector<RealType>& p) { domain.exchangeNodeHalo(p, dofMapPtr); });
    }

    pt.lap("CG setup (Matrix alloc + halo cb)");

    solver.solve(A, b, x);
    cudaDeviceSynchronize();
    pt.lap("CG solve");

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

    domain.exchangeNodeHalo(d_nodeSolution);

    d_errorPerElement = TetErrorIndicator<KeyType, RealType>::computeError(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
        std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
        d_nodeSolution.data(), d_x.data(), d_y.data(), d_z.data(),
        d_node_to_dof.data(), elementCount, blockSize);
    pt.lap("error indicator");

    cudaFree(d_matrix);

    auto& cstoneDom = domain.getDomain();
    RealType errNorm = ErrorIndicator<KeyType, RealType>::globalErrorNormOwned(
        d_errorPerElement, cstoneDom.startIndex(), cstoneDom.endIndex());

    auto x_ptr           = thrust::device_pointer_cast(x.data());
    RealType localSumSq  = thrust::inner_product(thrust::device, x_ptr, x_ptr + numOwnedDofs, x_ptr, RealType(0));
    RealType globalSumSq = 0;
    auto mpiTypeR        = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    MPI_Allreduce(&localSumSq, &globalSumSq, 1, mpiTypeR, MPI_SUM, MPI_COMM_WORLD);
    RealType solNorm = std::sqrt(globalSumSq);

    int globalDofs = 0;
    MPI_Allreduce(&numOwnedDofs, &globalDofs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0)
    {
        std::cout << "    Solve: " << globalDofs << " global DOFs (rank0=" << numOwnedDofs
                  << "), " << nnz << " NNZ, ||u||=" << std::scientific
                  << solNorm << ", ||err||=" << errNorm << "\n"
                  << std::defaultfloat;
    }

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

    std::string meshFile;
    int blockSize         = 256;
    int bucketSize        = 64;
    double sourceTerm     = 1.0;
    int maxIter           = 1000;
    double tolerance      = 1e-10;
    int amrLevels         = 3;
    double refineFraction = 0.3;
    double coarsenFraction = 0.03;

    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg.find("--mesh=") == 0)                 meshFile       = arg.substr(7);
        else if (arg.find("--block-size=") == 0)      blockSize      = std::stoi(arg.substr(13));
        else if (arg.find("--bucket-size=") == 0)     bucketSize     = std::stoi(arg.substr(14));
        else if (arg.find("--source=") == 0)          sourceTerm     = std::stod(arg.substr(9));
        else if (arg.find("--max-iter=") == 0)        maxIter        = std::stoi(arg.substr(11));
        else if (arg.find("--tol=") == 0)             tolerance      = std::stod(arg.substr(6));
        else if (arg.find("--amr-levels=") == 0)      amrLevels      = std::stoi(arg.substr(13));
        else if (arg.find("--refine-frac=") == 0)     refineFraction = std::stod(arg.substr(14));
        else if (arg.find("--coarsen-frac=") == 0)    coarsenFraction = std::stod(arg.substr(15));
        else if (arg[0] != '-' && meshFile.empty())   meshFile       = arg;
    }

    if (meshFile.empty())
    {
        if (rank == 0)
        {
            std::cout << "Usage: " << argv[0] << " --mesh=FILE [options]\n"
                      << "  --source=VALUE       source term f (default 1.0)\n"
                      << "  --max-iter=N         CG max iters (default 1000)\n"
                      << "  --tol=VALUE          CG tolerance (default 1e-10)\n"
                      << "  --amr-levels=N       max AMR levels (default 3)\n"
                      << "  --refine-frac=VALUE  refine fraction (default 0.3)\n"
                      << "  --coarsen-frac=VALUE coarsen fraction (default 0.03)\n"
                      << "  --bucket-size=N      cornerstone bucket size (default 64)\n"
                      << "  --block-size=N       CUDA block size (default 256)\n";
        }
        MPI_Finalize();
        return 1;
    }

    using KeyType  = uint64_t;
    using RealType = double;

    if (rank == 0)
    {
        std::cout << "\n========================================\n"
                  << "MARS AMR CVFEM Tet (Graph) Solver\n"
                  << "========================================\n"
                  << "Problem: -du = " << sourceTerm << ", u = 0 on boundary\n"
                  << "Mesh: " << meshFile << "\n"
                  << "MPI ranks: " << numRanks << "\n"
                  << "AMR levels: " << amrLevels << "\n"
                  << "Bucket size: " << bucketSize << "\n"
                  << "========================================\n\n";
    }

    AmrManager<TetTag, KeyType, RealType>::Config amrConfig;
    amrConfig.maxLevels       = amrLevels;
    amrConfig.refineFraction  = refineFraction;
    amrConfig.coarsenFraction = coarsenFraction;
    amrConfig.blockSize       = blockSize;
    amrConfig.bucketSize      = bucketSize;

    AmrManager<TetTag, KeyType, RealType> amr(amrConfig);
    amr.initialize(meshFile, rank, numRanks);
    auto initT = amr.initTimings();

    if (rank == 0)
    {
        std::cout << "Initial mesh:\n"
                  << "  Elements: " << amr.domain().getElementCount() << "\n"
                  << "  Nodes:    " << amr.domain().getNodeCount() << "\n"
                  << "  Init breakdown (ms):\n"
                  << "    domain (file read + sync):    " << std::fixed << initT.domainSyncTimeMs << "\n"
                  << "    halo + node topology:         " << initT.haloTopoTimeMs   << "\n"
                  << "    adjacency CSR (e2n + n2e):    " << initT.adjacencyTimeMs  << "\n"
                  << "    coord cache (SFC decode):     " << initT.coordCacheTimeMs << "\n"
                  << "    AMR octree state:             " << initT.octreeTimeMs     << "\n"
                  << "    TOTAL:                        " << initT.totalMs          << "\n\n";
    }

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
            std::cout << "--- AMR Level " << level << " ---\n"
                      << "  Elements: " << amr.domain().getElementCount()
                      << " (local: " << amr.domain().localElementCount() << ")\n"
                      << "  Nodes:    " << amr.domain().getNodeCount() << "\n";
        }

        solvePoissonGraphTet<KeyType, RealType>(amr.domain(), sourceTerm, blockSize, maxIter, tolerance,
                                                 rank, d_nodeSolution, d_errorPerElement);
        cudaDeviceSynchronize();
        MPI_Barrier(MPI_COMM_WORLD);

        auto levelEnd = std::chrono::high_resolution_clock::now();
        float levelMs = std::chrono::duration<float, std::milli>(levelEnd - levelStart).count();
        if (rank == 0) std::cout << "    Level time: " << std::fixed << levelMs << " ms\n";

        if (level >= amrLevels) break;

        cstone::DeviceVector<RealType> d_newSolution;
        auto amrStart = std::chrono::high_resolution_clock::now();
        stats = amr.adaptMesh(d_errorPerElement.data(), d_nodeSolution.data(), d_newSolution);
        cudaDeviceSynchronize();
        MPI_Barrier(MPI_COMM_WORLD);
        auto amrEnd = std::chrono::high_resolution_clock::now();
        float amrMs = std::chrono::duration<float, std::milli>(amrEnd - amrStart).count();
        AmrManager<TetTag, KeyType, RealType>::printStats(stats, rank);
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

        std::cout << "\n========================================\n"
                  << "Final Solution\n"
                  << "========================================\n"
                  << "  Min:     " << std::scientific << solMin << "\n"
                  << "  Max:     " << solMax << "\n"
                  << "  L2 norm: " << solNorm << "\n"
                  << "  AMR levels: " << amr.currentLevel() << "\n"
                  << "  Total time: " << std::fixed << totalMs << " ms\n"
                  << "========================================\n";
    }

    MPI_Finalize();
    return 0;
}
