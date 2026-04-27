// AMR Poisson Solver: -Δu = f with u = 0 on boundary
// Demonstrates adaptive mesh refinement with CVFEM assembly on cstone-based unstructured meshes.
// All mesh operations are GPU-native.

#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_fem.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_hex_kernel.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_assembler.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_utils.hpp"
#include "backend/distributed/unstructured/fem/mars_sparse_matrix.hpp"
#include "backend/distributed/unstructured/solvers/mars_cg_solver.hpp"
#include "backend/distributed/unstructured/amr/mars_amr.hpp"

#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/extrema.h>
#include <thrust/inner_product.h>
#include <thrust/transform.h>
#include <set>
#include <mpi.h>
#include <iomanip>
#include <chrono>
#include <cmath>

using namespace mars;
using namespace mars::fem;
using namespace mars::amr;

// Build DOF mapping: GPU kernel assigns DOF numbers to owned nodes
__global__ void buildDofMappingKernel(const uint8_t* ownership,
                                      int* nodeToDof,
                                      int* dofCounter,
                                      size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;

    if (ownership[i] == 1)
    {
        int dof      = atomicAdd(dofCounter, 1);
        nodeToDof[i] = dof;
    }
    else
    {
        nodeToDof[i] = -1;
    }
}

// Solve Poisson on a given domain, return solution per node
template<typename KeyType, typename RealType>
void solvePoisson(ElementDomain<HexTag, RealType, KeyType, cstone::GpuTag>& domain,
                  RealType sourceTerm,
                  CvfemKernelVariant kernelVariant,
                  int blockSize,
                  int maxIter,
                  RealType tolerance,
                  int rank,
                  // outputs
                  cstone::DeviceVector<RealType>& d_nodeSolution,
                  cstone::DeviceVector<int>& d_nodeToDof,
                  int& numOwnedDofs)
{
    size_t nodeCount    = domain.getNodeCount();
    size_t elementCount = domain.getElementCount();

    const auto& d_nodeOwnership = domain.getNodeOwnershipMap();
    const auto& d_conn          = domain.getConnectivity();

    // Build DOF mapping on GPU
    d_nodeToDof.resize(nodeCount);
    cstone::DeviceVector<int> d_counter(1, 0);

    int numBlocks = (nodeCount + blockSize - 1) / blockSize;
    buildDofMappingKernel<<<numBlocks, blockSize>>>(d_nodeOwnership.data(), d_nodeToDof.data(), d_counter.data(),
                                                     nodeCount);
    cudaDeviceSynchronize();
    cudaMemcpy(&numOwnedDofs, d_counter.data(), sizeof(int), cudaMemcpyDeviceToHost);

    if (numOwnedDofs == 0)
    {
        d_nodeSolution.resize(nodeCount);
        thrust::fill(d_nodeSolution.begin(), d_nodeSolution.end(), RealType(0));
        return;
    }

    // Build sparsity pattern (host-side, as in poisson example)
    std::vector<uint8_t> h_ownership(nodeCount);
    cudaMemcpy(h_ownership.data(), d_nodeOwnership.data(), nodeCount * sizeof(uint8_t), cudaMemcpyDeviceToHost);

    std::vector<int> h_node_to_dof(nodeCount);
    cudaMemcpy(h_node_to_dof.data(), d_nodeToDof.data(), nodeCount * sizeof(int), cudaMemcpyDeviceToHost);

    std::vector<KeyType> h_conn0(elementCount), h_conn1(elementCount), h_conn2(elementCount), h_conn3(elementCount);
    std::vector<KeyType> h_conn4(elementCount), h_conn5(elementCount), h_conn6(elementCount), h_conn7(elementCount);
    cudaMemcpy(h_conn0.data(), std::get<0>(d_conn).data(), elementCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_conn1.data(), std::get<1>(d_conn).data(), elementCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_conn2.data(), std::get<2>(d_conn).data(), elementCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_conn3.data(), std::get<3>(d_conn).data(), elementCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_conn4.data(), std::get<4>(d_conn).data(), elementCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_conn5.data(), std::get<5>(d_conn).data(), elementCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_conn6.data(), std::get<6>(d_conn).data(), elementCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_conn7.data(), std::get<7>(d_conn).data(), elementCount * sizeof(KeyType), cudaMemcpyDeviceToHost);

    std::vector<std::set<int>> adj(numOwnedDofs);
    for (size_t e = 0; e < elementCount; ++e)
    {
        KeyType nodes[8] = {h_conn0[e], h_conn1[e], h_conn2[e], h_conn3[e],
                            h_conn4[e], h_conn5[e], h_conn6[e], h_conn7[e]};
        int dofs[8];
        for (int i = 0; i < 8; ++i)
            dofs[i] = h_node_to_dof[nodes[i]];

        for (int i = 0; i < 8; ++i)
        {
            if (dofs[i] >= 0 && h_ownership[nodes[i]] == 1)
            {
                for (int j = 0; j < 8; ++j)
                {
                    if (dofs[j] >= 0) adj[dofs[i]].insert(dofs[j]);
                }
            }
        }
    }

    std::vector<int> rowPtr(numOwnedDofs + 1);
    std::vector<int> colInd;
    std::vector<int> diagPtr(numOwnedDofs);

    rowPtr[0] = 0;
    for (int d = 0; d < numOwnedDofs; ++d)
    {
        int diagFound = -1;
        for (int col : adj[d])
        {
            if (col == d) diagFound = colInd.size();
            colInd.push_back(col);
        }
        diagPtr[d]     = diagFound;
        rowPtr[d + 1] = colInd.size();
    }
    int nnz = colInd.size();

    // Allocate sparse matrix
    using Matrix = SparseMatrix<int, RealType, cstone::GpuTag>;
    Matrix A;
    A.allocate(numOwnedDofs, numOwnedDofs, nnz);

    thrust::copy(rowPtr.begin(), rowPtr.end(), thrust::device_pointer_cast(A.rowOffsetsPtr()));
    thrust::copy(colInd.begin(), colInd.end(), thrust::device_pointer_cast(A.colIndicesPtr()));
    thrust::fill(thrust::device_pointer_cast(A.valuesPtr()), thrust::device_pointer_cast(A.valuesPtr() + nnz),
                 RealType(0));

    cstone::DeviceVector<int> d_diagPtr(numOwnedDofs);
    cudaMemcpy(d_diagPtr.data(), diagPtr.data(), numOwnedDofs * sizeof(int), cudaMemcpyHostToDevice);
    CSRMatrix<RealType> assemblyMatrix;
    assemblyMatrix.numRows = numOwnedDofs;
    assemblyMatrix.nnz     = nnz;
    assemblyMatrix.rowPtr  = A.rowOffsetsPtr();
    assemblyMatrix.colInd  = A.colIndicesPtr();
    assemblyMatrix.values  = A.valuesPtr();
    assemblyMatrix.diagPtr = d_diagPtr.data();

    cstone::DeviceVector<RealType> d_rhs(numOwnedDofs, 0.0);

    // Initialize fields
    domain.cacheNodeCoordinates();
    const auto& d_x = domain.getNodeX();
    const auto& d_y = domain.getNodeY();
    const auto& d_z = domain.getNodeZ();

    cstone::DeviceVector<RealType> d_gamma(nodeCount, 1.0);
    cstone::DeviceVector<RealType> d_phi(nodeCount, 0.0);
    cstone::DeviceVector<RealType> d_beta(nodeCount, 0.0);
    cstone::DeviceVector<RealType> d_grad_phi_x(nodeCount, 0.0);
    cstone::DeviceVector<RealType> d_grad_phi_y(nodeCount, 0.0);
    cstone::DeviceVector<RealType> d_grad_phi_z(nodeCount, 0.0);

    cstone::DeviceVector<RealType> d_areaVec_x(elementCount * 12);
    cstone::DeviceVector<RealType> d_areaVec_y(elementCount * 12);
    cstone::DeviceVector<RealType> d_areaVec_z(elementCount * 12);

    precomputeAreaVectorsGpu<KeyType, RealType>(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(), std::get<2>(d_conn).data(),
        std::get<3>(d_conn).data(), std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
        std::get<6>(d_conn).data(), std::get<7>(d_conn).data(), elementCount, d_x.data(), d_y.data(), d_z.data(),
        d_areaVec_x.data(), d_areaVec_y.data(), d_areaVec_z.data());

    cstone::DeviceVector<RealType> d_mdot(elementCount * 12, 0.0);

    // Assemble
    typename CvfemHexAssembler<KeyType, RealType>::Config config;
    config.blockSize = blockSize;
    config.variant   = kernelVariant;

    CvfemHexAssembler<KeyType, RealType>::assembleFull(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(), std::get<2>(d_conn).data(),
        std::get<3>(d_conn).data(), std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
        std::get<6>(d_conn).data(), std::get<7>(d_conn).data(), elementCount, d_x.data(), d_y.data(), d_z.data(),
        d_gamma.data(), d_phi.data(), d_beta.data(), d_grad_phi_x.data(), d_grad_phi_y.data(), d_grad_phi_z.data(),
        d_mdot.data(), d_areaVec_x.data(), d_areaVec_y.data(), d_areaVec_z.data(), d_nodeToDof.data(),
        d_nodeOwnership.data(), &assemblyMatrix, d_rhs.data(), config);
    cudaDeviceSynchronize();

    // Source term
    thrust::transform(d_rhs.begin(), d_rhs.end(), d_rhs.begin(),
                       [sourceTerm] __device__(RealType x) { return -sourceTerm + x; });

    // Apply boundary conditions
    std::vector<RealType> h_x(nodeCount), h_y(nodeCount), h_z(nodeCount);
    cudaMemcpy(h_x.data(), d_x.data(), nodeCount * sizeof(RealType), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_y.data(), d_y.data(), nodeCount * sizeof(RealType), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_z.data(), d_z.data(), nodeCount * sizeof(RealType), cudaMemcpyDeviceToHost);

    RealType xmin = *std::min_element(h_x.begin(), h_x.end());
    RealType xmax = *std::max_element(h_x.begin(), h_x.end());
    RealType ymin = *std::min_element(h_y.begin(), h_y.end());
    RealType ymax = *std::max_element(h_y.begin(), h_y.end());
    RealType zmin = *std::min_element(h_z.begin(), h_z.end());
    RealType zmax = *std::max_element(h_z.begin(), h_z.end());
    RealType eps  = 1e-10 * std::max({xmax - xmin, ymax - ymin, zmax - zmin});

    std::vector<bool> isBoundary(numOwnedDofs, false);
    for (size_t i = 0; i < nodeCount; ++i)
    {
        if (h_ownership[i] == 1)
        {
            bool onBdry = (std::abs(h_x[i] - xmin) < eps || std::abs(h_x[i] - xmax) < eps ||
                           std::abs(h_y[i] - ymin) < eps || std::abs(h_y[i] - ymax) < eps ||
                           std::abs(h_z[i] - zmin) < eps || std::abs(h_z[i] - zmax) < eps);
            if (onBdry && h_node_to_dof[i] >= 0)
            {
                isBoundary[h_node_to_dof[i]] = true;
            }
        }
    }

    std::vector<RealType> h_values(nnz);
    std::vector<RealType> h_rhs(numOwnedDofs);
    cudaMemcpy(h_values.data(), A.valuesPtr(), nnz * sizeof(RealType), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_rhs.data(), d_rhs.data(), numOwnedDofs * sizeof(RealType), cudaMemcpyDeviceToHost);

    for (int d = 0; d < numOwnedDofs; ++d)
    {
        if (isBoundary[d])
        {
            for (int j = rowPtr[d]; j < rowPtr[d + 1]; ++j)
                h_values[j] = 0.0;
            if (diagPtr[d] >= 0) h_values[diagPtr[d]] = 1.0;
            h_rhs[d] = 0.0;
        }
    }

    cudaMemcpy(A.valuesPtr(), h_values.data(), nnz * sizeof(RealType), cudaMemcpyHostToDevice);
    cudaMemcpy(d_rhs.data(), h_rhs.data(), numOwnedDofs * sizeof(RealType), cudaMemcpyHostToDevice);

    // Solve
    using Vector = cstone::DeviceVector<RealType>;
    Vector b(numOwnedDofs), x(numOwnedDofs);
    thrust::copy(d_rhs.begin(), d_rhs.end(), b.begin());
    thrust::fill(x.begin(), x.end(), RealType(0));

    ConjugateGradientSolver<RealType, int, cstone::GpuTag> solver(maxIter, tolerance);
    solver.setVerbose(false);
    solver.solve(A, b, x);
    cudaDeviceSynchronize();

    // Convert DOF solution to node solution
    d_nodeSolution.resize(nodeCount);
    thrust::fill(d_nodeSolution.begin(), d_nodeSolution.end(), RealType(0));

    // Scatter DOF values back to nodes
    std::vector<RealType> h_x_sol(numOwnedDofs);
    cudaMemcpy(h_x_sol.data(), x.data(), numOwnedDofs * sizeof(RealType), cudaMemcpyDeviceToHost);

    std::vector<RealType> h_nodeSol(nodeCount, 0.0);
    for (size_t i = 0; i < nodeCount; ++i)
    {
        if (h_node_to_dof[i] >= 0) h_nodeSol[i] = h_x_sol[h_node_to_dof[i]];
    }
    cudaMemcpy(d_nodeSolution.data(), h_nodeSol.data(), nodeCount * sizeof(RealType), cudaMemcpyHostToDevice);

    if (rank == 0)
    {
        RealType solNorm = std::sqrt(thrust::inner_product(x.begin(), x.end(), x.begin(), RealType(0)));
        std::cout << "    Solve: " << numOwnedDofs << " DOFs, " << nnz << " NNZ, ||u||=" << std::scientific << solNorm
                  << "\n"
                  << std::defaultfloat;
    }
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
            else if (v == "shmem") kernelVariant = CvfemKernelVariant::Shmem;
            else if (v == "optimized") kernelVariant = CvfemKernelVariant::Optimized;
            else if (v == "original") kernelVariant = CvfemKernelVariant::Original;
        }
        else if (arg.find("--block-size=") == 0)
            blockSize = std::stoi(arg.substr(13));
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
            std::cout << "  --kernel=VARIANT     tensor, shmem, optimized, original (default: tensor)\n";
            std::cout << "  --source=VALUE       Source term f (default: 1.0)\n";
            std::cout << "  --max-iter=N         CG max iterations (default: 1000)\n";
            std::cout << "  --tol=VALUE          CG tolerance (default: 1e-10)\n";
            std::cout << "  --amr-levels=N       Max AMR refinement levels (default: 3)\n";
            std::cout << "  --refine-frac=VALUE  Fraction of max error to refine (default: 0.3)\n";
            std::cout << "  --coarsen-frac=VALUE Fraction of max error to coarsen (default: 0.03)\n";
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
        std::cout << "MARS AMR Poisson Solver\n";
        std::cout << "========================================\n";
        std::cout << "Problem: -du = " << sourceTerm << ", u = 0 on boundary\n";
        std::cout << "Mesh: " << meshFile << "\n";
        std::cout << "Kernel: " << CvfemHexAssembler<KeyType, RealType>::variantName(kernelVariant) << "\n";
        std::cout << "MPI ranks: " << numRanks << "\n";
        std::cout << "AMR levels: " << amrLevels << "\n";
        std::cout << "Refine fraction: " << refineFraction << "\n";
        std::cout << "========================================\n\n";
    }

    // Configure AMR
    AmrManager<HexTag, KeyType, RealType>::Config amrConfig;
    amrConfig.maxLevels       = amrLevels;
    amrConfig.refineFraction  = refineFraction;
    amrConfig.coarsenFraction = coarsenFraction;
    amrConfig.blockSize       = blockSize;
    amrConfig.bucketSize      = 64;

    AmrManager<HexTag, KeyType, RealType> amr(amrConfig);
    amr.initialize(meshFile, rank, numRanks);

    if (rank == 0)
    {
        std::cout << "Initial mesh:\n";
        std::cout << "  Elements: " << amr.domain().getElementCount() << "\n";
        std::cout << "  Nodes:    " << amr.domain().getNodeCount() << "\n\n";
    }

    // AMR loop
    cstone::DeviceVector<RealType> d_nodeSolution;
    cstone::DeviceVector<int> d_nodeToDof;
    AmrStats stats;
    stats.elementsRefined = 1; // force at least one iteration

    for (int level = 0; level <= amrLevels; ++level)
    {
        if (rank == 0)
        {
            std::cout << "--- AMR Level " << level << " ---\n";
            std::cout << "  Elements: " << amr.domain().getElementCount() << "\n";
            std::cout << "  Nodes:    " << amr.domain().getNodeCount() << "\n";
        }

        // Solve on current mesh
        int numOwnedDofs = 0;
        solvePoisson<KeyType, RealType>(amr.domain(), sourceTerm, kernelVariant, blockSize, maxIter, tolerance, rank,
                                         d_nodeSolution, d_nodeToDof, numOwnedDofs);

        if (level >= amrLevels) break; // last level: just solve, don't refine

        // Compute error indicator
        size_t numElements = amr.domain().getElementCount();
        const auto& d_conn = amr.domain().getConnectivity();
        const auto& d_x    = amr.domain().getNodeX();
        const auto& d_y    = amr.domain().getNodeY();
        const auto& d_z    = amr.domain().getNodeZ();

        auto d_error = HexErrorIndicator<KeyType, RealType>::computeError(
            std::get<0>(d_conn).data(), std::get<1>(d_conn).data(), std::get<2>(d_conn).data(),
            std::get<3>(d_conn).data(), std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
            std::get<6>(d_conn).data(), std::get<7>(d_conn).data(), d_nodeSolution.data(), d_x.data(), d_y.data(),
            d_z.data(), d_nodeToDof.data(), numElements, blockSize);

        // Adapt mesh
        cstone::DeviceVector<RealType> d_newSolution;
        stats = amr.adaptMesh(d_error.data(), d_nodeSolution.data(), d_newSolution);
        AmrManager<HexTag, KeyType, RealType>::printStats(stats, rank);

        d_nodeSolution = std::move(d_newSolution);

        if (!amr.shouldContinue(stats))
        {
            if (rank == 0) std::cout << "\n  AMR converged or no more refinement needed.\n";
            break;
        }

        if (rank == 0) std::cout << "\n";
    }

    // Final solution statistics
    if (rank == 0)
    {
        RealType solMin = thrust::reduce(d_nodeSolution.begin(), d_nodeSolution.end(),
                                          std::numeric_limits<RealType>::max(), thrust::minimum<RealType>());
        RealType solMax = thrust::reduce(d_nodeSolution.begin(), d_nodeSolution.end(),
                                          std::numeric_limits<RealType>::lowest(), thrust::maximum<RealType>());
        RealType solNorm = std::sqrt(
            thrust::inner_product(d_nodeSolution.begin(), d_nodeSolution.end(), d_nodeSolution.begin(), RealType(0)));

        std::cout << "\n========================================\n";
        std::cout << "Final Solution\n";
        std::cout << "========================================\n";
        std::cout << "  Min:     " << std::scientific << solMin << "\n";
        std::cout << "  Max:     " << solMax << "\n";
        std::cout << "  L2 norm: " << solNorm << "\n";
        std::cout << "  AMR levels completed: " << amr.currentLevel() << "\n";
        std::cout << "========================================\n";
    }

    MPI_Finalize();
    return 0;
}
