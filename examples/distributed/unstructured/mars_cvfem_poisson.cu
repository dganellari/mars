// CVFEM Poisson Solver: -Δu = f with u = 0 on boundary
// Validates MARS CVFEM against MFEM ex0/ex1

#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_fem.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_hex_kernel.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_assembler.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_utils.hpp"
#include "backend/distributed/unstructured/fem/mars_sparse_matrix.hpp"
#include "backend/distributed/unstructured/solvers/mars_cg_solver.hpp"
#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/extrema.h>
#include <thrust/inner_product.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/for_each.h>
#include <algorithm>
#include <limits>
#include <mpi.h>
#include <iomanip>
#include <chrono>
#include <cmath>

using namespace mars;
using namespace mars::fem;

// Lumped control volume per node: each hex adds 1/8 of its volume to every corner. The trilinear
// hex volume is integrated exactly by 2x2x2 Gauss (detJ is quadratic in each reference direction).
template<typename KeyType, typename RealType>
__global__ void lumpedNodeVolumeKernel(const KeyType* c0, const KeyType* c1, const KeyType* c2,
                                       const KeyType* c3, const KeyType* c4, const KeyType* c5,
                                       const KeyType* c6, const KeyType* c7,
                                       const RealType* x, const RealType* y, const RealType* z,
                                       RealType* nodeVol, size_t numElements)
{
    size_t e = blockIdx.x * size_t(blockDim.x) + threadIdx.x;
    if (e >= numElements) return;
    const KeyType n[8] = {c0[e], c1[e], c2[e], c3[e], c4[e], c5[e], c6[e], c7[e]};
    const RealType s[8][3] = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
                              {-1, -1, 1},  {1, -1, 1},  {1, 1, 1},  {-1, 1, 1}};
    const RealType g = RealType(0.5773502691896257);
    RealType vol = 0;
    for (int q = 0; q < 8; ++q)
    {
        RealType a = s[q][0] * g, b = s[q][1] * g, c = s[q][2] * g;
        RealType J[3][3] = {};
        for (int i = 0; i < 8; ++i)
        {
            RealType dN[3] = {s[i][0] * (1 + s[i][1] * b) * (1 + s[i][2] * c) / 8,
                              s[i][1] * (1 + s[i][0] * a) * (1 + s[i][2] * c) / 8,
                              s[i][2] * (1 + s[i][0] * a) * (1 + s[i][1] * b) / 8};
            RealType p[3] = {x[n[i]], y[n[i]], z[n[i]]};
            for (int r = 0; r < 3; ++r)
                for (int t = 0; t < 3; ++t) J[r][t] += p[r] * dN[t];
        }
        RealType det = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
                     - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
                     + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
        vol += fabs(det);
    }
    for (int i = 0; i < 8; ++i) atomicAdd(&nodeVol[n[i]], vol / 8);
}

int main(int argc, char** argv) {
    // Initialize MPI
    MPI_Init(&argc, &argv);

    int rank, numRanks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    // Set CUDA device based on local MPI rank
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount > 0) {
        int device = rank % deviceCount;
        cudaSetDevice(device);
    }

    // Parse command-line options
    std::string meshFile;
    CvfemKernelVariant kernelVariant = CvfemKernelVariant::Tensor;
    int blockSize = 256;
    double sourceTerm = 1.0;  // RHS: -Δu = f
    int maxIter = 1000;
    double tolerance = 1e-10;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.find("--mesh=") == 0) {
            meshFile = arg.substr(7);
        } else if (arg.find("--kernel=") == 0) {
            std::string v = arg.substr(9);
            if (v == "tensor") kernelVariant = CvfemKernelVariant::Tensor;
            else if (v == "shmem") kernelVariant = CvfemKernelVariant::Shmem;
            else if (v == "optimized") kernelVariant = CvfemKernelVariant::Optimized;
            else if (v == "original") kernelVariant = CvfemKernelVariant::Original;
            else if (v == "wmma_tensor")  kernelVariant = CvfemKernelVariant::WmmaTensor;
            else if (v == "wgmma_tensor") kernelVariant = CvfemKernelVariant::WgmmaTensor;
        } else if (arg.find("--block-size=") == 0) {
            blockSize = std::stoi(arg.substr(13));
        } else if (arg.find("--source=") == 0) {
            sourceTerm = std::stod(arg.substr(9));
        } else if (arg.find("--max-iter=") == 0) {
            maxIter = std::stoi(arg.substr(11));
        } else if (arg.find("--tol=") == 0) {
            tolerance = std::stod(arg.substr(6));
        } else if (arg[0] != '-' && meshFile.empty()) {
            meshFile = arg;
        }
    }

    if (meshFile.empty()) {
        if (rank == 0) {
            std::cout << "Usage: " << argv[0] << " [options]\n";
            std::cout << "\nOptions:\n";
            std::cout << "  --mesh=FILE         Mesh file (.mesh or .exo format) [REQUIRED]\n";
            std::cout << "  --kernel=VARIANT    tensor, shmem, optimized, original (default: tensor)\n";
            std::cout << "  --source=VALUE      Source term f (default: 1.0)\n";
            std::cout << "  --max-iter=N        CG max iterations (default: 1000)\n";
            std::cout << "  --tol=VALUE         CG tolerance (default: 1e-10)\n";
            std::cout << "  --block-size=N      CUDA block size (default: 256)\n";
        }
        MPI_Finalize();
        return 1;
    }

    using KeyType = uint64_t;
    using RealType = double;
    using ElemTag = HexTag;

    if (rank == 0) {
        std::cout << "\n========================================\n";
        std::cout << "MARS CVFEM Poisson Solver\n";
        std::cout << "========================================\n";
        std::cout << "Problem: -Δu = " << sourceTerm << ", u = 0 on boundary\n";
        std::cout << "Mesh: " << meshFile << "\n";
        std::cout << "Kernel: " << CvfemHexAssembler<KeyType, RealType>::variantName(kernelVariant) << "\n";
        std::cout << "MPI ranks: " << numRanks << "\n";
        std::cout << "========================================\n\n";
    }

    // Load mesh and create domain
    ElementDomain<ElemTag, RealType, KeyType, cstone::execution::Gpu> domain(meshFile, rank, numRanks, true);
    const auto& d_nodeOwnership = domain.getNodeOwnershipMap();

    size_t nodeCount = domain.getNodeCount();
    size_t elementCount = domain.getElementCount();
    const auto& d_conn = domain.getElementToNodeConnectivity();
    domain.cacheNodeCoordinates();
    const auto& d_x = domain.getNodeX();
    const auto& d_y = domain.getNodeY();
    const auto& d_z = domain.getNodeZ();

    if (rank == 0) {
        std::cout << "Mesh loaded (rank 0):\n";
        std::cout << "  Nodes:    " << nodeCount << "\n";
        std::cout << "  Elements: " << elementCount << "\n\n";
    }

    // Owned nodes first, then ghosts, from the domain's node ownership (as the Navier-Stokes solvers)
    cstone::DeviceVector<int> d_node_to_dof(nodeCount);
    int numOwnedDofs = buildDofMappingGpu<KeyType>(d_nodeOwnership.data(), d_node_to_dof.data(), nodeCount);
    int numTotalDofs = int(nodeCount);
    long numGlobalDofs = numOwnedDofs;
    MPI_Allreduce(MPI_IN_PLACE, &numGlobalDofs, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "DOF mapping:\n";
        std::cout << "  Global DOFs: " << numGlobalDofs << "\n\n";
    }

    // Full 8x8 sparsity for every local DOF, ghost columns included; the owned rows come first
    cstone::DeviceVector<int> d_rowPtr(numTotalDofs + 1);
    cstone::DeviceVector<int> d_diagPtr(numTotalDofs);
    int nnzAll = CvfemSparsityBuilder<KeyType>::buildFullSparsity(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(), std::get<2>(d_conn).data(),
        std::get<3>(d_conn).data(), std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
        std::get<6>(d_conn).data(), std::get<7>(d_conn).data(), elementCount, d_node_to_dof.data(), numTotalDofs,
        d_rowPtr.data(), nullptr, nullptr, 0);
    cstone::DeviceVector<int> d_colInd(nnzAll);
    CvfemSparsityBuilder<KeyType>::buildFullSparsity(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(), std::get<2>(d_conn).data(),
        std::get<3>(d_conn).data(), std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
        std::get<6>(d_conn).data(), std::get<7>(d_conn).data(), elementCount, d_node_to_dof.data(), numTotalDofs,
        d_rowPtr.data(), d_colInd.data(), d_diagPtr.data(), 0);
    int nnz = 0;  // owned rows only
    cudaMemcpy(&nnz, d_rowPtr.data() + numOwnedDofs, sizeof(int), cudaMemcpyDeviceToHost);

    if (rank == 0) {
        std::cout << "Sparsity pattern (rank 0 owned rows):\n";
        std::cout << "  NNZ: " << nnz << "\n";
        std::cout << "  Avg NNZ/row: " << (double)nnz / std::max(numOwnedDofs, 1) << "\n\n";
    }

    cstone::DeviceVector<RealType> d_values(nnzAll, RealType(0));
    cstone::DeviceVector<RealType> d_rhs(numTotalDofs, RealType(0));

    // The assembly kernels read the matrix descriptor on the device; ghost rows exist only so ghost
    // columns are addressable and are skipped
    CSRMatrix<RealType> h_matrix{d_rowPtr.data(), d_colInd.data(), d_values.data(), d_diagPtr.data(),
                                 numTotalDofs, nnzAll};
    h_matrix.numOwnedRows = numOwnedDofs;
    CSRMatrix<RealType>* d_matrix = nullptr;
    cudaMalloc(&d_matrix, sizeof(CSRMatrix<RealType>));
    cudaMemcpy(d_matrix, &h_matrix, sizeof(CSRMatrix<RealType>), cudaMemcpyHostToDevice);

    // Initialize fields for CVFEM assembly
    // For Poisson equation: -Δu = f
    // In CVFEM advection-diffusion form: ∂φ/∂t + ∇·(βφu - γ∇φ) = 0
    // Set: β = 0 (no advection), γ = 1 (diffusion), source = f

    cstone::DeviceVector<RealType> d_gamma(nodeCount, 1.0);  // Diffusion coefficient
    cstone::DeviceVector<RealType> d_phi(nodeCount, 0.0);    // Solution (initial guess)
    cstone::DeviceVector<RealType> d_beta(nodeCount, 0.0);   // No advection
    cstone::DeviceVector<RealType> d_grad_phi_x(nodeCount, 0.0);
    cstone::DeviceVector<RealType> d_grad_phi_y(nodeCount, 0.0);
    cstone::DeviceVector<RealType> d_grad_phi_z(nodeCount, 0.0);

    // Precompute area vectors
    cstone::DeviceVector<RealType> d_areaVec_x(elementCount * 12);
    cstone::DeviceVector<RealType> d_areaVec_y(elementCount * 12);
    cstone::DeviceVector<RealType> d_areaVec_z(elementCount * 12);

    precomputeAreaVectorsGpu<KeyType, RealType>(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
        std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
        std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
        std::get<6>(d_conn).data(), std::get<7>(d_conn).data(),
        elementCount,
        d_x.data(), d_y.data(), d_z.data(),
        d_areaVec_x.data(), d_areaVec_y.data(), d_areaVec_z.data()
    );

    // mdot = 0 (no mass flux for Poisson)
    cstone::DeviceVector<RealType> d_mdot(elementCount * 12, 0.0);

    if (rank == 0) {
        std::cout << "Assembling system...\n";
    }

    // Assemble the system
    CvfemHexAssembler<KeyType, RealType>::Config config;
    config.blockSize = blockSize;
    config.variant = kernelVariant;

    auto assemblyStart = std::chrono::high_resolution_clock::now();

    CvfemHexAssembler<KeyType, RealType>::assembleFull(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
        std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
        std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
        std::get<6>(d_conn).data(), std::get<7>(d_conn).data(),
        elementCount,
        d_x.data(), d_y.data(), d_z.data(),
        d_gamma.data(), d_phi.data(), d_beta.data(),
        d_grad_phi_x.data(), d_grad_phi_y.data(), d_grad_phi_z.data(),
        d_mdot.data(),
        d_areaVec_x.data(), d_areaVec_y.data(), d_areaVec_z.data(),
        d_node_to_dof.data(),
        d_nodeOwnership.data(),
        d_matrix,
        d_rhs.data(),
        config
    );

    cudaDeviceSynchronize();
    cudaFree(d_matrix);
    auto assemblyEnd = std::chrono::high_resolution_clock::now();
    float assemblyTime = std::chrono::duration<float, std::milli>(assemblyEnd - assemblyStart).count();

    // Source term: the assembled operator is the positive-definite -Δ integrated over control
    // volumes, so the load of node i is +f * V_i (lumped control volume), not a bare -f. Every held
    // element adds to the volumes, so owned nodes get their complete control volume.
    cstone::DeviceVector<RealType> d_nodeVol(nodeCount, RealType(0));
    {
        int nb = int((elementCount + blockSize - 1) / blockSize);
        lumpedNodeVolumeKernel<KeyType, RealType><<<nb, blockSize>>>(
            std::get<0>(d_conn).data(), std::get<1>(d_conn).data(), std::get<2>(d_conn).data(),
            std::get<3>(d_conn).data(), std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
            std::get<6>(d_conn).data(), std::get<7>(d_conn).data(),
            d_x.data(), d_y.data(), d_z.data(), d_nodeVol.data(), elementCount);
        cudaDeviceSynchronize();
    }
    thrust::for_each(thrust::device, thrust::counting_iterator<size_t>(0),
                     thrust::counting_iterator<size_t>(nodeCount),
                     [rhs = d_rhs.data(), vol = d_nodeVol.data(), n2d = d_node_to_dof.data(), sourceTerm,
                      numOwnedDofs] __device__(size_t i)
                     {
                         int dof = n2d[i];
                         if (dof >= 0 && dof < numOwnedDofs) rhs[dof] += RealType(sourceTerm) * vol[i];
                     });

    if (rank == 0) {
        std::cout << "Assembly completed in " << assemblyTime << " ms\n\n";
    }

    // Boundary conditions u = 0 on the faces of the global bounding box: owned boundary rows become
    // identity rows with zero right-hand side
    const RealType* coords[3] = {d_x.data(), d_y.data(), d_z.data()};
    RealType lo[3], hi[3];
    for (int d = 0; d < 3; ++d) {
        lo[d] = thrust::reduce(thrust::device, coords[d], coords[d] + nodeCount,
                               std::numeric_limits<RealType>::max(), thrust::minimum<RealType>());
        hi[d] = thrust::reduce(thrust::device, coords[d], coords[d] + nodeCount,
                               std::numeric_limits<RealType>::lowest(), thrust::maximum<RealType>());
    }
    MPI_Allreduce(MPI_IN_PLACE, lo, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, hi, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    RealType eps = 1e-10 * std::max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]});

    cstone::DeviceVector<uint8_t> d_isBoundary(numOwnedDofs, 0);
    thrust::for_each(thrust::device, thrust::counting_iterator<size_t>(0),
                     thrust::counting_iterator<size_t>(nodeCount),
                     [x = d_x.data(), y = d_y.data(), z = d_z.data(), n2d = d_node_to_dof.data(),
                      isB = d_isBoundary.data(), numOwnedDofs, eps, x0 = lo[0], x1 = hi[0], y0 = lo[1],
                      y1 = hi[1], z0 = lo[2], z1 = hi[2]] __device__(size_t i)
                     {
                         int dof = n2d[i];
                         if (dof < 0 || dof >= numOwnedDofs) return;
                         bool onBoundary = fabs(x[i] - x0) < eps || fabs(x[i] - x1) < eps ||
                                           fabs(y[i] - y0) < eps || fabs(y[i] - y1) < eps ||
                                           fabs(z[i] - z0) < eps || fabs(z[i] - z1) < eps;
                         if (onBoundary) isB[dof] = 1;
                     });
    thrust::for_each(thrust::device, thrust::counting_iterator<int>(0),
                     thrust::counting_iterator<int>(numOwnedDofs),
                     [rowPtr = d_rowPtr.data(), diagPtr = d_diagPtr.data(), values = d_values.data(),
                      rhs = d_rhs.data(), isB = d_isBoundary.data()] __device__(int dof)
                     {
                         if (!isB[dof]) return;
                         for (int j = rowPtr[dof]; j < rowPtr[dof + 1]; ++j) values[j] = RealType(0);
                         values[diagPtr[dof]] = RealType(1);
                         rhs[dof] = RealType(0);
                     });
    long numBoundaryDofs = thrust::reduce(thrust::device, d_isBoundary.data(), d_isBoundary.data() + numOwnedDofs, 0L);
    MPI_Allreduce(MPI_IN_PLACE, &numBoundaryDofs, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "Boundary conditions:\n";
        std::cout << "  Boundary DOFs: " << numBoundaryDofs << "\n\n";
    }

    // CG on the owned rows; columns include ghosts, refreshed through the node halo
    using Matrix = SparseMatrix<int, RealType, cstone::execution::Gpu>;
    Matrix A;
    A.allocate(numOwnedDofs, numTotalDofs, nnz);
    cudaMemcpy(A.rowOffsetsPtr(), d_rowPtr.data(), (numOwnedDofs + 1) * sizeof(int), cudaMemcpyDeviceToDevice);
    cudaMemcpy(A.colIndicesPtr(), d_colInd.data(), size_t(nnz) * sizeof(int), cudaMemcpyDeviceToDevice);
    cudaMemcpy(A.valuesPtr(), d_values.data(), size_t(nnz) * sizeof(RealType), cudaMemcpyDeviceToDevice);

    using Vector = cstone::DeviceVector<RealType>;
    Vector b(numOwnedDofs), x(numTotalDofs);
    cudaMemcpy(b.data(), d_rhs.data(), numOwnedDofs * sizeof(RealType), cudaMemcpyDeviceToDevice);
    thrust::fill(thrust::device, x.begin(), x.end(), RealType(0));

    if (rank == 0) {
        std::cout << "Solving with CG...\n";
    }

    auto solveStart = std::chrono::high_resolution_clock::now();

    ConjugateGradientSolver<RealType, int, cstone::execution::Gpu> solver(maxIter, tolerance);
    solver.setVerbose(rank == 0);  // Only rank 0 prints
    solver.setOwnedSize(numOwnedDofs);
    if (numRanks > 1) {
        const int* dofMap = d_node_to_dof.data();
        solver.setHaloExchangeCallback([&domain, dofMap](Vector& p) { domain.exchangeNodeHalo(p, dofMap); });
    }
    bool converged = solver.solve(A, b, x);

    cudaDeviceSynchronize();
    auto solveEnd = std::chrono::high_resolution_clock::now();
    float solveTime = std::chrono::duration<float, std::milli>(solveEnd - solveStart).count();

    if (rank == 0) {
        std::cout << "\n========================================\n";
        std::cout << "Solver Results\n";
        std::cout << "========================================\n";
        std::cout << "Converged: " << (converged ? "YES" : "NO") << "\n";
        std::cout << "Solve time: " << std::fixed << solveTime << " ms\n";
        std::cout << "========================================\n\n";
    }

    // Solution statistics over the owned DOFs of all ranks
    const RealType* xOwned = x.data();
    RealType solMin = thrust::reduce(thrust::device, xOwned, xOwned + numOwnedDofs,
                                     std::numeric_limits<RealType>::max(), thrust::minimum<RealType>());
    RealType solMax = thrust::reduce(thrust::device, xOwned, xOwned + numOwnedDofs,
                                     std::numeric_limits<RealType>::lowest(), thrust::maximum<RealType>());
    RealType solSq  = thrust::inner_product(thrust::device, xOwned, xOwned + numOwnedDofs, xOwned, RealType(0));
    MPI_Allreduce(MPI_IN_PLACE, &solMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &solMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &solSq, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "Solution statistics:\n";
        std::cout << "  Min:  " << std::scientific << solMin << "\n";
        std::cout << "  Max:  " << solMax << "\n";
        std::cout << "  L2 norm: " << std::sqrt(solSq) << "\n";
        std::cout << "========================================\n";
    }

    MPI_Finalize();
    return converged ? 0 : 1;
}
