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
#include <set>
#include <mpi.h>
#include <iomanip>
#include <chrono>
#include <cmath>

using namespace mars;
using namespace mars::fem;

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
    RealType sourceTerm = 1.0;  // RHS: -Δu = f
    int maxIter = 1000;
    RealType tolerance = 1e-10;

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
    ElementDomain<ElemTag, RealType, KeyType, cstone::GpuTag> domain(meshFile, rank, numRanks, true);
    const auto& d_nodeOwnership = domain.getNodeOwnershipMap();
    cudaDeviceSynchronize();

    size_t nodeCount = domain.getNodeCount();
    size_t elementCount = domain.getElementCount();
    const auto& d_conn = domain.getConnectivity();
    const auto& d_coords = domain.getCoordinates();

    if (rank == 0) {
        std::cout << "Mesh loaded:\n";
        std::cout << "  Nodes:    " << nodeCount << "\n";
        std::cout << "  Elements: " << elementCount << "\n\n";
    }

    // Create node-to-DOF mapping (1 DOF per node for scalar Poisson)
    // DOFs are numbered only for owned nodes
    std::vector<int> h_node_to_dof(nodeCount, -1);
    std::vector<uint8_t> h_ownership(nodeCount);
    cudaMemcpy(h_ownership.data(), d_nodeOwnership.data(), nodeCount * sizeof(uint8_t), cudaMemcpyDeviceToHost);

    int numOwnedDofs = 0;
    for (size_t i = 0; i < nodeCount; ++i) {
        if (h_ownership[i] == 1) {
            h_node_to_dof[i] = numOwnedDofs++;
        }
    }

    cstone::DeviceVector<int> d_node_to_dof(h_node_to_dof.begin(), h_node_to_dof.end());

    if (rank == 0) {
        std::cout << "DOF mapping:\n";
        std::cout << "  Owned DOFs: " << numOwnedDofs << "\n\n";
    }

    // Build sparsity pattern (full 8×8 connectivity for hex8)
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

    // Build adjacency for owned DOFs
    std::vector<std::set<int>> adj(numOwnedDofs);
    for (size_t e = 0; e < elementCount; ++e) {
        KeyType nodes[8] = {h_conn0[e], h_conn1[e], h_conn2[e], h_conn3[e],
                            h_conn4[e], h_conn5[e], h_conn6[e], h_conn7[e]};
        int dofs[8];
        for (int i = 0; i < 8; ++i) {
            dofs[i] = h_node_to_dof[nodes[i]];
        }
        // Add edges for owned DOFs
        for (int i = 0; i < 8; ++i) {
            if (dofs[i] >= 0 && h_ownership[nodes[i]] == 1) {
                for (int j = 0; j < 8; ++j) {
                    if (dofs[j] >= 0) {
                        adj[dofs[i]].insert(dofs[j]);
                    }
                }
            }
        }
    }

    // Build CSR matrix
    std::vector<int> rowPtr(numOwnedDofs + 1);
    std::vector<int> colInd;
    std::vector<int> diagPtr(numOwnedDofs);

    rowPtr[0] = 0;
    for (int d = 0; d < numOwnedDofs; ++d) {
        int diagFound = -1;
        for (int col : adj[d]) {
            if (col == d) {
                diagFound = colInd.size();
            }
            colInd.push_back(col);
        }
        diagPtr[d] = diagFound;
        rowPtr[d + 1] = colInd.size();
    }

    int nnz = colInd.size();
    if (rank == 0) {
        std::cout << "Sparsity pattern:\n";
        std::cout << "  NNZ: " << nnz << "\n";
        std::cout << "  Avg NNZ/row: " << (double)nnz / numOwnedDofs << "\n\n";
    }

    // Allocate matrix using SparseMatrix class
    using Matrix = SparseMatrix<int, RealType, cstone::GpuTag>;
    Matrix A;
    A.allocate(numOwnedDofs, numOwnedDofs, nnz);

    // Copy sparsity pattern
    thrust::copy(rowPtr.begin(), rowPtr.end(),
                 thrust::device_pointer_cast(A.rowOffsetsPtr()));
    thrust::copy(colInd.begin(), colInd.end(),
                 thrust::device_pointer_cast(A.colIndicesPtr()));
    thrust::fill(thrust::device_pointer_cast(A.valuesPtr()),
                 thrust::device_pointer_cast(A.valuesPtr() + nnz),
                 RealType(0));

    // Create temporary CsrMatrix for assembly (CVFEM kernels need this format)
    cstone::DeviceVector<int> d_diagPtr(diagPtr.begin(), diagPtr.end());
    CsrMatrix<RealType> assemblyMatrix;
    assemblyMatrix.numRows = numOwnedDofs;
    assemblyMatrix.numCols = numOwnedDofs;
    assemblyMatrix.nnz = nnz;
    assemblyMatrix.rowPtr = A.rowOffsetsPtr();
    assemblyMatrix.colInd = A.colIndicesPtr();
    assemblyMatrix.values = A.valuesPtr();
    assemblyMatrix.diagPtr = d_diagPtr.data();

    cstone::DeviceVector<RealType> d_rhs(numOwnedDofs, 0.0);

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
        std::get<0>(d_coords).data(), std::get<1>(d_coords).data(), std::get<2>(d_coords).data(),
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
        std::get<0>(d_coords).data(), std::get<1>(d_coords).data(), std::get<2>(d_coords).data(),
        d_gamma.data(), d_phi.data(), d_beta.data(),
        d_grad_phi_x.data(), d_grad_phi_y.data(), d_grad_phi_z.data(),
        d_mdot.data(),
        d_areaVec_x.data(), d_areaVec_y.data(), d_areaVec_z.data(),
        d_node_to_dof.data(),
        d_nodeOwnership.data(),
        &assemblyMatrix,
        d_rhs.data(),
        config
    );

    cudaDeviceSynchronize();
    auto assemblyEnd = std::chrono::high_resolution_clock::now();
    float assemblyTime = std::chrono::duration<float, std::milli>(assemblyEnd - assemblyStart).count();

    // Add source term to RHS: RHS = -f (since we have -Δu on LHS)
    thrust::transform(d_rhs.begin(), d_rhs.end(), d_rhs.begin(),
                      [sourceTerm] __device__ (RealType x) { return -sourceTerm + x; });

    if (rank == 0) {
        std::cout << "Assembly completed in " << assemblyTime << " ms\n\n";
    }

    // Apply boundary conditions: u = 0 on boundary
    // Detect boundary nodes geometrically
    std::vector<RealType> h_x(nodeCount), h_y(nodeCount), h_z(nodeCount);
    cudaMemcpy(h_x.data(), std::get<0>(d_coords).data(), nodeCount * sizeof(RealType), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_y.data(), std::get<1>(d_coords).data(), nodeCount * sizeof(RealType), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_z.data(), std::get<2>(d_coords).data(), nodeCount * sizeof(RealType), cudaMemcpyDeviceToHost);

    // Find bounding box
    RealType xmin = *std::min_element(h_x.begin(), h_x.end());
    RealType xmax = *std::max_element(h_x.begin(), h_x.end());
    RealType ymin = *std::min_element(h_y.begin(), h_y.end());
    RealType ymax = *std::max_element(h_y.begin(), h_y.end());
    RealType zmin = *std::min_element(h_z.begin(), h_z.end());
    RealType zmax = *std::max_element(h_z.begin(), h_z.end());

    RealType eps = 1e-10 * std::max({xmax - xmin, ymax - ymin, zmax - zmin});

    // Mark boundary DOFs
    std::vector<bool> isBoundary(numOwnedDofs, false);
    int numBoundaryDofs = 0;

    for (size_t i = 0; i < nodeCount; ++i) {
        if (h_ownership[i] == 1) {
            bool onBoundary = (std::abs(h_x[i] - xmin) < eps || std::abs(h_x[i] - xmax) < eps ||
                               std::abs(h_y[i] - ymin) < eps || std::abs(h_y[i] - ymax) < eps ||
                               std::abs(h_z[i] - zmin) < eps || std::abs(h_z[i] - zmax) < eps);
            if (onBoundary) {
                int dof = h_node_to_dof[i];
                if (dof >= 0) {
                    isBoundary[dof] = true;
                    numBoundaryDofs++;
                }
            }
        }
    }

    if (rank == 0) {
        std::cout << "Boundary conditions:\n";
        std::cout << "  Boundary DOFs: " << numBoundaryDofs << "\n\n";
    }

    // Apply BC by zeroing rows and setting diagonal to 1, RHS to 0
    std::vector<RealType> h_values(nnz);
    std::vector<RealType> h_rhs(numOwnedDofs);
    cudaMemcpy(h_values.data(), A.valuesPtr(), nnz * sizeof(RealType), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_rhs.data(), d_rhs.data(), numOwnedDofs * sizeof(RealType), cudaMemcpyDeviceToHost);

    for (int dof = 0; dof < numOwnedDofs; ++dof) {
        if (isBoundary[dof]) {
            // Zero row
            for (int j = rowPtr[dof]; j < rowPtr[dof + 1]; ++j) {
                h_values[j] = 0.0;
            }
            // Set diagonal to 1
            if (diagPtr[dof] >= 0) {
                h_values[diagPtr[dof]] = 1.0;
            }
            // Set RHS to 0
            h_rhs[dof] = 0.0;
        }
    }

    cudaMemcpy(A.valuesPtr(), h_values.data(), nnz * sizeof(RealType), cudaMemcpyHostToDevice);
    cudaMemcpy(d_rhs.data(), h_rhs.data(), numOwnedDofs * sizeof(RealType), cudaMemcpyHostToDevice);

    // Solve using CG
    using Vector = cstone::DeviceVector<RealType>;
    Vector b(numOwnedDofs), x(numOwnedDofs);
    thrust::copy(d_rhs.begin(), d_rhs.end(), b.begin());
    thrust::fill(x.begin(), x.end(), RealType(0));

    if (rank == 0) {
        std::cout << "Solving with CG...\n";
    }

    auto solveStart = std::chrono::high_resolution_clock::now();

    ConjugateGradientSolver<RealType, int, cstone::GpuTag> solver(maxIter, tolerance);
    solver.setVerbose(rank == 0);  // Only rank 0 prints
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

    // Compute solution statistics
    RealType solMin = thrust::reduce(x.begin(), x.end(),
                                     std::numeric_limits<RealType>::max(), thrust::minimum<RealType>());
    RealType solMax = thrust::reduce(x.begin(), x.end(),
                                     std::numeric_limits<RealType>::lowest(), thrust::maximum<RealType>());
    RealType solNorm = std::sqrt(thrust::inner_product(x.begin(), x.end(), x.begin(), 0.0));

    if (rank == 0) {
        std::cout << "Solution statistics:\n";
        std::cout << "  Min:  " << std::scientific << solMin << "\n";
        std::cout << "  Max:  " << solMax << "\n";
        std::cout << "  L2 norm: " << solNorm << "\n";
        std::cout << "========================================\n";
    }

    MPI_Finalize();
    return 0;
}
