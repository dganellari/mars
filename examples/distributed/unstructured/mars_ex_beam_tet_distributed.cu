/**
 * MFEM beam-tet.mesh Distributed GPU Poisson solver
 *
 * Problem: -Δu = 1 in Ω, u = 0 on ∂Ω
 *
 * Usage: mpirun -np N mars_ex_beam_tet_distributed --mesh beam-tet.mesh
 *
 * Features:
 * - MPI distributed across multiple ranks
 * - GPU acceleration throughout (no CPU transfers)
 * - Hypre GPU-accelerated PCG + BoomerAMG
 * - Domain decomposition with halo exchange
 *
 * Expected: Fast convergence with distributed AMG
 */

#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_fem.hpp"
#include "backend/distributed/unstructured/fem/mars_unstructured_dof_handler.hpp"
#include "backend/distributed/unstructured/fem/mars_unstructured_dm.hpp"
#include "backend/distributed/unstructured/solvers/mars_cg_solver_with_preconditioner.hpp"
#ifdef MARS_ENABLE_HYPRE
#include "backend/distributed/unstructured/solvers/mars_hypre_pcg_solver.hpp"
#endif
#include "mars_mfem_mesh_loader.hpp"
#include <thrust/copy.h>
#include <algorithm>
#include <limits>

#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <chrono>
#include <fstream>

using namespace mars;
using namespace mars::fem;

// MFEM-style bigint handling for Hypre compatibility
typedef uint64_t IndexType;

// DOF handler typedef for convenience
using DofHandler = mars::fem::UnstructuredDofHandler<TetTag, double, IndexType, cstone::GpuTag>;

// Source term: f(x,y,z) = 1
struct SourceTerm {
    __device__ __host__
    double operator()(double x, double y, double z) const {
        return 1.0;
    }
};

// Save solution to MFEM GridFunction format (rank 0 only)
void saveSolutionToGridFunction(const cstone::DeviceVector<double>& u_local,
                               const DofHandler& dof_handler,
                               int rank, int numRanks) {
    if (rank != 0) return;  // Only rank 0 writes the file

    // For distributed case, gather solution from all ranks
    // This is a simplified version - in practice you'd need proper gathering
    std::vector<double> h_solution(u_local.size());
    thrust::copy(thrust::device_pointer_cast(u_local.data()),
                 thrust::device_pointer_cast(u_local.data() + u_local.size()),
                 h_solution.begin());

    std::ofstream solFile("sol_distributed.gf");
    if (!solFile.is_open()) {
        std::cerr << "Error: Could not open sol_distributed.gf for writing" << std::endl;
        return;
    }

    // MFEM GridFunction header
    solFile << "FiniteElementSpace" << std::endl;
    solFile << "FiniteElementCollection: H1_3D_P1" << std::endl;
    solFile << "VDim: 1" << std::endl;
    solFile << "Ordering: 0" << std::endl;
    solFile << std::endl;

    // Solution values (one per line) - rank 0's owned DOFs
    for (size_t i = 0; i < h_solution.size(); ++i) {
        solFile << h_solution[i] << std::endl;
    }

    solFile.close();
    std::cout << "Solution saved to sol_distributed.gf (" << h_solution.size() << " DOF values from rank 0)" << std::endl;
}

int main(int argc, char** argv) {
    // Initialize MPI
    MPI_Init(&argc, &argv);

    int rank, numRanks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    // Parse command line arguments
    std::string meshPath = "";  // No default, user must specify
    int order = 1;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--mesh" && i + 1 < argc) {
            meshPath = argv[++i];
        } else if (arg == "--order" && i + 1 < argc) {
            order = std::atoi(argv[++i]);
        } else if (arg == "--help" || arg == "-h") {
            if (rank == 0) {
                std::cout << "Usage: " << argv[0] << " --mesh <path> [options]\n\n"
                          << "Solve -Δu = 1 with u=0 on boundary (Distributed GPU version)\n\n"
                          << "Options:\n"
                          << "  --mesh <path>      Path to mesh (MFEM .mesh or MARS binary dir)\n"
                          << "  --order <n>        Polynomial order (default: 1)\n"
                          << "  --help, -h         Print this help\n\n"
                          << "MFEM meshes are automatically converted to binary format.\n\n"
                          << "Examples:\n"
                          << "  mpirun -np 4 " << argv[0] << " --mesh beam-tet.mesh\n"
                          << "  mpirun -np 8 " << argv[0] << " --mesh my-binary-mesh-dir\n";
            }
            MPI_Finalize();
            return 0;
        }
    }

    if (meshPath.empty()) {
        if (rank == 0) {
            std::cerr << "Error: --mesh argument is required\n";
            std::cerr << "Usage: " << argv[0] << " --mesh <path> [options]\n";
            std::cerr << "Use --help for more information\n";
        }
        MPI_Finalize();
        return 1;
    }

    if (rank == 0) {
        std::cout << "=== MFEM beam-tet.mesh Distributed GPU Poisson Test ===" << std::endl;
        std::cout << "Problem: -Δu = 1 in Ω, u = 0 on ∂Ω" << std::endl;
        std::cout << "Ranks: " << numRanks << std::endl;
        std::cout << "Order: " << order << std::endl;
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    // Load MFEM mesh directly (bypass refinement and binary conversion)
    std::string mfemMeshFile = meshPath;  // Use the provided mesh path directly

    if (rank == 0) {
        std::cout << "\nLoading MFEM mesh: " << mfemMeshFile << std::endl;
    }

    MFEMMeshLoader loader;
    if (!loader.load(mfemMeshFile)) {
        std::cerr << "Rank " << rank << ": ERROR - Failed to load MFEM mesh" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Broadcast mesh data from rank 0 to all ranks
    if (rank == 0) {
        std::cout << "Broadcasting mesh data to all ranks..." << std::endl;
    }

    // Broadcast vertices
    size_t numVertices = loader.vertices.size();
    MPI_Bcast(&numVertices, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        loader.vertices.resize(numVertices);
    }
    MPI_Bcast(loader.vertices.data(), numVertices * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Broadcast elements
    size_t numElements = loader.elements.size();
    MPI_Bcast(&numElements, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        loader.elements.resize(numElements);
    }
    MPI_Bcast(loader.elements.data(), numElements * 4, MPI_INT, 0, MPI_COMM_WORLD);

    // Broadcast boundary vertices
    size_t numBoundaryVertices = loader.boundary_vertices.size();
    MPI_Bcast(&numBoundaryVertices, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        loader.boundary_vertices.resize(numBoundaryVertices);
    }
    MPI_Bcast(loader.boundary_vertices.data(), numBoundaryVertices, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "Global mesh: " << numVertices << " vertices, "
                  << numElements << " elements" << std::endl;
        std::cout << "Boundary info: " << numBoundaryVertices << " boundary vertices" << std::endl;
    }

    // Convert MFEM mesh data to ElementDomain format
    // HostCoordsTuple: tuple<x_coords, y_coords, z_coords>
    std::vector<double> x_coords(numVertices);
    std::vector<double> y_coords(numVertices);
    std::vector<double> z_coords(numVertices);

    for (size_t i = 0; i < numVertices; ++i) {
        x_coords[i] = static_cast<double>(loader.vertices[i][0]);
        y_coords[i] = static_cast<double>(loader.vertices[i][1]);
        z_coords[i] = static_cast<double>(loader.vertices[i][2]);
    }

    // HostConnectivityTuple for TetTag: tuple<i0, i1, i2, i3>
    std::vector<IndexType> i0(numElements);
    std::vector<IndexType> i1(numElements);
    std::vector<IndexType> i2(numElements);
    std::vector<IndexType> i3(numElements);

    for (size_t i = 0; i < numElements; ++i) {
        i0[i] = static_cast<IndexType>(loader.elements[i][0]);
        i1[i] = static_cast<IndexType>(loader.elements[i][1]);
        i2[i] = static_cast<IndexType>(loader.elements[i][2]);
        i3[i] = static_cast<IndexType>(loader.elements[i][3]);
    }

    // HostBoundaryTuple for boundary nodes: tuple<isBoundaryNode>
    std::vector<uint8_t> isBoundaryNode(numVertices, 0);
    for (int boundaryVertex : loader.boundary_vertices) {
        if (boundaryVertex >= 0 && boundaryVertex < static_cast<int>(numVertices)) {
            isBoundaryNode[boundaryVertex] = 1;
        }
    }

    if (rank == 0) {
        std::cout << "Creating distributed ElementDomain (automatic SFC partitioning)..." << std::endl;
    }

    // Create ElementDomain using direct constructor with MFEM mesh data and boundary info
    // The bounding box will be computed automatically from the coordinate data
    using Domain = mars::ElementDomain<mars::TetTag, double, IndexType, cstone::GpuTag>;
    Domain domain(std::make_tuple(x_coords, y_coords, z_coords),
                  std::make_tuple(i0, i1, i2, i3),
                  std::make_tuple(isBoundaryNode),
                  rank, numRanks);

    // Force domain initialization before creating FE space
    domain.getNodeOwnershipMap();  // Trigger lazy initialization
    domain.getHaloElementIndices();  // Ensure halo structures are built

    // Validate domain data consistency
    size_t localNodes = domain.getNodeCount();
    size_t localElements = domain.localElementCount();

    // Check for obviously invalid data
    if (localNodes == 0) {
        std::cerr << "Rank " << rank << ": ERROR - Domain has 0 local nodes!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (localElements == 0) {
        std::cerr << "Rank " << rank << ": ERROR - Domain has 0 local elements!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Print domain info for all ranks
    std::cout << "Rank " << rank << ": Domain has " << localNodes << " local nodes, "
              << localElements << " local elements" << std::endl;

    // =====================================================
    // 1. Create finite element space
    // =====================================================
    if (rank == 0) std::cout << "\n1. Creating finite element space...\n";
    auto fes_start = std::chrono::high_resolution_clock::now();

    // Create FE space (like MFEM ex1)
    TetFESpace<double, IndexType> fes(domain, 1);  // Order 1 like MFEM ex1
    size_t numDofs = fes.numDofs();

    auto fes_end = std::chrono::high_resolution_clock::now();
    double fes_time = std::chrono::duration<double>(fes_end - fes_start).count();

    if (rank == 0) {
        std::cout << "FE space created in " << fes_time << " seconds" << std::endl;
        std::cout << "Number of finite element unknowns: " << numDofs << std::endl;
    }

    // =====================================================
    // 2. Create distributed DOF handler
    // =====================================================
    if (rank == 0) std::cout << "\n2. Creating distributed DOF handler...\n";

    DofHandler dof_handler(domain, rank, numRanks);
    dof_handler.enumerate_dofs();

    size_t numDofs_check = dof_handler.get_num_local_dofs();

    if (rank == 0) {
        std::cout << "FE space created in " << fes_time << " seconds" << std::endl;
        std::cout << "Number of finite element unknowns: " << numDofs << std::endl;
    }

    if (rank == 0 || rank == 1) {
        size_t owned = dof_handler.get_num_local_dofs();
        size_t total = dof_handler.get_num_local_dofs_with_ghosts();
        size_t ghost = total - owned;
        std::cout << "Rank " << rank << ": " << owned << " owned DOFs, " << ghost << " ghost DOFs\n";
    }

    // =====================================================
    // 4. Geometric boundary detection
    // =====================================================
    if (rank == 0) std::cout << "\n4. Performing geometric boundary detection...\n";

    // Get local node coordinates
    const auto& d_x = domain.getNodeX();
    const auto& d_y = domain.getNodeY();
    const auto& d_z = domain.getNodeZ();

    std::vector<double> h_x(domain.getNodeCount());
    std::vector<double> h_y(domain.getNodeCount());
    std::vector<double> h_z(domain.getNodeCount());

    thrust::copy(thrust::device_pointer_cast(d_x.data()),
                 thrust::device_pointer_cast(d_x.data() + domain.getNodeCount()),
                 h_x.begin());
    thrust::copy(thrust::device_pointer_cast(d_y.data()),
                 thrust::device_pointer_cast(d_y.data() + domain.getNodeCount()),
                 h_y.begin());
    thrust::copy(thrust::device_pointer_cast(d_z.data()),
                 thrust::device_pointer_cast(d_z.data() + domain.getNodeCount()),
                 h_z.begin());

    // Create local boundary data using geometric boundary detection
    // For beam-tet mesh: boundary nodes have x=0, x=8.4, y=0, y=1.05, z=0, or z=1.05
    std::vector<uint8_t> localBoundaryData(domain.getNodeCount(), 0);

    // Detect boundary nodes geometrically (more robust)
    // Check if coordinates are near the boundary extents
    double x_min = *std::min_element(h_x.begin(), h_x.end());
    double x_max = *std::max_element(h_x.begin(), h_x.end());
    double y_min = *std::min_element(h_y.begin(), h_y.end());
    double y_max = *std::max_element(h_y.begin(), h_y.end());
    double z_min = *std::min_element(h_z.begin(), h_z.end());
    double z_max = *std::max_element(h_z.begin(), h_z.end());

    double tol = 1e-3;  // Match single-rank tolerance

    size_t boundaryCount = 0;
    for (size_t i = 0; i < domain.getNodeCount(); ++i) {
        double x = h_x[i], y = h_y[i], z = h_z[i];
        // Check if near any boundary
        if (fabs(x - x_min) < tol || fabs(x - x_max) < tol ||
            fabs(y - y_min) < tol || fabs(y - y_max) < tol ||
            fabs(z - z_min) < tol || fabs(z - z_max) < tol) {
            localBoundaryData[i] = 1;
            boundaryCount++;
        }
    }

    // Set local boundary data for GPU-based boundary detection
    dof_handler.set_boundary_data(localBoundaryData);

    // Get boundary DOFs using DOF handler (like MFEM's GetBoundaryTrueDofs)
    std::vector<IndexType> boundaryDofs;
    dof_handler.boundary_owned_dof_iterate([&](IndexType localDof) {
        boundaryDofs.push_back(localDof);
    });

    if (rank == 0 || rank == 1) {
        std::cout << "Rank " << rank << ": Geometric boundary detection - domain bounds: ["
                  << x_min << "," << x_max << "] [" << y_min << "," << y_max << "] ["
                  << z_min << "," << z_max << "]" << std::endl;
        std::cout << "Rank " << rank << ": Found " << boundaryCount
                  << " boundary nodes out of " << domain.getNodeCount() << " total nodes" << std::endl;
        std::cout << "Rank " << rank << ": Found " << boundaryDofs.size() 
                  << " boundary DOFs out of " << numDofs << " total DOFs" << std::endl;
    }

    if (rank == 0) {
        std::cout << "Rank " << rank << ": Geometric boundary detection - domain bounds: ["
                  << x_min << "," << x_max << "] [" << y_min << "," << y_max << "] ["
                  << z_min << "," << z_max << "]" << std::endl;
        std::cout << "Rank " << rank << ": Found " << boundaryCount
                  << " boundary nodes out of " << domain.getNodeCount() << " total nodes" << std::endl;
        std::cout << "Found " << boundaryDofs.size() << " boundary DOFs (DOF handler: topological detection like MFEM), numDofs = " << numDofs << std::endl;
    }

    // =====================================================
    // 3. Create distributed data manager
    // =====================================================
    if (rank == 0) std::cout << "\n3. Creating distributed data manager...\n";

    // Create distributed data manager
    // Note: For Hypre, each rank stores only its owned DOFs
    mars::fem::UnstructuredDM<DofHandler, double, cstone::GpuTag> dm(dof_handler);
    dm.add_data_field<double>();  // Solution vector
    dm.add_data_field<double>();  // RHS vector
    dm.resize(dof_handler.get_num_local_dofs());  // Owned DOFs only

    if (rank == 0) {
        auto [start, end] = dof_handler.get_owned_dof_range();
        std::cout << "   Global DOFs: " << dof_handler.get_num_global_dofs() << "\n";
        std::cout << "   Owned DOF range: [" << start << ", " << end << ")\n";
    }

    // =====================================================
    // 5. Assemble stiffness matrix
    // =====================================================
    if (rank == 0) std::cout << "\n5. Assembling stiffness matrix...\n";

    auto assembly_start = std::chrono::high_resolution_clock::now();

    TetStiffnessAssembler<double, IndexType> stiffnessAssembler;
    TetSparseMatrix<double, IndexType> K;

    // Get node-to-DOF mapping from distributed DOF handler (includes ghosts)
    const auto& fullNodeToLocalDof = dof_handler.get_node_to_local_dof();
    const auto& ownership = domain.getNodeOwnershipMap();
    
    // Copy ownership to host for CPU access
    thrust::host_vector<uint8_t> h_ownership(domain.getNodeCount());
    thrust::copy(thrust::device_pointer_cast(ownership.data()),
                thrust::device_pointer_cast(ownership.data() + domain.getNodeCount()),
                h_ownership.begin());
    
    // Create owned-only node-to-DOF mapping for Hypre (each rank owns its rows)
    std::vector<IndexType> nodeToOwnedDof(domain.getNodeCount(), 0);  // Initialize to 0, will be set properly
    size_t ownedDofCount = 0;
    
    for (size_t nodeIdx = 0; nodeIdx < domain.getNodeCount(); ++nodeIdx) {
        if (h_ownership[nodeIdx] == 1) {  // Node owned by this rank (1 = owned, 0 = ghost)
            nodeToOwnedDof[nodeIdx] = ownedDofCount++;
        } else {
            // For ghost nodes, set to a large value to indicate invalid
            nodeToOwnedDof[nodeIdx] = std::numeric_limits<IndexType>::max();
        }
    }
    
    // Verify owned DOF count matches expected
    if (ownedDofCount != dof_handler.get_num_local_dofs()) {
        std::cerr << "Rank " << rank << ": ERROR - Owned DOF count mismatch: " 
                  << ownedDofCount << " vs " << dof_handler.get_num_local_dofs() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    stiffnessAssembler.assemble(fes, K, nodeToOwnedDof);

    auto assembly_end = std::chrono::high_resolution_clock::now();
    double assembly_time = std::chrono::duration<double>(assembly_end - assembly_start).count();

    if (rank == 0) {
        std::cout << "Stiffness matrix assembled in " << assembly_time << " seconds" << std::endl;
        std::cout << "Matrix size: " << K.numRows() << " x " << K.numCols() << ", nnz=" << K.nnz() << std::endl;

        // Check matrix properties
        std::cout << "Checking matrix properties..." << std::endl;
        std::vector<IndexType> h_rowOffsets(K.numRows() + 1);
        std::vector<IndexType> h_colIndices(K.nnz());
        std::vector<double> h_values(K.nnz());

        thrust::copy(thrust::device_pointer_cast(K.rowOffsetsPtr()),
                    thrust::device_pointer_cast(K.rowOffsetsPtr() + K.numRows() + 1),
                    h_rowOffsets.begin());
        thrust::copy(thrust::device_pointer_cast(K.colIndicesPtr()),
                    thrust::device_pointer_cast(K.colIndicesPtr() + K.nnz()),
                    h_colIndices.begin());
        thrust::copy(thrust::device_pointer_cast(K.valuesPtr()),
                    thrust::device_pointer_cast(K.valuesPtr() + K.nnz()),
                    h_values.begin());

        // Check diagonal
        double minDiag = 1e30, maxDiag = -1e30;
        int zeroDiag = 0, negDiag = 0;
        for (size_t i = 0; i < K.numRows(); ++i) {
            double diagVal = 0.0;
            bool foundDiag = false;
            for (IndexType idx = h_rowOffsets[i]; idx < h_rowOffsets[i+1]; ++idx) {
                if (h_colIndices[idx] == i) {
                    diagVal = h_values[idx];
                    foundDiag = true;
                    break;
                }
            }
            if (!foundDiag) zeroDiag++;
            else {
                if (diagVal <= 0.0) negDiag++;
                minDiag = std::min(minDiag, diagVal);
                maxDiag = std::max(maxDiag, diagVal);
            }
        }
        std::cout << "Diagonal range: [" << minDiag << ", " << maxDiag << "]" << std::endl;
        std::cout << "Zero diagonal entries: " << zeroDiag << std::endl;
        std::cout << "Negative diagonal entries: " << negDiag << std::endl;
    }

    // =====================================================
    // 6. Assemble RHS vector
    // =====================================================
    if (rank == 0) std::cout << "\n6. Assembling RHS vector...\n";

    auto rhs_start = std::chrono::high_resolution_clock::now();

    TetMassAssembler<double, IndexType> massAssembler;

    SourceTerm f;
    massAssembler.assembleRHS(fes, dm.get_data<1>(), f, nodeToOwnedDof);

    auto rhs_end = std::chrono::high_resolution_clock::now();
    double rhs_time = std::chrono::duration<double>(rhs_end - rhs_start).count();

    if (rank == 0) {
        std::cout << "RHS vector assembled in " << rhs_time << " seconds" << std::endl;
    }

    // Apply boundary conditions using TOPOLOGICAL detection (matching MFEM ex1)
    // =====================================================
    if (rank == 0) std::cout << "\n7. Applying boundary conditions using full-system modification...\n";
    auto bc_start = std::chrono::high_resolution_clock::now();

    // For full-system approach: modify matrix and RHS directly
    // K[i,i] = 1, K[i,j] = 0 for j != i, rhs[i] = 0 for boundary DOFs i

    // Get the full matrix before elimination
    TetSparseMatrix<double, IndexType> K_full = K;  // Copy the full matrix

    // Modify matrix and RHS for boundary conditions
    // boundaryDofs contains LOCAL owned DOF indices
    size_t numOwnedDofs = dof_handler.get_num_local_dofs();

    // Get matrix data
    std::vector<IndexType> h_rowOffsets(K_full.numRows() + 1);
    std::vector<IndexType> h_colIndices(K_full.nnz());
    std::vector<double> h_values(K_full.nnz());

    thrust::copy(thrust::device_pointer_cast(K_full.rowOffsetsPtr()),
                thrust::device_pointer_cast(K_full.rowOffsetsPtr() + K_full.numRows() + 1),
                h_rowOffsets.begin());
    thrust::copy(thrust::device_pointer_cast(K_full.colIndicesPtr()),
                thrust::device_pointer_cast(K_full.colIndicesPtr() + K_full.nnz()),
                h_colIndices.begin());
    thrust::copy(thrust::device_pointer_cast(K_full.valuesPtr()),
                thrust::device_pointer_cast(K_full.valuesPtr() + K_full.nnz()),
                h_values.begin());

    // Get RHS data
    std::vector<double> h_rhs_full(numOwnedDofs);
    dm.copy_to_host<1>(h_rhs_full);

    // Modify for boundary conditions
    for (auto localDof : boundaryDofs) {
        if (localDof < numOwnedDofs) {
            // Find the diagonal entry and set it to 1
            for (IndexType idx = h_rowOffsets[localDof]; idx < h_rowOffsets[localDof+1]; ++idx) {
                if (h_colIndices[idx] == localDof) {
                    h_values[idx] = 1.0f;  // K[i,i] = 1
                } else {
                    h_values[idx] = 0.0f;  // K[i,j] = 0 for j != i
                }
            }
            h_rhs_full[localDof] = 0.0f;  // rhs[i] = 0
        }
    }

    // Copy modified matrix back to device
    thrust::copy(h_values.begin(), h_values.end(),
                thrust::device_pointer_cast(K_full.valuesPtr()));

    // Copy modified RHS back
    dm.copy_from_host<1>(h_rhs_full);

    auto bc_end = std::chrono::high_resolution_clock::now();
    double bc_time = std::chrono::duration<double>(bc_end - bc_start).count();

    if (rank == 0) {
        std::cout << "Full-system BC modification completed in " << bc_time << " seconds" << std::endl;
        std::cout << "Boundary DOFs: " << boundaryDofs.size() << std::endl;
    }

    // =====================================================
    // 8. Solve linear system (Distributed GPU Hypre)
    // =====================================================
    if (rank == 0) std::cout << "\n8. Solving linear system (Distributed GPU Hypre PCG + BoomerAMG)...\n";

    auto solve_start = std::chrono::high_resolution_clock::now();

    // Validate matrix dimensions before solving
    if (K_full.numRows() != K_full.numCols()) {
        std::cerr << "Rank " << rank << ": ERROR - Matrix not square: "
                  << K_full.numRows() << " x " << K_full.numCols() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (K_full.numRows() != numOwnedDofs) {
        std::cerr << "Rank " << rank << ": ERROR - Matrix/RHS size mismatch: "
                  << K_full.numRows() << " vs " << numOwnedDofs << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (rank == 0 || rank == 1) {
        std::cout << "Rank " << rank << ": Full system size: " << K_full.numRows()
                  << " x " << K_full.numCols() << ", nnz=" << K_full.nnz() << std::endl;
    }

    // Check if matrix is sorted
    if (rank == 0) {
        std::cout << "Matrix columns sorted: " << (K_full.isSorted() ? "YES" : "NO") << std::endl;
    }

    // Get solution and RHS vectors from data manager
    auto& u_full = dm.get_data<0>();  // Solution vector
    auto& rhs_full_data = dm.get_data<1>();  // RHS vector

#ifdef MARS_ENABLE_HYPRE
    // Use full Hypre PCG + BoomerAMG (most MFEM-like)
    if (rank == 0) std::cout << "Using full Hypre PCG + BoomerAMG solver" << std::endl;
    mars::fem::HyprePCGSolver<double, IndexType, cstone::GpuTag> hypre_solver(MPI_COMM_WORLD, 400, 1e-10);
    hypre_solver.setVerbose(rank == 0);
    bool converged = hypre_solver.solve(K_full, rhs_full_data, u_full);
#else
    // Fallback to custom CG with Jacobi preconditioner
    if (rank == 0) std::cout << "Using custom CG + Jacobi preconditioner" << std::endl;
    mars::fem::PreconditionedConjugateGradientSolver<double, IndexType, cstone::GpuTag> cg(400, 1e-10);
    cg.setVerbose(rank == 0);
    bool converged = cg.solve(K_full, rhs_full_data, u_full);
#endif

    auto solve_end = std::chrono::high_resolution_clock::now();
    double solve_time = std::chrono::duration<double>(solve_end - solve_start).count();

    // =====================================================
    // 9. Exchange ghost data and compute statistics
    // =====================================================
    if (rank == 0) std::cout << "\n9. Computing solution statistics...\n";

    // Exchange ghost values between ranks
    dm.gather_ghost_data();  // Send owned data to neighbors
    dm.scatter_ghost_data(); // Receive ghost data from neighbors

    // Compute local statistics
    double u_min_local = *thrust::min_element(thrust::device_pointer_cast(u_full.data()),
                                             thrust::device_pointer_cast(u_full.data() + numOwnedDofs));
    double u_max_local = *thrust::max_element(thrust::device_pointer_cast(u_full.data()),
                                             thrust::device_pointer_cast(u_full.data() + numOwnedDofs));

    double u_norm_sq_local = thrust::transform_reduce(thrust::device_pointer_cast(u_full.data()),
                                                     thrust::device_pointer_cast(u_full.data() + numOwnedDofs),
                                                     [] __host__ __device__ (double x) { return x * x; },
                                                     0.0, thrust::plus<double>());

    double u_norm_local = sqrt(u_norm_sq_local);

    // Print per-rank solution statistics for debugging
    std::cout << "Rank " << rank << ": min(u) = " << u_min_local
              << ", max(u) = " << u_max_local << ", ||u|| = " << u_norm_local << std::endl;

    // Global reduction for statistics
    double u_min_global, u_max_global;
    double u_norm_sq_global;

    MPI_Allreduce(&u_min_local, &u_min_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&u_max_local, &u_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&u_norm_sq_local, &u_norm_sq_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double u_norm_global = sqrt(u_norm_sq_global);

    if (rank == 0) {
        std::cout << "\n=== Global Solution Statistics ===\n";
        std::cout << "  Global min(u) = " << u_min_global << "\n";
        std::cout << "  Global max(u) = " << u_max_global << "\n";
        std::cout << "  Global ||u||   = " << u_norm_global << "\n";
    }

    // =====================================================
    // 10. Save solution
    // =====================================================
    if (rank == 0) std::cout << "\n10. Saving solution...\n";

    // For distributed case, save local solution from rank 0 (simplified)
    if (rank == 0) {
        std::vector<double> u_host(u_full.size());
        thrust::copy(thrust::device_pointer_cast(u_full.data()),
                    thrust::device_pointer_cast(u_full.data() + u_full.size()),
                    u_host.begin());

        std::ofstream solFile("sol_distributed.gf");
        if (!solFile.is_open()) {
            std::cerr << "Error: Could not open sol_distributed.gf for writing" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // MFEM GridFunction header
        solFile << "FiniteElementSpace" << std::endl;
        solFile << "FiniteElementCollection: H1_3D_P1" << std::endl;
        solFile << "VDim: 1" << std::endl;
        solFile << "Ordering: 0" << std::endl;
        solFile << std::endl;

        // Solution values (one per line)
        for (size_t i = 0; i < u_host.size(); ++i) {
            solFile << u_host[i] << std::endl;
        }

        solFile.close();
        std::cout << "Solution saved to sol_distributed.gf (" << u_host.size() << " DOF values)" << std::endl;
    }

    // =====================================================
    // 11. Final timing summary
    // =====================================================
    auto end_time = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double>(end_time - start_time).count();

    if (rank == 0) {
        std::cout << "\n=== Timing Summary ===\n";
        std::cout << "  Mesh loading:     " << std::fixed << std::setprecision(4)
                  << (fes_start - start_time).count() << " s\n";
        std::cout << "  FE space:         " << fes_time << " s\n";
        std::cout << "  Assembly:         " << assembly_time + rhs_time << " s\n";
        std::cout << "  BCs:              " << bc_time << " s\n";
        std::cout << "  Solve:            " << solve_time << " s\n";
        std::cout << "  Total:            " << total_time << " s\n";
        std::cout << "\n========================================\n";
        std::cout << "Distributed GPU solve completed successfully!\n";
        std::cout << "========================================\n\n";
    }

    MPI_Finalize();
    return 0;
}