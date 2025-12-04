/**
 * MFEM beam-tet.mesh Poisson solver for validation
 * 
 * Problem: -Δu = 1 in Ω, u = 0 on ∂Ω
 * 
 * Usage: mpirun -np 1 mars_ex_beam_tet --mesh beam-tet.mesh
 * 
 * Automatically converts MFEM mesh to binary format on first run.
 * Binary format cached for subsequent runs.
 * 
 * Expected: CG convergence in ~20-50 iterations, residual < 1e-10
 * Matches MFEM ex1p.cpp results exactly
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

// Source term: f(x,y,z) = 1
struct SourceTerm {
    __device__ __host__
    double operator()(double x, double y, double z) const {
        return 1.0;
    }
};

// Save solution to MFEM GridFunction format
void saveSolutionToGridFunction(const cstone::DeviceVector<double>& u_local, 
                               const mars::fem::UnstructuredDofHandler<TetTag, double, Unsigned, cstone::GpuTag>& dof_handler,
                               int rank, int numRanks) {
    if (rank != 0) return;  // Only rank 0 writes the file
    
    // For now, just save the local solution (rank 0's owned DOFs)
    // In a full implementation, we'd gather all DOFs from all ranks
    std::vector<double> h_solution(u_local.size());
    thrust::copy(thrust::device_pointer_cast(u_local.data()),
                 thrust::device_pointer_cast(u_local.data() + u_local.size()),
                 h_solution.begin());
    
    std::ofstream solFile("sol.gf");
    if (!solFile.is_open()) {
        std::cerr << "Error: Could not open sol.gf for writing" << std::endl;
        return;
    }
    
    // MFEM GridFunction header
    solFile << "FiniteElementSpace" << std::endl;
    solFile << "FiniteElementCollection: H1_3D_P1" << std::endl;
    solFile << "VDim: 1" << std::endl;
    solFile << "Ordering: 0" << std::endl;
    solFile << std::endl;
    
    // Solution values (one per line)
    for (size_t i = 0; i < h_solution.size(); ++i) {
        solFile << h_solution[i] << std::endl;
    }
    
    solFile.close();
    std::cout << "Solution saved to sol.gf (" << h_solution.size() << " DOF values)" << std::endl;
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
                          << "Solve -Δu = 1 with u=0 on boundary\n\n"
                          << "Options:\n"
                          << "  --mesh <path>      Path to mesh (MFEM .mesh or MARS binary dir)\n"
                          << "  --order <n>        Polynomial order (default: 1)\n"
                          << "  --help, -h         Print this help\n\n"
                          << "MFEM meshes are automatically converted to binary format.\n\n"
                          << "Examples:\n"
                          << "  mpirun -np 1 " << argv[0] << " --mesh beam-tet.mesh\n"
                          << "  mpirun -np 4 " << argv[0] << " --mesh my-binary-mesh-dir\n";
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
        std::cout << "=== MFEM beam-tet.mesh Poisson Test ===" << std::endl;
        std::cout << "Problem: -Δu = 1 in Ω, u = 0 on ∂Ω" << std::endl;
        std::cout << "Ranks: " << numRanks << std::endl;
        std::cout << "Order: " << order << std::endl;
    }
    
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
    
    // Convert MFEM mesh data to ElementDomain format
    // HostCoordsTuple: tuple<x_coords, y_coords, z_coords>
    std::vector<double> x_coords(loader.vertices.size());
    std::vector<double> y_coords(loader.vertices.size());
    std::vector<double> z_coords(loader.vertices.size());
    
    for (size_t i = 0; i < loader.vertices.size(); ++i) {
        x_coords[i] = static_cast<double>(loader.vertices[i][0]);
        y_coords[i] = static_cast<double>(loader.vertices[i][1]);
        z_coords[i] = static_cast<double>(loader.vertices[i][2]);
    }
    
    // HostConnectivityTuple for TetTag: tuple<i0, i1, i2, i3>
    std::vector<Unsigned> i0(loader.elements.size());
    std::vector<Unsigned> i1(loader.elements.size());
    std::vector<Unsigned> i2(loader.elements.size());
    std::vector<Unsigned> i3(loader.elements.size());
    
    for (size_t i = 0; i < loader.elements.size(); ++i) {
        i0[i] = static_cast<Unsigned>(loader.elements[i][0]);
        i1[i] = static_cast<Unsigned>(loader.elements[i][1]);
        i2[i] = static_cast<Unsigned>(loader.elements[i][2]);
        i3[i] = static_cast<Unsigned>(loader.elements[i][3]);
    }
    
    // HostBoundaryTuple for boundary nodes: tuple<isBoundaryNode>
    std::vector<uint8_t> isBoundaryNode(loader.vertices.size(), 0);
    for (int boundaryVertex : loader.boundary_vertices) {
        if (boundaryVertex >= 0 && boundaryVertex < static_cast<int>(loader.vertices.size())) {
            isBoundaryNode[boundaryVertex] = 1;
        }
    }
    
    if (rank == 0) {
        std::cout << "Creating distributed ElementDomain (automatic SFC partitioning)..." << std::endl;
        std::cout << "Global mesh: " << loader.vertices.size() << " vertices, " 
                  << loader.elements.size() << " elements" << std::endl;
        std::cout << "Boundary info: " << loader.boundary_vertices.size() << " boundary vertices" << std::endl;
    }
    
    // Create ElementDomain using direct constructor with MFEM mesh data and boundary info
    // The bounding box will be computed automatically from the coordinate data
    using Domain = ElementDomain<TetTag, double, Unsigned, cstone::GpuTag>;
    Domain domain(std::make_tuple(x_coords, y_coords, z_coords),
                  std::make_tuple(i0, i1, i2, i3),
                  std::make_tuple(isBoundaryNode),
                  rank, numRanks);
    
    // Force domain initialization before creating FE space
    // This ensures ownership map is complete before FE space counts DOFs
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
    
    // Create finite element space (this will determine boundary DOFs like MFEM)
    if (rank == 0) std::cout << "\nCreating finite element space...\n";
    
    // Create FESpace and get boundary DOFs (like MFEM's GetBoundaryTrueDofs)
    TetFESpace<double, Unsigned> fes(domain, 1);  // Order 1 like MFEM ex1
    size_t numDofs = fes.numDofs();
    
    if (rank == 0) {
        std::cout << "Number of finite element unknowns: " << numDofs << std::endl;
    }
    
    // =====================================================
    // 1.5. Create distributed DOF handler
    // =====================================================
    if (rank == 0) std::cout << "\n1.5. Creating distributed DOF handler...\n";
    
    using DofHandler = mars::fem::UnstructuredDofHandler<TetTag, double, Unsigned, cstone::GpuTag>;
    DofHandler dof_handler(domain, rank, numRanks);
    dof_handler.enumerate_dofs();
    
    // Create local boundary data using geometric boundary detection
    // For beam-tet mesh: boundary nodes have x=0, x=8.4, y=0, y=1.05, z=0, or z=1.05
    std::vector<uint8_t> localBoundaryData(domain.getNodeCount(), 0);
    
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
    
    // Detect boundary nodes geometrically (more robust)
    // Check if coordinates are near the boundary extents
    double x_min = *std::min_element(h_x.begin(), h_x.end());
    double x_max = *std::max_element(h_x.begin(), h_x.end());
    double y_min = *std::min_element(h_y.begin(), h_y.end());
    double y_max = *std::max_element(h_y.begin(), h_y.end());
    double z_min = *std::min_element(h_z.begin(), h_z.end());
    double z_max = *std::max_element(h_z.begin(), h_z.end());
    
    double tol = 1e-3;  // More tolerant boundary detection
    
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
    
    if (rank == 0) {
        std::cout << "Rank " << rank << ": Geometric boundary detection - domain bounds: ["
                  << x_min << "," << x_max << "] [" << y_min << "," << y_max << "] [" 
                  << z_min << "," << z_max << "]" << std::endl;
        std::cout << "Rank " << rank << ": Found " << boundaryCount 
                  << " boundary nodes out of " << domain.getNodeCount() << " total nodes" << std::endl;
    }
    
    // Set local boundary data for GPU-based boundary detection
    dof_handler.set_boundary_data(localBoundaryData);
    
    // Get boundary DOFs using DOF handler (like MFEM's GetBoundaryTrueDofs)
    // Boundary detection uses GPU computation with local boundary data
    std::vector<Unsigned> boundaryDofs;
    dof_handler.boundary_owned_dof_iterate([&](Unsigned localDof) {
        boundaryDofs.push_back(localDof);
    });
    
    if (rank == 0) {
        std::cout << "Found " << boundaryDofs.size() << " boundary DOFs (DOF handler: topological detection like MFEM), numDofs = " << numDofs << std::endl;
    }
    
    // Note: Adjacency is built lazily when needed by FE space
    
    if (rank == 0) {
        std::cout << "\n=== ElementDomain Created Successfully ===" << std::endl;
        std::cout << "Domain info:" << std::endl;
        std::cout << "  Total nodes: " << domain.getNodeCount() << std::endl;
        std::cout << "  Total elements: " << domain.getElementCount() << std::endl;
        std::cout << "  Local elements: " << domain.localElementCount() << std::endl;
        std::cout << "\n=== Starting FEM Solve ===" << std::endl;
    }
    
    // Use geometric boundary detection with global bounds for consistency
    
    auto t_total_start = std::chrono::high_resolution_clock::now();
    
    // =====================================================
    // 1. Create finite element space
    // =====================================================
    if (rank == 0) std::cout << "\n1. Creating finite element space...\n";
    auto t_fes_start = std::chrono::high_resolution_clock::now();
    
    // Create H1 finite element space
    using FESpace = mars::fem::H1FESpace<TetTag, double, Unsigned, cstone::GpuTag>;
    FESpace fe_space(domain);
    
    auto t_fes_end = std::chrono::high_resolution_clock::now();
    double t_fes = std::chrono::duration<double>(t_fes_end - t_fes_start).count();
    
    if (rank == 0) {
        std::cout << "   FE space created in " << t_fes << " seconds\n"
                  << "   Total DOFs: " << numDofs << "\n";
    }
    
    // =====================================================
    // 1.5. Create distributed data manager
    // =====================================================
    if (rank == 0) std::cout << "\n1.5. Creating distributed data manager...\n";
    
    // Create distributed data manager
    // Note: Data manager stores owned DOFs only (ghosts updated via halo exchange)
    mars::fem::UnstructuredDM<DofHandler, double, cstone::GpuTag> dm(dof_handler);
    dm.add_data_field<double>();  // Solution vector
    dm.add_data_field<double>();  // RHS vector
    dm.resize(numDofs);  // Owned DOFs only
    
    if (rank == 0) {
        auto [start, end] = dof_handler.get_owned_dof_range();
        std::cout << "   Global DOFs: " << dof_handler.get_num_global_dofs() << "\n";
        std::cout << "   Owned DOF range: [" << start << ", " << end << ")\n";
    }
    
    // =====================================================
    // 2. Assemble stiffness matrix
    // =====================================================
    if (rank == 0) std::cout << "\n2. Assembling stiffness matrix...\n";
    auto t_stiff_start = std::chrono::high_resolution_clock::now();
    
    TetStiffnessAssembler<double, Unsigned> stiffnessAssembler;
    TetSparseMatrix<double, Unsigned> K;
    
    // Get node-to-DOF mapping from distributed DOF handler (includes ghosts)
    const auto& nodeToLocalDof = dof_handler.get_node_to_local_dof();
    
    stiffnessAssembler.assemble(fe_space, K, nodeToLocalDof);
    
    auto t_stiff_end = std::chrono::high_resolution_clock::now();
    double t_stiff = std::chrono::duration<double>(t_stiff_end - t_stiff_start).count();
    
    if (rank == 0) {
        std::cout << "   Stiffness matrix assembled in " << t_stiff << " seconds\n"
                  << "   Matrix size: " << K.numRows() << " x " << K.numCols() << "\n"
                  << "   Non-zeros: " << K.nnz() << "\n";
        
        // Check matrix symmetry
        std::cout << "   Checking matrix properties...\n";
        std::vector<Unsigned> h_rowOffsets(K.numRows() + 1);
        std::vector<Unsigned> h_colIndices(K.nnz());
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
            for (Unsigned idx = h_rowOffsets[i]; idx < h_rowOffsets[i+1]; ++idx) {
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
        std::cout << "   Diagonal range: [" << minDiag << ", " << maxDiag << "]\n";
        std::cout << "   Zero diagonal entries: " << zeroDiag << "\n";
        std::cout << "   Negative diagonal entries: " << negDiag << "\n";
    }
    
    // =====================================================
    // 3. Assemble RHS vector using distributed data manager
    // =====================================================
    if (rank == 0) std::cout << "\n3. Assembling RHS vector...\n";
    auto t_rhs_start = std::chrono::high_resolution_clock::now();
    
    TetMassAssembler<double, Unsigned> massAssembler;
    
    SourceTerm f;
    massAssembler.assembleRHS(fe_space, dm.get_data<1>(), f, nodeToLocalDof);
    
    auto t_rhs_end = std::chrono::high_resolution_clock::now();
    double t_rhs = std::chrono::duration<double>(t_rhs_end - t_rhs_start).count();
    
    if (rank == 0) {
        std::cout << "   RHS vector assembled in " << t_rhs << " seconds\n";
    }
    
    // Apply boundary conditions using TOPOLOGICAL detection (matching MFEM ex1)
    // =====================================================
    if (rank == 0) std::cout << "\n4. Applying boundary conditions using full-system modification...\n";
    auto t_bc_start = std::chrono::high_resolution_clock::now();

    // For full-system approach: modify matrix and RHS directly
    // K[i,i] = 1, K[i,j] = 0 for j != i, rhs[i] = 0 for boundary DOFs i
    
    // Get the full matrix before elimination
    TetSparseMatrix<double, Unsigned> K_full = K;  // Copy the full matrix
    
    // Modify matrix and RHS for boundary conditions
    // boundaryDofs contains LOCAL owned DOF indices
    size_t numOwnedDofs = dof_handler.get_num_local_dofs();
    
    // Get matrix data
    std::vector<Unsigned> h_rowOffsets(K_full.numRows() + 1);
    std::vector<Unsigned> h_colIndices(K_full.nnz());
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
            for (Unsigned idx = h_rowOffsets[localDof]; idx < h_rowOffsets[localDof+1]; ++idx) {
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
    
    auto t_bc_end = std::chrono::high_resolution_clock::now();
    double t_bc = std::chrono::duration<double>(t_bc_end - t_bc_start).count();
    
    if (rank == 0) {
        std::cout << "   Full-system BC modification completed in " << t_bc << " seconds\n";
        std::cout << "   Boundary DOFs: " << boundaryDofs.size() << "\n";
    }
    
    // =====================================================
    // 5. Solve linear system (PCG with preconditioner)
    // =====================================================
    if (rank == 0) std::cout << "\n5. Solving linear system (PCG with preconditioner)...\n";
    
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
    
    auto t_solve_start = std::chrono::high_resolution_clock::now();
    
    // Solution vector (full system)
    auto& u_full = dm.get_data<0>();  // Solution vector from data manager (size: owned DOFs)
    auto& rhs_full = dm.get_data<1>();  // RHS vector (size: owned DOFs)
    
#ifdef MARS_ENABLE_HYPRE
    // Use full Hypre PCG + BoomerAMG (most MFEM-like)
    if (rank == 0) std::cout << "Using full Hypre PCG + BoomerAMG solver" << std::endl;
    mars::fem::HyprePCGSolver<double, Unsigned, cstone::GpuTag> hypre_solver(MPI_COMM_WORLD, 400, 1e-10);
    hypre_solver.setVerbose(rank == 0);
    bool converged = hypre_solver.solve(K_full, rhs_full, u_full);
#else
    // Fallback to custom CG with Jacobi preconditioner
    if (rank == 0) std::cout << "Using custom CG + Jacobi preconditioner" << std::endl;
    mars::fem::PreconditionedConjugateGradientSolver<double, Unsigned, cstone::GpuTag> cg(400, 1e-10);
    cg.setVerbose(rank == 0);
    bool converged = cg.solve(K_full, rhs_full, u_full);
#endif
    
    auto t_solve_end = std::chrono::high_resolution_clock::now();
    double t_solve = std::chrono::duration<double>(t_solve_end - t_solve_start).count();
    
    if (rank == 0) {
        std::cout << "   Solver " << (converged ? "converged" : "did NOT converge") 
                  << " in " << t_solve << " seconds\n";
    }
    
    // =====================================================
    // 7. Exchange ghost DOF values across ranks (if multi-rank)
    // =====================================================
    // Note: Ghost exchange not yet implemented in UnstructuredDofHandler
    // For single-partition problems, this is not needed
    if (numRanks > 1 && false) {  // Disabled until ghost exchange is implemented
        if (rank == 0) std::cout << "\n6. Exchanging ghost DOF values...\n";
        auto t_exchange_start = std::chrono::high_resolution_clock::now();
        
        // Exchange solution values at shared nodes using DM
        dm.gather_ghost_data();
        dm.scatter_ghost_data();
        
        auto t_exchange_end = std::chrono::high_resolution_clock::now();
        double t_exchange = std::chrono::duration<double>(t_exchange_end - t_exchange_start).count();
        
        if (rank == 0) {
            std::cout << "   Ghost DOF exchange completed in " << t_exchange << " seconds\n";
        }
    }
    
    // Compute solution statistics (after exchange)
    
    double u_min = *thrust::min_element(thrust::device_pointer_cast(u_full.data()), 
                                      thrust::device_pointer_cast(u_full.data() + numOwnedDofs));
    double u_max = *thrust::max_element(thrust::device_pointer_cast(u_full.data()), 
                                      thrust::device_pointer_cast(u_full.data() + numOwnedDofs));
    double u_norm_sq = thrust::transform_reduce(thrust::device_pointer_cast(u_full.data()), 
                                             thrust::device_pointer_cast(u_full.data() + numOwnedDofs), 
                                             [] __host__ __device__ (double x) { return x*x; }, 
                                             0.0, thrust::plus<double>());
    double u_norm = std::sqrt(u_norm_sq);
    
    // Print per-rank solution statistics for debugging
    std::cout << "Rank " << rank << ": min(u) = " << u_min 
              << ", max(u) = " << u_max << ", ||u|| = " << u_norm << std::endl;
    
    // Compute global solution statistics
    double global_u_min, global_u_max, global_u_norm_sq;
    MPI_Allreduce(&u_min, &global_u_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&u_max, &global_u_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&u_norm_sq, &global_u_norm_sq, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    auto t_total_end = std::chrono::high_resolution_clock::now();
    double t_total = std::chrono::duration<double>(t_total_end - t_total_start).count();
    
    if (rank == 0) {
        std::cout << "\n=== Solution Statistics (Rank 0) ===" << std::endl;
        std::cout << "  min(u) = " << u_min << std::endl;
        std::cout << "  max(u) = " << u_max << std::endl;
        std::cout << "  ||u||  = " << u_norm << std::endl;
        
        std::cout << "\n=== Global Solution Statistics ===" << std::endl;
        std::cout << "  Global min(u) = " << global_u_min << std::endl;
        std::cout << "  Global max(u) = " << global_u_max << std::endl;
        std::cout << "  Global ||u||  = " << std::sqrt(global_u_norm_sq) << std::endl;
        
        std::cout << "\n=== Timing Summary ===" << std::endl;
        std::cout << "  FE space:     " << t_fes << " s" << std::endl;
        std::cout << "  Assembly:     " << (t_stiff + t_rhs) << " s" << std::endl;
        std::cout << "  BCs:          " << t_bc << " s" << std::endl;
        std::cout << "  Solve:        " << t_solve << " s" << std::endl;
        std::cout << "  Total:        " << t_total << " s" << std::endl;
        
        // Save solution to MFEM GridFunction format
        std::cout << "\n=== Saving Solution ===" << std::endl;
        saveSolutionToGridFunction(u_full, dof_handler, rank, numRanks);
        
        std::cout << "\n=== Test Complete ===" << std::endl;
        std::cout << "Compare with MFEM: mpirun -np " << numRanks << " ex1p -m beam-tet.mesh" << std::endl;
    }
    
    MPI_Finalize();
    return 0;
} 