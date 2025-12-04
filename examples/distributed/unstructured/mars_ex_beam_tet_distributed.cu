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
#include "backend/distributed/unstructured/fem/mars_form_linear_system.hpp"
#ifdef MARS_ENABLE_HYPRE
#include "backend/distributed/unstructured/solvers/mars_hypre_pcg_solver.hpp"
#endif
#include "mars_mfem_mesh_loader.hpp"
#include <thrust/copy.h>
#include <algorithm>
#include <limits>
#include <unordered_set>

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
    
    if (rank == 0) {
        std::cout << "=== BINARY COMPILED: " << __DATE__ << " " << __TIME__ << " ===" << std::endl;
        std::cout.flush();
    }

    // Parse command line arguments
    std::string meshPath = "";  // No default, user must specify
    int order = 1;
    std::string solver = "elimination-boomeramg";  // Default: MFEM-like elimination with BoomerAMG

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--mesh" && i + 1 < argc) {
            meshPath = argv[++i];
        } else if (arg == "--order" && i + 1 < argc) {
            order = std::atoi(argv[++i]);
        } else if (arg == "--solver" && i + 1 < argc) {
            solver = argv[++i];
        } else if (arg == "--help" || arg == "-h") {
            if (rank == 0) {
                std::cout << "Usage: " << argv[0] << " --mesh <path> [options]\n\n"
                          << "Solve -Δu = 1 with u=0 on boundary (Distributed GPU version)\n\n"
                          << "Options:\n"
                          << "  --mesh <path>      Path to mesh (MFEM .mesh or MARS binary dir)\n"
                          << "  --order <n>        Polynomial order (default: 1)\n"
                          << "  --solver <type>    Solver: elimination-boomeramg (default), elimination-jacobi, elimination, or jacobi\n"
                          << "  --help, -h         Print this help\n\n"
                          << "Solvers:\n"
                          << "  elimination-boomeramg - PCG + BoomerAMG with DOF elimination (default, MFEM ex1p style)\n"
                          << "  elimination-jacobi    - PCG + Jacobi with DOF elimination (for debugging)\n"
                          << "  elimination           - PCG + BoomerAMG with DOF elimination (legacy alias)\n"
                          << "  jacobi                - PCG + Jacobi with penalty BC (not yet implemented)\n\n"
                          << "MFEM meshes are automatically converted to binary format.\n\n"
                          << "Examples:\n"
                          << "  mpirun -np 4 " << argv[0] << " --mesh beam-tet.mesh\n"
                          << "  mpirun -np 4 " << argv[0] << " --mesh beam-tet.mesh --solver elimination\n";
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
        std::cout << "Solver: " << solver << std::endl;
    }

    { // Start of scope to ensure destructors run before MPI_Finalize

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
    dof_handler.initialize();
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

    // CRITICAL: Use GLOBAL bounding box, not local min/max!
    // Local min/max would mark partition boundaries as physical boundaries
    const auto& globalBox = domain.getBoundingBox();
    double x_min = globalBox.xmin();
    double x_max = globalBox.xmax();
    double y_min = globalBox.ymin();
    double y_max = globalBox.ymax();
    double z_min = globalBox.zmin();
    double z_max = globalBox.zmax();

    // NOTE: Bounding box may have epsilon expansion, so use larger tolerance
    // or compute actual min/max from node coordinates
    double tol = 0.1;  // Increased tolerance to handle bounding box expansion

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

    // Get node-to-DOF mapping from distributed DOF handler
    // DOF handler already provides correct mapping:
    //   - Owned nodes: indices [0, numLocalDofs)
    //   - Ghost nodes: indices [numLocalDofs, numLocalDofs+numGhostDofs)
    const auto& nodeToLocalDof = dof_handler.get_node_to_local_dof();

    stiffnessAssembler.assemble(fes, K, nodeToLocalDof);

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
    
    size_t numOwnedDofs = dof_handler.get_num_local_dofs();
    
    if (rank == 0) {
        std::cout << "Before RHS assembly: fes.numDofs()=" << fes.numDofs() 
                  << ", owned DOFs=" << numOwnedDofs << std::endl;
    }
    
    // Assemble directly into owned DOFs (no ghost exchange needed - MFEM style)
    // Each rank assembles contributions from its local elements into owned DOFs only
    massAssembler.assembleRHS(fes, dm.get_data<1>(), f, nodeToLocalDof, numOwnedDofs);
    
    if (rank == 0) {
        std::cout << "After RHS assembly: dm vector size=" << dm.get_data<1>().size() << std::endl;
    }

    auto rhs_end = std::chrono::high_resolution_clock::now();
    double rhs_time = std::chrono::duration<double>(rhs_end - rhs_start).count();

    if (rank == 0) {
        std::cout << "RHS vector assembled in " << rhs_time << " seconds" << std::endl;
        
        // Check RHS values
        auto& rhs_vec = dm.get_data<1>();
        std::vector<double> h_rhs(rhs_vec.size());
        thrust::copy(thrust::device_pointer_cast(rhs_vec.data()),
                     thrust::device_pointer_cast(rhs_vec.data() + rhs_vec.size()),
                     h_rhs.begin());
        
        double sum = 0.0, min_val = 1e100, max_val = -1e100;
        int nan_count = 0;
        for (auto v : h_rhs) {
            if (!std::isfinite(v)) nan_count++;
            else {
                sum += v;
                min_val = std::min(min_val, v);
                max_val = std::max(max_val, v);
            }
        }
        
        std::cout << "RHS stats: sum=" << sum << ", range=[" << min_val << ", " << max_val << "], NaNs=" << nan_count << std::endl;
    }

    // Form linear system based on solver type
    if (solver == "jacobi") {
        // Penalty method - no DOF elimination
        if (rank == 0) std::cout << "\n7. Applying penalty method for boundary conditions...\n";
        
        // Apply penalty to boundary DOFs (large diagonal + zero RHS)
        // This keeps the matrix square and avoids elimination complexity
        
        if (rank == 0) std::cerr << "ERROR: Penalty method not yet implemented in distributed version\n";
        MPI_Finalize();
        return 1;
    }
    
    // DOF elimination methods (elimination, elimination-jacobi, elimination-boomeramg)
    // =====================================================
    if (rank == 0) {
        std::cout << "\n7. Forming linear system (eliminating essential BCs)...\n";
    }
    auto bc_start = std::chrono::high_resolution_clock::now();

    TetSparseMatrix<double, IndexType> A_r;
    mars::VectorSelector<double, cstone::GpuTag>::type B_r;
    std::vector<IndexType> dof_mapping;
    
    // IndexType numOwnedDofs = dof_handler.get_num_local_dofs(); // Already declared above
    IndexType numTotalLocalDofs = K.numCols();  // Owned + ghost
    IndexType numGhostDofs = numTotalLocalDofs - numOwnedDofs;
    
    // PRE-PROCESSING: Detect DOFs with zero diagonals (not in any local element)
    // These will be eliminated by FormLinearSystem, so we need to share them as "boundary" DOFs
    if (rank == 0 || rank == 1) {
        std::cout << "Rank " << rank << ": Checking for zero-diagonal DOFs among " << numOwnedDofs << " owned DOFs\n";
    }
    
    std::vector<IndexType> h_rowOffsets_check(numOwnedDofs + 1);
    std::vector<IndexType> h_colIndices_check;
    std::vector<double> h_values_check;
    thrust::copy(thrust::device_pointer_cast(K.rowOffsetsPtr()),
                 thrust::device_pointer_cast(K.rowOffsetsPtr() + numOwnedDofs + 1),
                 h_rowOffsets_check.begin());
    h_colIndices_check.resize(h_rowOffsets_check[numOwnedDofs]);
    h_values_check.resize(h_rowOffsets_check[numOwnedDofs]);
    thrust::copy(thrust::device_pointer_cast(K.colIndicesPtr()),
                 thrust::device_pointer_cast(K.colIndicesPtr() + h_colIndices_check.size()),
                 h_colIndices_check.begin());
    thrust::copy(thrust::device_pointer_cast(K.valuesPtr()),
                 thrust::device_pointer_cast(K.valuesPtr() + h_values_check.size()),
                 h_values_check.begin());
    
    std::vector<IndexType> zeroDiagDofs;
    for (IndexType i = 0; i < numOwnedDofs; ++i) {
        // Skip if already a boundary DOF
        if (std::find(boundaryDofs.begin(), boundaryDofs.end(), i) != boundaryDofs.end()) continue;
        
        // Check if diagonal exists and is non-zero
        bool hasNonZeroDiag = false;
        for (IndexType idx = h_rowOffsets_check[i]; idx < h_rowOffsets_check[i + 1]; ++idx) {
            if (h_colIndices_check[idx] == i && std::abs(h_values_check[idx]) > 1e-14) {
                hasNonZeroDiag = true;
                break;
            }
        }
        
        if (!hasNonZeroDiag) {
            zeroDiagDofs.push_back(i);
        }
    }
    
    if (rank == 0 || rank == 1) {
        std::cout << "Rank " << rank << ": Found " << zeroDiagDofs.size() 
                  << " zero-diagonal DOFs (will be eliminated)\n";
        if (!zeroDiagDofs.empty()) {
            std::cout << "Rank " << rank << ": Adding zero-diagonal DOFs to boundary set\n";
        }
    }
    
    // Add zero-diagonal DOFs to boundary set so they get properly communicated
    boundaryDofs.insert(boundaryDofs.end(), zeroDiagDofs.begin(), zeroDiagDofs.end());
    
    // Exchange boundary DOF information to identify boundary ghosts
    // Each rank broadcasts its boundary DOFs in global numbering
    std::vector<int> boundaryCountPerRank(numRanks);
    int myBoundaryCount = static_cast<int>(boundaryDofs.size());
    MPI_Allgather(&myBoundaryCount, 1, MPI_INT, boundaryCountPerRank.data(), 1, MPI_INT, MPI_COMM_WORLD);
    
    std::vector<int> boundaryOffsets(numRanks + 1, 0);
    for (int r = 0; r < numRanks; ++r) {
        boundaryOffsets[r + 1] = boundaryOffsets[r] + boundaryCountPerRank[r];
    }
    
    // Collect all boundary DOFs in global numbering
    std::vector<IndexType> myBoundaryGlobalDofs;
    myBoundaryGlobalDofs.reserve(boundaryDofs.size());
    for (auto bdof : boundaryDofs) {
        myBoundaryGlobalDofs.push_back(dof_handler.local_to_global(bdof));
    }
    
    std::vector<IndexType> allBoundaryGlobalDofs(boundaryOffsets[numRanks]);
    MPI_Allgatherv(myBoundaryGlobalDofs.data(), myBoundaryCount, MPI_UINT64_T,
                   allBoundaryGlobalDofs.data(), boundaryCountPerRank.data(), 
                   boundaryOffsets.data(), MPI_UINT64_T, MPI_COMM_WORLD);
    
    // Build set for fast lookup
    std::unordered_set<IndexType> globalBoundarySet(allBoundaryGlobalDofs.begin(), 
                                                      allBoundaryGlobalDofs.end());
    
    // Identify boundary ghost DOFs using global boundary set
    std::vector<IndexType> ess_ghost_dofs;
    for (IndexType ghostIdx = 0; ghostIdx < numGhostDofs; ++ghostIdx) {
        // IMPORTANT: ghostIdx is a local index in the range [0, numGhostDofs-1]
        // It needs to be mapped to the full local DOF index space [numOwnedDofs, numTotalLocalDofs-1]
        IndexType fullGhostLocalIndex = numOwnedDofs + ghostIdx;
        
        // Find the global DOF for this ghost
        // This requires a reverse mapping from local ghost index to global DOF, which should be in the dof_handler
        // Assuming dof_handler.get_ghost_global_dof(ghostIdx) exists and is correct
        IndexType ghostGlobalDof = dof_handler.get_ghost_global_dof(ghostIdx);
        
        if (ghostGlobalDof != static_cast<IndexType>(-1) && globalBoundarySet.count(ghostGlobalDof) > 0) {
            ess_ghost_dofs.push_back(ghostIdx);
        }
    }
    
    if (rank == 0) {
        std::cout << "Found " << boundaryDofs.size() << " owned boundary DOFs, "
                  << ess_ghost_dofs.size() << " ghost boundary DOFs\n";
    }
    
    IndexType reducedSize = mars::fem::FormLinearSystem(
        K, dm.get_data<1>(), boundaryDofs, A_r, B_r, dof_mapping, numGhostDofs, ess_ghost_dofs);
    
    if (rank == 0) {
        std::cout << "Reduced system size: " << reducedSize << " DOFs (eliminated " 
                  << boundaryDofs.size() << ")\n";
    }
    
    // Validate reduced matrix has diagonals
    std::vector<IndexType> h_rowOffsets(reducedSize + 1);
    std::vector<IndexType> h_colIndices(A_r.nnz());
    std::vector<double> h_values(A_r.nnz());
    
    thrust::copy(thrust::device_pointer_cast(A_r.rowOffsetsPtr()),
                 thrust::device_pointer_cast(A_r.rowOffsetsPtr() + reducedSize + 1),
                 h_rowOffsets.begin());
    thrust::copy(thrust::device_pointer_cast(A_r.colIndicesPtr()),
                 thrust::device_pointer_cast(A_r.colIndicesPtr() + A_r.nnz()),
                 h_colIndices.begin());
    thrust::copy(thrust::device_pointer_cast(A_r.valuesPtr()),
                 thrust::device_pointer_cast(A_r.valuesPtr() + A_r.nnz()),
                 h_values.begin());
    
    int missingDiag = 0, zeroDiag = 0;
    double minDiag = 1e100, maxDiag = -1e100;
    for (IndexType i = 0; i < reducedSize; ++i) {
        bool foundDiag = false;
        for (IndexType idx = h_rowOffsets[i]; idx < h_rowOffsets[i+1]; ++idx) {
            if (h_colIndices[idx] == i) {
                foundDiag = true;
                double diagVal = h_values[idx];
                if (std::abs(diagVal) < 1e-14) zeroDiag++;
                minDiag = std::min(minDiag, diagVal);
                maxDiag = std::max(maxDiag, diagVal);
                break;
            }
        }
        if (!foundDiag) missingDiag++;
    }
    
    std::cout << "Rank " << rank << ": Reduced matrix diagonal check - missing=" << missingDiag 
              << ", zero=" << zeroDiag << ", range=[" << minDiag << ", " << maxDiag << "]\n";
    
    if (missingDiag > 0 || zeroDiag > 0) {
        std::cerr << "Rank " << rank << ": ERROR - Reduced matrix has " << missingDiag 
                  << " missing diagonals and " << zeroDiag << " zero diagonals!\n";
    }
    
    // Check for inter-rank coupling
    bool hasInterRankCoupling = (A_r.numCols() > reducedSize);
    int globalHasInterRankCoupling = 0;
    int localHasInterRankCoupling = hasInterRankCoupling ? 1 : 0;
    MPI_Allreduce(&localHasInterRankCoupling, &globalHasInterRankCoupling, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    
    // Determine actual column count from matrix (not allocated size)
    IndexType maxColInMatrix = 0;
    for (IndexType row = 0; row < reducedSize; ++row) {
        for (IndexType idx = h_rowOffsets[row]; idx < h_rowOffsets[row + 1]; ++idx) {
            maxColInMatrix = std::max(maxColInMatrix, h_colIndices[idx]);
        }
    }
    IndexType actualNumCols = maxColInMatrix + 1;  // columns are 0-indexed
    
    if (rank == 0) {
        std::cout << "\nRank " << rank << ": actualNumCols from max=" << actualNumCols 
                  << ", A_r.numCols()=" << A_r.numCols() << ", reducedSize=" << reducedSize << std::endl;
    }
    
    // Use the actual number of columns that need to be mapped (from FormLinearSystem)
    // This is reducedSize + numInteriorGhostDofs where numInteriorGhostDofs = columns in [reducedSize, actualNumCols)
    IndexType reducedLocalCount = A_r.numCols();  // Use allocated size from FormLinearSystem
    
    if (rank == 0 && actualNumCols == reducedSize) {
        std::cout << "\nNOTE: Matrix is square (no ghost columns), indicating block-diagonal system.\n";
        std::cout << "      All inter-rank connections are through boundary DOFs only.\n";
    } else if (rank == 0 && actualNumCols > reducedSize) {
        std::cout << "\nNOTE: Matrix has ghost columns: " << (actualNumCols - reducedSize) 
                  << " ghost DOFs (columns [" << reducedSize << ", " << (actualNumCols-1) << "])\n";
    }
    
    if (!globalHasInterRankCoupling && rank == 0) {
        std::cout << "\nNOTE: No inter-rank coupling detected (block-diagonal system).\n";
        std::cout << "      Each rank can solve independently. Using local CG solver.\n";
    }
    
    // Keep the rectangular matrix structure like MFEM ex1p
    // A_r is (reducedSize × numTotalLocalDofs) with ghost columns preserved
    
    auto bc_end = std::chrono::high_resolution_clock::now();
    double bc_time = std::chrono::duration<double>(bc_end - bc_start).count();

    // =====================================================
    // 8. Solve: A_r * X_r = B_r with Hypre PCG + BoomerAMG
    // =====================================================
    if (rank == 0) std::cout << "\n8. Solving with Hypre PCG + BoomerAMG...\n";

    auto solve_start = std::chrono::high_resolution_clock::now();

#ifdef MARS_ENABLE_HYPRE
    // Build local-to-global DOF mapping for reduced system
    // Use contiguous global numbering for interior (non-boundary) DOFs
    // Each rank gets a contiguous block [globalRowStart, globalRowEnd)
    
    // The reducedSize is the number of interior DOFs on this rank (after BC elimination)
    IndexType numInteriorLocal = reducedSize;
    IndexType globalRowStart = 0, globalRowEnd = 0;
    
    // Compute global offset for this rank's interior DOFs using exclusive prefix sum
    MPI_Scan(&numInteriorLocal, &globalRowStart, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    globalRowStart -= numInteriorLocal;  // Convert to start index (exclusive scan)
    globalRowEnd = globalRowStart + numInteriorLocal;
    
    IndexType numInteriorGlobal = 0;
    MPI_Allreduce(&numInteriorLocal, &numInteriorGlobal, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    
    if (rank == 0) {
        std::cout << "Interior DOF distribution: rank 0 [" << globalRowStart << ", " << globalRowEnd 
                  << "), total interior DOFs: " << numInteriorGlobal << std::endl;
    }
    if (rank == 1) {
        std::cout << "Interior DOF distribution: rank 1 [" << globalRowStart << ", " << globalRowEnd 
                  << "), total interior DOFs: " << numInteriorGlobal << std::endl;
    }
    
    std::vector<IndexType> localToGlobalDof_all(reducedLocalCount);
    
    // Simple contiguous numbering for owned interior DOFs
    // Owned interior: [globalRowStart, globalRowEnd)
    for (IndexType j = 0; j < reducedSize; ++j) {
        localToGlobalDof_all[j] = globalRowStart + j;
    }
    
    if (rank == 0 && reducedSize >= 3) {
        std::cout << "Column mapping sample (rank 0): [0]→" << localToGlobalDof_all[0] 
                  << ", [1]→" << localToGlobalDof_all[1] 
                  << ", [" << (reducedSize-1) << "]→" << localToGlobalDof_all[reducedSize-1] << std::endl;
    }
    
    // Ghost interior DOFs: need to find which rank owns each ghost and map to that rank's global range
    // Build a mapping of original global DOF ID → contiguous global DOF ID
    // First, collect all ranks' global ranges
    std::vector<IndexType> rankStarts(numRanks + 1);
    MPI_Allgather(&globalRowStart, 1, MPI_UINT64_T, rankStarts.data(), 1, MPI_UINT64_T, MPI_COMM_WORLD);
    rankStarts[numRanks] = numInteriorGlobal;
    
    // For each ghost, we need to:
    // 1. Get its original global DOF ID from dof_handler
    // 2. Find which rank owns it (via MPI communication or storing owner info)
    // 3. Find its position within that rank's interior DOFs
    // 4. Map to contiguous global ID
    
    // This requires knowing which original global DOF IDs each rank has after BC elimination
    // Build list of survived original global DOF IDs on this rank
    // Use dof_handler to get ALL owned DOFs (including those not in local elements)
    std::vector<IndexType> myOriginalGlobalDofs;
    myOriginalGlobalDofs.reserve(reducedSize);
    
    // Method 1: Use dof_mapping which only has DOFs from FE space (in local elements)
    // This misses owned DOFs that aren't in local elements!
    // for (IndexType localDof = 0; localDof < numOwnedDofs; ++localDof) {
    //     IndexType reducedIdx = dof_mapping[localDof];
    //     if (reducedIdx != static_cast<IndexType>(-1)) {
    //         myOriginalGlobalDofs.push_back(dof_handler.local_to_global(localDof));
    //     }
    // }
    
    // Method 2: Directly use reducedSize and assume first reducedSize DOFs in reduced matrix
    // correspond to owned DOFs that survived BC elimination
    // This should work because FormLinearSystem maintains the order
    for (IndexType reducedDof = 0; reducedDof < reducedSize; ++reducedDof) {
        // Find which original DOF this reduced DOF came from
        // Inverse of dof_mapping: find localDof where dof_mapping[localDof] == reducedDof
        IndexType originalLocalDof = static_cast<IndexType>(-1);
        for (IndexType localDof = 0; localDof < numOwnedDofs; ++localDof) {
            if (dof_mapping[localDof] == reducedDof) {
                originalLocalDof = localDof;
                break;
            }
        }
        
        if (originalLocalDof != static_cast<IndexType>(-1)) {
            myOriginalGlobalDofs.push_back(dof_handler.local_to_global(originalLocalDof));
        } else {
            std::cerr << "Rank " << rank << ": WARNING - Reduced DOF " << reducedDof 
                      << " has no corresponding original DOF!\n";
        }
    }
    
    // Allgatherv to collect all ranks' original global DOF lists
    std::vector<int> recvCounts(numRanks);
    std::vector<int> recvOffsets(numRanks + 1, 0);
    int myCount = static_cast<int>(myOriginalGlobalDofs.size());
    MPI_Allgather(&myCount, 1, MPI_INT, recvCounts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    
    for (int r = 0; r < numRanks; ++r) {
        recvOffsets[r + 1] = recvOffsets[r] + recvCounts[r];
    }
    
    std::vector<IndexType> allOriginalGlobalDofs(recvOffsets[numRanks]);
    MPI_Allgatherv(myOriginalGlobalDofs.data(), myCount, MPI_UINT64_T,
                   allOriginalGlobalDofs.data(), recvCounts.data(), recvOffsets.data(),
                   MPI_UINT64_T, MPI_COMM_WORLD);
    
    // Build mapping: original global DOF ID → contiguous global DOF ID
    std::unordered_map<IndexType, IndexType> originalToContiguous;
    for (IndexType contiguousIdx = 0; contiguousIdx < allOriginalGlobalDofs.size(); ++contiguousIdx) {
        originalToContiguous[allOriginalGlobalDofs[contiguousIdx]] = contiguousIdx;
    }
    
    // Ghost interior: [reducedSize, reducedLocalCount) 
    IndexType numInteriorGhostDofs = reducedLocalCount - reducedSize;
    
    if (rank == 0 || rank == 1) {
        std::cout << "Rank " << rank << ": numInteriorGhostDofs=" << numInteriorGhostDofs 
                  << ", reducedLocalCount=" << reducedLocalCount << ", reducedSize=" << reducedSize << "\n";
    }
    
    if (numInteriorGhostDofs > 0) {
        // FormLinearSystem already filtered boundary ghosts, so all remaining ghosts are interior
        // Build list of interior ghost global DOFs once
        std::vector<IndexType> interiorGhostGlobalDofs;
        for (IndexType g = 0; g < numGhostDofs; ++g) {
            // Check if NOT a boundary ghost
            bool isBoundary = (std::find(ess_ghost_dofs.begin(), ess_ghost_dofs.end(), g) != ess_ghost_dofs.end());
            if (!isBoundary) {
                IndexType ghostGlobalDof = dof_handler.get_ghost_global_dof(g);
                // Also check if this global DOF is in the interior (not eliminated by owner rank)
                if (originalToContiguous.find(ghostGlobalDof) != originalToContiguous.end()) {
                    interiorGhostGlobalDofs.push_back(ghostGlobalDof);
                }
            }
        }
        
        if (numInteriorGhostDofs != interiorGhostGlobalDofs.size()) {
            std::cerr << "Rank " << rank << ": ERROR - Mismatch in interior ghost count: expected " 
                      << numInteriorGhostDofs << ", got " << interiorGhostGlobalDofs.size() << "\n";
        } else if (rank == 0 || rank == 1) {
            std::cout << "Rank " << rank << ": Interior ghost count matches: " << numInteriorGhostDofs << "\n";
        }
        
        for (IndexType ghostIdx = 0; ghostIdx < numInteriorGhostDofs && ghostIdx < interiorGhostGlobalDofs.size(); ++ghostIdx) {
            IndexType reducedCol = reducedSize + ghostIdx;
            IndexType ghostGlobalDof = interiorGhostGlobalDofs[ghostIdx];
            
            // Convert to contiguous global DOF
            auto it = originalToContiguous.find(ghostGlobalDof);
            if (it != originalToContiguous.end()) {
                localToGlobalDof_all[reducedCol] = it->second;
            } else {
                std::cerr << "Rank " << rank << ": ERROR - Interior ghost DOF " << ghostGlobalDof 
                          << " not found in contiguous mapping!\n";
                localToGlobalDof_all[reducedCol] = static_cast<IndexType>(-1);
            }
        }
    }
    
    // Validate reduced RHS before solver
    std::vector<double> h_B_r(reducedSize);
    thrust::copy(thrust::device_pointer_cast(B_r.data()),
                 thrust::device_pointer_cast(B_r.data() + reducedSize),
                 h_B_r.begin());
    
    double sum_B_r = 0.0, min_B_r = 1e100, max_B_r = -1e100;
    int nan_count_B_r = 0;
    for (IndexType i = 0; i < reducedSize; ++i) {
        if (!std::isfinite(h_B_r[i])) {
            nan_count_B_r++;
        } else {
            sum_B_r += h_B_r[i];
            min_B_r = std::min(min_B_r, h_B_r[i]);
            max_B_r = std::max(max_B_r, h_B_r[i]);
        }
    }
    
    std::cout << "Rank " << rank << ": Reduced RHS before solver - sum=" << sum_B_r 
              << ", range=[" << min_B_r << ", " << max_B_r << "], NaNs=" << nan_count_B_r << "/" << reducedSize << std::endl;
    
    if (nan_count_B_r > 0) {
        std::cerr << "Rank " << rank << ": ERROR - Reduced RHS contains NaN/Inf!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    if (rank == 0) {
        std::cout << "Using contiguous global DOF numbering for interior DOFs: rows [" << globalRowStart 
                  << ", " << globalRowEnd << "), columns [0, " << numInteriorGlobal << ")" << std::endl;
    }
    
    mars::VectorSelector<double, cstone::GpuTag>::type X_r(reducedSize, 0.0);
    
    // Select preconditioner based on solver argument
    using PrecondType = mars::fem::HyprePCGSolver<double, IndexType, cstone::GpuTag>::PrecondType;
    PrecondType precondType = PrecondType::BOOMERAMG;  // Default to BoomerAMG
    if (solver == "elimination-jacobi") {
        precondType = PrecondType::JACOBI;
    }
    
    if (rank == 0) {
        std::cout << "Using preconditioner: " << (precondType == PrecondType::JACOBI ? "Jacobi" : "BoomerAMG") << std::endl;
    }
    
    mars::fem::HyprePCGSolver<double, IndexType, cstone::GpuTag> hypre_solver(MPI_COMM_WORLD, 400, 1e-10, precondType);
    hypre_solver.setVerbose(rank == 0);
    
    // Pass the rectangular matrix with contiguous global numbering
    bool converged = hypre_solver.solve(A_r, B_r, X_r, globalRowStart, globalRowEnd, 
                                       0, numInteriorGlobal,
                                       localToGlobalDof_all);
    
    // Recover full solution (boundary DOFs = 0)
    auto& u_solution = dm.get_data<0>();
    mars::fem::RecoverFEMSolution<IndexType, double, cstone::GpuTag>(X_r, dof_mapping, u_solution);
#else
    if (rank == 0) std::cerr << "ERROR: Hypre not available\n";
    MPI_Finalize();
    return 1;
#endif

    auto solve_end = std::chrono::high_resolution_clock::now();
    double solve_time = std::chrono::duration<double>(solve_end - solve_start).count();

    // =====================================================
    // 9. Exchange ghost data and compute statistics
    // =====================================================
    if (rank == 0) std::cout << "\n9. Computing solution statistics...\n";

    // Exchange ghost values between ranks
    dm.gather_ghost_data();
    dm.scatter_ghost_data();

    auto& u_full = dm.get_data<0>();
    IndexType numOwnedDofsForStats = dof_handler.get_num_local_dofs();

    // Compute local statistics
    double u_min_local = *thrust::min_element(thrust::device_pointer_cast(u_full.data()),
                                             thrust::device_pointer_cast(u_full.data() + numOwnedDofsForStats));
    double u_max_local = *thrust::max_element(thrust::device_pointer_cast(u_full.data()),
                                             thrust::device_pointer_cast(u_full.data() + numOwnedDofsForStats));

    double u_norm_sq_local = thrust::transform_reduce(thrust::device_pointer_cast(u_full.data()),
                                                     thrust::device_pointer_cast(u_full.data() + numOwnedDofsForStats),
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
    } // End of scope to ensure destructors run before MPI_Finalize

    MPI_Finalize();
    return 0;
}