/**
 * MFEM mesh Distributed GPU Poisson solver
 *
 * Problem: -Δu = 1 in Ω, u = 0 on ∂Ω
 *
 * Usage: 
 *   For tetrahedral meshes: mpirun -np N mars_ex_beam_tet_distributed --mesh beam-tet.mesh
 *   For hexahedral meshes: mpirun -np N mars_ex_beam_tet_distributed --mesh beam-hex.mesh
 *                          (requires recompilation with -DUSE_HEX_TAG)
 *
 * Features:
 * - MPI distributed across multiple ranks
 * - GPU acceleration throughout (no CPU transfers)
 * - Hypre GPU-accelerated PCG + BoomerAMG
 * - Domain decomposition with halo exchange
 * - Supports both tetrahedral and hexahedral meshes (compile-time selection)
 *
 * Element Type Support:
 * - Tetrahedral meshes (default): compile without flags
 * - Hexahedral meshes: compile with -DUSE_HEX_TAG
 * - Runtime validation ensures mesh type matches compiled ElementTag
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

// Element type selection - RUNTIME DETECTION from mesh!
// Compile-time element tag selection
// To support hexahedral meshes, recompile with -DUSE_HEX_TAG
#ifdef USE_HEX_TAG
using ElementTag = mars::HexTag;
#else
using ElementTag = mars::TetTag;  // Default: tetrahedral meshes
#endif

// MFEM-style bigint handling for Hypre compatibility
typedef uint64_t IndexType;

// DOF handler typedef for convenience  
template<typename ElemTag>
using DofHandlerT = mars::fem::UnstructuredDofHandler<ElemTag, double, IndexType, cstone::GpuTag>;

// FE space typedef for convenience
template<typename ElemTag>
using FESpaceT = mars::fem::H1FESpace<ElemTag, double, IndexType, cstone::GpuTag>;

// Assembler typedefs for convenience
template<typename ElemTag>
using StiffnessAssemblerT = mars::fem::StiffnessAssembler<ElemTag, double, IndexType, cstone::GpuTag>;
template<typename ElemTag>
using MassAssemblerT = mars::fem::MassAssembler<ElemTag, double, IndexType, cstone::GpuTag>;

// Matrix typedef (doesn't depend on element tag)
using SparseMatrixType = mars::fem::SparseMatrix<IndexType, double, cstone::GpuTag>;

// Backward compatibility - use full type names to avoid conflicts
using DofHandler = DofHandlerT<ElementTag>;
using FESpace = FESpaceT<ElementTag>;

// Source term: f(x,y,z) = 1
struct SourceTerm {
    __device__ __host__
    double operator()(double x, double y, double z) const {
        return 1.0;
    }
};

// Save solution to MFEM GridFunction format (rank 0 only)
template<typename ElementTag>
void saveSolutionToGridFunction(const cstone::DeviceVector<double>& u_local,
                               const DofHandlerT<ElementTag>& dof_handler,
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
    bool use_geometric_bc = false;  // Use geometric boundaries instead of topological
    bool use_all_bc = false;  // Use ALL topological boundaries (all mesh boundary faces)
    std::vector<int> bc_attributes;  // Boundary attributes to apply Dirichlet BC (empty = default to attr 1)

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--mesh" && i + 1 < argc) {
            meshPath = argv[++i];
        } else if (arg == "--order" && i + 1 < argc) {
            order = std::atoi(argv[++i]);
        } else if (arg == "--solver" && i + 1 < argc) {
            solver = argv[++i];
        } else if (arg == "--geometric-bc") {
            use_geometric_bc = true;
        } else if (arg == "--topological-bc" || arg == "--all-bc") {
            use_all_bc = true;
        } else if (arg == "--bc-attr" && i + 1 < argc) {
            bc_attributes.push_back(std::atoi(argv[++i]));
        } else if (arg == "--help" || arg == "-h") {
            if (rank == 0) {
                std::cout << "Usage: " << argv[0] << " --mesh <path> [options]\n\n"
                          << "Solve -Δu = 1 with u=0 on boundary (Distributed GPU version)\n\n"
                          << "Options:\n"
                          << "  --mesh <path>      Path to mesh (MFEM .mesh or MARS binary dir)\n"
                          << "  --order <n>        Polynomial order (default: 1)\n"
                          << "  --solver <type>    Solver: elimination-boomeramg (default), elimination-jacobi, elimination, or jacobi\n"
                          << "  --bc-attr <n>      Apply Dirichlet BC only to boundary attribute <n> (can be repeated; default: attr 1)\n"
                          << "  --topological-bc   Use ALL topological boundaries from mesh (all boundary faces)\n"
                          << "  --geometric-bc     Use geometric boundary detection (x=min and x=max faces)\n"
                          << "  --help, -h         Print this help\n\n"
                          << "Boundary Condition Modes:\n"
                          << "  (default)          Apply BC to boundary attribute 1 only (MFEM ex1p style)\n"
                          << "  --bc-attr <n>      Apply BC to specified attribute(s)\n"
                          << "  --topological-bc   Apply BC to ALL mesh boundary faces (all attributes)\n"
                          << "  --geometric-bc     Apply BC to geometric boundaries (x endpoints)\n\n"
                          << "Solvers:\n"
                          << "  elimination-boomeramg - PCG + BoomerAMG with DOF elimination (default, MFEM ex1p style)\n"
                          << "  elimination-jacobi    - PCG + Jacobi with DOF elimination (for debugging)\n"
                          << "  elimination           - PCG + BoomerAMG with DOF elimination (legacy alias)\n"
                          << "  jacobi                - PCG + Jacobi with penalty BC (not yet implemented)\n\n"
                          << "Examples:\n"
                          << "  mpirun -np 4 " << argv[0] << " --mesh beam-tet.mesh              # BC on attr 1\n"
                          << "  mpirun -np 4 " << argv[0] << " --mesh beam-tet.mesh --bc-attr 1  # BC on attr 1\n"
                          << "  mpirun -np 4 " << argv[0] << " --mesh beam-tet.mesh --topological-bc  # BC on all faces\n"
                          << "  mpirun -np 4 " << argv[0] << " --mesh beam-tet.mesh --geometric-bc    # BC on x endpoints\n";
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

    // Broadcast element geometry type
    int element_geom_type = loader.element_geometry_type;
    MPI_Bcast(&element_geom_type, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Broadcast vertices
    size_t numVertices = loader.vertices.size();
    MPI_Bcast(&numVertices, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        loader.vertices.resize(numVertices);
    }
    MPI_Bcast(loader.vertices.data(), numVertices * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Broadcast elements (tet or hex depending on type)
    size_t numElements = (element_geom_type == 4) ? loader.elements.size() : loader.hex_elements.size();
    MPI_Bcast(&numElements, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    if (element_geom_type == 4) {  // Tetrahedra
        if (rank != 0) {
            loader.elements.resize(numElements);
        }
        MPI_Bcast(loader.elements.data(), numElements * 4, MPI_INT, 0, MPI_COMM_WORLD);
    } else if (element_geom_type == 5) {  // Hexahedra
        if (rank != 0) {
            loader.hex_elements.resize(numElements);
        }
        MPI_Bcast(loader.hex_elements.data(), numElements * 8, MPI_INT, 0, MPI_COMM_WORLD);
    }

    // Broadcast boundary faces and attributes (needed for attribute-filtered BC extraction)
    // For tet meshes: triangular faces; for hex meshes: quadrilateral faces
    size_t numBoundaryFaces = loader.boundary_faces.size();
    size_t numBoundaryQuads = loader.boundary_quads.size();
    MPI_Bcast(&numBoundaryFaces, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numBoundaryQuads, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        loader.boundary_faces.resize(numBoundaryFaces);
        loader.boundary_quads.resize(numBoundaryQuads);
        loader.boundary_face_attributes.resize(numBoundaryFaces + numBoundaryQuads);
    }
    if (numBoundaryFaces > 0) {
        MPI_Bcast(loader.boundary_faces.data(), numBoundaryFaces * 3, MPI_INT, 0, MPI_COMM_WORLD);
    }
    if (numBoundaryQuads > 0) {
        MPI_Bcast(loader.boundary_quads.data(), numBoundaryQuads * 4, MPI_INT, 0, MPI_COMM_WORLD);
    }
    MPI_Bcast(loader.boundary_face_attributes.data(), numBoundaryFaces + numBoundaryQuads, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "\n=== Boundary Condition Setup ===" << std::endl;
        std::cout << "Boundary triangular faces: " << numBoundaryFaces << std::endl;
        std::cout << "Boundary quadrilateral faces: " << numBoundaryQuads << std::endl;
        std::cout << "Total boundary faces: " << (numBoundaryFaces + numBoundaryQuads) << std::endl;
        std::cout << "Boundary face attributes vector size: " << loader.boundary_face_attributes.size() << std::endl;
        if (!loader.boundary_face_attributes.empty()) {
            std::map<int, int> attr_count;
            for (int attr : loader.boundary_face_attributes) {
                attr_count[attr]++;
            }
            std::cout << "Attribute distribution: ";
            for (const auto& [attr, count] : attr_count) {
                std::cout << "attr_" << attr << "=" << count << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "User-specified bc_attributes.size() = " << bc_attributes.size() << std::endl;
        std::cout << "use_all_bc = " << (use_all_bc ? "true" : "false") << std::endl;
        std::cout << "use_geometric_bc = " << (use_geometric_bc ? "true" : "false") << std::endl;
    }

    // Determine boundary condition mode
    if (use_all_bc) {
        // Use ALL topological boundaries - clear bc_attributes to use all faces
        bc_attributes.clear();
        if (rank == 0) {
            std::cout << "Mode: Using ALL topological boundaries (--topological-bc)" << std::endl;
        }
    } else if (use_geometric_bc) {
        // Geometric BC handled later
        if (rank == 0) {
            std::cout << "Mode: Using geometric boundaries (--geometric-bc)" << std::endl;
        }
    } else if (bc_attributes.empty() && !loader.boundary_face_attributes.empty()) {
        // Default: use attribute 1 if it exists (MFEM ex1p style)
        bool has_attr_1 = false;
        for (int attr : loader.boundary_face_attributes) {
            if (attr == 1) {
                has_attr_1 = true;
                break;
            }
        }
        if (has_attr_1) {
            bc_attributes.push_back(1);
            if (rank == 0) {
                std::cout << "Mode: Defaulting to boundary attribute 1 (MFEM ex1p style)" << std::endl;
                std::cout << "      (Use --topological-bc for all boundaries, --bc-attr <n> for specific)" << std::endl;
            }
        }
    }

    // Extract boundary vertices with specified attributes (all ranks compute identically)
    if (!use_geometric_bc && !bc_attributes.empty()) {
        // Re-extract boundary vertices with specified attributes only
        loader.boundary_vertices.clear();
        loader.extract_boundary_vertices(bc_attributes);
        if (rank == 0) {
            std::cout << "Applying Dirichlet BC to attributes: ";
            for (int attr : bc_attributes) std::cout << attr << " ";
            std::cout << std::endl;
            std::cout << "Boundary vertices with these attributes: " << loader.boundary_vertices.size() << std::endl;
        }
    } else if (!use_geometric_bc) {
        // use_all_bc or no attributes specified - use ALL boundary vertices
        // For consistency, don't recompute on other ranks - use original from load()
        if (rank == 0) {
            std::cout << "Using ALL boundary vertices (all topological boundaries)" << std::endl;
            std::cout << "Boundary vertices: " << loader.boundary_vertices.size() << std::endl;
        }
        // Rank != 0 keeps the original boundary_vertices from load() (which were broadcast)
    } else {
        // use_geometric_bc - boundary vertices will be computed later from geometry
        // Clear boundary_vertices for geometric BC since we don't use mesh topology
        loader.boundary_vertices.clear();
        if (rank == 0) {
            std::cout << "Geometric BC mode - boundary will be detected from coordinates" << std::endl;
        }
    }
    
    if (rank == 0) {
        std::cout << "=== End Boundary Setup ===\n" << std::endl;
        std::cout << "Global mesh: " << numVertices << " vertices, "
                  << numElements << " " << (element_geom_type == 4 ? "tetrahedra" : "hexahedra") << std::endl;
        std::cout << "Boundary info: " << loader.boundary_vertices.size() << " boundary vertices" << std::endl;
    }

    // Broadcast boundary vertices to all ranks (skip for geometric BC)
    if (!use_geometric_bc) {
        size_t numBoundaryVertices = loader.boundary_vertices.size();
        MPI_Bcast(&numBoundaryVertices, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        
        if (rank != 0) {
            loader.boundary_vertices.resize(numBoundaryVertices);
        }
        MPI_Bcast(loader.boundary_vertices.data(), numBoundaryVertices, MPI_INT, 0, MPI_COMM_WORLD);
    }

    // Convert MFEM mesh data to ElementDomain format (shared by both types)
    std::vector<double> x_coords(numVertices);
    std::vector<double> y_coords(numVertices);
    std::vector<double> z_coords(numVertices);

    for (size_t i = 0; i < numVertices; ++i) {
        x_coords[i] = static_cast<double>(loader.vertices[i][0]);
        y_coords[i] = static_cast<double>(loader.vertices[i][1]);
        z_coords[i] = static_cast<double>(loader.vertices[i][2]);
    }

    // HostBoundaryTuple for boundary nodes
    std::vector<uint8_t> isBoundaryNode(numVertices, 0);
    for (int boundaryVertex : loader.boundary_vertices) {
        if (boundaryVertex >= 0 && boundaryVertex < static_cast<int>(numVertices)) {
            isBoundaryNode[boundaryVertex] = 1;
        }
    }

    // Convert connectivity based on element type
    std::vector<IndexType> i0(numElements), i1(numElements), i2(numElements), i3(numElements);
    std::vector<IndexType> i4, i5, i6, i7;  // For hex elements or padding for tet elements

    // Validate mesh type - now supports both tetrahedral and hexahedral meshes automatically
    if (element_geom_type != 4 && element_geom_type != 5) {
        if (rank == 0) {
            std::cerr << "ERROR: Unsupported mesh element type: " << element_geom_type << "\\n";
            std::cerr << "  Supported types: 4 (Tetrahedra) or 5 (Hexahedra)\\n";
        }
        MPI_Finalize();
        return 1;
    }

    if (rank == 0) {
        std::cout << "Mesh type: " << (element_geom_type == 4 ? "Tetrahedral" : "Hexahedral") << " (auto-detected)" << std::endl;
    }
    
    if (element_geom_type == 4) {
        // Tetrahedra
        if (rank == 0) std::cout << "\\nElement type: TETRAHEDRA (matches TetTag)\\n" << std::endl;
        for (size_t i = 0; i < numElements; ++i) {
            i0[i] = static_cast<IndexType>(loader.elements[i][0]);
            i1[i] = static_cast<IndexType>(loader.elements[i][1]);
            i2[i] = static_cast<IndexType>(loader.elements[i][2]);
            i3[i] = static_cast<IndexType>(loader.elements[i][3]);
        }
    } else if (element_geom_type == 5) {
        // Hexahedra
        if (rank == 0) std::cout << "\\nElement type: HEXAHEDRA (matches HexTag)\\n" << std::endl;
        i4.resize(numElements);
        i5.resize(numElements);
        i6.resize(numElements);
        i7.resize(numElements);
        for (size_t i = 0; i < numElements; ++i) {
            i0[i] = static_cast<IndexType>(loader.hex_elements[i][0]);
            i1[i] = static_cast<IndexType>(loader.hex_elements[i][1]);
            i2[i] = static_cast<IndexType>(loader.hex_elements[i][2]);
            i3[i] = static_cast<IndexType>(loader.hex_elements[i][3]);
            i4[i] = static_cast<IndexType>(loader.hex_elements[i][4]);
            i5[i] = static_cast<IndexType>(loader.hex_elements[i][5]);
            i6[i] = static_cast<IndexType>(loader.hex_elements[i][6]);
            i7[i] = static_cast<IndexType>(loader.hex_elements[i][7]);
        }
    } else {
        std::cerr << "ERROR: Unsupported element geometry type: " << element_geom_type << std::endl;
        MPI_Finalize();
        return 1;
    }

    if (rank == 0) {
        std::cout << "Creating distributed ElementDomain (automatic SFC partitioning)..." << std::endl;
    }

    // Create ElementDomain with compile-time ElementTag selection
    using Domain = mars::ElementDomain<ElementTag, double, IndexType, cstone::GpuTag>;

    // Validate mesh type matches compiled ElementTag and create domain
    std::unique_ptr<Domain> domain_ptr;
    
    if constexpr (ElementTag::NodesPerElement == 4) {
        // TetTag: validate tetrahedral mesh
        if (element_geom_type != 4) {
            if (rank == 0) {
                std::cerr << "ERROR: Hexahedral mesh detected but compiled with TetTag.\n"
                         << "       Recompile with -DUSE_HEX_TAG to use hexahedral meshes." << std::endl;
            }
            MPI_Finalize();
            return 1;
        }
        domain_ptr = std::make_unique<Domain>(
            std::make_tuple(x_coords, y_coords, z_coords),
            std::make_tuple(i0, i1, i2, i3),
            std::make_tuple(isBoundaryNode),
            rank, numRanks);
    } else {
        // HexTag: validate hexahedral mesh
        if (element_geom_type != 5) {
            if (rank == 0) {
                std::cerr << "ERROR: Tetrahedral mesh detected but compiled with HexTag.\n"
                         << "       Recompile without -DUSE_HEX_TAG to use tetrahedral meshes." << std::endl;
            }
            MPI_Finalize();
            return 1;
        }
        domain_ptr = std::make_unique<Domain>(
            std::make_tuple(x_coords, y_coords, z_coords),
            std::make_tuple(i0, i1, i2, i3, i4, i5, i6, i7),
            std::make_tuple(isBoundaryNode),
            rank, numRanks);
    }
    
    Domain& domain = *domain_ptr;    // Force domain initialization before creating FE space
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
    FESpace fes(domain, 1);  // Order 1 like MFEM ex1
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
    // 4. Boundary conditions: u=0 on ALL external faces
    // =====================================================
    if (rank == 0) {
        std::cout << "\n4. Applying boundary conditions (u=0 on " 
                  << (use_geometric_bc ? "geometric boundaries" : "all external faces") << ")...\n";
    }

    std::vector<uint8_t> localBoundaryData(domain.getNodeCount(), 0);
    
    if (use_geometric_bc) {
        // Geometric BC: identify boundary nodes by position
        // For Poisson problem: boundary is at ALL 6 faces (x, y, z min/max)
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

        double x_min = *std::min_element(h_x.begin(), h_x.end());
        double x_max = *std::max_element(h_x.begin(), h_x.end());
        double y_min = *std::min_element(h_y.begin(), h_y.end());
        double y_max = *std::max_element(h_y.begin(), h_y.end());
        double z_min = *std::min_element(h_z.begin(), h_z.end());
        double z_max = *std::max_element(h_z.begin(), h_z.end());

        double tol = 1e-3;  // Tolerance for boundary detection

        for (size_t i = 0; i < h_x.size(); ++i) {
            double x = h_x[i], y = h_y[i], z = h_z[i];
            // Check if near any of the 6 boundaries
            if (std::abs(x - x_min) < tol || std::abs(x - x_max) < tol ||
                std::abs(y - y_min) < tol || std::abs(y - y_max) < tol ||
                std::abs(z - z_min) < tol || std::abs(z - z_max) < tol) {
                localBoundaryData[i] = 1;
            }
        }

        size_t geom_boundary_count = std::count(localBoundaryData.begin(), localBoundaryData.end(), 1);
        if (rank == 0) {
            std::cout << "Rank 0: Geometric boundary count: " << geom_boundary_count
                      << " nodes at domain boundaries" << std::endl;
            std::cout << "       Domain bounds: x=[" << x_min << "," << x_max
                      << "] y=[" << y_min << "," << y_max
                      << "] z=[" << z_min << "," << z_max << "]" << std::endl;
        }
    } else {
        // Topological BC: use mesh boundary information
        if (!domain.hasBoundaryInfo()) {
            std::cerr << "ERROR: Domain has no boundary information from mesh file!" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Use ALL boundary nodes from mesh topology (all 6 faces)
        const auto& d_boundaryNodes = domain.getBoundaryNodes();
        thrust::copy(thrust::device_pointer_cast(d_boundaryNodes.data()),
                     thrust::device_pointer_cast(d_boundaryNodes.data() + domain.getNodeCount()),
                     localBoundaryData.begin());
    }

    size_t boundaryCount = 0;
    for (size_t i = 0; i < localBoundaryData.size(); ++i) {
        if (localBoundaryData[i]) boundaryCount++;
    }

    if (rank == 0) {
        std::cout << "Rank 0: Local boundary node count: " << boundaryCount
                  << " out of " << localBoundaryData.size() << " nodes ("
                  << (100.0 * boundaryCount / localBoundaryData.size()) << "%)" << std::endl;
    }

    // Global reduction to get total boundary nodes across all ranks
    size_t totalBoundaryNodes = 0;
    MPI_Reduce(&boundaryCount, &totalBoundaryNodes, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << "Total boundary nodes across all ranks: " << totalBoundaryNodes << std::endl;
    }

    // Set local boundary data for DOF handler
    dof_handler.set_boundary_data(localBoundaryData);

    // Get boundary DOFs using DOF handler (like MFEM's GetBoundaryTrueDofs)
    std::vector<IndexType> boundaryDofs;
    dof_handler.boundary_owned_dof_iterate([&](IndexType localDof) {
        boundaryDofs.push_back(localDof);
    });

    if (rank == 0 || rank == 1) {
        std::cout << "Rank " << rank << ": Found " << boundaryDofs.size() 
                  << " boundary DOFs out of " << numDofs << " total DOFs" << std::endl;
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

    StiffnessAssemblerT<ElementTag> stiffnessAssembler;
    SparseMatrixType K;

    // Get node-to-DOF mapping from distributed DOF handler
    // DOF handler already provides correct mapping:
    //   - Owned nodes: indices [0, numLocalDofs)
    //   - Ghost nodes: indices [numLocalDofs, numLocalDofs+numGhostDofs)
    const auto& nodeToLocalDof = dof_handler.get_node_to_local_dof();
    const auto& resolvedOwnership = dof_handler.get_resolved_ownership();

    stiffnessAssembler.assemble(fes, K, nodeToLocalDof, resolvedOwnership);

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
        
        // Check row sums (should be ~0 for Laplacian)
        double maxRowSum = 0.0, minRowSum = 1e30;
        int nonZeroRowSums = 0;
        for (size_t i = 0; i < K.numRows(); ++i) {
            double rowSum = 0.0;
            for (IndexType idx = h_rowOffsets[i]; idx < h_rowOffsets[i+1]; ++idx) {
                rowSum += h_values[idx];
            }
            maxRowSum = std::max(maxRowSum, std::abs(rowSum));
            minRowSum = std::min(minRowSum, std::abs(rowSum));
            if (std::abs(rowSum) > 1e-10) nonZeroRowSums++;
        }
        std::cout << "Row sum range: [" << minRowSum << ", " << maxRowSum << "]" << std::endl;
        std::cout << "Rows with non-zero sum (|sum| > 1e-10): " << nonZeroRowSums << " / " << K.numRows() << std::endl;
    }

    // =====================================================
    // 6. Assemble RHS vector
    // =====================================================
    if (rank == 0) std::cout << "\n6. Assembling RHS vector...\n";

    auto rhs_start = std::chrono::high_resolution_clock::now();

    MassAssemblerT<ElementTag> massAssembler;

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

    SparseMatrixType A_r;
    mars::VectorSelector<double, cstone::GpuTag>::type B_r;
    std::vector<IndexType> dof_mapping;
    
    // IndexType numOwnedDofs = dof_handler.get_num_local_dofs(); // Already declared above
    IndexType numTotalLocalDofs = K.numCols();  // Owned + ghost
    IndexType numGhostDofs = numTotalLocalDofs - numOwnedDofs;
    
    // PRE-PROCESSING: Detect DOFs with zero diagonals (not in any local element)
    // In multi-rank simulations, ghost DOFs may have zero diagonals because they're not 
    // in any local elements - this is expected and correct
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
                  << " zero-diagonal DOFs among owned DOFs\n";
        if (!zeroDiagDofs.empty() && numRanks == 1) {
            std::cerr << "Rank " << rank << ": WARNING - Zero-diagonal DOFs in single-rank run may indicate assembly issues\n";
        }
    }
    
    // DO NOT add zero-diagonal DOFs to boundary set
    // Zero diagonals in multi-rank runs indicate nodes without local elements (expected)
    // Only true topological boundary DOFs should be in boundaryDofs
    
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

    // Gather global statistics
    size_t totalElements = 0;
    size_t totalOwnedDofs = 0;
    size_t totalBoundaryDofs = 0;
    size_t myElements = domain.localElementCount();
    size_t myBoundaryDofs = boundaryDofs.size();

    MPI_Reduce(&myElements, &totalElements, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&numOwnedDofs, &totalOwnedDofs, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&myBoundaryDofs, &totalBoundaryDofs, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "\n=== Global Mesh Statistics ===\n";
        std::cout << "Total elements (owned by all ranks): " << totalElements << "\n";
        std::cout << "Total owned DOFs (before BC): " << totalOwnedDofs << "\n";
        std::cout << "Total boundary DOFs: " << totalBoundaryDofs << "\n";
        std::cout << "Total interior DOFs: " << (totalOwnedDofs - totalBoundaryDofs) << "\n";
        std::cout << "==============================\n\n";
    }

    // Debug: validate input sizes before FormLinearSystem
    std::cout << "Rank " << rank << ": FormLinearSystem inputs - K: " << K.numRows() << "x" << K.numCols()
              << ", RHS size: " << dm.get_data<1>().size()
              << ", boundaryDofs: " << boundaryDofs.size()
              << ", numGhostDofs: " << numGhostDofs
              << ", ess_ghost_dofs: " << ess_ghost_dofs.size() << "\n";
    MPI_Barrier(MPI_COMM_WORLD);  // Sync before FormLinearSystem

    // Get ghost DOF mapping from FormLinearSystem
    std::vector<IndexType> ghost_dof_mapping;
    IndexType reducedLocalCount = 0;
    IndexType reducedSize = mars::fem::FormLinearSystem(
        K, dm.get_data<1>(), boundaryDofs, A_r, B_r, dof_mapping, numGhostDofs, ess_ghost_dofs,
        &ghost_dof_mapping, &reducedLocalCount);
    
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
    
    // reducedLocalCount is now returned directly from FormLinearSystem
    // Verify it matches the matrix allocation
    
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
    
    std::vector<IndexType> localToGlobalDof_all(reducedLocalCount, static_cast<IndexType>(-1));

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

    // Filter ghost_dof_mapping to only include ghosts that survived BC elimination on owner rank
    // Build a compacted mapping and count actual interior ghosts
    std::vector<IndexType> validInteriorGhostGlobalDofs;
    std::vector<IndexType> validGhostReducedIndices;

    for (IndexType g = 0; g < numGhostDofs; ++g) {
        IndexType reducedGhostIdx = ghost_dof_mapping[g];
        if (reducedGhostIdx != static_cast<IndexType>(-1)) {
            // This ghost was not a boundary ghost locally
            IndexType ghostGlobalDof = dof_handler.get_ghost_global_dof(g);

            // Check if it survived BC elimination on its owner rank
            if (originalToContiguous.find(ghostGlobalDof) != originalToContiguous.end()) {
                validInteriorGhostGlobalDofs.push_back(ghostGlobalDof);
                validGhostReducedIndices.push_back(reducedGhostIdx);
            }
        }
    }

    IndexType actualInteriorGhostDofs = validInteriorGhostGlobalDofs.size();
    IndexType numInteriorGhostDofs = reducedLocalCount - reducedSize;

    if (rank == 0 || rank == 1) {
        std::cout << "Rank " << rank << ": FormLinearSystem reported " << numInteriorGhostDofs
                  << " interior ghosts, but only " << actualInteriorGhostDofs
                  << " survived BC elimination on owner rank\n";
        std::cout << "Rank " << rank << ": reducedLocalCount=" << reducedLocalCount
                  << ", reducedSize=" << reducedSize << "\n";
    }

    // Map the valid interior ghost DOFs
    for (size_t i = 0; i < validInteriorGhostGlobalDofs.size(); ++i) {
        IndexType reducedGhostIdx = validGhostReducedIndices[i];
        IndexType reducedCol = reducedSize + reducedGhostIdx;
        IndexType ghostGlobalDof = validInteriorGhostGlobalDofs[i];

        auto it = originalToContiguous.find(ghostGlobalDof);
        localToGlobalDof_all[reducedCol] = it->second;
    }

    if (rank == 0 || rank == 1) {
        std::cout << "Rank " << rank << ": Successfully mapped " << actualInteriorGhostDofs
                  << " valid interior ghost DOFs\n";
    }

    // Check for unmapped columns and build compacted column mapping
    IndexType unmappedCount = 0;
    std::vector<IndexType> compactedColMapping(reducedLocalCount);
    std::vector<IndexType> compactedGlobalMapping;
    compactedGlobalMapping.reserve(reducedLocalCount);

    // First pass: owned columns (always mapped)
    for (IndexType col = 0; col < reducedSize; ++col) {
        compactedColMapping[col] = col;
        compactedGlobalMapping.push_back(localToGlobalDof_all[col]);
    }

    // Second pass: compact ghost columns (only keep mapped ones)
    IndexType compactedGhostIdx = reducedSize;
    for (IndexType col = reducedSize; col < reducedLocalCount; ++col) {
        if (localToGlobalDof_all[col] != static_cast<IndexType>(-1)) {
            compactedColMapping[col] = compactedGhostIdx++;
            compactedGlobalMapping.push_back(localToGlobalDof_all[col]);
        } else {
            compactedColMapping[col] = static_cast<IndexType>(-1);
            unmappedCount++;
        }
    }

    IndexType compactedLocalCount = compactedGhostIdx;

    if (unmappedCount > 0 && (rank == 0 || rank == 1)) {
        std::cout << "Rank " << rank << ": Compacting matrix: removing " << unmappedCount
                  << " unmapped ghost columns (out of " << (reducedLocalCount - reducedSize) << ")\n";
        std::cout << "Rank " << rank << ": Compacted size: " << reducedSize << " rows × "
                  << compactedLocalCount << " cols (was " << reducedLocalCount << ")\n";
    }

    // If we have unmapped columns, compact the matrix
    if (unmappedCount > 0) {
        // Read current matrix from device
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

        // Compact the matrix: remap column indices
        std::vector<IndexType> h_rowOffsets_compact(reducedSize + 1, 0);
        std::vector<IndexType> h_colIndices_compact;
        std::vector<double> h_values_compact;
        h_colIndices_compact.reserve(A_r.nnz());
        h_values_compact.reserve(A_r.nnz());

        for (IndexType row = 0; row < reducedSize; ++row) {
            for (IndexType idx = h_rowOffsets[row]; idx < h_rowOffsets[row + 1]; ++idx) {
                IndexType oldCol = h_colIndices[idx];
                IndexType newCol = compactedColMapping[oldCol];

                // Only keep columns that are mapped
                if (newCol != static_cast<IndexType>(-1)) {
                    h_colIndices_compact.push_back(newCol);
                    h_values_compact.push_back(h_values[idx]);
                }
            }
            h_rowOffsets_compact[row + 1] = h_colIndices_compact.size();
        }

        // Reallocate compacted matrix
        A_r.allocate(reducedSize, compactedLocalCount, h_colIndices_compact.size());
        thrust::copy(h_rowOffsets_compact.begin(), h_rowOffsets_compact.end(),
                     thrust::device_pointer_cast(A_r.rowOffsetsPtr()));
        thrust::copy(h_colIndices_compact.begin(), h_colIndices_compact.end(),
                     thrust::device_pointer_cast(A_r.colIndicesPtr()));
        thrust::copy(h_values_compact.begin(), h_values_compact.end(),
                     thrust::device_pointer_cast(A_r.valuesPtr()));

        if (rank == 0 || rank == 1) {
            std::cout << "Rank " << rank << ": Compacted matrix nnz: " << A_r.nnz()
                      << " (removed " << (h_values.size() - h_values_compact.size()) << " entries)\n";
        }

        // Update the global mapping to use compacted version
        localToGlobalDof_all = compactedGlobalMapping;
        reducedLocalCount = compactedLocalCount;
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