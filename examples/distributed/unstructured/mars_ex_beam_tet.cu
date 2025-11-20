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
#include "backend/distributed/unstructured/utils/mars_read_mfem_mesh.hpp"

#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <chrono>

using namespace mars;
using namespace mars::fem;

// Source term: f(x,y,z) = 1
struct SourceTerm {
    __device__ __host__
    float operator()(float x, float y, float z) const {
        return 1.0f;
    }
};

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
    
    // Auto-convert MFEM mesh to binary if needed
    std::string binaryMeshPath = meshPath;
    if (meshPath.find(".mesh") != std::string::npos) {
        // This is an MFEM mesh file - convert to binary first
        binaryMeshPath = meshPath.substr(0, meshPath.find_last_of(".")) + "-binary";
        
        if (rank == 0) {
            std::cout << "\nDetected MFEM format: " << meshPath << std::endl;
            
            // Check if binary version already exists
            std::ifstream check(binaryMeshPath + "/x.float32");
            if (!check.good()) {
                std::cout << "Converting to binary format (3 refinement levels)..." << std::endl;
                
                // Read MFEM mesh WITHOUT partitioning - just get the raw global mesh
                // We'll let ElementDomain handle partitioning when it reads the binary
                std::ifstream file(meshPath);
                if (!file) {
                    std::cerr << "Failed to open MFEM mesh file: " << meshPath << std::endl;
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                
                // Parse MFEM header
                std::string line;
                int dimension = 0;
                size_t num_elements = 0;
                size_t num_vertices = 0;
                
                while (std::getline(file, line)) {
                    if (line.find("dimension") != std::string::npos) {
                        file >> dimension;
                        std::getline(file, line);
                        break;
                    }
                }
                
                while (std::getline(file, line)) {
                    if (line.find("elements") != std::string::npos) {
                        file >> num_elements;
                        std::getline(file, line);
                        break;
                    }
                }
                
                // Read elements
                std::vector<unsigned> i0, i1, i2, i3;
                i0.reserve(num_elements);
                i1.reserve(num_elements);
                i2.reserve(num_elements);
                i3.reserve(num_elements);
                
                for (size_t i = 0; i < num_elements; i++) {
                    int attr, geom_type;
                    unsigned n0, n1, n2, n3;
                    file >> attr >> geom_type >> n0 >> n1 >> n2 >> n3;
                    i0.push_back(n0);
                    i1.push_back(n1);
                    i2.push_back(n2);
                    i3.push_back(n3);
                }
                
                // Skip boundary section
                while (std::getline(file, line)) {
                    if (line.find("boundary") != std::string::npos) {
                        size_t num_boundary;
                        file >> num_boundary;
                        std::getline(file, line);
                        for (size_t i = 0; i < num_boundary; i++) {
                            std::getline(file, line);
                        }
                        break;
                    }
                }
                
                // Read vertices
                while (std::getline(file, line)) {
                    if (line.find("vertices") != std::string::npos) {
                        file >> num_vertices;
                        std::getline(file, line);
                        break;
                    }
                }
                
                int vertex_dim;
                file >> vertex_dim;
                std::getline(file, line);
                
                std::vector<float> x_data, y_data, z_data;
                x_data.reserve(num_vertices);
                y_data.reserve(num_vertices);
                z_data.reserve(num_vertices);
                
                for (size_t i = 0; i < num_vertices; i++) {
                    double x, y, z = 0.0;
                    file >> x >> y;
                    if (dimension == 3) file >> z;
                    x_data.push_back(static_cast<float>(x));
                    y_data.push_back(static_cast<float>(y));
                    z_data.push_back(static_cast<float>(z));
                }
                file.close();
                
                auto connectivity = std::make_tuple(std::move(i0), std::move(i1), std::move(i2), std::move(i3));
                size_t nodeCount = x_data.size();
                size_t elementCount = std::get<0>(connectivity).size();
                
                std::cout << "Initial mesh: " << elementCount << " elements, " << nodeCount << " nodes" << std::endl;
                
                // Apply 3 refinement levels
                for (int l = 0; l < 3; l++) {
                    uniformRefineMFEMMesh<4, float, unsigned>(x_data, y_data, z_data, connectivity, nodeCount, elementCount);
                }
                std::cout << "Refined mesh: " << elementCount << " elements, " << nodeCount << " nodes" << std::endl;
                
                // Write binary format
                std::system(("mkdir -p " + binaryMeshPath).c_str());
                
                auto writeFloatBinary = [](const std::string& filename, const std::vector<float>& data) {
                    std::ofstream file(filename, std::ios::binary);
                    file.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(float));
                };
                
                auto writeIntBinary = [](const std::string& filename, const std::vector<unsigned>& data) {
                    std::ofstream file(filename, std::ios::binary);
                    std::vector<int32_t> int32_data(data.begin(), data.end());
                    file.write(reinterpret_cast<const char*>(int32_data.data()), int32_data.size() * sizeof(int32_t));
                };
                
                writeFloatBinary(binaryMeshPath + "/x.float32", x_data);
                writeFloatBinary(binaryMeshPath + "/y.float32", y_data);
                writeFloatBinary(binaryMeshPath + "/z.float32", z_data);
                writeIntBinary(binaryMeshPath + "/i0.int32", std::get<0>(connectivity));
                writeIntBinary(binaryMeshPath + "/i1.int32", std::get<1>(connectivity));
                writeIntBinary(binaryMeshPath + "/i2.int32", std::get<2>(connectivity));
                writeIntBinary(binaryMeshPath + "/i3.int32", std::get<3>(connectivity));
                
                std::cout << "Conversion complete: " << binaryMeshPath << std::endl;
            } else {
                std::cout << "Using cached binary mesh: " << binaryMeshPath << std::endl;
            }
        }
        
        // Wait for rank 0
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    if (rank == 0) {
        std::cout << "\nLoading mesh from: " << binaryMeshPath << std::endl;
    }
    
    // Create ElementDomain using standard file-based constructor
    using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
    Domain domain(binaryMeshPath, rank, numRanks);
    // Note: Adjacency is built lazily when needed by FE space
    
    if (rank == 0) {
        std::cout << "\n=== ElementDomain Created Successfully ===" << std::endl;
        std::cout << "Domain info:" << std::endl;
        std::cout << "  Total nodes: " << domain.getNodeCount() << std::endl;
        std::cout << "  Total elements: " << domain.getElementCount() << std::endl;
        std::cout << "  Local elements: " << domain.localElementCount() << std::endl;
        std::cout << "\n=== Starting FEM Solve ===" << std::endl;
    }
    
    auto t_total_start = std::chrono::high_resolution_clock::now();
    
    // =====================================================
    // 1. Create finite element space
    // =====================================================
    if (rank == 0) std::cout << "\n1. Creating finite element space...\n";
    auto t_fes_start = std::chrono::high_resolution_clock::now();
    
    TetFESpace<float, unsigned> fes(domain, order);
    size_t numDofs = fes.numDofs();
    
    auto t_fes_end = std::chrono::high_resolution_clock::now();
    double t_fes = std::chrono::duration<double>(t_fes_end - t_fes_start).count();
    
    if (rank == 0) {
        std::cout << "   FE space created in " << t_fes << " seconds\n"
                  << "   Total DOFs: " << numDofs << "\n";
    }
    
    // Debug: Check element connectivity
    if (rank == 0) {
        std::cout << "\nDebug: Checking element connectivity...\n";
        std::cout << "   Domain nodes: " << domain.getNodeCount() << "\n";
        std::cout << "   Domain elements: " << domain.getElementCount() << "\n";
        std::cout << "   Local elements: " << domain.localElementCount() << "\n";
    }
    
    // =====================================================
    // 2. Assemble stiffness matrix
    // =====================================================
    if (rank == 0) std::cout << "\n2. Assembling stiffness matrix...\n";
    auto t_stiff_start = std::chrono::high_resolution_clock::now();
    
    TetStiffnessAssembler<float, unsigned> stiffnessAssembler;
    TetSparseMatrix<float, unsigned> K;
    stiffnessAssembler.assemble(fes, K);
    
    auto t_stiff_end = std::chrono::high_resolution_clock::now();
    double t_stiff = std::chrono::duration<double>(t_stiff_end - t_stiff_start).count();
    
    if (rank == 0) {
        std::cout << "   Stiffness matrix assembled in " << t_stiff << " seconds\n"
                  << "   Matrix size: " << K.numRows() << " x " << K.numCols() << "\n"
                  << "   Non-zeros: " << K.nnz() << "\n";
        
        // Check matrix symmetry
        std::cout << "   Checking matrix properties...\n";
        std::vector<unsigned> h_rowOffsets(K.numRows() + 1);
        std::vector<unsigned> h_colIndices(K.nnz());
        std::vector<float> h_values(K.nnz());
        
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
        float minDiag = 1e30f, maxDiag = -1e30f;
        int zeroDiag = 0, negDiag = 0;
        for (size_t i = 0; i < K.numRows(); ++i) {
            float diagVal = 0.0f;
            bool foundDiag = false;
            for (unsigned idx = h_rowOffsets[i]; idx < h_rowOffsets[i+1]; ++idx) {
                if (h_colIndices[idx] == i) {
                    diagVal = h_values[idx];
                    foundDiag = true;
                    break;
                }
            }
            if (!foundDiag) zeroDiag++;
            else {
                if (diagVal <= 0.0f) negDiag++;
                minDiag = std::min(minDiag, diagVal);
                maxDiag = std::max(maxDiag, diagVal);
            }
        }
        std::cout << "   Diagonal range: [" << minDiag << ", " << maxDiag << "]\n";
        std::cout << "   Zero diagonal entries: " << zeroDiag << "\n";
        std::cout << "   Negative diagonal entries: " << negDiag << "\n";
    }
    
    // =====================================================
    // 3. Assemble RHS vector
    // =====================================================
    if (rank == 0) std::cout << "\n3. Assembling RHS vector...\n";
    auto t_rhs_start = std::chrono::high_resolution_clock::now();
    
    TetMassAssembler<float, unsigned> massAssembler;
    cstone::DeviceVector<float> b;
    
    SourceTerm f;
    massAssembler.assembleRHS(fes, b, f);
    
    auto t_rhs_end = std::chrono::high_resolution_clock::now();
    double t_rhs = std::chrono::duration<double>(t_rhs_end - t_rhs_start).count();
    
    if (rank == 0) {
        std::cout << "   RHS vector assembled in " << t_rhs << " seconds\n";
    }
    
    // =====================================================
    // 4. Apply boundary conditions (GEOMETRIC detection)
    // =====================================================
    if (rank == 0) std::cout << "\n4. Applying boundary conditions (geometric detection)...\n";
    auto t_bc_start = std::chrono::high_resolution_clock::now();
    
    // Use geometric boundary detection for beam-tet mesh
    // Apply Dirichlet BCs: u = 0 on all 6 faces
    std::vector<bool> isBoundaryDOF(fes.numDofs(), false);
    std::vector<float> boundaryValues(fes.numDofs(), 0.0f);
    
    // Get ownership map to identify owned vs ghost nodes
    const auto& ownership = domain.getNodeOwnershipMap();
    thrust::host_vector<uint8_t> h_ownership(domain.getNodeCount());
    thrust::copy(thrust::device_pointer_cast(ownership.data()),
                 thrust::device_pointer_cast(ownership.data() + domain.getNodeCount()),
                 h_ownership.begin());
    
    // Create mapping from owned node indices to DOF indices
    std::vector<size_t> ownedNodeToDof(domain.getNodeCount(), SIZE_MAX);
    size_t dofIdx = 0;
    for (size_t nodeIdx = 0; nodeIdx < domain.getNodeCount(); ++nodeIdx) {
        if (h_ownership[nodeIdx] == 1) {
            ownedNodeToDof[nodeIdx] = dofIdx++;
        }
    }
    
    // Get node coordinates from domain
    const auto& d_x = domain.getNodeX();
    const auto& d_y = domain.getNodeY();
    const auto& d_z = domain.getNodeZ();
    
    thrust::host_vector<float> h_x(domain.getNodeCount());
    thrust::host_vector<float> h_y(domain.getNodeCount());
    thrust::host_vector<float> h_z(domain.getNodeCount());
    
    thrust::copy(thrust::device_pointer_cast(d_x.data()),
                thrust::device_pointer_cast(d_x.data() + d_x.size()),
                h_x.begin());
    thrust::copy(thrust::device_pointer_cast(d_y.data()),
                thrust::device_pointer_cast(d_y.data() + d_y.size()),
                h_y.begin());
    thrust::copy(thrust::device_pointer_cast(d_z.data()),
                thrust::device_pointer_cast(d_z.data() + d_z.size()),
                h_z.begin());
    
    // First, find coordinate bounds on this rank
    float x_min = 1e30f, x_max = -1e30f;
    float y_min = 1e30f, y_max = -1e30f;
    float z_min = 1e30f, z_max = -1e30f;
    for (size_t i = 0; i < domain.getNodeCount(); ++i) {
        x_min = std::min(x_min, h_x[i]);
        x_max = std::max(x_max, h_x[i]);
        y_min = std::min(y_min, h_y[i]);
        y_max = std::max(y_max, h_y[i]);
        z_min = std::min(z_min, h_z[i]);
        z_max = std::max(z_max, h_z[i]);
    }
    
    // Get global bounds across all ranks for consistent boundary detection
    float global_x_min, global_x_max, global_y_min, global_y_max, global_z_min, global_z_max;
    MPI_Allreduce(&x_min, &global_x_min, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&x_max, &global_x_max, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&y_min, &global_y_min, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&y_max, &global_y_max, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&z_min, &global_z_min, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&z_max, &global_z_max, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    
    if (rank == 0) {
        std::cout << "   Global mesh bounds: x=[" << global_x_min << ", " << global_x_max << "], "
                  << "y=[" << global_y_min << ", " << global_y_max << "], "
                  << "z=[" << global_z_min << ", " << global_z_max << "]\n";
    }
    
    const float tol = 1e-4f;
    size_t numBoundary = 0;
    
    // Iterate over all local nodes, mark boundary conditions for owned nodes only
    for (size_t nodeIdx = 0; nodeIdx < domain.getNodeCount(); ++nodeIdx) {
        if (h_ownership[nodeIdx] != 1) continue;  // Skip ghost nodes
        
        float x = h_x[nodeIdx];
        float y = h_y[nodeIdx];
        float z = h_z[nodeIdx];
        
        // Check if on any face of the domain (using global bounds for consistency)
        bool onBoundary = (std::abs(x - global_x_min) < tol) || (std::abs(x - global_x_max) < tol) ||
                         (std::abs(y - global_y_min) < tol) || (std::abs(y - global_y_max) < tol) ||
                         (std::abs(z - global_z_min) < tol) || (std::abs(z - global_z_max) < tol);
        
        if (onBoundary) {
            size_t dof = ownedNodeToDof[nodeIdx];
            isBoundaryDOF[dof] = true;
            boundaryValues[dof] = 0.0f;
            numBoundary++;
        }
    }
    
    if (rank == 0) {
        std::cout << "   Total DOFs: " << fes.numDofs() << "\n";
        std::cout << "   Boundary DOFs (geometric): " << numBoundary << "\n";
        std::cout << "   Interior DOFs: " << (fes.numDofs() - numBoundary) << "\n";
    }
    
    // Eliminate boundary DOFs
    fem::DOFElimination<float, unsigned, cstone::GpuTag> eliminator;
    fem::SparseMatrix<unsigned, float, cstone::GpuTag> K_int;
    cstone::DeviceVector<float> b_int;
    
    eliminator.buildInteriorSystem(K, b, isBoundaryDOF, boundaryValues, K_int, b_int);
    
    auto t_bc_end = std::chrono::high_resolution_clock::now();
    double t_bc = std::chrono::duration<double>(t_bc_end - t_bc_start).count();
    
    if (rank == 0) {
        std::cout << "   DOF elimination completed in " << t_bc << " seconds\n";
        std::cout << "   Interior DOFs: " << K_int.numRows() << "\n";
    }
    
    // =====================================================
    // 5. Solve linear system (PCG with Gauss-Seidel preconditioner - matching MFEM ex1)
    // =====================================================
    if (rank == 0) std::cout << "\n5. Solving linear system (PCG + Jacobi)...\n";
    auto t_solve_start = std::chrono::high_resolution_clock::now();
    
    cstone::DeviceVector<float> u_int(K_int.numRows(), 0.0f);
    
    // CG solver with tolerance and max iterations matching MFEM ex1.cpp
    ConjugateGradientSolver<float, unsigned, cstone::GpuTag> cg(400, 1e-6f);
    cg.setVerbose(rank == 0);  // Only rank 0 prints
    
    bool converged = cg.solve(K_int, b_int, u_int);
    
    auto t_solve_end = std::chrono::high_resolution_clock::now();
    double t_solve = std::chrono::duration<double>(t_solve_end - t_solve_start).count();
    
    if (rank == 0) {
        std::cout << "   Solver " << (converged ? "converged" : "did NOT converge") 
                  << " in " << t_solve << " seconds\n";
    }
    
    // =====================================================
    // 6. Reconstruct full solution from interior solution
    // =====================================================
    cstone::DeviceVector<float> u_full(fes.numDofs(), 0.0f);
    eliminator.reconstructFullSolution(u_int, isBoundaryDOF, boundaryValues, u_full);
    
    // =====================================================
    // 7. Exchange ghost DOF values across ranks
    // =====================================================
    if (numRanks > 1) {
        if (rank == 0) std::cout << "\n6. Exchanging ghost DOF values...\n";
        auto t_exchange_start = std::chrono::high_resolution_clock::now();
        
        // Create DOF handler for communication
        UnstructuredDofHandler<TetTag, float, unsigned, cstone::GpuTag> dofHandler(domain, rank, numRanks);
        dofHandler.initialize();
        
        // Exchange solution values at shared nodes
        dofHandler.exchangeNodeData(u_full);
        
        auto t_exchange_end = std::chrono::high_resolution_clock::now();
        double t_exchange = std::chrono::duration<double>(t_exchange_end - t_exchange_start).count();
        
        if (rank == 0) {
            std::cout << "   Ghost DOF exchange completed in " << t_exchange << " seconds\n";
        }
    }
    
    // Compute solution statistics (after exchange)
    std::vector<float> h_u(u_full.size());
    thrust::copy(thrust::device_pointer_cast(u_full.data()),
                thrust::device_pointer_cast(u_full.data() + u_full.size()),
                h_u.begin());
    
    float u_min = *std::min_element(h_u.begin(), h_u.end());
    float u_max = *std::max_element(h_u.begin(), h_u.end());
    float u_norm = 0.0f;
    for (float val : h_u) u_norm += val * val;
    u_norm = std::sqrt(u_norm);
    
    auto t_total_end = std::chrono::high_resolution_clock::now();
    double t_total = std::chrono::duration<double>(t_total_end - t_total_start).count();
    
    if (rank == 0) {
        std::cout << "\n=== Solution Statistics ===" << std::endl;
        std::cout << "  min(u) = " << u_min << std::endl;
        std::cout << "  max(u) = " << u_max << std::endl;
        std::cout << "  ||u||  = " << u_norm << std::endl;
        
        std::cout << "\n=== Timing Summary ===" << std::endl;
        std::cout << "  FE space:     " << t_fes << " s" << std::endl;
        std::cout << "  Assembly:     " << (t_stiff + t_rhs) << " s" << std::endl;
        std::cout << "  BCs:          " << t_bc << " s" << std::endl;
        std::cout << "  Solve:        " << t_solve << " s" << std::endl;
        std::cout << "  Total:        " << t_total << " s" << std::endl;
        std::cout << "\n=== Test Complete ===" << std::endl;
        std::cout << "Compare with MFEM: mpirun -np " << numRanks << " ex1p -m beam-tet.mesh" << std::endl;
    }
    
    MPI_Finalize();
    return 0;
} 