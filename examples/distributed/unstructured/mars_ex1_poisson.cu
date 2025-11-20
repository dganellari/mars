// MARS Poisson Example - GPU-native finite element solver
// Solves: -Δu = f in Ω, u = 0 on ∂Ω
// 
// Equivalent to MFEM examples/ex1.cpp but fully GPU-accelerated
//
// Compile:
//   nvcc -std=c++17 mars_ex1_poisson.cu -o mars_ex1_poisson \
//        -I/path/to/mars -I/path/to/cornerstone \
//        -lcusparse -lcublas
//
// Run:
//   mpirun -np 4 ./mars_ex1_poisson --mesh mesh_parts --order 1

#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_fem.hpp"
#include <mpi.h>
#include <iostream>
#include <chrono>
#include <cmath>

using namespace mars;
using namespace mars::fem;

// Source term: f(x,y,z) for RHS
struct SourceTerm {
    __device__ __host__
    float operator()(float x, float y, float z) const {
        // Constant source: f = 1
        return 1.0f;
    }
};

// Exact solution for verification (if known)
struct ExactSolution {
    __device__ __host__
    float operator()(float x, float y, float z) const {
        // For f=1 on unit cube with u=0 on boundary
        // Exact solution depends on domain geometry
        return 0.0f;  // Placeholder
    }
};

int main(int argc, char* argv[]) {
    // Initialize MPI
    MPI_Init(&argc, &argv);
    
    int rank, numRanks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
    
    // Parse command line arguments
    std::string meshPath = "mesh_parts";
    int order = 1;
    int maxIter = 500;
    float tolerance = 1e-6f;
    
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--mesh" && i + 1 < argc) {
            meshPath = argv[++i];
        } else if (arg == "--order" && i + 1 < argc) {
            order = std::atoi(argv[++i]);
        } else if (arg == "--max-iter" && i + 1 < argc) {
            maxIter = std::atoi(argv[++i]);
        } else if (arg == "--tol" && i + 1 < argc) {
            tolerance = std::atof(argv[++i]);
        } else if (arg == "--help" || arg == "-h") {
            if (rank == 0) {
                std::cout << "Usage: " << argv[0] << " [options]\n"
                          << "Options:\n"
                          << "  --mesh <path>      Path to mesh files (default: mesh_parts)\n"
                          << "  --order <n>        Polynomial order (default: 1)\n"
                          << "  --max-iter <n>     Max CG iterations (default: 1000)\n"
                          << "  --tol <val>        CG tolerance (default: 1e-10)\n"
                          << "  --help, -h         Print this help message\n";
            }
            MPI_Finalize();
            return 0;
        }
    }
    
    if (rank == 0) {
        std::cout << "========================================\n"
                  << "   MARS Poisson Example (GPU-native)\n"
                  << "========================================\n"
                  << "Problem: -Δu = f in Ω, u = 0 on ∂Ω\n"
                  << "Mesh: " << meshPath << "\n"
                  << "MPI ranks: " << numRanks << "\n"
                  << "Order: " << order << "\n"
                  << "========================================\n\n";
    }
    
    try {
        auto t_total_start = std::chrono::high_resolution_clock::now();
        
        // =====================================================
        // 1. Load mesh and create domain
        // =====================================================
        if (rank == 0) std::cout << "1. Loading mesh...\n";
        auto t_mesh_start = std::chrono::high_resolution_clock::now();
        
        using Domain = ElementDomain<TetTag, float, uint64_t, cstone::GpuTag>;
        Domain domain(meshPath, rank, numRanks);
        
        auto t_mesh_end = std::chrono::high_resolution_clock::now();
        double t_mesh = std::chrono::duration<double>(t_mesh_end - t_mesh_start).count();
        
        size_t localElements = domain.localElementCount();
        size_t localNodes = domain.getNodeCount();
        
        if (rank == 0) {
            std::cout << "   Mesh loaded in " << t_mesh << " seconds\n"
                      << "   Local elements: " << localElements << "\n"
                      << "   Local nodes: " << localNodes << "\n\n";
        }
        
        // =====================================================
        // 2. Create finite element space
        // =====================================================
        if (rank == 0) std::cout << "2. Creating finite element space...\n";
        auto t_fes_start = std::chrono::high_resolution_clock::now();
        
        TetFESpace<float, uint64_t> fes(domain, order);
        size_t numDofs = fes.numDofs();
        
        auto t_fes_end = std::chrono::high_resolution_clock::now();
        double t_fes = std::chrono::duration<double>(t_fes_end - t_fes_start).count();
        
        if (rank == 0) {
            std::cout << "   FE space created in " << t_fes << " seconds\n"
                      << "   Total DOFs: " << numDofs << "\n\n";
        }
        
        // =====================================================
        // 3. Check mesh quality
        // =====================================================
        if (rank == 0) std::cout << "3. Checking mesh quality...\n";
        TetStiffnessAssembler<float, uint64_t> stiffnessAssembler;
        if (rank == 0) {
            stiffnessAssembler.checkMeshQuality(fes);
            std::cout << "\n";
        }
        
        // =====================================================
        // 4. Assemble stiffness matrix
        // =====================================================
        if (rank == 0) std::cout << "4. Assembling stiffness matrix...\n";
        auto t_stiff_start = std::chrono::high_resolution_clock::now();
        
        TetSparseMatrix<float, uint64_t> K;
        stiffnessAssembler.assemble(fes, K);
        
        auto t_stiff_end = std::chrono::high_resolution_clock::now();
        double t_stiff = std::chrono::duration<double>(t_stiff_end - t_stiff_start).count();
        
        if (rank == 0) {
            std::cout << "   Stiffness matrix assembled in " << t_stiff << " seconds\n"
                      << "   Matrix size: " << K.numRows() << " x " << K.numCols() << "\n"
                      << "   Non-zeros: " << K.nnz() << "\n\n";
        }
        
        // =====================================================
        // 5. Assemble RHS vector
        // =====================================================
        if (rank == 0) std::cout << "5. Assembling RHS vector...\n";
        auto t_rhs_start = std::chrono::high_resolution_clock::now();
        
        TetMassAssembler<float, uint64_t> massAssembler;
        cstone::DeviceVector<float> b;
        
        SourceTerm f;
        massAssembler.assembleRHS(fes, b, f);
        
        auto t_rhs_end = std::chrono::high_resolution_clock::now();
        double t_rhs = std::chrono::duration<double>(t_rhs_end - t_rhs_start).count();
        
        if (rank == 0) {
            std::cout << "   RHS vector assembled in " << t_rhs << " seconds\n";
            
            // Check RHS norm before BCs
            std::vector<float> h_b_before(b.size());
            thrust::copy(thrust::device_pointer_cast(b.data()),
                        thrust::device_pointer_cast(b.data() + b.size()),
                        h_b_before.begin());
            float b_norm_before = 0.0f;
            for (float val : h_b_before) b_norm_before += val * val;
            b_norm_before = std::sqrt(b_norm_before);
            std::cout << "   ||b|| before BCs: " << b_norm_before << "\n\n";
        }
        
        debug::checkVector("RHS vector b", b);
        
        // =====================================================
        // 6. Apply boundary conditions and form linear system
        // =====================================================
        if (rank == 0) std::cout << "6. Applying boundary conditions...\n";
        auto t_bc_start = std::chrono::high_resolution_clock::now();
        
        // Use geometric boundary detection for beam-tet mesh
        // Mesh domain: [0, 8] x [0, 1] x [0, 1]
        std::vector<bool> isBoundaryDOF(fes.numDofs(), false);
        std::vector<float> boundaryValues(fes.numDofs(), 0.0f);
        
        // Get node coordinates from domain (coordinate cache is initialized lazily)
        auto& domain_ref = fes.domain();
        
        std::vector<float> h_x(domain_ref.getNodeCount());
        std::vector<float> h_y(domain_ref.getNodeCount());
        std::vector<float> h_z(domain_ref.getNodeCount());
        
        const auto& d_x = domain_ref.getNodeX();
        const auto& d_y = domain_ref.getNodeY();
        const auto& d_z = domain_ref.getNodeZ();
        
        thrust::copy(thrust::device_pointer_cast(d_x.data()),
                    thrust::device_pointer_cast(d_x.data() + d_x.size()),
                    h_x.begin());
        thrust::copy(thrust::device_pointer_cast(d_y.data()),
                    thrust::device_pointer_cast(d_y.data() + d_y.size()),
                    h_y.begin());
        thrust::copy(thrust::device_pointer_cast(d_z.data()),
                    thrust::device_pointer_cast(d_z.data() + d_z.size()),
                    h_z.begin());
        
        // Mark boundary DOFs geometrically (same as MFEM's ess_tdof_list)
        size_t numBoundary = 0;
        const float tol = 1e-6f;
        for (size_t i = 0; i < fes.numDofs(); ++i) {
            float x = h_x[i];
            float y = h_y[i];
            float z = h_z[i];
            
            // Check if on any face of [0,8] x [0,1] x [0,1] domain
            bool onBoundary = (std::abs(x - 0.0f) < tol) || (std::abs(x - 8.0f) < tol) ||
                             (std::abs(y - 0.0f) < tol) || (std::abs(y - 1.0f) < tol) ||
                             (std::abs(z - 0.0f) < tol) || (std::abs(z - 1.0f) < tol);
            
            if (onBoundary) {
                isBoundaryDOF[i] = true;
                boundaryValues[i] = 0.0f;  // u = 0 on boundary
                numBoundary++;
            }
        }
        
        if (rank == 0) {
            std::cout << "   Total DOFs: " << fes.numDofs() << "\n";
            std::cout << "   Boundary DOFs (geometric): " << numBoundary << "\n";
            std::cout << "   Interior DOFs: " << (fes.numDofs() - numBoundary) << "\n";
        }
        
        // Form linear system with BCs (like MFEM's FormLinearSystem)
        // This eliminates boundary DOFs: A_full, b_full -> A_int, b_int
        fem::DOFElimination<float, uint64_t, cstone::GpuTag> eliminator;
        fem::SparseMatrix<uint64_t, float, cstone::GpuTag> K_int;
        cstone::DeviceVector<float> b_int;
        
        eliminator.buildInteriorSystem(K, b, isBoundaryDOF, boundaryValues, K_int, b_int);
        
        auto t_bc_end = std::chrono::high_resolution_clock::now();
        double t_bc = std::chrono::duration<double>(t_bc_end - t_bc_start).count();
        
        if (rank == 0) {
            std::cout << "   DOF elimination completed in " << t_bc << " seconds\n\n";
        }
        
        debug::checkMatrix("Interior stiffness matrix K_int", K_int);
        debug::checkVector("Interior RHS vector b_int", b_int);
        
        // Check matrix symmetry
        if (rank == 0 && K_int.numRows() > 0) {
            std::cout << "   Checking matrix symmetry...\n";
            
            // Test: create random vector, compute Av, check if max|A[i,j] - A[j,i]| is small
            std::vector<uint64_t> h_rowOffsets(K_int.numRows() + 1);
            std::vector<uint64_t> h_colIndices(K_int.nnz());
            std::vector<float> h_values(K_int.nnz());
            
            thrust::copy(thrust::device_pointer_cast(K_int.rowOffsetsPtr()),
                        thrust::device_pointer_cast(K_int.rowOffsetsPtr() + K_int.numRows() + 1),
                        h_rowOffsets.begin());
            thrust::copy(thrust::device_pointer_cast(K_int.colIndicesPtr()),
                        thrust::device_pointer_cast(K_int.colIndicesPtr() + K_int.nnz()),
                        h_colIndices.begin());
            thrust::copy(thrust::device_pointer_cast(K_int.valuesPtr()),
                        thrust::device_pointer_cast(K_int.valuesPtr() + K_int.nnz()),
                        h_values.begin());
            
            // Check a few entries for symmetry
            int asymmetries = 0;
            float maxAsymmetry = 0.0f;
            for (size_t i = 0; i < std::min<size_t>(100, K_int.numRows()); ++i) {
                for (uint64_t idx = h_rowOffsets[i]; idx < h_rowOffsets[i + 1]; ++idx) {
                    uint64_t j = h_colIndices[idx];
                    float a_ij = h_values[idx];
                    
                    // Find A[j,i]
                    float a_ji = 0.0f;
                    bool found = false;
                    for (uint64_t idx2 = h_rowOffsets[j]; idx2 < h_rowOffsets[j + 1]; ++idx2) {
                        if (h_colIndices[idx2] == i) {
                            a_ji = h_values[idx2];
                            found = true;
                            break;
                        }
                    }
                    
                    if (found && std::abs(a_ij - a_ji) > 1e-6) {
                        asymmetries++;
                        maxAsymmetry = std::max(maxAsymmetry, std::abs(a_ij - a_ji));
                    }
                }
            }
            
            std::cout << "   Asymmetries found: " << asymmetries 
                      << ", max |A[i,j] - A[j,i]| = " << maxAsymmetry << "\n\n";
        }
        
        // =====================================================
        // 7. Solve reduced interior system
        // =====================================================
        if (rank == 0) std::cout << "7. Solving reduced linear system (GMRES)...\n";
        auto t_solve_start = std::chrono::high_resolution_clock::now();
        
        // Use GMRES with larger restart and more iterations
        TetGMRESSolver<float, uint64_t> solver(5000, 1e-4f, 100);
        solver.setVerbose(rank == 0);
        
        // Solve interior system (without preconditioner first)
        cstone::DeviceVector<float> u_int;
        bool converged = solver.solve(K_int, b_int, u_int, false);
        
        auto t_solve_end = std::chrono::high_resolution_clock::now();
        double t_solve = std::chrono::duration<double>(t_solve_end - t_solve_start).count();
        
        if (rank == 0) {
            std::cout << "   System solved in " << t_solve << " seconds\n"
                      << "   Converged: " << (converged ? "Yes" : "No") << "\n\n";
        }
        
        // =====================================================
        // 8. Reconstruct full solution
        // =====================================================
        if (rank == 0) std::cout << "8. Reconstructing full solution...\n";
        
        cstone::DeviceVector<float> u;
        eliminator.reconstructFullSolution(u_int, isBoundaryDOF, boundaryValues, u);
        
        if (rank == 0) std::cout << "   Full solution reconstructed\n\n";
        
        debug::checkVector("Full solution u", u);
        
        // =====================================================
        // 9. Compute solution statistics
        // =====================================================
        if (rank == 0) std::cout << "9. Computing solution statistics...\n";
        
        // Copy solution to host for analysis
        std::vector<float> h_u(u.size());
        thrust::copy(thrust::device_pointer_cast(u.data()), 
                    thrust::device_pointer_cast(u.data() + u.size()), 
                    h_u.begin());
        
        float u_min = *std::min_element(h_u.begin(), h_u.end());
        float u_max = *std::max_element(h_u.begin(), h_u.end());
        float u_sum = std::accumulate(h_u.begin(), h_u.end(), 0.0f);
        float u_mean = u_sum / h_u.size();
        
        // Compute L2 norm
        float l2_norm = 0.0f;
        for (float val : h_u) {
            l2_norm += val * val;
        }
        l2_norm = std::sqrt(l2_norm);
        
        if (rank == 0) {
            std::cout << "   L2 norm: " << l2_norm << "\n\n";
        }
        
        // =====================================================
        // 8. Write VTK output
        // =====================================================
        // if (rank == 0) std::cout << "8. Writing VTK output...\n";
        // 
        // VTKWriter writer;
        // writer.writeVTU("mars_poisson_solution.vtu", domain, h_u, "solution");
        // 
        // if (rank == 0) {
        //     std::cout << "   VTK file written to mars_poisson_solution.vtu\n\n";
        // }
        
        // =====================================================
        // 9. Timing summary
        // =====================================================
        
        if (rank == 0) {
            std::cout << "   Solution statistics:\n"
                      << "     Min: " << u_min << "\n"
                      << "     Max: " << u_max << "\n"
                      << "     Mean: " << u_mean << "\n"
                      << "     L2 norm: " << l2_norm << "\n\n";
        }
        
        // =====================================================
        // 8. Write VTK output
        // =====================================================
        // if (rank == 0) std::cout << "8. Writing VTK output...\n";
        // 
        // VTKWriter vtkWriter;
        // vtkWriter.writeVTU("mars_poisson_solution.vtu", domain, h_u, "solution");
        // 
        // if (rank == 0) {
        //     std::cout << "   VTK file written to mars_poisson_solution.vtu\n\n";
        // }
        
        // =====================================================
        // 10. Timing summary
        // =====================================================
        auto t_total_end = std::chrono::high_resolution_clock::now();
        double t_total = std::chrono::duration<double>(t_total_end - t_total_start).count();
        
        if (rank == 0) {
            std::cout << "========================================\n"
                      << "   Timing Summary\n"
                      << "========================================\n"
                      << "Mesh loading:     " << t_mesh << " s\n"
                      << "FE space setup:   " << t_fes << " s\n"
                      << "Stiffness assembly: " << t_stiff << " s\n"
                      << "RHS assembly:     " << t_rhs << " s\n"
                      << "BC application:   " << t_bc << " s\n"
                      << "Linear solve:     " << t_solve << " s\n"
                      << "----------------------------------------\n"
                      << "Total time:       " << t_total << " s\n"
                      << "========================================\n\n";
            
            std::cout << "MARS Poisson example completed successfully!\n";
        }
        
    } catch (const std::exception& e) {
        if (rank == 0) {
            std::cerr << "Error: " << e.what() << std::endl;
        }
        MPI_Finalize();
        return 1;
    }
    
    MPI_Finalize();
    return 0;
}
