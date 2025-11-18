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
    int maxIter = 1000;
    float tolerance = 1e-10f;
    
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
        // 3. Assemble stiffness matrix
        // =====================================================
        if (rank == 0) std::cout << "3. Assembling stiffness matrix...\n";
        auto t_stiff_start = std::chrono::high_resolution_clock::now();
        
        TetStiffnessAssembler<float, uint64_t> stiffnessAssembler;
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
        // 4. Assemble RHS vector
        // =====================================================
        if (rank == 0) std::cout << "4. Assembling RHS vector...\n";
        auto t_rhs_start = std::chrono::high_resolution_clock::now();
        
        TetMassAssembler<float, uint64_t> massAssembler;
        cstone::DeviceVector<float> b;
        
        SourceTerm f;
        massAssembler.assembleRHS(fes, b, f);
        
        auto t_rhs_end = std::chrono::high_resolution_clock::now();
        double t_rhs = std::chrono::duration<double>(t_rhs_end - t_rhs_start).count();
        
        if (rank == 0) {
            std::cout << "   RHS vector assembled in " << t_rhs << " seconds\n\n";
        }
        
        // =====================================================
        // 5. Apply boundary conditions
        // =====================================================
        if (rank == 0) std::cout << "5. Applying boundary conditions...\n";
        auto t_bc_start = std::chrono::high_resolution_clock::now();
        
        TetBCHandler<float, uint64_t> bcHandler;
        auto boundaryDofs = fes.getBoundaryDofs();
        
        if (rank == 0) {
            std::cout << "   Boundary DOFs: " << boundaryDofs.size() << "\n";
        }
        
        bcHandler.applyDirichlet(fes, K, b, boundaryDofs, 0.0f);
        
        auto t_bc_end = std::chrono::high_resolution_clock::now();
        double t_bc = std::chrono::duration<double>(t_bc_end - t_bc_start).count();
        
        if (rank == 0) {
            std::cout << "   BCs applied in " << t_bc << " seconds\n\n";
        }
        
        // =====================================================
        // 6. Solve linear system with CG
        // =====================================================
        if (rank == 0) std::cout << "6. Solving linear system (CG)...\n";
        auto t_solve_start = std::chrono::high_resolution_clock::now();
        
        TetCGSolver<float, uint64_t> solver(maxIter, tolerance);
        solver.setVerbose(rank == 0);
        
        cstone::DeviceVector<float> u;
        bool converged = solver.solve(K, b, u);
        
        auto t_solve_end = std::chrono::high_resolution_clock::now();
        double t_solve = std::chrono::duration<double>(t_solve_end - t_solve_start).count();
        
        if (rank == 0) {
            std::cout << "   System solved in " << t_solve << " seconds\n"
                      << "   Converged: " << (converged ? "Yes" : "No") << "\n\n";
        }
        
        // =====================================================
        // 7. Compute solution statistics
        // =====================================================
        if (rank == 0) std::cout << "7. Computing solution statistics...\n";
        
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
        float u_l2_sq = 0.0f;
        for (float val : h_u) {
            u_l2_sq += val * val;
        }
        float u_l2 = std::sqrt(u_l2_sq);
        
        if (rank == 0) {
            std::cout << "   Solution statistics:\n"
                      << "     Min: " << u_min << "\n"
                      << "     Max: " << u_max << "\n"
                      << "     Mean: " << u_mean << "\n"
                      << "     L2 norm: " << u_l2 << "\n\n";
        }
        
        // =====================================================
        // 8. Export solution (VTK) - TODO
        // =====================================================
        // if (rank == 0) std::cout << "8. Exporting solution to VTK...\n";
        // VTKWriter writer;
        // writer.write("solution.vtu", domain, fes, u);
        
        // =====================================================
        // Summary
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
