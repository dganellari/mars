// Simple Poisson test matching MFEM ex0p: -Delta u = 1, u = 0 on boundary
// This is for direct comparison with MFEM to validate our implementation

#include "mars.hpp"
#include "mars_distributed_staggered_mesh_block.hpp"
#include "mars_distributed_dof_management.hpp"
#include "mars_distributed_fe_octant_kokkos.hpp"
#include "mars_err.hpp"

#include "mars_fem.hpp"

#include <mpi.h>
#include <chrono>
#include <iostream>
#include <cmath>

using namespace mars;

int main(int argc, char* argv[]) {
    // Initialize MPI
    MPI_Init(&argc, &argv);
    
    int rank, numRanks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
    
    // Problem: -Delta u = 1 with u = 0 on boundary
    // Exact solution on unit cube [0,1]^3 is not simple, but we know:
    // - u = 0 on boundary
    // - u > 0 in interior
    // - max(u) should be at center, roughly O(0.01) for unit cube
    
    if (rank == 0) {
        std::cout << "\n========================================\n";
        std::cout << "MARS Poisson Test (matching MFEM ex0p)\n";
        std::cout << "========================================\n";
        std::cout << "Problem: -Delta u = 1, u = 0 on boundary\n";
        std::cout << "MPI ranks: " << numRanks << "\n";
        std::cout << "========================================\n\n";
    }
    
    try {
        // Load mesh
        std::string meshPath = "mesh_parts";
        if (argc > 1) {
            meshPath = argv[1];
        }
        
        if (rank == 0) std::cout << "1. Loading mesh from: " << meshPath << "\n";
        auto t_mesh_start = std::chrono::high_resolution_clock::now();
        
        using Mesh = DistributedDMesh<double>;
        Mesh mesh(MPI_COMM_WORLD);
        mesh.load(meshPath, 1);
        
        auto t_mesh_end = std::chrono::high_resolution_clock::now();
        double t_mesh = std::chrono::duration<double>(t_mesh_end - t_mesh_start).count();
        
        if (rank == 0) {
            std::cout << "   Mesh loaded in " << t_mesh << " seconds\n";
            std::cout << "   Elements: " << mesh.get_mesh().get_elem_size() << "\n";
            std::cout << "   Nodes: " << mesh.get_view_size() << "\n\n";
        }
        
        // Create FE space
        if (rank == 0) std::cout << "2. Creating finite element space...\n";
        auto t_fes_start = std::chrono::high_resolution_clock::now();
        
        fem::TetFESpace<float, uint64_t> fes(mesh);
        
        auto t_fes_end = std::chrono::high_resolution_clock::now();
        double t_fes = std::chrono::duration<double>(t_fes_end - t_fes_start).count();
        
        if (rank == 0) {
            std::cout << "   FE space created in " << t_fes << " seconds\n";
            std::cout << "   DOFs: " << fes.numDofs() << "\n\n";
        }
        
        // Check mesh quality
        if (rank == 0) std::cout << "3. Checking mesh quality...\n";
        fem::TetStiffnessAssembler<float, uint64_t> tempAssembler;
        tempAssembler.checkMeshQuality(mesh);
        if (rank == 0) std::cout << "\n";
        
        // Assemble stiffness matrix
        if (rank == 0) std::cout << "4. Assembling stiffness matrix...\n";
        auto t_stiff_start = std::chrono::high_resolution_clock::now();
        
        fem::SparseMatrix<uint64_t, float, cstone::GpuTag> K;
        fem::TetStiffnessAssembler<float, uint64_t> stiffAssembler;
        stiffAssembler.assemble(fes, K);
        
        auto t_stiff_end = std::chrono::high_resolution_clock::now();
        double t_stiff = std::chrono::duration<double>(t_stiff_end - t_stiff_start).count();
        
        if (rank == 0) {
            std::cout << "   Stiffness matrix assembled in " << t_stiff << " seconds\n";
            std::cout << "   Matrix size: " << K.numRows() << " x " << K.numCols() << "\n";
            std::cout << "   Non-zeros: " << K.nnz() << "\n\n";
        }
        
        // Assemble RHS vector
        if (rank == 0) std::cout << "5. Assembling RHS vector (f=1)...\n";
        auto t_rhs_start = std::chrono::high_resolution_clock::now();
        
        cstone::DeviceVector<float> b(fes.numDofs());
        fem::TetMassAssembler<float, uint64_t> massAssembler;
        massAssembler.assembleRHS(fes, b, 1.0f);  // source term f = 1
        
        auto t_rhs_end = std::chrono::high_resolution_clock::now();
        double t_rhs = std::chrono::duration<double>(t_rhs_end - t_rhs_start).count();
        
        // Check RHS norm
        if (rank == 0) {
            std::vector<float> h_b_before(b.size());
            thrust::copy(thrust::device_pointer_cast(b.data()),
                        thrust::device_pointer_cast(b.data() + b.size()),
                        h_b_before.begin());
            float b_norm_before = 0.0f;
            for (float val : h_b_before) b_norm_before += val * val;
            b_norm_before = std::sqrt(b_norm_before);
            
            std::cout << "   RHS assembled in " << t_rhs << " seconds\n";
            std::cout << "   ||b|| = " << b_norm_before << "\n\n";
        }
        
        // Get boundary DOFs
        auto boundaryDofs = fes.getBoundaryDofs();
        
        if (rank == 0) {
            std::cout << "6. Boundary conditions...\n";
            std::cout << "   Total DOFs: " << fes.numDofs() << "\n";
            std::cout << "   Boundary DOFs: " << boundaryDofs.size() << "\n";
            std::cout << "   Interior DOFs: " << (fes.numDofs() - boundaryDofs.size()) << "\n\n";
        }
        
        // DOF elimination
        if (rank == 0) std::cout << "7. Eliminating boundary DOFs...\n";
        auto t_elim_start = std::chrono::high_resolution_clock::now();
        
        std::vector<bool> isBoundaryDOF(fes.numDofs(), false);
        std::vector<float> boundaryValues(fes.numDofs(), 0.0f);
        
        for (uint64_t dof : boundaryDofs) {
            isBoundaryDOF[dof] = true;
            boundaryValues[dof] = 0.0f;  // u = 0 on boundary
        }
        
        fem::DOFElimination<float, uint64_t, cstone::GpuTag> eliminator;
        fem::SparseMatrix<uint64_t, float, cstone::GpuTag> K_int;
        cstone::DeviceVector<float> b_int;
        
        eliminator.buildInteriorSystem(K, b, isBoundaryDOF, boundaryValues, K_int, b_int);
        
        auto t_elim_end = std::chrono::high_resolution_clock::now();
        double t_elim = std::chrono::duration<double>(t_elim_end - t_elim_start).count();
        
        if (rank == 0) {
            std::cout << "   Elimination completed in " << t_elim << " seconds\n\n";
        }
        
        // Solve with GMRES
        if (rank == 0) std::cout << "8. Solving with GMRES...\n";
        auto t_solve_start = std::chrono::high_resolution_clock::now();
        
        fem::TetGMRESSolver<float, uint64_t> solver(1000, 1e-6f, 50);
        solver.setVerbose(rank == 0);
        
        cstone::DeviceVector<float> u_int;
        bool converged = solver.solve(K_int, b_int, u_int, false);
        
        auto t_solve_end = std::chrono::high_resolution_clock::now();
        double t_solve = std::chrono::duration<double>(t_solve_end - t_solve_start).count();
        
        if (rank == 0) {
            std::cout << "   Solved in " << t_solve << " seconds\n";
            std::cout << "   Converged: " << (converged ? "YES" : "NO") << "\n\n";
        }
        
        // Reconstruct full solution
        if (rank == 0) std::cout << "9. Reconstructing solution...\n";
        cstone::DeviceVector<float> u;
        eliminator.reconstructFullSolution(u_int, isBoundaryDOF, boundaryValues, u);
        
        // Compute solution statistics
        if (rank == 0) std::cout << "10. Solution statistics...\n";
        std::vector<float> h_u(u.size());
        thrust::copy(thrust::device_pointer_cast(u.data()),
                    thrust::device_pointer_cast(u.data() + u.size()),
                    h_u.begin());
        
        float u_min = *std::min_element(h_u.begin(), h_u.end());
        float u_max = *std::max_element(h_u.begin(), h_u.end());
        float u_sum = std::accumulate(h_u.begin(), h_u.end(), 0.0f);
        float u_mean = u_sum / h_u.size();
        
        float l2_norm = 0.0f;
        for (float val : h_u) l2_norm += val * val;
        l2_norm = std::sqrt(l2_norm);
        
        if (rank == 0) {
            std::cout << "   Min: " << u_min << "\n";
            std::cout << "   Max: " << u_max << "\n";
            std::cout << "   Mean: " << u_mean << "\n";
            std::cout << "   L2 norm: " << l2_norm << "\n\n";
            
            std::cout << "========================================\n";
            std::cout << "Timing Summary\n";
            std::cout << "========================================\n";
            std::cout << "Mesh loading:     " << t_mesh << " s\n";
            std::cout << "FE space:         " << t_fes << " s\n";
            std::cout << "Stiffness matrix: " << t_stiff << " s\n";
            std::cout << "RHS assembly:     " << t_rhs << " s\n";
            std::cout << "BC elimination:   " << t_elim << " s\n";
            std::cout << "Linear solve:     " << t_solve << " s\n";
            std::cout << "========================================\n\n";
            
            // Expected results for unit cube:
            // - u_min should be exactly 0 (boundary)
            // - u_max should be positive, at interior (roughly 0.01-0.02 for refined mesh)
            // - If u_max >> 0.1 or u_min < 0, something is wrong
            
            std::cout << "VALIDATION:\n";
            if (std::abs(u_min) < 1e-6) {
                std::cout << "✓ Boundary condition satisfied (u_min ≈ 0)\n";
            } else {
                std::cout << "✗ WARNING: Boundary condition violated (u_min = " << u_min << ")\n";
            }
            
            if (u_max > 0 && u_max < 1.0) {
                std::cout << "✓ Solution in reasonable range\n";
            } else {
                std::cout << "✗ WARNING: Solution out of expected range (u_max = " << u_max << ")\n";
            }
            
            if (converged) {
                std::cout << "✓ Solver converged\n";
            } else {
                std::cout << "✗ WARNING: Solver did not converge\n";
            }
            std::cout << "\n";
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
