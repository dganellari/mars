// MARS Poisson Example - GPU-native finite element solver
// Solves: -Δu = f in Ω, u = 0 on ∂Ω, for a box-shaped Ω (the boundary is found as the faces of the
// mesh's global bounding box)
//
// Galerkin P1 on tetrahedra, like MFEM examples/ex1.cpp. On several ranks each rank assembles the rows
// of the nodes it owns (columns include its ghost nodes), and CG refreshes the ghost values through the
// domain's node halo: the same DOF numbering and solve as the Navier-Stokes solvers.
//
// Run:
//   mpirun -np 4 ./mars_ex1_poisson --mesh mesh_parts

#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_fem.hpp"
#include "backend/distributed/unstructured/solvers/mars_cg_solver.hpp"
#include <thrust/fill.h>
#include <thrust/for_each.h>
#include <thrust/functional.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/reduce.h>
#include <thrust/transform_reduce.h>
#include <mpi.h>
#include <algorithm>
#include <iostream>
#include <chrono>
#include <cmath>
#include <limits>
#include <string>

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

int main(int argc, char* argv[]) {
    // Initialize MPI
    MPI_Init(&argc, &argv);

    int rank, numRanks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    // Tet drivers must bind a GPU per rank before the domain touches the device.
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount > 0) cudaSetDevice(rank % deviceCount);

    // Parse command line arguments
    std::string meshPath = "mesh_parts";
    int order = 1;
    int maxIter = 500;
    float tolerance = 1e-6f;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.rfind("--mesh=", 0) == 0) {
            meshPath = arg.substr(7);
        } else if (arg == "--mesh" && i + 1 < argc) {
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
                          << "  --mesh <path>      Path to a box-shaped tet mesh (default: mesh_parts)\n"
                          << "  --order <n>        Polynomial order (default: 1)\n"
                          << "  --max-iter <n>     Max CG iterations (default: 500)\n"
                          << "  --tol <val>        CG tolerance (default: 1e-6)\n"
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

    bool converged = false;
    try {
        auto t_total_start = std::chrono::high_resolution_clock::now();

        // =====================================================
        // 1. Load mesh and create domain
        // =====================================================
        if (rank == 0) std::cout << "1. Loading mesh...\n";
        auto t_mesh_start = std::chrono::high_resolution_clock::now();

        using Domain = ElementDomain<TetTag, float, uint64_t, cstone::GpuTag>;
        Domain domain(meshPath, rank, numRanks);
        const auto& d_ownership = domain.getNodeOwnershipMap();
        size_t nodeCount        = domain.getNodeCount();

        auto t_mesh_end = std::chrono::high_resolution_clock::now();
        double t_mesh = std::chrono::duration<double>(t_mesh_end - t_mesh_start).count();

        if (rank == 0) {
            std::cout << "   Mesh loaded in " << t_mesh << " seconds\n"
                      << "   Local elements: " << domain.localElementCount() << "\n"
                      << "   Local nodes: " << nodeCount << "\n\n";
        }

        // =====================================================
        // 2. Create finite element space
        // =====================================================
        if (rank == 0) std::cout << "2. Creating finite element space...\n";
        auto t_fes_start = std::chrono::high_resolution_clock::now();

        TetFESpace<float, uint64_t> fes(domain, order);

        // Owned nodes first, then ghosts, from the domain's node ownership
        cstone::DeviceVector<int> d_nodeToDof(nodeCount);
        int numOwnedDofs = buildDofMappingGpu<uint64_t>(d_ownership.data(), d_nodeToDof.data(), nodeCount);
        long numGlobalDofs = numOwnedDofs;
        MPI_Allreduce(MPI_IN_PLACE, &numGlobalDofs, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

        auto t_fes_end = std::chrono::high_resolution_clock::now();
        double t_fes = std::chrono::duration<double>(t_fes_end - t_fes_start).count();

        if (rank == 0) {
            std::cout << "   FE space and DOF numbering in " << t_fes << " seconds\n"
                      << "   Global DOFs: " << numGlobalDofs << "\n\n";
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
        // 4. Assemble stiffness matrix: owned rows, owned + ghost columns
        // =====================================================
        if (rank == 0) std::cout << "4. Assembling stiffness matrix...\n";
        auto t_stiff_start = std::chrono::high_resolution_clock::now();

        SparseMatrix<int, float, cstone::GpuTag> A;
        stiffnessAssembler.assemble(domain, A, d_nodeToDof.data(), numOwnedDofs);

        auto t_stiff_end = std::chrono::high_resolution_clock::now();
        double t_stiff = std::chrono::duration<double>(t_stiff_end - t_stiff_start).count();

        if (rank == 0) {
            std::cout << "   Stiffness matrix assembled in " << t_stiff << " seconds\n"
                      << "   Rank 0 rows: " << A.numRows() << ", non-zeros: " << A.nnz() << "\n\n";
        }

        // =====================================================
        // 5. Assemble RHS vector on the owned DOFs
        // =====================================================
        if (rank == 0) std::cout << "5. Assembling RHS vector...\n";
        auto t_rhs_start = std::chrono::high_resolution_clock::now();

        TetMassAssembler<float, uint64_t> massAssembler;
        cstone::DeviceVector<float> b(numOwnedDofs);
        massAssembler.assembleRHS(domain, b, SourceTerm{}, d_nodeToDof.data(), numOwnedDofs);

        auto t_rhs_end = std::chrono::high_resolution_clock::now();
        double t_rhs = std::chrono::duration<double>(t_rhs_end - t_rhs_start).count();

        if (rank == 0) std::cout << "   RHS vector assembled in " << t_rhs << " seconds\n\n";

        // =====================================================
        // 6. Apply boundary conditions
        // =====================================================
        // Boundary rows become identity rows with zero right-hand side. CG starts from zero, so the
        // boundary entries of every iterate stay zero and the columns need no change; ghost copies of
        // boundary nodes stay zero through the halo exchange.
        if (rank == 0) std::cout << "6. Applying boundary conditions...\n";
        auto t_bc_start = std::chrono::high_resolution_clock::now();

        const auto& d_x = domain.getNodeX();
        const auto& d_y = domain.getNodeY();
        const auto& d_z = domain.getNodeZ();
        const float* coords[3] = {d_x.data(), d_y.data(), d_z.data()};
        float lo[3], hi[3];
        for (int d = 0; d < 3; ++d) {
            lo[d] = thrust::reduce(thrust::device, coords[d], coords[d] + nodeCount,
                                   std::numeric_limits<float>::max(), thrust::minimum<float>());
            hi[d] = thrust::reduce(thrust::device, coords[d], coords[d] + nodeCount,
                                   std::numeric_limits<float>::lowest(), thrust::maximum<float>());
        }
        MPI_Allreduce(MPI_IN_PLACE, lo, 3, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, hi, 3, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
        const float tol = 1e-6f * std::max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]});

        cstone::DeviceVector<uint8_t> d_isBoundaryDof(numOwnedDofs, 0);
        {
            const float *x = d_x.data(), *y = d_y.data(), *z = d_z.data();
            const int* nodeToDof = d_nodeToDof.data();
            uint8_t* isBoundary  = d_isBoundaryDof.data();
            float x0 = lo[0], x1 = hi[0], y0 = lo[1], y1 = hi[1], z0 = lo[2], z1 = hi[2];
            thrust::for_each(thrust::device, thrust::counting_iterator<size_t>(0),
                             thrust::counting_iterator<size_t>(nodeCount),
                             [=] __device__(size_t n) {
                                 int dof = nodeToDof[n];
                                 if (dof >= numOwnedDofs) return;
                                 bool onBoundary = fabsf(x[n] - x0) < tol || fabsf(x[n] - x1) < tol ||
                                                   fabsf(y[n] - y0) < tol || fabsf(y[n] - y1) < tol ||
                                                   fabsf(z[n] - z0) < tol || fabsf(z[n] - z1) < tol;
                                 if (onBoundary) isBoundary[dof] = 1;
                             });
        }

        {
            const int* rowPtr         = A.rowOffsetsPtr();
            const int* colInd         = A.colIndicesPtr();
            float* values             = A.valuesPtr();
            float* rhs                = b.data();
            const uint8_t* isBoundary = d_isBoundaryDof.data();
            thrust::for_each(thrust::device, thrust::counting_iterator<int>(0),
                             thrust::counting_iterator<int>(numOwnedDofs),
                             [=] __device__(int i) {
                                 if (!isBoundary[i]) return;
                                 for (int j = rowPtr[i]; j < rowPtr[i + 1]; ++j)
                                     values[j] = colInd[j] == i ? 1.0f : 0.0f;
                                 rhs[i] = 0.0f;
                             });
        }
        long numBoundary = thrust::reduce(thrust::device, d_isBoundaryDof.data(),
                                          d_isBoundaryDof.data() + numOwnedDofs, 0L);
        MPI_Allreduce(MPI_IN_PLACE, &numBoundary, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

        auto t_bc_end = std::chrono::high_resolution_clock::now();
        double t_bc = std::chrono::duration<double>(t_bc_end - t_bc_start).count();

        if (rank == 0) {
            std::cout << "   Boundary DOFs: " << numBoundary << " of " << numGlobalDofs << "\n"
                      << "   BCs applied in " << t_bc << " seconds\n\n";
        }

        // =====================================================
        // 7. Solve with CG
        // =====================================================
        if (rank == 0) std::cout << "7. Solving linear system (CG)...\n";
        auto t_solve_start = std::chrono::high_resolution_clock::now();

        ConjugateGradientSolver<float, int, cstone::GpuTag> solver(maxIter, tolerance);
        solver.setVerbose(false);
        solver.setOwnedSize(numOwnedDofs);
        if (numRanks > 1) {
            const int* dofMap = d_nodeToDof.data();
            solver.setHaloExchangeCallback(
                [&domain, dofMap](cstone::DeviceVector<float>& p) { domain.exchangeNodeHalo(p, dofMap); });
        }

        cstone::DeviceVector<float> u(nodeCount);
        thrust::fill(thrust::device, u.data(), u.data() + nodeCount, 0.0f);
        converged = solver.solve(A, b, u);

        auto t_solve_end = std::chrono::high_resolution_clock::now();
        double t_solve = std::chrono::duration<double>(t_solve_end - t_solve_start).count();

        if (rank == 0) {
            std::cout << "   System solved in " << t_solve << " seconds\n"
                      << "   Converged: " << (converged ? "Yes" : "No") << " (" << solver.getIterations()
                      << " iterations)\n\n";
        }

        // =====================================================
        // 8. Solution statistics over all owned DOFs of all ranks
        // =====================================================
        if (rank == 0) std::cout << "8. Computing solution statistics...\n";

        const float* uOwned = u.data();
        float uMin = thrust::reduce(thrust::device, uOwned, uOwned + numOwnedDofs,
                                    std::numeric_limits<float>::max(), thrust::minimum<float>());
        float uMax = thrust::reduce(thrust::device, uOwned, uOwned + numOwnedDofs,
                                    std::numeric_limits<float>::lowest(), thrust::maximum<float>());
        double sums[2] = {thrust::reduce(thrust::device, uOwned, uOwned + numOwnedDofs, 0.0),
                          thrust::transform_reduce(thrust::device, uOwned, uOwned + numOwnedDofs,
                                                   [] __device__(float v) -> double { return double(v) * v; }, 0.0,
                                                   thrust::plus<double>())};
        MPI_Allreduce(MPI_IN_PLACE, &uMin, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &uMax, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, sums, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        if (rank == 0) {
            std::cout << "   Solution statistics:\n"
                      << "     Min: " << uMin << "\n"
                      << "     Max: " << uMax << "\n"
                      << "     Mean: " << sums[0] / numGlobalDofs << "\n"
                      << "     L2 norm: " << std::sqrt(sums[1]) << "\n\n";
        }

        // =====================================================
        // 9. Timing summary
        // =====================================================
        auto t_total_end = std::chrono::high_resolution_clock::now();
        double t_total = std::chrono::duration<double>(t_total_end - t_total_start).count();

        if (rank == 0) {
            std::cout << "========================================\n"
                      << "   Timing Summary\n"
                      << "========================================\n"
                      << "Mesh loading:     " << t_mesh << " s\n"
                      << "FE space + DOFs:  " << t_fes << " s\n"
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
    return converged ? 0 : 1;
}
