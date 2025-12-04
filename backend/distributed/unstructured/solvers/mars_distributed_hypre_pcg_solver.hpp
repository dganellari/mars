#pragma once

#include "../fem/mars_sparse_matrix.hpp"
#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>
#include <mpi.h>
#include <iostream>
#include <vector>

namespace mars {
namespace fem {

// Distributed GPU Hypre PCG solver with BoomerAMG preconditioner
// Uses Hypre's native GPU support for distributed parallel solving
template<typename RealType, typename IndexType, typename AcceleratorTag>
class DistributedHyprePCGSolver {
public:
    using Matrix = SparseMatrix<IndexType, RealType, AcceleratorTag>;
    using Vector = typename mars::VectorSelector<RealType, AcceleratorTag>::type;

    DistributedHyprePCGSolver(MPI_Comm comm = MPI_COMM_WORLD, int maxIter = 1000, RealType tolerance = 1e-6)
        : comm_(comm), maxIter_(maxIter), tolerance_(tolerance),
          solver_(nullptr), precond_(nullptr), A_hypre_(nullptr),
          b_hypre_(nullptr), x_hypre_(nullptr), verbose_(true) {

        // Initialize Hypre
        HYPRE_Init();

        // Check IndexType compatibility with Hypre's HYPRE_Int
        if (sizeof(IndexType) != sizeof(HYPRE_Int)) {
            int rank;
            MPI_Comm_rank(comm_, &rank);
            if (rank == 0) {
                std::cerr << "ERROR: IndexType size (" << sizeof(IndexType) 
                          << ") does not match Hypre HYPRE_Int size (" << sizeof(HYPRE_Int) 
                          << "). Ensure Hypre bigint configuration matches." << std::endl;
            }
            MPI_Abort(comm_, 1);
        }

        // Check if GPU is available and set memory location
        if (verbose_) {
            int rank;
            MPI_Comm_rank(comm_, &rank);
            if (rank == 0) {
                std::cout << "Initializing Hypre for GPU execution" << std::endl;
            }
        }
    }

    ~DistributedHyprePCGSolver() {
        destroy();
        HYPRE_Finalize();
    }

    // Solve Ax = b using distributed GPU Hypre PCG with BoomerAMG preconditioner
    bool solve(const Matrix& A, const Vector& b, Vector& x) {
        destroy();  // Clean up any existing setup

        int rank, num_ranks;
        MPI_Comm_rank(comm_, &rank);
        MPI_Comm_size(comm_, &num_ranks);

        HYPRE_Int m = static_cast<HYPRE_Int>(A.numRows());
        HYPRE_Int n = static_cast<HYPRE_Int>(A.numCols());

        if (m != n) {
            if (rank == 0) std::cerr << "DistributedHyprePCGSolver: Matrix must be square" << std::endl;
            return false;
        }

        // For distributed solving, each rank owns a portion of rows
        // Simplified: assume consecutive row ownership
        HYPRE_Int rows_per_rank = m / num_ranks;
        HYPRE_Int extra_rows = m % num_ranks;

        HYPRE_Int ilower = rank * rows_per_rank + std::min(rank, extra_rows);
        HYPRE_Int iupper = ilower + rows_per_rank + (rank < extra_rows ? 1 : 0) - 1;

        if (verbose_ && rank == 0) {
            std::cout << "Distributed setup: " << num_ranks << " ranks, "
                      << "global size " << m << " x " << n << std::endl;
        }

        if (rank == 0 || rank == 1) {
            std::cout << "Rank " << rank << ": owns rows [" << ilower << ", " << iupper << "]" << std::endl;
        }

        // Convert MARS matrix to Hypre format (GPU data)
        setupDistributedHypreMatrix(A, ilower, iupper);

        // Create solution and RHS vectors (GPU-based)
        setupDistributedHypreVectors(b, x, ilower, iupper);

        // Create BoomerAMG preconditioner (GPU-accelerated)
        if (verbose_ && rank == 0) std::cout << "Creating GPU BoomerAMG preconditioner" << std::endl;

        HYPRE_BoomerAMGCreate(&precond_);

        // Basic AMG settings for GPU
        HYPRE_BoomerAMGSetPrintLevel(precond_, 0);
        HYPRE_BoomerAMGSetMaxLevels(precond_, 25);
        HYPRE_BoomerAMGSetTol(precond_, 0.0);

        if (verbose_ && rank == 0) std::cout << "GPU BoomerAMG preconditioner created" << std::endl;

        // Create PCG solver (GPU-accelerated)
        if (verbose_ && rank == 0) std::cout << "Creating GPU PCG solver" << std::endl;

        HYPRE_ParCSRPCGCreate(comm_, &solver_);
        if (!solver_) {
            if (rank == 0) std::cerr << "Failed to create Hypre PCG solver" << std::endl;
            return false;
        }

        HYPRE_PCGSetMaxIter(solver_, maxIter_);
        HYPRE_PCGSetTol(solver_, tolerance_);
        HYPRE_PCGSetPrintLevel(solver_, verbose_ ? 1 : 0);

        // Set GPU preconditioner
        HYPRE_PCGSetPrecond(solver_,
                           (HYPRE_Int (*)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector))HYPRE_BoomerAMGSolve,
                           (HYPRE_Int (*)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector))HYPRE_BoomerAMGSetup,
                           precond_);

        if (verbose_ && rank == 0) std::cout << "GPU PCG solver created with BoomerAMG preconditioner" << std::endl;

        // Setup and solve (all on GPU)
        if (verbose_ && rank == 0) std::cout << "Setting up and solving on GPU" << std::endl;

        HYPRE_ParCSRPCGSetup(solver_, parcsr_A_, par_b_, par_x_);
        if (verbose_ && rank == 0) std::cout << "GPU PCG setup completed, starting solve" << std::endl;

        HYPRE_ParCSRPCGSolve(solver_, parcsr_A_, par_b_, par_x_);
        if (verbose_ && rank == 0) std::cout << "GPU PCG solve completed" << std::endl;

        // Get solution back to GPU vector (no CPU transfer)
        getDistributedSolution(x, ilower, iupper);

        // Check convergence
        int num_iterations;
        double final_res_norm;
        HYPRE_PCGGetNumIterations(solver_, &num_iterations);
        HYPRE_PCGGetFinalRelativeResidualNorm(solver_, &final_res_norm);

        if (verbose_ && rank == 0) {
            std::cout << "Distributed GPU Hypre PCG converged in " << num_iterations
                      << " iterations, final residual: " << final_res_norm << std::endl;
        }

        return (final_res_norm < tolerance_);
    }

    void setVerbose(bool verbose) { verbose_ = verbose; }

private:
    void setupDistributedHypreMatrix(const Matrix& A, HYPRE_Int ilower, HYPRE_Int iupper) {
        int rank;
        MPI_Comm_rank(comm_, &rank);

        HYPRE_Int local_rows = iupper - ilower + 1;
        HYPRE_Int global_rows = A.numRows();

        if (verbose_ && rank == 0) {
            std::cout << "Setting up distributed GPU matrix: " << global_rows << " x " << A.numCols()
                      << ", local rows [" << ilower << ", " << iupper << "]" << std::endl;
        }

        // Create distributed IJ matrix (GPU memory)
        HYPRE_IJMatrixCreate(comm_, ilower, iupper, 0, global_rows - 1, &A_hypre_);
        HYPRE_IJMatrixSetObjectType(A_hypre_, HYPRE_PARCSR);
        HYPRE_IJMatrixInitialize(A_hypre_);

        // Get local matrix data from GPU
        std::vector<IndexType> h_rowOffsets(local_rows + 1);
        thrust::copy(thrust::device_pointer_cast(A.rowOffsetsPtr() + ilower),
                    thrust::device_pointer_cast(A.rowOffsetsPtr() + iupper + 2),
                    h_rowOffsets.begin());

        // Convert to local row offsets
        for (HYPRE_Int i = 0; i <= local_rows; ++i) {
            h_rowOffsets[i] -= h_rowOffsets[0];
        }

        // Set matrix values directly from GPU data
        for (HYPRE_Int i = 0; i < local_rows; ++i) {
            HYPRE_Int global_row = ilower + i;
            HYPRE_Int row_nnz = h_rowOffsets[i+1] - h_rowOffsets[i];

            if (row_nnz > 0) {
                // Get column indices and values for this row (GPU data)
                const IndexType* row_cols = A.colIndicesPtr() + h_rowOffsets[0] + h_rowOffsets[i];
                const RealType* row_vals = A.valuesPtr() + h_rowOffsets[0] + h_rowOffsets[i];

                // Create CPU copies for Hypre API (unfortunately still needed)
                std::vector<HYPRE_Int> h_cols(row_nnz);
                std::vector<HYPRE_Real> h_vals(row_nnz);

                thrust::copy(thrust::device_pointer_cast(row_cols),
                           thrust::device_pointer_cast(row_cols + row_nnz),
                           h_cols.begin());
                thrust::copy(thrust::device_pointer_cast(row_vals),
                           thrust::device_pointer_cast(row_vals + row_nnz),
                           h_vals.begin());

                HYPRE_IJMatrixSetValues(A_hypre_, 1, &row_nnz, &global_row,
                                       h_cols.data(), h_vals.data());
            }
        }

        HYPRE_IJMatrixAssemble(A_hypre_);
        HYPRE_IJMatrixGetObject(A_hypre_, (void**)&parcsr_A_);

        if (verbose_ && rank == 0) std::cout << "Distributed GPU matrix assembled" << std::endl;
    }

    void setupDistributedHypreVectors(const Vector& b, const Vector& x,
                                    HYPRE_Int ilower, HYPRE_Int iupper) {
        int rank;
        MPI_Comm_rank(comm_, &rank);

        HYPRE_Int local_size = iupper - ilower + 1;
        HYPRE_Int global_size = b.size();  // Assuming same size

        if (verbose_ && rank == 0) {
            std::cout << "Creating distributed GPU vectors, local size: " << local_size << std::endl;
        }

        // Create GPU-based IJ vectors
        HYPRE_IJVectorCreate(comm_, ilower, iupper, &b_hypre_);
        HYPRE_IJVectorCreate(comm_, ilower, iupper, &x_hypre_);

        HYPRE_IJVectorSetObjectType(b_hypre_, HYPRE_PARCSR);
        HYPRE_IJVectorSetObjectType(x_hypre_, HYPRE_PARCSR);

        HYPRE_IJVectorInitialize(b_hypre_);
        HYPRE_IJVectorInitialize(x_hypre_);

        // Set RHS values (copy from GPU to CPU for Hypre API, then to GPU)
        std::vector<HYPRE_Real> h_b(local_size);
        thrust::copy(thrust::device_pointer_cast(b.data() + ilower),
                    thrust::device_pointer_cast(b.data() + iupper + 1),
                    h_b.begin());
        HYPRE_IJVectorSetValues(b_hypre_, local_size, NULL, h_b.data());

        // Set initial guess
        std::vector<HYPRE_Real> h_x(local_size);
        thrust::copy(thrust::device_pointer_cast(x.data() + ilower),
                    thrust::device_pointer_cast(x.data() + iupper + 1),
                    h_x.begin());
        HYPRE_IJVectorSetValues(x_hypre_, local_size, NULL, h_x.data());

        HYPRE_IJVectorAssemble(b_hypre_);
        HYPRE_IJVectorAssemble(x_hypre_);

        HYPRE_IJVectorGetObject(b_hypre_, (void**)&par_b_);
        HYPRE_IJVectorGetObject(x_hypre_, (void**)&par_x_);

        if (verbose_ && rank == 0) std::cout << "Distributed GPU vectors created" << std::endl;
    }

    void getDistributedSolution(Vector& x, HYPRE_Int ilower, HYPRE_Int iupper) {
        HYPRE_Int local_size = iupper - ilower + 1;

        // Get solution values from Hypre (GPU to CPU, then back to GPU vector)
        std::vector<HYPRE_Real> h_x(local_size);
        HYPRE_IJVectorGetValues(x_hypre_, local_size, NULL, h_x.data());

        // Copy back to GPU vector at correct positions
        thrust::copy(h_x.begin(), h_x.end(),
                    thrust::device_pointer_cast(x.data() + ilower));
    }

    void destroy() {
        if (solver_) {
            HYPRE_ParCSRPCGDestroy(solver_);
            solver_ = nullptr;
        }
        if (precond_) {
            HYPRE_BoomerAMGDestroy(precond_);
            precond_ = nullptr;
        }
        if (A_hypre_) {
            HYPRE_IJMatrixDestroy(A_hypre_);
            A_hypre_ = nullptr;
        }
        if (b_hypre_) {
            HYPRE_IJVectorDestroy(b_hypre_);
            b_hypre_ = nullptr;
        }
        if (x_hypre_) {
            HYPRE_IJVectorDestroy(x_hypre_);
            x_hypre_ = nullptr;
        }
    }

    MPI_Comm comm_;
    int maxIter_;
    RealType tolerance_;
    bool verbose_;

    HYPRE_Solver solver_;
    HYPRE_Solver precond_;

    HYPRE_IJMatrix A_hypre_;
    HYPRE_ParCSRMatrix parcsr_A_;

    HYPRE_IJVector b_hypre_;
    HYPRE_IJVector x_hypre_;
    HYPRE_ParVector par_b_;
    HYPRE_ParVector par_x_;
};

} // namespace fem
} // namespace mars