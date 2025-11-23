#pragma once

#include "../fem/mars_sparse_matrix.hpp"
#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>
#include <mpi.h>
#include <iostream>
#include <vector>

namespace mars {
namespace fem {

// Full Hypre PCG solver with BoomerAMG preconditioner
// This replaces both preconditioner and solver for maximum Hypre integration
template<typename RealType, typename IndexType, typename AcceleratorTag>
class HyprePCGSolver {
public:
    using Matrix = SparseMatrix<IndexType, RealType, AcceleratorTag>;
    using Vector = typename mars::VectorSelector<RealType, AcceleratorTag>::type;

    HyprePCGSolver(MPI_Comm comm = MPI_COMM_WORLD, int maxIter = 1000, RealType tolerance = 1e-6)
        : comm_(comm), maxIter_(maxIter), tolerance_(tolerance),
          solver_(nullptr), precond_(nullptr), A_hypre_(nullptr),
          b_hypre_(nullptr), x_hypre_(nullptr), verbose_(true) {
    }

    ~HyprePCGSolver() {
        destroy();
    }

    // Solve Ax = b using Hypre PCG with BoomerAMG preconditioner
    bool solve(const Matrix& A, const Vector& b, Vector& x) {
        destroy();  // Clean up any existing setup

        HYPRE_Int m = static_cast<HYPRE_Int>(A.numRows());
        HYPRE_Int n = static_cast<HYPRE_Int>(A.numCols());

        if (m != n) {
            std::cerr << "HyprePCGSolver: Matrix must be square" << std::endl;
            return false;
        }

        // Convert MARS matrix to Hypre format
        setupHypreMatrix(A);

        // Create solution and RHS vectors
        HYPRE_Int ilower = 0, iupper = m - 1;  // Simplified for now

        if (verbose_) std::cout << "Creating vectors with range [" << ilower << ", " << iupper << "]" << std::endl;
        HYPRE_IJVectorCreate(comm_, ilower, iupper, &b_hypre_);
        HYPRE_IJVectorCreate(comm_, ilower, iupper, &x_hypre_);
        HYPRE_IJVectorSetObjectType(b_hypre_, HYPRE_PARCSR);
        HYPRE_IJVectorSetObjectType(x_hypre_, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(b_hypre_);
        HYPRE_IJVectorInitialize(x_hypre_);

        if (verbose_) std::cout << "Vectors created and initialized" << std::endl;

        // Set RHS values
        std::vector<HYPRE_Real> h_b(m);
        thrust::copy(thrust::device_pointer_cast(b.data()),
                    thrust::device_pointer_cast(b.data() + m),
                    h_b.begin());
        HYPRE_IJVectorSetValues(b_hypre_, m, NULL, h_b.data());

        // Set initial guess
        std::vector<HYPRE_Real> h_x(m);
        thrust::copy(thrust::device_pointer_cast(x.data()),
                    thrust::device_pointer_cast(x.data() + m),
                    h_x.begin());
        HYPRE_IJVectorSetValues(x_hypre_, m, NULL, h_x.data());

        HYPRE_IJVectorAssemble(b_hypre_);
        HYPRE_IJVectorAssemble(x_hypre_);

        HYPRE_IJVectorGetObject(b_hypre_, (void**)&par_b_);
        HYPRE_IJVectorGetObject(x_hypre_, (void**)&par_x_);
        
        if (!par_b_ || !par_x_) {
            std::cerr << "Failed to get Hypre ParVector objects" << std::endl;
            return false;
        }

        // Use MFEM ex1p style BoomerAMG settings (default/simple)
        if (verbose_) std::cout << "Creating BoomerAMG preconditioner (MFEM ex1p style)" << std::endl;
        HYPRE_BoomerAMGCreate(&precond_);
        HYPRE_BoomerAMGSetPrintLevel(precond_, 0);
        // Use default MFEM settings - no special tuning
        if (verbose_) std::cout << "BoomerAMG preconditioner created with default settings" << std::endl;

        // Create PCG solver
        if (verbose_) std::cout << "Creating PCG solver" << std::endl;
        HYPRE_ParCSRPCGCreate(comm_, &solver_);
        if (!solver_) {
            std::cerr << "Failed to create Hypre PCG solver" << std::endl;
            return false;
        }
        if (verbose_) std::cout << "PCG solver created" << std::endl;
        
        HYPRE_PCGSetMaxIter(solver_, maxIter_);
        HYPRE_PCGSetTol(solver_, tolerance_);
        HYPRE_PCGSetPrintLevel(solver_, verbose_ ? 2 : 0);  // Print convergence info
        
        // Set preconditioner (BoomerAMG)
        if (precond_) {
            HYPRE_PCGSetPrecond(solver_,
                               (HYPRE_Int (*)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector))HYPRE_BoomerAMGSolve,
                               (HYPRE_Int (*)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector))HYPRE_BoomerAMGSetup,
                               precond_);
            if (verbose_) std::cout << "PCG BoomerAMG preconditioner set" << std::endl;
        } else {
            if (verbose_) std::cout << "Using unpreconditioned CG" << std::endl;
        }

        // Setup and solve
        if (verbose_) std::cout << "Setting up PCG solver" << std::endl;
        HYPRE_ParCSRPCGSetup(solver_, parcsr_A_, par_b_, par_x_);
        if (verbose_) std::cout << "PCG setup completed, starting solve" << std::endl;
        HYPRE_ParCSRPCGSolve(solver_, parcsr_A_, par_b_, par_x_);
        if (verbose_) std::cout << "PCG solve completed" << std::endl;

        // Get solution
        HYPRE_IJVectorGetValues(x_hypre_, m, NULL, h_x.data());
        thrust::copy(h_x.begin(), h_x.end(),
                    thrust::device_pointer_cast(x.data()));

        // Check convergence
        int num_iterations;
        double final_res_norm;
        HYPRE_PCGGetNumIterations(solver_, &num_iterations);
        HYPRE_PCGGetFinalRelativeResidualNorm(solver_, &final_res_norm);

        if (verbose_) {
            std::cout << "Hypre PCG converged in " << num_iterations
                      << " iterations, final residual: " << final_res_norm << std::endl;
        }

        return (final_res_norm < tolerance_);
    }

    void setVerbose(bool verbose) { verbose_ = verbose; }

private:
    void setupHypreMatrix(const Matrix& A) {
        HYPRE_Int m = static_cast<HYPRE_Int>(A.numRows());
        HYPRE_Int ilower = 0, iupper = m - 1;

        if (verbose_) {
            std::cout << "Setting up Hypre matrix: " << m << " x " << A.numCols() 
                      << ", nnz = " << A.nnz() << std::endl;
        }

        // Get matrix data
        std::vector<IndexType> h_rowOffsets(m + 1);
        std::vector<IndexType> h_colIndices(A.nnz());
        std::vector<RealType> h_values(A.nnz());

        thrust::copy(thrust::device_pointer_cast(A.rowOffsetsPtr()),
                    thrust::device_pointer_cast(A.rowOffsetsPtr() + m + 1),
                    h_rowOffsets.begin());
        thrust::copy(thrust::device_pointer_cast(A.colIndicesPtr()),
                    thrust::device_pointer_cast(A.colIndicesPtr() + A.nnz()),
                    h_colIndices.begin());
        thrust::copy(thrust::device_pointer_cast(A.valuesPtr()),
                    thrust::device_pointer_cast(A.valuesPtr() + A.nnz()),
                    h_values.begin());

        // Debug: Check row offsets
        if (verbose_) {
            std::cout << "Row offsets: " << h_rowOffsets[0] << " ... " << h_rowOffsets[m] << std::endl;
            std::cout << "Last row offset should be nnz: " << h_rowOffsets[m] << " == " << A.nnz() << std::endl;
            
            // Check for empty rows or invalid structure
            bool hasEmptyRows = false;
            bool hasInvalidOffsets = false;
            for (HYPRE_Int i = 0; i < m; ++i) {
                if (h_rowOffsets[i] > h_rowOffsets[i+1]) {
                    hasInvalidOffsets = true;
                    std::cout << "Invalid offset at row " << i << ": " << h_rowOffsets[i] << " > " << h_rowOffsets[i+1] << std::endl;
                }
                if (h_rowOffsets[i] == h_rowOffsets[i+1]) {
                    hasEmptyRows = true;
                    std::cout << "Empty row at " << i << std::endl;
                }
            }
            if (hasEmptyRows) std::cout << "Warning: Matrix has empty rows" << std::endl;
            if (hasInvalidOffsets) std::cout << "Error: Matrix has invalid row offsets" << std::endl;
            
            // Check column indices range
            bool hasInvalidCols = false;
            for (size_t i = 0; i < h_colIndices.size(); ++i) {
                if (h_colIndices[i] >= static_cast<IndexType>(A.numCols())) {
                    hasInvalidCols = true;
                    std::cout << "Invalid column index " << h_colIndices[i] << " >= " << A.numCols() << " at position " << i << std::endl;
                    break;
                }
            }
            if (hasInvalidCols) std::cout << "Error: Matrix has invalid column indices" << std::endl;
        }

        // Create Hypre matrix
        HYPRE_IJMatrixCreate(comm_, ilower, iupper, ilower, iupper, &A_hypre_);
        HYPRE_IJMatrixSetObjectType(A_hypre_, HYPRE_PARCSR);
        HYPRE_IJMatrixInitialize(A_hypre_);

        // Set matrix values
        for (HYPRE_Int i = 0; i < m; ++i) {
            HYPRE_Int row = ilower + i;
            HYPRE_Int ncols = static_cast<HYPRE_Int>(h_rowOffsets[i+1] - h_rowOffsets[i]);
            std::vector<HYPRE_Int> cols(ncols);
            std::vector<HYPRE_Real> vals(ncols);

            for (HYPRE_Int j = 0; j < ncols; ++j) {
                cols[j] = static_cast<HYPRE_Int>(h_colIndices[h_rowOffsets[i] + j]);
                vals[j] = static_cast<HYPRE_Real>(h_values[h_rowOffsets[i] + j]);
            }

            HYPRE_IJMatrixSetValues(A_hypre_, 1, &ncols, &row, cols.data(), vals.data());
        }

        HYPRE_IJMatrixAssemble(A_hypre_);
        if (verbose_) std::cout << "Matrix assembled successfully" << std::endl;
        HYPRE_IJMatrixGetObject(A_hypre_, (void**)&parcsr_A_);
        
        if (!parcsr_A_) {
            std::cerr << "Failed to get Hypre ParCSR matrix object" << std::endl;
            return;
        }
        if (verbose_) std::cout << "Got ParCSR matrix object" << std::endl;
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