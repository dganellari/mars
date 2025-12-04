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

    enum PrecondType { BOOMERAMG, JACOBI };

    HyprePCGSolver(MPI_Comm comm = MPI_COMM_WORLD, int maxIter = 1000, RealType tolerance = 1e-6, PrecondType precondType = BOOMERAMG)
        : comm_(comm), maxIter_(maxIter), tolerance_(tolerance),
          solver_(nullptr), precond_(nullptr), A_hypre_(nullptr),
          b_hypre_(nullptr), x_hypre_(nullptr), verbose_(true), precondType_(precondType) {
    }

    ~HyprePCGSolver() {
        destroy();
    }

    // Solve Ax = b using Hypre PCG with BoomerAMG preconditioner (simple version)
    bool solve(const Matrix& A, const Vector& b, Vector& x, 
               IndexType globalDofStart = 0, IndexType globalDofEnd = 0) {
        // For simple case, assume square matrix (column range = row range)
        return solve<IndexType>(A, b, x, globalDofStart, globalDofEnd, 
                               globalDofStart, globalDofEnd, {});
    }

    // Solve Ax = b using Hypre PCG with BoomerAMG preconditioner
    // globalDofStart and globalDofEnd define this rank's owned global DOF range [start, end)
    // globalColStart and globalColEnd define the global column range [start, end)
    // localToGlobalDof maps local column indices to global DOF indices
    template<typename KeyType>
    bool solve(const Matrix& A, const Vector& b, Vector& x, 
               IndexType globalDofStart, IndexType globalDofEnd,
               IndexType globalColStart, IndexType globalColEnd,
               const std::vector<KeyType>& localToGlobalDof) {
        int rank;
        MPI_Comm_rank(comm_, &rank);
        std::cout << "Rank " << rank << ": Entering Hypre solve with globalDofRange [" 
                  << globalDofStart << ", " << globalDofEnd << ")" << std::endl;
        
        destroy();  // Clean up any existing setup

        HYPRE_Int m = static_cast<HYPRE_Int>(A.numRows());
        HYPRE_Int n = static_cast<HYPRE_Int>(A.numCols());

        // Note: For distributed FEM, matrix may be rectangular (owned rows × owned+ghost columns)
        // Hypre will handle the parallel structure using global column indices

        // If global DOF range not provided, assume single rank or local numbering
        if (globalDofEnd == 0) {
            globalDofEnd = globalDofStart + m;
        }
        
        globalDofStart_ = globalDofStart;
        globalDofEnd_ = globalDofEnd;
        
        // Convert the provided local-to-global DOF mapping to HYPRE format
        std::vector<HYPRE_BigInt> hypreLocalToGlobalDof(localToGlobalDof.size());
        for (size_t i = 0; i < localToGlobalDof.size(); ++i) {
            hypreLocalToGlobalDof[i] = static_cast<HYPRE_BigInt>(localToGlobalDof[i]);
        }
        
        // Convert MARS matrix to Hypre format
        setupHypreMatrix(A, hypreLocalToGlobalDof, globalColStart, globalColEnd);
        
        // Create solution and RHS vectors
        HYPRE_Int ilower = static_cast<HYPRE_Int>(globalDofStart_);
        HYPRE_Int iupper = static_cast<HYPRE_Int>(globalDofEnd_ - 1);

        HYPRE_IJVectorCreate(comm_, ilower, iupper, &b_hypre_);
        HYPRE_IJVectorCreate(comm_, ilower, iupper, &x_hypre_);
        HYPRE_IJVectorSetObjectType(b_hypre_, HYPRE_PARCSR);
        HYPRE_IJVectorSetObjectType(x_hypre_, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(b_hypre_);
        HYPRE_IJVectorInitialize(x_hypre_);

        // Set RHS values - need to provide global indices for each value
        std::vector<HYPRE_Real> h_b(m);
        std::vector<HYPRE_Int> indices(m);
        
        thrust::copy(thrust::device_pointer_cast(b.data()),
                    thrust::device_pointer_cast(b.data() + m),
                    h_b.begin());
        
        // Validate RHS for NaNs/Infs
        bool hasNaN = false;
        double minVal = 1e100, maxVal = -1e100, sumVal = 0.0;
        int nanCount = 0;
        for (HYPRE_Int i = 0; i < m; ++i) {
            if (std::isnan(h_b[i]) || std::isinf(h_b[i])) {
                hasNaN = true;
                nanCount++;
            }
            if (std::isfinite(h_b[i])) {
                minVal = std::min(minVal, static_cast<double>(h_b[i]));
                maxVal = std::max(maxVal, static_cast<double>(h_b[i]));
                sumVal += h_b[i];
            }
        }
        
        if (verbose_ || hasNaN) {
            std::cout << "Rank " << rank << ": Before Hypre - RHS sum=" << sumVal 
                      << ", range=[" << minVal << ", " << maxVal << "], NaNs=" << nanCount << "/" << m << std::endl;
            if (hasNaN) {
                std::cout << "Rank " << rank << " ERROR: RHS contains NaN/Inf!" << std::endl;
            }
        }
        
        // Create global index array for this rank's DOFs
        for (HYPRE_Int i = 0; i < m; ++i) {
            indices[i] = ilower + i;
        }
        
        HYPRE_IJVectorSetValues(b_hypre_, m, indices.data(), h_b.data());

        // Set initial guess with global indices
        std::vector<HYPRE_Real> h_x(m);
        thrust::copy(thrust::device_pointer_cast(x.data()),
                    thrust::device_pointer_cast(x.data() + m),
                    h_x.begin());
        HYPRE_IJVectorSetValues(x_hypre_, m, indices.data(), h_x.data());

        HYPRE_IJVectorAssemble(b_hypre_);
        HYPRE_IJVectorAssemble(x_hypre_);

        std::cout << "Rank " << rank << ": Vectors assembled, getting ParVector objects..." << std::endl;
        HYPRE_IJVectorGetObject(b_hypre_, (void**)&par_b_);
        HYPRE_IJVectorGetObject(x_hypre_, (void**)&par_x_);
        
        if (!par_b_ || !par_x_) {
            std::cerr << "Rank " << rank << ": Failed to get Hypre ParVector objects" << std::endl;
            return false;
        }
        std::cout << "Rank " << rank << ": ParVector objects obtained successfully" << std::endl;

        if (verbose_ && rank == 0) std::cout << "Creating preconditioner..." << std::endl;
        
        // Create preconditioner based on type
        if (precondType_ == BOOMERAMG) {
            if (verbose_ && rank == 0) std::cout << "Using BoomerAMG preconditioner" << std::endl;
            HYPRE_BoomerAMGCreate(&precond_);
            HYPRE_BoomerAMGSetPrintLevel(precond_, 0);
        } else if (precondType_ == JACOBI) {
            if (verbose_ && rank == 0) std::cout << "Using Jacobi preconditioner" << std::endl;
            // Create Jacobi (diagonal scaling) preconditioner
            precond_ = (HYPRE_Solver) parcsr_A_;  // Hypre uses matrix for diagonal preconditioner
        }

        // Create PCG solver
        if (verbose_ && rank == 0) std::cout << "Creating PCG solver..." << std::endl;
        HYPRE_ParCSRPCGCreate(comm_, &solver_);
        if (!solver_) {
            std::cerr << "Failed to create Hypre PCG solver" << std::endl;
            return false;
        }
        
        HYPRE_PCGSetMaxIter(solver_, maxIter_);
        HYPRE_PCGSetTol(solver_, tolerance_);
        HYPRE_PCGSetPrintLevel(solver_, verbose_ ? 2 : 0);
        
        // Set preconditioner
        if (precondType_ == BOOMERAMG && precond_) {
            if (verbose_ && rank == 0) std::cout << "Setting BoomerAMG preconditioner..." << std::endl;
            HYPRE_PCGSetPrecond(solver_,
                               (HYPRE_Int (*)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector))HYPRE_BoomerAMGSolve,
                               (HYPRE_Int (*)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector))HYPRE_BoomerAMGSetup,
                               precond_);
        } else if (precondType_ == JACOBI) {
            if (verbose_ && rank == 0) std::cout << "Setting Jacobi (diagonal) preconditioner..." << std::endl;
            HYPRE_PCGSetPrecond(solver_,
                               (HYPRE_Int (*)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector))HYPRE_ParCSRDiagScale,
                               (HYPRE_Int (*)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector))HYPRE_ParCSRDiagScaleSetup,
                               (HYPRE_Solver) parcsr_A_);
        }

        // Setup and solve
        if (verbose_ && rank == 0) std::cout << "Setting up PCG solver..." << std::endl;
        std::cout << "Rank " << rank << ": Waiting at barrier before PCG setup..." << std::endl;
        MPI_Barrier(comm_);
        std::cout << "Rank " << rank << ": Starting PCG setup..." << std::endl;
        HYPRE_ParCSRPCGSetup(solver_, parcsr_A_, par_b_, par_x_);
        std::cout << "Rank " << rank << ": PCG setup complete" << std::endl;
        if (verbose_ && rank == 0) std::cout << "PCG setup complete, starting solve..." << std::endl;
        MPI_Barrier(comm_);
        HYPRE_ParCSRPCGSolve(solver_, parcsr_A_, par_b_, par_x_);
        if (verbose_ && rank == 0) std::cout << "PCG solve complete." << std::endl;

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
    void setupHypreMatrix(const Matrix& A, const std::vector<HYPRE_BigInt>& localToGlobalDof,
                          IndexType globalColStart, IndexType globalColEnd) {
        HYPRE_Int m = static_cast<HYPRE_Int>(A.numRows());
        HYPRE_Int ilower = static_cast<HYPRE_Int>(globalDofStart_);
        HYPRE_Int iupper = static_cast<HYPRE_Int>(globalDofEnd_ - 1);
        
        // CRITICAL: For BoomerAMG, column range must match the global interior DOF space [0, N)
        // NOT the local column range which includes ghosts
        // Hypre ParCSR will automatically handle off-processor columns
        HYPRE_Int jlower = 0;
        HYPRE_Int jupper = static_cast<HYPRE_Int>(globalColEnd - 1);  // globalColEnd should be total interior DOFs

        // Matrix setup for Hypre

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

        // Validate matrix structure
        bool hasNaN = false, hasInf = false;
        for (size_t i = 0; i < h_values.size(); ++i) {
            if (std::isnan(h_values[i]) || std::isinf(h_values[i])) {
                if (verbose_) std::cout << "ERROR: Matrix contains NaN/Inf values!" << std::endl;
                return;
            }
        }

        // Create Hypre matrix
        // For distributed FEM, we want the vector x to be distributed exactly like b.
        // So the column partitioning (jlower, jupper) must match the row partitioning (ilower, iupper).
        // Hypre handles off-processor columns (ghosts) automatically in the ParCSR format.
        
        HYPRE_IJMatrixCreate(comm_, ilower, iupper, ilower, iupper, &A_hypre_);
        HYPRE_IJMatrixSetObjectType(A_hypre_, HYPRE_PARCSR);
        HYPRE_IJMatrixInitialize(A_hypre_);

        // Set matrix values row by row
        // Note: Column indices must be in GLOBAL numbering for ParCSR format
        int emptyRows = 0, rowsMissingDiagonal = 0;
        size_t totalFiltered = 0, totalOriginal = 0;
        
        HYPRE_BigInt minGlobalCol = HYPRE_BigInt(1e18), maxGlobalCol = -1;
        bool foundInvalidIndex = false;
        
        for (HYPRE_Int i = 0; i < m; ++i) {
            HYPRE_Int row = ilower + i;  // Global row index
            HYPRE_Int ncols = static_cast<HYPRE_Int>(h_rowOffsets[i+1] - h_rowOffsets[i]);
            totalOriginal += ncols;
            std::vector<HYPRE_Int> cols(ncols);
            std::vector<HYPRE_Real> vals(ncols);

            // Filter out columns with invalid global indices
            std::vector<HYPRE_Int> validCols;
            std::vector<HYPRE_Real> validVals;
            validCols.reserve(ncols);
            validVals.reserve(ncols);
            
            bool hasDiagonal = false;
            HYPRE_Int globalRow = ilower + i;
            
            for (HYPRE_Int j = 0; j < ncols; ++j) {
                IndexType localCol = h_colIndices[h_rowOffsets[i] + j];
                HYPRE_BigInt globalCol = -1;
                
                // Convert local column index to global using provided mapping
                if (!localToGlobalDof.empty()) {
                    if (localCol >= localToGlobalDof.size()) {
                        if (i == 0 && !foundInvalidIndex) {
                            std::cerr << "Rank " << rank << ": WARNING - local column " << localCol 
                                      << " >= mapping size " << localToGlobalDof.size() << std::endl;
                            foundInvalidIndex = true;
                        }
                        continue;  // Skip out of range columns
                    }
                    globalCol = localToGlobalDof[localCol];
                } else {
                    // Fallback: assume local indices are in owned range
                    globalCol = ilower + static_cast<HYPRE_Int>(localCol);
                }
                
                // Skip columns with invalid global indices
                if (globalCol < 0) {
                    continue;
                }
                
                minGlobalCol = std::min(minGlobalCol, globalCol);
                maxGlobalCol = std::max(maxGlobalCol, globalCol);
                
                if (globalCol == globalRow) {
                    hasDiagonal = true;
                }
                
                validCols.push_back(static_cast<HYPRE_Int>(globalCol));
                validVals.push_back(static_cast<HYPRE_Real>(h_values[h_rowOffsets[i] + j]));
            }

            totalFiltered += validCols.size();
            if (!hasDiagonal) rowsMissingDiagonal++;
            
            // Only set row if it has valid columns
            if (!validCols.empty()) {
                HYPRE_Int validNcols = static_cast<HYPRE_Int>(validCols.size());
                HYPRE_IJMatrixSetValues(A_hypre_, 1, &validNcols, &row, validCols.data(), validVals.data());
            } else {
                // Empty row - this will cause NaN in Hypre!
                emptyRows++;
                int mpiRank;
                MPI_Comm_rank(comm_, &mpiRank);
                std::cout << "Rank " << mpiRank << ": ERROR - Row " << i << " (global " << globalRow 
                          << ") has no valid columns after filtering! Original ncols=" << ncols << std::endl;
            }
        }

        // Print filtering statistics
        int mpiRank;
        MPI_Comm_rank(comm_, &mpiRank);
        std::cout << "Rank " << mpiRank << ": Matrix filtering: " << totalOriginal 
                  << " original entries → " << totalFiltered << " valid entries ("
                  << (totalOriginal - totalFiltered) << " filtered)" << std::endl;
        std::cout << "Rank " << mpiRank << ": Global column range: [" << minGlobalCol 
                  << ", " << maxGlobalCol << "], global row range: [" << ilower 
                  << ", " << iupper << "]" << std::endl;
        if (emptyRows > 0) {
            std::cout << "Rank " << mpiRank << ": WARNING - " << emptyRows << " rows became empty after filtering!" << std::endl;
        }
        if (rowsMissingDiagonal > 0) {
            std::cout << "Rank " << mpiRank << ": WARNING - " << rowsMissingDiagonal 
                      << " rows missing diagonal after filtering!" << std::endl;
        }

        HYPRE_IJMatrixAssemble(A_hypre_);
        std::cout << "Rank " << mpiRank << ": Matrix assembled, getting ParCSR object..." << std::endl;
        HYPRE_IJMatrixGetObject(A_hypre_, (void**)&parcsr_A_);
        
        if (!parcsr_A_) {
            std::cerr << "Rank " << mpiRank << ": Failed to get Hypre ParCSR matrix object" << std::endl;
            return;
        }
        std::cout << "Rank " << mpiRank << ": ParCSR matrix object obtained successfully" << std::endl;
    }

    void destroy() {
        if (solver_) {
            HYPRE_ParCSRPCGDestroy(solver_);
            solver_ = nullptr;
        }
        if (precond_ && precondType_ == BOOMERAMG) {
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
    PrecondType precondType_;

    HYPRE_Solver solver_;
    HYPRE_Solver precond_;

    HYPRE_IJMatrix A_hypre_;
    HYPRE_ParCSRMatrix parcsr_A_;

    HYPRE_IJVector b_hypre_;
    HYPRE_IJVector x_hypre_;
    HYPRE_ParVector par_b_;
    HYPRE_ParVector par_x_;
    
    IndexType globalDofStart_;  // Global DOF range start (inclusive)
    IndexType globalDofEnd_;    // Global DOF range end (exclusive)
};

} // namespace fem
} // namespace mars