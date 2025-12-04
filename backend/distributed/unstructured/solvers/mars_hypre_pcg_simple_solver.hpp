#pragma once

#include "../fem/mars_sparse_matrix.hpp"
#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>
#include <mpi.h>
#include <iostream>
#include <vector>

namespace mars {
namespace fem {

// Hypre PCG solver with Jacobi (diagonal) preconditioner
// Works well with full system matrix (zero rows/columns, diagonal=1)
// More robust and better conditioned than penalty method
template<typename RealType, typename IndexType, typename AcceleratorTag>
class HyprePCGJacobiSolver {
public:
    using Matrix = SparseMatrix<IndexType, RealType, AcceleratorTag>;
    using Vector = typename mars::VectorSelector<RealType, AcceleratorTag>::type;

    HyprePCGJacobiSolver(MPI_Comm comm = MPI_COMM_WORLD, int maxIter = 1000, RealType tolerance = 1e-6)
        : comm_(comm), maxIter_(maxIter), tolerance_(tolerance),
          solver_(nullptr), A_hypre_(nullptr),
          b_hypre_(nullptr), x_hypre_(nullptr), verbose_(true) {
        MPI_Comm_rank(comm_, &rank_);
    }

    ~HyprePCGJacobiSolver() {
        destroy();
    }

    // Solve Ax = b using Hypre PCG with Jacobi preconditioner
    // Works with full system matrix (boundary rows/cols zeroed, diagonal=1)
    template<typename KeyType>
    bool solve(const Matrix& A, const Vector& b, Vector& x, 
               IndexType globalDofStart, IndexType globalDofEnd,
               const std::vector<KeyType>& nodeToLocalDof,
               const std::vector<KeyType>& nodeToGlobalDof) {
        
        destroy();  // Clean up any existing setup

        HYPRE_Int m = static_cast<HYPRE_Int>(A.numRows());
        
        if (globalDofEnd == 0) {
            globalDofEnd = globalDofStart + m;
        }
        
        // Build local DOF to global DOF mapping
        std::vector<HYPRE_BigInt> localToGlobalDof;
        if (!nodeToLocalDof.empty() && !nodeToGlobalDof.empty()) {
            size_t maxLocalDof = 0;
            for (auto localDof : nodeToLocalDof) {
                if (localDof != static_cast<decltype(localDof)>(-1)) {
                    maxLocalDof = std::max(maxLocalDof, static_cast<size_t>(localDof));
                }
            }
            localToGlobalDof.resize(maxLocalDof + 1, -1);
            
            for (size_t nodeIdx = 0; nodeIdx < nodeToLocalDof.size() && nodeIdx < nodeToGlobalDof.size(); ++nodeIdx) {
                auto localDof = nodeToLocalDof[nodeIdx];
                auto globalDof = nodeToGlobalDof[nodeIdx];
                if (localDof != static_cast<decltype(localDof)>(-1) && 
                    globalDof != static_cast<decltype(globalDof)>(-1) &&
                    localDof <= maxLocalDof) {
                    localToGlobalDof[localDof] = static_cast<HYPRE_BigInt>(globalDof);
                }
            }
        }
        
        // Setup Hypre matrix
        setupHypreMatrix(A, localToGlobalDof, globalDofStart, globalDofEnd);
        
        // Create RHS and solution vectors
        HYPRE_Int ilower = static_cast<HYPRE_Int>(globalDofStart);
        HYPRE_Int iupper = static_cast<HYPRE_Int>(globalDofEnd - 1);

        HYPRE_IJVectorCreate(comm_, ilower, iupper, &b_hypre_);
        HYPRE_IJVectorCreate(comm_, ilower, iupper, &x_hypre_);
        HYPRE_IJVectorSetObjectType(b_hypre_, HYPRE_PARCSR);
        HYPRE_IJVectorSetObjectType(x_hypre_, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(b_hypre_);
        HYPRE_IJVectorInitialize(x_hypre_);

        // Set RHS and initial guess
        std::vector<HYPRE_Real> h_b(m);
        std::vector<HYPRE_Real> h_x(m);
        std::vector<HYPRE_Int> indices(m);
        
        thrust::copy(thrust::device_pointer_cast(b.data()),
                    thrust::device_pointer_cast(b.data() + m),
                    h_b.begin());
        thrust::copy(thrust::device_pointer_cast(x.data()),
                    thrust::device_pointer_cast(x.data() + m),
                    h_x.begin());
        
        // Validate input vectors for NaN/Inf
        bool bHasNaN = false, bHasInf = false, xHasNaN = false, xHasInf = false;
        for (HYPRE_Int i = 0; i < m; ++i) {
            if (std::isnan(h_b[i])) bHasNaN = true;
            if (std::isinf(h_b[i])) bHasInf = true;
            if (std::isnan(h_x[i])) xHasNaN = true;
            if (std::isinf(h_x[i])) xHasInf = true;
        }
        if ((bHasNaN || bHasInf || xHasNaN || xHasInf) && verbose_ && rank_ == 0) {
            if (bHasNaN) std::cerr << "ERROR: RHS vector contains NaN!" << std::endl;
            if (bHasInf) std::cerr << "ERROR: RHS vector contains Inf!" << std::endl;
            if (xHasNaN) std::cerr << "ERROR: Initial guess contains NaN!" << std::endl;
            if (xHasInf) std::cerr << "ERROR: Initial guess contains Inf!" << std::endl;
        }
        
        for (HYPRE_Int i = 0; i < m; ++i) {
            indices[i] = ilower + i;
        }
        
        HYPRE_IJVectorSetValues(b_hypre_, m, indices.data(), h_b.data());
        HYPRE_IJVectorSetValues(x_hypre_, m, indices.data(), h_x.data());

        HYPRE_IJVectorAssemble(b_hypre_);
        HYPRE_IJVectorAssemble(x_hypre_);

        HYPRE_IJVectorGetObject(b_hypre_, (void**)&par_b_);
        HYPRE_IJVectorGetObject(x_hypre_, (void**)&par_x_);
        
        // Create PCG solver with BoomerAMG preconditioner
        // Full system matrix (diagonal=1 for BCs) works well with BoomerAMG
        HYPRE_ParCSRPCGCreate(comm_, &solver_);
        HYPRE_PCGSetMaxIter(solver_, maxIter_);
        HYPRE_PCGSetTol(solver_, tolerance_);
        HYPRE_PCGSetPrintLevel(solver_, verbose_ ? 2 : 0);
        HYPRE_PCGSetTwoNorm(solver_, 1);  // Use 2-norm for residual
        
        // Setup BoomerAMG preconditioner
        HYPRE_Solver precond;
        HYPRE_BoomerAMGCreate(&precond);
        HYPRE_BoomerAMGSetCoarsenType(precond, 10);  // HMIS
        HYPRE_BoomerAMGSetRelaxType(precond, 6);     // Hybrid Gauss-Seidel
        HYPRE_BoomerAMGSetNumSweeps(precond, 1);
        HYPRE_BoomerAMGSetMaxLevels(precond, 25);
        HYPRE_BoomerAMGSetTol(precond, 0.0);
        HYPRE_BoomerAMGSetPrintLevel(precond, 0);
        
        HYPRE_PCGSetPrecond(solver_,
                           (HYPRE_Int (*)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector))HYPRE_BoomerAMGSolve,
                           (HYPRE_Int (*)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector))HYPRE_BoomerAMGSetup,
                           precond);
        
        // Setup and solve
        HYPRE_ParCSRPCGSetup(solver_, parcsr_A_, par_b_, par_x_);
        HYPRE_ParCSRPCGSolve(solver_, parcsr_A_, par_b_, par_x_);
        
        // Clean up
        HYPRE_BoomerAMGDestroy(precond);

        // Get solution
        HYPRE_IJVectorGetValues(x_hypre_, m, indices.data(), h_x.data());
        thrust::copy(h_x.begin(), h_x.end(),
                    thrust::device_pointer_cast(x.data()));

        // Check convergence
        int num_iterations;
        double final_res_norm;
        HYPRE_PCGGetNumIterations(solver_, &num_iterations);
        HYPRE_PCGGetFinalRelativeResidualNorm(solver_, &final_res_norm);

        if (verbose_ && rank_ == 0) {
            std::cout << "Hypre PCG+Jacobi converged in " << num_iterations
                      << " iterations, final residual: " << final_res_norm << std::endl;
        }

        return (final_res_norm < tolerance_);
    }

    void setVerbose(bool verbose) { verbose_ = verbose; }

private:
    void setupHypreMatrix(const Matrix& A, const std::vector<HYPRE_BigInt>& localToGlobalDof,
                         IndexType globalDofStart, IndexType globalDofEnd) {
        HYPRE_Int m = static_cast<HYPRE_Int>(A.numRows());
        HYPRE_Int ilower = static_cast<HYPRE_Int>(globalDofStart);
        HYPRE_Int iupper = static_cast<HYPRE_Int>(globalDofEnd - 1);

        // Get matrix data from device
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
        
        // Validate input matrix for NaN/Inf
        bool hasNaN = false, hasInf = false;
        for (size_t i = 0; i < h_values.size(); ++i) {
            if (std::isnan(h_values[i])) hasNaN = true;
            if (std::isinf(h_values[i])) hasInf = true;
        }
        if (hasNaN && verbose_ && rank_ == 0) {
            std::cerr << "ERROR: Input matrix contains NaN values!" << std::endl;
        }
        if (hasInf && verbose_ && rank_ == 0) {
            std::cerr << "ERROR: Input matrix contains Inf values!" << std::endl;
        }

        // Create Hypre matrix
        HYPRE_IJMatrixCreate(comm_, ilower, iupper, ilower, iupper, &A_hypre_);
        HYPRE_IJMatrixSetObjectType(A_hypre_, HYPRE_PARCSR);
        HYPRE_IJMatrixInitialize(A_hypre_);

        // Set matrix values row by row
        int totalValidEntries = 0;
        int totalSkippedEntries = 0;
        for (HYPRE_Int i = 0; i < m; ++i) {
            HYPRE_Int row = ilower + i;
            HYPRE_Int ncols = static_cast<HYPRE_Int>(h_rowOffsets[i+1] - h_rowOffsets[i]);
            std::vector<HYPRE_Int> cols;
            std::vector<HYPRE_Real> vals;
            cols.reserve(ncols);
            vals.reserve(ncols);
            
            for (HYPRE_Int j = 0; j < ncols; ++j) {
                IndexType localCol = h_colIndices[h_rowOffsets[i] + j];
                HYPRE_BigInt globalCol = -1;
                RealType val = h_values[h_rowOffsets[i] + j];
                
                // Validate value before adding
                if (std::isnan(val) || std::isinf(val)) {
                    if (verbose_ && rank_ == 0) {
                        std::cerr << "ERROR: Found invalid value at row " << i 
                                  << ", local col " << localCol << ": " << val << std::endl;
                    }
                    totalSkippedEntries++;
                    continue;  // Skip invalid entries
                }
                
                if (!localToGlobalDof.empty() && localCol < localToGlobalDof.size()) {
                    globalCol = localToGlobalDof[localCol];
                } else if (localCol < m) {
                    globalCol = ilower + static_cast<HYPRE_Int>(localCol);
                } else {
                    // Column index out of range - skip it
                    totalSkippedEntries++;
                    continue;
                }
                
                if (globalCol >= 0) {
                    cols.push_back(static_cast<HYPRE_Int>(globalCol));
                    vals.push_back(static_cast<HYPRE_Real>(val));
                    totalValidEntries++;
                }
            }
            
            if (!cols.empty()) {
                HYPRE_Int validNcols = static_cast<HYPRE_Int>(cols.size());
                HYPRE_IJMatrixSetValues(A_hypre_, 1, &validNcols, &row, cols.data(), vals.data());
            }
        }
        
        if (verbose_ && rank_ == 0) {
            std::cout << "Hypre matrix setup: " << totalValidEntries << " valid entries, " 
                      << totalSkippedEntries << " skipped entries" << std::endl;
        }

        HYPRE_IJMatrixAssemble(A_hypre_);
        HYPRE_IJMatrixGetObject(A_hypre_, (void**)&parcsr_A_);
    }

    void destroy() {
        if (solver_) {
            HYPRE_ParCSRPCGDestroy(solver_);
            solver_ = nullptr;
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
    int rank_;
    int maxIter_;
    RealType tolerance_;
    bool verbose_;

    HYPRE_Solver solver_;

    HYPRE_IJMatrix A_hypre_;
    HYPRE_ParCSRMatrix parcsr_A_;

    HYPRE_IJVector b_hypre_;
    HYPRE_IJVector x_hypre_;
    HYPRE_ParVector par_b_;
    HYPRE_ParVector par_x_;
};

} // namespace fem
} // namespace mars
