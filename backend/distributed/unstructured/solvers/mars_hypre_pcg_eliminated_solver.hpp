#pragma once

#include "../fem/mars_sparse_matrix.hpp"
#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>
#include <mpi.h>
#include <iostream>
#include <vector>
#include <unordered_set>

namespace mars {
namespace fem {

// Hypre PCG solver with TRUE DOF ELIMINATION (like MFEM)
// This solver eliminates boundary DOFs from the system entirely, creating a reduced system.
// The reduced system is then solved with Hypre PCG + BoomerAMG.
template<typename RealType, typename IndexType, typename AcceleratorTag>
class HyprePCGEliminatedSolver {
public:
    using Matrix = SparseMatrix<IndexType, RealType, AcceleratorTag>;
    using Vector = typename mars::VectorSelector<RealType, AcceleratorTag>::type;

    HyprePCGEliminatedSolver(MPI_Comm comm = MPI_COMM_WORLD, int maxIter = 1000, RealType tolerance = 1e-6)
        : comm_(comm), maxIter_(maxIter), tolerance_(tolerance),
          solver_(nullptr), precond_(nullptr), A_hypre_(nullptr),
          b_hypre_(nullptr), x_hypre_(nullptr), verbose_(true) {
        MPI_Comm_rank(comm_, &rank_);
    }

    ~HyprePCGEliminatedSolver() {
        destroy();
    }

    // Solve Ax = b with boundary condition elimination
    // boundaryDofs: list of LOCAL DOF indices that are on the boundary (to be eliminated)
    // boundaryValues: values for boundary DOFs (typically 0 for homogeneous Dirichlet)
    template<typename KeyType>
    bool solve(const Matrix& A_full, const Vector& b_full, Vector& x_full,
               IndexType globalDofStart, IndexType globalDofEnd,
               const std::vector<IndexType>& boundaryDofs,
               const std::vector<RealType>& boundaryValues,
               const std::vector<KeyType>& nodeToLocalDof,
               const std::vector<KeyType>& nodeToGlobalDof) {
        
        destroy();  // Clean up any existing setup

        HYPRE_Int m_full = static_cast<HYPRE_Int>(A_full.numRows());
        
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
        
        // Create boundary DOF set for quick lookup
        std::unordered_set<IndexType> boundarySet(boundaryDofs.begin(), boundaryDofs.end());
        
        if (verbose_ && rank_ == 0) {
            std::cout << "HyprePCGEliminatedSolver: Eliminating " << boundaryDofs.size() 
                      << " boundary DOFs from " << m_full << " total DOFs..." << std::endl;
        }
        
        // Create mapping from full DOF index to reduced DOF index
        // -1 means this DOF is eliminated (boundary)
        std::vector<IndexType> fullToReduced(m_full, -1);
        std::vector<IndexType> reducedToFull;
        reducedToFull.reserve(m_full - boundaryDofs.size());
        
        IndexType reducedIdx = 0;
        for (IndexType i = 0; i < m_full; ++i) {
            if (boundarySet.count(i) == 0) {
                // This is a free (non-boundary) DOF
                fullToReduced[i] = reducedIdx;
                reducedToFull.push_back(i);
                reducedIdx++;
            }
        }
        
        HYPRE_Int m_reduced = static_cast<HYPRE_Int>(reducedToFull.size());
        
        if (verbose_ && rank_ == 0) {
            std::cout << "  Reduced system size: " << m_reduced << " (was " << m_full << ")" << std::endl;
        }
        
        // Get full matrix data from device
        std::vector<IndexType> h_rowOffsets_full(m_full + 1);
        std::vector<IndexType> h_colIndices_full(A_full.nnz());
        std::vector<RealType> h_values_full(A_full.nnz());

        thrust::copy(thrust::device_pointer_cast(A_full.rowOffsetsPtr()),
                    thrust::device_pointer_cast(A_full.rowOffsetsPtr() + m_full + 1),
                    h_rowOffsets_full.begin());
        thrust::copy(thrust::device_pointer_cast(A_full.colIndicesPtr()),
                    thrust::device_pointer_cast(A_full.colIndicesPtr() + A_full.nnz()),
                    h_colIndices_full.begin());
        thrust::copy(thrust::device_pointer_cast(A_full.valuesPtr()),
                    thrust::device_pointer_cast(A_full.valuesPtr() + A_full.nnz()),
                    h_values_full.begin());
        
        // Get full RHS and solution vectors
        std::vector<RealType> h_b_full(m_full);
        thrust::copy(thrust::device_pointer_cast(b_full.data()),
                    thrust::device_pointer_cast(b_full.data() + m_full),
                    h_b_full.begin());
        
        // Extract reduced system A_red and b_red
        // A_red only contains rows and columns for non-boundary DOFs
        // b_red is adjusted for known boundary values: b_red = b_free - A_fb * x_boundary
        
        std::vector<RealType> h_b_reduced(m_reduced);
        std::vector<std::vector<IndexType>> reducedCols(m_reduced);
        std::vector<std::vector<RealType>> reducedVals(m_reduced);
        
        // Build reduced RHS and matrix row by row
        for (HYPRE_Int i_red = 0; i_red < m_reduced; ++i_red) {
            IndexType i_full = reducedToFull[i_red];
            
            // Start with original RHS value
            RealType rhs = h_b_full[i_full];
            
            // Process this row: separate free vs boundary columns
            for (IndexType idx = h_rowOffsets_full[i_full]; idx < h_rowOffsets_full[i_full + 1]; ++idx) {
                IndexType j_full_local = h_colIndices_full[idx];
                RealType value = h_values_full[idx];
                
                // Get global column index
                HYPRE_BigInt j_global = -1;
                if (!localToGlobalDof.empty() && j_full_local < localToGlobalDof.size()) {
                    j_global = localToGlobalDof[j_full_local];
                } else if (j_full_local < m_full) {
                    // Owned DOF without mapping
                    j_global = globalDofStart + static_cast<HYPRE_BigInt>(j_full_local);
                }
                
                // Skip invalid columns
                if (j_global < 0) {
                    continue;
                }
                
                // Check if this column corresponds to a boundary DOF
                bool isColBoundary = (j_full_local < m_full) && (boundarySet.count(j_full_local) > 0);
                
                if (!isColBoundary) {
                    // Free DOF column - keep it in reduced matrix
                    // Store global column index for now
                    reducedCols[i_red].push_back(j_global);
                    reducedVals[i_red].push_back(value);
                } else {
                    // Boundary DOF column - eliminate it
                    // Adjust RHS: b_red[i] -= A[i,j] * x_boundary[j]
                    RealType boundaryValue = 0.0;
                    if (!boundaryValues.empty() && j_full_local < boundaryValues.size()) {
                        boundaryValue = boundaryValues[j_full_local];
                    }
                    rhs -= value * boundaryValue;
                }
            }
            
            h_b_reduced[i_red] = rhs;
        }
        
        // Now build Hypre matrix from reduced system
        // The column indices are already in GLOBAL numbering
        HYPRE_Int ilower = static_cast<HYPRE_Int>(globalDofStart);
        HYPRE_Int iupper = static_cast<HYPRE_Int>(globalDofStart + m_reduced - 1);
        
        // Create Hypre matrix
        HYPRE_IJMatrixCreate(comm_, ilower, iupper, ilower, iupper, &A_hypre_);
        HYPRE_IJMatrixSetObjectType(A_hypre_, HYPRE_PARCSR);
        HYPRE_IJMatrixInitialize(A_hypre_);
        
        // Set matrix values
        for (HYPRE_Int i = 0; i < m_reduced; ++i) {
            HYPRE_Int row = ilower + i;
            HYPRE_Int ncols = static_cast<HYPRE_Int>(reducedCols[i].size());
            
            if (ncols > 0) {
                // Column indices are already global
                std::vector<HYPRE_Int> cols_global(ncols);
                for (HYPRE_Int j = 0; j < ncols; ++j) {
                    cols_global[j] = static_cast<HYPRE_Int>(reducedCols[i][j]);
                }
                
                HYPRE_IJMatrixSetValues(A_hypre_, 1, &ncols, &row, cols_global.data(), 
                                       reinterpret_cast<HYPRE_Real*>(reducedVals[i].data()));
            }
        }
        
        HYPRE_IJMatrixAssemble(A_hypre_);
        HYPRE_IJMatrixGetObject(A_hypre_, (void**)&parcsr_A_);
        
        // Create RHS and solution vectors
        HYPRE_IJVectorCreate(comm_, ilower, iupper, &b_hypre_);
        HYPRE_IJVectorCreate(comm_, ilower, iupper, &x_hypre_);
        HYPRE_IJVectorSetObjectType(b_hypre_, HYPRE_PARCSR);
        HYPRE_IJVectorSetObjectType(x_hypre_, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(b_hypre_);
        HYPRE_IJVectorInitialize(x_hypre_);
        
        // Set RHS values
        std::vector<HYPRE_Int> indices(m_reduced);
        for (HYPRE_Int i = 0; i < m_reduced; ++i) {
            indices[i] = ilower + i;
        }
        HYPRE_IJVectorSetValues(b_hypre_, m_reduced, indices.data(), 
                               reinterpret_cast<HYPRE_Real*>(h_b_reduced.data()));
        
        // Set initial guess (zero)
        std::vector<HYPRE_Real> h_x_reduced(m_reduced, 0.0);
        HYPRE_IJVectorSetValues(x_hypre_, m_reduced, indices.data(), h_x_reduced.data());
        
        HYPRE_IJVectorAssemble(b_hypre_);
        HYPRE_IJVectorAssemble(x_hypre_);
        
        HYPRE_IJVectorGetObject(b_hypre_, (void**)&par_b_);
        HYPRE_IJVectorGetObject(x_hypre_, (void**)&par_x_);
        
        // Create BoomerAMG preconditioner
        HYPRE_BoomerAMGCreate(&precond_);
        HYPRE_BoomerAMGSetPrintLevel(precond_, 0);
        HYPRE_BoomerAMGSetCoarsenType(precond_, 6); // Falgout coarsening
        HYPRE_BoomerAMGSetRelaxType(precond_, 3);   // Hybrid Gauss-Seidel
        HYPRE_BoomerAMGSetNumSweeps(precond_, 1);
        HYPRE_BoomerAMGSetMaxLevels(precond_, 25);
        HYPRE_BoomerAMGSetTol(precond_, 0.0);       // No convergence check in preconditioner
        
        // Create PCG solver
        HYPRE_ParCSRPCGCreate(comm_, &solver_);
        HYPRE_PCGSetMaxIter(solver_, maxIter_);
        HYPRE_PCGSetTol(solver_, tolerance_);
        HYPRE_PCGSetPrintLevel(solver_, verbose_ ? 2 : 0);
        
        // Set preconditioner
        HYPRE_PCGSetPrecond(solver_,
                           (HYPRE_Int (*)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector))HYPRE_BoomerAMGSolve,
                           (HYPRE_Int (*)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector))HYPRE_BoomerAMGSetup,
                           precond_);
        
        if (verbose_ && rank_ == 0) {
            std::cout << "  Solving reduced system with Hypre PCG + BoomerAMG..." << std::endl;
        }
        
        // Setup and solve
        HYPRE_ParCSRPCGSetup(solver_, parcsr_A_, par_b_, par_x_);
        HYPRE_ParCSRPCGSolve(solver_, parcsr_A_, par_b_, par_x_);
        
        // Get solution for reduced system
        HYPRE_IJVectorGetValues(x_hypre_, m_reduced, indices.data(), h_x_reduced.data());
        
        // Recover full solution by inserting boundary values
        std::vector<RealType> h_x_full(m_full);
        for (HYPRE_Int i_red = 0; i_red < m_reduced; ++i_red) {
            IndexType i_full = reducedToFull[i_red];
            h_x_full[i_full] = h_x_reduced[i_red];
        }
        
        // Set boundary DOF values
        for (size_t i = 0; i < boundaryDofs.size(); ++i) {
            IndexType dof = boundaryDofs[i];
            RealType value = 0.0;
            if (!boundaryValues.empty() && i < boundaryValues.size()) {
                value = boundaryValues[i];
            }
            h_x_full[dof] = value;
        }
        
        // Copy solution back to device
        thrust::copy(h_x_full.begin(), h_x_full.end(),
                    thrust::device_pointer_cast(x_full.data()));
        
        // Check convergence
        int num_iterations;
        double final_res_norm;
        HYPRE_PCGGetNumIterations(solver_, &num_iterations);
        HYPRE_PCGGetFinalRelativeResidualNorm(solver_, &final_res_norm);

        if (verbose_ && rank_ == 0) {
            std::cout << "  Hypre PCG converged in " << num_iterations
                      << " iterations, final residual: " << final_res_norm << std::endl;
        }

        return (final_res_norm < tolerance_);
    }

    void setVerbose(bool verbose) { verbose_ = verbose; }

private:
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
    int rank_;
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
