#pragma once

#include "../fem/mars_sparse_matrix.hpp"
#include "mars_cg_solver_with_preconditioner.hpp"
#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>
#include <mpi.h>
#include <iostream>
#include <vector>

namespace mars {
namespace fem {

// Hypre AMG Preconditioner wrapper for MARS sparse matrices
template<typename RealType, typename IndexType, typename AcceleratorTag>
class HypreAMGPreconditioner : public Preconditioner<RealType, IndexType, AcceleratorTag> {
public:
    using Matrix = SparseMatrix<IndexType, RealType, AcceleratorTag>;
    using Vector = typename mars::VectorSelector<RealType, AcceleratorTag>::type;

    HypreAMGPreconditioner(MPI_Comm comm = MPI_COMM_WORLD)
        : comm_(comm), precond_(nullptr), A_hypre_(nullptr), b_hypre_(nullptr), x_hypre_(nullptr),
          ilower_(0), iupper_(0), m_(0) {
    }

    ~HypreAMGPreconditioner() {
        destroy();
    }

    // Setup the preconditioner with the matrix
    void setup(const Matrix& A) override {
        destroy();  // Clean up any existing setup

        m_ = A.numRows();  // Local rows (owned DOFs)
        int n = A.numCols();  // Local columns (owned + ghost DOFs)

        // For now, assume the matrix is square and we have all DOFs
        // In a full distributed implementation, we'd need to handle partitioning properly
        if (m_ != n) {
            std::cerr << "HypreAMGPreconditioner: Matrix must be square for now" << std::endl;
            return;
        }

        // Get matrix data from device
        std::vector<IndexType> h_rowOffsets(m_ + 1);
        std::vector<IndexType> h_colIndices(A.nnz());
        std::vector<RealType> h_values(A.nnz());

        // Use local variable for pointer arithmetic to avoid template issues
        auto row_ptr = A.rowOffsetsPtr();
        auto col_ptr = A.colIndicesPtr();
        auto val_ptr = A.valuesPtr();
        
        thrust::copy(thrust::device_pointer_cast(row_ptr),
                    thrust::device_pointer_cast(row_ptr + m_ + 1),
                    h_rowOffsets.begin());
        thrust::copy(thrust::device_pointer_cast(col_ptr),
                    thrust::device_pointer_cast(col_ptr + A.nnz()),
                    h_colIndices.begin());
        thrust::copy(thrust::device_pointer_cast(val_ptr),
                    thrust::device_pointer_cast(val_ptr + A.nnz()),
                    h_values.begin());

        // Create Hypre matrix
        // For simplicity, assume single rank for now (global partitioning)
        // In full implementation, need proper distributed matrix setup
        ilower_ = 0;  // Global row start (would be rank * m_ in distributed)
        iupper_ = m_ - 1;  // Global row end

        HYPRE_IJMatrixCreate(comm_, ilower_, iupper_, ilower_, iupper_, &A_hypre_);
        HYPRE_IJMatrixSetObjectType(A_hypre_, HYPRE_PARCSR);
        HYPRE_IJMatrixInitialize(A_hypre_);

        // Set matrix values
        for (int i = 0; i < m_; ++i) {
            int row = ilower + i;
            int ncols = h_rowOffsets[i+1] - h_rowOffsets[i];
            std::vector<int> cols(ncols);
            std::vector<double> vals(ncols);

            for (int j = 0; j < ncols; ++j) {
                cols[j] = h_colIndices[h_rowOffsets[i] + j];
                vals[j] = h_values[h_rowOffsets[i] + j];
            }

            HYPRE_IJMatrixSetValues(A_hypre_, 1, &ncols, &row, cols.data(), vals.data());
        }

        HYPRE_IJMatrixAssemble(A_hypre_);
        HYPRE_IJMatrixGetObject(A_hypre_, (void**)&parcsr_A_);

        // Create AMG preconditioner
        HYPRE_BoomerAMGCreate(&precond_);
        HYPRE_BoomerAMGSetPrintLevel(precond_, 0);  // No output
        HYPRE_BoomerAMGSetCoarsenType(precond_, 6);  // Falgout coarsening
        HYPRE_BoomerAMGSetRelaxType(precond_, 6);    // Sym GS/Hybrid relaxation
        HYPRE_BoomerAMGSetNumSweeps(precond_, 1);    // 1 sweep
        HYPRE_BoomerAMGSetMaxLevels(precond_, 20);   // Max levels
        HYPRE_BoomerAMGSetTol(precond_, 0.0);        // Tolerance for convergence

        HYPRE_BoomerAMGSetup(precond_, parcsr_A_, NULL, NULL);
    }

    // Apply preconditioner: solve M*z = r
    void apply(const Vector& r, Vector& z) override {
        if (!precond_) {
            std::cerr << "HypreAMGPreconditioner: Preconditioner not set up" << std::endl;
            return;
        }

        int m = r.size();
        std::vector<double> h_r(m_), h_z(m_);

        // Copy r to host
        thrust::copy(thrust::device_pointer_cast(r.data()),
                    thrust::device_pointer_cast(r.data() + m_),
                    h_r.begin());

        // Create Hypre vectors (temporary - in full implementation, reuse)
        HYPRE_IJVectorCreate(comm_, ilower_, iupper_, &b_hypre_);
        HYPRE_IJVectorCreate(comm_, ilower_, iupper_, &x_hypre_);
        HYPRE_IJVectorSetObjectType(b_hypre_, HYPRE_PARCSR);
        HYPRE_IJVectorSetObjectType(x_hypre_, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(b_hypre_);
        HYPRE_IJVectorInitialize(x_hypre_);

        // Set vector values
        HYPRE_IJVectorSetValues(b_hypre_, m_, NULL, h_r.data());
        HYPRE_IJVectorAssemble(b_hypre_);

        HYPRE_IJVectorGetObject(b_hypre_, (void**)&par_b_);
        HYPRE_IJVectorGetObject(x_hypre_, (void**)&par_x_);

        // Apply preconditioner
        HYPRE_BoomerAMGSolve(precond_, parcsr_A_, par_b_, par_x_);

        // Get solution
        HYPRE_IJVectorGetValues(x_hypre_, m_, NULL, h_z.data());

        // Copy z back to device
        thrust::copy(h_z.begin(), h_z.end(),
                    thrust::device_pointer_cast(z.data()));

        // Clean up temporary vectors
        HYPRE_IJVectorDestroy(b_hypre_);
        HYPRE_IJVectorDestroy(x_hypre_);
        b_hypre_ = nullptr;
        x_hypre_ = nullptr;
    }

private:
    void destroy() {
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
    HYPRE_Solver precond_;
    HYPRE_Solver solver_;  // For potential future use

    HYPRE_IJMatrix A_hypre_;
    HYPRE_ParCSRMatrix parcsr_A_;

    HYPRE_IJVector b_hypre_;
    HYPRE_IJVector x_hypre_;
    HYPRE_ParVector par_b_;
    HYPRE_ParVector par_x_;

    // Matrix dimensions
    int ilower_, iupper_, m_;
};

} // namespace fem
} // namespace mars