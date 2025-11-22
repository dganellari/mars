#pragma once

#include "mars_cg_solver.hpp"
#include "../fem/mars_unstructured_dof_handler.hpp"

namespace mars {
namespace fem {

/**
 * @brief Distributed Conjugate Gradient solver with automatic halo exchange
 * 
 * Wraps the standard CG solver to handle ghost DOF updates via MPI communication
 * before each matrix-vector multiplication.
 */
template<typename ElementTag, typename RealType, typename IndexType, typename AcceleratorTag>
class DistributedCGSolver {
public:
    using Matrix = SparseMatrix<IndexType, RealType, AcceleratorTag>;
    using Vector = typename mars::VectorSelector<RealType, AcceleratorTag>::type;
    using DofHandler = UnstructuredDofHandler<ElementTag, RealType, IndexType, AcceleratorTag>;
    using CGSolver = ConjugateGradientSolver<RealType, IndexType, AcceleratorTag>;
    
    DistributedCGSolver(DofHandler& dofHandler, int maxIter = 1000, RealType tolerance = 1e-10)
        : dofHandler_(dofHandler)
        , cgSolver_(maxIter, tolerance)
    {}
    
    void setVerbose(bool verbose) {
        cgSolver_.setVerbose(verbose);
    }
    
    /**
     * @brief Solve Ax = b with automatic halo exchange
     * 
     * Note: This currently assumes square matrix (owned DOFs only).
     * For rectangular matrix (owned + ghost columns), we need to:
     * 1. Extend solution vector to include ghost DOF slots
     * 2. Call halo exchange before each SpMV
     * 3. Extract owned portion after solve
     */
    bool solve(const Matrix& A, const Vector& b, Vector& x) {
        // For now, delegate to standard CG solver
        // TODO: Modify CG to handle rectangular matrices with halo exchange
        return cgSolver_.solve(A, b, x);
    }
    
    /**
     * @brief Solve with pre-allocated extended vectors (owned + ghost DOFs)
     * 
     * This version properly handles rectangular matrices by maintaining
     * ghost DOF values through halo exchange.
     * 
     * @param A Matrix (numLocalDofs × numLocalDofsWithGhosts)
     * @param b RHS vector (size numLocalDofs - owned DOFs only)
     * @param x Solution vector (size numLocalDofs - owned DOFs only)
     */
    bool solveWithHaloExchange(const Matrix& A, const Vector& b, Vector& x) {
        size_t numLocalDofs = dofHandler_.get_num_local_dofs();
        size_t numLocalDofsWithGhosts = dofHandler_.get_num_local_dofs_with_ghosts();
        
        // Extend solution vector to include ghost DOF slots
        Vector x_extended(numLocalDofsWithGhosts, 0.0);
        
        // Copy initial guess (if any) to owned portion
        if (x.size() >= numLocalDofs) {
            thrust::copy(thrust::device_pointer_cast(x.data()),
                        thrust::device_pointer_cast(x.data() + numLocalDofs),
                        thrust::device_pointer_cast(x_extended.data()));
        }
        
        // Initialize ghost values via halo exchange
        dofHandler_.updateGhostDofValues(x_extended);
        
        // TODO: Implement CG with halo exchange callbacks
        // For now, this is a placeholder
        std::cerr << "Error: solveWithHaloExchange not yet fully implemented\n";
        std::cerr << "       Need to modify CG solver to call halo exchange before each SpMV\n";
        
        // Extract owned portion back to output
        thrust::copy(thrust::device_pointer_cast(x_extended.data()),
                    thrust::device_pointer_cast(x_extended.data() + numLocalDofs),
                    thrust::device_pointer_cast(x.data()));
        
        return false;
    }
    
private:
    DofHandler& dofHandler_;
    CGSolver cgSolver_;
};

} // namespace fem
} // namespace mars
