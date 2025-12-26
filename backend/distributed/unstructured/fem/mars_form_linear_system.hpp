#ifndef MARS_FORM_LINEAR_SYSTEM_HPP
#define MARS_FORM_LINEAR_SYSTEM_HPP

#include <vector>
#include <unordered_set>
#include "mars_base.hpp"
#include "mars_sparse_matrix.hpp"

namespace mars {
namespace fem {

/**
 * @brief Form reduced linear system by eliminating essential (boundary) DOFs.
 * 
 * This function mimics MFEM's FormLinearSystem operation:
 * - Extracts non-essential DOF rows/columns from matrix A
 * - Modifies RHS: b_interior -= A_interior_boundary * u_boundary (u_boundary assumed 0 for homogeneous BCs)
 * - Creates reduced system A_r * x_r = b_r with only interior DOFs
 * - Preserves original local column indexing for parallel coupling (needed for BoomerAMG)
 * 
 * @param A_full Full stiffness matrix (owned DOFs x all local DOFs including ghosts)
 * @param b_full Full RHS vector (owned DOFs only)
 * @param ess_dofs List of essential (boundary) DOF indices to eliminate
 * @param A_reduced Output: reduced matrix (interior owned x all interior+ghost)
 * @param b_reduced Output: reduced RHS vector (interior owned DOFs only)
 * @param dof_mapping Output: mapping from full DOF indices to reduced indices (owned DOFs only)
 * @param numGhostDofs Number of ghost DOFs in the full system
 * @param ghost_dof_mapping Output: mapping from full ghost DOF indices to reduced ghost indices (ghost DOFs only)
 * @param reducedLocalCount Output: total number of columns in reduced system (owned + interior ghosts)
 * @return Number of interior owned DOFs in reduced system
 */
template<typename KeyType, typename RealType, typename AcceleratorTag>
KeyType FormLinearSystem(
    const SparseMatrix<KeyType, RealType, AcceleratorTag>& A_full,
    const typename mars::VectorSelector<RealType, AcceleratorTag>::type& b_full,
    const std::vector<KeyType>& ess_dofs,
    SparseMatrix<KeyType, RealType, AcceleratorTag>& A_reduced,
    typename mars::VectorSelector<RealType, AcceleratorTag>::type& b_reduced,
    std::vector<KeyType>& dof_mapping,
    KeyType numGhostDofs = 0,
    const std::vector<KeyType>& ess_ghost_dofs = {},
    std::vector<KeyType>* ghost_dof_mapping = nullptr,
    KeyType* reducedLocalCount = nullptr)
{
    using Vector = typename mars::VectorSelector<RealType, AcceleratorTag>::type;
    
    KeyType numOwnedDofs = b_full.size();  // RHS only has owned DOFs
    KeyType numLocalDofs = A_full.numCols();  // Matrix has owned + ghost columns
    
    // Mark essential DOFs (only in owned range)
    std::vector<int> ess_marker(numOwnedDofs, 0);
    for (auto dof : ess_dofs) {
        if (dof < numOwnedDofs) {
            ess_marker[dof] = 1;
        }
    }
    
    // Copy matrix and RHS to host for processing
    std::vector<KeyType> h_rowOffsets(A_full.numRows() + 1);
    std::vector<KeyType> h_colIndices(A_full.nnz());
    std::vector<RealType> h_values(A_full.nnz());
    std::vector<RealType> h_b_full(numOwnedDofs);
    
    thrust::copy(thrust::device_pointer_cast(A_full.rowOffsetsPtr()),
                 thrust::device_pointer_cast(A_full.rowOffsetsPtr() + A_full.numRows() + 1),
                 h_rowOffsets.begin());
    thrust::copy(thrust::device_pointer_cast(A_full.colIndicesPtr()),
                 thrust::device_pointer_cast(A_full.colIndicesPtr() + A_full.nnz()),
                 h_colIndices.begin());
    thrust::copy(thrust::device_pointer_cast(A_full.valuesPtr()),
                 thrust::device_pointer_cast(A_full.valuesPtr() + A_full.nnz()),
                 h_values.begin());
    thrust::copy(thrust::device_pointer_cast(b_full.data()),
                 thrust::device_pointer_cast(b_full.data() + numOwnedDofs),
                 h_b_full.begin());
    
    // First pass: identify rows with zero diagonals (these are DOFs that are not in any local element)
    // This is a symptom of an incorrect ownership definition. A proper fix is in the DofHandler/ElementDomain.
    // This check is kept as a safeguard.
    std::vector<int> zerodiag_marker(numOwnedDofs, 0);
    int zeroDiagsFound = 0;
    for (KeyType i = 0; i < numOwnedDofs; ++i) {
        if (ess_marker[i]) continue;  // Boundary DOFs already marked
        
        // Check if diagonal exists and is non-zero
        bool hasNonZeroDiag = false;
        for (KeyType idx = h_rowOffsets[i]; idx < h_rowOffsets[i + 1]; ++idx) {
            if (h_colIndices[idx] == i) {
                if (std::abs(h_values[idx]) > 1e-14) {
                    hasNonZeroDiag = true;
                }
                break;
            }
        }
        
        if (!hasNonZeroDiag) {
            zerodiag_marker[i] = 1;
            zeroDiagsFound++;
        }
    }
    
    if (zeroDiagsFound > 0) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        std::cout << "Rank " << rank << ": WARNING - Found " << zeroDiagsFound 
                  << " owned DOFs with zero/missing diagonals. These are likely ghost nodes incorrectly marked as owned." << std::endl;
    }
    
    // Create mapping: full owned DOF index -> reduced owned DOF index
    // Ghost DOFs will be remapped separately
    dof_mapping.resize(numOwnedDofs);
    KeyType reducedOwnedCount = 0;
    for (KeyType i = 0; i < numOwnedDofs; ++i) {
        if (!ess_marker[i] && !zerodiag_marker[i]) {
            dof_mapping[i] = reducedOwnedCount++;
        } else {
            dof_mapping[i] = static_cast<KeyType>(-1);
        }
    }
    
    // Mark boundary ghost DOFs
    std::vector<int> ess_ghost_marker(numGhostDofs, 0);
    for (auto ghostDof : ess_ghost_dofs) {
        if (ghostDof < numGhostDofs) {
            ess_ghost_marker[ghostDof] = 1;
        }
    }
    
    // Count interior ghost DOFs
    KeyType numInteriorGhostDofs = 0;
    for (KeyType g = 0; g < numGhostDofs; ++g) {
        if (!ess_ghost_marker[g]) {
            numInteriorGhostDofs++;
        }
    }
    
    // Build reduced system (rectangular: owned rows × owned+interior_ghost columns)
    KeyType totalReducedCols = reducedOwnedCount + numInteriorGhostDofs;

    std::vector<KeyType> h_rowOffsets_r(reducedOwnedCount + 1, 0);
    std::vector<KeyType> h_colIndices_r;
    std::vector<RealType> h_values_r;
    std::vector<RealType> h_b_r(reducedOwnedCount);

    h_colIndices_r.reserve(A_full.nnz());
    h_values_r.reserve(A_full.nnz());
    
    // Create column index mapping: old local -> new local
    std::vector<KeyType> colMapping(numLocalDofs, static_cast<KeyType>(-1));
    for (KeyType j = 0; j < numOwnedDofs; ++j) {
        if (!ess_marker[j]) {
            colMapping[j] = dof_mapping[j];  // Interior owned DOF
        }
    }
    // Map interior ghost DOFs only (skip boundary ghosts)
    KeyType reducedGhostIdx = 0;
    std::vector<KeyType> ghostMapping(numGhostDofs, static_cast<KeyType>(-1));
    for (KeyType g = 0; g < numGhostDofs; ++g) {
        if (!ess_ghost_marker[g]) {
            KeyType j = numOwnedDofs + g;
            colMapping[j] = reducedOwnedCount + reducedGhostIdx;
            ghostMapping[g] = reducedGhostIdx;  // Store mapping for ghost DOFs
            reducedGhostIdx++;
        }
    }
    
    KeyType newRow = 0;
    int ghostColsKept = 0;
    
    for (KeyType i = 0; i < numOwnedDofs; ++i) {
        if (dof_mapping[i] == static_cast<KeyType>(-1)) continue;  // Skip eliminated DOFs (boundary or zero-diag)
        
        RealType rhs_i = h_b_full[i];
        
        // Extract non-essential columns
        for (KeyType idx = h_rowOffsets[i]; idx < h_rowOffsets[i + 1]; ++idx) {
            KeyType j = h_colIndices[idx];
            RealType A_ij = h_values[idx];
            
            if (j < numOwnedDofs) {
                // Owned column
                if (dof_mapping[j] != static_cast<KeyType>(-1)) {
                    // Interior-interior: keep in reduced matrix
                    h_colIndices_r.push_back(colMapping[j]);
                    h_values_r.push_back(A_ij);
                }
                // Boundary columns are eliminated (modify RHS if non-homogeneous)
            } else {
                // Ghost column: keep only if interior (skip boundary ghosts)
                if (colMapping[j] != static_cast<KeyType>(-1)) {
                    h_colIndices_r.push_back(colMapping[j]);
                    h_values_r.push_back(A_ij);
                    ghostColsKept++;
                }
            }
        }
        
        h_b_r[newRow] = rhs_i;
        h_rowOffsets_r[newRow + 1] = h_colIndices_r.size();
        newRow++;
    }
    
    if (ghostColsKept == 0 && numGhostDofs > 0) {
        std::cout << "WARNING: Expected " << numGhostDofs << " ghost DOFs but found NO ghost columns in reduced matrix!\n";
        std::cout << "         This indicates no inter-rank coupling to interior DOFs.\n";
    } else if (ghostColsKept > 0) {
        // Ghost columns preserved for inter-rank coupling
    }
    
    // Copy reduced system to device
    A_reduced.allocate(reducedOwnedCount, totalReducedCols, h_values_r.size());

    // Validate column indices
    KeyType maxCol = 0, minCol = totalReducedCols;
    int outOfRange = 0, missingDiag = 0;
    std::vector<bool> hasDiag(reducedOwnedCount, false);

    for (KeyType row = 0; row < reducedOwnedCount; ++row) {
        for (KeyType idx = h_rowOffsets_r[row]; idx < h_rowOffsets_r[row + 1]; ++idx) {
            auto col = h_colIndices_r[idx];
            if (col >= totalReducedCols || col < 0) outOfRange++;
            if (col == row) hasDiag[row] = true;
            maxCol = std::max(maxCol, col);
            minCol = std::min(minCol, col);
        }
    }
    for (auto d : hasDiag) if (!d) missingDiag++;

    std::cout << "Reduced matrix: " << reducedOwnedCount << "x" << totalReducedCols
              << " (owned x local), nnz=" << h_values_r.size()
              << ", col range=[" << minCol << ", " << maxCol << "]"
              << ", out_of_range=" << outOfRange
              << ", missing_diag=" << missingDiag << std::endl;
    
    thrust::copy(h_rowOffsets_r.begin(), h_rowOffsets_r.end(),
                 thrust::device_pointer_cast(A_reduced.rowOffsetsPtr()));
    thrust::copy(h_colIndices_r.begin(), h_colIndices_r.end(),
                 thrust::device_pointer_cast(A_reduced.colIndicesPtr()));
    thrust::copy(h_values_r.begin(), h_values_r.end(),
                 thrust::device_pointer_cast(A_reduced.valuesPtr()));
    
    b_reduced.resize(reducedOwnedCount);
    thrust::copy(h_b_r.begin(), h_b_r.end(),
                 thrust::device_pointer_cast(b_reduced.data()));
    
    // Diagnostic: check reduced RHS
    double sum_b = 0.0, min_b = 1e100, max_b = -1e100;
    int nan_count = 0;
    for (KeyType i = 0; i < reducedOwnedCount; ++i) {
        if (!std::isfinite(h_b_r[i])) {
            nan_count++;
        } else {
            sum_b += h_b_r[i];
            min_b = std::min(min_b, static_cast<double>(h_b_r[i]));
            max_b = std::max(max_b, static_cast<double>(h_b_r[i]));
        }
    }
    std::cout << "Reduced RHS: sum=" << sum_b << ", range=[" << min_b << ", " << max_b
              << "], NaNs=" << nan_count << "/" << reducedOwnedCount << std::endl;

    // Return optional outputs
    if (ghost_dof_mapping) {
        *ghost_dof_mapping = ghostMapping;
    }
    if (reducedLocalCount) {
        *reducedLocalCount = reducedOwnedCount + numInteriorGhostDofs;
    }

    return reducedOwnedCount;
}

/**
 * @brief Recover full FEM solution from reduced solution.
 * 
 * Expands the solution of the reduced system back to the full DOF space,
 * setting essential DOFs to their prescribed values (0 for homogeneous BCs).
 * 
 * @param x_reduced Solution of reduced system (interior DOFs only)
 * @param dof_mapping Mapping from full DOF indices to reduced indices
 * @param x_full Output: full solution vector (essential DOFs = 0)
 */
template<typename KeyType, typename RealType, typename AcceleratorTag>
void RecoverFEMSolution(
    const typename mars::VectorSelector<RealType, AcceleratorTag>::type& x_reduced,
    const std::vector<KeyType>& dof_mapping,
    typename mars::VectorSelector<RealType, AcceleratorTag>::type& x_full)
{
    using Vector = typename mars::VectorSelector<RealType, AcceleratorTag>::type;
    
    KeyType numOwnedDofs = dof_mapping.size();
    KeyType reducedSize = x_reduced.size();
    
    // Copy reduced solution to host
    std::vector<RealType> h_x_reduced(reducedSize);
    thrust::copy(thrust::device_pointer_cast(x_reduced.data()),
                 thrust::device_pointer_cast(x_reduced.data() + reducedSize),
                 h_x_reduced.begin());
    
    // Expand to full size
    std::vector<RealType> h_x_full(numOwnedDofs, 0.0);
    for (KeyType i = 0; i < numOwnedDofs; ++i) {
        KeyType reduced_idx = dof_mapping[i];
        if (reduced_idx != static_cast<KeyType>(-1)) {
            h_x_full[i] = h_x_reduced[reduced_idx];
        }
        // else: essential DOF, already 0
    }
    
    // Copy back to device
    x_full.resize(numOwnedDofs);
    thrust::copy(h_x_full.begin(), h_x_full.end(),
                 thrust::device_pointer_cast(x_full.data()));
}

} // namespace fem
} // namespace mars

#endif // MARS_FORM_LINEAR_SYSTEM_HPP
