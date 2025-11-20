#pragma once

#include "mars_sparse_matrix.hpp"
#include <thrust/device_vector.h>
#include <thrust/sequence.h>
#include <thrust/remove.h>
#include <thrust/copy.h>
#include <iostream>
#include <vector>

namespace mars
{
namespace fem
{

// DOF elimination for essential boundary conditions
// Removes boundary DOFs from the system before solving, then reconstructs full solution
template<typename RealType, typename IndexType, typename AcceleratorTag>
class DOFElimination
{
public:
    using Matrix = SparseMatrix<IndexType, RealType, AcceleratorTag>;
    using Vector = typename mars::VectorSelector<RealType, AcceleratorTag>::type;
    using IntVector = typename mars::VectorSelector<IndexType, AcceleratorTag>::type;

    // Build interior-only system from full system
    // Input: Full matrix A (n x n), full RHS b (n), boundary DOF flags, boundary values
    // Output: Reduced matrix A_int (n_int x n_int), reduced RHS b_int (n_int)
    void buildInteriorSystem(const Matrix& A_full,
                            const Vector& b_full,
                            const std::vector<bool>& isBoundaryDOF,
                            const std::vector<RealType>& boundaryValues,
                            Matrix& A_int,
                            Vector& b_int)
    {
        size_t n = A_full.numRows();
        
        // Build global-to-interior DOF map
        std::vector<IndexType> globalToInterior(n);
        std::vector<IndexType> interiorToGlobal;
        interiorToGlobal.reserve(n);
        
        IndexType interiorIdx = 0;
        for (size_t i = 0; i < n; ++i) {
            if (!isBoundaryDOF[i]) {
                globalToInterior[i] = interiorIdx++;
                interiorToGlobal.push_back(i);
            } else {
                globalToInterior[i] = static_cast<IndexType>(-1);  // Mark as boundary
            }
        }
        
        size_t n_int = interiorToGlobal.size();
        
        if (n_int == 0) {
            std::cout << "Warning: All DOFs are boundary DOFs - trivial system\n";
            A_int.allocate(0, 0, 0);
            b_int.resize(0);
            return;
        }
        
        std::cout << "   Eliminating BCs: " << n << " DOFs -> " << n_int << " interior DOFs\n";
        
        // Diagnostics
        RealType minBoundaryValue = std::numeric_limits<RealType>::max();
        RealType maxBoundaryValue = std::numeric_limits<RealType>::lowest();
        for (size_t i = 0; i < n; ++i) {
            if (isBoundaryDOF[i]) {
                minBoundaryValue = std::min(minBoundaryValue, boundaryValues[i]);
                maxBoundaryValue = std::max(maxBoundaryValue, boundaryValues[i]);
            }
        }
        std::cout << "   Boundary value range: [" << minBoundaryValue << ", " << maxBoundaryValue << "]\n";
        
        // Copy full matrix to host for processing
        std::vector<IndexType> h_rowOffsets(n + 1);
        std::vector<IndexType> h_colIndices(A_full.nnz());
        std::vector<RealType> h_values(A_full.nnz());
        
        thrust::copy(thrust::device_pointer_cast(A_full.rowOffsetsPtr()),
                    thrust::device_pointer_cast(A_full.rowOffsetsPtr() + n + 1),
                    h_rowOffsets.begin());
        thrust::copy(thrust::device_pointer_cast(A_full.colIndicesPtr()),
                    thrust::device_pointer_cast(A_full.colIndicesPtr() + A_full.nnz()),
                    h_colIndices.begin());
        thrust::copy(thrust::device_pointer_cast(A_full.valuesPtr()),
                    thrust::device_pointer_cast(A_full.valuesPtr() + A_full.nnz()),
                    h_values.begin());
        
        // Copy RHS to host
        std::vector<RealType> h_b_full(n);
        thrust::copy(thrust::device_pointer_cast(b_full.data()),
                    thrust::device_pointer_cast(b_full.data() + n),
                    h_b_full.begin());
        
        // Build interior-only matrix
        std::vector<IndexType> h_rowOffsets_int(n_int + 1, 0);
        std::vector<IndexType> h_colIndices_int;
        std::vector<RealType> h_values_int;
        std::vector<RealType> h_b_int(n_int, 0.0);
        
        h_colIndices_int.reserve(A_full.nnz() / 2);  // Estimate
        h_values_int.reserve(A_full.nnz() / 2);
        
        IndexType nnz_int = 0;
        
        for (size_t i_int = 0; i_int < n_int; ++i_int) {
            IndexType i_global = interiorToGlobal[i_int];
            
            // Start of this row
            h_rowOffsets_int[i_int] = nnz_int;
            
            // Initialize RHS from full system
            h_b_int[i_int] = h_b_full[i_global];
            
            // Process row i_global
            for (IndexType idx = h_rowOffsets[i_global]; idx < h_rowOffsets[i_global + 1]; ++idx) {
                IndexType j_global = h_colIndices[idx];
                RealType val = h_values[idx];
                
                if (!isBoundaryDOF[j_global]) {
                    // Interior-interior coupling: add to A_int
                    IndexType j_int = globalToInterior[j_global];
                    h_colIndices_int.push_back(j_int);
                    h_values_int.push_back(val);
                    ++nnz_int;
                } else {
                    // Interior-boundary coupling: move to RHS
                    // b_int[i] -= A[i,j] * u_boundary[j]
                    h_b_int[i_int] -= val * boundaryValues[j_global];
                }
            }
        }
        h_rowOffsets_int[n_int] = nnz_int;
        
        std::cout << "   Interior matrix: " << n_int << " x " << n_int 
                  << ", nnz = " << nnz_int << "\n";
        
        // Check interior RHS statistics
        RealType b_int_min = *std::min_element(h_b_int.begin(), h_b_int.end());
        RealType b_int_max = *std::max_element(h_b_int.begin(), h_b_int.end());
        RealType b_int_norm = 0.0;
        for (RealType val : h_b_int) b_int_norm += val * val;
        b_int_norm = std::sqrt(b_int_norm);
        std::cout << "   Interior RHS: ||b_int|| = " << b_int_norm 
                  << ", range [" << b_int_min << ", " << b_int_max << "]\n";
        
        // Check interior matrix diagonal
        std::vector<RealType> diag_int(n_int, 0.0);
        for (size_t i_int = 0; i_int < n_int; ++i_int) {
            for (IndexType idx = h_rowOffsets_int[i_int]; idx < h_rowOffsets_int[i_int + 1]; ++idx) {
                if (h_colIndices_int[idx] == static_cast<IndexType>(i_int)) {
                    diag_int[i_int] = h_values_int[idx];
                    break;
                }
            }
        }
        RealType diag_min = *std::min_element(diag_int.begin(), diag_int.end());
        RealType diag_max = *std::max_element(diag_int.begin(), diag_int.end());
        int num_zero_diag = std::count(diag_int.begin(), diag_int.end(), RealType(0));
        std::cout << "   Interior diagonal: range [" << diag_min << ", " << diag_max 
                  << "], zeros: " << num_zero_diag << "\n";
        
        // Allocate and populate interior system matrix
        A_int.allocate(n_int, n_int, nnz_int);
        
        // Copy data directly to the matrix vectors
        thrust::copy(h_rowOffsets_int.begin(), h_rowOffsets_int.end(),
                    thrust::device_pointer_cast(A_int.rowOffsetsPtr()));
        thrust::copy(h_colIndices_int.begin(), h_colIndices_int.end(),
                    thrust::device_pointer_cast(A_int.colIndicesPtr()));
        thrust::copy(h_values_int.begin(), h_values_int.end(),
                    thrust::device_pointer_cast(A_int.valuesPtr()));
        
        // Copy RHS to device
        b_int.resize(n_int);
        thrust::copy(h_b_int.begin(), h_b_int.end(),
                    thrust::device_pointer_cast(b_int.data()));
        
        // Store mapping for reconstruction
        interiorToGlobal_ = interiorToGlobal;
    }
    
    // Reconstruct full solution from interior solution
    void reconstructFullSolution(const Vector& u_int,
                                const std::vector<bool>& isBoundaryDOF,
                                const std::vector<RealType>& boundaryValues,
                                Vector& u_full)
    {
        size_t n = isBoundaryDOF.size();
        
        // Handle trivial case: all DOFs are boundary
        if (interiorToGlobal_.size() == 0) {
            u_full.resize(n);
            std::vector<RealType> h_u_full(n);
            for (size_t i = 0; i < n; ++i) {
                h_u_full[i] = boundaryValues[i];
            }
            cudaMemcpy(u_full.data(), h_u_full.data(), 
                      n * sizeof(RealType), cudaMemcpyHostToDevice);
            return;
        }
        
        u_full.resize(n);
        std::vector<RealType> h_u_full(n, 0.0);
        
        // Copy interior solution to host
        std::vector<RealType> h_u_int(u_int.size());
        cudaMemcpy(h_u_int.data(), u_int.data(), 
                  u_int.size() * sizeof(RealType), cudaMemcpyDeviceToHost);
        
        for (size_t i_int = 0; i_int < interiorToGlobal_.size(); ++i_int) {
            IndexType i_global = interiorToGlobal_[i_int];
            h_u_full[i_global] = h_u_int[i_int];
        }
        
        // Set boundary values
        for (size_t i = 0; i < n; ++i) {
            if (isBoundaryDOF[i]) {
                h_u_full[i] = boundaryValues[i];
            }
        }
        
        // Copy back to device
        cudaMemcpy(u_full.data(), h_u_full.data(), 
                  n * sizeof(RealType), cudaMemcpyHostToDevice);
    }

private:
    std::vector<IndexType> interiorToGlobal_;
};

} // namespace fem
} // namespace mars
