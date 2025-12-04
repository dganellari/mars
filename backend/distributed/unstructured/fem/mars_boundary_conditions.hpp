#pragma once

#include "mars_sparse_matrix.hpp"
#include "mars_h1_fe_space.hpp"
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/binary_search.h>

namespace mars
{
namespace fem
{

// Penalty method kernel - adds scaled penalty to diagonal
template<typename RealType, typename IndexType>
__global__ void applyPenaltyBCKernel(const IndexType* boundaryDofs,
                                     IndexType numBoundaryDofs,
                                     RealType bcValue,
                                     RealType penalty,
                                     const IndexType* rowOffsets,
                                     const IndexType* colIndices,
                                     RealType* values,
                                     RealType* rhs,
                                     IndexType numRows)
{
    IndexType idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numBoundaryDofs) return;

    IndexType row = boundaryDofs[idx];

    // Find diagonal entry
    IndexType rowStart = rowOffsets[row];
    IndexType rowEnd   = rowOffsets[row + 1];

    RealType originalDiag = 0.0;
    IndexType diagIdx = -1;
    
    for (IndexType i = rowStart; i < rowEnd; ++i)
    {
        if (colIndices[i] == row)
        {
            originalDiag = values[i];
            diagIdx = i;
            break;
        }
    }
    
    if (diagIdx >= 0) {
        // Use scaled penalty: multiply diagonal by large factor
        values[diagIdx] = originalDiag * penalty;
        
        // Scale RHS accordingly: b[row] = originalDiag * penalty * bcValue
        rhs[row] = originalDiag * penalty * bcValue;
    }
}

// CUDA kernel to apply Dirichlet boundary conditions
// Sets matrix rows for boundary DOFs to identity, modifies RHS
template<typename RealType, typename IndexType>
__global__ void applyDirichletKernel(const IndexType* boundaryDofs,
                                     IndexType numBoundaryDofs,
                                     RealType bcValue,
                                     const IndexType* rowOffsets,
                                     const IndexType* colIndices,
                                     RealType* values,
                                     RealType* rhs,
                                     IndexType numRows)
{
    IndexType idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numBoundaryDofs) return;

    IndexType row = boundaryDofs[idx];

    // Set RHS to boundary value
    rhs[row] = bcValue;

    // Zero out row and set diagonal to 1
    IndexType rowStart = rowOffsets[row];
    IndexType rowEnd   = rowOffsets[row + 1];

    for (IndexType i = rowStart; i < rowEnd; ++i)
    {
        if (colIndices[i] == row)
        {
            values[i] = 1.0; // Diagonal
        }
        else
        {
            values[i] = 0.0; // Off-diagonal
        }
    }
}

// Kernel 1: Modify RHS for boundary column contributions (before zeroing)
template<typename RealType, typename IndexType>
__global__ void eliminateBCStep1Kernel(const IndexType* boundaryDofs,
                                       IndexType numBoundaryDofs,
                                       RealType bcValue,
                                       const IndexType* rowOffsets,
                                       const IndexType* colIndices,
                                       const RealType* values,
                                       RealType* rhs,
                                       IndexType numRows)
{
    IndexType row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= numRows) return;

    // Check if this row is a boundary DOF
    IndexType left = 0, right = numBoundaryDofs - 1;
    bool isBoundaryRow = false;
    while (left <= right)
    {
        IndexType mid = (left + right) / 2;
        if (boundaryDofs[mid] == row)
        {
            isBoundaryRow = true;
            break;
        }
        else if (boundaryDofs[mid] < row)
            left = mid + 1;
        else
            right = mid - 1;
    }

    if (!isBoundaryRow)
    {
        // Interior row: subtract contributions from boundary columns
        IndexType rowStart = rowOffsets[row];
        IndexType rowEnd = rowOffsets[row + 1];

        for (IndexType i = rowStart; i < rowEnd; ++i)
        {
            IndexType col = colIndices[i];

            // Check if column is a boundary DOF
            left = 0;
            right = numBoundaryDofs - 1;
            while (left <= right)
            {
                IndexType mid = (left + right) / 2;
                if (boundaryDofs[mid] == col)
                {
                    // Subtract contribution: b[row] -= A[row,col] * bcValue
                    rhs[row] -= values[i] * bcValue;
                    break;
                }
                else if (boundaryDofs[mid] < col)
                    left = mid + 1;
                else
                    right = mid - 1;
            }
        }
    }
}

// Kernel 2: Zero boundary rows/columns and set identity
template<typename RealType, typename IndexType>
__global__ void eliminateBCStep2Kernel(const IndexType* boundaryDofs,
                                       IndexType numBoundaryDofs,
                                       RealType bcValue,
                                       const IndexType* rowOffsets,
                                       const IndexType* colIndices,
                                       RealType* values,
                                       RealType* rhs,
                                       IndexType numRows)
{
    IndexType row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= numRows) return;

    IndexType rowStart = rowOffsets[row];
    IndexType rowEnd = rowOffsets[row + 1];

    // Check if this row is a boundary DOF
    IndexType left = 0, right = numBoundaryDofs - 1;
    bool isBoundaryRow = false;
    while (left <= right)
    {
        IndexType mid = (left + right) / 2;
        if (boundaryDofs[mid] == row)
        {
            isBoundaryRow = true;
            break;
        }
        else if (boundaryDofs[mid] < row)
            left = mid + 1;
        else
            right = mid - 1;
    }

    if (isBoundaryRow)
    {
        // Boundary row: set to identity, RHS to bcValue
        for (IndexType i = rowStart; i < rowEnd; ++i)
        {
            IndexType col = colIndices[i];
            values[i] = (col == row) ? 1.0 : 0.0;
        }
        rhs[row] = bcValue;
    }
    else
    {
        // Interior row: zero columns corresponding to boundary DOFs
        for (IndexType i = rowStart; i < rowEnd; ++i)
        {
            IndexType col = colIndices[i];

            // Check if column is a boundary DOF
            left = 0;
            right = numBoundaryDofs - 1;
            while (left <= right)
            {
                IndexType mid = (left + right) / 2;
                if (boundaryDofs[mid] == col)
                {
                    values[i] = 0.0;
                    break;
                }
                else if (boundaryDofs[mid] < col)
                    left = mid + 1;
                else
                    right = mid - 1;
            }
        }
    }
}

// Old unused kernel - keeping for reference
// Simple kernel to apply Dirichlet BCs - SYMMETRIC penalty method
template<typename RealType, typename IndexType>
__global__ void applySymmetricDirichletKernel_OLD(const IndexType* boundaryDofs,
                                               IndexType numBoundaryDofs,
                                               RealType bcValue,
                                               const IndexType* rowOffsets,
                                               const IndexType* colIndices,
                                               RealType* values,
                                               RealType* rhs,
                                               IndexType numRows)
{
    IndexType row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= numRows) return;

    IndexType rowStart = rowOffsets[row];
    IndexType rowEnd   = rowOffsets[row + 1];

    // Check if this row is a boundary DOF (binary search on sorted array)
    IndexType left = 0, right = numBoundaryDofs - 1;
    bool isBoundaryRow = false;
    while (left <= right)
    {
        IndexType mid = (left + right) / 2;
        if (boundaryDofs[mid] == row)
        {
            isBoundaryRow = true;
            break;
        }
        else if (boundaryDofs[mid] < row)
            left = mid + 1;
        else
            right = mid - 1;
    }

    if (isBoundaryRow)
    {
        // Penalty method: Set diagonal to large value, RHS to large*bcValue
        // This maintains symmetry and positive definiteness
        // Use penalty = max_diagonal * 1e7 for good conditioning
        RealType penalty = 1e7; // Moderate penalty value
        
        for (IndexType i = rowStart; i < rowEnd; ++i)
        {
            IndexType col = colIndices[i];
            if (col == row)
            {
                values[i] += penalty; // Add penalty to existing diagonal
            }
            // Don't zero off-diagonals - keep original stiffness
        }
        rhs[row] += penalty * bcValue;
    }
    else
    {
        // Interior rows: no modification needed with penalty method
        // The penalty on boundary DOFs naturally enforces the BC
    }
}

// Boundary condition handler class
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
class BoundaryConditionHandler
{
public:
    using FESpace         = H1FESpace<ElementTag, RealType, KeyType, AcceleratorTag>;
    using Matrix          = SparseMatrix<KeyType, RealType, AcceleratorTag>;
    using Vector          = typename mars::VectorSelector<RealType, AcceleratorTag>::type;
    using DeviceVectorKey = typename mars::VectorSelector<KeyType, AcceleratorTag>::type;

    // Apply homogeneous Dirichlet boundary conditions (u = value on boundary)
    void
    applyDirichlet(FESpace& fes, Matrix& A, Vector& b, const std::vector<KeyType>& boundaryDofs, RealType bcValue = 0.0)
    {
        if (boundaryDofs.empty()) return;

        // Copy boundary DOFs to device
        DeviceVectorKey d_boundaryDofs(boundaryDofs.size());
        thrust::copy(boundaryDofs.begin(), boundaryDofs.end(), thrust::device_pointer_cast(d_boundaryDofs.data()));

        // Sort boundary DOFs for efficient searching
        thrust::sort(thrust::device_pointer_cast(d_boundaryDofs.data()),
                     thrust::device_pointer_cast(d_boundaryDofs.data() + d_boundaryDofs.size()));

        // Remove duplicates
        auto newEnd = thrust::unique(thrust::device_pointer_cast(d_boundaryDofs.data()),
                                     thrust::device_pointer_cast(d_boundaryDofs.data() + d_boundaryDofs.size()));
        d_boundaryDofs.resize(newEnd - thrust::device_pointer_cast(d_boundaryDofs.data()));

        size_t numBoundaryDofs = d_boundaryDofs.size();
        size_t numRows         = A.numRows();

        const int blockSize = 256;
        int gridSize = (numBoundaryDofs + blockSize - 1) / blockSize;

        // Simple row replacement for now
        // TODO: Implement proper DOF elimination for better conditioning
        applyDirichletKernel<RealType, KeyType><<<gridSize, blockSize>>>(
            thrust::raw_pointer_cast(d_boundaryDofs.data()), numBoundaryDofs, bcValue, A.rowOffsetsPtr(),
            A.colIndicesPtr(), A.valuesPtr(), thrust::raw_pointer_cast(b.data()), numRows);

        cudaDeviceSynchronize();

        // Check for errors
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess)
        {
            throw std::runtime_error(std::string("CUDA error in BC application: ") + cudaGetErrorString(err));
        }
    }

    // Apply inhomogeneous Dirichlet BC with function
    template<typename BCFunc>
    void applyDirichletFunction(
        FESpace& fes, Matrix& A, Vector& b, const std::vector<KeyType>& boundaryDofs, BCFunc bcFunction)
    {
        // For inhomogeneous BC, need to evaluate function at boundary nodes
        // This is more complex - for now, implement constant BC version above
        // TODO: Implement with coordinate lookups
        throw std::runtime_error("Inhomogeneous Dirichlet BC not yet implemented");
    }
};

} // namespace fem
} // namespace mars
