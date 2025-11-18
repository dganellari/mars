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

// CUDA kernel to modify RHS for columns corresponding to boundary DOFs
template<typename RealType, typename IndexType>
__global__ void modifyRHSForBCKernel(const IndexType* boundaryDofs,
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
    // (Could use thrust::binary_search but keeping it simple)
    bool isBoundaryRow = false;
    for (IndexType i = 0; i < numBoundaryDofs; ++i)
    {
        if (row == boundaryDofs[i])
        {
            isBoundaryRow = true;
            break;
        }
    }

    if (isBoundaryRow) return; // Skip boundary rows

    // For each column that's a boundary DOF, subtract contribution
    IndexType rowStart = rowOffsets[row];
    IndexType rowEnd   = rowOffsets[row + 1];

    for (IndexType i = rowStart; i < rowEnd; ++i)
    {
        IndexType col = colIndices[i];

        // Check if column is a boundary DOF
        for (IndexType j = 0; j < numBoundaryDofs; ++j)
        {
            if (col == boundaryDofs[j])
            {
                rhs[row] -= values[i] * bcValue;
                break;
            }
        }
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

        // Step 1: Modify RHS for boundary column contributions
        const int blockSize = 256;
        int gridSize        = (numRows + blockSize - 1) / blockSize;

        modifyRHSForBCKernel<RealType, KeyType><<<gridSize, blockSize>>>(
            thrust::raw_pointer_cast(d_boundaryDofs.data()), numBoundaryDofs, bcValue, A.rowOffsetsPtr(),
            A.colIndicesPtr(), A.valuesPtr(), thrust::raw_pointer_cast(b.data()), numRows);

        cudaDeviceSynchronize();

        // Step 2: Set boundary rows to identity
        gridSize = (numBoundaryDofs + blockSize - 1) / blockSize;

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
