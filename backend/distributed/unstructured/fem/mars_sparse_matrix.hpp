#pragma once

#include <memory>
#include <vector>

namespace mars {
namespace fem {

// Forward declare kernel
template<typename RealType, typename IndexType>
__global__ void extractDiagonalKernel(const IndexType* rowOffsets,
                                      const IndexType* colIndices,
                                      const RealType* values,
                                      RealType* diagonal,
                                      IndexType numRows);

// CSR Sparse Matrix for GPU
template<typename IndexType, typename RealType, typename AcceleratorTag>
class SparseMatrix {
public:
    using DeviceVectorIndex = typename mars::VectorSelector<IndexType, AcceleratorTag>::type;
    using DeviceVectorReal = typename mars::VectorSelector<RealType, AcceleratorTag>::type;
    
    SparseMatrix() : numRows_(0), numCols_(0), nnz_(0) {}
    
    SparseMatrix(IndexType numRows, IndexType numCols)
        : numRows_(numRows), numCols_(numCols), nnz_(0) {
        // Initialize with empty CSR structure
        d_rowOffsets_.resize(numRows_ + 1, 0);
    }
    
    // Allocate matrix with known sparsity pattern
    void allocate(IndexType numRows, IndexType numCols, IndexType nnz) {
        numRows_ = numRows;
        numCols_ = numCols;
        nnz_ = nnz;
        
        d_rowOffsets_.resize(numRows_ + 1);
        d_colIndices_.resize(nnz_);
        d_values_.resize(nnz_);
    }
    
    // Build sparsity pattern from connectivity (before assembly)
    void buildSparsityPattern(const DeviceVectorIndex& nodeToElemOffsets,
                             const DeviceVectorIndex& nodeToElemList,
                             const DeviceVectorIndex& connectivity,
                             IndexType numNodesPerElem);
    
    // Zero out matrix values (keep sparsity pattern)
    void zero() {
        thrust::fill(thrust::device_pointer_cast(d_values_.data()),
                     thrust::device_pointer_cast(d_values_.data() + d_values_.size()),
                     RealType(0));
    }
    
    // Accessors
    IndexType numRows() const { return numRows_; }
    IndexType numCols() const { return numCols_; }
    IndexType nnz() const { return nnz_; }
    
    DeviceVectorIndex& rowOffsets() { return d_rowOffsets_; }
    const DeviceVectorIndex& rowOffsets() const { return d_rowOffsets_; }
    
    DeviceVectorReal& values() { return d_values_; }
    const DeviceVectorReal& values() const { return d_values_; }
    
    DeviceVectorIndex& colIndices() { return d_colIndices_; }
    const DeviceVectorIndex& colIndices() const { return d_colIndices_; }
    
    // Raw pointers for kernel access
    IndexType* rowOffsetsPtr() { return thrust::raw_pointer_cast(d_rowOffsets_.data()); }
    const IndexType* rowOffsetsPtr() const { return thrust::raw_pointer_cast(d_rowOffsets_.data()); }
    
    RealType* valuesPtr() { return thrust::raw_pointer_cast(d_values_.data()); }
    const RealType* valuesPtr() const { return thrust::raw_pointer_cast(d_values_.data()); }
    
    IndexType* colIndicesPtr() { return thrust::raw_pointer_cast(d_colIndices_.data()); }
    const IndexType* colIndicesPtr() const { return thrust::raw_pointer_cast(d_colIndices_.data()); }
    
    // Extract diagonal values
    DeviceVectorReal getDiagonal() const {
        DeviceVectorReal diag(numRows_);
        thrust::fill(thrust::device_pointer_cast(diag.data()),
                     thrust::device_pointer_cast(diag.data() + diag.size()),
                     RealType(0));
        
        // Kernel to extract diagonal
        const int blockSize = 256;
        const int gridSize = (numRows_ + blockSize - 1) / blockSize;
        
        extractDiagonalKernel<<<gridSize, blockSize>>>(
            thrust::raw_pointer_cast(d_rowOffsets_.data()),
            thrust::raw_pointer_cast(d_colIndices_.data()),
            thrust::raw_pointer_cast(d_values_.data()),
            thrust::raw_pointer_cast(diag.data()),
            numRows_);
        
        cudaDeviceSynchronize();
        return diag;
    }
    
private:
    IndexType numRows_;
    IndexType numCols_;
    IndexType nnz_;
    
    DeviceVectorIndex d_rowOffsets_;   // Size: numRows + 1
    DeviceVectorIndex d_colIndices_;   // Size: nnz
    DeviceVectorReal d_values_;        // Size: nnz
};

// CUDA kernel to extract diagonal
template<typename RealType, typename IndexType>
__global__ void extractDiagonalKernel(const IndexType* rowOffsets,
                                      const IndexType* colIndices,
                                      const RealType* values,
                                      RealType* diagonal,
                                      IndexType numRows)
{
    IndexType row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= numRows) return;
    
    IndexType rowStart = rowOffsets[row];
    IndexType rowEnd = rowOffsets[row + 1];
    
    for (IndexType i = rowStart; i < rowEnd; ++i) {
        if (colIndices[i] == row) {
            diagonal[row] = values[i];
            return;
        }
    }
}

// Device function to find column index in CSR row
template<typename IndexType>
__device__ inline int findColumnIndex(const IndexType* colIndices, 
                                     IndexType rowStart, 
                                     IndexType rowEnd, 
                                     IndexType targetCol) {
    for (IndexType i = rowStart; i < rowEnd; ++i) {
        if (colIndices[i] == targetCol) {
            return i;
        }
    }
    return -1;  // Not found
}

// Atomic add to sparse matrix entry
template<typename RealType, typename IndexType>
__device__ inline void atomicAddSparseEntry(RealType* values,
                                           const IndexType* colIndices,
                                           IndexType rowStart,
                                           IndexType rowEnd,
                                           IndexType col,
                                           RealType value) {
    int idx = findColumnIndex(colIndices, rowStart, rowEnd, col);
    if (idx >= 0) {
        atomicAdd(&values[idx], value);
    }
}

} // namespace fem
} // namespace mars
