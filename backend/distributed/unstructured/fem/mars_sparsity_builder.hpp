#pragma once

#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/execution_policy.h>
#include "mars_cvfem_utils.hpp"

namespace mars {
namespace fem {

// Forward declare the kernel (defined below as free function)
template<typename KeyType>
__global__ void buildSparsityEdgeListKernel(
    const KeyType* __restrict__ conn0,
    const KeyType* __restrict__ conn1,
    const KeyType* __restrict__ conn2,
    const KeyType* __restrict__ conn3,
    const KeyType* __restrict__ conn4,
    const KeyType* __restrict__ conn5,
    const KeyType* __restrict__ conn6,
    const KeyType* __restrict__ conn7,
    const int* __restrict__ nodeToDof,
    size_t numElements,
    int* edgeListRow,
    int* edgeListCol);

// GPU-based CSR sparsity pattern builder for CVFEM hex elements
template<typename KeyType>
class CvfemSparsityBuilder {
public:
    // Build reduced sparsity pattern (7 NNZ/row) from element connectivity
    // Returns number of non-zeros
    // If d_diagPtr is provided, also computes diagonal positions for each row
    static int buildGraphSparsity(
        const KeyType* d_conn0,
        const KeyType* d_conn1,
        const KeyType* d_conn2,
        const KeyType* d_conn3,
        const KeyType* d_conn4,
        const KeyType* d_conn5,
        const KeyType* d_conn6,
        const KeyType* d_conn7,
        size_t numElements,
        const int* d_nodeToDof,
        int numDofs,
        int* d_rowPtr,
        int* d_colInd,
        int* d_diagPtr = nullptr,  // Output: diagonal position for each row
        cudaStream_t stream = 0);
};

// GPU kernel implementation (free function)
template<typename KeyType>
__global__ void buildSparsityEdgeListKernel(
    const KeyType* __restrict__ conn0,
    const KeyType* __restrict__ conn1,
    const KeyType* __restrict__ conn2,
    const KeyType* __restrict__ conn3,
    const KeyType* __restrict__ conn4,
    const KeyType* __restrict__ conn5,
    const KeyType* __restrict__ conn6,
    const KeyType* __restrict__ conn7,
    const int* __restrict__ nodeToDof,
    size_t numElements,
    int* edgeListRow,
    int* edgeListCol)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements) return;

    KeyType nodes[8];
    nodes[0] = conn0[e];
    nodes[1] = conn1[e];
    nodes[2] = conn2[e];
    nodes[3] = conn3[e];
    nodes[4] = conn4[e];
    nodes[5] = conn5[e];
    nodes[6] = conn6[e];
    nodes[7] = conn7[e];

    // Each element contributes 12 SCS edges
    for (int scs = 0; scs < 12; ++scs) {
        int nodeL = d_hexLRSCV[scs * 2];
        int nodeR = d_hexLRSCV[scs * 2 + 1];

        int dofL = nodeToDof[nodes[nodeL]];
        int dofR = nodeToDof[nodes[nodeR]];

        size_t edgeIdx = e * 12 + scs;

        // Store both directions (L->R and R->L)
        if (dofL >= 0 && dofR >= 0) {
            edgeListRow[edgeIdx * 2] = dofL;
            edgeListCol[edgeIdx * 2] = dofR;
            edgeListRow[edgeIdx * 2 + 1] = dofR;
            edgeListCol[edgeIdx * 2 + 1] = dofL;
        } else {
            edgeListRow[edgeIdx * 2] = -1;
            edgeListCol[edgeIdx * 2] = -1;
            edgeListRow[edgeIdx * 2 + 1] = -1;
            edgeListCol[edgeIdx * 2 + 1] = -1;
        }
    }
}

template<typename KeyType>
int CvfemSparsityBuilder<KeyType>::buildGraphSparsity(
    const KeyType* d_conn0,
    const KeyType* d_conn1,
    const KeyType* d_conn2,
    const KeyType* d_conn3,
    const KeyType* d_conn4,
    const KeyType* d_conn5,
    const KeyType* d_conn6,
    const KeyType* d_conn7,
    size_t numElements,
    const int* d_nodeToDof,
    int numDofs,
    int* d_rowPtr,
    int* d_colInd,
    int* d_diagPtr,
    cudaStream_t stream)
{
    // Build edge list (each element contributes 24 edges: 12 SCS * 2 directions)
    size_t numEdges = numElements * 12 * 2;
    thrust::device_vector<int> d_edgeRow(numEdges);
    thrust::device_vector<int> d_edgeCol(numEdges);

    int blockSize = 256;
    int numBlocks = (numElements + blockSize - 1) / blockSize;
    buildSparsityEdgeListKernel<<<numBlocks, blockSize, 0, stream>>>(
        d_conn0, d_conn1, d_conn2, d_conn3,
        d_conn4, d_conn5, d_conn6, d_conn7,
        d_nodeToDof,
        numElements,
        thrust::raw_pointer_cast(d_edgeRow.data()),
        thrust::raw_pointer_cast(d_edgeCol.data())
    );

    // Add diagonal entries
    thrust::device_vector<int> d_diagRow(numDofs);
    thrust::device_vector<int> d_diagCol(numDofs);
    thrust::sequence(thrust::device, d_diagRow.begin(), d_diagRow.end());
    thrust::sequence(thrust::device, d_diagCol.begin(), d_diagCol.end());

    // Merge edges and diagonals
    thrust::device_vector<int> d_allRow(numEdges + numDofs);
    thrust::device_vector<int> d_allCol(numEdges + numDofs);
    thrust::copy(thrust::device, d_edgeRow.begin(), d_edgeRow.end(), d_allRow.begin());
    thrust::copy(thrust::device, d_diagRow.begin(), d_diagRow.end(), d_allRow.begin() + numEdges);
    thrust::copy(thrust::device, d_edgeCol.begin(), d_edgeCol.end(), d_allCol.begin());
    thrust::copy(thrust::device, d_diagCol.begin(), d_diagCol.end(), d_allCol.begin() + numEdges);

    // Remove invalid entries (-1)
    auto new_end = thrust::remove_if(
        thrust::device,
        thrust::make_zip_iterator(thrust::make_tuple(d_allRow.begin(), d_allCol.begin())),
        thrust::make_zip_iterator(thrust::make_tuple(d_allRow.end(), d_allCol.end())),
        [] __device__ (const thrust::tuple<int, int>& t) {
            return thrust::get<0>(t) < 0;
        }
    );

    size_t validEntries = thrust::distance(
        thrust::make_zip_iterator(thrust::make_tuple(d_allRow.begin(), d_allCol.begin())),
        new_end
    );
    d_allRow.resize(validEntries);
    d_allCol.resize(validEntries);

    // Sort by (row, col)
    thrust::sort(
        thrust::device,
        thrust::make_zip_iterator(thrust::make_tuple(d_allRow.begin(), d_allCol.begin())),
        thrust::make_zip_iterator(thrust::make_tuple(d_allRow.end(), d_allCol.end()))
    );

    // Remove duplicates
    auto unique_end = thrust::unique(
        thrust::device,
        thrust::make_zip_iterator(thrust::make_tuple(d_allRow.begin(), d_allCol.begin())),
        thrust::make_zip_iterator(thrust::make_tuple(d_allRow.end(), d_allCol.end()))
    );

    int nnz = thrust::distance(
        thrust::make_zip_iterator(thrust::make_tuple(d_allRow.begin(), d_allCol.begin())),
        unique_end
    );
    d_allRow.resize(nnz);
    d_allCol.resize(nnz);

    // Build CSR rowPtr using histogram approach
    // First count entries per row
    thrust::device_vector<int> d_rowCounts(numDofs, 0);
    thrust::for_each(
        thrust::device,
        d_allRow.begin(), d_allRow.end(),
        [counts = thrust::raw_pointer_cast(d_rowCounts.data())] __device__ (int row) {
            atomicAdd(&counts[row], 1);
        }
    );

    // Exclusive prefix sum to get row pointers
    thrust::device_vector<int> d_rowPtrTemp(numDofs + 1);
    thrust::exclusive_scan(thrust::device, d_rowCounts.begin(), d_rowCounts.end(), d_rowPtrTemp.begin());
    // Set last element
    int lastCount, lastPtr;
    cudaMemcpy(&lastCount, thrust::raw_pointer_cast(d_rowCounts.data() + numDofs - 1), sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&lastPtr, thrust::raw_pointer_cast(d_rowPtrTemp.data() + numDofs - 1), sizeof(int), cudaMemcpyDeviceToHost);
    int finalVal = lastPtr + lastCount;
    cudaMemcpy(thrust::raw_pointer_cast(d_rowPtrTemp.data() + numDofs), &finalVal, sizeof(int), cudaMemcpyHostToDevice);

    // Copy results
    cudaMemcpy(d_rowPtr, thrust::raw_pointer_cast(d_rowPtrTemp.data()),
               (numDofs + 1) * sizeof(int), cudaMemcpyDeviceToDevice);

    if (d_colInd != nullptr) {
        cudaMemcpy(d_colInd, thrust::raw_pointer_cast(d_allCol.data()),
                   nnz * sizeof(int), cudaMemcpyDeviceToDevice);
    }

    // Compute diagonal positions if requested
    // Since diagonal is always first entry in sorted row (col[row] == row),
    // diagPtr[row] = rowPtr[row]
    if (d_diagPtr != nullptr) {
        cudaMemcpy(d_diagPtr, d_rowPtr, numDofs * sizeof(int), cudaMemcpyDeviceToDevice);
    }

    return nnz;
}

} // namespace fem
} // namespace mars
