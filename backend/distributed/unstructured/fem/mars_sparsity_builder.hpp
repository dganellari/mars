#pragma once

#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/execution_policy.h>
#include <cub/cub.cuh>
#include <cstdlib>
#include <string>
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
    //
    // At runtime, MARS_SPARSITY_CUB=1 selects the CUB-backed implementation
    // (radix sort on packed (row,col) keys; faster on large nnz). Default is
    // the original thrust-based path until CUB is validated bit-exact.
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
        int* d_diagPtr = nullptr,
        cudaStream_t stream = 0);

    // CUB-backed implementation. Same I/O contract as buildGraphSparsity.
    // Internal differences:
    //   - Packs (row, col) into uint64_t (row in high 32, col in low 32)
    //   - cub::DeviceRadixSort::SortKeys on uint64_t (bulk-typed, faster than
    //     thrust comparator sort over zip<int,int>)
    //   - cub::DeviceSelect::Unique to dedupe
    //   - Histogram via thrust::for_each + atomicAdd (fastest in practice;
    //     CUB DeviceHistogram has worse perf for sparse-row distributions)
    //   - cub::DeviceScan::ExclusiveSum for rowPtr
    static int buildGraphSparsityCub(
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
        int* d_diagPtr = nullptr,
        cudaStream_t stream = 0);

    // Full sparsity pattern (27 NNZ/row): every (i,j) pair within each
    // hex8's 8 nodes contributes a NNZ entry. Required when assembling
    // with assembleFull (full element-local 8x8 stiffness).
    static int buildFullSparsity(
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
        int* d_diagPtr = nullptr,
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

// Runtime dispatch: env MARS_SPARSITY_CUB=1 selects CUB path.
namespace mars_sparsity_detail {
inline bool useCubSparsity()
{
    static const bool on = []() {
        const char* e = std::getenv("MARS_SPARSITY_CUB");
        return e != nullptr && std::string(e) != "0";
    }();
    return on;
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
    if (mars_sparsity_detail::useCubSparsity()) {
        return buildGraphSparsityCub(
            d_conn0, d_conn1, d_conn2, d_conn3,
            d_conn4, d_conn5, d_conn6, d_conn7,
            numElements, d_nodeToDof, numDofs,
            d_rowPtr, d_colInd, d_diagPtr, stream);
    }

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
    // For each row r, find position k where colInd[k] == r
    if (d_diagPtr != nullptr && d_colInd != nullptr) {
        // Launch kernel to find diagonal positions
        thrust::for_each(
            thrust::device,
            thrust::make_counting_iterator(0),
            thrust::make_counting_iterator(numDofs),
            [rowPtr = d_rowPtr, colInd = d_colInd, diagPtr = d_diagPtr] __device__ (int row) {
                int start = rowPtr[row];
                int end = rowPtr[row + 1];
                // Binary search for diagonal (columns are sorted)
                int lo = start, hi = end;
                while (lo < hi) {
                    int mid = (lo + hi) / 2;
                    if (colInd[mid] < row) {
                        lo = mid + 1;
                    } else {
                        hi = mid;
                    }
                }
                // lo should now point to the diagonal
                diagPtr[row] = lo;
            }
        );
    }

    return nnz;
}

// =====================================================================
// CUB-backed buildGraphSparsityCub
// =====================================================================
// Same I/O contract as buildGraphSparsity. Internal pipeline:
//   1. Build edge list with the existing kernel.
//   2. Append diagonals (sequence).
//   3. Pack (row, col) → uint64_t key (row in high 32, col in low 32).
//      Invalid entries (row<0 or col<0) packed as MAX so they sort to the end
//      and get clipped before unique.
//   4. cub::DeviceRadixSort::SortKeys over uint64_t (single-pass radix; faster
//      than thrust comparator sort over zip<int,int>).
//   5. cub::DeviceSelect::Unique to dedupe.
//   6. Decode keys back to row[], col[]. Histogram rows via atomicAdd.
//   7. cub::DeviceScan::ExclusiveSum for rowPtr.
//   8. Diagonal lookup: same as thrust path.
//
// Expected speedup at 217M nonzeros (cube256/16 L1 case): the sort step
// drops from ~150ms (thrust) to ~30-50ms (CUB radix). End-to-end win
// roughly 2-3x on this phase.
template<typename KeyType>
int CvfemSparsityBuilder<KeyType>::buildGraphSparsityCub(
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
    using PackedKey = unsigned long long;
    constexpr PackedKey INVALID_PACKED = ~PackedKey(0);

    // Step 1: edge list as before (12 SCS × 2 directions = 24 edges/element)
    size_t numEdges = numElements * 12 * 2;
    thrust::device_vector<int> d_edgeRow(numEdges);
    thrust::device_vector<int> d_edgeCol(numEdges);

    int blockSize = 256;
    int numBlocks = (numElements + blockSize - 1) / blockSize;
    buildSparsityEdgeListKernel<<<numBlocks, blockSize, 0, stream>>>(
        d_conn0, d_conn1, d_conn2, d_conn3,
        d_conn4, d_conn5, d_conn6, d_conn7,
        d_nodeToDof, numElements,
        thrust::raw_pointer_cast(d_edgeRow.data()),
        thrust::raw_pointer_cast(d_edgeCol.data()));

    // Step 2-3: pack edges + diagonals into uint64_t keys.
    size_t totalEntries = numEdges + numDofs;
    thrust::device_vector<PackedKey> d_keysIn(totalEntries);

    // Pack edges. Threads with row<0 || col<0 emit INVALID_PACKED.
    {
        int nB = (numEdges + blockSize - 1) / blockSize;
        const int* rowPtr_ = thrust::raw_pointer_cast(d_edgeRow.data());
        const int* colPtr_ = thrust::raw_pointer_cast(d_edgeCol.data());
        PackedKey* outPtr  = thrust::raw_pointer_cast(d_keysIn.data());
        thrust::for_each(thrust::cuda::par.on(stream),
            thrust::counting_iterator<size_t>(0),
            thrust::counting_iterator<size_t>(numEdges),
            [rowPtr_, colPtr_, outPtr] __device__ (size_t i) {
                int r = rowPtr_[i];
                int c = colPtr_[i];
                if (r < 0 || c < 0) {
                    outPtr[i] = INVALID_PACKED;
                } else {
                    outPtr[i] = (static_cast<PackedKey>(static_cast<unsigned int>(r)) << 32)
                              |  static_cast<PackedKey>(static_cast<unsigned int>(c));
                }
            });
    }
    // Pack diagonals starting at offset numEdges. Diagonal (i,i) for i in [0, numDofs).
    {
        PackedKey* outPtr = thrust::raw_pointer_cast(d_keysIn.data()) + numEdges;
        thrust::for_each(thrust::cuda::par.on(stream),
            thrust::counting_iterator<int>(0),
            thrust::counting_iterator<int>(numDofs),
            [outPtr] __device__ (int i) {
                outPtr[i] = (static_cast<PackedKey>(static_cast<unsigned int>(i)) << 32)
                          |  static_cast<PackedKey>(static_cast<unsigned int>(i));
            });
    }

    // Free the int arrays now; CUB sort is on packed keys only.
    d_edgeRow.clear(); d_edgeRow.shrink_to_fit();
    d_edgeCol.clear(); d_edgeCol.shrink_to_fit();

    // Step 4: CUB radix sort.
    thrust::device_vector<PackedKey> d_keysOut(totalEntries);
    {
        size_t tempBytes = 0;
        cub::DeviceRadixSort::SortKeys(
            nullptr, tempBytes,
            thrust::raw_pointer_cast(d_keysIn.data()),
            thrust::raw_pointer_cast(d_keysOut.data()),
            int(totalEntries), 0, sizeof(PackedKey) * 8, stream);
        thrust::device_vector<unsigned char> d_temp(tempBytes);
        cub::DeviceRadixSort::SortKeys(
            thrust::raw_pointer_cast(d_temp.data()), tempBytes,
            thrust::raw_pointer_cast(d_keysIn.data()),
            thrust::raw_pointer_cast(d_keysOut.data()),
            int(totalEntries), 0, sizeof(PackedKey) * 8, stream);
    }

    // Free the input buffer; we keep d_keysOut.
    d_keysIn.clear(); d_keysIn.shrink_to_fit();

    // Step 5: CUB unique. Output count via DeviceSelect::Unique.
    thrust::device_vector<PackedKey> d_keysUniq(totalEntries);
    int  h_uniqCount = 0;
    {
        thrust::device_vector<int> d_uniqCount(1);
        size_t tempBytes = 0;
        cub::DeviceSelect::Unique(
            nullptr, tempBytes,
            thrust::raw_pointer_cast(d_keysOut.data()),
            thrust::raw_pointer_cast(d_keysUniq.data()),
            thrust::raw_pointer_cast(d_uniqCount.data()),
            int(totalEntries), stream);
        thrust::device_vector<unsigned char> d_temp(tempBytes);
        cub::DeviceSelect::Unique(
            thrust::raw_pointer_cast(d_temp.data()), tempBytes,
            thrust::raw_pointer_cast(d_keysOut.data()),
            thrust::raw_pointer_cast(d_keysUniq.data()),
            thrust::raw_pointer_cast(d_uniqCount.data()),
            int(totalEntries), stream);
        cudaMemcpyAsync(&h_uniqCount, thrust::raw_pointer_cast(d_uniqCount.data()),
                        sizeof(int), cudaMemcpyDeviceToHost, stream);
        cudaStreamSynchronize(stream);
    }

    // Strip trailing INVALID_PACKED entries (sort puts them at the end).
    // Binary-search for the first INVALID. Cheap one-shot device call.
    int nnz = h_uniqCount;
    {
        // Pull the last entry; if it's INVALID, scan back. In practice INVALID
        // is a contiguous tail block; one binary search on host suffices via a
        // single d2h of the last few entries. Simplest: pull last entry, if
        // INVALID, do a host-side binary search via a small d2h.
        if (h_uniqCount > 0) {
            PackedKey lastKey = 0;
            cudaMemcpyAsync(&lastKey,
                thrust::raw_pointer_cast(d_keysUniq.data()) + h_uniqCount - 1,
                sizeof(PackedKey), cudaMemcpyDeviceToHost, stream);
            cudaStreamSynchronize(stream);
            if (lastKey == INVALID_PACKED) {
                // Linear scan from the end (worst case ~few k entries are INVALID).
                // Using thrust::lower_bound to keep all on device:
                auto it = thrust::lower_bound(
                    thrust::cuda::par.on(stream),
                    d_keysUniq.begin(), d_keysUniq.begin() + h_uniqCount,
                    INVALID_PACKED);
                nnz = int(it - d_keysUniq.begin());
            }
        }
    }
    d_keysOut.clear(); d_keysOut.shrink_to_fit();

    // Step 6: decode keys back → row[], col[]. Histogram rows.
    thrust::device_vector<int> d_rowCounts(numDofs, 0);
    thrust::device_vector<int> d_decodedCol(nnz);
    {
        const PackedKey* keysPtr = thrust::raw_pointer_cast(d_keysUniq.data());
        int* colOutPtr = thrust::raw_pointer_cast(d_decodedCol.data());
        int* rowCntPtr = thrust::raw_pointer_cast(d_rowCounts.data());
        thrust::for_each(thrust::cuda::par.on(stream),
            thrust::counting_iterator<int>(0),
            thrust::counting_iterator<int>(nnz),
            [keysPtr, colOutPtr, rowCntPtr] __device__ (int i) {
                PackedKey k = keysPtr[i];
                int row = static_cast<int>(static_cast<unsigned int>(k >> 32));
                int col = static_cast<int>(static_cast<unsigned int>(k & 0xFFFFFFFFull));
                colOutPtr[i] = col;
                atomicAdd(&rowCntPtr[row], 1);
            });
    }

    // Step 7: rowPtr via cub::DeviceScan::ExclusiveSum.
    {
        size_t tempBytes = 0;
        cub::DeviceScan::ExclusiveSum(
            nullptr, tempBytes,
            thrust::raw_pointer_cast(d_rowCounts.data()),
            d_rowPtr,
            numDofs + 1, stream);
        thrust::device_vector<unsigned char> d_temp(tempBytes);
        cub::DeviceScan::ExclusiveSum(
            thrust::raw_pointer_cast(d_temp.data()), tempBytes,
            thrust::raw_pointer_cast(d_rowCounts.data()),
            d_rowPtr,
            numDofs + 1, stream);
        // CUB ExclusiveSum writes numDofs+1 entries reading numDofs+1 inputs.
        // Last input is undefined (we wrote only numDofs counts); the scan's
        // last output equals sum of first numDofs counts = nnz. To be safe,
        // explicitly write rowPtr[numDofs] = nnz.
        cudaMemcpyAsync(d_rowPtr + numDofs, &nnz, sizeof(int),
                        cudaMemcpyHostToDevice, stream);
    }

    // Copy decoded col indices to caller's d_colInd (if requested).
    if (d_colInd != nullptr) {
        cudaMemcpyAsync(d_colInd,
                        thrust::raw_pointer_cast(d_decodedCol.data()),
                        nnz * sizeof(int),
                        cudaMemcpyDeviceToDevice, stream);
    }

    // Diagonal positions: same logic as thrust path.
    if (d_diagPtr != nullptr && d_colInd != nullptr) {
        thrust::for_each(thrust::cuda::par.on(stream),
            thrust::make_counting_iterator(0),
            thrust::make_counting_iterator(numDofs),
            [rowPtr = d_rowPtr, colInd = d_colInd, diagPtr = d_diagPtr]
            __device__ (int row) {
                int start = rowPtr[row];
                int end   = rowPtr[row + 1];
                int lo = start, hi = end;
                while (lo < hi) {
                    int mid = (lo + hi) / 2;
                    if (colInd[mid] < row) lo = mid + 1;
                    else                   hi = mid;
                }
                diagPtr[row] = lo;
            });
    }
    cudaStreamSynchronize(stream);
    return nnz;
}

// Full element-local sparsity: emit all 8x8 = 64 (dofI, dofJ) pairs per hex.
template<typename KeyType>
__global__ void buildFullEdgeListKernel(
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
    nodes[0] = conn0[e]; nodes[1] = conn1[e]; nodes[2] = conn2[e]; nodes[3] = conn3[e];
    nodes[4] = conn4[e]; nodes[5] = conn5[e]; nodes[6] = conn6[e]; nodes[7] = conn7[e];

    int dofs[8];
    for (int n = 0; n < 8; ++n) dofs[n] = nodeToDof[nodes[n]];

    size_t base = e * 64;
    for (int i = 0; i < 8; ++i)
    {
        for (int j = 0; j < 8; ++j)
        {
            size_t idx = base + i * 8 + j;
            if (dofs[i] >= 0 && dofs[j] >= 0)
            {
                edgeListRow[idx] = dofs[i];
                edgeListCol[idx] = dofs[j];
            }
            else
            {
                edgeListRow[idx] = -1;
                edgeListCol[idx] = -1;
            }
        }
    }
}

template<typename KeyType>
int CvfemSparsityBuilder<KeyType>::buildFullSparsity(
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
    size_t numEdges = numElements * 64;
    thrust::device_vector<int> d_edgeRow(numEdges);
    thrust::device_vector<int> d_edgeCol(numEdges);

    int blockSize = 256;
    int numBlocks = (numElements + blockSize - 1) / blockSize;
    buildFullEdgeListKernel<<<numBlocks, blockSize, 0, stream>>>(
        d_conn0, d_conn1, d_conn2, d_conn3,
        d_conn4, d_conn5, d_conn6, d_conn7,
        d_nodeToDof, numElements,
        thrust::raw_pointer_cast(d_edgeRow.data()),
        thrust::raw_pointer_cast(d_edgeCol.data()));

    auto new_end = thrust::remove_if(
        thrust::device,
        thrust::make_zip_iterator(thrust::make_tuple(d_edgeRow.begin(), d_edgeCol.begin())),
        thrust::make_zip_iterator(thrust::make_tuple(d_edgeRow.end(), d_edgeCol.end())),
        [] __device__(const thrust::tuple<int, int>& t) { return thrust::get<0>(t) < 0; });

    size_t validEntries = thrust::distance(
        thrust::make_zip_iterator(thrust::make_tuple(d_edgeRow.begin(), d_edgeCol.begin())), new_end);
    d_edgeRow.resize(validEntries);
    d_edgeCol.resize(validEntries);

    thrust::sort(
        thrust::device,
        thrust::make_zip_iterator(thrust::make_tuple(d_edgeRow.begin(), d_edgeCol.begin())),
        thrust::make_zip_iterator(thrust::make_tuple(d_edgeRow.end(), d_edgeCol.end())));

    auto unique_end = thrust::unique(
        thrust::device,
        thrust::make_zip_iterator(thrust::make_tuple(d_edgeRow.begin(), d_edgeCol.begin())),
        thrust::make_zip_iterator(thrust::make_tuple(d_edgeRow.end(), d_edgeCol.end())));

    int nnz = thrust::distance(
        thrust::make_zip_iterator(thrust::make_tuple(d_edgeRow.begin(), d_edgeCol.begin())), unique_end);
    d_edgeRow.resize(nnz);
    d_edgeCol.resize(nnz);

    thrust::device_vector<int> d_rowCounts(numDofs, 0);
    thrust::for_each(
        thrust::device, d_edgeRow.begin(), d_edgeRow.end(),
        [counts = thrust::raw_pointer_cast(d_rowCounts.data())] __device__(int row) { atomicAdd(&counts[row], 1); });

    thrust::device_vector<int> d_rowPtrTemp(numDofs + 1);
    thrust::exclusive_scan(thrust::device, d_rowCounts.begin(), d_rowCounts.end(), d_rowPtrTemp.begin());
    int lastCount, lastPtr;
    cudaMemcpy(&lastCount, thrust::raw_pointer_cast(d_rowCounts.data() + numDofs - 1), sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&lastPtr, thrust::raw_pointer_cast(d_rowPtrTemp.data() + numDofs - 1), sizeof(int), cudaMemcpyDeviceToHost);
    int finalVal = lastPtr + lastCount;
    cudaMemcpy(thrust::raw_pointer_cast(d_rowPtrTemp.data() + numDofs), &finalVal, sizeof(int), cudaMemcpyHostToDevice);

    cudaMemcpy(d_rowPtr, thrust::raw_pointer_cast(d_rowPtrTemp.data()), (numDofs + 1) * sizeof(int),
                cudaMemcpyDeviceToDevice);

    if (d_colInd != nullptr)
    {
        cudaMemcpy(d_colInd, thrust::raw_pointer_cast(d_edgeCol.data()), nnz * sizeof(int), cudaMemcpyDeviceToDevice);
    }

    if (d_diagPtr != nullptr && d_colInd != nullptr)
    {
        thrust::for_each(
            thrust::device, thrust::make_counting_iterator(0), thrust::make_counting_iterator(numDofs),
            [rowPtr = d_rowPtr, colInd = d_colInd, diagPtr = d_diagPtr] __device__(int row)
            {
                int start = rowPtr[row];
                int end   = rowPtr[row + 1];
                int lo = start, hi = end;
                while (lo < hi)
                {
                    int mid = lo + (hi - lo) / 2;
                    if (colInd[mid] < row) lo = mid + 1;
                    else hi = mid;
                }
                diagPtr[row] = lo;
            });
    }

    return nnz;
}

} // namespace fem
} // namespace mars
