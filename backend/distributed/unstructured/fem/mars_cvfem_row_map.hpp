#pragma once

#include "mars_cvfem_hex_kernel_tensor_rows.hpp"
#include "../utils/mars_domain_utils.hpp"
#include "cstone/cuda/device_vector.h"
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/scan.h>
#include <climits>
#include <cstdint>

namespace mars {
namespace fem {

// Which assembly a row map is for: a column outside the pattern goes to the diagonal (graph, lumped) or is
// dropped (full).
enum class CvfemRowMapPattern
{
    Graph,
    Full
};

// The row map of the tensor assemblies (mars_cvfem_row_map_format.hpp), built once on the GPU for a mesh, a
// DOF map and a sparsity pattern: 72 bytes per element (40 if compact), 6 per block and 6 per block row.
struct CvfemRowMap
{
    size_t numElements = 0;
    int numBlocks      = 0;
    int numBlockRows   = 0; // the rows of all blocks together
    int numDirectBlocks = 0; // blocks whose rows do not fit in shared memory
    int numSearchRows  = 0; // element rows longer than rowMapMaxLen
    bool lumped        = false;
    bool compact       = false; // 4-bit elemCols (rowMapCompactMaxLen)
    cstone::DeviceVector<uint8_t> elemRow;
    cstone::DeviceVector<uint8_t> elemCols;
    // uint16_t values, 2 bytes each: DeviceVector has no uint16_t instantiation
    cstone::DeviceVector<uint8_t> blockSize;
    cstone::DeviceVector<uint8_t> rowBase;
    cstone::DeviceVector<int> rowDof;
    // a thrust vector: scanned with the sparsity builder's own thrust call, so no new scan kernels compile
    thrust::device_vector<int> blockRowStart;
    TensorRowsMap view{};

    CvfemRowMap() = default;
    CvfemRowMap(const CvfemRowMap&)            = delete;
    CvfemRowMap& operator=(const CvfemRowMap&) = delete;

    size_t bytes() const
    {
        return elemRow.size() + elemCols.size() + blockSize.size() + rowBase.size() +
               (rowDof.size() + blockRowStart.size()) * sizeof(int);
    }
};

template<typename KeyType>
struct RowMapConn
{
    const KeyType* c[8];
};

// the DOF of element e's node i, or -1 if its row is not assembled here: the same test as the assembly's
template<typename KeyType>
__device__ __forceinline__ int rowMapDof(const RowMapConn<KeyType>& conn,
                                         size_t e,
                                         int i,
                                         const int* nodeToDof,
                                         const uint8_t* ownership,
                                         int numOwnedRows)
{
    const KeyType n = conn.c[i][e];
    const int dof   = nodeToDof[n];
    return (ownership[n] != 0 && dof >= 0 && dof < numOwnedRows) ? dof : -1;
}

// One thread per element node of a block of 32 elements: sorts their DOFs, numbers the distinct ones in
// increasing order, and writes the block's rows, where each starts in shared memory, and every element
// node's place among the rows. A block whose rows do not fit in capacity doubles, and a row longer than
// rowMapMaxLen, go straight to the matrix (rowMapDirect). Without Write it only counts the rows and their
// doubles, the direct blocks (counts[0]) and the blocks with a row too long for a compact map (counts[2]).
template<typename KeyType, bool Write>
__global__ void __launch_bounds__(256) cvfem_row_map_block_kernel(RowMapConn<KeyType> conn,
                                                                   size_t numElements,
                                                                   const int* nodeToDof,
                                                                   const uint8_t* ownership,
                                                                   int numOwnedRows,
                                                                   const int* rowPtr,
                                                                   int capacity,
                                                                   int* blockCount,
                                                                   uint16_t* blockSize,
                                                                   const int* blockRowStart,
                                                                   int* rowDof,
                                                                   uint16_t* rowBase,
                                                                   uint8_t* elemRow,
                                                                   int* counts)
{
    static_assert(tensorRowsElems * 8 == 256, "one thread per element node");
    __shared__ int key[256], rank[256], firstOf[256], total, wide;
    const int t     = threadIdx.x;
    const size_t e  = size_t(blockIdx.x) * tensorRowsElems + t / 8;
    const int dof   = e < numElements ? rowMapDof(conn, e, t % 8, nodeToDof, ownership, numOwnedRows) : -1;
    const int len   = dof >= 0 ? rowPtr[dof + 1] - rowPtr[dof] : 0;
    const bool tooLong = len > rowMapMaxLen;
    key[t]          = (dof >= 0 && !tooLong) ? dof : INT_MAX;
    if (t == 0) total = wide = 0;
    __syncthreads();

    // bitonic sort, ascending; the element nodes without a row in shared memory go to the end
    for (int k = 2; k <= 256; k <<= 1)
        for (int j = k >> 1; j > 0; j >>= 1)
        {
            const int partner = t ^ j;
            if (partner > t)
            {
                const int a = key[t], b = key[partner];
                if ((a > b) == ((t & k) == 0))
                {
                    key[t]       = b;
                    key[partner] = a;
                }
            }
            __syncthreads();
        }

    const bool first = key[t] != INT_MAX && (t == 0 || key[t - 1] != key[t]);
    rank[t]          = first;
    __syncthreads();
    for (int d = 1; d < 256; d <<= 1)
    {
        const int v = t >= d ? rank[t - d] : 0;
        __syncthreads();
        rank[t] += v;
        __syncthreads();
    }
    // each row's doubles: its right-hand side and its CSR row
    if (first) atomicAdd(&total, rowPtr[key[t] + 1] - rowPtr[key[t]] + 1);
    // rows too long for a compact map, counted in the block, so that one thread per block adds to counts
    if (!Write && len > rowMapCompactMaxLen) atomicAdd(&wide, 1);
    __syncthreads();
    const int count   = rank[255];
    const bool direct = count > tensorRowsMaxRows || total > capacity;

    if constexpr (!Write)
    {
        if (t == 0)
        {
            blockCount[blockIdx.x] = direct ? 0 : count;
            blockSize[blockIdx.x]  = direct ? 0 : uint16_t(total);
            if (direct) atomicAdd(&counts[0], 1);
            if (wide > 0) atomicAdd(&counts[2], 1);
        }
    }
    else
    {
        if (first && !direct) firstOf[rank[t] - 1] = t;
        __syncthreads();
        // where the rows start: a scan of their doubles, in rank, which is free now
        const int row   = t < count && !direct ? key[firstOf[t]] : -1;
        const int slots = row >= 0 ? rowPtr[row + 1] - rowPtr[row] + 1 : 0;
        rank[t]         = slots;
        __syncthreads();
        for (int d = 1; d < 256; d <<= 1)
        {
            const int v = t >= d ? rank[t - d] : 0;
            __syncthreads();
            rank[t] += v;
            __syncthreads();
        }
        if (row >= 0)
        {
            rowDof[blockRowStart[blockIdx.x] + t]  = row;
            rowBase[blockRowStart[blockIdx.x] + t] = uint16_t(rank[t] - slots);
        }
        if (e < numElements)
        {
            int place = rowMapNone;
            if (dof >= 0 && (direct || tooLong)) place = rowMapDirect;
            else if (dof >= 0)
            {
                int lo = 0, hi = count - 1;
                while (lo < hi)
                {
                    const int mid = (lo + hi) >> 1;
                    if (key[firstOf[mid]] < dof) lo = mid + 1;
                    else hi = mid;
                }
                place = lo;
            }
            elemRow[e * 8 + t % 8] = uint8_t(place);
        }
    }
}

// For each element row, where each of the element's 8 columns goes in the CSR row: the same rules as the
// assembly (a column without a DOF is dropped; one outside the pattern goes to the diagonal if lumped, else
// it is dropped). A row longer than rowMapMaxLen gets rowMapSearch. Compact: 4 bits per column, one word per
// element row.
template<typename KeyType, bool Compact>
__global__ void cvfem_row_map_cols_kernel(RowMapConn<KeyType> conn,
                                          size_t numElements,
                                          const int* nodeToDof,
                                          const uint8_t* ownership,
                                          int numOwnedRows,
                                          const int* rowPtr,
                                          const int* colInd,
                                          const int* diagPtr,
                                          bool lumped,
                                          uint8_t* elemCols,
                                          int* counts)
{
    const size_t t = blockIdx.x * size_t(blockDim.x) + threadIdx.x;
    if (t >= numElements * 8) return;
    const size_t e = t / 8;
    const int row  = rowMapDof(conn, e, int(t % 8), nodeToDof, ownership, numOwnedRows);
    unsigned word[2] = {~0u, ~0u};
    if (row >= 0)
    {
        const int start = rowPtr[row], len = rowPtr[row + 1] - start;
        if (len > rowMapMaxLen)
        {
            word[0] = word[1] = 0x01010101u * rowMapSearch;
            atomicAdd(&counts[1], 1);
        }
        else
        {
            const unsigned diag = unsigned(diagPtr[row] - start);
            word[0] = word[1] = 0;
            for (int c = 0; c < 8; ++c)
            {
                const int col = nodeToDof[conn.c[c][e]];
                unsigned o    = rowMapNone;
                if (col >= 0)
                {
                    o = lumped ? diag : rowMapNone;
                    for (int k = 0; k < len; ++k)
                        if (colInd[start + k] == col) o = unsigned(k);
                }
                word[c / 4] |= o << (8 * (c % 4));
            }
        }
    }
    if constexpr (Compact)
    {
        unsigned packed = 0;
        for (int c = 0; c < 8; ++c)
        {
            const unsigned o = (word[c / 4] >> (8 * (c % 4))) & 0xFF;
            packed |= (o == rowMapNone ? rowMapCompactNone : o) << (4 * c);
        }
        reinterpret_cast<unsigned*>(elemCols)[t] = packed;
    }
    else
    {
        unsigned* out = reinterpret_cast<unsigned*>(elemCols + t * 8);
        out[0]        = word[0];
        out[1]        = word[1];
    }
}

// Builds the row map on the GPU. matrix is the host copy of the device CSR struct (the pattern of the graph
// or the full assembly), and nodeToDof and ownership are those of the assembly. Every mesh gets a map: the
// blocks and rows that do not fit go straight to the matrix (numDirectBlocks, numSearchRows).
template<typename KeyType, typename RealType>
void buildCvfemRowMap(const KeyType* d_conn0,
                      const KeyType* d_conn1,
                      const KeyType* d_conn2,
                      const KeyType* d_conn3,
                      const KeyType* d_conn4,
                      const KeyType* d_conn5,
                      const KeyType* d_conn6,
                      const KeyType* d_conn7,
                      size_t numElements,
                      const int* d_nodeToDof,
                      const uint8_t* d_ownership,
                      const CSRMatrix<RealType>& matrix,
                      CvfemRowMapPattern pattern,
                      CvfemRowMap& map)
{
    const RowMapConn<KeyType> conn{{d_conn0, d_conn1, d_conn2, d_conn3, d_conn4, d_conn5, d_conn6, d_conn7}};
    const int numOwned      = matrix.numOwnedRows >= 0 ? matrix.numOwnedRows : matrix.numRows;
    const size_t incidences = numElements * 8;
    constexpr int Threads   = 256;
    const unsigned grid     = unsigned((incidences + Threads - 1) / Threads);
    map.numElements         = numElements;
    map.numBlocks           = int((numElements + tensorRowsElems - 1) / tensorRowsElems);
    map.lumped              = pattern == CvfemRowMapPattern::Graph;

    // direct blocks, searched element rows, a row too long for a compact map
    cstone::DeviceVector<int> counts(3, 0);
    // one count more than blocks, left 0, so that the scan's last entry is the total
    thrust::device_vector<int> blockCount(map.numBlocks + 1, 0);
    map.blockRowStart.resize(map.numBlocks + 1);
    map.blockSize.resize(size_t(map.numBlocks) * 2);
    if (map.numBlocks > 0)
        cvfem_row_map_block_kernel<KeyType, false><<<map.numBlocks, 256>>>(
            conn, numElements, d_nodeToDof, d_ownership, numOwned, matrix.rowPtr, tensorRowsCapacity<RealType>,
            thrust::raw_pointer_cast(blockCount.data()), reinterpret_cast<uint16_t*>(map.blockSize.data()), nullptr,
            nullptr, nullptr, nullptr, counts.data());
    cudaCheckError();
    thrust::exclusive_scan(thrust::device, blockCount.begin(), blockCount.end(), map.blockRowStart.begin());
    cudaCheckCall(cudaMemcpy(&map.numBlockRows, thrust::raw_pointer_cast(map.blockRowStart.data()) + map.numBlocks,
                             sizeof(int), cudaMemcpyDeviceToHost));
    // the format of elemCols, which sizes it
    int h_counts[3];
    cudaCheckCall(cudaMemcpy(h_counts, counts.data(), sizeof(h_counts), cudaMemcpyDeviceToHost));
    map.numDirectBlocks = h_counts[0];
    map.compact         = h_counts[0] == 0 && h_counts[2] == 0;

    map.rowDof.resize(map.numBlockRows);
    map.rowBase.resize(size_t(map.numBlockRows) * 2);
    map.elemRow.resize(incidences);
    map.elemCols.resize(incidences * (map.compact ? 4 : 8));
    if (map.numBlocks > 0)
        cvfem_row_map_block_kernel<KeyType, true><<<map.numBlocks, 256>>>(
            conn, numElements, d_nodeToDof, d_ownership, numOwned, matrix.rowPtr, tensorRowsCapacity<RealType>,
            nullptr, nullptr, thrust::raw_pointer_cast(map.blockRowStart.data()), map.rowDof.data(),
            reinterpret_cast<uint16_t*>(map.rowBase.data()), map.elemRow.data(), counts.data());
    cudaCheckError();
    if (incidences > 0)
    {
        auto cols = map.compact ? cvfem_row_map_cols_kernel<KeyType, true> : cvfem_row_map_cols_kernel<KeyType, false>;
        cols<<<grid, Threads>>>(conn, numElements, d_nodeToDof, d_ownership, numOwned, matrix.rowPtr, matrix.colInd,
                                matrix.diagPtr, map.lumped, map.elemCols.data(), counts.data());
    }
    cudaCheckError();
    // a map rebuilt in place must not keep the count of its last mesh; a compact map has no searched rows
    map.numSearchRows = 0;
    if (!map.compact)
        cudaCheckCall(cudaMemcpy(&map.numSearchRows, counts.data() + 1, sizeof(int), cudaMemcpyDeviceToHost));
    map.view = {map.elemRow.data(),
                map.elemCols.data(),
                thrust::raw_pointer_cast(map.blockRowStart.data()),
                reinterpret_cast<const uint16_t*>(map.blockSize.data()),
                reinterpret_cast<const uint16_t*>(map.rowBase.data()),
                map.rowDof.data(),
                map.lumped};
}

} // namespace fem
} // namespace mars
