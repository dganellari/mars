#pragma once

#include "mars_cvfem_hex_kernel_tensor_rows.hpp"
#include "../utils/mars_domain_utils.hpp"
#include "cstone/cuda/device_vector.h"
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/scan.h>
#include <climits>
#include <cstdint>
#include <stdexcept>

namespace mars {
namespace fem {

// The row map of the graph assembly (mars_cvfem_hex_kernel_tensor_rows.hpp, and the column offsets alone
// for cvfem_hex_assembly_kernel_tensor), built once on the GPU for a mesh, a DOF map and a graph sparsity
// pattern: 40 bytes per element and 4 per block row.
struct CvfemRowMap
{
    size_t numElements = 0;
    int numBlocks      = 0;
    int numBlockRows   = 0; // the rows of all blocks together
    cstone::DeviceVector<uint8_t> elemRow;
    cstone::DeviceVector<unsigned> elemCols;
    cstone::DeviceVector<int> rowDof;
    // a thrust vector: scanned with the sparsity builder's own thrust call, so no new scan kernels compile
    thrust::device_vector<int> blockRowStart;
    TensorRowsMap view{};

    CvfemRowMap() = default;
    CvfemRowMap(const CvfemRowMap&)            = delete;
    CvfemRowMap& operator=(const CvfemRowMap&) = delete;

    size_t bytes() const
    {
        return elemRow.size() + elemCols.size() * sizeof(unsigned) + (rowDof.size() + blockRowStart.size()) * sizeof(int);
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
// increasing order, and writes the block's rows and every element node's place among them. Without Write
// it only counts the rows.
template<typename KeyType, bool Write>
__global__ void __launch_bounds__(256) cvfem_row_map_block_kernel(RowMapConn<KeyType> conn,
                                                                   size_t numElements,
                                                                   const int* nodeToDof,
                                                                   const uint8_t* ownership,
                                                                   int numOwnedRows,
                                                                   int* blockCount,
                                                                   const int* blockRowStart,
                                                                   int* rowDof,
                                                                   uint8_t* elemRow,
                                                                   int* error)
{
    static_assert(tensorRowsElems * 8 == 256, "one thread per element node");
    __shared__ int key[256], rank[256], firstOf[256];
    const int t    = threadIdx.x;
    const size_t e = size_t(blockIdx.x) * tensorRowsElems + t / 8;
    const int dof  = e < numElements ? rowMapDof(conn, e, t % 8, nodeToDof, ownership, numOwnedRows) : -1;
    key[t]         = dof >= 0 ? dof : INT_MAX;
    __syncthreads();

    // bitonic sort, ascending; the element nodes without a row go to the end
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
    const int count = rank[255];

    if constexpr (!Write)
    {
        if (t == 0)
        {
            blockCount[blockIdx.x] = count;
            if (count > tensorRowsMaxRows) atomicMax(error, 1);
        }
    }
    else
    {
        if (first) firstOf[rank[t] - 1] = t;
        __syncthreads();
        if (t < count) rowDof[blockRowStart[blockIdx.x] + t] = key[firstOf[t]];
        if (e < numElements)
        {
            int place = 0xFF;
            if (dof >= 0)
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
// lumped graph assembly (a column without a DOF is dropped, one outside the pattern goes to the diagonal).
template<typename KeyType>
__global__ void cvfem_row_map_cols_kernel(RowMapConn<KeyType> conn,
                                          size_t numElements,
                                          const int* nodeToDof,
                                          const uint8_t* ownership,
                                          int numOwnedRows,
                                          const int* rowPtr,
                                          const int* colInd,
                                          const int* diagPtr,
                                          unsigned* elemCols,
                                          int* error)
{
    const size_t t = blockIdx.x * size_t(blockDim.x) + threadIdx.x;
    if (t >= numElements * 8) return;
    const size_t e = t / 8;
    const int row  = rowMapDof(conn, e, int(t % 8), nodeToDof, ownership, numOwnedRows);
    unsigned word  = ~0u;
    if (row >= 0)
    {
        const int start = rowPtr[row], len = rowPtr[row + 1] - start;
        if (len > tensorRowsStride - 1) atomicMax(error, 2);
        else
        {
            const unsigned diag = unsigned(diagPtr[row] - start);
            word                = 0;
            for (int c = 0; c < 8; ++c)
            {
                const int col = nodeToDof[conn.c[c][e]];
                unsigned o    = 0xF;
                if (col >= 0)
                {
                    o = diag;
                    for (int k = 0; k < len; ++k)
                        if (colInd[start + k] == col) o = unsigned(k);
                }
                word |= o << (4 * c);
            }
        }
    }
    elemCols[t] = word;
}

// Builds the row map on the GPU. matrix is the host copy of the device CSR struct (the graph pattern), and
// nodeToDof and ownership are those of the assembly. Throws when a block of 32 elements touches more than
// 255 rows or a row has more than 8 entries: the kernels that use the map do not apply then.
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
                      CvfemRowMap& map)
{
    const RowMapConn<KeyType> conn{{d_conn0, d_conn1, d_conn2, d_conn3, d_conn4, d_conn5, d_conn6, d_conn7}};
    const int numOwned      = matrix.numOwnedRows >= 0 ? matrix.numOwnedRows : matrix.numRows;
    const size_t incidences = numElements * 8;
    constexpr int Threads   = 256;
    const unsigned grid     = unsigned((incidences + Threads - 1) / Threads);
    map.numElements         = numElements;
    map.numBlocks           = int((numElements + tensorRowsElems - 1) / tensorRowsElems);

    auto failed = [](const cstone::DeviceVector<int>& error) {
        int code = 0;
        cudaCheckCall(cudaMemcpy(&code, error.data(), sizeof(int), cudaMemcpyDeviceToHost));
        return code;
    };
    cstone::DeviceVector<int> error(1, 0);

    // one count more than blocks, left 0, so that the scan's last entry is the total
    thrust::device_vector<int> blockCount(map.numBlocks + 1, 0);
    map.blockRowStart.resize(map.numBlocks + 1);
    if (map.numBlocks > 0)
        cvfem_row_map_block_kernel<KeyType, false><<<map.numBlocks, 256>>>(
            conn, numElements, d_nodeToDof, d_ownership, numOwned, thrust::raw_pointer_cast(blockCount.data()),
            nullptr, nullptr, nullptr, error.data());
    cudaCheckError();
    if (failed(error) == 1) throw std::runtime_error("row map: a block of 32 elements touches more than 255 rows");
    thrust::exclusive_scan(thrust::device, blockCount.begin(), blockCount.end(), map.blockRowStart.begin());
    cudaCheckCall(cudaMemcpy(&map.numBlockRows, thrust::raw_pointer_cast(map.blockRowStart.data()) + map.numBlocks,
                             sizeof(int), cudaMemcpyDeviceToHost));

    map.rowDof.resize(map.numBlockRows);
    map.elemRow.resize(incidences);
    map.elemCols.resize(incidences);
    if (map.numBlocks > 0)
        cvfem_row_map_block_kernel<KeyType, true><<<map.numBlocks, 256>>>(
            conn, numElements, d_nodeToDof, d_ownership, numOwned, nullptr,
            thrust::raw_pointer_cast(map.blockRowStart.data()), map.rowDof.data(), map.elemRow.data(), error.data());
    cudaCheckError();
    if (incidences > 0)
        cvfem_row_map_cols_kernel<<<grid, Threads>>>(conn, numElements, d_nodeToDof, d_ownership, numOwned,
                                                     matrix.rowPtr, matrix.colInd, matrix.diagPtr,
                                                     map.elemCols.data(), error.data());
    cudaCheckError();
    if (failed(error) == 2) throw std::runtime_error("row map: a CSR row has more than 8 entries");

    cudaCheckCall(cudaDeviceSynchronize());
    map.view = {map.elemRow.data(), map.elemCols.data(), thrust::raw_pointer_cast(map.blockRowStart.data()),
                map.rowDof.data()};
}

} // namespace fem
} // namespace mars
