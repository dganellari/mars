#pragma once

#include "mars_cvfem_hex_kernel_tensor_faces.hpp"
#include "mars_cvfem_row_map_format.hpp"

namespace mars {
namespace fem {

// rows per block: the place of a row is a byte, and rowMapDirect and rowMapNone are not places
constexpr int tensorRowsMaxRows = 0xFE;
constexpr int tensorRowsElems   = 32;
// the doubles of a block's rows in shared memory: the element storage, less each row's DOF, CSR start and
// base, and one byte per 8 doubles for the flush (tensorRowsChunk). A block whose rows need more adds them
// straight to the matrix.
constexpr int tensorRowsChunk = 8;
template<typename RealType>
constexpr int tensorRowsCapacity =
    int((sizeof(TensorFacesElement<RealType>) * tensorRowsElems - (3 * tensorRowsMaxRows + 1) * sizeof(int)) *
        tensorRowsChunk / (tensorRowsChunk * sizeof(RealType) + 1)) /
    tensorRowsChunk * tensorRowsChunk;

template<typename KeyType>
__device__ __forceinline__ KeyType tensorRowsNode(int i,
                                                  size_t e,
                                                  const KeyType* c0,
                                                  const KeyType* c1,
                                                  const KeyType* c2,
                                                  const KeyType* c3,
                                                  const KeyType* c4,
                                                  const KeyType* c5,
                                                  const KeyType* c6,
                                                  const KeyType* c7)
{
    switch (i)
    {
        case 0: return c0[e];
        case 1: return c1[e];
        case 2: return c2[e];
        case 3: return c3[e];
        case 4: return c4[e];
        case 5: return c5[e];
        case 6: return c6[e];
        default: return c7[e];
    }
}

// the place of column col in a sorted CSR row of len entries, or miss
__device__ __forceinline__ int tensorRowsFind(const int* colInd, int len, int col, int miss)
{
    int lo = 0, hi = len - 1;
    while (lo <= hi)
    {
        const int mid = (lo + hi) >> 1, c = tensorFacesLoadGlobal(colInd + mid);
        if (c == col) return mid;
        if (c < col) lo = mid + 1;
        else hi = mid - 1;
    }
    return miss;
}

// Rows (32 elements per block, graph or full assembly): the block adds up its elements' contributions in
// shared memory, laid out as the CSR rows of the DOFs it touches, then adds each row whole to the matrix, in
// DOF order. Unlike the aggregated tensor_faces kernel it has no hash table and no CSR search: the row map
// gives every contribution its place. Consecutive threads add consecutive entries of a row, and on MI300A,
// where every global atomic goes past the L2, that cuts the atomic requests 1.35e8 to 3.3e7 (hex192, graph):
// 8.2 to 4.5 ms. Storing the rows only one block touches, or replacing the atomics by partial rows and
// a reduction kernel, was slower (tests/performance/AMD-assembly.md). The rows the map marks rowMapDirect
// go straight to the matrix, with their own atomics. Compact must be the map's (CvfemRowMap::compact): its
// 4-bit columns, and no direct or searched rows.
template<typename KeyType, typename RealType, bool Compact>
__global__ void __launch_bounds__(512) cvfem_hex_assembly_kernel_tensor_rows(
    const KeyType* __restrict__ d_conn0,
    const KeyType* __restrict__ d_conn1,
    const KeyType* __restrict__ d_conn2,
    const KeyType* __restrict__ d_conn3,
    const KeyType* __restrict__ d_conn4,
    const KeyType* __restrict__ d_conn5,
    const KeyType* __restrict__ d_conn6,
    const KeyType* __restrict__ d_conn7,
    size_t numElements,
    const RealType* __restrict__ d_x,
    const RealType* __restrict__ d_y,
    const RealType* __restrict__ d_z,
    const RealType* __restrict__ d_gamma,
    const RealType* __restrict__ d_phi,
    const RealType* __restrict__ d_beta,
    const RealType* __restrict__ d_grad_phi_x,
    const RealType* __restrict__ d_grad_phi_y,
    const RealType* __restrict__ d_grad_phi_z,
    const RealType* __restrict__ d_mdot,
    const RealType* __restrict__ d_areaVec_x,
    const RealType* __restrict__ d_areaVec_y,
    const RealType* __restrict__ d_areaVec_z,
    const int* __restrict__ d_node_to_dof,
    CSRMatrix<RealType>* matrix,
    RealType* __restrict__ d_rhs,
    TensorRowsMap map)
{
    constexpr int BlockSize = 512, MaxRows = tensorRowsMaxRows, Capacity = tensorRowsCapacity<RealType>;
    constexpr int Chunk = tensorRowsChunk;
    static_assert(BlockSize == tensorRowsElems * 16, "16 threads per element");
    __shared__ TensorFacesElement<RealType> shared[tensorRowsElems];
    static_assert(Capacity * sizeof(RealType) + (3 * MaxRows + 1) * sizeof(int) + Capacity / Chunk <= sizeof(shared),
                  "the rows fit in the element storage");

    const int lane       = threadIdx.x % 16;
    const size_t elemIdx = size_t(blockIdx.x) * tensorRowsElems + threadIdx.x / 16;
    const bool active    = elemIdx < numElements;

    // The block's rows, one per thread, and this lane's map entries: loaded first, so that the element
    // work hides their latency.
    const int rowBegin = map.blockRowStart[blockIdx.x];
    const int numRows  = map.blockRowStart[blockIdx.x + 1] - rowBegin;
    const int size     = map.blockSize[blockIdx.x];
    int rowDof = 0, rowStart = 0, rowBase = 0, rowEnd = 0;
    if (int(threadIdx.x) < numRows)
    {
        rowDof   = map.rowDof[rowBegin + threadIdx.x];
        rowStart = tensorFacesLoadGlobal(&matrix->rowPtr[rowDof]);
        rowBase  = map.rowBase[rowBegin + threadIdx.x];
        rowEnd   = int(threadIdx.x) + 1 < numRows ? map.rowBase[rowBegin + threadIdx.x + 1] : size;
    }
    const int i     = lane & 7;
    const int col0  = (lane >> 3) * 4;
    const int place = active ? map.elemRow[elemIdx * 8 + i] : rowMapNone;
    // this lane's 4 columns: one aligned word of bytes, or the 4-bit word of the row, shifted to them
    unsigned cols = ~0u;
    if (place != rowMapNone)
        cols = Compact ? reinterpret_cast<const unsigned*>(map.elemCols)[elemIdx * 8 + i] >> (4 * col0)
                       : *reinterpret_cast<const unsigned*>(map.elemCols + elemIdx * 64 + i * 8 + col0);
    int dof = 0, start = 0;
    if (!Compact && place == rowMapDirect)
    {
        dof   = d_node_to_dof[tensorRowsNode(i, elemIdx, d_conn0, d_conn1, d_conn2, d_conn3, d_conn4, d_conn5,
                                             d_conn6, d_conn7)];
        start = tensorFacesLoadGlobal(&matrix->rowPtr[dof]);
    }

    RealType lhs[4], rhs;
    tensorFacesRow<KeyType, RealType, false>(shared[threadIdx.x / 16], lane, elemIdx, active, d_conn0, d_conn1,
                                             d_conn2, d_conn3, d_conn4, d_conn5, d_conn6, d_conn7, d_x, d_y, d_z,
                                             d_gamma, d_phi, d_beta, d_grad_phi_x, d_grad_phi_y, d_grad_phi_z,
                                             d_mdot, d_areaVec_x, d_areaVec_y, d_areaVec_z, nullptr, nullptr, lhs,
                                             rhs);

    // The element storage is dead now: it holds the rows, then each row's DOF, CSR start and base (one more:
    // the block's size), then for every 8 doubles the row of the first one.
    RealType* acc     = reinterpret_cast<RealType*>(shared);
    int* infoDof      = reinterpret_cast<int*>(acc + Capacity);
    int* infoStart    = infoDof + MaxRows;
    int* infoBase     = infoStart + MaxRows;
    uint8_t* chunkRow = reinterpret_cast<uint8_t*>(infoBase + MaxRows + 1);
    __syncthreads();
    for (int t = threadIdx.x; t < size; t += BlockSize)
        acc[t] = 0;
    if (int(threadIdx.x) < numRows)
    {
        infoDof[threadIdx.x]   = rowDof;
        infoStart[threadIdx.x] = rowStart;
        infoBase[threadIdx.x]  = rowBase;
        for (int c = (rowBase + Chunk - 1) / Chunk; c * Chunk < rowEnd; ++c)
            chunkRow[c] = uint8_t(threadIdx.x);
    }
    if (threadIdx.x == 0) infoBase[numRows] = size;
    __syncthreads();

    if (place != rowMapNone)
    {
        int offset[4]; // -1: dropped
#pragma unroll
        for (int j = 0; j < 4; ++j)
        {
            const unsigned b = Compact ? (cols >> (4 * j)) & 0xF : (cols >> (8 * j)) & 0xFF;
            offset[j]        = b == (Compact ? rowMapCompactNone : rowMapNone) ? -1 : int(b);
        }
        if (!Compact && offset[0] == rowMapSearch)
        {
            const int len  = tensorFacesLoadGlobal(&matrix->rowPtr[dof + 1]) - start;
            const int miss = map.lumped ? tensorFacesLoadGlobal(&matrix->diagPtr[dof]) - start : -1;
#pragma unroll
            for (int j = 0; j < 4; ++j)
            {
                const int col = d_node_to_dof[tensorRowsNode(col0 + j, elemIdx, d_conn0, d_conn1, d_conn2, d_conn3,
                                                             d_conn4, d_conn5, d_conn6, d_conn7)];
                offset[j]     = col < 0 ? -1 : tensorRowsFind(matrix->colInd + start, len, col, miss);
            }
        }
        // the columns lumped into the diagonal share its offset: one atomic for all of them
#pragma unroll
        for (int j = 1; j < 4; ++j)
#pragma unroll
            for (int k = 0; k < j; ++k)
                if (offset[j] >= 0 && offset[k] == offset[j])
                {
                    lhs[k] += lhs[j];
                    offset[j] = -1;
                }
        if (!Compact && place == rowMapDirect)
        {
            RealType* row = matrix->values + start;
#pragma unroll
            for (int j = 0; j < 4; ++j)
                if (offset[j] >= 0 && lhs[j] != 0) tensorFacesAtomicAddGlobal(row + offset[j], lhs[j]);
            if (col0 == 0 && rhs != 0) tensorFacesAtomicAddGlobal(d_rhs + dof, rhs);
        }
        else
        {
            RealType* r = acc + infoBase[place];
#pragma unroll
            for (int j = 0; j < 4; ++j)
                if (offset[j] >= 0) atomicAdd(&r[1 + offset[j]], lhs[j]);
            if (col0 == 0) atomicAdd(r, rhs);
        }
    }
    __syncthreads();

    RealType* values = matrix->values;
    for (int t = threadIdx.x; t < size; t += BlockSize)
    {
        const RealType v = acc[t];
        if (v == 0) continue;
        // every row has at least one double, so this steps at most 7 rows
        int r = chunkRow[t / Chunk];
        while (t >= infoBase[r + 1])
            ++r;
        const int e = t - infoBase[r];
        tensorFacesAtomicAddGlobal(e == 0 ? d_rhs + infoDof[r] : values + infoStart[r] + e - 1, v);
    }
}

} // namespace fem
} // namespace mars
