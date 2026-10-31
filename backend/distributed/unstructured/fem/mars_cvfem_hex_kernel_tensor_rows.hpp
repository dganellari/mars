#pragma once

#include "mars_cvfem_hex_kernel_tensor_faces.hpp"

namespace mars {
namespace fem {

// The row map (mars_cvfem_row_map.hpp) as the kernels read it. Block b, 32 elements, touches the rows
// rowDof[blockRowStart[b] .. blockRowStart[b + 1]), in increasing DOF order. elemRow gives each element
// node's place in that list (0xFF: the row is not assembled here); elemCols, for each element row, the
// CSR offset in the row that each of the element's 8 columns adds to, 4 bits per column (0xF: dropped),
// and ~0u for a row not assembled here.
struct TensorRowsMap
{
    const uint8_t* elemRow;
    const unsigned* elemCols;
    const int* blockRowStart;
    const int* rowDof;
};

// rows per block: the place of a row is a byte, and 0xFF marks the rows that are not assembled
constexpr int tensorRowsMaxRows = 255;
// a row in shared memory: the CSR row, at most 8 entries, then the right-hand side. 9 and not 8 doubles
// also spreads the rows over the shared-memory banks.
constexpr int tensorRowsStride = 9;
constexpr int tensorRowsElems  = 32;

// Rows (graph mode, 32 elements per block): the block adds up its elements' contributions in shared memory,
// laid out as the CSR rows of the DOFs it touches, then adds each row whole to the matrix, in DOF order.
// Unlike the aggregated tensor_faces kernel it has no hash table and no CSR search: the row map gives
// every contribution its place. Consecutive threads add consecutive entries of a row, and on MI300A,
// where every global atomic goes past the L2, that cuts the atomic requests 1.35e8 to 3.3e7 (hex192):
// 8.2 to 4.5 ms. Storing the rows only one block touches, or replacing the atomics by partial rows and
// a reduction kernel, was slower (tests/performance/AMD-assembly.md).
template<typename KeyType, typename RealType>
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
    CSRMatrix<RealType>* matrix,
    RealType* __restrict__ d_rhs,
    TensorRowsMap map)
{
    constexpr int BlockSize = 512, Stride = tensorRowsStride, MaxRows = tensorRowsMaxRows;
    static_assert(BlockSize == tensorRowsElems * 16, "16 threads per element");
    __shared__ TensorFacesElement<RealType> shared[tensorRowsElems];
    static_assert(sizeof(shared) >= MaxRows * (Stride * sizeof(RealType) + 3 * sizeof(int)),
                  "the rows reuse the element storage");

    const int lane       = threadIdx.x % 16;
    const size_t elemIdx = size_t(blockIdx.x) * tensorRowsElems + threadIdx.x / 16;
    const bool active    = elemIdx < numElements;

    // The block's rows, one per thread, and this lane's map entries: loaded first, so that the element
    // work hides their latency.
    const int rowBegin = map.blockRowStart[blockIdx.x];
    const int numRows  = map.blockRowStart[blockIdx.x + 1] - rowBegin;
    int rowDof = 0, rowStart = 0, rowLen = 0;
    if (int(threadIdx.x) < numRows)
    {
        rowDof   = map.rowDof[rowBegin + threadIdx.x];
        rowStart = tensorFacesLoadGlobal(&matrix->rowPtr[rowDof]);
        rowLen   = tensorFacesLoadGlobal(&matrix->rowPtr[rowDof + 1]) - rowStart;
    }
    const int i         = lane & 7;
    const int row       = active ? map.elemRow[elemIdx * 8 + i] : 0xFF;
    const unsigned cols = (row != 0xFF) ? map.elemCols[elemIdx * 8 + i] : 0u;

    RealType lhs[4], rhs;
    tensorFacesRow<KeyType, RealType, false>(shared[threadIdx.x / 16], lane, elemIdx, active, d_conn0, d_conn1,
                                             d_conn2, d_conn3, d_conn4, d_conn5, d_conn6, d_conn7, d_x, d_y, d_z,
                                             d_gamma, d_phi, d_beta, d_grad_phi_x, d_grad_phi_y, d_grad_phi_z,
                                             d_mdot, d_areaVec_x, d_areaVec_y, d_areaVec_z, nullptr, nullptr, lhs,
                                             rhs);

    // The element storage is dead now: it holds the rows, then each row's DOF, CSR start and length.
    RealType* acc  = reinterpret_cast<RealType*>(shared);
    int* infoDof   = reinterpret_cast<int*>(acc + MaxRows * Stride);
    int* infoStart = infoDof + MaxRows;
    int* infoLen   = infoStart + MaxRows;
    __syncthreads();
    for (int t = threadIdx.x; t < numRows * Stride; t += BlockSize)
        acc[t] = 0;
    if (int(threadIdx.x) < numRows)
    {
        infoDof[threadIdx.x]   = rowDof;
        infoStart[threadIdx.x] = rowStart;
        infoLen[threadIdx.x]   = rowLen;
    }
    __syncthreads();

    if (row != 0xFF)
    {
        const int col0 = (lane >> 3) * 4;
        int offset[4];
#pragma unroll
        for (int j = 0; j < 4; ++j)
            offset[j] = (cols >> (4 * (col0 + j))) & 0xF;
        // the columns lumped into the diagonal share its offset: one shared atomic for all of them
#pragma unroll
        for (int j = 1; j < 4; ++j)
#pragma unroll
            for (int k = 0; k < j; ++k)
                if (offset[j] != 0xF && offset[k] == offset[j])
                {
                    lhs[k] += lhs[j];
                    offset[j] = 0xF;
                }
        RealType* r = acc + row * Stride;
#pragma unroll
        for (int j = 0; j < 4; ++j)
            if (offset[j] != 0xF) atomicAdd(&r[offset[j]], lhs[j]);
        if (col0 == 0) atomicAdd(&r[Stride - 1], rhs);
    }
    __syncthreads();

    RealType* values = matrix->values;
    for (int t = threadIdx.x; t < numRows * Stride; t += BlockSize)
    {
        const int r = t / Stride, e = t - r * Stride;
        if (e < Stride - 1 && e >= infoLen[r]) continue;
        const RealType v = acc[t];
        if (v != 0) tensorFacesAtomicAddGlobal((e == Stride - 1) ? d_rhs + infoDof[r] : values + infoStart[r] + e, v);
    }
}

} // namespace fem
} // namespace mars
