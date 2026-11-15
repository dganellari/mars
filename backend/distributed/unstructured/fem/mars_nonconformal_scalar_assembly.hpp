#pragma once
#if defined(__CUDACC__) || defined(__HIPCC__)
#include "../utils/mars_gpu_runtime.hpp"
#endif
#include "mars_nonconformal_scalar.hpp"

namespace mars::fem
{
// Sorted owned global row IDs, sorted global columns within each CSR row. Remote columns remain
// global until the linear-system adapter builds its column/halo map.
template<class RealType>
struct InterfaceCsrView
{
    int rows, entries;
    const int64_t* rowIds;
    const int* offsets;
    const int64_t* columns;
    RealType* values;
};

namespace nonconformal_detail
{
MARS_FS_FN inline int columnSlot(const int64_t* columns, int begin, int end, int64_t column)
{
    int stop = end;
    while (begin < end)
    {
        int mid = begin + (end - begin) / 2;
        if (columns[mid] < column)
            begin = mid + 1;
        else
            end = mid;
    }
    return begin < stop && columns[begin] == column ? begin : -1;
}

// Each row has one writer. Point lists are sorted and contain every stencil exactly once;
// building this inverse map once avoids atomics and preserves the accumulation order.
template<class RealType>
MARS_FS_FN bool assembleScalarRow(int row,
                                  InterfaceCsrView<RealType> matrix,
                                  int count,
                                  const InterfaceStencil<RealType>* stencils,
                                  const int* pointOffsets,
                                  const int* pointIds)
{
    if (row < 0 || row >= matrix.rows || matrix.entries < 0 || count < 0 || !matrix.rowIds || !matrix.offsets ||
        !pointOffsets || (matrix.entries && (!matrix.columns || !matrix.values)) || (count && (!stencils || !pointIds)))
        return false;
    if (matrix.offsets[0] != 0 || matrix.offsets[matrix.rows] != matrix.entries || pointOffsets[0] != 0 ||
        pointOffsets[matrix.rows] != count || matrix.rowIds[row] < 0 ||
        (row && matrix.rowIds[row - 1] >= matrix.rowIds[row]))
        return false;
    int begin = matrix.offsets[row], end = matrix.offsets[row + 1];
    int first = pointOffsets[row], last = pointOffsets[row + 1];
    if (begin < 0 || end < begin || end > matrix.entries || first < 0 || last < first || last > count) return false;
    for (int k = begin; k < end; ++k)
        if (matrix.columns[k] < 0 || (k > begin && matrix.columns[k - 1] >= matrix.columns[k]) ||
            !finite(matrix.values[k]))
            return false;
    // Validate every column before changing this row, so a missing remote entry cannot be dropped.
    for (int i = first; i < last; ++i)
    {
        int p = pointIds[i];
        if (p < 0 || p >= count || (i > first && pointIds[i - 1] >= p)) return false;
        const auto& s = stencils[p];
        if (s.status != InterfaceStatus::valid || s.row != matrix.rowIds[row]) return false;
        for (int j = 0; j < 8; ++j)
            if (!finite(s.value[j]) || columnSlot(matrix.columns, begin, end, s.column[j]) < 0) return false;
    }
    for (int i = first; i < last; ++i)
    {
        const auto& s = stencils[pointIds[i]];
        for (int j = 0; j < 8; ++j)
        {
            int k = columnSlot(matrix.columns, begin, end, s.column[j]);
            matrix.values[k] += s.value[j];
            if (!finite(matrix.values[k])) return false;
        }
    }
    return true;
}
} // namespace nonconformal_detail

template<class RealType>
bool assembleNonconformalScalarRowsHost(InterfaceCsrView<RealType> matrix,
                                        int count,
                                        const InterfaceStencil<RealType>* stencils,
                                        const int* pointOffsets,
                                        const int* pointIds)
{
    if (matrix.rows < 0 || matrix.entries < 0 || count < 0) return false;
    if (matrix.rows == 0) return count == 0 && matrix.entries == 0;
    bool valid = true;
    for (int row = 0; row < matrix.rows; ++row)
        valid = nonconformal_detail::assembleScalarRow(row, matrix, count, stencils, pointOffsets, pointIds) && valid;
    return valid;
}

#if defined(__CUDACC__) || defined(__HIPCC__)
namespace nonconformal_detail
{
template<class RealType>
__global__ void scalarRowsKernel(InterfaceCsrView<RealType> matrix,
                                 int count,
                                 const InterfaceStencil<RealType>* stencils,
                                 const int* offsets,
                                 const int* ids,
                                 int* valid)
{
    size_t row = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (row < size_t(matrix.rows)) valid[row] = assembleScalarRow(int(row), matrix, count, stencils, offsets, ids);
}
} // namespace nonconformal_detail
// Reduce rowValid collectively before using the matrix. On any rejection, discard the assembly;
// other valid rows may already have been added. Storage lengths must match the view and counts.
template<class RealType>
cudaError_t assembleNonconformalScalarRows(InterfaceCsrView<RealType> matrix,
                                           int count,
                                           const InterfaceStencil<RealType>* stencils,
                                           const int* pointOffsets,
                                           const int* pointIds,
                                           int* rowValid,
                                           cudaStream_t stream = 0)
{
    if (matrix.rows < 0 || matrix.entries < 0 || count < 0) return cudaErrorInvalidValue;
    if (!matrix.rows) return count == 0 && matrix.entries == 0 ? cudaSuccess : cudaErrorInvalidValue;
    if (!rowValid) return cudaErrorInvalidValue;
    nonconformal_detail::scalarRowsKernel<<<(matrix.rows - 1) / 128 + 1, 128, 0, stream>>>(
        matrix, count, stencils, pointOffsets, pointIds, rowValid);
    return cudaGetLastError();
}
#endif
} // namespace mars::fem
