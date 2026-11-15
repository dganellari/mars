#pragma once
#include "mars_nonconformal_scalar_assembly.hpp"
#include <algorithm>
#include <climits>
#include <numeric>
#include <vector>
#if defined(__CUDACC__) || defined(__HIPCC__)
#include "../utils/mars_gpu_thrust.hpp"
#include <stdexcept>
#include <thrust/copy.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/binary_search.h>
#include <thrust/count.h>
#include <thrust/for_each.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/sort.h>
#include <thrust/transform.h>
#include <thrust/unique.h>
#endif

namespace mars::fem
{
namespace nonconformal_detail
{
struct GraphEntry
{
    int row;
    int64_t column;
    MARS_FS_FN bool operator<(const GraphEntry& b) const { return row < b.row || (row == b.row && column < b.column); }
    MARS_FS_FN bool operator==(const GraphEntry& b) const { return row == b.row && column == b.column; }
};
struct GraphRow
{
    MARS_FS_FN int operator()(GraphEntry e) const { return e.row; }
};
struct GraphColumn
{
    MARS_FS_FN int64_t operator()(GraphEntry e) const { return e.column; }
};
struct RowKey
{
    MARS_FS_FN GraphEntry operator()(int r) const { return {r, -1}; }
};
template<class T>
bool graphSizes(InterfaceCsrView<T> base, int count, const InterfaceStencil<T>* stencils)
{
    return base.rows >= 0 && base.rows < INT_MAX && base.entries >= 0 && count >= 0 &&
           int64_t(base.rows) + base.entries + int64_t(8) * count <= INT_MAX &&
           (base.rows || (!base.entries && !count)) && base.offsets && (!base.rows || base.rowIds) &&
           (!base.entries || (base.columns && base.values)) && (!count || stencils);
}
template<class T>
struct InvalidGraphRow
{
    InterfaceCsrView<T> base;
    MARS_FS_FN bool operator()(int r) const
    {
        if (base.offsets[0] != 0 || base.offsets[base.rows] != base.entries) return true;
        if (r == base.rows) return false;
        int first = base.offsets[r], last = base.offsets[r + 1];
        if (base.rowIds[r] < 0 || (r && base.rowIds[r - 1] >= base.rowIds[r]) || first < 0 || last < first ||
            last > base.entries)
            return true;
        for (int k = first; k < last; ++k)
            if (base.columns[k] < 0 || (k > first && base.columns[k - 1] >= base.columns[k]) || !finite(base.values[k]))
                return true;
        return false;
    }
};
template<class T>
struct InvalidGraphPoint
{
    InterfaceCsrView<T> base;
    const InterfaceStencil<T>* stencils;
    MARS_FS_FN bool operator()(int p) const
    {
        const auto& s = stencils[p];
        if (s.status != InterfaceStatus::valid || columnSlot(base.rowIds, 0, base.rows, s.row) < 0) return true;
        for (int j = 0; j < 8; ++j)
            if (s.column[j] < 0 || !finite(s.value[j])) return true;
        return false;
    }
};
template<class T>
struct BaseGraphEntries
{
    InterfaceCsrView<T> base;
    GraphEntry* entries;
    MARS_FS_FN void operator()(int r) const
    {
        entries[r] = {r, base.rowIds[r]};
        for (int k = base.offsets[r]; k < base.offsets[r + 1]; ++k)
            entries[base.rows + k] = {r, base.columns[k]};
    }
};
template<class T>
struct PointGraphEntries
{
    InterfaceCsrView<T> base;
    const InterfaceStencil<T>* stencils;
    GraphEntry* entries;
    GraphEntry* points;
    MARS_FS_FN void operator()(int p) const
    {
        int r     = columnSlot(base.rowIds, 0, base.rows, stencils[p].row);
        points[p] = {r, p};
        for (int j = 0; j < 8; ++j)
            entries[base.rows + base.entries + 8 * p + j] = {r, stencils[p].column[j]};
    }
};
template<class T>
struct CopyBaseValues
{
    InterfaceCsrView<T> base, output;
    MARS_FS_FN void operator()(int r) const
    {
        for (int k = base.offsets[r]; k < base.offsets[r + 1]; ++k)
        {
            int slot            = columnSlot(output.columns, output.offsets[r], output.offsets[r + 1], base.columns[k]);
            output.values[slot] = base.values[k];
        }
    }
};
} // namespace nonconformal_detail

template<class T>
struct NonconformalScalarGraphHost
{
    std::vector<int64_t> rowIds, columns;
    std::vector<int> offsets, pointOffsets, pointIds;
    std::vector<T> values;
    InterfaceCsrView<T> view()
    {
        return {int(rowIds.size()), int(columns.size()), rowIds.data(), offsets.data(), columns.data(), values.data()};
    }
};
// The base graph contains owned rows and global columns. Add all donor columns and a diagonal
// for each row, retaining base coefficients. Rebuild after pairing changes, not numeric refreshes.
// Input storage must match its declared lengths. Rejection leaves output unchanged.
template<class T>
bool buildNonconformalScalarGraphHost(InterfaceCsrView<T> base,
                                      int count,
                                      const InterfaceStencil<T>* stencils,
                                      NonconformalScalarGraphHost<T>& output)
{
    using namespace nonconformal_detail;
    if (!graphSizes(base, count, stencils)) return false;
    for (int r = 0; r <= base.rows; ++r)
        if (InvalidGraphRow<T>{base}(r)) return false;
    for (int p = 0; p < count; ++p)
        if (InvalidGraphPoint<T>{base, stencils}(p)) return false;
    std::vector<GraphEntry> entries(base.rows + base.entries + 8 * count), points(count);
    for (int r = 0; r < base.rows; ++r)
        BaseGraphEntries<T>{base, entries.data()}(r);
    for (int p = 0; p < count; ++p)
        PointGraphEntries<T>{base, stencils, entries.data(), points.data()}(p);
    std::sort(entries.begin(), entries.end());
    entries.erase(std::unique(entries.begin(), entries.end()), entries.end());
    std::sort(points.begin(), points.end());
    NonconformalScalarGraphHost<T> next;
    if (base.rows) next.rowIds.assign(base.rowIds, base.rowIds + base.rows);
    for (auto e : entries)
        next.columns.push_back(e.column);
    for (auto p : points)
        next.pointIds.push_back(int(p.column));
    for (int r = 0; r <= base.rows; ++r)
    {
        next.offsets.push_back(
            int(std::lower_bound(entries.begin(), entries.end(), GraphEntry{r, -1}) - entries.begin()));
        next.pointOffsets.push_back(
            int(std::lower_bound(points.begin(), points.end(), GraphEntry{r, -1}) - points.begin()));
    }
    next.values.assign(entries.size(), T(0));
    for (int r = 0; r < base.rows; ++r)
        CopyBaseValues<T>{base, next.view()}(r);
    output = std::move(next);
    return true;
}

#if defined(__CUDACC__) || defined(__HIPCC__)
template<class T>
struct NonconformalScalarGraph
{
    thrust::device_vector<int64_t> rowIds, columns;
    thrust::device_vector<int> offsets, pointOffsets, pointIds;
    thrust::device_vector<T> values;
    InterfaceCsrView<T> view()
    {
        return {int(rowIds.size()),
                int(columns.size()),
                thrust::raw_pointer_cast(rowIds.data()),
                thrust::raw_pointer_cast(offsets.data()),
                thrust::raw_pointer_cast(columns.data()),
                thrust::raw_pointer_cast(values.data())};
    }
};
// Setup uses device sort/unique and binary searches; only validity/count results reach the host.
// No MPI here: the caller must handle failure collectively before entering a distributed solve.
// Backend/allocation errors throw; invalid inputs return false without replacing output.
template<class T>
bool buildNonconformalScalarGraph(InterfaceCsrView<T> base,
                                  int count,
                                  const InterfaceStencil<T>* stencils,
                                  NonconformalScalarGraph<T>& output,
                                  cudaStream_t stream = 0)
{
    using namespace nonconformal_detail;
    if (!graphSizes(base, count, stencils)) return false;
    auto policy = mars::gpuPar.on(stream);
    auto start  = thrust::make_counting_iterator(0);
    if (thrust::count_if(policy, start, start + base.rows + 1, InvalidGraphRow<T>{base}) ||
        thrust::count_if(policy, start, start + count, InvalidGraphPoint<T>{base, stencils}))
        return false;
    thrust::device_vector<GraphEntry> entries(base.rows + base.entries + 8 * count), points(count);
    thrust::for_each(policy, start, start + base.rows,
                     BaseGraphEntries<T>{base, thrust::raw_pointer_cast(entries.data())});
    thrust::for_each(policy, start, start + count,
                     PointGraphEntries<T>{base, stencils, thrust::raw_pointer_cast(entries.data()),
                                          thrust::raw_pointer_cast(points.data())});
    thrust::sort(policy, entries.begin(), entries.end());
    entries.erase(thrust::unique(policy, entries.begin(), entries.end()), entries.end());
    thrust::sort(policy, points.begin(), points.end());
    NonconformalScalarGraph<T> next;
    next.rowIds.resize(base.rows);
    if (base.rows) thrust::copy_n(policy, base.rowIds, base.rows, next.rowIds.begin());
    next.columns.resize(entries.size());
    next.pointIds.resize(count);
    thrust::transform(policy, entries.begin(), entries.end(), next.columns.begin(), GraphColumn{});
    thrust::transform(policy, points.begin(), points.end(), next.pointIds.begin(), GraphColumn{});
    next.offsets.resize(base.rows + 1);
    next.pointOffsets.resize(base.rows + 1);
    auto keys = thrust::make_transform_iterator(start, RowKey{});
    thrust::lower_bound(policy, entries.begin(), entries.end(), keys, keys + base.rows + 1, next.offsets.begin());
    thrust::lower_bound(policy, points.begin(), points.end(), keys, keys + base.rows + 1, next.pointOffsets.begin());
    next.values.assign(entries.size(), T(0));
    thrust::for_each(policy, start, start + base.rows, CopyBaseValues<T>{base, next.view()});
    // Setup buffers are temporary; completion also catches asynchronous copy/kernel failures.
    auto error = cudaStreamSynchronize(stream);
    if (error != cudaSuccess) throw std::runtime_error(cudaGetErrorString(error));
    output = std::move(next);
    return true;
}
#endif
} // namespace mars::fem
