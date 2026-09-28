#pragma once
// CPU emulation of storage and kernel launch primitives for exercising the actual
// wrapper against a sequential Hypre build. This does not validate CUDA or MPI.
#include <HYPRE.h>
#include <HYPRE_IJ_mv.h>
#include <HYPRE_parcsr_ls.h>
#include <HYPRE_utilities.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
#ifndef HYPRE_SEQUENTIAL
#error "This host emulation requires a sequential CPU Hypre build"
#endif
#define __host__
#define __device__
#define __global__
using MPI_Comm = int;
constexpr int MPI_COMM_WORLD = 0, MPI_INT = 0, MPI_DOUBLE = 1, MPI_MAX = 0;
inline int host_barriers = 0;
inline int host_reductions = 0;
inline void MPI_Comm_rank(MPI_Comm, int* rank) { *rank = 0; }
inline void MPI_Allreduce(const void* src, void* dst, int n, int type, int, MPI_Comm) {
    ++host_reductions;
    std::memcpy(dst, src, n * (type == MPI_DOUBLE ? sizeof(double) : sizeof(int)));
}
inline void MPI_Barrier(MPI_Comm) { ++host_barriers; }
[[noreturn]] inline void MPI_Abort(MPI_Comm, int) { throw std::runtime_error("collective failure"); }
inline double MPI_Wtime() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
constexpr int cudaSuccess = 0, cudaMemcpyDeviceToHost = 0;
inline int cudaGetLastError() { return 0; }
inline void cudaMemcpy(void* dst, const void* src, size_t bytes, int) { std::memcpy(dst, src, bytes); }
struct HostDimension { int x = 0; };
inline HostDimension blockIdx, blockDim, threadIdx;
template<class Function> void launch_host(int grid, int block, Function function) {
    blockDim.x = block;
    for (blockIdx.x = 0; blockIdx.x < grid; ++blockIdx.x)
        for (threadIdx.x = 0; threadIdx.x < block; ++threadIdx.x) function();
}
namespace thrust {
template<class T> using device_vector = std::vector<T>;
template<class T> T* raw_pointer_cast(T* p) { return p; }
template<class T> T* device_pointer_cast(T* p) { return p; }
using std::any_of;
using std::copy;
using std::count;
using std::minmax_element;
using std::plus;
template<class T> struct maximum { T operator()(T a, T b) const { return std::max(a, b); } };
template<class T> struct minimum { T operator()(T a, T b) const { return std::min(a, b); } };
template<class In, class Out> void exclusive_scan(In first, In last, Out out) {
    std::exclusive_scan(first, last, out, 0);
}
template<class It, class T> T reduce(It first, It last, T initial) { return std::accumulate(first, last, initial); }
template<class It, class Map, class T, class Reduce>
T transform_reduce(It first, It last, Map map, T initial, Reduce reduce_op) {
    for (; first != last; ++first) initial = reduce_op(initial, map(*first));
    return initial;
}
template<class Slots, class In, class Out> void gather(Slots first, Slots last, In input, Out output) {
    for (; first != last; ++first, ++output) *output = input[*first];
}
struct CountingIterator {
    int value;
    int operator*() const { return value; }
    CountingIterator& operator++() { ++value; return *this; }
    bool operator!=(CountingIterator other) const { return value != other.value; }
};
inline CountingIterator make_counting_iterator(int value) { return {value}; }
}
namespace mars {
struct HostTestTag {};
template<class T, class Tag> struct VectorSelector { using type = std::vector<T>; };
namespace fem {
template<class Index, class Real, class Tag> class SparseMatrix {
public:
    std::vector<Index> offsets, columns;
    std::vector<Real> values;
    Index column_count = 0;
    Index numRows() const { return offsets.size() - 1; }
    Index numCols() const { return column_count; }
    Index nnz() const { return values.size(); }
    const Index* rowOffsetsPtr() const { return offsets.data(); }
    const Index* colIndicesPtr() const { return columns.data(); }
    const Real* valuesPtr() const { return values.data(); }
};
struct HypreInitGuard {
    HypreInitGuard() { HYPRE_Init(); }
    ~HypreInitGuard() { HYPRE_Finalize(); }
};
struct SolverProfile {
    enum Phase { Prepare, Setup, Solve, Finish };
    double stamp() const { return 0; }
    double lap(Phase, double) { return 0; }
};
}
}
