#pragma once

// Device reductions with the arguments of thrust::reduce, transform_reduce, count, count_if and any_of,
// minus the execution policy.
//
// On NVIDIA these are the Thrust calls themselves. On AMD every Thrust reduction becomes about ten
// rocPRIM kernels for each of 13 GPU architectures, and the Navier-Stokes drivers have more than a
// hundred of them, so one driver took over an hour to compile. There a call site adds one kernel instead:
// each block folds its part into a partial, and a shared second kernel folds the partials in a fixed
// order. -DMARS_HIP_ROCPRIM_REDUCTIONS=ON goes back to rocPRIM on AMD.
#include "mars_gpu_runtime.hpp"
#include <thrust/functional.h>
#include <cstddef>
#include <iterator>

#if (defined(__HIP_PLATFORM_AMD__) || defined(__HIPCC__)) && !defined(MARS_HIP_ROCPRIM_REDUCTIONS)
#define MARS_OWN_REDUCTIONS 1
#endif

#ifdef MARS_OWN_REDUCTIONS

#include "mars_domain_utils.hpp"
#include <algorithm>
#include <new>

namespace mars
{
namespace detail
{

constexpr int reduceBlockSize = 256;

// Per-process scratch for the partials and the result, grown on demand and never freed: the GPU context
// goes away with the process anyway. Reductions run on the default stream, one after the other.
inline char* reduceScratch(std::size_t bytes)
{
    static char* buffer         = nullptr;
    static std::size_t capacity = 0;
    if (bytes > capacity)
    {
        if (buffer) cudaCheckCall(cudaFree(buffer));
        cudaCheckCall(cudaMalloc(&buffer, bytes));
        capacity = bytes;
    }
    return buffer;
}

inline int reduceMaxBlocks()
{
    static int maxBlocks = 0;
    if (maxBlocks == 0)
    {
        int device = 0, sms = 0;
        cudaCheckCall(cudaGetDevice(&device));
        cudaCheckCall(cudaDeviceGetAttribute(&sms, cudaDevAttrMultiProcessorCount, device));
        maxBlocks = 4 * sms;
    }
    return maxBlocks;
}

// Folds s[0..reduceBlockSize) into s[0]. h[i] says whether s[i] holds a value: the operator needs no
// identity element (there is none for thrust::maximum that is safe for every input).
// Shared memory for reduceBlockSize values of T, left unconstructed: a __shared__ array of a type with a
// constructor (thrust::tuple) is not allowed. The values are placed with placement new.
template<class T>
struct alignas(T) SharedValues
{
    unsigned char bytes[sizeof(T) * reduceBlockSize];
    __device__ T* data() { return reinterpret_cast<T*>(bytes); }
};

template<class T, class Op>
__device__ void blockFold(T* s, unsigned char* h, Op op)
{
    for (int width = reduceBlockSize / 2; width > 0; width /= 2)
    {
        const int t = threadIdx.x;
        if (t < width && h[t + width])
        {
            s[t] = h[t] ? op(s[t], s[t + width]) : s[t + width];
            h[t] = 1;
        }
        __syncthreads();
    }
}

template<class T, class It, class F, class Op>
__global__ void __launch_bounds__(reduceBlockSize)
    transformReducePartials(It first, std::size_t n, F f, Op op, T* partial, unsigned char* partialHas)
{
    __shared__ SharedValues<T> values;
    __shared__ unsigned char h[reduceBlockSize];
    T* s = values.data();

    T acc{};
    bool has               = false;
    const std::size_t step = std::size_t(gridDim.x) * reduceBlockSize;
    for (std::size_t i = std::size_t(blockIdx.x) * reduceBlockSize + threadIdx.x; i < n; i += step)
    {
        const T v = f(first[i]);
        acc       = has ? op(acc, v) : v;
        has       = true;
    }
    new (s + threadIdx.x) T(acc);
    h[threadIdx.x] = has;
    __syncthreads();
    blockFold(s, h, op);
    if (threadIdx.x == 0)
    {
        partial[blockIdx.x]    = s[0];
        partialHas[blockIdx.x] = h[0];
    }
}

// One block: result = op(init, partial[0], ..., partial[numPartials - 1]).
template<class T, class Op>
__global__ void __launch_bounds__(reduceBlockSize)
    foldPartials(const T* partial, const unsigned char* partialHas, int numPartials, T init, Op op, T* result)
{
    __shared__ SharedValues<T> values;
    __shared__ unsigned char h[reduceBlockSize];
    T* s = values.data();

    T acc{};
    bool has = false;
    for (int b = threadIdx.x; b < numPartials; b += reduceBlockSize)
    {
        if (!partialHas[b]) continue;
        acc = has ? op(acc, partial[b]) : partial[b];
        has = true;
    }
    new (s + threadIdx.x) T(acc);
    h[threadIdx.x] = has;
    __syncthreads();
    blockFold(s, h, op);
    if (threadIdx.x == 0) *result = h[0] ? op(init, s[0]) : init;
}

template<class T>
struct ConvertTo
{
    template<class V>
    __host__ __device__ T operator()(const V& v) const
    {
        return T(v);
    }
};

template<class D, class V>
struct CountEqual
{
    V value;
    template<class U>
    __host__ __device__ D operator()(const U& u) const
    {
        return u == value ? D(1) : D(0);
    }
};

// device-only: the predicates are often __device__ lambdas
template<class D, class Pred>
struct CountTrue
{
    Pred pred;
    template<class U>
    __device__ D operator()(const U& u) const
    {
        return pred(u) ? D(1) : D(0);
    }
};

} // namespace detail

// op(init, f(first[0]), ..., f(first[n-1])) on the device, for any random-access device iterator: raw and
// thrust::device_ptr pointers, device_vector iterators, counting iterators. Returns to the host, like Thrust.
template<class It, class F, class T, class Op>
T transformReduce(It first, It last, F f, T init, Op op)
{
    const std::size_t n = std::size_t(last - first);
    if (n == 0) return init;

    const int blocks = int(std::min<std::size_t>((n + detail::reduceBlockSize - 1) / detail::reduceBlockSize,
                                                 std::size_t(detail::reduceMaxBlocks())));
    // result, partials, their flags; each part 16-byte aligned
    const std::size_t slot         = (sizeof(T) + 15) / 16 * 16;
    const std::size_t partialBytes = (blocks * sizeof(T) + 15) / 16 * 16;
    char* scratch                  = detail::reduceScratch(slot + partialBytes + blocks);
    T* result                      = reinterpret_cast<T*>(scratch);
    T* partial                     = reinterpret_cast<T*>(scratch + slot);
    unsigned char* hasValue        = reinterpret_cast<unsigned char*>(scratch + slot + partialBytes);

    detail::transformReducePartials<T><<<blocks, detail::reduceBlockSize>>>(first, n, f, op, partial, hasValue);
    detail::foldPartials<T><<<1, detail::reduceBlockSize>>>(partial, hasValue, blocks, init, op, result);
    cudaCheckCall(cudaGetLastError());
    T value;
    cudaCheckCall(cudaMemcpy(&value, result, sizeof(T), cudaMemcpyDeviceToHost));
    return value;
}

template<class It, class T, class Op>
T reduce(It first, It last, T init, Op op)
{
    return mars::transformReduce(first, last, detail::ConvertTo<T>{}, init, op);
}

template<class It, class V>
typename std::iterator_traits<It>::difference_type count(It first, It last, const V& value)
{
    using D = typename std::iterator_traits<It>::difference_type;
    return mars::transformReduce(first, last, detail::CountEqual<D, V>{value}, D(0), thrust::plus<D>());
}

template<class It, class Pred>
typename std::iterator_traits<It>::difference_type countIf(It first, It last, Pred pred)
{
    using D = typename std::iterator_traits<It>::difference_type;
    return mars::transformReduce(first, last, detail::CountTrue<D, Pred>{pred}, D(0), thrust::plus<D>());
}

template<class It, class Pred>
bool anyOf(It first, It last, Pred pred)
{
    return mars::countIf(first, last, pred) > 0;
}

} // namespace mars

#else // Thrust: on NVIDIA, and on AMD with MARS_HIP_ROCPRIM_REDUCTIONS

#include <thrust/count.h>
#include <thrust/execution_policy.h>
#include <thrust/logical.h>
#include <thrust/reduce.h>
#include <thrust/transform_reduce.h>

namespace mars
{

template<class It, class F, class T, class Op>
T transformReduce(It first, It last, F f, T init, Op op)
{
    return thrust::transform_reduce(thrust::device, first, last, f, init, op);
}

template<class It, class T, class Op>
T reduce(It first, It last, T init, Op op)
{
    return thrust::reduce(thrust::device, first, last, init, op);
}

template<class It, class V>
typename std::iterator_traits<It>::difference_type count(It first, It last, const V& value)
{
    return thrust::count(thrust::device, first, last, value);
}

template<class It, class Pred>
typename std::iterator_traits<It>::difference_type countIf(It first, It last, Pred pred)
{
    return thrust::count_if(thrust::device, first, last, pred);
}

template<class It, class Pred>
bool anyOf(It first, It last, Pred pred)
{
    return thrust::any_of(thrust::device, first, last, pred);
}

} // namespace mars

#endif

namespace mars
{

// mars:: qualified: an unqualified call would also find thrust::reduce through the argument types
template<class It, class T>
T reduce(It first, It last, T init)
{
    return mars::reduce(first, last, init, thrust::plus<T>());
}

template<class It>
typename std::iterator_traits<It>::value_type reduce(It first, It last)
{
    using T = typename std::iterator_traits<It>::value_type;
    return mars::reduce(first, last, T{}, thrust::plus<T>());
}

} // namespace mars
