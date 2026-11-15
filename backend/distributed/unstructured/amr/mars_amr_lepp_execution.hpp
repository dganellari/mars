#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <type_traits>
#include <vector>

namespace mars::amr
{
struct LeppHostExecution
{
    template<class T>
    using Vector = std::vector<T>;
    template<class V>
    static auto data(V& values)
    {
        return values.data();
    }
    template<class F>
    static void for_each(size_t n, F function)
    {
        for (size_t i = 0; i < n; ++i)
            function(i);
    }
    template<class T, class F>
    static void sort(T* values, size_t n, F less)
    {
        if (n) std::sort(values, values + n, less);
    }
    template<class T>
    static void copy(const T* source, size_t n, T* target)
    {
        if (n) std::copy_n(source, n, target);
    }
    template<class T>
    static void download(const T* source, size_t n, T* target)
    {
        copy(source, n, target);
    }
    template<class T>
    static void upload(const T* source, size_t n, T* target)
    {
        copy(source, n, target);
    }
    template<class T>
    static T read(const T* source)
    {
        return *source;
    }
    static void scan(const uint64_t* counts, size_t n, uint64_t* offsets)
    {
        if (n) std::exclusive_scan(counts, counts + n, offsets, uint64_t(0));
    }
    static uint64_t sum(const uint64_t* values, size_t n)
    {
        return n ? std::accumulate(values, values + n, uint64_t(0)) : 0;
    }
    static void synchronize() {}
};
} // namespace mars::amr

#if defined(__CUDACC__) || defined(__HIPCC__)
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/for_each.h>
#include <thrust/iterator/counting_iterator.h>
#include "backend/distributed/unstructured/utils/mars_gpu_reduce.hpp"
#include <thrust/scan.h>
#include <thrust/sort.h>
#include "cstone/cuda/cuda_utils.hpp"
#include "cstone/cuda/device_vector.h"
#include "backend/distributed/unstructured/utils/mars_domain_utils.hpp"

namespace mars::amr
{
struct LeppDeviceExecution
{
    // Keep coordinate/connectivity buffers directly movable into ElementDomain.
    template<class T>
    using Vector = std::conditional_t<std::is_arithmetic_v<T>, cstone::DeviceVector<T>, thrust::device_vector<T>>;
    template<class V>
    static auto data(V& values)
    {
        return thrust::raw_pointer_cast(values.data());
    }
    template<class F>
    static void for_each(size_t n, F function)
    {
        if (n) thrust::for_each_n(thrust::device, thrust::counting_iterator<size_t>(0), n, function);
    }
    template<class T, class F>
    static void sort(T* values, size_t n, F less)
    {
        if (n)
            thrust::sort(thrust::device, thrust::device_pointer_cast(values), thrust::device_pointer_cast(values) + n,
                         less);
    }
    template<class T>
    static void copy(const T* source, size_t n, T* target)
    {
        if (n) cudaCheckCall(cudaMemcpyAsync(target, source, n * sizeof(T), cudaMemcpyDeviceToDevice));
    }
    template<class T>
    static void download(const T* source, size_t n, T* target)
    {
        if (n) cudaCheckCall(cudaMemcpy(target, source, n * sizeof(T), cudaMemcpyDeviceToHost));
    }
    template<class T>
    static void upload(const T* source, size_t n, T* target)
    {
        if (n) cudaCheckCall(cudaMemcpy(target, source, n * sizeof(T), cudaMemcpyHostToDevice));
    }
    template<class T>
    static T read(const T* source)
    {
        T value;
        download(source, 1, &value);
        return value;
    }
    static void scan(const uint64_t* counts, size_t n, uint64_t* offsets)
    {
        if (n)
            thrust::exclusive_scan(thrust::device, thrust::device_pointer_cast(counts),
                                   thrust::device_pointer_cast(counts) + n, thrust::device_pointer_cast(offsets));
    }
    static uint64_t sum(const uint64_t* values, size_t n)
    {
        return n ? mars::reduce(thrust::device_pointer_cast(values), thrust::device_pointer_cast(values) + n,
                                uint64_t(0))
                 : 0;
    }
    static void synchronize() { cudaCheckCall(cudaStreamSynchronize(0)); }
};
} // namespace mars::amr
#endif
