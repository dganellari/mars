#pragma once

#include <algorithm>
#include <cstddef>
#include <tuple>
#include <utility>
#include <vector>

#if defined(__CUDACC__) || defined(__HIPCC__)
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/copy.h>
#endif

namespace mars
{

template<int N, typename KeyType = unsigned>
auto createNVectors(size_t size)
{
    if constexpr (N == 3)
    {
        return std::tuple<std::vector<KeyType>, std::vector<KeyType>, std::vector<KeyType>>(
            std::vector<KeyType>(size), std::vector<KeyType>(size), std::vector<KeyType>(size));
    }
    else if constexpr (N == 4)
    {
        return std::tuple<std::vector<KeyType>, std::vector<KeyType>, std::vector<KeyType>, std::vector<KeyType>>(
            std::vector<KeyType>(size), std::vector<KeyType>(size), std::vector<KeyType>(size), std::vector<KeyType>(size));
    }
    else if constexpr (N == 8)
    {
        return std::tuple<std::vector<KeyType>, std::vector<KeyType>, std::vector<KeyType>, std::vector<KeyType>,
                          std::vector<KeyType>, std::vector<KeyType>, std::vector<KeyType>, std::vector<KeyType>>(
            std::vector<KeyType>(size), std::vector<KeyType>(size), std::vector<KeyType>(size), std::vector<KeyType>(size),
            std::vector<KeyType>(size), std::vector<KeyType>(size), std::vector<KeyType>(size), std::vector<KeyType>(size));
    }
}

// Large local node lists use the GPU; small lists avoid the transfer overhead.
template<typename KeyType>
inline std::vector<KeyType> dedupNodes(std::vector<KeyType>& allNodes)
{
#if defined(__CUDACC__) || defined(__HIPCC__)
    constexpr size_t cpuThreshold = 1ull << 20;  // ~1M keys: H2D round-trip not worth it
    if (allNodes.size() >= cpuThreshold)
    {
        thrust::device_vector<KeyType> d_keys(allNodes.size());
        thrust::copy(allNodes.begin(), allNodes.end(), d_keys.begin());
        thrust::sort(d_keys.begin(), d_keys.end());
        auto newEnd = thrust::unique(d_keys.begin(), d_keys.end());
        size_t uniqueCount = newEnd - d_keys.begin();

        std::vector<KeyType> result(uniqueCount);
        thrust::copy(d_keys.begin(), d_keys.begin() + uniqueCount, result.begin());
        return result;
    }
#endif
    std::sort(allNodes.begin(), allNodes.end());
    auto last = std::unique(allNodes.begin(), allNodes.end());
    allNodes.erase(last, allNodes.end());
    return std::move(allNodes);
}

} // namespace mars
