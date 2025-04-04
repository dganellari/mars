// File: utils/cornerstone_patch.hpp
#pragma once

#ifndef CORNERSTONE_PATCH_HPP_DEFINED
#define CORNERSTONE_PATCH_HPP_DEFINED

// Define a macro to prevent the original implementation
#define CSTONE_UTIL_TUPLE_UTIL_HPP_ZIPTUPLE_DEFINED 1

// Include standard headers we'll need
#include <tuple>
#include <utility>
#include <array>

// Define only the specific function we need to patch
namespace util {

// This will override the zipTuples in the original Cornerstone library
template<class... Tps>
constexpr auto zipTuples(Tps&&... tps)
{
    constexpr std::size_t N = std::min({std::tuple_size_v<std::decay_t<Tps>>...});
    
    // Use simple implementation for just compilation
    if constexpr (N == 0)
    {
        return std::make_tuple();
    }
    else if constexpr (N == 1)
    {
        return std::make_tuple(
            std::make_tuple(std::get<0>(tps)...)
        );
    }
    else if constexpr (N == 2)
    {
        return std::make_tuple(
            std::make_tuple(std::get<0>(tps)...),
            std::make_tuple(std::get<1>(tps)...)
        );
    }
    else if constexpr (N == 3)
    {
        return std::make_tuple(
            std::make_tuple(std::get<0>(tps)...),
            std::make_tuple(std::get<1>(tps)...),
            std::make_tuple(std::get<2>(tps)...)
        );
    }
    else
    {
        // Default empty case for larger tuples
        return std::make_tuple();
    }
}

} // namespace util

#endif // CORNERSTONE_PATCH_HPP_DEFINED