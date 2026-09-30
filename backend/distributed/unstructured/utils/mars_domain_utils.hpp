#include "cstone/cuda/cuda_utils.cuh" // For IsDeviceVector and the memcpy helpers
#include <tuple>
#include <type_traits>
#include <algorithm>
#include <stdexcept>

namespace mars
{
#define cudaCheckError()                                                                                               \
    {                                                                                                                  \
        cudaError_t e = cudaGetLastError();                                                                            \
        if (e != cudaSuccess)                                                                                          \
        {                                                                                                              \
            printf("CUDA error %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(e));                               \
            exit(EXIT_FAILURE);                                                                                        \
        }                                                                                                              \
    }

// Add this macro near the top after includes (around line 25)
#ifdef __CUDACC__
#define MARS_HOST_DEVICE __host__ __device__
#else
#define MARS_HOST_DEVICE
#endif

// In mars_cuda_utils.hpp - Simple implementation using only cornerstone utilities
template<typename DstTuple, typename SrcTuple, std::size_t I = 0>
void copyTupleElements(DstTuple& dst, const SrcTuple& src)
{
    if constexpr (I < std::tuple_size_v<DstTuple>)
    {
        auto& dstVec       = std::get<I>(dst);
        const auto& srcVec = std::get<I>(src);

        // Resize destination to match source
        dstVec.resize(srcVec.size());

        // Determine copy direction based on vector types; the copies are complete on return
        using cstone::IsDeviceVector;
        constexpr auto stream = cstone::execution::gpuDefaultStream;
        if constexpr (IsDeviceVector<std::decay_t<decltype(srcVec)>>::value &&
                      IsDeviceVector<std::decay_t<decltype(dstVec)>>::value)
        {
            cstone::memcpyD2DAsync(stream, srcVec.data(), srcVec.size(), dstVec.data());
            cstone::syncGpu(stream);
        }
        else if constexpr (IsDeviceVector<std::decay_t<decltype(dstVec)>>::value)
        {
            cstone::memcpyH2DAsync(stream, srcVec.data(), srcVec.size(), dstVec.data());
            cstone::syncGpu(stream);
        }
        else if constexpr (IsDeviceVector<std::decay_t<decltype(srcVec)>>::value)
        {
            cstone::memcpyD2HAsync(stream, srcVec.data(), srcVec.size(), dstVec.data());
            cstone::syncGpu(stream);
        }
        else
        {
            // Host to Host
            std::copy(srcVec.begin(), srcVec.end(), dstVec.begin());
        }

        // Recursive call for next tuple element
        copyTupleElements<DstTuple, SrcTuple, I + 1>(dst, src);
    }
}

// Helper for accessing tuple elements by runtime index
template<typename Tuple, std::size_t... Is>
decltype(auto) get_tuple_element_impl(Tuple& t, std::size_t i, std::index_sequence<Is...>)
{
    decltype(auto) result = std::get<0>(t); // Default value (will be overwritten)

    bool found = false;
    (void)std::initializer_list<int>{(i == Is ? (found = true, result = std::get<Is>(t), 0) : 0)...};

    if (!found) { throw std::out_of_range("Tuple index out of range"); }

    return result;
}

// Public interface that creates the right index sequence
template<typename Tuple>
decltype(auto) getTupleElement(Tuple& t, std::size_t i)
{
    constexpr std::size_t size = std::tuple_size_v<std::remove_reference_t<Tuple>>;
    return get_tuple_element_impl(t, i, std::make_index_sequence<size>{});
}
} // namespace mars
