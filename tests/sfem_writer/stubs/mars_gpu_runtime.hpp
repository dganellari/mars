#pragma once
#include <cstring>
#include <vector>
inline bool fail_next_copy           = false;
inline size_t copied_field_bytes     = 0;
constexpr int cudaMemcpyDeviceToHost = 1, cudaMemcpyHostToDevice = 2, cudaSuccess = 0;
inline int cudaMemcpy(void* destination, const void* source, size_t bytes, int kind)
{
    if (fail_next_copy)
    {
        fail_next_copy = false;
        return 1;
    }
    if (kind == cudaMemcpyDeviceToHost) copied_field_bytes += bytes;
    std::memcpy(destination, source, bytes);
    return cudaSuccess;
}
namespace cstone
{
template<class T>
using DeviceVector = std::vector<T>;
namespace execution
{
struct Gpu
{
};
} // namespace execution
} // namespace cstone
