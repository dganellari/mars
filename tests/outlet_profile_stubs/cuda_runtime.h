#pragma once

// Test-only include directory: never add it to a production target.
inline int outlet_profile_test_syncs = 0;
using cudaError_t = int;
inline constexpr cudaError_t cudaSuccess = 0;
inline cudaError_t cudaDeviceSynchronize()
{
    ++outlet_profile_test_syncs;
    return cudaSuccess;
}
inline const char* cudaGetErrorString(cudaError_t) { return "stubbed CUDA"; }
