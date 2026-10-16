#pragma once

// Thrust's stream execution policy and CUB for CUDA and HIP builds: rocThrust names the policy
// thrust::hip::par, and hipCUB provides CUB's device algorithms on AMD.
#include "mars_gpu_runtime.hpp"

#if defined(__HIP_PLATFORM_AMD__) || defined(__HIPCC__)

#include <hipcub/hipcub.hpp>
#include <thrust/system/hip/execution_policy.h>
namespace cub = hipcub;

namespace mars
{
inline constexpr auto& gpuPar = thrust::hip::par;
}

#else

#include <cub/cub.cuh>
#include <thrust/system/cuda/execution_policy.h>

namespace mars
{
inline constexpr auto& gpuPar = thrust::cuda::par;
}

#endif
