#pragma once

// Warp-level helpers for CUDA and HIP, on top of cornerstone's warp width (cstone::GpuConfig: 32 lanes
// on NVIDIA, 64 on AMD). MARS kernels work on groups of at most 32 lanes; with an explicit group width
// they compute the same on a 32-wide CUDA warp and on a 64-wide AMD wavefront. Device code only.
#include "mars_gpu_runtime.hpp"
#include "cstone/cuda/gpu_config.cuh"

namespace mars
{

//! Lanes in the groups MARS kernels are written for; a divisor of every supported warp width
inline constexpr int laneGroup = 32;
static_assert(cstone::GpuConfig::warpSize % laneGroup == 0);

//! value of srcLane within this lane's group of `width` lanes (width a power of 2, at most laneGroup)
template<class T>
__device__ __forceinline__ T shflGroup(T value, int srcLane, int width)
{
#if defined(__HIP_PLATFORM_AMD__) || defined(__HIPCC__)
    return __shfl(value, srcLane, width);
#else
    return __shfl_sync(0xFFFFFFFFu, value, srcLane, width);
#endif
}

//! Orders the shared-memory writes of a lane group before the group reads them. A wavefront runs in
//! lockstep, but the compiler may still move a read above another lane's write without a fence, so
//! unlike cstone::syncWarp this is not empty on AMD.
__device__ __forceinline__ void syncWarp()
{
#if defined(__HIP_PLATFORM_AMD__) || defined(__HIPCC__)
    __builtin_amdgcn_fence(__ATOMIC_RELEASE, "wavefront");
    __builtin_amdgcn_wave_barrier();
    __builtin_amdgcn_fence(__ATOMIC_ACQUIRE, "wavefront");
#else
    __syncwarp();
#endif
}

} // namespace mars
