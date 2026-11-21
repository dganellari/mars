#pragma once

// The driver API the MARSIR loaders use, for CUDA and HIP builds. The loaders call the CUDA
// names; on AMD they map to HIP's module API, which takes an AMDGPU code object (.hsaco) where
// CUDA takes PTX. Runtime calls (cudaFree, cudaStream_t, ...) come from mars_gpu_runtime.hpp.
#include "backend/distributed/unstructured/utils/mars_gpu_runtime.hpp"

#if defined(__HIP_PLATFORM_AMD__) || defined(__HIPCC__)

#include <hip/hip_runtime.h>

using CUmodule = hipModule_t;
using CUfunction = hipFunction_t;
using CUdeviceptr = hipDeviceptr_t;
using CUstream = hipStream_t;
using CUresult = hipError_t;
#define CUDA_SUCCESS hipSuccess
#define cuInit hipInit
#define cuModuleLoadData hipModuleLoadData
#define cuModuleGetFunction hipModuleGetFunction
#define cuModuleGetGlobal hipModuleGetGlobal
#define cuModuleUnload hipModuleUnload
#define cuLaunchKernel hipModuleLaunchKernel
#define cuMemcpyDtoH hipMemcpyDtoH
#define cuMemcpyDtoDAsync hipMemcpyDtoDAsync
#define cuOccupancyMaxActiveBlocksPerMultiprocessor hipModuleOccupancyMaxActiveBlocksPerMultiprocessor
#define cuFuncGetAttribute hipFuncGetAttribute
#define CU_FUNC_ATTRIBUTE_NUM_REGS HIP_FUNC_ATTRIBUTE_NUM_REGS
#define CU_FUNC_ATTRIBUTE_SHARED_SIZE_BYTES HIP_FUNC_ATTRIBUTE_SHARED_SIZE_BYTES
inline hipError_t cuMemAlloc(hipDeviceptr_t* p, size_t bytes) { return hipMalloc(p, bytes); }
inline hipError_t cuMemFree(hipDeviceptr_t p) { return hipFree(p); }
inline hipError_t cuGetErrorString(hipError_t e, const char** s)
{
    *s = hipGetErrorString(e);
    return hipSuccess;
}
// The generated kernel files: PTX for sm_90 on NVIDIA, a gfx942 code object on AMD.
#define MARS_MARSIR_KERNEL_SUFFIX "_gfx942.hsaco"

#else

#include <cuda.h>
#include <cuda_runtime.h>
#define MARS_MARSIR_KERNEL_SUFFIX "_sm90.ptx"

#endif
