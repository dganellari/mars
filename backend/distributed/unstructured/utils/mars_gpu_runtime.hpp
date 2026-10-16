#pragma once

// The GPU runtime for CUDA and HIP builds. cornerstone maps the common CUDA runtime calls to HIP
// (cstone/cuda/cuda_runtime.hpp); this adds the ones MARS uses beyond its list, so MARS code calls the
// CUDA names on both platforms. Host-safe: no device code here.
#include "cstone/cuda/cuda_runtime.hpp"

#if defined(__HIP_PLATFORM_AMD__) || defined(__HIPCC__)

#define cudaDevAttrMaxSharedMemoryPerBlockOptin hipDeviceAttributeSharedMemPerBlockOptin
#define cudaDeviceGetAttribute hipDeviceGetAttribute
#define cudaErrorIllegalAddress hipErrorIllegalAddress
#define cudaErrorInvalidConfiguration hipErrorInvalidConfiguration
#define cudaErrorInvalidDevice hipErrorInvalidDevice
#define cudaFuncAttributeMaxDynamicSharedMemorySize hipFuncAttributeMaxDynamicSharedMemorySize
#define cudaFuncAttributePreferredSharedMemoryCarveout hipFuncAttributePreferredSharedMemoryCarveout
#define cudaFuncSetAttribute ::mars::funcSetAttribute
#define cudaGraphDestroy hipGraphDestroy
#define cudaGraphExecDestroy hipGraphExecDestroy
#define cudaGraphExec_t hipGraphExec_t
#define cudaGraphInstantiate hipGraphInstantiate
#define cudaGraphInstantiateWithFlags hipGraphInstantiateWithFlags
#define cudaGraphLaunch hipGraphLaunch
#define cudaGraph_t hipGraph_t
#define cudaMallocAsync hipMallocAsync
#define cudaMemcpyAsync hipMemcpyAsync
#define cudaMemsetAsync hipMemsetAsync
#define cudaStreamBeginCapture hipStreamBeginCapture
#define cudaStreamCaptureModeThreadLocal hipStreamCaptureModeThreadLocal
#define cudaStreamEndCapture hipStreamEndCapture

namespace mars
{
// HIP takes the kernel only as const void*; CUDA also has a template taking the kernel pointer itself
template<class Kernel>
inline hipError_t funcSetAttribute(Kernel* kernel, hipFuncAttribute attr, int value)
{
    return hipFuncSetAttribute(reinterpret_cast<const void*>(kernel), attr, value);
}
} // namespace mars

#endif
