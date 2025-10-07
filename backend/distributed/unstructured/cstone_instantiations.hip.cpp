#ifdef MARS_ENABLE_HIP

// Include HIP runtime first - this provides the CUDA compatibility layer
#include <hip/hip_runtime.h>
// Make sure platform is defined
#ifndef __HIP_PLATFORM_AMD__
#define __HIP_PLATFORM_AMD__ 1
#endif

namespace cstone {
    template class GlobalAssignmentGpu<unsigned int, float>;
    template class GlobalAssignmentGpu<unsigned int, double>;
    template class GlobalAssignmentGpu<uint64_t, float>;
    template class GlobalAssignmentGpu<uint64_t, double>;
}
#endif