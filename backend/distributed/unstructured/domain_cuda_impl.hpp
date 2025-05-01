#pragma once

#include <cstddef>  // For size_t definition

// Handle both CUDA and HIP includes
#ifdef MARS_ENABLE_CUDA
#include <cuda_runtime.h>
#endif

#ifdef MARS_ENABLE_HIP
#include <hip/hip_runtime.h>
#endif

namespace cstone {
    // Forward declaration of Domain class
    template<class KeyType, class T, class AcceleratorType>
    class Domain;

    // Forward declaration of GpuTag
    struct GpuTag;

    // Forward declaration of DeviceVector
    template<typename T>
    class DeviceVector;
    
    // Forward declarations of GPU-specific classes
    template<class KeyType, class T>
    class GlobalAssignmentGpu;
}

namespace mars {
    // Forward declaration of sync implementation function
    template<typename KeyType, typename RealType>
    void syncDomainImpl(cstone::Domain<KeyType, RealType, cstone::GpuTag>* domain,
                      cstone::DeviceVector<KeyType>& elemSfcCodes,
                      cstone::DeviceVector<RealType>& elemX,
                      cstone::DeviceVector<RealType>& elemY,
                      cstone::DeviceVector<RealType>& elemZ,
                      cstone::DeviceVector<RealType>& elemH,
                      size_t elementCount);
}