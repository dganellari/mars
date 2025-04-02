// domain_cuda_impl.hpp
#pragma once

#include <cstddef>  // For size_t definition

namespace cstone {
    // Forward declaration of Domain class
    template<class KeyType, class T, class AcceleratorType>
    class Domain;

    // Forward declaration of GpuTag
    struct GpuTag;

    // Forward declaration of DeviceVector
    template<typename T>
    class DeviceVector;
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