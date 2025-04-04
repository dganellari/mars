// domain_cuda_impl.cpp
#include "domain_cuda_impl.hpp"
#include <cstone/domain/domain.hpp>  // Include complex templates here
#include <tuple>

namespace mars {
    // Implementation of syncDomainImpl function
    template<typename KeyType, typename RealType>
    void syncDomainImpl(cstone::Domain<KeyType, RealType, cstone::GpuTag>* domain,
                      cstone::DeviceVector<KeyType>& elemSfcCodes,
                      cstone::DeviceVector<RealType>& elemX,
                      cstone::DeviceVector<RealType>& elemY,
                      cstone::DeviceVector<RealType>& elemZ,
                      cstone::DeviceVector<RealType>& elemH,
                      size_t elementCount)
    {
        // Create scratch vectors required by Cornerstone
        cstone::DeviceVector<RealType> s1(elementCount);
        cstone::DeviceVector<RealType> s2(elementCount);
        cstone::DeviceVector<RealType> s3(elementCount);

        // Create a single property
        cstone::DeviceVector<RealType> d_m(elementCount, 1.0f);

        // Call sync on the domain
        domain->sync(elemSfcCodes, elemX, elemY, elemZ, elemH, std::tie(d_m), std::tie(s1, s2, s3));
    }

    // Explicit template instantiations
    template void syncDomainImpl<unsigned int, float>(
        cstone::Domain<unsigned int, float, cstone::GpuTag>* domain,
        cstone::DeviceVector<unsigned int>& elemSfcCodes,
        cstone::DeviceVector<float>& elemX,
        cstone::DeviceVector<float>& elemY,
        cstone::DeviceVector<float>& elemZ,
        cstone::DeviceVector<float>& elemH,
        size_t elementCount);

    template void syncDomainImpl<unsigned int, double>(
        cstone::Domain<unsigned int, double, cstone::GpuTag>* domain,
        cstone::DeviceVector<unsigned int>& elemSfcCodes,
        cstone::DeviceVector<double>& elemX,
        cstone::DeviceVector<double>& elemY,
        cstone::DeviceVector<double>& elemZ,
        cstone::DeviceVector<double>& elemH,
        size_t elementCount);

    template void syncDomainImpl<uint64_t, float>(
        cstone::Domain<uint64_t, float, cstone::GpuTag>* domain,
        cstone::DeviceVector<uint64_t>& elemSfcCodes,
        cstone::DeviceVector<float>& elemX,
        cstone::DeviceVector<float>& elemY,
        cstone::DeviceVector<float>& elemZ,
        cstone::DeviceVector<float>& elemH,
        size_t elementCount);

    template void syncDomainImpl<uint64_t, double>(
        cstone::Domain<uint64_t, double, cstone::GpuTag>* domain,
        cstone::DeviceVector<uint64_t>& elemSfcCodes,
        cstone::DeviceVector<double>& elemX,
        cstone::DeviceVector<double>& elemY,
        cstone::DeviceVector<double>& elemZ,
        cstone::DeviceVector<double>& elemH,
        size_t elementCount);
}