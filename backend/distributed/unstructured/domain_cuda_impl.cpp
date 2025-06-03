// domain_cuda_impl.cpp
#include "domain_cuda_impl.hpp"

// GPU-specific includes
#ifdef MARS_ENABLE_CUDA
#include <cuda_runtime.h>
#include "cstone/cuda/cuda_utils.hpp"
#endif

#ifdef MARS_ENABLE_HIP
#include <hip/hip_runtime.h>
// Include both for compatibility - HIP can use CUDA headers too
#include "cstone/cuda/cuda_utils.hpp"
#endif

#include <cstone/domain/domain.hpp>
#include <cstone/domain/assignment_gpu.cuh>
#include <cstone/domain/assignment.hpp>
#include <tuple>

namespace mars
{

template<typename KeyType, typename RealType, typename SfcConnTuple>
void syncDomainImpl(cstone::Domain<KeyType, RealType, cstone::GpuTag>* domain,
                    cstone::DeviceVector<KeyType>& elemSfcCodes,
                    cstone::DeviceVector<RealType>& elemX,
                    cstone::DeviceVector<RealType>& elemY,
                    cstone::DeviceVector<RealType>& elemZ,
                    cstone::DeviceVector<RealType>& elemH,
                    size_t& elementCount,
                    SfcConnTuple& d_conn_keys_)
{
    // Create scratch buffers to match the data types being synced:
    // - First 3 for coordinates (x, y, z) -> RealType
    // - Next 4+ for SFC connectivity keys -> KeyType
    cstone::DeviceVector<RealType> s1(elementCount);  // For x, y, z coordinates
    cstone::DeviceVector<RealType> s2(elementCount);  
    cstone::DeviceVector<RealType> s3(elementCount);  
    cstone::DeviceVector<KeyType> s4(elementCount);   // For SFC keys (larger type)
    cstone::DeviceVector<KeyType> s5(elementCount);   
    cstone::DeviceVector<KeyType> s6(elementCount);   
    cstone::DeviceVector<KeyType> s7(elementCount);   

    // Create a tuple of references to match the expected signature
    auto properties_refs = std::tie(std::get<0>(d_conn_keys_), std::get<1>(d_conn_keys_), 
                                    std::get<2>(d_conn_keys_), std::get<3>(d_conn_keys_));

    // Call sync with mixed-type scratch buffers
    domain->sync(elemSfcCodes, elemX, elemY, elemZ, elemH,
                 properties_refs,                                    // 4 properties of KeyType
                 std::tie(s1, s2, s3, s4, s5, s6, s7));             // Mixed scratch types

    // Update element count after sync
    elementCount = domain->endIndex() - domain->startIndex();
}

// Explicit template instantiations for syncDomainWithSfcConn
// For tet elements with unsigned int keys and float coordinates
template void syncDomainImpl<unsigned int,
                             float,
                             std::tuple<cstone::DeviceVector<unsigned int>,
                                        cstone::DeviceVector<unsigned int>,
                                        cstone::DeviceVector<unsigned int>,
                                        cstone::DeviceVector<unsigned int>>>(
    cstone::Domain<unsigned int, float, cstone::GpuTag>* domain,
    cstone::DeviceVector<unsigned int>& elemSfcCodes,
    cstone::DeviceVector<float>& elemX,
    cstone::DeviceVector<float>& elemY,
    cstone::DeviceVector<float>& elemZ,
    cstone::DeviceVector<float>& elemH,
    size_t& elementCount,
    std::tuple<cstone::DeviceVector<unsigned int>,
               cstone::DeviceVector<unsigned int>,
               cstone::DeviceVector<unsigned int>,
               cstone::DeviceVector<unsigned int>>& d_conn_keys_);

// For tet elements with unsigned int keys and double coordinates
template void syncDomainImpl<unsigned int,
                             double,
                             std::tuple<cstone::DeviceVector<unsigned int>,
                                        cstone::DeviceVector<unsigned int>,
                                        cstone::DeviceVector<unsigned int>,
                                        cstone::DeviceVector<unsigned int>>>(
    cstone::Domain<unsigned int, double, cstone::GpuTag>* domain,
    cstone::DeviceVector<unsigned int>& elemSfcCodes,
    cstone::DeviceVector<double>& elemX,
    cstone::DeviceVector<double>& elemY,
    cstone::DeviceVector<double>& elemZ,
    cstone::DeviceVector<double>& elemH,
    size_t& elementCount,
    std::tuple<cstone::DeviceVector<unsigned int>,
               cstone::DeviceVector<unsigned int>,
               cstone::DeviceVector<unsigned int>,
               cstone::DeviceVector<unsigned int>>& d_conn_keys_);

// For tet elements with uint64_t keys and float coordinates
template void
syncDomainImpl<uint64_t,
               float,
               std::tuple<cstone::DeviceVector<uint64_t>,
                          cstone::DeviceVector<uint64_t>,
                          cstone::DeviceVector<uint64_t>,
                          cstone::DeviceVector<uint64_t>>>(cstone::Domain<uint64_t, float, cstone::GpuTag>* domain,
                                                           cstone::DeviceVector<uint64_t>& elemSfcCodes,
                                                           cstone::DeviceVector<float>& elemX,
                                                           cstone::DeviceVector<float>& elemY,
                                                           cstone::DeviceVector<float>& elemZ,
                                                           cstone::DeviceVector<float>& elemH,
                                                           size_t& elementCount,
                                                           std::tuple<cstone::DeviceVector<uint64_t>,
                                                                      cstone::DeviceVector<uint64_t>,
                                                                      cstone::DeviceVector<uint64_t>,
                                                                      cstone::DeviceVector<uint64_t>>& d_conn_keys_);

// For tet elements with uint64_t keys and double coordinates
template void
syncDomainImpl<uint64_t,
               double,
               std::tuple<cstone::DeviceVector<uint64_t>,
                          cstone::DeviceVector<uint64_t>,
                          cstone::DeviceVector<uint64_t>,
                          cstone::DeviceVector<uint64_t>>>(cstone::Domain<uint64_t, double, cstone::GpuTag>* domain,
                                                           cstone::DeviceVector<uint64_t>& elemSfcCodes,
                                                           cstone::DeviceVector<double>& elemX,
                                                           cstone::DeviceVector<double>& elemY,
                                                           cstone::DeviceVector<double>& elemZ,
                                                           cstone::DeviceVector<double>& elemH,
                                                           size_t& elementCount,
                                                           std::tuple<cstone::DeviceVector<uint64_t>,
                                                                      cstone::DeviceVector<uint64_t>,
                                                                      cstone::DeviceVector<uint64_t>,
                                                                      cstone::DeviceVector<uint64_t>>& d_conn_keys_);

} // namespace mars
