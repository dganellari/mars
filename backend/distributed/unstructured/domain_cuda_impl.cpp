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
#include <cstone/domain/assignment.hpp>
#include <tuple>

namespace mars
{

// Implementation of syncDomainImpl for various KeyType and RealType combinations
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
    int numRanks;
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
    // Early check for insufficient elements
    if (elementCount < numRanks)
    {
        throw std::runtime_error("Mesh has fewer elements (" + std::to_string(elementCount) + 
                               ") than MPI ranks (" + std::to_string(numRanks) +
                               "). Each rank must get at least one element for domain decomposition.");
    }

    // Create scratch buffers to match the data types being synced:
    // - First 3 for coordinates (x, y, z) -> RealType
    // - Next N for SFC connectivity keys -> KeyType (N = tuple size)
    cstone::DeviceVector<RealType> s1(elementCount);  // For x, y, z coordinates
    cstone::DeviceVector<RealType> s2(elementCount);  
    cstone::DeviceVector<RealType> s3(elementCount);  
    cstone::DeviceVector<KeyType> s4(elementCount);   // For SFC keys (larger type)
    cstone::DeviceVector<KeyType> s5(elementCount);   
    cstone::DeviceVector<KeyType> s6(elementCount);   
    cstone::DeviceVector<KeyType> s7(elementCount);
    cstone::DeviceVector<KeyType> s8(elementCount);
    cstone::DeviceVector<KeyType> s9(elementCount);
    cstone::DeviceVector<KeyType> s10(elementCount);
    cstone::DeviceVector<KeyType> s11(elementCount);

    // Create a tuple of references based on tuple size
    constexpr size_t tupleSize = std::tuple_size<SfcConnTuple>::value;
    
    if constexpr (tupleSize == 4) {
        // 4-tuple for Tet/Quad elements
        auto properties_refs = std::tie(std::get<0>(d_conn_keys_), std::get<1>(d_conn_keys_), 
                                      std::get<2>(d_conn_keys_), std::get<3>(d_conn_keys_));
        
        try {
            // Call sync with 4 properties + 7 scratch buffers
            domain->sync(elemSfcCodes, elemX, elemY, elemZ, elemH,
                       properties_refs,                                  // 4 properties of KeyType
                       std::tie(s1, s2, s3, s4, s5, s6, s7));           // Mixed scratch types
            domain->exchangeHalos(properties_refs, s4, s5);
        } catch (const std::exception& e) {
            std::string errorMsg = e.what();
            // If there's any sync error and element count is low, provide a more helpful message
            if (errorMsg.find("invalid device ordinal") != std::string::npos)
            {
                throw std::runtime_error("Domain decomposition failed. This may be due to insufficient elements (" + 
                                       std::to_string(elementCount) + ") for " + std::to_string(numRanks) +
                                       " ranks. Possibly causing a size-0 CUDA kernel launch \n" +
                                       "Original error: " + e.what());
            } else {
                // Re-throw the original exception
                throw;
            }
        }
    } else if constexpr (tupleSize == 8) {
        // 8-tuple for Hex elements
        auto properties_refs = std::tie(std::get<0>(d_conn_keys_), std::get<1>(d_conn_keys_), 
                                      std::get<2>(d_conn_keys_), std::get<3>(d_conn_keys_),
                                      std::get<4>(d_conn_keys_), std::get<5>(d_conn_keys_),
                                      std::get<6>(d_conn_keys_), std::get<7>(d_conn_keys_));

        try {
            // Call sync with 8 properties + 11 scratch buffers
            domain->sync(elemSfcCodes, elemX, elemY, elemZ, elemH,
                       properties_refs,                                              // 8 properties of KeyType
                       std::tie(s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11));     // Mixed scratch types
            domain->exchangeHalos(properties_refs, s4, s5);
        } catch (const std::exception& e) {
            std::string errorMsg = e.what();
            // If there's any sync error and element count is low, provide a more helpful message
            if (errorMsg.find("invalid device ordinal") != std::string::npos)
            {
                throw std::runtime_error("Domain decomposition failed. This may be due to insufficient elements (" + 
                                       std::to_string(elementCount) + ") for " + std::to_string(numRanks) +
                                       " ranks. Possibly causing a size-0 CUDA kernel launch \n" +
                                       "Original error: " + e.what());
            } else {
                // Re-throw the original exception
                throw;
            }
        }
    } else {
        throw std::runtime_error("Unsupported tuple size: " + std::to_string(tupleSize));
    }
                              

    // Update element count after sync
    elementCount = domain->nParticlesWithHalos();
}

// Explicit template instantiations for syncDomainWithSfcConn
// For elements with unsigned int keys and float coordinates
template void syncDomainImpl<unsigned int,
                             float,
                             std::tuple<cstone::DeviceVector<unsigned int>,
                                        cstone::DeviceVector<unsigned int>,
                                        cstone::DeviceVector<unsigned int>,
                                        cstone::DeviceVector<unsigned int>,
                                        cstone::DeviceVector<unsigned int>,
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
               cstone::DeviceVector<unsigned int>,
               cstone::DeviceVector<unsigned int>,
               cstone::DeviceVector<unsigned int>,
               cstone::DeviceVector<unsigned int>,
               cstone::DeviceVector<unsigned int>>& d_conn_keys_);

// For elements with unsigned int keys and double coordinates
template void syncDomainImpl<unsigned int,
                             double,
                             std::tuple<cstone::DeviceVector<unsigned int>,
                                        cstone::DeviceVector<unsigned int>,
                                        cstone::DeviceVector<unsigned int>,
                                        cstone::DeviceVector<unsigned int>,
                                        cstone::DeviceVector<unsigned int>,
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
               cstone::DeviceVector<unsigned int>,
               cstone::DeviceVector<unsigned int>,
               cstone::DeviceVector<unsigned int>,
               cstone::DeviceVector<unsigned int>,
               cstone::DeviceVector<unsigned int>>& d_conn_keys_);

// For elements with uint64_t keys and float coordinates
template void
syncDomainImpl<uint64_t,
               float,
               std::tuple<cstone::DeviceVector<uint64_t>,
                          cstone::DeviceVector<uint64_t>,
                          cstone::DeviceVector<uint64_t>,
                          cstone::DeviceVector<uint64_t>,
                          cstone::DeviceVector<uint64_t>,
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
                                                                      cstone::DeviceVector<uint64_t>,
                                                                      cstone::DeviceVector<uint64_t>,
                                                                      cstone::DeviceVector<uint64_t>,
                                                                      cstone::DeviceVector<uint64_t>,
                                                                      cstone::DeviceVector<uint64_t>>& d_conn_keys_);

// For elements with uint64_t keys and double coordinates
template void
syncDomainImpl<uint64_t,
               double,
               std::tuple<cstone::DeviceVector<uint64_t>,
                          cstone::DeviceVector<uint64_t>,
                          cstone::DeviceVector<uint64_t>,
                          cstone::DeviceVector<uint64_t>,
                          cstone::DeviceVector<uint64_t>,
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
                                                                      cstone::DeviceVector<uint64_t>,
                                                                      cstone::DeviceVector<uint64_t>,
                                                                      cstone::DeviceVector<uint64_t>,
                                                                      cstone::DeviceVector<uint64_t>,
                                                                      cstone::DeviceVector<uint64_t>>& d_conn_keys_);

// Add 4-tuple instantiations for TetTag/QuadTag (4 nodes per element)
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

template void syncDomainImpl<uint64_t,
                             float,
                             std::tuple<cstone::DeviceVector<uint64_t>,
                                        cstone::DeviceVector<uint64_t>,
                                        cstone::DeviceVector<uint64_t>,
                                        cstone::DeviceVector<uint64_t>>>(
    cstone::Domain<uint64_t, float, cstone::GpuTag>* domain,
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

template void syncDomainImpl<uint64_t,
                             double,
                             std::tuple<cstone::DeviceVector<uint64_t>,
                                        cstone::DeviceVector<uint64_t>,
                                        cstone::DeviceVector<uint64_t>,
                                        cstone::DeviceVector<uint64_t>>>(
    cstone::Domain<uint64_t, double, cstone::GpuTag>* domain,
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
