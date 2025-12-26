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

// Implementation of syncDomainImplWithOrigCoords for syncing with original coordinates
template<typename KeyType, typename RealType, typename SfcConnTuple, typename OrigCoordsTuple>
void syncDomainImplWithOrigCoords(cstone::Domain<KeyType, RealType, cstone::GpuTag>* domain,
                                   cstone::DeviceVector<KeyType>& elemSfcCodes,
                                   cstone::DeviceVector<RealType>& elemX,
                                   cstone::DeviceVector<RealType>& elemY,
                                   cstone::DeviceVector<RealType>& elemZ,
                                   cstone::DeviceVector<RealType>& elemH,
                                   size_t& elementCount,
                                   SfcConnTuple& d_conn_keys_,
                                   OrigCoordsTuple& d_orig_coords_)
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

    // Create scratch buffers - need MORE for original coordinates
    // We need: 3 for x,y,z coords + 8 for SFC keys + 24 for original coords = 35 total
    cstone::DeviceVector<RealType> s1(elementCount);   // For x, y, z coordinates
    cstone::DeviceVector<RealType> s2(elementCount);
    cstone::DeviceVector<RealType> s3(elementCount);
    cstone::DeviceVector<KeyType> s4(elementCount);    // For SFC keys
    cstone::DeviceVector<KeyType> s5(elementCount);
    cstone::DeviceVector<KeyType> s6(elementCount);
    cstone::DeviceVector<KeyType> s7(elementCount);
    cstone::DeviceVector<KeyType> s8(elementCount);
    cstone::DeviceVector<KeyType> s9(elementCount);
    cstone::DeviceVector<KeyType> s10(elementCount);
    cstone::DeviceVector<KeyType> s11(elementCount);
    // Additional scratch buffers for original coordinates (24 RealType properties)
    cstone::DeviceVector<RealType> s12(elementCount);
    cstone::DeviceVector<RealType> s13(elementCount);
    cstone::DeviceVector<RealType> s14(elementCount);
    cstone::DeviceVector<RealType> s15(elementCount);
    cstone::DeviceVector<RealType> s16(elementCount);
    cstone::DeviceVector<RealType> s17(elementCount);
    cstone::DeviceVector<RealType> s18(elementCount);
    cstone::DeviceVector<RealType> s19(elementCount);
    cstone::DeviceVector<RealType> s20(elementCount);
    cstone::DeviceVector<RealType> s21(elementCount);
    cstone::DeviceVector<RealType> s22(elementCount);
    cstone::DeviceVector<RealType> s23(elementCount);
    cstone::DeviceVector<RealType> s24(elementCount);
    cstone::DeviceVector<RealType> s25(elementCount);
    cstone::DeviceVector<RealType> s26(elementCount);
    cstone::DeviceVector<RealType> s27(elementCount);
    cstone::DeviceVector<RealType> s28(elementCount);
    cstone::DeviceVector<RealType> s29(elementCount);
    cstone::DeviceVector<RealType> s30(elementCount);
    cstone::DeviceVector<RealType> s31(elementCount);
    cstone::DeviceVector<RealType> s32(elementCount);
    cstone::DeviceVector<RealType> s33(elementCount);
    cstone::DeviceVector<RealType> s34(elementCount);
    cstone::DeviceVector<RealType> s35(elementCount);

    constexpr size_t sfcTupleSize = std::tuple_size<SfcConnTuple>::value;
    constexpr size_t origCoordsTupleSize = std::tuple_size<OrigCoordsTuple>::value;

    // Original coordinates should be 24 properties (8 nodes × 3 coords)
    static_assert(origCoordsTupleSize == 24, "OrigCoordsTuple must have 24 elements for hex8");

    if constexpr (sfcTupleSize == 8) {
        // 8-tuple for Hex elements + 24 original coordinate properties
        auto sfc_refs = std::tie(std::get<0>(d_conn_keys_), std::get<1>(d_conn_keys_),
                                std::get<2>(d_conn_keys_), std::get<3>(d_conn_keys_),
                                std::get<4>(d_conn_keys_), std::get<5>(d_conn_keys_),
                                std::get<6>(d_conn_keys_), std::get<7>(d_conn_keys_));

        auto orig_refs = std::tie(std::get<0>(d_orig_coords_), std::get<1>(d_orig_coords_),
                                 std::get<2>(d_orig_coords_), std::get<3>(d_orig_coords_),
                                 std::get<4>(d_orig_coords_), std::get<5>(d_orig_coords_),
                                 std::get<6>(d_orig_coords_), std::get<7>(d_orig_coords_),
                                 std::get<8>(d_orig_coords_), std::get<9>(d_orig_coords_),
                                 std::get<10>(d_orig_coords_), std::get<11>(d_orig_coords_),
                                 std::get<12>(d_orig_coords_), std::get<13>(d_orig_coords_),
                                 std::get<14>(d_orig_coords_), std::get<15>(d_orig_coords_),
                                 std::get<16>(d_orig_coords_), std::get<17>(d_orig_coords_),
                                 std::get<18>(d_orig_coords_), std::get<19>(d_orig_coords_),
                                 std::get<20>(d_orig_coords_), std::get<21>(d_orig_coords_),
                                 std::get<22>(d_orig_coords_), std::get<23>(d_orig_coords_));

        // Combine SFC connectivity + original coords into one tuple
        auto all_properties = std::tuple_cat(sfc_refs, orig_refs);

        try {
            // Call sync with 8 SFC properties + 24 original coordinate properties = 32 total
            domain->sync(elemSfcCodes, elemX, elemY, elemZ, elemH,
                       all_properties,  // 32 properties total
                       std::tie(s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11,
                               s12, s13, s14, s15, s16, s17, s18, s19, s20, s21, s22, s23,
                               s24, s25, s26, s27, s28, s29, s30, s31, s32, s33, s34, s35));
            domain->exchangeHalos(all_properties, s4, s5);
        } catch (const std::exception& e) {
            std::string errorMsg = e.what();
            if (errorMsg.find("invalid device ordinal") != std::string::npos)
            {
                throw std::runtime_error("Domain decomposition failed with original coordinates. " +
                                       std::string("Element count: ") + std::to_string(elementCount) +
                                       ", ranks: " + std::to_string(numRanks) +
                                       "\nOriginal error: " + e.what());
            } else {
                throw;
            }
        }
    } else {
        throw std::runtime_error("syncDomainImplWithOrigCoords only supports 8-tuple SFC connectivity (hex8)");
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

// Explicit template instantiations for syncDomainImplWithOrigCoords
// For hex8 elements with original coordinates (8 SFC keys + 24 coordinate arrays)

// Helper type alias for 24-tuple of original coordinates (references, as created by std::tie)
template<typename T>
using OrigCoords24 = std::tuple<
    cstone::DeviceVector<T>&, cstone::DeviceVector<T>&, cstone::DeviceVector<T>&,
    cstone::DeviceVector<T>&, cstone::DeviceVector<T>&, cstone::DeviceVector<T>&,
    cstone::DeviceVector<T>&, cstone::DeviceVector<T>&, cstone::DeviceVector<T>&,
    cstone::DeviceVector<T>&, cstone::DeviceVector<T>&, cstone::DeviceVector<T>&,
    cstone::DeviceVector<T>&, cstone::DeviceVector<T>&, cstone::DeviceVector<T>&,
    cstone::DeviceVector<T>&, cstone::DeviceVector<T>&, cstone::DeviceVector<T>&,
    cstone::DeviceVector<T>&, cstone::DeviceVector<T>&, cstone::DeviceVector<T>&,
    cstone::DeviceVector<T>&, cstone::DeviceVector<T>&, cstone::DeviceVector<T>&
>;

// unsigned int keys, float coords
template void syncDomainImplWithOrigCoords<
    unsigned int, float,
    std::tuple<cstone::DeviceVector<unsigned int>, cstone::DeviceVector<unsigned int>,
               cstone::DeviceVector<unsigned int>, cstone::DeviceVector<unsigned int>,
               cstone::DeviceVector<unsigned int>, cstone::DeviceVector<unsigned int>,
               cstone::DeviceVector<unsigned int>, cstone::DeviceVector<unsigned int>>,
    OrigCoords24<float>>(
    cstone::Domain<unsigned int, float, cstone::GpuTag>* domain,
    cstone::DeviceVector<unsigned int>& elemSfcCodes,
    cstone::DeviceVector<float>& elemX,
    cstone::DeviceVector<float>& elemY,
    cstone::DeviceVector<float>& elemZ,
    cstone::DeviceVector<float>& elemH,
    size_t& elementCount,
    std::tuple<cstone::DeviceVector<unsigned int>, cstone::DeviceVector<unsigned int>,
               cstone::DeviceVector<unsigned int>, cstone::DeviceVector<unsigned int>,
               cstone::DeviceVector<unsigned int>, cstone::DeviceVector<unsigned int>,
               cstone::DeviceVector<unsigned int>, cstone::DeviceVector<unsigned int>>& d_conn_keys_,
    OrigCoords24<float>& d_orig_coords_);

// unsigned int keys, double coords
template void syncDomainImplWithOrigCoords<
    unsigned int, double,
    std::tuple<cstone::DeviceVector<unsigned int>, cstone::DeviceVector<unsigned int>,
               cstone::DeviceVector<unsigned int>, cstone::DeviceVector<unsigned int>,
               cstone::DeviceVector<unsigned int>, cstone::DeviceVector<unsigned int>,
               cstone::DeviceVector<unsigned int>, cstone::DeviceVector<unsigned int>>,
    OrigCoords24<double>>(
    cstone::Domain<unsigned int, double, cstone::GpuTag>* domain,
    cstone::DeviceVector<unsigned int>& elemSfcCodes,
    cstone::DeviceVector<double>& elemX,
    cstone::DeviceVector<double>& elemY,
    cstone::DeviceVector<double>& elemZ,
    cstone::DeviceVector<double>& elemH,
    size_t& elementCount,
    std::tuple<cstone::DeviceVector<unsigned int>, cstone::DeviceVector<unsigned int>,
               cstone::DeviceVector<unsigned int>, cstone::DeviceVector<unsigned int>,
               cstone::DeviceVector<unsigned int>, cstone::DeviceVector<unsigned int>,
               cstone::DeviceVector<unsigned int>, cstone::DeviceVector<unsigned int>>& d_conn_keys_,
    OrigCoords24<double>& d_orig_coords_);

// uint64_t keys, float coords
template void syncDomainImplWithOrigCoords<
    uint64_t, float,
    std::tuple<cstone::DeviceVector<uint64_t>, cstone::DeviceVector<uint64_t>,
               cstone::DeviceVector<uint64_t>, cstone::DeviceVector<uint64_t>,
               cstone::DeviceVector<uint64_t>, cstone::DeviceVector<uint64_t>,
               cstone::DeviceVector<uint64_t>, cstone::DeviceVector<uint64_t>>,
    OrigCoords24<float>>(
    cstone::Domain<uint64_t, float, cstone::GpuTag>* domain,
    cstone::DeviceVector<uint64_t>& elemSfcCodes,
    cstone::DeviceVector<float>& elemX,
    cstone::DeviceVector<float>& elemY,
    cstone::DeviceVector<float>& elemZ,
    cstone::DeviceVector<float>& elemH,
    size_t& elementCount,
    std::tuple<cstone::DeviceVector<uint64_t>, cstone::DeviceVector<uint64_t>,
               cstone::DeviceVector<uint64_t>, cstone::DeviceVector<uint64_t>,
               cstone::DeviceVector<uint64_t>, cstone::DeviceVector<uint64_t>,
               cstone::DeviceVector<uint64_t>, cstone::DeviceVector<uint64_t>>& d_conn_keys_,
    OrigCoords24<float>& d_orig_coords_);

// uint64_t keys, double coords
template void syncDomainImplWithOrigCoords<
    uint64_t, double,
    std::tuple<cstone::DeviceVector<uint64_t>, cstone::DeviceVector<uint64_t>,
               cstone::DeviceVector<uint64_t>, cstone::DeviceVector<uint64_t>,
               cstone::DeviceVector<uint64_t>, cstone::DeviceVector<uint64_t>,
               cstone::DeviceVector<uint64_t>, cstone::DeviceVector<uint64_t>>,
    OrigCoords24<double>>(
    cstone::Domain<uint64_t, double, cstone::GpuTag>* domain,
    cstone::DeviceVector<uint64_t>& elemSfcCodes,
    cstone::DeviceVector<double>& elemX,
    cstone::DeviceVector<double>& elemY,
    cstone::DeviceVector<double>& elemZ,
    cstone::DeviceVector<double>& elemH,
    size_t& elementCount,
    std::tuple<cstone::DeviceVector<uint64_t>, cstone::DeviceVector<uint64_t>,
               cstone::DeviceVector<uint64_t>, cstone::DeviceVector<uint64_t>,
               cstone::DeviceVector<uint64_t>, cstone::DeviceVector<uint64_t>,
               cstone::DeviceVector<uint64_t>, cstone::DeviceVector<uint64_t>>& d_conn_keys_,
    OrigCoords24<double>& d_orig_coords_);

} // namespace mars
