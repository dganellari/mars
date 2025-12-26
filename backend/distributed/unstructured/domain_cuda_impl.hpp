#pragma once

#include <cstddef> // For size_t definition

// Handle both CUDA and HIP includes
#ifdef MARS_ENABLE_CUDA
#include <cuda_runtime.h>
#endif

#ifdef MARS_ENABLE_HIP
#include <hip/hip_runtime.h>
#endif

namespace cstone
{
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
} // namespace cstone

namespace mars
{
// Forward declaration of sync implementation function
template<typename KeyType, typename RealType, typename SfcConnTuple>
void syncDomainImpl(cstone::Domain<KeyType, RealType, cstone::GpuTag>* domain,
                    cstone::DeviceVector<KeyType>& elemSfcCodes,
                    cstone::DeviceVector<RealType>& elemX,
                    cstone::DeviceVector<RealType>& elemY,
                    cstone::DeviceVector<RealType>& elemZ,
                    cstone::DeviceVector<RealType>& elemH,
                    size_t& elementCount,
                    SfcConnTuple& d_conn_keys_);

// Overload for syncing with original coordinates (24 additional properties for hex8)
template<typename KeyType, typename RealType, typename SfcConnTuple, typename OrigCoordsTuple>
void syncDomainImplWithOrigCoords(cstone::Domain<KeyType, RealType, cstone::GpuTag>* domain,
                                   cstone::DeviceVector<KeyType>& elemSfcCodes,
                                   cstone::DeviceVector<RealType>& elemX,
                                   cstone::DeviceVector<RealType>& elemY,
                                   cstone::DeviceVector<RealType>& elemZ,
                                   cstone::DeviceVector<RealType>& elemH,
                                   size_t& elementCount,
                                   SfcConnTuple& d_conn_keys_,
                                   OrigCoordsTuple& d_orig_coords_);

// Forward declarations of CUDA kernels
template<typename KeyType, typename SfcConnTuple>
__global__ void extractAllTupleComponentsKernel(const SfcConnTuple conn_keys,
                                                KeyType* flattenedKeys,
                                                int numElements);

template<typename KeyType>
__global__ void connectivityRebuildKernel(
    const KeyType* sfcKeys, int* connectivity, const KeyType* uniqueSfcKeys, int numElements, int numUniqueKeys);

} // namespace mars
