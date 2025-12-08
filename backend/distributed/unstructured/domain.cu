#include "domain.hpp"
#include "thrust/sort.h"
#include "thrust/unique.h"
#include "thrust/device_vector.h"
#include "cub/cub.cuh"

namespace mars
{

// Single source of truth for SFC-to-physical conversion
template<typename KeyType, typename RealType>
__device__ __host__ std::tuple<RealType, RealType, RealType> decodeSfcToPhysical(KeyType sfcKey, const cstone::Box<RealType>& box) {
    // Convert raw key to SfcKind strong type
    auto sfcKindKey = cstone::SfcKind<KeyType>(sfcKey);
    auto [ix, iy, iz] = cstone::decodeSfc(sfcKindKey);
    
    // Use SfcKind for maxTreeLevel
    constexpr unsigned maxCoord = (1u << cstone::maxTreeLevel<cstone::SfcKind<KeyType>>{}) - 1;
    RealType invMaxCoord = RealType(1.0) / maxCoord;
    
    RealType x = box.xmin() + ix * invMaxCoord * (box.xmax() - box.xmin());
    RealType y = box.ymin() + iy * invMaxCoord * (box.ymax() - box.ymin());
    RealType z = box.zmin() + iz * invMaxCoord * (box.zmax() - box.zmin());
    
    return std::make_tuple(x, y, z);
}

// CUDA kernels with RealType template parameter instead of Real
template<typename RealType>
__global__ void
transformCharacteristicSizesKernel(RealType* d_h, size_t size, RealType meshFactor, RealType minH, RealType maxH)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < size)
    {
        RealType val    = d_h[idx];
        RealType result = val * meshFactor;

        if constexpr (std::is_same_v<RealType, float>) {
            d_h[idx] = fmaxf(minH, fminf(maxH, result));
        } else {
            d_h[idx] = fmax(minH, fmin(maxH, result));  // For double
        }
    }
}

template<typename RealType>
__global__ void fillCharacteristicSizesKernel(RealType* d_h, size_t size, RealType value)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < size) { d_h[idx] = value; }
}

// CUDA kernel to calculate element characteristic sizes
template<typename ElementTag, typename KeyType, typename RealType>
__global__ void computeCharacteristicSizesKernel(const RealType* x,
                                                 const RealType* y,
                                                 const RealType* z,
                                                 const KeyType* indices0,
                                                 const KeyType* indices1,
                                                 const KeyType* indices2,
                                                 const KeyType* indices3,
                                                 const KeyType* indices4,
                                                 const KeyType* indices5,
                                                 const KeyType* indices6,
                                                 const KeyType* indices7,
                                                 int* nodeTetCount,
                                                 RealType* h,
                                                 int numElements)
{
    int elemIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (elemIdx < numElements)
    {
        if constexpr (std::is_same_v<ElementTag, TetTag>)
        {
            // Get the four nodes of this tetrahedron
            auto n0 = indices0[elemIdx];
            auto n1 = indices1[elemIdx];
            auto n2 = indices2[elemIdx];
            auto n3 = indices3[elemIdx];

            // Calculate edge lengths and contribute to characteristic size
            // Edge n0-n1
            RealType dx         = x[n0] - x[n1];
            RealType dy         = y[n0] - y[n1];
            RealType dz         = z[n0] - z[n1];
            RealType edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n0], edgeLength);
            atomicAdd(&h[n1], edgeLength);
            atomicAdd(&nodeTetCount[n0], 1);
            atomicAdd(&nodeTetCount[n1], 1);

            // Edge n0-n2
            dx         = x[n0] - x[n2];
            dy         = y[n0] - y[n2];
            dz         = z[n0] - z[n2];
            edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n0], edgeLength);
            atomicAdd(&h[n2], edgeLength);
            atomicAdd(&nodeTetCount[n0], 1);
            atomicAdd(&nodeTetCount[n2], 1);

            // Edge n0-n3
            dx         = x[n0] - x[n3];
            dy         = y[n0] - y[n3];
            dz         = z[n0] - z[n3];
            edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n0], edgeLength);
            atomicAdd(&h[n3], edgeLength);
            atomicAdd(&nodeTetCount[n0], 1);
            atomicAdd(&nodeTetCount[n3], 1);

            // Edge n1-n2
            dx         = x[n1] - x[n2];
            dy         = y[n1] - y[n2];
            dz         = z[n1] - z[n2];
            edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n1], edgeLength);
            atomicAdd(&h[n2], edgeLength);
            atomicAdd(&nodeTetCount[n1], 1);
            atomicAdd(&nodeTetCount[n2], 1);

            // Edge n1-n3
            dx         = x[n1] - x[n3];
            dy         = y[n1] - y[n3];
            dz         = z[n1] - z[n3];
            edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n1], edgeLength);
            atomicAdd(&h[n3], edgeLength);
            atomicAdd(&nodeTetCount[n1], 1);
            atomicAdd(&nodeTetCount[n3], 1);

            // Edge n2-n3
            dx         = x[n2] - x[n3];
            dy         = y[n2] - y[n3];
            dz         = z[n2] - z[n3];
            edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n2], edgeLength);
            atomicAdd(&h[n3], edgeLength);
            atomicAdd(&nodeTetCount[n2], 1);
            atomicAdd(&nodeTetCount[n3], 1);
        }
        else if constexpr (std::is_same_v<ElementTag, HexTag>)
        {
            // Get all eight nodes of this hexahedron
            auto n0 = indices0[elemIdx];
            auto n1 = indices1[elemIdx];
            auto n2 = indices2[elemIdx];
            auto n3 = indices3[elemIdx];
            auto n4 = indices4[elemIdx];
            auto n5 = indices5[elemIdx];
            auto n6 = indices6[elemIdx];
            auto n7 = indices7[elemIdx];
            
            // Calculate edge lengths for all 12 edges of the hexahedron
            RealType dx, dy, dz, edgeLength;
            
            // Bottom face edges (0-1, 1-2, 2-3, 3-0)
            dx = x[n0] - x[n1]; dy = y[n0] - y[n1]; dz = z[n0] - z[n1];
            edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n0], edgeLength); atomicAdd(&h[n1], edgeLength);
            atomicAdd(&nodeTetCount[n0], 1); atomicAdd(&nodeTetCount[n1], 1);

            dx = x[n1] - x[n2]; dy = y[n1] - y[n2]; dz = z[n1] - z[n2];
            edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n1], edgeLength); atomicAdd(&h[n2], edgeLength);
            atomicAdd(&nodeTetCount[n1], 1); atomicAdd(&nodeTetCount[n2], 1);

            dx = x[n2] - x[n3]; dy = y[n2] - y[n3]; dz = z[n2] - z[n3];
            edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n2], edgeLength); atomicAdd(&h[n3], edgeLength);
            atomicAdd(&nodeTetCount[n2], 1); atomicAdd(&nodeTetCount[n3], 1);

            dx = x[n3] - x[n0]; dy = y[n3] - y[n0]; dz = z[n3] - z[n0];
            edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n3], edgeLength); atomicAdd(&h[n0], edgeLength);
            atomicAdd(&nodeTetCount[n3], 1); atomicAdd(&nodeTetCount[n0], 1);

            // Top face edges (4-5, 5-6, 6-7, 7-4)
            dx = x[n4] - x[n5]; dy = y[n4] - y[n5]; dz = z[n4] - z[n5];
            edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n4], edgeLength); atomicAdd(&h[n5], edgeLength);
            atomicAdd(&nodeTetCount[n4], 1); atomicAdd(&nodeTetCount[n5], 1);

            dx = x[n5] - x[n6]; dy = y[n5] - y[n6]; dz = z[n5] - z[n6];
            edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n5], edgeLength); atomicAdd(&h[n6], edgeLength);
            atomicAdd(&nodeTetCount[n5], 1); atomicAdd(&nodeTetCount[n6], 1);

            dx = x[n6] - x[n7]; dy = y[n6] - y[n7]; dz = z[n6] - z[n7];
            edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n6], edgeLength); atomicAdd(&h[n7], edgeLength);
            atomicAdd(&nodeTetCount[n6], 1); atomicAdd(&nodeTetCount[n7], 1);

            dx = x[n7] - x[n4]; dy = y[n7] - y[n4]; dz = z[n7] - z[n4];
            edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n7], edgeLength); atomicAdd(&h[n4], edgeLength);
            atomicAdd(&nodeTetCount[n7], 1); atomicAdd(&nodeTetCount[n4], 1);

            // Vertical edges (0-4, 1-5, 2-6, 3-7)
            dx = x[n0] - x[n4]; dy = y[n0] - y[n4]; dz = z[n0] - z[n4];
            edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n0], edgeLength); atomicAdd(&h[n4], edgeLength);
            atomicAdd(&nodeTetCount[n0], 1); atomicAdd(&nodeTetCount[n4], 1);

            dx = x[n1] - x[n5]; dy = y[n1] - y[n5]; dz = z[n1] - z[n5];
            edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n1], edgeLength); atomicAdd(&h[n5], edgeLength);
            atomicAdd(&nodeTetCount[n1], 1); atomicAdd(&nodeTetCount[n5], 1);

            dx = x[n2] - x[n6]; dy = y[n2] - y[n6]; dz = z[n2] - z[n6];
            edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n2], edgeLength); atomicAdd(&h[n6], edgeLength);
            atomicAdd(&nodeTetCount[n2], 1); atomicAdd(&nodeTetCount[n6], 1);

            dx = x[n3] - x[n7]; dy = y[n3] - y[n7]; dz = z[n3] - z[n7];
            edgeLength = sqrt(dx * dx + dy * dy + dz * dz);
            atomicAdd(&h[n3], edgeLength); atomicAdd(&h[n7], edgeLength);
            atomicAdd(&nodeTetCount[n3], 1); atomicAdd(&nodeTetCount[n7], 1);
        }
    }
}

template<typename KeyType, typename RealType>
__global__ void finalizeCharacteristicSizesKernel(RealType* h, int* nodeTetCount, int numNodes)
{
    int nodeIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (nodeIdx < numNodes)
    {
        if (nodeTetCount[nodeIdx] > 0) { h[nodeIdx] /= nodeTetCount[nodeIdx]; }
        else
        {
            h[nodeIdx] = 0.01; // Default for isolated nodes
        }
    }
}

// Generic kernel for finding representative nodes
template<typename ElementTag, typename KeyType, typename RealType>
__global__ void findRepresentativeNodesKernel(const KeyType* indices0,
                                              const KeyType* indices1,
                                              const KeyType* indices2,
                                              const KeyType* indices3,
                                              const KeyType* indices4,
                                              const KeyType* indices5,
                                              const KeyType* indices6,
                                              const KeyType* indices7,
                                              const KeyType* sfcCodes,
                                              KeyType* elemToNodeMap,
                                              int numElements)
{
    int elemIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (elemIdx < numElements)
    {
        if constexpr (std::is_same_v<ElementTag, TetTag>)
        {
            // Get the four nodes of this tetrahedron
            auto node0 = indices0[elemIdx];
            auto node1 = indices1[elemIdx];
            auto node2 = indices2[elemIdx];
            auto node3 = indices3[elemIdx];

            // Get SFC codes
            auto sfc0 = sfcCodes[node0];
            auto sfc1 = sfcCodes[node1];
            auto sfc2 = sfcCodes[node2];
            auto sfc3 = sfcCodes[node3];

            // Find minimum SFC code
            auto repNode     = node0;
            auto minSfc = sfc0;

            if (sfc1 < minSfc)
            {
                minSfc  = sfc1;
                repNode = node1;
            }

            if (sfc2 < minSfc)
            {
                minSfc  = sfc2;
                repNode = node2;
            }

            if (sfc3 < minSfc)
            {
                minSfc  = sfc3;
                repNode = node3;
            }

            // Store the representative node
            elemToNodeMap[elemIdx] = repNode;
        }
        else if constexpr (std::is_same_v<ElementTag, HexTag>)
        {
            // For hexahedra - find node with minimum SFC code among all 8 nodes
            auto node0 = indices0[elemIdx];
            auto node1 = indices1[elemIdx];
            auto node2 = indices2[elemIdx];
            auto node3 = indices3[elemIdx];
            auto node4 = indices4[elemIdx];
            auto node5 = indices5[elemIdx];
            auto node6 = indices6[elemIdx];
            auto node7 = indices7[elemIdx];

            auto sfc0 = sfcCodes[node0];
            auto sfc1 = sfcCodes[node1];
            auto sfc2 = sfcCodes[node2];
            auto sfc3 = sfcCodes[node3];
            auto sfc4 = sfcCodes[node4];
            auto sfc5 = sfcCodes[node5];
            auto sfc6 = sfcCodes[node6];
            auto sfc7 = sfcCodes[node7];

            auto repNode = node0;
            auto minSfc = sfc0;

            if (sfc1 < minSfc) { minSfc = sfc1; repNode = node1; }
            if (sfc2 < minSfc) { minSfc = sfc2; repNode = node2; }
            if (sfc3 < minSfc) { minSfc = sfc3; repNode = node3; }
            if (sfc4 < minSfc) { minSfc = sfc4; repNode = node4; }
            if (sfc5 < minSfc) { minSfc = sfc5; repNode = node5; }
            if (sfc6 < minSfc) { minSfc = sfc6; repNode = node6; }
            if (sfc7 < minSfc) { minSfc = sfc7; repNode = node7; }

            elemToNodeMap[elemIdx] = repNode;
        }
        else if constexpr (std::is_same_v<ElementTag, TriTag>)
        {
            // For triangles - would use only indices0-2
            auto node0 = indices0[elemIdx];
            auto node1 = indices1[elemIdx];
            auto node2 = indices2[elemIdx];

            auto sfc0 = sfcCodes[node0];
            auto sfc1 = sfcCodes[node1];
            auto sfc2 = sfcCodes[node2];

            auto repNode     = node0;
            auto minSfc = sfc0;

            if (sfc1 < minSfc)
            {
                minSfc  = sfc1;
                repNode = node1;
            }

            if (sfc2 < minSfc)
            {
                minSfc  = sfc2;
                repNode = node2;
            }

            elemToNodeMap[elemIdx] = repNode;
        }
    }
}

template<typename KeyType, typename RealType>
void generateSfcKeys(const RealType* x,
                     const RealType* y,
                     const RealType* z,
                     KeyType* keys,
                     size_t numKeys,
                     const cstone::Box<RealType>& box)
{
    // Use sfcKindPointer to match cornerstone's template instantiation
    cstone::computeSfcKeysGpu(x, y, z, cstone::sfcKindPointer(keys), numKeys, box);
    cudaCheckError();
}

// Kernel to extract representative node coordinates and compute SFC keys
template<typename ElementTag, typename KeyType, typename RealType>
__global__ void extractRepCoordinatesKernel(const RealType* x,
                                            const RealType* y,
                                            const RealType* z,
                                            const RealType* h,
                                            const KeyType* elemToNodeMap,
                                            RealType* elemX,
                                            RealType* elemY,
                                            RealType* elemZ,
                                            RealType* elemH,
                                            int numElements)
{
    int elemIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (elemIdx < numElements)
    {
        auto repNodeIdx = elemToNodeMap[elemIdx];

        // Extract coordinates and properties
        elemX[elemIdx] = x[repNodeIdx];
        elemY[elemIdx] = y[repNodeIdx];
        elemZ[elemIdx] = z[repNodeIdx];
        elemH[elemIdx] = h[repNodeIdx];
    }
}

template<typename KeyType>
__global__ void extractAllTupleComponentsKernel(const KeyType* conn_key0,
                                               const KeyType* conn_key1,
                                               const KeyType* conn_key2,
                                               const KeyType* conn_key3,
                                               KeyType* flattenedKeys,
                                               int numElements,
                                               int nodeCount)
{
    int elemIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (elemIdx < numElements)
    {
        flattenedKeys[elemIdx * nodeCount + 0] = conn_key0[elemIdx];
        if (nodeCount >= 2) flattenedKeys[elemIdx * nodeCount + 1] = conn_key1[elemIdx];
        if (nodeCount >= 3) flattenedKeys[elemIdx * nodeCount + 2] = conn_key2[elemIdx];
        if (nodeCount >= 4) flattenedKeys[elemIdx * nodeCount + 3] = conn_key3[elemIdx];
    }
}

template<typename KeyType>
__global__ void connectivityRebuildKernel(
    const KeyType* sfcKeys, KeyType* connectivity, const KeyType* uniqueSfcKeys, int numElements, int numUniqueKeys)
{
    int elemIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (elemIdx < numElements)
    {
        KeyType key = sfcKeys[elemIdx];
        // Use CUB's device-side binary search
        auto localId = cub::LowerBound(uniqueSfcKeys, numUniqueKeys, key);
        // If the key is not found, localId will be numUniqueKeys
        connectivity[elemIdx] = localId;
        // Ensure the key is valid
        // This assertion is for debugging purposes, it will fail if the key is not found
        assert(localId < numUniqueKeys && uniqueSfcKeys[localId] == key);
    }
}

// Generic kernel for finding representative nodes
template<typename ElementTag, typename KeyType, typename RealType>
__global__ void buildSfcConnectivity(const KeyType* indices0,
                                     const KeyType* indices1,
                                     const KeyType* indices2,
                                     const KeyType* indices3,
                                     const KeyType* indices4,
                                     const KeyType* indices5,
                                     const KeyType* indices6,
                                     const KeyType* indices7,
                                     const KeyType* sfcCodes,
                                     KeyType* keys0,
                                     KeyType* keys1,
                                     KeyType* keys2,
                                     KeyType* keys3,
                                     KeyType* keys4,
                                     KeyType* keys5,
                                     KeyType* keys6,
                                     KeyType* keys7,
                                     int numElements)
{
    int elemIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (elemIdx < numElements)
    {
        if constexpr (std::is_same_v<ElementTag, TetTag>)
        {
            // Get the four nodes of this tetrahedron
            auto node0 = indices0[elemIdx];
            auto node1 = indices1[elemIdx];
            auto node2 = indices2[elemIdx];
            auto node3 = indices3[elemIdx];

            keys0[elemIdx] = sfcCodes[node0];
            keys1[elemIdx] = sfcCodes[node1];
            keys2[elemIdx] = sfcCodes[node2];
            keys3[elemIdx] = sfcCodes[node3];
        }
        else if constexpr (std::is_same_v<ElementTag, HexTag>)
        {
            // For hexahedra - store SFC keys for all 8 nodes
            auto node0 = indices0[elemIdx];
            auto node1 = indices1[elemIdx];
            auto node2 = indices2[elemIdx];
            auto node3 = indices3[elemIdx];
            auto node4 = indices4[elemIdx];
            auto node5 = indices5[elemIdx];
            auto node6 = indices6[elemIdx];
            auto node7 = indices7[elemIdx];

            keys0[elemIdx] = sfcCodes[node0];
            keys1[elemIdx] = sfcCodes[node1];
            keys2[elemIdx] = sfcCodes[node2];
            keys3[elemIdx] = sfcCodes[node3];
            keys4[elemIdx] = sfcCodes[node4];
            keys5[elemIdx] = sfcCodes[node5];
            keys6[elemIdx] = sfcCodes[node6];
            keys7[elemIdx] = sfcCodes[node7];
        }
        else if constexpr (std::is_same_v<ElementTag, TriTag>)
        {
            // For triangles - would use only indices0-2
            auto node0 = indices0[elemIdx];
            auto node1 = indices1[elemIdx];
            auto node2 = indices2[elemIdx];

            keys0[elemIdx] = sfcCodes[node0];
            keys1[elemIdx] = sfcCodes[node1];
            keys2[elemIdx] = sfcCodes[node2];
        }
    }
}

template<typename KeyType, typename SfcConnTuple, typename ConnTuple>
void rebuildElementConnectivity(
    SfcConnTuple& d_conn_keys_,
    ConnTuple& d_conn_,
    size_t newElementCount)
{
    constexpr size_t NodesPerElement = std::tuple_size_v<SfcConnTuple>;
    size_t totalNodeEntries = newElementCount * NodesPerElement;

    thrust::device_vector<KeyType> d_flattenedKeys(totalNodeEntries);

    int blockSize = 256;
    int numBlocks = (newElementCount + blockSize - 1) / blockSize;

    extractAllTupleComponentsKernel<KeyType><<<numBlocks, blockSize>>>(
        thrust::raw_pointer_cast(std::get<0>(d_conn_keys_).data()),
        NodesPerElement >= 2 ? thrust::raw_pointer_cast(std::get<1>(d_conn_keys_).data()) : nullptr,
        NodesPerElement >= 3 ? thrust::raw_pointer_cast(std::get<2>(d_conn_keys_).data()) : nullptr,
        NodesPerElement >= 4 ? thrust::raw_pointer_cast(std::get<3>(d_conn_keys_).data()) : nullptr,
        thrust::raw_pointer_cast(d_flattenedKeys.data()),
        newElementCount, NodesPerElement);
    cudaCheckError();

    thrust::sort(d_flattenedKeys.begin(), d_flattenedKeys.end());

    thrust::device_vector<KeyType> d_uniqueSfcKeys(totalNodeEntries);
    auto endIter = thrust::unique_copy(d_flattenedKeys.begin(), d_flattenedKeys.end(), d_uniqueSfcKeys.begin());

    size_t uniqueKeyCount = endIter - d_uniqueSfcKeys.begin();
    d_uniqueSfcKeys.resize(uniqueKeyCount);

    ConnTuple newConn;
    if constexpr (NodesPerElement >= 1) { std::get<0>(newConn).resize(newElementCount); }
    if constexpr (NodesPerElement >= 2) { std::get<1>(newConn).resize(newElementCount); }
    if constexpr (NodesPerElement >= 3) { std::get<2>(newConn).resize(newElementCount); }
    if constexpr (NodesPerElement >= 4) { std::get<3>(newConn).resize(newElementCount); }

    // Use thrust's optimized binary search
    if constexpr (NodesPerElement >= 1)
    {
        connectivityRebuildKernel<<<numBlocks, blockSize>>>(
            thrust::raw_pointer_cast(std::get<0>(d_conn_keys_).data()),
            thrust::raw_pointer_cast(std::get<0>(newConn).data()), 
            thrust::raw_pointer_cast(d_uniqueSfcKeys.data()),
            newElementCount, uniqueKeyCount);
        cudaCheckError();
    }

    if constexpr (NodesPerElement >= 2)
    {
        connectivityRebuildKernel<<<numBlocks, blockSize>>>(
            thrust::raw_pointer_cast(std::get<1>(d_conn_keys_).data()),
            thrust::raw_pointer_cast(std::get<1>(newConn).data()), 
            thrust::raw_pointer_cast(d_uniqueSfcKeys.data()),
            newElementCount, uniqueKeyCount);
        cudaCheckError();
    }

    if constexpr (NodesPerElement >= 3)
    {
        connectivityRebuildKernel<<<numBlocks, blockSize>>>(
            thrust::raw_pointer_cast(std::get<2>(d_conn_keys_).data()),
            thrust::raw_pointer_cast(std::get<2>(newConn).data()), 
            thrust::raw_pointer_cast(d_uniqueSfcKeys.data()),
            newElementCount, uniqueKeyCount);
        cudaCheckError();
    }

    if constexpr (NodesPerElement >= 4)
    {
        connectivityRebuildKernel<<<numBlocks, blockSize>>>(
            thrust::raw_pointer_cast(std::get<3>(d_conn_keys_).data()),
            thrust::raw_pointer_cast(std::get<3>(newConn).data()), 
            thrust::raw_pointer_cast(d_uniqueSfcKeys.data()),
            newElementCount, uniqueKeyCount);
        cudaCheckError();
    }

    std::swap(d_conn_, newConn);
}

template __global__ void connectivityRebuildKernel<unsigned int>
    (const unsigned int*, unsigned int*, const unsigned int*, int, int);

template __global__ void connectivityRebuildKernel<uint64_t>
    (const uint64_t*, uint64_t*, const uint64_t*, int, int);

template void rebuildElementConnectivity<unsigned int, 
    std::tuple<cstone::DeviceVector<unsigned int>, cstone::DeviceVector<unsigned int>, 
               cstone::DeviceVector<unsigned int>, cstone::DeviceVector<unsigned int>>,
    std::tuple<cstone::DeviceVector<unsigned int>, cstone::DeviceVector<unsigned int>, 
               cstone::DeviceVector<unsigned int>, cstone::DeviceVector<unsigned int>>>(
    std::tuple<cstone::DeviceVector<unsigned int>, cstone::DeviceVector<unsigned int>, 
               cstone::DeviceVector<unsigned int>, cstone::DeviceVector<unsigned int>>& d_conn_keys_,
    std::tuple<cstone::DeviceVector<unsigned int>, cstone::DeviceVector<unsigned int>, 
               cstone::DeviceVector<unsigned int>, cstone::DeviceVector<unsigned int>>& d_conn_,
    size_t newElementCount);

template void rebuildElementConnectivity<uint64_t, 
    std::tuple<cstone::DeviceVector<uint64_t>, cstone::DeviceVector<uint64_t>, 
               cstone::DeviceVector<uint64_t>, cstone::DeviceVector<uint64_t>>,
    std::tuple<cstone::DeviceVector<uint64_t>, cstone::DeviceVector<uint64_t>,    
               cstone::DeviceVector<uint64_t>, cstone::DeviceVector<uint64_t>>>(  
    std::tuple<cstone::DeviceVector<uint64_t>, cstone::DeviceVector<uint64_t>, 
               cstone::DeviceVector<uint64_t>, cstone::DeviceVector<uint64_t>>& d_conn_keys_,
    std::tuple<cstone::DeviceVector<uint64_t>, cstone::DeviceVector<uint64_t>,    
               cstone::DeviceVector<uint64_t>, cstone::DeviceVector<uint64_t>>& d_conn_,  
    size_t newElementCount);

// Compute tetrahedron volumes
template<typename ElementTag, typename RealType>
__global__ void computeElementVolumesKernel(const RealType* x,
                                            const RealType* y,
                                            const RealType* z,
                                            const int* indices0,
                                            const int* indices1,
                                            const int* indices2,
                                            const int* indices3,
                                            RealType* volumes,
                                            int numElements)
{
    int elemIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (elemIdx < numElements)
    {
        if constexpr (std::is_same_v<ElementTag, TetTag>)
        {
            // Get node indices
            int n0 = indices0[elemIdx];
            int n1 = indices1[elemIdx];
            int n2 = indices2[elemIdx];
            int n3 = indices3[elemIdx];

            // Get node coordinates
            RealType x0 = x[n0], y0 = y[n0], z0 = z[n0];
            RealType x1 = x[n1], y1 = y[n1], z1 = z[n1];
            RealType x2 = x[n2], y2 = y[n2], z2 = z[n2];
            RealType x3 = x[n3], y3 = y[n3], z3 = z[n3];

            // Compute vectors for volume calculation
            RealType v1x = x1 - x0, v1y = y1 - y0, v1z = z1 - z0;
            RealType v2x = x2 - x0, v2y = y2 - y0, v2z = z2 - z0;
            RealType v3x = x3 - x0, v3y = y3 - y0, v3z = z3 - z0;

            // Compute volume using the scalar triple product
            volumes[elemIdx] =
                fabs(v1x * (v2y * v3z - v2z * v3y) + v1y * (v2z * v3x - v2x * v3z) + v1z * (v2x * v3y - v2y * v3x)) /
                6.0;
        }
        else if constexpr (std::is_same_v<ElementTag, HexTag>)
        {
            // Hexahedron volume calculation would go here
            // Would need additional indices parameters
        }
    }
}

template<typename KeyType, size_t NodesPerElement>
__global__ void flattenConnectivityKernel(
    ConnPtrs<KeyType, NodesPerElement> conn,
    KeyType* flat_keys, 
    size_t numElements)
{
    size_t elementIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (elementIdx >= numElements) return;

    size_t baseIdx = elementIdx * NodesPerElement;
    
    // Clean loop - compiler will unroll automatically
    #pragma unroll
    for (int i = 0; i < NodesPerElement; ++i) {
        flat_keys[baseIdx + i] = conn.ptrs[i][elementIdx];
    }
}

template<typename KeyType>
__global__ void mapSfcToLocalIdKernel(const KeyType* sfc_conn, 
                                      KeyType* local_conn, 
                                      const KeyType* sorted_sfc, 
                                      size_t num_elements, 
                                      size_t num_nodes)
{
    size_t elemIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (elemIdx >= num_elements) return;

    KeyType sfc_key = sfc_conn[elemIdx];
    
    // Binary search using CUB for optimal performance
    KeyType local_id = cub::LowerBound(sorted_sfc, static_cast<int>(num_nodes), sfc_key);
    
    // Verify the key was found (debug check)
    assert(local_id < num_nodes && sorted_sfc[local_id] == sfc_key);
    
    local_conn[elemIdx] = local_id;
}

template<typename KeyType, typename RealType>
__global__ void decodeAllNodesKernel(const KeyType* sfcKeys, RealType* x, RealType* y, RealType* z,
                                     size_t numNodes, cstone::Box<RealType> box) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numNodes) return;
    
    auto [xi, yi, zi] = decodeSfcToPhysical(sfcKeys[idx], box);
    x[idx] = xi; y[idx] = yi; z[idx] = zi;
}

template <typename KeyType, int NodesPerElement>
__global__ void markNodesInElementRangeKernel(
    unsigned int* nodeFlags,
    const KeyType* conn0,
    const KeyType* conn1,
    const KeyType* conn2,
    const KeyType* conn3,
    const KeyType* conn4,
    const KeyType* conn5,
    const KeyType* conn6,
    const KeyType* conn7,
    const KeyType* localToGlobalSfcMap,
    size_t elementStart,
    size_t elementEnd,
    size_t nodeCount,
    unsigned int flagValue)
{
    size_t elemIdx = blockIdx.x * blockDim.x + threadIdx.x + elementStart;
    if (elemIdx < elementEnd)
    {
        // For each node in this element, find its local node ID and set flag
        KeyType sfcKeys[8];
        sfcKeys[0] = conn0[elemIdx];
        sfcKeys[1] = conn1[elemIdx];
        sfcKeys[2] = conn2[elemIdx];
        sfcKeys[3] = conn3[elemIdx];
        if constexpr (NodesPerElement > 4) {
            sfcKeys[4] = conn4[elemIdx];
            sfcKeys[5] = conn5[elemIdx];
            sfcKeys[6] = conn6[elemIdx];
            sfcKeys[7] = conn7[elemIdx];
        }
        
        for (int i = 0; i < NodesPerElement; ++i) {
            KeyType targetSfc = sfcKeys[i];
            if (targetSfc == 0 && NodesPerElement > 4 && i >= 4) continue; // Skip padded zeros for tet in hex domain
            // Binary search to find this SFC key in localToGlobalSfcMap
            int left = 0, right = nodeCount - 1;
            while (left <= right) {
                int mid = left + (right - left) / 2;
                KeyType midSfc = localToGlobalSfcMap[mid];
                if (midSfc == targetSfc) {
                    atomicOr(&nodeFlags[mid], flagValue);
                    break;
                }
                if (midSfc < targetSfc) left = mid + 1;
                else right = mid - 1;
            }
        }
    }
}

template <typename KeyType>
__global__ void determineOwnershipDirectKernel(
    uint8_t* nodeOwnership, 
    const KeyType* nodeSfcKeys,
    const KeyType* assignment, 
    int numRanks, 
    int myRank, 
    size_t numNodes)
{
    size_t nodeIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (nodeIdx < numNodes)
    {
        KeyType sfc = nodeSfcKeys[nodeIdx];
        
        // Find owner rank using binary search on assignment array
        int owner = -1;
        int left = 0;
        int right = numRanks - 1;
        
        while (left <= right) {
            int mid = left + (right - left) / 2;
            if (sfc >= assignment[mid] && sfc < assignment[mid+1]) {
                owner = mid;
                break;
            }
            if (sfc < assignment[mid]) {
                right = mid - 1;
            } else {
                left = mid + 1;
            }
        }
        
        // Mark ownership: 1 if owned by this rank, 0 if ghost
        nodeOwnership[nodeIdx] = (owner == myRank) ? 1 : 0;
    }
}

template <typename KeyType>
__global__ void detectSharedNodesKernel(
    uint8_t* nodeOwnership,
    const unsigned int* nodeFlags,
    size_t numNodes)
{
    size_t nodeIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (nodeIdx < numNodes)
    {
        unsigned int flags = nodeFlags[nodeIdx];
        uint8_t ownership = nodeOwnership[nodeIdx];
        
        // Upgrade owned nodes to shared if they appear in both local and halo elements
        // Note: Some owned nodes may not appear in local elements (only in halo)
        // These will have zero diagonals but we keep them as owned to avoid orphaning them
        
        if (ownership == 1 && (flags & 0x03) == 0x03) {
            // In both local and halo → shared boundary
            nodeOwnership[nodeIdx] = 2;
        }
        // Keep SFC-based ownership otherwise (don't downgrade based on element participation)
    }
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void HaloData<ElementTag, RealType, KeyType, AcceleratorTag>::buildNodeOwnership(
    const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain)
{
    std::cout << "Rank " << domain.rank() << ": Building node ownership..." << std::endl;
    std::cout.flush();
    
    size_t nodeCount = domain.getNodeCount();
    d_nodeOwnership_.resize(nodeCount);
    
    if (domain.numRanks() == 1)
    {
        thrust::fill(thrust::device, 
                    thrust::device_pointer_cast(d_nodeOwnership_.data()),
                    thrust::device_pointer_cast(d_nodeOwnership_.data() + nodeCount), 
                    1); // All owned
        return;
    }

    
    // NEW APPROACH: SFC-key-based ownership (matches cornerstone's assignment)
    // This ensures neighbor detection works correctly with findRank()
    // - Nodes with SFC keys in [assignmentStart, assignmentEnd) → owned (1)
    // - Shared boundaries will be detected by checking if nodes appear in both local and halo elements
    
    const auto& nodeSfcKeys = domain.getLocalToGlobalSfcMap();
    const auto& conn_sfc = domain.getConnectivity();
    
    if (nodeSfcKeys.size() != nodeCount) {
        std::cerr << "Rank " << domain.rank() << ": ERROR - nodeSfcKeys.size()=" << nodeSfcKeys.size() 
                  << " != nodeCount=" << nodeCount << std::endl;
        return;
    }

    // Get SFC assignment range for this rank
    KeyType assignmentStart = domain.getDomain().assignmentStart();
    
    int blockSize = 256;
    
    // Initialize all nodes as ghost (0)
    thrust::fill(thrust::device, 
                thrust::device_pointer_cast(d_nodeOwnership_.data()),
                thrust::device_pointer_cast(d_nodeOwnership_.data() + nodeCount), 
                0);
    
    // Phase 1: Gather all rank assignments for SFC-key-based ownership
    cstone::DeviceVector<KeyType> d_assignment(domain.numRanks() + 1);
    std::vector<KeyType> h_assignment(domain.numRanks() + 1);
    
    MPI_Datatype mpiKeyType = (sizeof(KeyType) == 8) ? MPI_UNSIGNED_LONG : MPI_UNSIGNED;
    MPI_Allgather(&assignmentStart, 1, mpiKeyType,
                  h_assignment.data(), 1, mpiKeyType,
                  MPI_COMM_WORLD);
    h_assignment[domain.numRanks()] = std::numeric_limits<KeyType>::max();
    
    thrust::copy(h_assignment.begin(), h_assignment.end(), 
                thrust::device_pointer_cast(d_assignment.data()));
    
    // Phase 2: Assign ownership based on SFC key ranges
    if (nodeCount > 0) {
        int numBlocks = (nodeCount + blockSize - 1) / blockSize;
        determineOwnershipDirectKernel<KeyType><<<numBlocks, blockSize>>>(
            thrust::raw_pointer_cast(d_nodeOwnership_.data()),
            thrust::raw_pointer_cast(nodeSfcKeys.data()),
            thrust::raw_pointer_cast(d_assignment.data()),
            domain.numRanks(),
            domain.rank(),
            nodeCount);
        cudaCheckError();
    }
    
    // Phase 3: Mark nodes in LOCAL elements to detect potential shared boundaries
    cstone::DeviceVector<unsigned int> d_nodeFlags(nodeCount, 0);
    
    constexpr int NodesPerElem = ElementTag::NodesPerElement;
    auto conn0 = std::get<0>(conn_sfc);
    auto conn1 = std::get<1>(conn_sfc);
    auto conn2 = std::get<2>(conn_sfc);
    auto conn3 = std::get<3>(conn_sfc);
    
    // For HexTag, get all 8 connectivity vectors; for TetTag, use nullptr placeholders
    const KeyType* conn4_ptr = nullptr;
    const KeyType* conn5_ptr = nullptr;
    const KeyType* conn6_ptr = nullptr;
    const KeyType* conn7_ptr = nullptr;
    if constexpr (NodesPerElem == 8) {
        conn4_ptr = thrust::raw_pointer_cast(std::get<4>(conn_sfc).data());
        conn5_ptr = thrust::raw_pointer_cast(std::get<5>(conn_sfc).data());
        conn6_ptr = thrust::raw_pointer_cast(std::get<6>(conn_sfc).data());
        conn7_ptr = thrust::raw_pointer_cast(std::get<7>(conn_sfc).data());
    }
    
    size_t localElementCount = domain.localElementCount();
    
    if (localElementCount > 0) {
        int numBlocks = (localElementCount + blockSize - 1) / blockSize;
        markNodesInElementRangeKernel<KeyType, NodesPerElem><<<numBlocks, blockSize>>>(
            thrust::raw_pointer_cast(d_nodeFlags.data()),
            thrust::raw_pointer_cast(conn0.data()),
            thrust::raw_pointer_cast(conn1.data()),
            thrust::raw_pointer_cast(conn2.data()),
            thrust::raw_pointer_cast(conn3.data()),
            conn4_ptr, conn5_ptr, conn6_ptr, conn7_ptr,
            thrust::raw_pointer_cast(nodeSfcKeys.data()),
            0, localElementCount, nodeCount, 0x01);
        cudaCheckError();
    }
    
    // Phase 4: Mark nodes in HALO elements
    size_t totalElementCount = domain.getElementCount();
    
    if (totalElementCount > localElementCount) {
        int numBlocks = (totalElementCount - localElementCount + blockSize - 1) / blockSize;
        markNodesInElementRangeKernel<KeyType, NodesPerElem><<<numBlocks, blockSize>>>(
            thrust::raw_pointer_cast(d_nodeFlags.data()),
            thrust::raw_pointer_cast(conn0.data()),
            thrust::raw_pointer_cast(conn1.data()),
            thrust::raw_pointer_cast(conn2.data()),
            thrust::raw_pointer_cast(conn3.data()),
            conn4_ptr, conn5_ptr, conn6_ptr, conn7_ptr,
            thrust::raw_pointer_cast(nodeSfcKeys.data()),
            localElementCount, totalElementCount, nodeCount, 0x02);
        cudaCheckError();
    }
    
    // Phase 5: Upgrade owned nodes (state 1) to shared (state 2) if they appear in both local and halo
    if (nodeCount > 0) {
        int numBlocks = (nodeCount + blockSize - 1) / blockSize;
        detectSharedNodesKernel<KeyType><<<numBlocks, blockSize>>>(
            thrust::raw_pointer_cast(d_nodeOwnership_.data()),
            thrust::raw_pointer_cast(d_nodeFlags.data()),
            nodeCount);
        cudaCheckError();
    }
    
    // Debug: Count ownership states and verify element-node associations
    if (domain.rank() == 0 || domain.rank() == 1) {
        thrust::host_vector<uint8_t> h_ownership(nodeCount);
        thrust::host_vector<unsigned int> h_flags(nodeCount);
        thrust::copy(thrust::device_pointer_cast(d_nodeOwnership_.data()),
                    thrust::device_pointer_cast(d_nodeOwnership_.data() + nodeCount),
                    h_ownership.begin());
        thrust::copy(thrust::device_pointer_cast(d_nodeFlags.data()),
                    thrust::device_pointer_cast(d_nodeFlags.data() + nodeCount),
                    h_flags.begin());
        
        size_t cnt_ghost = 0, cnt_owned = 0, cnt_shared = 0;
        size_t cnt_in_local = 0, cnt_in_halo = 0, cnt_in_both = 0, cnt_in_neither = 0;
        size_t cnt_owned_not_in_local = 0;  // Problematic: owned but not in local elements
        for (size_t i = 0; i < nodeCount; ++i) {
            if (h_ownership[i] == 0) cnt_ghost++;
            else if (h_ownership[i] == 1) cnt_owned++;
            else if (h_ownership[i] == 2) cnt_shared++;
            
            if ((h_flags[i] & 0x01) && (h_flags[i] & 0x02)) cnt_in_both++;
            else if (h_flags[i] & 0x01) cnt_in_local++;
            else if (h_flags[i] & 0x02) cnt_in_halo++;
            else cnt_in_neither++;
            
            // Check for owned nodes not in local elements (will have zero diagonal)
            if ((h_ownership[i] == 1 || h_ownership[i] == 2) && !(h_flags[i] & 0x01)) {
                cnt_owned_not_in_local++;
            }
        }
        std::cout << "Rank " << domain.rank() << ": Ownership states after 3-phase - " 
                  << cnt_ghost << " ghost, " << cnt_owned << " pure owned, " 
                  << cnt_shared << " shared\n";
        std::cout << "Rank " << domain.rank() << ": Element flags - " 
                  << cnt_in_local << " local only, " << cnt_in_halo << " halo only, "
                  << cnt_in_both << " both, " << cnt_in_neither << " neither\n";
        if (cnt_owned_not_in_local > 0) {
            std::cout << "Rank " << domain.rank() << ": WARNING - " << cnt_owned_not_in_local 
                      << " owned nodes NOT in local elements (will have zero diagonals!)\n";
        }
    }
}

// ===== For unsigned KeyType =====
// Float combinations
template __global__ void computeCharacteristicSizesKernel<TetTag, unsigned, float>(
    const float* x, const float* y, const float* z,
    const unsigned* indices0, const unsigned* indices1, const unsigned* indices2, const unsigned* indices3,
    const unsigned* indices4, const unsigned* indices5, const unsigned* indices6, const unsigned* indices7,
    int* nodeTetCount, float* h, int numElements);

template __global__ void extractRepCoordinatesKernel<TetTag, unsigned, float>(
    const float* x, const float* y, const float* z, const float* h,
    const unsigned* elemToNodeMap, float* elemX, float* elemY, float* elemZ, float* elemH, int numElements);

template __global__ void findRepresentativeNodesKernel<TetTag, unsigned, float>(
    const unsigned* indices0, const unsigned* indices1, const unsigned* indices2, const unsigned* indices3,
    const unsigned* indices4, const unsigned* indices5, const unsigned* indices6, const unsigned* indices7,
    const unsigned* sfcCodes, unsigned* elemToNodeMap, int numElements);

template __global__ void finalizeCharacteristicSizesKernel<unsigned, float>(float* h, int* nodeTetCount, int numNodes);

// Double combinations
template __global__ void computeCharacteristicSizesKernel<TetTag, unsigned, double>(
    const double* x, const double* y, const double* z,
    const unsigned* indices0, const unsigned* indices1, const unsigned* indices2, const unsigned* indices3,
    const unsigned* indices4, const unsigned* indices5, const unsigned* indices6, const unsigned* indices7,
    int* nodeTetCount, double* h, int numElements);

template __global__ void extractRepCoordinatesKernel<TetTag, unsigned, double>(
    const double* x, const double* y, const double* z, const double* h,
    const unsigned* elemToNodeMap, double* elemX, double* elemY, double* elemZ, double* elemH, int numElements);

template __global__ void findRepresentativeNodesKernel<TetTag, unsigned, double>(
    const unsigned* indices0, const unsigned* indices1, const unsigned* indices2, const unsigned* indices3,
    const unsigned* indices4, const unsigned* indices5, const unsigned* indices6, const unsigned* indices7,
    const unsigned* sfcCodes, unsigned* elemToNodeMap, int numElements);

template __global__ void finalizeCharacteristicSizesKernel<unsigned, double>(double* h, int* nodeTetCount, int numNodes);

// ===== For uint64_t KeyType =====
// Float combinations
template __global__ void computeCharacteristicSizesKernel<TetTag, uint64_t, float>(
    const float* x, const float* y, const float* z,
    const uint64_t* indices0, const uint64_t* indices1, const uint64_t* indices2, const uint64_t* indices3,
    const uint64_t* indices4, const uint64_t* indices5, const uint64_t* indices6, const uint64_t* indices7,
    int* nodeTetCount, float* h, int numElements);

template __global__ void extractRepCoordinatesKernel<TetTag, uint64_t, float>(
    const float* x, const float* y, const float* z, const float* h,
    const uint64_t* elemToNodeMap, float* elemX, float* elemY, float* elemZ, float* elemH, int numElements);

template __global__ void findRepresentativeNodesKernel<TetTag, uint64_t, float>(
    const uint64_t* indices0, const uint64_t* indices1, const uint64_t* indices2, const uint64_t* indices3,
    const uint64_t* indices4, const uint64_t* indices5, const uint64_t* indices6, const uint64_t* indices7,
    const uint64_t* sfcCodes, uint64_t* elemToNodeMap, int numElements);

template __global__ void finalizeCharacteristicSizesKernel<uint64_t, float>(float* h, int* nodeTetCount, int numNodes);

// Double combinations
template __global__ void computeCharacteristicSizesKernel<TetTag, uint64_t, double>(
    const double* x, const double* y, const double* z,
    const uint64_t* indices0, const uint64_t* indices1, const uint64_t* indices2, const uint64_t* indices3,
    const uint64_t* indices4, const uint64_t* indices5, const uint64_t* indices6, const uint64_t* indices7,
    int* nodeTetCount, double* h, int numElements);

template __global__ void extractRepCoordinatesKernel<TetTag, uint64_t, double>(
    const double* x, const double* y, const double* z, const double* h,
    const uint64_t* elemToNodeMap, double* elemX, double* elemY, double* elemZ, double* elemH, int numElements);

template __global__ void findRepresentativeNodesKernel<TetTag, uint64_t, double>(
    const uint64_t* indices0, const uint64_t* indices1, const uint64_t* indices2, const uint64_t* indices3,
    const uint64_t* indices4, const uint64_t* indices5, const uint64_t* indices6, const uint64_t* indices7,
    const uint64_t* sfcCodes, uint64_t* elemToNodeMap, int numElements);

template __global__ void finalizeCharacteristicSizesKernel<uint64_t, double>(double* h, int* nodeTetCount, int numNodes);

// ===== Hex8 instantiations =====
// Float combinations
template __global__ void computeCharacteristicSizesKernel<HexTag, unsigned, float>(
    const float* x, const float* y, const float* z,
    const unsigned* indices0, const unsigned* indices1, const unsigned* indices2, const unsigned* indices3,
    const unsigned* indices4, const unsigned* indices5, const unsigned* indices6, const unsigned* indices7,
    int* nodeTetCount, float* h, int numElements);

template __global__ void findRepresentativeNodesKernel<HexTag, unsigned, float>(
    const unsigned* indices0, const unsigned* indices1, const unsigned* indices2, const unsigned* indices3,
    const unsigned* indices4, const unsigned* indices5, const unsigned* indices6, const unsigned* indices7,
    const unsigned* sfcCodes, unsigned* elemToNodeMap, int numElements);

template __global__ void buildSfcConnectivity<HexTag, unsigned, float>(
    const unsigned* indices0, const unsigned* indices1, const unsigned* indices2, const unsigned* indices3,
    const unsigned* indices4, const unsigned* indices5, const unsigned* indices6, const unsigned* indices7,
    const unsigned* sfcCodes, unsigned* keys0, unsigned* keys1, unsigned* keys2, unsigned* keys3,
    unsigned* keys4, unsigned* keys5, unsigned* keys6, unsigned* keys7, int numElements);

// Double combinations
template __global__ void computeCharacteristicSizesKernel<HexTag, unsigned, double>(
    const double* x, const double* y, const double* z,
    const unsigned* indices0, const unsigned* indices1, const unsigned* indices2, const unsigned* indices3,
    const unsigned* indices4, const unsigned* indices5, const unsigned* indices6, const unsigned* indices7,
    int* nodeTetCount, double* h, int numElements);

template __global__ void findRepresentativeNodesKernel<HexTag, unsigned, double>(
    const unsigned* indices0, const unsigned* indices1, const unsigned* indices2, const unsigned* indices3,
    const unsigned* indices4, const unsigned* indices5, const unsigned* indices6, const unsigned* indices7,
    const unsigned* sfcCodes, unsigned* elemToNodeMap, int numElements);

template __global__ void buildSfcConnectivity<HexTag, unsigned, double>(
    const unsigned* indices0, const unsigned* indices1, const unsigned* indices2, const unsigned* indices3,
    const unsigned* indices4, const unsigned* indices5, const unsigned* indices6, const unsigned* indices7,
    const unsigned* sfcCodes, unsigned* keys0, unsigned* keys1, unsigned* keys2, unsigned* keys3,
    unsigned* keys4, unsigned* keys5, unsigned* keys6, unsigned* keys7, int numElements);

// uint64_t combinations
template __global__ void computeCharacteristicSizesKernel<HexTag, uint64_t, float>(
    const float* x, const float* y, const float* z,
    const uint64_t* indices0, const uint64_t* indices1, const uint64_t* indices2, const uint64_t* indices3,
    const uint64_t* indices4, const uint64_t* indices5, const uint64_t* indices6, const uint64_t* indices7,
    int* nodeTetCount, float* h, int numElements);

template __global__ void findRepresentativeNodesKernel<HexTag, uint64_t, float>(
    const uint64_t* indices0, const uint64_t* indices1, const uint64_t* indices2, const uint64_t* indices3,
    const uint64_t* indices4, const uint64_t* indices5, const uint64_t* indices6, const uint64_t* indices7,
    const uint64_t* sfcCodes, uint64_t* elemToNodeMap, int numElements);

template __global__ void buildSfcConnectivity<HexTag, uint64_t, float>(
    const uint64_t* indices0, const uint64_t* indices1, const uint64_t* indices2, const uint64_t* indices3,
    const uint64_t* indices4, const uint64_t* indices5, const uint64_t* indices6, const uint64_t* indices7,
    const uint64_t* sfcCodes, uint64_t* keys0, uint64_t* keys1, uint64_t* keys2, uint64_t* keys3,
    uint64_t* keys4, uint64_t* keys5, uint64_t* keys6, uint64_t* keys7, int numElements);

template __global__ void computeCharacteristicSizesKernel<HexTag, uint64_t, double>(
    const double* x, const double* y, const double* z,
    const uint64_t* indices0, const uint64_t* indices1, const uint64_t* indices2, const uint64_t* indices3,
    const uint64_t* indices4, const uint64_t* indices5, const uint64_t* indices6, const uint64_t* indices7,
    int* nodeTetCount, double* h, int numElements);

template __global__ void findRepresentativeNodesKernel<HexTag, uint64_t, double>(
    const uint64_t* indices0, const uint64_t* indices1, const uint64_t* indices2, const uint64_t* indices3,
    const uint64_t* indices4, const uint64_t* indices5, const uint64_t* indices6, const uint64_t* indices7,
    const uint64_t* sfcCodes, uint64_t* elemToNodeMap, int numElements);

template __global__ void buildSfcConnectivity<HexTag, uint64_t, double>(
    const uint64_t* indices0, const uint64_t* indices1, const uint64_t* indices2, const uint64_t* indices3,
    const uint64_t* indices4, const uint64_t* indices5, const uint64_t* indices6, const uint64_t* indices7,
    const uint64_t* sfcCodes, uint64_t* keys0, uint64_t* keys1, uint64_t* keys2, uint64_t* keys3,
    uint64_t* keys4, uint64_t* keys5, uint64_t* keys6, uint64_t* keys7, int numElements);

// Add extractRepCoordinatesKernel instantiations for HexTag
template __global__ void extractRepCoordinatesKernel<HexTag, unsigned, float>(
    const float* x, const float* y, const float* z, const float* h,
    const unsigned* elemToNodeMap, float* elemX, float* elemY, float* elemZ, float* elemH, int numElements);

template __global__ void extractRepCoordinatesKernel<HexTag, unsigned, double>(
    const double* x, const double* y, const double* z, const double* h,
    const unsigned* elemToNodeMap, double* elemX, double* elemY, double* elemZ, double* elemH, int numElements);

template __global__ void extractRepCoordinatesKernel<HexTag, uint64_t, float>(
    const float* x, const float* y, const float* z, const float* h,
    const uint64_t* elemToNodeMap, float* elemX, float* elemY, float* elemZ, float* elemH, int numElements);

template __global__ void extractRepCoordinatesKernel<HexTag, uint64_t, double>(
    const double* x, const double* y, const double* z, const double* h,
    const uint64_t* elemToNodeMap, double* elemX, double* elemY, double* elemZ, double* elemH, int numElements);

template __global__ void transformCharacteristicSizesKernel<float>(float* d_h, size_t size, float meshFactor, float minH, float maxH);
template __global__ void transformCharacteristicSizesKernel<double>(double* d_h, size_t size, double meshFactor, double minH, double maxH);
template __global__ void fillCharacteristicSizesKernel<float>(float* d_h, size_t size, float value);
template __global__ void fillCharacteristicSizesKernel<double>(double* d_h, size_t size, double value);

template __global__ void computeElementVolumesKernel<TetTag, float>(const float* x, const float* y, const float* z,
    const int* indices0, const int* indices1, const int* indices2, const int* indices3, float* volumes, int numElements);

template __global__ void computeElementVolumesKernel<TetTag, double>(const double* x, const double* y, const double* z,
    const int* indices0, const int* indices1, const int* indices2, const int* indices3, double* volumes, int numElements);

// Explicit instantiation for computeSfcKeysGpu with common combinations
template void generateSfcKeys<unsigned, float>(
    const float* x, const float* y, const float* z, unsigned* keys, size_t numKeys, const cstone::Box<float>& box);
template void generateSfcKeys<unsigned, double>(
    const double* x, const double* y, const double* z, unsigned* keys, size_t numKeys, const cstone::Box<double>& box);
template void generateSfcKeys<uint64_t, float>(
    const float* x, const float* y, const float* z, uint64_t* keys, size_t numKeys, const cstone::Box<float>& box);
template void generateSfcKeys<uint64_t, double>(
    const double* x, const double* y, const double* z, uint64_t* keys, size_t numKeys, const cstone::Box<double>& box);


// Instantiations for buildSfcConnectivity kernel
// For tet elements with unsigned int keys
template __global__ void buildSfcConnectivity<TetTag, unsigned int, float>
    (const unsigned int*, const unsigned int*, const unsigned int*, const unsigned int*,
     const unsigned int*, const unsigned int*, const unsigned int*, const unsigned int*,
     const unsigned int*, unsigned int*, unsigned int*, 
     unsigned int*, unsigned int*, unsigned int*, unsigned int*, unsigned int*, unsigned int*, int);

template __global__ void buildSfcConnectivity<TetTag, unsigned int, double>
    (const unsigned int*, const unsigned int*, const unsigned int*, const unsigned int*,
     const unsigned int*, const unsigned int*, const unsigned int*, const unsigned int*,
     const unsigned int*, unsigned int*, unsigned int*, 
     unsigned int*, unsigned int*, unsigned int*, unsigned int*, unsigned int*, unsigned int*, int);

// For tet elements with uint64_t keys
template __global__ void buildSfcConnectivity<TetTag, uint64_t, float>
    (const uint64_t*, const uint64_t*, const uint64_t*, const uint64_t*,
     const uint64_t*, const uint64_t*, const uint64_t*, const uint64_t*,
     const uint64_t*, uint64_t*, uint64_t*, 
     uint64_t*, uint64_t*, uint64_t*, uint64_t*, uint64_t*, uint64_t*, int);

template __global__ void buildSfcConnectivity<TetTag, uint64_t, double>
    (const uint64_t*, const uint64_t*, const uint64_t*, const uint64_t*,
     const uint64_t*, const uint64_t*, const uint64_t*, const uint64_t*,
     const uint64_t*, uint64_t*, uint64_t*, 
     uint64_t*, uint64_t*, uint64_t*, uint64_t*, uint64_t*, uint64_t*, int);

// Update instantiations
template __global__ void flattenConnectivityKernel<unsigned int, 3>(
    ConnPtrs<unsigned int, 3>, unsigned int*, size_t);

template __global__ void flattenConnectivityKernel<unsigned int, 4>(
    ConnPtrs<unsigned int, 4>, unsigned int*, size_t);

template __global__ void flattenConnectivityKernel<unsigned int, 8>(
    ConnPtrs<unsigned int, 8>, unsigned int*, size_t);

template __global__ void flattenConnectivityKernel<uint64_t, 3>(
    ConnPtrs<uint64_t, 3>, uint64_t*, size_t);

template __global__ void flattenConnectivityKernel<uint64_t, 4>(
    ConnPtrs<uint64_t, 4>, uint64_t*, size_t);

template __global__ void flattenConnectivityKernel<uint64_t, 8>(
    ConnPtrs<uint64_t, 8>, uint64_t*, size_t);

template __global__ void mapSfcToLocalIdKernel<unsigned int>(
    const unsigned int*, unsigned int*, const unsigned int*, size_t, size_t);

template __global__ void mapSfcToLocalIdKernel<uint64_t>(
    const uint64_t*, uint64_t*, const uint64_t*, size_t, size_t);

// Explicit instantiations for decodeSfcToPhysical
template __device__ __host__ std::tuple<float, float, float> decodeSfcToPhysical<unsigned, float>(unsigned, const cstone::Box<float>&);
template __device__ __host__ std::tuple<double, double, double> decodeSfcToPhysical<unsigned, double>(unsigned, const cstone::Box<double>&);
template __device__ __host__ std::tuple<float, float, float> decodeSfcToPhysical<uint64_t, float>(uint64_t, const cstone::Box<float>&);
template __device__ __host__ std::tuple<double, double, double> decodeSfcToPhysical<uint64_t, double>(uint64_t, const cstone::Box<double>&);

template __global__ void decodeAllNodesKernel<unsigned int, float>(
    const unsigned int*, float*, float*, float*, size_t, cstone::Box<float>);
template __global__ void decodeAllNodesKernel<unsigned int, double>(
    const unsigned int*, double*, double*, double*, size_t, cstone::Box<double>);
template __global__ void decodeAllNodesKernel<uint64_t, float>(
    const uint64_t*, float*, float*, float*, size_t, cstone::Box<float>);
template __global__ void decodeAllNodesKernel<uint64_t, double>(
    const uint64_t*, double*, double*, double*, size_t, cstone::Box<double>);

// Explicit instantiations for HaloData
template struct HaloData<TetTag, float, unsigned, cstone::GpuTag>;
template struct HaloData<TetTag, double, unsigned, cstone::GpuTag>;
template struct HaloData<TetTag, float, uint64_t, cstone::GpuTag>;
template struct HaloData<TetTag, double, uint64_t, cstone::GpuTag>;

template struct HaloData<HexTag, float, unsigned, cstone::GpuTag>;
template struct HaloData<HexTag, double, unsigned, cstone::GpuTag>;
template struct HaloData<HexTag, float, uint64_t, cstone::GpuTag>;
template struct HaloData<HexTag, double, uint64_t, cstone::GpuTag>;

} // namespace mars
