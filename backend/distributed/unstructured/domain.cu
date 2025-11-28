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
            // Hexahedral element implementation
            // This would need all 8 node indices and 12 edges
            // ...
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
            // For hexahedra - would need additional parameters for indices4-7
            // Implementation would be similar but with 8 nodes
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
                                     const KeyType* sfcCodes,
                                     KeyType* keys0,
                                     KeyType* keys1,
                                     KeyType* keys2,
                                     KeyType* keys3,
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
            // For hexahedra - would need additional parameters for indices4-7
            // Implementation would be similar but with 8 nodes
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

template<typename KeyType>
__global__ void markHaloNodesKernel(uint8_t* nodeOwnership,
                                    const KeyType* nodeToElementList,
                                    const KeyType* nodeToElementOffsets,
                                    size_t nodeCount,
                                    size_t localStartIdx,
                                    size_t localEndIdx)
{
    size_t nodeIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (nodeIdx >= nodeCount) return;
    
    // Get the range of elements that reference this node
    KeyType startOffset = nodeToElementOffsets[nodeIdx];
    KeyType endOffset = nodeToElementOffsets[nodeIdx + 1];
    
    // Check if ANY element referencing this node is local (owned)
    uint8_t hasLocalElement = 0;
    for (KeyType i = startOffset; i < endOffset; ++i) {
        KeyType elemIdx = nodeToElementList[i];
        
        // Element is local if it's in [localStartIdx, localEndIdx)
        if (elemIdx >= localStartIdx && elemIdx < localEndIdx) {
            hasLocalElement = 1;
            break;
        }
    }
    
    // Node is owned (local) if it has at least one local element
    // Otherwise, it's a halo node (appears only in halo elements)
    nodeOwnership[nodeIdx] = hasLocalElement;
}

// Phase 1: Initialize all nodes as halo (pure ghost)
__global__ void initNodeOwnershipKernel(uint8_t* nodeOwnership, size_t nodeCount) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < nodeCount) {
        nodeOwnership[idx] = 0;  // All ghost initially
    }
}

// Phase 2: Mark nodes of local elements as owned (optimized for tets: 4 nodes)
template<typename KeyType, int NodesPerElement>
__global__ void markLocalElementNodesKernel(uint8_t* nodeOwnership,
                                            const KeyType* conn0,
                                            const KeyType* conn1,
                                            const KeyType* conn2,
                                            const KeyType* conn3,
                                            size_t numLocalElements,
                                            size_t localStartIdx)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numLocalElements) return;
    
    size_t elemIdx = idx + localStartIdx;
    
    // Mark all 4 nodes of this local tetrahedron as owned (1 = owned)
    nodeOwnership[conn0[elemIdx]] = 1;
    nodeOwnership[conn1[elemIdx]] = 1;
    nodeOwnership[conn2[elemIdx]] = 1;
    nodeOwnership[conn3[elemIdx]] = 1;
}

// Phase 3: Mark nodes that appear in BOTH local and halo elements as shared (2)
// These are partition boundary nodes
template<typename KeyType, int NodesPerElement>
__global__ void markSharedNodesKernel(uint8_t* nodeOwnership,
                                      const KeyType* conn0,
                                      const KeyType* conn1,
                                      const KeyType* conn2,
                                      const KeyType* conn3,
                                      size_t numHaloElements,
                                      const KeyType* haloIndices)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numHaloElements) return;
    
    size_t elemIdx = haloIndices[idx];
    
    // Check all 4 nodes of this halo element
    // If a node is marked as owned (1), it's also in a halo element, so mark as shared (2)
    KeyType n0 = conn0[elemIdx];
    KeyType n1 = conn1[elemIdx];
    KeyType n2 = conn2[elemIdx];
    KeyType n3 = conn3[elemIdx];
    
    if (nodeOwnership[n0] == 1) nodeOwnership[n0] = 2;  // Owned + halo = shared
    if (nodeOwnership[n1] == 1) nodeOwnership[n1] = 2;
    if (nodeOwnership[n2] == 1) nodeOwnership[n2] = 2;
    if (nodeOwnership[n3] == 1) nodeOwnership[n3] = 2;
}

template __global__ void markLocalElementNodesKernel<unsigned int, 4>(
    uint8_t*, const unsigned int*, const unsigned int*, const unsigned int*, const unsigned int*, size_t, size_t);
template __global__ void markLocalElementNodesKernel<uint64_t, 4>(
    uint8_t*, const uint64_t*, const uint64_t*, const uint64_t*, const uint64_t*, size_t, size_t);

template __global__ void markSharedNodesKernel<unsigned int, 4>(
    uint8_t*, const unsigned int*, const unsigned int*, const unsigned int*, const unsigned int*, size_t, const unsigned int*);
template __global__ void markSharedNodesKernel<uint64_t, 4>(
    uint8_t*, const uint64_t*, const uint64_t*, const uint64_t*, const uint64_t*, size_t, const uint64_t*);

// ===== For unsigned KeyType =====
// Float combinations
template __global__ void computeCharacteristicSizesKernel<TetTag, unsigned, float>(
    const float* x, const float* y, const float* z,
    const unsigned* indices0, const unsigned* indices1, const unsigned* indices2, const unsigned* indices3,
    int* nodeTetCount, float* h, int numElements);

template __global__ void extractRepCoordinatesKernel<TetTag, unsigned, float>(
    const float* x, const float* y, const float* z, const float* h,
    const unsigned* elemToNodeMap, float* elemX, float* elemY, float* elemZ, float* elemH, int numElements);

template __global__ void findRepresentativeNodesKernel<TetTag, unsigned, float>(
    const unsigned* indices0, const unsigned* indices1, const unsigned* indices2, const unsigned* indices3,
    const unsigned* sfcCodes, unsigned* elemToNodeMap, int numElements);

template __global__ void finalizeCharacteristicSizesKernel<unsigned, float>(float* h, int* nodeTetCount, int numNodes);

// Double combinations
template __global__ void computeCharacteristicSizesKernel<TetTag, unsigned, double>(
    const double* x, const double* y, const double* z,
    const unsigned* indices0, const unsigned* indices1, const unsigned* indices2, const unsigned* indices3,
    int* nodeTetCount, double* h, int numElements);

template __global__ void extractRepCoordinatesKernel<TetTag, unsigned, double>(
    const double* x, const double* y, const double* z, const double* h,
    const unsigned* elemToNodeMap, double* elemX, double* elemY, double* elemZ, double* elemH, int numElements);

template __global__ void findRepresentativeNodesKernel<TetTag, unsigned, double>(
    const unsigned* indices0, const unsigned* indices1, const unsigned* indices2, const unsigned* indices3,
    const unsigned* sfcCodes, unsigned* elemToNodeMap, int numElements);

template __global__ void finalizeCharacteristicSizesKernel<unsigned, double>(double* h, int* nodeTetCount, int numNodes);

// ===== For uint64_t KeyType =====
// Float combinations
template __global__ void computeCharacteristicSizesKernel<TetTag, uint64_t, float>(
    const float* x, const float* y, const float* z,
    const uint64_t* indices0, const uint64_t* indices1, const uint64_t* indices2, const uint64_t* indices3,
    int* nodeTetCount, float* h, int numElements);

template __global__ void extractRepCoordinatesKernel<TetTag, uint64_t, float>(
    const float* x, const float* y, const float* z, const float* h,
    const uint64_t* elemToNodeMap, float* elemX, float* elemY, float* elemZ, float* elemH, int numElements);

template __global__ void findRepresentativeNodesKernel<TetTag, uint64_t, float>(
    const uint64_t* indices0, const uint64_t* indices1, const uint64_t* indices2, const uint64_t* indices3,
    const uint64_t* sfcCodes, uint64_t* elemToNodeMap, int numElements);

template __global__ void finalizeCharacteristicSizesKernel<uint64_t, float>(float* h, int* nodeTetCount, int numNodes);

// Double combinations
template __global__ void computeCharacteristicSizesKernel<TetTag, uint64_t, double>(
    const double* x, const double* y, const double* z,
    const uint64_t* indices0, const uint64_t* indices1, const uint64_t* indices2, const uint64_t* indices3,
    int* nodeTetCount, double* h, int numElements);

template __global__ void extractRepCoordinatesKernel<TetTag, uint64_t, double>(
    const double* x, const double* y, const double* z, const double* h,
    const uint64_t* elemToNodeMap, double* elemX, double* elemY, double* elemZ, double* elemH, int numElements);

template __global__ void findRepresentativeNodesKernel<TetTag, uint64_t, double>(
    const uint64_t* indices0, const uint64_t* indices1, const uint64_t* indices2, const uint64_t* indices3,
    const uint64_t* sfcCodes, uint64_t* elemToNodeMap, int numElements);

template __global__ void finalizeCharacteristicSizesKernel<uint64_t, double>(double* h, int* nodeTetCount, int numNodes);

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
     const unsigned int*, unsigned int*, unsigned int*, 
     unsigned int*, unsigned int*, int);

template __global__ void buildSfcConnectivity<TetTag, unsigned int, double>
    (const unsigned int*, const unsigned int*, const unsigned int*, const unsigned int*,
     const unsigned int*, unsigned int*, unsigned int*, 
     unsigned int*, unsigned int*, int);

// For tet elements with uint64_t keys
template __global__ void buildSfcConnectivity<TetTag, uint64_t, float>
    (const uint64_t*, const uint64_t*, const uint64_t*, const uint64_t*,
     const uint64_t*, uint64_t*, uint64_t*, 
     uint64_t*, uint64_t*, int);

template __global__ void buildSfcConnectivity<TetTag, uint64_t, double>
    (const uint64_t*, const uint64_t*, const uint64_t*, const uint64_t*,
     const uint64_t*, uint64_t*, uint64_t*, 
     uint64_t*, uint64_t*, int);

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

template __global__ void markHaloNodesKernel<unsigned int>(
    uint8_t*, const unsigned int*, const unsigned int*, size_t, size_t, size_t);
template __global__ void markHaloNodesKernel<uint64_t>(
    uint8_t*, const uint64_t*, const uint64_t*, size_t, size_t, size_t);

} // namespace mars