#include "domain.hpp"
#include "thrust/sort.h"
#include "thrust/unique.h"
#include "thrust/binary_search.h"
#include "thrust/device_vector.h"
#include "cub/cub.cuh"
#include <unordered_map>
#include <algorithm>
#include <cstdlib>
#include <climits>
#include <map>
#include <string>

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

// Kernel to extract original coordinates of all 8 nodes for each hex8 element
// This extracts 24 per-element properties (8 nodes × 3 coords)
template<typename ElementTag, typename KeyType, typename RealType>
__global__ void extractElementNodeCoordsKernel(const RealType* x,
                                               const RealType* y,
                                               const RealType* z,
                                               const KeyType* idx0,
                                               const KeyType* idx1,
                                               const KeyType* idx2,
                                               const KeyType* idx3,
                                               const KeyType* idx4,
                                               const KeyType* idx5,
                                               const KeyType* idx6,
                                               const KeyType* idx7,
                                               RealType* elemOrigX0, RealType* elemOrigY0, RealType* elemOrigZ0,
                                               RealType* elemOrigX1, RealType* elemOrigY1, RealType* elemOrigZ1,
                                               RealType* elemOrigX2, RealType* elemOrigY2, RealType* elemOrigZ2,
                                               RealType* elemOrigX3, RealType* elemOrigY3, RealType* elemOrigZ3,
                                               RealType* elemOrigX4, RealType* elemOrigY4, RealType* elemOrigZ4,
                                               RealType* elemOrigX5, RealType* elemOrigY5, RealType* elemOrigZ5,
                                               RealType* elemOrigX6, RealType* elemOrigY6, RealType* elemOrigZ6,
                                               RealType* elemOrigX7, RealType* elemOrigY7, RealType* elemOrigZ7,
                                               int numElements)
{
    int elemIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (elemIdx < numElements)
    {
        // Extract coordinates for all 8 nodes of this element
        elemOrigX0[elemIdx] = x[idx0[elemIdx]];
        elemOrigY0[elemIdx] = y[idx0[elemIdx]];
        elemOrigZ0[elemIdx] = z[idx0[elemIdx]];

        elemOrigX1[elemIdx] = x[idx1[elemIdx]];
        elemOrigY1[elemIdx] = y[idx1[elemIdx]];
        elemOrigZ1[elemIdx] = z[idx1[elemIdx]];

        elemOrigX2[elemIdx] = x[idx2[elemIdx]];
        elemOrigY2[elemIdx] = y[idx2[elemIdx]];
        elemOrigZ2[elemIdx] = z[idx2[elemIdx]];

        elemOrigX3[elemIdx] = x[idx3[elemIdx]];
        elemOrigY3[elemIdx] = y[idx3[elemIdx]];
        elemOrigZ3[elemIdx] = z[idx3[elemIdx]];

        elemOrigX4[elemIdx] = x[idx4[elemIdx]];
        elemOrigY4[elemIdx] = y[idx4[elemIdx]];
        elemOrigZ4[elemIdx] = z[idx4[elemIdx]];

        elemOrigX5[elemIdx] = x[idx5[elemIdx]];
        elemOrigY5[elemIdx] = y[idx5[elemIdx]];
        elemOrigZ5[elemIdx] = z[idx5[elemIdx]];

        elemOrigX6[elemIdx] = x[idx6[elemIdx]];
        elemOrigY6[elemIdx] = y[idx6[elemIdx]];
        elemOrigZ6[elemIdx] = z[idx6[elemIdx]];

        elemOrigX7[elemIdx] = x[idx7[elemIdx]];
        elemOrigY7[elemIdx] = y[idx7[elemIdx]];
        elemOrigZ7[elemIdx] = z[idx7[elemIdx]];
    }
}

// Kernel to rebuild node coordinate arrays from element properties after sync
// Inverse of extractElementNodeCoordsKernel: elements → nodes
template<typename ElementTag, typename KeyType, typename RealType>
__global__ void rebuildNodeCoordsFromElementsKernel(
    const KeyType* connKey0, const KeyType* connKey1, const KeyType* connKey2, const KeyType* connKey3,
    const KeyType* connKey4, const KeyType* connKey5, const KeyType* connKey6, const KeyType* connKey7,
    const RealType* elemOrigX0, const RealType* elemOrigY0, const RealType* elemOrigZ0,
    const RealType* elemOrigX1, const RealType* elemOrigY1, const RealType* elemOrigZ1,
    const RealType* elemOrigX2, const RealType* elemOrigY2, const RealType* elemOrigZ2,
    const RealType* elemOrigX3, const RealType* elemOrigY3, const RealType* elemOrigZ3,
    const RealType* elemOrigX4, const RealType* elemOrigY4, const RealType* elemOrigZ4,
    const RealType* elemOrigX5, const RealType* elemOrigY5, const RealType* elemOrigZ5,
    const RealType* elemOrigX6, const RealType* elemOrigY6, const RealType* elemOrigZ6,
    const RealType* elemOrigX7, const RealType* elemOrigY7, const RealType* elemOrigZ7,
    RealType* nodeX, RealType* nodeY, RealType* nodeZ,
    int numElements)
{
    int elemIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (elemIdx < numElements)
    {
        // For each element, write the original coordinates to the node arrays
        // using the post-sync local node IDs from d_conn_keys_

        // Node 0
        KeyType nodeId0 = connKey0[elemIdx];
        nodeX[nodeId0] = elemOrigX0[elemIdx];
        nodeY[nodeId0] = elemOrigY0[elemIdx];
        nodeZ[nodeId0] = elemOrigZ0[elemIdx];

        // Node 1
        KeyType nodeId1 = connKey1[elemIdx];
        nodeX[nodeId1] = elemOrigX1[elemIdx];
        nodeY[nodeId1] = elemOrigY1[elemIdx];
        nodeZ[nodeId1] = elemOrigZ1[elemIdx];

        // Node 2
        KeyType nodeId2 = connKey2[elemIdx];
        nodeX[nodeId2] = elemOrigX2[elemIdx];
        nodeY[nodeId2] = elemOrigY2[elemIdx];
        nodeZ[nodeId2] = elemOrigZ2[elemIdx];

        // Node 3
        KeyType nodeId3 = connKey3[elemIdx];
        nodeX[nodeId3] = elemOrigX3[elemIdx];
        nodeY[nodeId3] = elemOrigY3[elemIdx];
        nodeZ[nodeId3] = elemOrigZ3[elemIdx];

        // Node 4
        KeyType nodeId4 = connKey4[elemIdx];
        nodeX[nodeId4] = elemOrigX4[elemIdx];
        nodeY[nodeId4] = elemOrigY4[elemIdx];
        nodeZ[nodeId4] = elemOrigZ4[elemIdx];

        // Node 5
        KeyType nodeId5 = connKey5[elemIdx];
        nodeX[nodeId5] = elemOrigX5[elemIdx];
        nodeY[nodeId5] = elemOrigY5[elemIdx];
        nodeZ[nodeId5] = elemOrigZ5[elemIdx];

        // Node 6
        KeyType nodeId6 = connKey6[elemIdx];
        nodeX[nodeId6] = elemOrigX6[elemIdx];
        nodeY[nodeId6] = elemOrigY6[elemIdx];
        nodeZ[nodeId6] = elemOrigZ6[elemIdx];

        // Node 7
        KeyType nodeId7 = connKey7[elemIdx];
        nodeX[nodeId7] = elemOrigX7[elemIdx];
        nodeY[nodeId7] = elemOrigY7[elemIdx];
        nodeZ[nodeId7] = elemOrigZ7[elemIdx];
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

// Generic kernel: mark nodes of elements with a given ownership value.
// connPtrs is a device array of NodesPerElem pointers to connectivity arrays.
template <typename KeyType, int NodesPerElem>
__global__ void markElementNodesOwnership(
    uint8_t* nodeOwnership,
    const KeyType* const* connPtrs,
    uint8_t value,
    size_t startElem, size_t numElems)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numElems) return;
    size_t e = startElem + idx;

    #pragma unroll
    for (int n = 0; n < NodesPerElem; ++n)
        nodeOwnership[connPtrs[n][e]] = value;
}


// Helper to extract raw pointers from a tuple of device vectors
template <typename KeyType, typename Tuple, std::size_t... Is>
void extractConnPtrs(const Tuple& t, const KeyType* ptrs[], std::index_sequence<Is...>)
{
    ((ptrs[Is] = thrust::raw_pointer_cast(std::get<Is>(t).data())), ...);
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

    constexpr int NPC = ElementTag::NodesPerElement;
    int blockSize = 256;
    
    // Initialize all nodes as ghost (0)
    thrust::fill(thrust::device, 
                thrust::device_pointer_cast(d_nodeOwnership_.data()),
                thrust::device_pointer_cast(d_nodeOwnership_.data() + nodeCount), 
                0);

    // Get element connectivity (dense local node IDs)
    const auto& d_conn = domain.getElementToNodeConnectivity();
    size_t startIdx = domain.startIndex();
    size_t endIdx   = domain.endIndex();

    // Build device pointer array for connectivity columns
    const KeyType* h_ptrs[NPC];
    extractConnPtrs<KeyType>(d_conn, h_ptrs, std::make_index_sequence<NPC>{});

    const KeyType** d_ptrs;
    cudaMalloc(&d_ptrs, NPC * sizeof(const KeyType*));
    cudaMemcpy(d_ptrs, h_ptrs, NPC * sizeof(const KeyType*), cudaMemcpyHostToDevice);

    // Step 1: Mark all nodes in LOCAL elements as owned
    size_t numLocal = endIdx - startIdx;
    if (numLocal > 0) {
        int numBlocks = (numLocal + blockSize - 1) / blockSize;
        markElementNodesOwnership<KeyType, NPC><<<numBlocks, blockSize>>>(
            thrust::raw_pointer_cast(d_nodeOwnership_.data()),
            d_ptrs, 1, startIdx, numLocal);
        cudaCheckError();
    }

    // Step 2: Yield boundary nodes to lower-SFC-key ranks.
    // Elements at indices [0, startIndex) are halo from lower-ranked partitions.
    // Those elements are LOCAL on lower-ranked ranks, so those ranks can fully
    // assemble the boundary nodes. Yielding ensures unique ownership: among all
    // ranks that have a node in their local elements, the lowest rank wins.
    if (startIdx > 0) {
        int numBlocks = (startIdx + blockSize - 1) / blockSize;
        markElementNodesOwnership<KeyType, NPC><<<numBlocks, blockSize>>>(
            thrust::raw_pointer_cast(d_nodeOwnership_.data()),
            d_ptrs, 0, 0, startIdx);
        cudaCheckError();
    }

    // Known limitation: a small number of corner-shared nodes (e.g. ~9 on
    // cube16/4-rank) end up doubly-owned because cstone's halo width does not
    // include the OPPOSITE-rank elements that touch corner-only contact nodes.
    // An attempted Step 3 (atomicMin claim across element halos) is ineffective
    // for the same reason: those elements aren't in any peer's halo set. A real
    // fix needs an MPI_Allgatherv of (sfc_key, owner_rank) pairs and a global
    // tiebreaker. Deferred until per-node halo topology is in place (#2).

    cudaFree(d_ptrs);
}

// ============================================================================
// NodeHaloTopology constructor implementation
// ============================================================================
//
// Builds per-peer send/recv lists of LOCAL node IDs for direct MPI exchange of
// node-DOF arrays. Two-stage protocol at init time:
//
// Stage 1: MPI_Allgatherv of (sfc_key, owner_rank) for OWNED nodes only.
//          Resulting global table has every (key, owner) pair across all ranks.
//          Apply lowest-rank-wins rule to detect & fix duplicates (the 9-node
//          ownership overlap from cstone halo-width gaps).
//
// Stage 2: For each ghost node N on this rank, look up its owner from the
//          global table -> append N's local id to recvNodeIds_[ownerRank].
//          For each owner R of mine, send R the list of N's sfc_keys I want
//          (Alltoallv). R receives the keys and resolves each to a local node
//          id via its sfc-map -> populates sendNodeIds_[requestingRank].
//
// Note: Stage 1 size = sum_r |owned_nodes_r| ~ globalNodes; tiny.
//       Stage 2 size = sum_r |ghost_nodes_r|; bounded by partition surface.

namespace {

// True when env var is set to a truthy value (1/true/on/yes).
inline bool envFlag(const char* name)
{
    const char* env = std::getenv(name);
    if (!env) return false;
    std::string s(env);
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s == "1" || s == "true" || s == "on" || s == "yes";
}

// v2 device-resident path is the default since cube1024/256 was validated
// bit-exact against host. MARS_NODEHALO_HOST=1 opts out (host O(global)
// fallback). MARS_NODEHALO_VALIDATE=1 runs both and aborts on mismatch.
inline bool nodeHaloHostFallbackEnabled() { return envFlag("MARS_NODEHALO_HOST"); }
inline bool nodeHaloValidateEnabled()     { return envFlag("MARS_NODEHALO_VALIDATE"); }

} // anonymous namespace

// Forward declaration of the host-driven path (the body that used to be the
// constructor). Now invocable both as the default and from validation.
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
static void buildNodeHaloTopologyHostPath(
    NodeHaloTopology<ElementTag, RealType, KeyType, AcceleratorTag>& topo,
    const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain);

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
NodeHaloTopology<ElementTag, RealType, KeyType, AcceleratorTag>::NodeHaloTopology(
    const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain)
{
    int rank     = domain.rank();
    int numRanks = domain.numRanks();

    if (numRanks == 1)
    {
        sendOffsets_.assign(1, 0);
        recvOffsets_.assign(1, 0);
        return;
    }

    // Runtime path selection. v2 is the default since cube1024/256 was
    // validated bit-exact (matrix=2.071203e+01, rhs=1.555497e-03 match host
    // reference; NodeHaloTopo phase 130 ms vs ~70 s host at cube512/32).
    // MARS_NODEHALO_HOST=1 opts out to the original O(global) Allgatherv path.
    // MARS_NODEHALO_VALIDATE=1 runs both and aborts on per-peer node-id mismatch.
    bool useHost  = nodeHaloHostFallbackEnabled();
    bool validate = nodeHaloValidateEnabled();

    if (validate)
    {
        bool ok = NodeHaloTopology::validateAgainstHost(domain);
        if (!ok)
        {
            std::cerr << "[Rank " << rank << "] NodeHaloTopology v2 vs host MISMATCH — aborting" << std::endl;
            std::cerr.flush();
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        this->buildFromCstoneHalos(domain);
    }
    else if (useHost) { buildNodeHaloTopologyHostPath(*this, domain); }
    else              { this->buildFromCstoneHalos(domain); }

    // Pre-size hot-path staging buffers (Gate 4 fix for v1 race).
    size_t sendTotal = sendOffsets_.empty() ? 0 : size_t(sendOffsets_.back());
    size_t recvTotal = recvOffsets_.empty() ? 0 : size_t(recvOffsets_.back());
    sendBuf_.resize(sendTotal);
    recvBuf_.resize(recvTotal);
}

// =====================================================================
// HOST FALLBACK PATH (NOT THE DEFAULT)
// =====================================================================
// O(global) construction using MPI_Allgatherv + MPI_Alltoallv. Active
// only when MARS_NODEHALO_HOST=1, or transiently inside
// MARS_NODEHALO_VALIDATE=1 for cross-checking.
//
// All MPI_Allgather*, MPI_Alltoall* in this file live in this function.
// The default path (buildFromCstoneHalos) uses only point-to-point
// MPI_Isend/Irecv in the peer subgroup — no collectives.
//
// Body is unchanged from commit 1d90501; refactored into a free
// function so the new constructor can call it as a fallback.
// =====================================================================
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
static void buildNodeHaloTopologyHostPath(
    NodeHaloTopology<ElementTag, RealType, KeyType, AcceleratorTag>& topo,
    const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain)
{
    using TopoT = NodeHaloTopology<ElementTag, RealType, KeyType, AcceleratorTag>;
    auto& peers_       = topo.peers_;
    auto& sendOffsets_ = topo.sendOffsets_;
    auto& recvOffsets_ = topo.recvOffsets_;
    auto& sendNodeIds_ = topo.sendNodeIds_;
    auto& recvNodeIds_ = topo.recvNodeIds_;

    int rank     = domain.rank();
    int numRanks = domain.numRanks();

    if (numRanks == 1)
    {
        sendOffsets_.assign(1, 0);
        recvOffsets_.assign(1, 0);
        return;
    }

    size_t nodeCount = domain.getNodeCount();

    // ----- Pull node ownership and SFC keys to host -----
    std::vector<uint8_t> h_owned(nodeCount);
    {
        const auto& d_own = domain.getNodeOwnershipMap();
        cudaMemcpy(h_owned.data(), thrust::raw_pointer_cast(d_own.data()),
                   nodeCount * sizeof(uint8_t), cudaMemcpyDeviceToHost);
    }
    std::vector<KeyType> h_sfc(nodeCount);
    {
        const auto& d_map = domain.getLocalToGlobalSfcMap();
        cudaMemcpy(h_sfc.data(), thrust::raw_pointer_cast(d_map.data()),
                   nodeCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
    }

    // ----- Build my local owned-keys list: (sfc, localNodeId) -----
    std::vector<KeyType> myOwnedKeys;
    std::vector<int>     myOwnedLocalIds;
    myOwnedKeys.reserve(nodeCount);
    myOwnedLocalIds.reserve(nodeCount);
    for (size_t n = 0; n < nodeCount; ++n)
    {
        if (h_owned[n] == 1)
        {
            myOwnedKeys.push_back(h_sfc[n]);
            myOwnedLocalIds.push_back(int(n));
        }
    }

    // ----- Stage 1: MPI_Allgatherv of OWNED sfc keys + ranks -----
    int myOwnedCount = int(myOwnedKeys.size());
    std::vector<int> ownedCounts(numRanks), ownedDispls(numRanks);
    MPI_Allgather(&myOwnedCount, 1, MPI_INT, ownedCounts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    int totalOwned = 0;
    for (int r = 0; r < numRanks; ++r) { ownedDispls[r] = totalOwned; totalOwned += ownedCounts[r]; }

    std::vector<KeyType> globalKeys(totalOwned);
    auto mpiKeyType = (sizeof(KeyType) == 8) ? MPI_UINT64_T : MPI_UNSIGNED;
    MPI_Allgatherv(myOwnedKeys.data(), myOwnedCount, mpiKeyType,
                   globalKeys.data(), ownedCounts.data(), ownedDispls.data(), mpiKeyType,
                   MPI_COMM_WORLD);

    // Build sfc -> first-claiming-rank map (lowest-rank wins for duplicates)
    std::unordered_map<KeyType, int> keyOwner;
    keyOwner.reserve(totalOwned);
    for (int r = 0; r < numRanks; ++r)
    {
        int beg = ownedDispls[r], end = beg + ownedCounts[r];
        for (int i = beg; i < end; ++i)
        {
            auto it = keyOwner.find(globalKeys[i]);
            if (it == keyOwner.end() || r < it->second) keyOwner[globalKeys[i]] = r;
        }
    }

    // Fix duplicates on this rank: if keyOwner[myKey] != myrank, yield ownership
    int yielded = 0;
    {
        auto& d_own = const_cast<typename HaloData<ElementTag, RealType, KeyType, AcceleratorTag>::template DeviceVector<uint8_t>&>(
            domain.halo_->d_nodeOwnership_);
        std::vector<uint8_t> h_ownNew = h_owned;
        for (size_t n = 0; n < nodeCount; ++n)
        {
            if (h_ownNew[n] == 1)
            {
                auto it = keyOwner.find(h_sfc[n]);
                if (it != keyOwner.end() && it->second != rank)
                {
                    h_ownNew[n] = 0;  // yield
                    ++yielded;
                }
            }
        }
        if (yielded > 0)
        {
            cudaMemcpy(thrust::raw_pointer_cast(d_own.data()), h_ownNew.data(),
                       nodeCount * sizeof(uint8_t), cudaMemcpyHostToDevice);
            h_owned = std::move(h_ownNew);
        }
        std::cout << "Rank " << rank << ": NodeHaloTopo yielded " << yielded
                  << " duplicate-owned nodes" << std::endl;
        std::cout.flush();
    }

    // ----- Stage 2: build recv lists (my ghosts grouped by owner) -----
    // Per-peer: keys I want (to send as request) + my local node ids (to populate recvNodeIds_)
    std::vector<std::vector<KeyType>> reqKeysPerPeer(numRanks);
    std::vector<std::vector<int>>     reqLocalIdsPerPeer(numRanks);
    for (size_t n = 0; n < nodeCount; ++n)
    {
        if (h_owned[n] == 0)
        {
            auto it = keyOwner.find(h_sfc[n]);
            if (it == keyOwner.end()) continue;  // unowned globally? skip
            int owner = it->second;
            reqKeysPerPeer[owner].push_back(h_sfc[n]);
            reqLocalIdsPerPeer[owner].push_back(int(n));
        }
    }

    // Stage 2a: Alltoall to learn how many keys each peer wants from me
    std::vector<int> sendCounts(numRanks), recvCounts(numRanks);
    for (int p = 0; p < numRanks; ++p) sendCounts[p] = int(reqKeysPerPeer[p].size());
    MPI_Alltoall(sendCounts.data(), 1, MPI_INT, recvCounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    std::vector<int> sendDispls(numRanks), recvDispls(numRanks);
    int totalReqSend = 0, totalReqRecv = 0;
    for (int p = 0; p < numRanks; ++p)
    {
        sendDispls[p] = totalReqSend; totalReqSend += sendCounts[p];
        recvDispls[p] = totalReqRecv; totalReqRecv += recvCounts[p];
    }

    // Flatten request keys to send
    std::vector<KeyType> reqKeysFlat(totalReqSend);
    for (int p = 0; p < numRanks; ++p)
        std::copy(reqKeysPerPeer[p].begin(), reqKeysPerPeer[p].end(),
                  reqKeysFlat.begin() + sendDispls[p]);

    // Stage 2b: Alltoallv to exchange the actual sfc keys
    std::vector<KeyType> recvReqKeys(totalReqRecv);
    MPI_Alltoallv(reqKeysFlat.data(),  sendCounts.data(), sendDispls.data(), mpiKeyType,
                  recvReqKeys.data(),  recvCounts.data(), recvDispls.data(), mpiKeyType,
                  MPI_COMM_WORLD);

    // Build sfc -> myLocalId map for fast lookup
    std::unordered_map<KeyType, int> myKeyToLocal;
    myKeyToLocal.reserve(nodeCount);
    for (size_t n = 0; n < nodeCount; ++n) myKeyToLocal[h_sfc[n]] = int(n);

    // Resolve received-keys -> local node ids; group by requester
    // recvReqKeys is already grouped by source rank via recvDispls
    std::vector<std::vector<int>> sendIdsPerPeer(numRanks);
    for (int p = 0; p < numRanks; ++p)
    {
        sendIdsPerPeer[p].reserve(recvCounts[p]);
        for (int i = 0; i < recvCounts[p]; ++i)
        {
            auto it = myKeyToLocal.find(recvReqKeys[recvDispls[p] + i]);
            // Owner side MUST find the key locally; if not, MPI bookkeeping bug.
            sendIdsPerPeer[p].push_back(it != myKeyToLocal.end() ? it->second : -1);
        }
    }

    // ----- Compact into peer-list / CSR structure -----
    peers_.clear();
    sendOffsets_.assign(1, 0);
    recvOffsets_.assign(1, 0);
    std::vector<int> sendNodesHost, recvNodesHost;

    for (int p = 0; p < numRanks; ++p)
    {
        int sCnt = int(sendIdsPerPeer[p].size());
        int rCnt = int(reqLocalIdsPerPeer[p].size());
        if (sCnt == 0 && rCnt == 0) continue;
        peers_.push_back(p);
        for (int v : sendIdsPerPeer[p]) sendNodesHost.push_back(v);
        for (int v : reqLocalIdsPerPeer[p]) recvNodesHost.push_back(v);
        sendOffsets_.push_back(int(sendNodesHost.size()));
        recvOffsets_.push_back(int(recvNodesHost.size()));
    }

    // Upload node-id arrays to device
    sendNodeIds_.resize(sendNodesHost.size());
    recvNodeIds_.resize(recvNodesHost.size());
    if (!sendNodesHost.empty())
        cudaMemcpy(thrust::raw_pointer_cast(sendNodeIds_.data()), sendNodesHost.data(),
                   sendNodesHost.size() * sizeof(int), cudaMemcpyHostToDevice);
    if (!recvNodesHost.empty())
        cudaMemcpy(thrust::raw_pointer_cast(recvNodeIds_.data()), recvNodesHost.data(),
                   recvNodesHost.size() * sizeof(int), cudaMemcpyHostToDevice);

    std::cout << "Rank " << rank << ": NodeHaloTopo " << peers_.size() << " peers, "
              << sendNodesHost.size() << " send nodes, "
              << recvNodesHost.size() << " recv nodes" << std::endl;
    std::cout.flush();
}

// ============================================================================
// Gate 1-4: Device-resident NodeHaloTopology construction
//
// Builds per-peer (sendNodeIds_, recvNodeIds_) CSRs by consuming cstone's
// element-halo bookkeeping (incomingHaloIndices, outgoingHaloIndices) and
// filtering by node ownership. O(local boundary), no MPI_Allgatherv.
//
// Gate 1 contract: same fields populated as buildNodeHaloTopologyHostPath,
// just via a different algorithm. Validation harness (validateAgainstHost)
// diffs per-peer sorted node-id lists; bit-exact match expected.
//
// UNTESTED on hardware. Build, then on the cluster:
//   MARS_NODEHALO_VALIDATE=1 mpirun ... ./mars_cvfem_graph ...
//      → runs both paths and aborts on first mismatch
//   MARS_NODEHALO_V2=1 mpirun ... ./mars_cvfem_graph ...
//      → runs only v2 path (use after Gate 1 passes)
//   (no env)            → host O(global) path (default, current production)
//
// KNOWN GAP: this builder consumes d_nodeOwnership_ as produced by
// HaloData::buildNodeOwnership() Steps 1+2 only. It does NOT run a
// global tiebreaker exchange to resolve corner-shared duplicate ownership
// (the host path does that via MPI_Allgatherv at L1085 and the keyOwner
// map at L1090).
//
// Consequence: on cube/structured meshes where Steps 1+2 already produce
// unambiguous ownership, v2 should match host bit-exactly. On meshes
// where multiple ranks initially own the same SFC key (corner-shared
// nodes on >=4-rank partitions), v2 will diverge from host. Gate 1's
// validateAgainstHost is precisely how we detect this — if it fires,
// we add a tiebreaker pass before flipping v2 on by default.
// ============================================================================

// Send-side per-node kernel: for each node N where I'm the authoritative
// owner, emit (claimant_rank, N) for every other rank that claimed N's
// SFC key (i.e. every peer that touches N and therefore needs it as a
// ghost). Using the merged peer-claim CSR avoids relying on cstone's
// outgoing-element-halo to peer P, which can miss corner-only contacts.
template<typename KeyType>
__global__ void emitSendNodesKernel(
    const KeyType* d_localToGlobalSfc,    // [nodeCount]
    const int*     d_authoritativeOwner,  // [nodeCount]
    const int*     d_claimRanks,          // [numClaims], sorted by key
    const int*     d_claimKeyOffsets,     // [numUniqueKeys+1]
    const KeyType* d_uniqueClaimKeys,     // [numUniqueKeys] sorted
    int            numUniqueKeys,
    int            myRank,
    int            outStride,             // bound on number of claimants per node
    int*           d_outPeer,             // [nodeCount * outStride]
    int*           d_outNodeId,           // [nodeCount * outStride]
    int            nodeCount)
{
    int n = blockIdx.x * blockDim.x + threadIdx.x;
    if (n >= nodeCount) return;

    int outBase = n * outStride;
    // Default-fill with sentinels.
    for (int s = 0; s < outStride; ++s)
    {
        d_outPeer[outBase + s]   = -1;
        d_outNodeId[outBase + s] = -1;
    }

    if (d_authoritativeOwner[n] != myRank) return;  // I don't own → don't send

    // Look up my key in unique-claims.
    KeyType myKey = d_localToGlobalSfc[n];
    int lo = 0, hi = numUniqueKeys;
    while (lo < hi)
    {
        int mid = lo + (hi - lo) / 2;
        if (d_uniqueClaimKeys[mid] < myKey) lo = mid + 1;
        else                                 hi = mid;
    }
    if (lo >= numUniqueKeys || d_uniqueClaimKeys[lo] != myKey) return;

    int beg = d_claimKeyOffsets[lo];
    int end = d_claimKeyOffsets[lo + 1];
    int slot = 0;
    constexpr int TOUCH_FLAG_LOCAL = (int)0x80000000;
    for (int i = beg; i < end && slot < outStride; ++i)
    {
        int rLabel = d_claimRanks[i];
        int r = rLabel & ~TOUCH_FLAG_LOCAL;  // strip touch flag
        if (r != myRank)  // skip self
        {
            d_outPeer[outBase + slot]   = r;
            d_outNodeId[outBase + slot] = n;
            ++slot;
        }
    }
}

// Recv-side per-node kernel: for each node N where I'm not the owner,
// emit (owner_rank, N) so I receive N from owner.
__global__ void emitRecvNodesKernel(
    const int* d_authoritativeOwner,      // [nodeCount]
    int        myRank,
    int*       d_outPeer,                  // [nodeCount]
    int*       d_outNodeId,                // [nodeCount]
    int        nodeCount)
{
    int n = blockIdx.x * blockDim.x + threadIdx.x;
    if (n >= nodeCount) return;

    int owner = d_authoritativeOwner[n];
    bool keep = (owner != myRank) && (owner >= 0);
    d_outPeer[n]   = keep ? owner : -1;
    d_outNodeId[n] = keep ? n     : -1;
}

// Per-node kernel: collect the SFC keys of all "boundary candidate" nodes —
// nodes that appear in any halo element corner (in or out). Boundary
// candidates are exactly the nodes whose ownership might be contested with
// peers; non-candidates are interior nodes owned exclusively by us.
template<typename KeyType, int NPC>
__global__ void markBoundaryCandidatesKernel(
    const KeyType* const* d_conn_local_ids,
    const int*     d_rangeStart,
    const int*     d_rangeEnd,
    int            numRanges,
    uint8_t*       d_isBoundary,               // [nodeCount]
    int            totalElems)
{
    int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= totalElems) return;

    int elemIdx = -1;
    int cumulative = 0;
    for (int r = 0; r < numRanges; ++r)
    {
        int len = d_rangeEnd[r] - d_rangeStart[r];
        if (e < cumulative + len)
        {
            elemIdx = d_rangeStart[r] + (e - cumulative);
            break;
        }
        cumulative += len;
    }
    if (elemIdx < 0) return;

#pragma unroll
    for (int c = 0; c < NPC; ++c)
    {
        int n = static_cast<int>(d_conn_local_ids[c][elemIdx]);
        if (n >= 0) d_isBoundary[n] = 1;
    }
}

// Per-element kernel: mark every node touched by any LOCAL element [startIdx,endIdx)
// as claimable by this rank. This is the "I claim because I have a local
// element with this corner" flag — the right starting point for the
// tiebreaker, before Step 2's halo-direction bias is applied.
template<typename KeyType, int NPC>
__global__ void markLocallyClaimedKernel(
    const KeyType* const* d_conn_local_ids,
    size_t startIdx, size_t endIdx,
    uint8_t* d_locallyClaimed)        // [nodeCount]
{
    size_t off = blockIdx.x * blockDim.x + threadIdx.x;
    size_t numLocal = endIdx - startIdx;
    if (off >= numLocal) return;
    size_t e = startIdx + off;

#pragma unroll
    for (int c = 0; c < NPC; ++c)
    {
        int n = static_cast<int>(d_conn_local_ids[c][e]);
        if (n >= 0) d_locallyClaimed[n] = 1;
    }
}

// Reduce per-node SFC key claims using min-rank-wins. For each boundary
// node N: the authoritative owner is the lowest-ranked rank that claimed
// the node's SFC key. "Claim" includes my own only if I currently flag
// the node as owned (d_nodeOwnership ∈ {1,2}); otherwise I'm just an
// observer and don't compete.
//
// Cases per boundary node N:
//   I claim (own=1/2):
//     - peer_min exists & < myRank → peer wins
//     - else                         → I win
//   I don't claim (own=0):
//     - peer_min exists              → peer with lowest rank wins
//     - peer_min absent              → no claimant; mark INT_MAX (will be
//                                       filtered out by emit kernel — neither
//                                       "owner==myRank" nor "owner==peer" hits)
template<typename KeyType>
__global__ void resolveOwnershipKernel(
    const KeyType* d_localToGlobalSfc,    // [nodeCount], sorted ascending
    const KeyType* d_claimKeys,           // [numClaims], sorted by key (unused, kept for API stability)
    const int*     d_claimRanks,          // [numClaims], same order
    const int*     d_claimKeyOffsets,     // [numUniqueKeys+1] CSR
    const KeyType* d_uniqueClaimKeys,     // [numUniqueKeys] unique keys sorted
    int            numUniqueKeys,
    const uint8_t* d_isBoundary,          // [nodeCount] (unused, kept for API stability)
    const uint8_t* d_locallyClaimed,      // [nodeCount], 1 if I have a local element with this corner
    int*           d_authoritativeOwner,  // [nodeCount] output
    int            myRank,
    int            nodeCount)
{
    int n = blockIdx.x * blockDim.x + threadIdx.x;
    if (n >= nodeCount) return;
    (void)d_claimKeys; (void)d_isBoundary;  // suppress unused warnings

    bool iClaim = (d_locallyClaimed[n] != 0);

    // Look up our key in the unique claim list. Don't short-circuit on
    // d_isBoundary — corner-shared nodes (the bug we're chasing) are NOT
    // in cstone halo elements but ARE in peer claim lists; the lookup
    // must happen for them too.
    KeyType myKey = d_localToGlobalSfc[n];
    int lo = 0, hi = numUniqueKeys;
    while (lo < hi)
    {
        int mid = lo + (hi - lo) / 2;
        if (d_uniqueClaimKeys[mid] < myKey) lo = mid + 1;
        else                                 hi = mid;
    }
    bool found = (lo < numUniqueKeys && d_uniqueClaimKeys[lo] == myKey);

    // Ownership competes only among true claimants (rank label without
    // touch-flag set). Touchers don't compete for ownership but their
    // presence in the claim list signals they need the value.
    constexpr int TOUCH_FLAG_LOCAL = (int)0x80000000;
    int peerMin = INT_MAX;
    if (found)
    {
        int beg = d_claimKeyOffsets[lo];
        int end = d_claimKeyOffsets[lo + 1];
        for (int i = beg; i < end; ++i)
        {
            int rLabel = d_claimRanks[i];
            if ((rLabel & TOUCH_FLAG_LOCAL) != 0) continue;  // skip touchers
            if (rLabel < peerMin) peerMin = rLabel;
        }
    }

    if (iClaim)
    {
        d_authoritativeOwner[n] = (peerMin < myRank) ? peerMin : myRank;
    }
    else
    {
        // I don't claim it locally. If a true claimant peer exists, they own.
        d_authoritativeOwner[n] = (peerMin != INT_MAX) ? peerMin : -1;
    }
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void NodeHaloTopology<ElementTag, RealType, KeyType, AcceleratorTag>::buildFromCstoneHalos(
    const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain)
{
    constexpr int NPC = ElementTag::NodesPerElement;
    int rank     = domain.rank();
    int numRanks = domain.numRanks();

    if (numRanks == 1)
    {
        sendOffsets_.assign(1, 0);
        recvOffsets_.assign(1, 0);
        return;
    }

    // ----- Pull cstone halo bookkeeping (host-side) -----
    const auto& cstoneDom    = domain.getDomain();
    const auto& incomingHalo = cstoneDom.incomingHaloIndices();
    const auto& outgoingHalo = cstoneDom.outgoingHaloIndices();

    // Derive peer list: any rank with non-empty incoming or outgoing halo.
    std::vector<int> peers;
    peers.reserve(numRanks);
    for (int r = 0; r < numRanks; ++r)
    {
        if (r == rank) continue;
        bool inHas  = (size_t(r) < incomingHalo.size()) &&
                      (incomingHalo[r].count() > 0);
        bool outHas = (size_t(r) < outgoingHalo.size()) &&
                      (outgoingHalo[r].totalCount() > 0);
        if (inHas || outHas) peers.push_back(r);
    }

    // ----- Build flat range tables for incoming (recv) and outgoing (send) -----
    // incomingHalo[peer] is one contiguous range; one entry per peer.
    // outgoingHalo[peer] is a SendManifest with potentially multiple ranges.
    std::vector<int> h_inStart, h_inEnd, h_inPeer;
    std::vector<int> h_outStart, h_outEnd, h_outPeer;
    int totalInElems = 0, totalOutElems = 0;
    for (int p : peers)
    {
        if (size_t(p) < incomingHalo.size() && incomingHalo[p].count() > 0)
        {
            h_inStart.push_back(int(incomingHalo[p].start()));
            h_inEnd.push_back(int(incomingHalo[p].end()));
            h_inPeer.push_back(p);
            totalInElems += int(incomingHalo[p].count());
        }
        if (size_t(p) < outgoingHalo.size() && outgoingHalo[p].totalCount() > 0)
        {
            const auto& mfst = outgoingHalo[p];
            for (size_t r = 0; r < mfst.nRanges(); ++r)
            {
                h_outStart.push_back(int(mfst.rangeStart(r)));
                h_outEnd.push_back(int(mfst.rangeEnd(r)));
                h_outPeer.push_back(p);
                totalOutElems += int(mfst.count(r));
            }
        }
    }

    // ----- Pull connectivity (NPC device pointers) and ownership pointer -----
    const auto& d_conn_tuple = domain.getElementToNodeConnectivity();
    const KeyType* h_connPtrs[NPC];
    extractConnPtrs<KeyType>(d_conn_tuple, h_connPtrs, std::make_index_sequence<NPC>{});
    // Upload pointer array to device
    thrust::device_vector<const KeyType*> d_connPtrs(h_connPtrs, h_connPtrs + NPC);

    size_t nodeCount = domain.getNodeCount();
    const KeyType* d_sfcMap = thrust::raw_pointer_cast(domain.getLocalToGlobalSfcMap().data());

    // ====================================================================
    // TIEBREAKER STAGE A: identify boundary candidate nodes (corners of any
    // halo element) so we don't publish claims for interior nodes.
    // ====================================================================
    thrust::device_vector<uint8_t> d_isBoundary(nodeCount, 0u);
    {
        // Combine in+out element ranges into one pass so a node touched by
        // either side gets marked.
        std::vector<int> hAllStart  = h_inStart;  hAllStart.insert (hAllStart.end(),  h_outStart.begin(),  h_outStart.end());
        std::vector<int> hAllEnd    = h_inEnd;    hAllEnd.insert   (hAllEnd.end(),    h_outEnd.begin(),    h_outEnd.end());
        int              totalAll   = totalInElems + totalOutElems;
        if (totalAll > 0)
        {
            thrust::device_vector<int> d_aStart(hAllStart.begin(), hAllStart.end());
            thrust::device_vector<int> d_aEnd  (hAllEnd.begin(),   hAllEnd.end());
            const int blockSize = 256;
            int blocks = (totalAll + blockSize - 1) / blockSize;
            markBoundaryCandidatesKernel<KeyType, NPC><<<blocks, blockSize>>>(
                thrust::raw_pointer_cast(d_connPtrs.data()),
                thrust::raw_pointer_cast(d_aStart.data()),
                thrust::raw_pointer_cast(d_aEnd.data()),
                int(hAllStart.size()),
                thrust::raw_pointer_cast(d_isBoundary.data()),
                totalAll);
            cudaDeviceSynchronize();
        }
    }

    // ====================================================================
    // TIEBREAKER STAGE B: gather (sfc_key) of locally-CLAIMED boundary
    // nodes. The "claim" predicate is "I have a local element [startIdx,
    // endIdx) whose corner is N" — NOT d_nodeOwnership, which has Step 2's
    // direction-biased yield baked in. Using d_nodeOwnership here is wrong
    // because Step 2 yields too aggressively: any node touched by a halo
    // element from a lower-SFC peer is yielded even if the global lowest-
    // rank claimant is actually us. We re-derive the unbiased claim flag
    // here directly from local element corners.
    // ====================================================================
    thrust::device_vector<uint8_t> d_locallyClaimed(nodeCount, 0u);
    {
        size_t startIdx = domain.startIndex();
        size_t endIdx   = domain.endIndex();
        size_t numLocal = endIdx - startIdx;
        if (numLocal > 0)
        {
            const int blockSize = 256;
            int blocks = (numLocal + blockSize - 1) / blockSize;
            markLocallyClaimedKernel<KeyType, NPC><<<blocks, blockSize>>>(
                thrust::raw_pointer_cast(d_connPtrs.data()),
                startIdx, endIdx,
                thrust::raw_pointer_cast(d_locallyClaimed.data()));
            cudaDeviceSynchronize();
        }
    }
    const uint8_t* d_claimedPtr  = thrust::raw_pointer_cast(d_locallyClaimed.data());

    // Build my publish list on device: keys of all nodes I "touch" — i.e.
    // any node N that's a corner of any element I see, whether the element
    // is local [startIdx, endIdx) OR a halo element on this rank. This is
    // a TOUCH list, not an OWNERSHIP CLAIM list:
    //
    //   - For ownership resolution, what matters is whether I have a local
    //     element with N (we use d_locallyClaimed for that, downstream).
    //   - For the "who needs N as a ghost" question (send-side derivation),
    //     we need to know which peers touch N, even if they don't own it.
    //     A peer that has N only in a halo element from me still needs me
    //     to send N's value at every CG iteration.
    //
    // Publishing TOUCH (= d_isBoundary, includes halo-element corners on
    // this rank) lets every rank reconstruct the full set of ranks that
    // need each node, without a request-back protocol.
    //
    // We tag each published rank with a flag bit:
    //   bit 31 set       → "I touch this key but don't own it (halo only)"
    //   bit 31 not set   → "I own this key (have a local element with this corner)"
    //
    // Ownership resolution looks only at unflagged entries (true owners
    // compete for min-rank).
    // Send-side derivation looks at ALL entries (touchers and other owners
    // both need values).
    static constexpr int TOUCH_FLAG = (int)0x80000000;  // high bit; -ve when interpreted signed

    thrust::device_vector<KeyType> d_myClaimKeys(nodeCount);
    thrust::device_vector<int>     d_myClaimRanks(nodeCount);
    int myClaimCount = 0;
    {
        thrust::device_vector<KeyType> d_keysAll(nodeCount);
        thrust::device_vector<int>     d_ranksAll(nodeCount);
        thrust::device_vector<uint8_t> d_flagAll(nodeCount);
        const uint8_t* dIsB = thrust::raw_pointer_cast(d_isBoundary.data());
        thrust::for_each(
            thrust::counting_iterator<int>(0),
            thrust::counting_iterator<int>(int(nodeCount)),
            [d_keysAll_ptr = thrust::raw_pointer_cast(d_keysAll.data()),
             d_ranksAll_ptr = thrust::raw_pointer_cast(d_ranksAll.data()),
             d_flagAll_ptr = thrust::raw_pointer_cast(d_flagAll.data()),
             dIsB, d_claimedPtr, d_sfcMap, rank]
            __device__ (int n) {
                d_keysAll_ptr[n] = d_sfcMap[n];
                bool touch = (dIsB[n] != 0u);
                bool claim = (d_claimedPtr[n] != 0u);
                bool publish = (touch || claim);
                int  rankLabel = claim ? rank : (rank | TOUCH_FLAG);
                d_ranksAll_ptr[n] = rankLabel;
                d_flagAll_ptr[n]  = publish ? 1u : 0u;
            });

        // Compact keys + ranks together using zip iterators.
        auto begIn  = thrust::make_zip_iterator(thrust::make_tuple(d_keysAll.begin(),  d_ranksAll.begin()));
        auto endIn  = thrust::make_zip_iterator(thrust::make_tuple(d_keysAll.end(),    d_ranksAll.end()));
        auto begOut = thrust::make_zip_iterator(thrust::make_tuple(d_myClaimKeys.begin(), d_myClaimRanks.begin()));
        auto endOut = thrust::copy_if(begIn, endIn, d_flagAll.begin(), begOut,
            [] __device__ (uint8_t f) { return f != 0; });
        myClaimCount = endOut - begOut;
        d_myClaimKeys.resize(myClaimCount);
        d_myClaimRanks.resize(myClaimCount);

        // Sort by key (zip).
        if (myClaimCount > 1)
        {
            thrust::sort_by_key(d_myClaimKeys.begin(), d_myClaimKeys.end(),
                                d_myClaimRanks.begin());
        }
    }

    // Exchange claim counts with each peer (point-to-point Isend/Irecv on
    // host counters; tiny — int per peer).
    int numPeers = int(peers.size());
    std::vector<int> peerSendCount(numPeers, myClaimCount);  // we send same list to each peer
    std::vector<int> peerRecvCount(numPeers, 0);
    {
        std::vector<MPI_Request> reqs;
        reqs.reserve(2 * numPeers);
        constexpr int tagCount = 0x4d54;  // "MT" — disambiguate
        for (int i = 0; i < numPeers; ++i)
        {
            MPI_Request r;
            MPI_Irecv(&peerRecvCount[i], 1, MPI_INT, peers[i], tagCount, MPI_COMM_WORLD, &r);
            reqs.push_back(r);
        }
        for (int i = 0; i < numPeers; ++i)
        {
            MPI_Request r;
            MPI_Isend(&peerSendCount[i], 1, MPI_INT, peers[i], tagCount, MPI_COMM_WORLD, &r);
            reqs.push_back(r);
        }
        MPI_Waitall(int(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
    }

    // Exchange the actual claim keys AND flagged rank labels via CUDA-aware MPI.
    int totalRecvKeys = 0;
    std::vector<int> peerRecvDispl(numPeers + 1, 0);
    for (int i = 0; i < numPeers; ++i)
    {
        peerRecvDispl[i + 1] = peerRecvDispl[i] + peerRecvCount[i];
    }
    totalRecvKeys = peerRecvDispl.back();

    thrust::device_vector<KeyType> d_recvClaimKeys(totalRecvKeys);
    thrust::device_vector<int>     d_recvClaimRanks(totalRecvKeys);
    auto mpiKeyType = (sizeof(KeyType) == 8) ? MPI_UINT64_T : MPI_UNSIGNED;
    {
        std::vector<MPI_Request> reqs;
        reqs.reserve(4 * numPeers);
        constexpr int tagKeys  = 0x4d55;  // "MU"
        constexpr int tagRanks = 0x4d56;  // "MV" — flagged rank labels
        for (int i = 0; i < numPeers; ++i)
        {
            int rcnt = peerRecvCount[i];
            if (rcnt == 0) continue;
            MPI_Request rk, rr;
            MPI_Irecv(thrust::raw_pointer_cast(d_recvClaimKeys.data()) + peerRecvDispl[i],
                      rcnt, mpiKeyType, peers[i], tagKeys, MPI_COMM_WORLD, &rk);
            MPI_Irecv(thrust::raw_pointer_cast(d_recvClaimRanks.data()) + peerRecvDispl[i],
                      rcnt, MPI_INT, peers[i], tagRanks, MPI_COMM_WORLD, &rr);
            reqs.push_back(rk);
            reqs.push_back(rr);
        }
        for (int i = 0; i < numPeers; ++i)
        {
            int scnt = myClaimCount;
            if (scnt == 0) continue;
            MPI_Request sk, sr;
            MPI_Isend(thrust::raw_pointer_cast(d_myClaimKeys.data()),
                      scnt, mpiKeyType, peers[i], tagKeys, MPI_COMM_WORLD, &sk);
            MPI_Isend(thrust::raw_pointer_cast(d_myClaimRanks.data()),
                      scnt, MPI_INT, peers[i], tagRanks, MPI_COMM_WORLD, &sr);
            reqs.push_back(sk);
            reqs.push_back(sr);
        }
        MPI_Waitall(int(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
    }

    // Sort received claims by (key, rank) so per-key min-rank reduces correctly.
    if (totalRecvKeys > 0)
    {
        thrust::sort_by_key(d_recvClaimKeys.begin(), d_recvClaimKeys.end(),
                            d_recvClaimRanks.begin());
    }

    // Build CSR of (unique_key → claim range) so resolveOwnershipKernel can
    // look up by binary search.
    thrust::device_vector<KeyType> d_uniqueKeys;
    thrust::device_vector<int>     d_uniqueOffsets;
    int numUnique = 0;
    if (totalRecvKeys > 0)
    {
        d_uniqueKeys.resize(totalRecvKeys);
        thrust::copy(d_recvClaimKeys.begin(), d_recvClaimKeys.end(), d_uniqueKeys.begin());
        auto uend = thrust::unique(d_uniqueKeys.begin(), d_uniqueKeys.end());
        numUnique = int(uend - d_uniqueKeys.begin());
        d_uniqueKeys.resize(numUnique);
        // For each unique key, lower_bound gives its first index in d_recvClaimKeys.
        // We need numUnique+1 offsets; last is totalRecvKeys.
        d_uniqueOffsets.resize(numUnique + 1);
        thrust::lower_bound(d_recvClaimKeys.begin(), d_recvClaimKeys.end(),
                            d_uniqueKeys.begin(), d_uniqueKeys.end(),
                            d_uniqueOffsets.begin());
        // Set the sentinel offset[numUnique] = totalRecvKeys.
        int lastOff = totalRecvKeys;
        cudaMemcpy(thrust::raw_pointer_cast(d_uniqueOffsets.data()) + numUnique,
                   &lastOff, sizeof(int), cudaMemcpyHostToDevice);
    }

    // ====================================================================
    // TIEBREAKER STAGE C: per-node compute authoritative owner.
    //   non-boundary owned → myRank
    //   boundary node not in unique-key claims → myRank (no peer claimed)
    //   boundary node in claims → min(myRank, all claimant ranks) if I own,
    //                              else min(claimant ranks)
    // ====================================================================
    thrust::device_vector<int> d_authoritativeOwner(nodeCount, -1);
    {
        const int blockSize = 256;
        int blocks = (int(nodeCount) + blockSize - 1) / blockSize;
        resolveOwnershipKernel<KeyType><<<blocks, blockSize>>>(
            d_sfcMap,
            thrust::raw_pointer_cast(d_recvClaimKeys.data()),
            thrust::raw_pointer_cast(d_recvClaimRanks.data()),
            thrust::raw_pointer_cast(d_uniqueOffsets.data()),
            thrust::raw_pointer_cast(d_uniqueKeys.data()),
            numUnique,
            thrust::raw_pointer_cast(d_isBoundary.data()),
            d_claimedPtr,
            thrust::raw_pointer_cast(d_authoritativeOwner.data()),
            rank,
            int(nodeCount));
        cudaDeviceSynchronize();
    }
    const int* d_authPtr = thrust::raw_pointer_cast(d_authoritativeOwner.data());

    // Mutate d_nodeOwnership_ to match the post-tiebreaker state, just as
    // host path does at L1212-1213. Required for downstream FEM code (DOF
    // handler, RHS assembly) to see the correct "owned vs ghost" flag.
    // Without this, ranks that lost the tiebreaker still report ownership=1
    // for nodes they no longer own, causing duplicate RHS contributions
    // and the matrix/rhs norm drift between host and v2 paths.
    {
        auto& dOwnMutable = const_cast<typename HaloData<ElementTag, RealType, KeyType, AcceleratorTag>::
                                       template DeviceVector<uint8_t>&>(domain.getNodeOwnershipMap());
        thrust::for_each(
            thrust::counting_iterator<int>(0),
            thrust::counting_iterator<int>(int(nodeCount)),
            [d_own_ptr = thrust::raw_pointer_cast(dOwnMutable.data()), d_authPtr, rank]
            __device__ (int n) {
                int owner = d_authPtr[n];
                // -1 = no owner found (interior unowned, shouldn't happen in
                // practice). Otherwise set to 1 iff this rank is the owner.
                d_own_ptr[n] = (owner == rank) ? uint8_t(1) : uint8_t(0);
            });
        cudaDeviceSynchronize();
    }

    // ----- Helper: emit per-peer (peer, node) tuples then compact -----
    // For send: emit one tuple per (owned node, peer claimant) using
    // emitSendNodesKernel. Output stride bounds claimants/node — use the
    // largest unique-claim range as upper bound.
    // For recv: emit at most one tuple per node (owner != me).
    auto compactSortBin = [&](
        thrust::device_vector<int>& d_outPeer,
        thrust::device_vector<int>& d_outNode,
        std::vector<std::vector<int>>& outIdsPerPeer)
    {
        outIdsPerPeer.assign(peers.size(), {});

        // Stream-compact: drop entries with peer < 0
        thrust::device_vector<int> d_keepPeer(d_outPeer.size());
        thrust::device_vector<int> d_keepNode(d_outNode.size());
        auto end = thrust::copy_if(
            thrust::make_zip_iterator(thrust::make_tuple(d_outPeer.begin(), d_outNode.begin())),
            thrust::make_zip_iterator(thrust::make_tuple(d_outPeer.end(),   d_outNode.end())),
            thrust::make_zip_iterator(thrust::make_tuple(d_keepPeer.begin(), d_keepNode.begin())),
            [] __device__ (const thrust::tuple<int,int>& t) {
                return thrust::get<0>(t) >= 0;
            });
        size_t kept = end - thrust::make_zip_iterator(thrust::make_tuple(d_keepPeer.begin(), d_keepNode.begin()));
        d_keepPeer.resize(kept);
        d_keepNode.resize(kept);

        // Sort by (peer, node) so per-peer dedup with one unique pass works.
        thrust::sort(
            thrust::make_zip_iterator(thrust::make_tuple(d_keepPeer.begin(), d_keepNode.begin())),
            thrust::make_zip_iterator(thrust::make_tuple(d_keepPeer.end(),   d_keepNode.end())));
        auto uend = thrust::unique(
            thrust::make_zip_iterator(thrust::make_tuple(d_keepPeer.begin(), d_keepNode.begin())),
            thrust::make_zip_iterator(thrust::make_tuple(d_keepPeer.end(),   d_keepNode.end())));
        size_t uniq = uend - thrust::make_zip_iterator(thrust::make_tuple(d_keepPeer.begin(), d_keepNode.begin()));

        std::vector<int> hPeerOut(uniq), hNodeOut(uniq);
        if (uniq > 0)
        {
            cudaMemcpy(hPeerOut.data(), thrust::raw_pointer_cast(d_keepPeer.data()),
                       uniq * sizeof(int), cudaMemcpyDeviceToHost);
            cudaMemcpy(hNodeOut.data(), thrust::raw_pointer_cast(d_keepNode.data()),
                       uniq * sizeof(int), cudaMemcpyDeviceToHost);
        }

        std::unordered_map<int, size_t> peerIdx;
        for (size_t i = 0; i < peers.size(); ++i) peerIdx[peers[i]] = i;
        for (size_t i = 0; i < uniq; ++i)
        {
            auto it = peerIdx.find(hPeerOut[i]);
            if (it != peerIdx.end()) outIdsPerPeer[it->second].push_back(hNodeOut[i]);
        }
    };

    // ===== SEND side =====
    // Emit owned nodes to every claimant peer (excluding self). Bound output
    // stride at numPeers (max possible claimants per node). At cube
    // partitions worst case is 7 (8-way corner with 8 ranks: 7 peer claimants).
    std::vector<std::vector<int>> sendIdsPerPeer;
    {
        int outStride = std::max(1, std::min(numPeers, 32));
        thrust::device_vector<int> d_outPeer(size_t(nodeCount) * outStride, -1);
        thrust::device_vector<int> d_outNode(size_t(nodeCount) * outStride, -1);
        const int blockSize = 256;
        int blocks = (int(nodeCount) + blockSize - 1) / blockSize;
        emitSendNodesKernel<KeyType><<<blocks, blockSize>>>(
            d_sfcMap,
            d_authPtr,
            thrust::raw_pointer_cast(d_recvClaimRanks.data()),
            thrust::raw_pointer_cast(d_uniqueOffsets.data()),
            thrust::raw_pointer_cast(d_uniqueKeys.data()),
            numUnique,
            rank,
            outStride,
            thrust::raw_pointer_cast(d_outPeer.data()),
            thrust::raw_pointer_cast(d_outNode.data()),
            int(nodeCount));
        cudaDeviceSynchronize();
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess)
        {
            std::cerr << "[Rank " << rank << "] emitSendNodesKernel failed: "
                      << cudaGetErrorString(err) << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        compactSortBin(d_outPeer, d_outNode, sendIdsPerPeer);
    }

    // ===== RECV side =====
    std::vector<std::vector<int>> recvIdsPerPeer;
    {
        thrust::device_vector<int> d_outPeer(size_t(nodeCount), -1);
        thrust::device_vector<int> d_outNode(size_t(nodeCount), -1);
        const int blockSize = 256;
        int blocks = (int(nodeCount) + blockSize - 1) / blockSize;
        emitRecvNodesKernel<<<blocks, blockSize>>>(
            d_authPtr, rank,
            thrust::raw_pointer_cast(d_outPeer.data()),
            thrust::raw_pointer_cast(d_outNode.data()),
            int(nodeCount));
        cudaDeviceSynchronize();
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess)
        {
            std::cerr << "[Rank " << rank << "] emitRecvNodesKernel failed: "
                      << cudaGetErrorString(err) << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        compactSortBin(d_outPeer, d_outNode, recvIdsPerPeer);
    }

    // ----- Compact into peer-list / CSR structure -----
    peers_.clear();
    sendOffsets_.assign(1, 0);
    recvOffsets_.assign(1, 0);
    std::vector<int> sendNodesHost, recvNodesHost;
    for (size_t i = 0; i < peers.size(); ++i)
    {
        int sCnt = int(sendIdsPerPeer[i].size());
        int rCnt = int(recvIdsPerPeer[i].size());
        if (sCnt == 0 && rCnt == 0) continue;
        peers_.push_back(peers[i]);
        for (int v : sendIdsPerPeer[i]) sendNodesHost.push_back(v);
        for (int v : recvIdsPerPeer[i]) recvNodesHost.push_back(v);
        sendOffsets_.push_back(int(sendNodesHost.size()));
        recvOffsets_.push_back(int(recvNodesHost.size()));
    }

    // Upload node-id arrays to device
    sendNodeIds_.resize(sendNodesHost.size());
    recvNodeIds_.resize(recvNodesHost.size());
    if (!sendNodesHost.empty())
        cudaMemcpy(thrust::raw_pointer_cast(sendNodeIds_.data()), sendNodesHost.data(),
                   sendNodesHost.size() * sizeof(int), cudaMemcpyHostToDevice);
    if (!recvNodesHost.empty())
        cudaMemcpy(thrust::raw_pointer_cast(recvNodeIds_.data()), recvNodesHost.data(),
                   recvNodesHost.size() * sizeof(int), cudaMemcpyHostToDevice);

    std::cout << "Rank " << rank << ": NodeHaloTopo[v2] " << peers_.size() << " peers, "
              << sendNodesHost.size() << " send nodes, "
              << recvNodesHost.size() << " recv nodes" << std::endl;
    std::cout.flush();
}

// ============================================================================
// validateAgainstHost: build both paths into separate topos, diff per-peer
// sorted node-id lists. Returns true on bit-exact match.
//
// Implementation note: we cannot reuse the constructor because it would
// recurse via the env flag. Build host path manually, then v2 path
// manually, then compare on host.
//
// Subtle: the host path mutates d_nodeOwnership_ during its tiebreaker
// yield step (host-path L1102-1129). We snapshot ownership before the
// host build and restore before the v2 build. This means v2 runs
// against the *original* (pre-yield) ownership. If host's tiebreaker
// changed any ownership flags, v2's send/recv lists will differ — that's
// the expected signal that v2 needs its own tiebreaker pass before
// becoming default.
// ============================================================================
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
bool NodeHaloTopology<ElementTag, RealType, KeyType, AcceleratorTag>::validateAgainstHost(
    const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain)
{
    int rank     = domain.rank();
    int numRanks = domain.numRanks();
    if (numRanks == 1) return true;

    // Snapshot d_nodeOwnership_ before running host path, since the host
    // path mutates it (yields duplicate ownership to lowest rank). v2 must
    // see the same starting ownership state to make the comparison fair.
    auto& dOwn = const_cast<typename HaloData<ElementTag, RealType, KeyType, AcceleratorTag>::
                            template DeviceVector<uint8_t>&>(domain.getNodeOwnershipMap());
    std::vector<uint8_t> ownSnapshot(dOwn.size());
    if (!ownSnapshot.empty())
        cudaMemcpy(ownSnapshot.data(), thrust::raw_pointer_cast(dOwn.data()),
                   ownSnapshot.size() * sizeof(uint8_t), cudaMemcpyDeviceToHost);

    NodeHaloTopology<ElementTag, RealType, KeyType, AcceleratorTag> hostTopo{};
    NodeHaloTopology<ElementTag, RealType, KeyType, AcceleratorTag> v2Topo{};
    buildNodeHaloTopologyHostPath(hostTopo, domain);

    // Restore ownership before v2 build so it sees pre-yield state.
    if (!ownSnapshot.empty())
        cudaMemcpy(thrust::raw_pointer_cast(dOwn.data()), ownSnapshot.data(),
                   ownSnapshot.size() * sizeof(uint8_t), cudaMemcpyHostToDevice);

    v2Topo.buildFromCstoneHalos(domain);

    auto pullToHost = [](const auto& d_vec) {
        std::vector<int> h(d_vec.size());
        if (!h.empty())
            cudaMemcpy(h.data(), thrust::raw_pointer_cast(d_vec.data()),
                       h.size() * sizeof(int), cudaMemcpyDeviceToHost);
        return h;
    };

    auto buildPerPeerSorted = [&](const auto& topo) {
        std::vector<std::vector<int>> sendByPeer, recvByPeer;
        std::vector<int> hSend = pullToHost(topo.sendNodeIds_);
        std::vector<int> hRecv = pullToHost(topo.recvNodeIds_);
        sendByPeer.resize(topo.peers_.size());
        recvByPeer.resize(topo.peers_.size());
        for (size_t i = 0; i < topo.peers_.size(); ++i)
        {
            int sb = topo.sendOffsets_[i],  se = topo.sendOffsets_[i+1];
            int rb = topo.recvOffsets_[i],  re = topo.recvOffsets_[i+1];
            sendByPeer[i].assign(hSend.begin() + sb, hSend.begin() + se);
            recvByPeer[i].assign(hRecv.begin() + rb, hRecv.begin() + re);
            std::sort(sendByPeer[i].begin(), sendByPeer[i].end());
            std::sort(recvByPeer[i].begin(), recvByPeer[i].end());
        }
        return std::make_tuple(std::move(sendByPeer), std::move(recvByPeer));
    };

    auto [hSend, hRecv] = buildPerPeerSorted(hostTopo);
    auto [vSend, vRecv] = buildPerPeerSorted(v2Topo);

    // Build peer-rank -> sorted lists maps to compare order-independently.
    auto byPeerMap = [](const auto& topo, const auto& byPeerVec) {
        std::map<int, std::vector<int>> m;
        for (size_t i = 0; i < topo.peers_.size(); ++i) m[topo.peers_[i]] = byPeerVec[i];
        return m;
    };
    auto hSendMap = byPeerMap(hostTopo, hSend);
    auto hRecvMap = byPeerMap(hostTopo, hRecv);
    auto vSendMap = byPeerMap(v2Topo,   vSend);
    auto vRecvMap = byPeerMap(v2Topo,   vRecv);

    bool ok = true;
    auto diffMap = [&](const char* label, const auto& a, const auto& b) {
        if (a.size() != b.size())
        {
            std::cerr << "[Rank " << rank << "] " << label
                      << " peer-count mismatch: host=" << a.size()
                      << " v2=" << b.size() << std::endl;
            ok = false;
        }
        for (const auto& [peer, hList] : a)
        {
            auto it = b.find(peer);
            if (it == b.end())
            {
                std::cerr << "[Rank " << rank << "] " << label
                          << " peer " << peer << " present in host but not v2" << std::endl;
                ok = false;
                continue;
            }
            const auto& vList = it->second;
            if (hList.size() != vList.size())
            {
                std::cerr << "[Rank " << rank << "] " << label
                          << " peer " << peer << " size: host=" << hList.size()
                          << " v2=" << vList.size() << std::endl;
                ok = false;
                continue;
            }
            for (size_t i = 0; i < hList.size(); ++i)
            {
                if (hList[i] != vList[i])
                {
                    std::cerr << "[Rank " << rank << "] " << label
                              << " peer " << peer << " idx " << i
                              << " host=" << hList[i] << " v2=" << vList[i] << std::endl;
                    ok = false;
                    break;
                }
            }
        }
        // Check for v2-only peers
        for (const auto& [peer, vList] : b)
        {
            if (a.find(peer) == a.end())
            {
                std::cerr << "[Rank " << rank << "] " << label
                          << " peer " << peer << " present in v2 but not host" << std::endl;
                ok = false;
            }
        }
    };
    diffMap("send", hSendMap, vSendMap);
    diffMap("recv", hRecvMap, vRecvMap);
    if (ok)
        std::cout << "[Rank " << rank << "] NodeHaloTopo v2 == host (Gate 1 pass)" << std::endl;
    std::cout.flush();
    std::cerr.flush();
    return ok;
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

// Template instantiations for extractElementNodeCoordsKernel (HexTag only)
template __global__ void extractElementNodeCoordsKernel<HexTag, unsigned, float>(
    const float* x, const float* y, const float* z,
    const unsigned* idx0, const unsigned* idx1, const unsigned* idx2, const unsigned* idx3,
    const unsigned* idx4, const unsigned* idx5, const unsigned* idx6, const unsigned* idx7,
    float* elemOrigX0, float* elemOrigY0, float* elemOrigZ0,
    float* elemOrigX1, float* elemOrigY1, float* elemOrigZ1,
    float* elemOrigX2, float* elemOrigY2, float* elemOrigZ2,
    float* elemOrigX3, float* elemOrigY3, float* elemOrigZ3,
    float* elemOrigX4, float* elemOrigY4, float* elemOrigZ4,
    float* elemOrigX5, float* elemOrigY5, float* elemOrigZ5,
    float* elemOrigX6, float* elemOrigY6, float* elemOrigZ6,
    float* elemOrigX7, float* elemOrigY7, float* elemOrigZ7,
    int numElements);

template __global__ void extractElementNodeCoordsKernel<HexTag, unsigned, double>(
    const double* x, const double* y, const double* z,
    const unsigned* idx0, const unsigned* idx1, const unsigned* idx2, const unsigned* idx3,
    const unsigned* idx4, const unsigned* idx5, const unsigned* idx6, const unsigned* idx7,
    double* elemOrigX0, double* elemOrigY0, double* elemOrigZ0,
    double* elemOrigX1, double* elemOrigY1, double* elemOrigZ1,
    double* elemOrigX2, double* elemOrigY2, double* elemOrigZ2,
    double* elemOrigX3, double* elemOrigY3, double* elemOrigZ3,
    double* elemOrigX4, double* elemOrigY4, double* elemOrigZ4,
    double* elemOrigX5, double* elemOrigY5, double* elemOrigZ5,
    double* elemOrigX6, double* elemOrigY6, double* elemOrigZ6,
    double* elemOrigX7, double* elemOrigY7, double* elemOrigZ7,
    int numElements);

template __global__ void extractElementNodeCoordsKernel<HexTag, uint64_t, float>(
    const float* x, const float* y, const float* z,
    const uint64_t* idx0, const uint64_t* idx1, const uint64_t* idx2, const uint64_t* idx3,
    const uint64_t* idx4, const uint64_t* idx5, const uint64_t* idx6, const uint64_t* idx7,
    float* elemOrigX0, float* elemOrigY0, float* elemOrigZ0,
    float* elemOrigX1, float* elemOrigY1, float* elemOrigZ1,
    float* elemOrigX2, float* elemOrigY2, float* elemOrigZ2,
    float* elemOrigX3, float* elemOrigY3, float* elemOrigZ3,
    float* elemOrigX4, float* elemOrigY4, float* elemOrigZ4,
    float* elemOrigX5, float* elemOrigY5, float* elemOrigZ5,
    float* elemOrigX6, float* elemOrigY6, float* elemOrigZ6,
    float* elemOrigX7, float* elemOrigY7, float* elemOrigZ7,
    int numElements);

template __global__ void extractElementNodeCoordsKernel<HexTag, uint64_t, double>(
    const double* x, const double* y, const double* z,
    const uint64_t* idx0, const uint64_t* idx1, const uint64_t* idx2, const uint64_t* idx3,
    const uint64_t* idx4, const uint64_t* idx5, const uint64_t* idx6, const uint64_t* idx7,
    double* elemOrigX0, double* elemOrigY0, double* elemOrigZ0,
    double* elemOrigX1, double* elemOrigY1, double* elemOrigZ1,
    double* elemOrigX2, double* elemOrigY2, double* elemOrigZ2,
    double* elemOrigX3, double* elemOrigY3, double* elemOrigZ3,
    double* elemOrigX4, double* elemOrigY4, double* elemOrigZ4,
    double* elemOrigX5, double* elemOrigY5, double* elemOrigZ5,
    double* elemOrigX6, double* elemOrigY6, double* elemOrigZ6,
    double* elemOrigX7, double* elemOrigY7, double* elemOrigZ7,
    int numElements);

// Template instantiations for rebuildNodeCoordsFromElementsKernel (HexTag only)
template __global__ void rebuildNodeCoordsFromElementsKernel<HexTag, unsigned, float>(
    const unsigned* connKey0, const unsigned* connKey1, const unsigned* connKey2, const unsigned* connKey3,
    const unsigned* connKey4, const unsigned* connKey5, const unsigned* connKey6, const unsigned* connKey7,
    const float* elemOrigX0, const float* elemOrigY0, const float* elemOrigZ0,
    const float* elemOrigX1, const float* elemOrigY1, const float* elemOrigZ1,
    const float* elemOrigX2, const float* elemOrigY2, const float* elemOrigZ2,
    const float* elemOrigX3, const float* elemOrigY3, const float* elemOrigZ3,
    const float* elemOrigX4, const float* elemOrigY4, const float* elemOrigZ4,
    const float* elemOrigX5, const float* elemOrigY5, const float* elemOrigZ5,
    const float* elemOrigX6, const float* elemOrigY6, const float* elemOrigZ6,
    const float* elemOrigX7, const float* elemOrigY7, const float* elemOrigZ7,
    float* nodeX, float* nodeY, float* nodeZ,
    int numElements);

template __global__ void rebuildNodeCoordsFromElementsKernel<HexTag, unsigned, double>(
    const unsigned* connKey0, const unsigned* connKey1, const unsigned* connKey2, const unsigned* connKey3,
    const unsigned* connKey4, const unsigned* connKey5, const unsigned* connKey6, const unsigned* connKey7,
    const double* elemOrigX0, const double* elemOrigY0, const double* elemOrigZ0,
    const double* elemOrigX1, const double* elemOrigY1, const double* elemOrigZ1,
    const double* elemOrigX2, const double* elemOrigY2, const double* elemOrigZ2,
    const double* elemOrigX3, const double* elemOrigY3, const double* elemOrigZ3,
    const double* elemOrigX4, const double* elemOrigY4, const double* elemOrigZ4,
    const double* elemOrigX5, const double* elemOrigY5, const double* elemOrigZ5,
    const double* elemOrigX6, const double* elemOrigY6, const double* elemOrigZ6,
    const double* elemOrigX7, const double* elemOrigY7, const double* elemOrigZ7,
    double* nodeX, double* nodeY, double* nodeZ,
    int numElements);

template __global__ void rebuildNodeCoordsFromElementsKernel<HexTag, uint64_t, float>(
    const uint64_t* connKey0, const uint64_t* connKey1, const uint64_t* connKey2, const uint64_t* connKey3,
    const uint64_t* connKey4, const uint64_t* connKey5, const uint64_t* connKey6, const uint64_t* connKey7,
    const float* elemOrigX0, const float* elemOrigY0, const float* elemOrigZ0,
    const float* elemOrigX1, const float* elemOrigY1, const float* elemOrigZ1,
    const float* elemOrigX2, const float* elemOrigY2, const float* elemOrigZ2,
    const float* elemOrigX3, const float* elemOrigY3, const float* elemOrigZ3,
    const float* elemOrigX4, const float* elemOrigY4, const float* elemOrigZ4,
    const float* elemOrigX5, const float* elemOrigY5, const float* elemOrigZ5,
    const float* elemOrigX6, const float* elemOrigY6, const float* elemOrigZ6,
    const float* elemOrigX7, const float* elemOrigY7, const float* elemOrigZ7,
    float* nodeX, float* nodeY, float* nodeZ,
    int numElements);

template __global__ void rebuildNodeCoordsFromElementsKernel<HexTag, uint64_t, double>(
    const uint64_t* connKey0, const uint64_t* connKey1, const uint64_t* connKey2, const uint64_t* connKey3,
    const uint64_t* connKey4, const uint64_t* connKey5, const uint64_t* connKey6, const uint64_t* connKey7,
    const double* elemOrigX0, const double* elemOrigY0, const double* elemOrigZ0,
    const double* elemOrigX1, const double* elemOrigY1, const double* elemOrigZ1,
    const double* elemOrigX2, const double* elemOrigY2, const double* elemOrigZ2,
    const double* elemOrigX3, const double* elemOrigY3, const double* elemOrigZ3,
    const double* elemOrigX4, const double* elemOrigY4, const double* elemOrigZ4,
    const double* elemOrigX5, const double* elemOrigY5, const double* elemOrigZ5,
    const double* elemOrigX6, const double* elemOrigY6, const double* elemOrigZ6,
    const double* elemOrigX7, const double* elemOrigY7, const double* elemOrigZ7,
    double* nodeX, double* nodeY, double* nodeZ,
    int numElements);

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

// Explicit instantiations for NodeHaloTopology
template struct NodeHaloTopology<TetTag, float, unsigned, cstone::GpuTag>;
template struct NodeHaloTopology<TetTag, double, unsigned, cstone::GpuTag>;
template struct NodeHaloTopology<TetTag, float, uint64_t, cstone::GpuTag>;
template struct NodeHaloTopology<TetTag, double, uint64_t, cstone::GpuTag>;
template struct NodeHaloTopology<HexTag, float, unsigned, cstone::GpuTag>;
template struct NodeHaloTopology<HexTag, double, unsigned, cstone::GpuTag>;
template struct NodeHaloTopology<HexTag, float, uint64_t, cstone::GpuTag>;
template struct NodeHaloTopology<HexTag, double, uint64_t, cstone::GpuTag>;

} // namespace mars
