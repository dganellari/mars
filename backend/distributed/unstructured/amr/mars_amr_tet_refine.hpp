#pragma once

#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/transform.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>
#include <thrust/count.h>

#include "cstone/cuda/cuda_utils.hpp"
#include "mars_amr_hex_refine.hpp" // reuse EdgeKey

namespace mars
{
namespace amr
{

// Tet4 red refinement: 1 tet -> 8 children via 6 edge midpoints.
//
//      3
//     /|\
//    / | \
//   /  |  \
//  0---+---2
//   \  |  /
//    \ | /
//     \|/
//      1
//
// 6 edges: 0-1, 0-2, 0-3, 1-2, 1-3, 2-3
// 6 midpoints: m01, m02, m03, m12, m13, m23
// 8 child tets (Bey's red refinement):
//   4 corner tets (each keeps one original vertex + 3 adjacent midpoints)
//   4 interior tets (octahedron triangulation, choice affects quality)

__global__ void countChildTetsKernel(const uint8_t* marks, uint64_t* childCounts, size_t numElements)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements) return;
    childCounts[e] = (marks[e] > 0) ? 8 : 1;
}

// Emit 6 edge keys per marked tet
template<typename KeyType>
__global__ void emitTetEdgeKeysKernel(const KeyType* conn0,
                                      const KeyType* conn1,
                                      const KeyType* conn2,
                                      const KeyType* conn3,
                                      const uint8_t* marks,
                                      const uint64_t* markedPrefix,
                                      uint64_t* edgeKeys,
                                      size_t numElements)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements || marks[e] <= 0) return;

    uint64_t base = markedPrefix[e] * 6;
    KeyType n[4]  = {conn0[e], conn1[e], conn2[e], conn3[e]};

    // 6 edges of tet4
    edgeKeys[base + 0] = EdgeKey::encode(n[0], n[1]);
    edgeKeys[base + 1] = EdgeKey::encode(n[0], n[2]);
    edgeKeys[base + 2] = EdgeKey::encode(n[0], n[3]);
    edgeKeys[base + 3] = EdgeKey::encode(n[1], n[2]);
    edgeKeys[base + 4] = EdgeKey::encode(n[1], n[3]);
    edgeKeys[base + 5] = EdgeKey::encode(n[2], n[3]);
}

template<typename KeyType, typename RealType>
__global__ void computeTetEdgeMidpointsKernel(const uint64_t* uniqueEdges,
                                              const RealType* nodeX,
                                              const RealType* nodeY,
                                              const RealType* nodeZ,
                                              RealType* newX,
                                              RealType* newY,
                                              RealType* newZ,
                                              KeyType baseNodeId,
                                              size_t numEdges)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numEdges) return;

    uint64_t key = uniqueEdges[i];
    KeyType a    = key & 0xFFFFFFFF;
    KeyType b    = key >> 32;

    KeyType idx  = baseNodeId + i;
    newX[idx]    = RealType(0.5) * (nodeX[a] + nodeX[b]);
    newY[idx]    = RealType(0.5) * (nodeY[a] + nodeY[b]);
    newZ[idx]    = RealType(0.5) * (nodeZ[a] + nodeZ[b]);
}

// Build 8 child tets per refined element using Bey's red refinement.
// Unrefinied elements are copied as-is.
template<typename KeyType>
__global__ void buildChildTetConnectivityKernel(const KeyType* conn0,
                                                const KeyType* conn1,
                                                const KeyType* conn2,
                                                const KeyType* conn3,
                                                const uint8_t* marks,
                                                const uint64_t* elemPrefix,
                                                const uint64_t* sortedEdgeKeys,
                                                size_t numUniqueEdges,
                                                KeyType edgeBaseNode,
                                                KeyType* outConn0,
                                                KeyType* outConn1,
                                                KeyType* outConn2,
                                                KeyType* outConn3,
                                                size_t numElements)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements) return;

    uint64_t outBase = elemPrefix[e];
    KeyType n[4]     = {conn0[e], conn1[e], conn2[e], conn3[e]};

    if (marks[e] <= 0)
    {
        outConn0[outBase] = n[0];
        outConn1[outBase] = n[1];
        outConn2[outBase] = n[2];
        outConn3[outBase] = n[3];
        return;
    }

    auto findEdge = [&](KeyType a, KeyType b) -> KeyType
    {
        uint64_t key = EdgeKey::encode(a, b);
        size_t lo = 0, hi = numUniqueEdges;
        while (lo < hi)
        {
            size_t mid = (lo + hi) / 2;
            if (sortedEdgeKeys[mid] < key) lo = mid + 1;
            else hi = mid;
        }
        return edgeBaseNode + lo;
    };

    // 6 edge midpoints
    KeyType m01 = findEdge(n[0], n[1]);
    KeyType m02 = findEdge(n[0], n[2]);
    KeyType m03 = findEdge(n[0], n[3]);
    KeyType m12 = findEdge(n[1], n[2]);
    KeyType m13 = findEdge(n[1], n[3]);
    KeyType m23 = findEdge(n[2], n[3]);

    auto writeChild = [&](int child, KeyType c0, KeyType c1, KeyType c2, KeyType c3)
    {
        uint64_t idx     = outBase + child;
        outConn0[idx] = c0;
        outConn1[idx] = c1;
        outConn2[idx] = c2;
        outConn3[idx] = c3;
    };

    // Bey's red refinement: 4 corner tets + 4 octahedron tets
    // Corner tets: each original vertex + its 3 adjacent midpoints
    writeChild(0, n[0], m01, m02, m03);
    writeChild(1, m01, n[1], m12, m13);
    writeChild(2, m02, m12, n[2], m23);
    writeChild(3, m03, m13, m23, n[3]);

    // Interior octahedron: 4 tets (consistent diagonal choice: m01-m23)
    // This diagonal choice keeps quality bounded across refinement levels
    writeChild(4, m01, m02, m03, m13);
    writeChild(5, m01, m02, m13, m12);
    writeChild(6, m02, m03, m13, m23);
    writeChild(7, m02, m12, m13, m23);
}

template<typename KeyType, typename RealType>
class TetRefiner
{
    using DeviceVector    = cstone::DeviceVector<KeyType>;
    using RealDeviceVector = cstone::DeviceVector<RealType>;

public:
    struct Result
    {
        RealDeviceVector d_x, d_y, d_z;
        DeviceVector d_conn0, d_conn1, d_conn2, d_conn3;
        size_t numNodes    = 0;
        size_t numElements = 0;

        // Kept for parentage-based solution transfer
        cstone::DeviceVector<uint64_t> d_edgeKeys;
        cstone::DeviceVector<uint64_t> d_markedPrefix;
        cstone::DeviceVector<uint8_t> d_marks;
        size_t numUniqueEdges = 0;
        size_t numMarked      = 0;
        size_t oldNumNodes    = 0;

        KeyType edgeBaseNode() const { return oldNumNodes; }
    };

    static Result refine(const KeyType* conn0,
                         const KeyType* conn1,
                         const KeyType* conn2,
                         const KeyType* conn3,
                         const uint8_t* marks,
                         size_t numElements,
                         const RealType* nodeX,
                         const RealType* nodeY,
                         const RealType* nodeZ,
                         size_t numNodes,
                         int blockSize = 256)
    {
        int numBlocks = (numElements + blockSize - 1) / blockSize;

        // Step 1: count children per element
        cstone::DeviceVector<uint64_t> d_childCounts(numElements);
        countChildTetsKernel<<<numBlocks, blockSize>>>(marks, d_childCounts.data(), numElements);

        cstone::DeviceVector<uint64_t> d_elemPrefix(numElements);
        auto cc_b = thrust::device_pointer_cast(d_childCounts.data());
        auto cc_e = thrust::device_pointer_cast(d_childCounts.data() + numElements);
        auto ep_b = thrust::device_pointer_cast(d_elemPrefix.data());
        thrust::exclusive_scan(thrust::device, cc_b, cc_e, ep_b);

        uint64_t totalNewElements;
        {
            uint64_t lastCount;
            cudaMemcpy(&totalNewElements, d_elemPrefix.data() + numElements - 1, sizeof(uint64_t), cudaMemcpyDeviceToHost);
            cudaMemcpy(&lastCount, d_childCounts.data() + numElements - 1, sizeof(uint64_t), cudaMemcpyDeviceToHost);
            totalNewElements += lastCount;
        }

        // Step 2: count marked elements
        cstone::DeviceVector<uint64_t> d_isMarked(numElements);
        auto im_b = thrust::device_pointer_cast(d_isMarked.data());
        auto im_e = thrust::device_pointer_cast(d_isMarked.data() + numElements);
        thrust::transform(thrust::device, marks, marks + numElements, im_b,
                           [] __device__(uint8_t m) -> uint64_t { return (m > 0) ? 1 : 0; });

        cstone::DeviceVector<uint64_t> d_markedPrefix(numElements);
        auto mp_b = thrust::device_pointer_cast(d_markedPrefix.data());
        thrust::exclusive_scan(thrust::device, im_b, im_e, mp_b);

        uint64_t numMarked;
        {
            uint64_t last, lastIs;
            cudaMemcpy(&last, d_markedPrefix.data() + numElements - 1, sizeof(uint64_t), cudaMemcpyDeviceToHost);
            cudaMemcpy(&lastIs, d_isMarked.data() + numElements - 1, sizeof(uint64_t), cudaMemcpyDeviceToHost);
            numMarked = last + lastIs;
        }

        if (numMarked == 0)
        {
            Result res;
            res.numNodes    = numNodes;
            res.numElements = numElements;
            res.d_x.resize(numNodes); res.d_y.resize(numNodes); res.d_z.resize(numNodes);
            cudaMemcpy(res.d_x.data(), nodeX, numNodes * sizeof(RealType), cudaMemcpyDeviceToDevice);
            cudaMemcpy(res.d_y.data(), nodeY, numNodes * sizeof(RealType), cudaMemcpyDeviceToDevice);
            cudaMemcpy(res.d_z.data(), nodeZ, numNodes * sizeof(RealType), cudaMemcpyDeviceToDevice);
            res.d_conn0.resize(numElements); res.d_conn1.resize(numElements);
            res.d_conn2.resize(numElements); res.d_conn3.resize(numElements);
            cudaMemcpy(res.d_conn0.data(), conn0, numElements * sizeof(KeyType), cudaMemcpyDeviceToDevice);
            cudaMemcpy(res.d_conn1.data(), conn1, numElements * sizeof(KeyType), cudaMemcpyDeviceToDevice);
            cudaMemcpy(res.d_conn2.data(), conn2, numElements * sizeof(KeyType), cudaMemcpyDeviceToDevice);
            cudaMemcpy(res.d_conn3.data(), conn3, numElements * sizeof(KeyType), cudaMemcpyDeviceToDevice);
            res.oldNumNodes = numNodes;
            return res;
        }

        // Step 3: emit and dedup edge keys
        size_t totalEdgeKeys = numMarked * 6;
        cstone::DeviceVector<uint64_t> d_edgeKeys(totalEdgeKeys);
        emitTetEdgeKeysKernel<<<numBlocks, blockSize>>>(conn0, conn1, conn2, conn3, marks,
                                                         d_markedPrefix.data(), d_edgeKeys.data(), numElements);

        auto ek_b = thrust::device_pointer_cast(d_edgeKeys.data());
        auto ek_e = thrust::device_pointer_cast(d_edgeKeys.data() + totalEdgeKeys);
        thrust::sort(thrust::device, ek_b, ek_e);
        auto edgeEnd          = thrust::unique(thrust::device, ek_b, ek_e);
        size_t numUniqueEdges = edgeEnd - ek_b;
        d_edgeKeys.resize(numUniqueEdges);

        // Step 4: allocate and compute new coordinates
        // Tets only produce edge midpoints (no face/body centers)
        size_t numNewNodes   = numNodes + numUniqueEdges;
        KeyType edgeBaseNode = numNodes;

        Result res;
        res.numNodes    = numNewNodes;
        res.numElements = totalNewElements;

        res.d_x.resize(numNewNodes); res.d_y.resize(numNewNodes); res.d_z.resize(numNewNodes);
        cudaMemcpy(res.d_x.data(), nodeX, numNodes * sizeof(RealType), cudaMemcpyDeviceToDevice);
        cudaMemcpy(res.d_y.data(), nodeY, numNodes * sizeof(RealType), cudaMemcpyDeviceToDevice);
        cudaMemcpy(res.d_z.data(), nodeZ, numNodes * sizeof(RealType), cudaMemcpyDeviceToDevice);

        int edgeBlocks = (numUniqueEdges + blockSize - 1) / blockSize;
        computeTetEdgeMidpointsKernel<KeyType, RealType><<<edgeBlocks, blockSize>>>(
            d_edgeKeys.data(), nodeX, nodeY, nodeZ,
            res.d_x.data(), res.d_y.data(), res.d_z.data(), edgeBaseNode, numUniqueEdges);

        // Step 5: build child connectivity
        res.d_conn0.resize(totalNewElements); res.d_conn1.resize(totalNewElements);
        res.d_conn2.resize(totalNewElements); res.d_conn3.resize(totalNewElements);

        buildChildTetConnectivityKernel<<<numBlocks, blockSize>>>(
            conn0, conn1, conn2, conn3, marks, d_elemPrefix.data(),
            d_edgeKeys.data(), numUniqueEdges, edgeBaseNode,
            res.d_conn0.data(), res.d_conn1.data(), res.d_conn2.data(), res.d_conn3.data(),
            numElements);

        cudaDeviceSynchronize();

        res.d_edgeKeys     = std::move(d_edgeKeys);
        res.d_markedPrefix = std::move(d_markedPrefix);
        res.numUniqueEdges = numUniqueEdges;
        res.numMarked      = numMarked;
        res.oldNumNodes    = numNodes;

        res.d_marks.resize(numElements);
        cudaMemcpy(res.d_marks.data(), marks, numElements * sizeof(uint8_t), cudaMemcpyDeviceToDevice);

        return res;
    }
};

} // namespace amr
} // namespace mars
