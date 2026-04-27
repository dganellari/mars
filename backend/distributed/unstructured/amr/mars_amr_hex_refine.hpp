#pragma once

#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/copy.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/sequence.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>
#include <thrust/reduce.h>
#include <thrust/scatter.h>
#include <thrust/gather.h>
#include <thrust/binary_search.h>
#include <thrust/transform.h>

#include "cstone/cuda/cuda_utils.hpp"

namespace mars
{
namespace amr
{

// Node numbering convention (hex8, same as MFEM/Exodus):
//
//        7--------6
//       /|       /|
//      / |      / |
//     4--------5  |
//     |  3-----|--2
//     | /      | /
//     |/       |/
//     0--------1
//
// Refinement inserts up to 19 new nodes per element (shared via dedup):
//   - 12 edge midpoints
//   - 6 face centers
//   - 1 body center
// Producing 8 child hexes per refined element.

// Ordered pair packed into uint64: lo in bits [0:31], hi in bits [32:63].
// Ordering ensures the same edge always produces the same key regardless of direction.
struct EdgeKey
{
    __host__ __device__ static uint64_t encode(uint64_t a, uint64_t b)
    {
        uint64_t lo = (a < b) ? a : b;
        uint64_t hi = (a < b) ? b : a;
        // Breaks for >4B nodes; would need a different encoding
        return (hi << 32) | (lo & 0xFFFFFFFF);
    }
};

// Face key: 4 node IDs sorted then packed into 2 uint64s.
// Sorting makes the key independent of face winding order.
struct FaceKey
{
    __host__ __device__ static void sort4(uint64_t& a, uint64_t& b, uint64_t& c, uint64_t& d)
    {
        if (a > b) { uint64_t t = a; a = b; b = t; }
        if (c > d) { uint64_t t = c; c = d; d = t; }
        if (a > c) { uint64_t t = a; a = c; c = t; }
        if (b > d) { uint64_t t = b; b = d; d = t; }
        if (b > c) { uint64_t t = b; b = c; c = t; }
    }

    __host__ __device__ static uint64_t encodeLo(uint64_t a, uint64_t b, uint64_t c, uint64_t d)
    {
        sort4(a, b, c, d);
        return (a << 32) | (b & 0xFFFFFFFF);
    }

    __host__ __device__ static uint64_t encodeHi(uint64_t a, uint64_t b, uint64_t c, uint64_t d)
    {
        sort4(a, b, c, d);
        return (c << 32) | (d & 0xFFFFFFFF);
    }
};

__global__ void countChildElementsKernel(const uint8_t* marks, uint64_t* childCounts, size_t numElements)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements) return;
    childCounts[e] = (marks[e] > 0) ? 8 : 1;
}

template<typename KeyType>
__global__ void emitEdgeKeysKernel(const KeyType* conn0,
                                   const KeyType* conn1,
                                   const KeyType* conn2,
                                   const KeyType* conn3,
                                   const KeyType* conn4,
                                   const KeyType* conn5,
                                   const KeyType* conn6,
                                   const KeyType* conn7,
                                   const uint8_t* marks,
                                   const uint64_t* markedPrefix,
                                   uint64_t* edgeKeys,
                                   size_t numElements)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements || marks[e] <= 0) return;

    uint64_t base = markedPrefix[e] * 12; // 12 edges per marked element
    KeyType n[8]  = {conn0[e], conn1[e], conn2[e], conn3[e], conn4[e], conn5[e], conn6[e], conn7[e]};

    // Bottom face
    edgeKeys[base + 0]  = EdgeKey::encode(n[0], n[1]);
    edgeKeys[base + 1]  = EdgeKey::encode(n[1], n[2]);
    edgeKeys[base + 2]  = EdgeKey::encode(n[2], n[3]);
    edgeKeys[base + 3]  = EdgeKey::encode(n[3], n[0]);
    // Top face
    edgeKeys[base + 4]  = EdgeKey::encode(n[4], n[5]);
    edgeKeys[base + 5]  = EdgeKey::encode(n[5], n[6]);
    edgeKeys[base + 6]  = EdgeKey::encode(n[6], n[7]);
    edgeKeys[base + 7]  = EdgeKey::encode(n[7], n[4]);
    // Vertical
    edgeKeys[base + 8]  = EdgeKey::encode(n[0], n[4]);
    edgeKeys[base + 9]  = EdgeKey::encode(n[1], n[5]);
    edgeKeys[base + 10] = EdgeKey::encode(n[2], n[6]);
    edgeKeys[base + 11] = EdgeKey::encode(n[3], n[7]);
}

template<typename KeyType>
__global__ void emitFaceKeysKernel(const KeyType* conn0,
                                   const KeyType* conn1,
                                   const KeyType* conn2,
                                   const KeyType* conn3,
                                   const KeyType* conn4,
                                   const KeyType* conn5,
                                   const KeyType* conn6,
                                   const KeyType* conn7,
                                   const uint8_t* marks,
                                   const uint64_t* markedPrefix,
                                   uint64_t* faceKeysLo,
                                   uint64_t* faceKeysHi,
                                   size_t numElements)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements || marks[e] <= 0) return;

    uint64_t base = markedPrefix[e] * 6;
    KeyType n[8]  = {conn0[e], conn1[e], conn2[e], conn3[e], conn4[e], conn5[e], conn6[e], conn7[e]};

    // 6 faces of hex8
    // Bottom: 0,1,2,3
    faceKeysLo[base + 0] = FaceKey::encodeLo(n[0], n[1], n[2], n[3]);
    faceKeysHi[base + 0] = FaceKey::encodeHi(n[0], n[1], n[2], n[3]);
    // Top: 4,5,6,7
    faceKeysLo[base + 1] = FaceKey::encodeLo(n[4], n[5], n[6], n[7]);
    faceKeysHi[base + 1] = FaceKey::encodeHi(n[4], n[5], n[6], n[7]);
    // Front: 0,1,5,4
    faceKeysLo[base + 2] = FaceKey::encodeLo(n[0], n[1], n[5], n[4]);
    faceKeysHi[base + 2] = FaceKey::encodeHi(n[0], n[1], n[5], n[4]);
    // Back: 3,2,6,7
    faceKeysLo[base + 3] = FaceKey::encodeLo(n[3], n[2], n[6], n[7]);
    faceKeysHi[base + 3] = FaceKey::encodeHi(n[3], n[2], n[6], n[7]);
    // Left: 0,3,7,4
    faceKeysLo[base + 4] = FaceKey::encodeLo(n[0], n[3], n[7], n[4]);
    faceKeysHi[base + 4] = FaceKey::encodeHi(n[0], n[3], n[7], n[4]);
    // Right: 1,2,6,5
    faceKeysLo[base + 5] = FaceKey::encodeLo(n[1], n[2], n[6], n[5]);
    faceKeysHi[base + 5] = FaceKey::encodeHi(n[1], n[2], n[6], n[5]);
}

template<typename KeyType, typename RealType>
__global__ void computeEdgeMidpointsKernel(const uint64_t* uniqueEdges,
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

    KeyType newIdx        = baseNodeId + i;
    newX[newIdx] = RealType(0.5) * (nodeX[a] + nodeX[b]);
    newY[newIdx] = RealType(0.5) * (nodeY[a] + nodeY[b]);
    newZ[newIdx] = RealType(0.5) * (nodeZ[a] + nodeZ[b]);
}

template<typename KeyType, typename RealType>
__global__ void computeFaceCentersKernel(const uint64_t* uniqueFaceLo,
                                         const uint64_t* uniqueFaceHi,
                                         const RealType* nodeX,
                                         const RealType* nodeY,
                                         const RealType* nodeZ,
                                         RealType* newX,
                                         RealType* newY,
                                         RealType* newZ,
                                         KeyType baseNodeId,
                                         size_t numFaces)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numFaces) return;

    KeyType a = uniqueFaceLo[i] >> 32;
    KeyType b = uniqueFaceLo[i] & 0xFFFFFFFF;
    KeyType c = uniqueFaceHi[i] >> 32;
    KeyType d = uniqueFaceHi[i] & 0xFFFFFFFF;

    KeyType newIdx        = baseNodeId + i;
    newX[newIdx] = RealType(0.25) * (nodeX[a] + nodeX[b] + nodeX[c] + nodeX[d]);
    newY[newIdx] = RealType(0.25) * (nodeY[a] + nodeY[b] + nodeY[c] + nodeY[d]);
    newZ[newIdx] = RealType(0.25) * (nodeZ[a] + nodeZ[b] + nodeZ[c] + nodeZ[d]);
}

template<typename KeyType, typename RealType>
__global__ void computeBodyCentersKernel(const KeyType* conn0,
                                         const KeyType* conn1,
                                         const KeyType* conn2,
                                         const KeyType* conn3,
                                         const KeyType* conn4,
                                         const KeyType* conn5,
                                         const KeyType* conn6,
                                         const KeyType* conn7,
                                         const uint8_t* marks,
                                         const uint64_t* markedPrefix,
                                         const RealType* nodeX,
                                         const RealType* nodeY,
                                         const RealType* nodeZ,
                                         RealType* newX,
                                         RealType* newY,
                                         RealType* newZ,
                                         KeyType baseNodeId,
                                         size_t numElements)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements || marks[e] <= 0) return;

    KeyType n[8] = {conn0[e], conn1[e], conn2[e], conn3[e], conn4[e], conn5[e], conn6[e], conn7[e]};

    KeyType newIdx = baseNodeId + markedPrefix[e];
    RealType inv8  = RealType(0.125);
    newX[newIdx]   = inv8 * (nodeX[n[0]] + nodeX[n[1]] + nodeX[n[2]] + nodeX[n[3]] + nodeX[n[4]] + nodeX[n[5]] +
                           nodeX[n[6]] + nodeX[n[7]]);
    newY[newIdx]   = inv8 * (nodeY[n[0]] + nodeY[n[1]] + nodeY[n[2]] + nodeY[n[3]] + nodeY[n[4]] + nodeY[n[5]] +
                           nodeY[n[6]] + nodeY[n[7]]);
    newZ[newIdx]   = inv8 * (nodeZ[n[0]] + nodeZ[n[1]] + nodeZ[n[2]] + nodeZ[n[3]] + nodeZ[n[4]] + nodeZ[n[5]] +
                           nodeZ[n[6]] + nodeZ[n[7]]);
}

template<typename KeyType>
__global__ void buildChildConnectivityKernel(const KeyType* conn0,
                                             const KeyType* conn1,
                                             const KeyType* conn2,
                                             const KeyType* conn3,
                                             const KeyType* conn4,
                                             const KeyType* conn5,
                                             const KeyType* conn6,
                                             const KeyType* conn7,
                                             const uint8_t* marks,
                                             const uint64_t* elemPrefix,
                                             const uint64_t* markedPrefix,
                                             const uint64_t* sortedEdgeKeys,
                                             size_t numUniqueEdges,
                                             KeyType edgeBaseNode,
                                             const uint64_t* sortedFaceKeysLo,
                                             const uint64_t* sortedFaceKeysHi,
                                             size_t numUniqueFaces,
                                             KeyType faceBaseNode,
                                             KeyType bodyBaseNode,
                                             KeyType* outConn0,
                                             KeyType* outConn1,
                                             KeyType* outConn2,
                                             KeyType* outConn3,
                                             KeyType* outConn4,
                                             KeyType* outConn5,
                                             KeyType* outConn6,
                                             KeyType* outConn7,
                                             size_t numElements)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements) return;

    uint64_t outBase = elemPrefix[e];
    KeyType n[8] = {conn0[e], conn1[e], conn2[e], conn3[e], conn4[e], conn5[e], conn6[e], conn7[e]};

    if (marks[e] <= 0)
    {
        outConn0[outBase] = n[0];
        outConn1[outBase] = n[1];
        outConn2[outBase] = n[2];
        outConn3[outBase] = n[3];
        outConn4[outBase] = n[4];
        outConn5[outBase] = n[5];
        outConn6[outBase] = n[6];
        outConn7[outBase] = n[7];
        return;
    }

    auto findEdge = [&](KeyType a, KeyType b) -> KeyType
    {
        uint64_t key = EdgeKey::encode(a, b);
        size_t lo = 0, hi = numUniqueEdges;
        while (lo < hi)
        {
            size_t mid = (lo + hi) / 2;
            if (sortedEdgeKeys[mid] < key)
                lo = mid + 1;
            else
                hi = mid;
        }
        return edgeBaseNode + lo;
    };

    auto findFace = [&](KeyType a, KeyType b, KeyType c, KeyType d) -> KeyType
    {
        uint64_t kLo = FaceKey::encodeLo(a, b, c, d);
        uint64_t kHi = FaceKey::encodeHi(a, b, c, d);
        size_t lo = 0, hi = numUniqueFaces;
        while (lo < hi)
        {
            size_t mid = (lo + hi) / 2;
            if (sortedFaceKeysLo[mid] < kLo || (sortedFaceKeysLo[mid] == kLo && sortedFaceKeysHi[mid] < kHi))
                lo = mid + 1;
            else
                hi = mid;
        }
        return faceBaseNode + lo;
    };

    KeyType e01 = findEdge(n[0], n[1]);
    KeyType e12 = findEdge(n[1], n[2]);
    KeyType e23 = findEdge(n[2], n[3]);
    KeyType e30 = findEdge(n[3], n[0]);
    KeyType e45 = findEdge(n[4], n[5]);
    KeyType e56 = findEdge(n[5], n[6]);
    KeyType e67 = findEdge(n[6], n[7]);
    KeyType e74 = findEdge(n[7], n[4]);
    KeyType e04 = findEdge(n[0], n[4]);
    KeyType e15 = findEdge(n[1], n[5]);
    KeyType e26 = findEdge(n[2], n[6]);
    KeyType e37 = findEdge(n[3], n[7]);

    KeyType fBot  = findFace(n[0], n[1], n[2], n[3]);
    KeyType fTop  = findFace(n[4], n[5], n[6], n[7]);
    KeyType fFrnt = findFace(n[0], n[1], n[5], n[4]);
    KeyType fBack = findFace(n[3], n[2], n[6], n[7]);
    KeyType fLeft = findFace(n[0], n[3], n[7], n[4]);
    KeyType fRght = findFace(n[1], n[2], n[6], n[5]);

    KeyType bCtr = bodyBaseNode + markedPrefix[e];

    auto writeChild = [&](int child, KeyType c0, KeyType c1, KeyType c2, KeyType c3, KeyType c4, KeyType c5,
                          KeyType c6, KeyType c7)
    {
        uint64_t idx     = outBase + child;
        outConn0[idx] = c0;
        outConn1[idx] = c1;
        outConn2[idx] = c2;
        outConn3[idx] = c3;
        outConn4[idx] = c4;
        outConn5[idx] = c5;
        outConn6[idx] = c6;
        outConn7[idx] = c7;
    };

    writeChild(0, n[0], e01, fBot, e30, e04, fFrnt, bCtr, fLeft);
    writeChild(1, e01, n[1], e12, fBot, fFrnt, e15, fRght, bCtr);
    writeChild(2, fBot, e12, n[2], e23, bCtr, fRght, e26, fBack);
    writeChild(3, e30, fBot, e23, n[3], fLeft, bCtr, fBack, e37);
    writeChild(4, e04, fFrnt, bCtr, fLeft, n[4], e45, fTop, e74);
    writeChild(5, fFrnt, e15, fRght, bCtr, e45, n[5], e56, fTop);
    writeChild(6, bCtr, fRght, e26, fBack, fTop, e56, n[6], e67);
    writeChild(7, fLeft, bCtr, fBack, e37, e74, fTop, e67, n[7]);
}

template<typename KeyType, typename RealType>
class HexRefiner
{
    using DeviceVector = cstone::DeviceVector<KeyType>;
    using RealDeviceVector = cstone::DeviceVector<RealType>;

public:
    struct Result
    {
        RealDeviceVector d_x, d_y, d_z;
        DeviceVector d_conn0, d_conn1, d_conn2, d_conn3;
        DeviceVector d_conn4, d_conn5, d_conn6, d_conn7;
        size_t numNodes    = 0;
        size_t numElements = 0;

        // Kept for parentage-based solution transfer after refinement
        cstone::DeviceVector<uint64_t> d_edgeKeys;
        cstone::DeviceVector<uint64_t> d_faceKeysLo;
        cstone::DeviceVector<uint64_t> d_faceKeysHi;
        cstone::DeviceVector<uint64_t> d_markedPrefix;
        cstone::DeviceVector<uint8_t> d_marks;
        size_t numUniqueEdges = 0;
        size_t numUniqueFaces = 0;
        size_t numMarked      = 0;
        size_t oldNumNodes    = 0; // node count before refinement

        // Base node IDs for each type of new node
        KeyType edgeBaseNode() const { return oldNumNodes; }
        KeyType faceBaseNode() const { return oldNumNodes + numUniqueEdges; }
        KeyType bodyBaseNode() const { return oldNumNodes + numUniqueEdges + numUniqueFaces; }
    };

    static Result refine(const KeyType* conn0,
                         const KeyType* conn1,
                         const KeyType* conn2,
                         const KeyType* conn3,
                         const KeyType* conn4,
                         const KeyType* conn5,
                         const KeyType* conn6,
                         const KeyType* conn7,
                         const uint8_t* marks,
                         size_t numElements,
                         const RealType* nodeX,
                         const RealType* nodeY,
                         const RealType* nodeZ,
                         size_t numNodes,
                         int blockSize = 256)
    {
        int numBlocks = (numElements + blockSize - 1) / blockSize;

        // ---- Step 1: Count children per element ----
        cstone::DeviceVector<uint64_t> d_childCounts(numElements);
        countChildElementsKernel<<<numBlocks, blockSize>>>(marks, d_childCounts.data(), numElements);

        cstone::DeviceVector<uint64_t> d_elemPrefix(numElements);
        auto cc_b  = thrust::device_pointer_cast(d_childCounts.data());
        auto cc_e  = thrust::device_pointer_cast(d_childCounts.data() + numElements);
        auto ep_b  = thrust::device_pointer_cast(d_elemPrefix.data());
        thrust::exclusive_scan(thrust::device, cc_b, cc_e, ep_b);

        uint64_t totalNewElements;
        cudaMemcpy(&totalNewElements, d_elemPrefix.data() + numElements - 1, sizeof(uint64_t),
                    cudaMemcpyDeviceToHost);
        uint64_t lastCount;
        cudaMemcpy(&lastCount, d_childCounts.data() + numElements - 1, sizeof(uint64_t), cudaMemcpyDeviceToHost);
        totalNewElements += lastCount;

        // ---- Step 2: Count marked elements for indexing ----
        cstone::DeviceVector<uint64_t> d_isMarked(numElements);
        auto im_b  = thrust::device_pointer_cast(d_isMarked.data());
        auto im_e  = thrust::device_pointer_cast(d_isMarked.data() + numElements);
        thrust::transform(thrust::device, marks, marks + numElements, im_b,
                           [] __device__(uint8_t m) -> uint64_t { return (m > 0) ? 1 : 0; });

        cstone::DeviceVector<uint64_t> d_markedPrefix(numElements);
        auto mp_b  = thrust::device_pointer_cast(d_markedPrefix.data());
        thrust::exclusive_scan(thrust::device, im_b, im_e, mp_b);

        uint64_t numMarked;
        {
            uint64_t lastMarked, lastIsMarked;
            cudaMemcpy(&lastMarked, d_markedPrefix.data() + numElements - 1, sizeof(uint64_t),
                        cudaMemcpyDeviceToHost);
            cudaMemcpy(&lastIsMarked, d_isMarked.data() + numElements - 1, sizeof(uint64_t), cudaMemcpyDeviceToHost);
            numMarked = lastMarked + lastIsMarked;
        }

        if (numMarked == 0)
        {
            // Nothing to refine
            Result res;
            res.numNodes    = numNodes;
            res.numElements = numElements;
            res.d_x.resize(numNodes);
            res.d_y.resize(numNodes);
            res.d_z.resize(numNodes);
            cudaMemcpy(res.d_x.data(), nodeX, numNodes * sizeof(RealType), cudaMemcpyDeviceToDevice);
            cudaMemcpy(res.d_y.data(), nodeY, numNodes * sizeof(RealType), cudaMemcpyDeviceToDevice);
            cudaMemcpy(res.d_z.data(), nodeZ, numNodes * sizeof(RealType), cudaMemcpyDeviceToDevice);
            res.d_conn0.resize(numElements);
            res.d_conn1.resize(numElements);
            res.d_conn2.resize(numElements);
            res.d_conn3.resize(numElements);
            res.d_conn4.resize(numElements);
            res.d_conn5.resize(numElements);
            res.d_conn6.resize(numElements);
            res.d_conn7.resize(numElements);
            cudaMemcpy(res.d_conn0.data(), conn0, numElements * sizeof(KeyType), cudaMemcpyDeviceToDevice);
            cudaMemcpy(res.d_conn1.data(), conn1, numElements * sizeof(KeyType), cudaMemcpyDeviceToDevice);
            cudaMemcpy(res.d_conn2.data(), conn2, numElements * sizeof(KeyType), cudaMemcpyDeviceToDevice);
            cudaMemcpy(res.d_conn3.data(), conn3, numElements * sizeof(KeyType), cudaMemcpyDeviceToDevice);
            cudaMemcpy(res.d_conn4.data(), conn4, numElements * sizeof(KeyType), cudaMemcpyDeviceToDevice);
            cudaMemcpy(res.d_conn5.data(), conn5, numElements * sizeof(KeyType), cudaMemcpyDeviceToDevice);
            cudaMemcpy(res.d_conn6.data(), conn6, numElements * sizeof(KeyType), cudaMemcpyDeviceToDevice);
            cudaMemcpy(res.d_conn7.data(), conn7, numElements * sizeof(KeyType), cudaMemcpyDeviceToDevice);
            return res;
        }

        // ---- Step 3: Emit and deduplicate edge keys ----
        size_t totalEdgeKeys = numMarked * 12;
        cstone::DeviceVector<uint64_t> d_edgeKeys(totalEdgeKeys);
        emitEdgeKeysKernel<<<numBlocks, blockSize>>>(conn0, conn1, conn2, conn3, conn4, conn5, conn6, conn7, marks,
                                                      d_markedPrefix.data(), d_edgeKeys.data(), numElements);

        auto ek_b = thrust::device_pointer_cast(d_edgeKeys.data());
        auto ek_e = thrust::device_pointer_cast(d_edgeKeys.data() + totalEdgeKeys);
        thrust::sort(thrust::device, ek_b, ek_e);
        auto edgeEnd        = thrust::unique(thrust::device, ek_b, ek_e);
        size_t numUniqueEdges = edgeEnd - ek_b;
        d_edgeKeys.resize(numUniqueEdges);

        // ---- Step 4: Emit and deduplicate face keys ----
        size_t totalFaceKeys = numMarked * 6;
        cstone::DeviceVector<uint64_t> d_faceKeysLo(totalFaceKeys);
        cstone::DeviceVector<uint64_t> d_faceKeysHi(totalFaceKeys);
        emitFaceKeysKernel<<<numBlocks, blockSize>>>(conn0, conn1, conn2, conn3, conn4, conn5, conn6, conn7, marks,
                                                      d_markedPrefix.data(), d_faceKeysLo.data(), d_faceKeysHi.data(),
                                                      numElements);

        auto fl_b = thrust::device_pointer_cast(d_faceKeysLo.data());
        auto fl_e = thrust::device_pointer_cast(d_faceKeysLo.data() + totalFaceKeys);
        auto fh_b = thrust::device_pointer_cast(d_faceKeysHi.data());
        auto fh_e = thrust::device_pointer_cast(d_faceKeysHi.data() + totalFaceKeys);
        thrust::sort(thrust::device,
            thrust::make_zip_iterator(thrust::make_tuple(fl_b, fh_b)),
            thrust::make_zip_iterator(thrust::make_tuple(fl_e, fh_e)));

        auto faceEnd = thrust::unique(thrust::device,
            thrust::make_zip_iterator(thrust::make_tuple(fl_b, fh_b)),
            thrust::make_zip_iterator(thrust::make_tuple(fl_e, fh_e)));

        size_t numUniqueFaces = faceEnd - thrust::make_zip_iterator(thrust::make_tuple(fl_b, fh_b));
        d_faceKeysLo.resize(numUniqueFaces);
        d_faceKeysHi.resize(numUniqueFaces);

        // ---- Step 5: Allocate new coordinate arrays ----
        size_t numNewNodes = numNodes + numUniqueEdges + numUniqueFaces + numMarked;
        KeyType edgeBaseNode = numNodes;
        KeyType faceBaseNode = numNodes + numUniqueEdges;
        KeyType bodyBaseNode = numNodes + numUniqueEdges + numUniqueFaces;

        Result res;
        res.numNodes    = numNewNodes;
        res.numElements = totalNewElements;

        res.d_x.resize(numNewNodes);
        res.d_y.resize(numNewNodes);
        res.d_z.resize(numNewNodes);

        cudaMemcpy(res.d_x.data(), nodeX, numNodes * sizeof(RealType), cudaMemcpyDeviceToDevice);
        cudaMemcpy(res.d_y.data(), nodeY, numNodes * sizeof(RealType), cudaMemcpyDeviceToDevice);
        cudaMemcpy(res.d_z.data(), nodeZ, numNodes * sizeof(RealType), cudaMemcpyDeviceToDevice);

        // ---- Step 6: Compute new node coordinates ----
        int edgeBlocks = (numUniqueEdges + blockSize - 1) / blockSize;
        computeEdgeMidpointsKernel<KeyType, RealType><<<edgeBlocks, blockSize>>>(
            d_edgeKeys.data(), nodeX, nodeY, nodeZ, res.d_x.data(), res.d_y.data(), res.d_z.data(), edgeBaseNode,
            numUniqueEdges);

        int faceBlocks = (numUniqueFaces + blockSize - 1) / blockSize;
        computeFaceCentersKernel<KeyType, RealType><<<faceBlocks, blockSize>>>(
            d_faceKeysLo.data(), d_faceKeysHi.data(), nodeX, nodeY, nodeZ, res.d_x.data(), res.d_y.data(),
            res.d_z.data(), faceBaseNode, numUniqueFaces);

        computeBodyCentersKernel<<<numBlocks, blockSize>>>(conn0, conn1, conn2, conn3, conn4, conn5, conn6, conn7,
                                                            marks, d_markedPrefix.data(), nodeX, nodeY, nodeZ,
                                                            res.d_x.data(), res.d_y.data(), res.d_z.data(),
                                                            bodyBaseNode, numElements);

        // ---- Step 7: Build child element connectivity ----
        res.d_conn0.resize(totalNewElements);
        res.d_conn1.resize(totalNewElements);
        res.d_conn2.resize(totalNewElements);
        res.d_conn3.resize(totalNewElements);
        res.d_conn4.resize(totalNewElements);
        res.d_conn5.resize(totalNewElements);
        res.d_conn6.resize(totalNewElements);
        res.d_conn7.resize(totalNewElements);

        buildChildConnectivityKernel<<<numBlocks, blockSize>>>(
            conn0, conn1, conn2, conn3, conn4, conn5, conn6, conn7, marks, d_elemPrefix.data(),
            d_markedPrefix.data(), d_edgeKeys.data(), numUniqueEdges, edgeBaseNode, d_faceKeysLo.data(),
            d_faceKeysHi.data(), numUniqueFaces, faceBaseNode, bodyBaseNode, res.d_conn0.data(), res.d_conn1.data(),
            res.d_conn2.data(), res.d_conn3.data(), res.d_conn4.data(), res.d_conn5.data(), res.d_conn6.data(),
            res.d_conn7.data(), numElements);

        cudaDeviceSynchronize();

        res.d_edgeKeys    = std::move(d_edgeKeys);
        res.d_faceKeysLo  = std::move(d_faceKeysLo);
        res.d_faceKeysHi  = std::move(d_faceKeysHi);
        res.d_markedPrefix = std::move(d_markedPrefix);
        res.numUniqueEdges = numUniqueEdges;
        res.numUniqueFaces = numUniqueFaces;
        res.numMarked      = numMarked;
        res.oldNumNodes    = numNodes;

        res.d_marks.resize(numElements);
        cudaMemcpy(res.d_marks.data(), marks, numElements * sizeof(uint8_t), cudaMemcpyDeviceToDevice);

        return res;
    }
};

} // namespace amr
} // namespace mars
