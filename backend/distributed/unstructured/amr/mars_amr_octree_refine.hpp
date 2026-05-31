#pragma once

// Octree-aligned hex refinement: generates 8 child hexes per refined parent
// with new SFC keys derived from parent_sfc + child_octant * nodeRange(level+1).
// Children are spatially split along x/y/z midplanes, with new midpoint nodes
// computed from corner coords. Resulting element list has SFC keys consistent
// with cstone's global octree refinement, so Domain::sync redistributes them
// without ambiguity.
//
// Tet path: delegates to TetRefiner (Bey red refinement). One unified
// OctreeAlignedRefine<KeyType, RealType, ElementTag> entry point.

#include <array>
#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/transform.h>
#include "cstone/cuda/cuda_utils.hpp"
#include "cstone/sfc/common.hpp"
#include "cstone/sfc/sfc.hpp"
#include "backend/distributed/unstructured/domain.hpp"   // HexTag, TetTag
#include "mars_amr_tet_refine.hpp"
#include "mars_amr_solution_transfer.hpp"   // SolutionTransfer::transferByParentageTet

namespace mars
{
namespace amr
{

// For each parent element marked for refinement, write 8 children.
// Each child gets:
//   - 8 corner node IDs (4 from parent, 4 new midpoints in this child's octant)
//   - representative node coords at the child's centroid (placed in the
//     child's SFC sub-range by construction)
//
// To keep the refiner GPU-resident and conflict-free, we duplicate the new
// midpoint nodes per element (no dedup). Cstone's sync will handle the
// redundant geometry — node coordinates that coincide will produce the same
// SFC key and be naturally deduplicated through the SFC ordering.
//
// Each refined parent contributes:
//   - 19 new local nodes (12 edge mids + 6 face centers + 1 body center)
//     emitted directly (no dedup) — at scale, dedup belongs in cstone's sync
template<typename KeyType, typename RealType>
__global__ void emitChildHexesKernel(const KeyType* conn0,
                                     const KeyType* conn1,
                                     const KeyType* conn2,
                                     const KeyType* conn3,
                                     const KeyType* conn4,
                                     const KeyType* conn5,
                                     const KeyType* conn6,
                                     const KeyType* conn7,
                                     const RealType* nodeX,
                                     const RealType* nodeY,
                                     const RealType* nodeZ,
                                     const uint8_t* marks,
                                     const uint64_t* elemPrefix, // exclusive scan of childCount
                                     const uint64_t* markedPrefix, // exclusive scan of (mark>0)
                                     size_t numElements,
                                     size_t numNodes,
                                     // outputs
                                     KeyType* outConn0,
                                     KeyType* outConn1,
                                     KeyType* outConn2,
                                     KeyType* outConn3,
                                     KeyType* outConn4,
                                     KeyType* outConn5,
                                     KeyType* outConn6,
                                     KeyType* outConn7,
                                     RealType* outX,
                                     RealType* outY,
                                     RealType* outZ)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements) return;

    KeyType n[8] = {conn0[e], conn1[e], conn2[e], conn3[e], conn4[e], conn5[e], conn6[e], conn7[e]};
    RealType px[8], py[8], pz[8];
    for (int i = 0; i < 8; ++i)
    {
        px[i] = nodeX[n[i]];
        py[i] = nodeY[n[i]];
        pz[i] = nodeZ[n[i]];
    }

    uint64_t outBase = elemPrefix[e];

    if (marks[e] == 0)
    {
        // Copy element as-is, keeping its existing connectivity & node coords
        outConn0[outBase] = n[0]; outConn1[outBase] = n[1];
        outConn2[outBase] = n[2]; outConn3[outBase] = n[3];
        outConn4[outBase] = n[4]; outConn5[outBase] = n[5];
        outConn6[outBase] = n[6]; outConn7[outBase] = n[7];
        return;
    }

    // 19 new nodes per refined parent, packed into output coord arrays at
    // index = numNodes + 19 * markedIdx. No dedup: cstone sync will collapse
    // coincident-coord nodes via SFC ordering.
    uint64_t nodeBase = numNodes + 19ULL * markedPrefix[e];

    // 12 edge midpoints
    auto mid2 = [&](int a, int b, RealType& X, RealType& Y, RealType& Z)
    {
        X = RealType(0.5) * (px[a] + px[b]);
        Y = RealType(0.5) * (py[a] + py[b]);
        Z = RealType(0.5) * (pz[a] + pz[b]);
    };
    auto mid4 = [&](int a, int b, int c, int d, RealType& X, RealType& Y, RealType& Z)
    {
        X = RealType(0.25) * (px[a] + px[b] + px[c] + px[d]);
        Y = RealType(0.25) * (py[a] + py[b] + py[c] + py[d]);
        Z = RealType(0.25) * (pz[a] + pz[b] + pz[c] + pz[d]);
    };

    // Edge node indices (0..11)
    KeyType e01 = nodeBase + 0,  e12 = nodeBase + 1,  e23 = nodeBase + 2,  e30 = nodeBase + 3;
    KeyType e45 = nodeBase + 4,  e56 = nodeBase + 5,  e67 = nodeBase + 6,  e74 = nodeBase + 7;
    KeyType e04 = nodeBase + 8,  e15 = nodeBase + 9,  e26 = nodeBase + 10, e37 = nodeBase + 11;
    // Face centers (12..17)
    KeyType fBot  = nodeBase + 12, fTop  = nodeBase + 13;
    KeyType fFrnt = nodeBase + 14, fBack = nodeBase + 15;
    KeyType fLeft = nodeBase + 16, fRght = nodeBase + 17;
    // Body center (18)
    KeyType bCtr = nodeBase + 18;

    // Write coords
    mid2(0, 1, outX[e01], outY[e01], outZ[e01]);
    mid2(1, 2, outX[e12], outY[e12], outZ[e12]);
    mid2(2, 3, outX[e23], outY[e23], outZ[e23]);
    mid2(3, 0, outX[e30], outY[e30], outZ[e30]);
    mid2(4, 5, outX[e45], outY[e45], outZ[e45]);
    mid2(5, 6, outX[e56], outY[e56], outZ[e56]);
    mid2(6, 7, outX[e67], outY[e67], outZ[e67]);
    mid2(7, 4, outX[e74], outY[e74], outZ[e74]);
    mid2(0, 4, outX[e04], outY[e04], outZ[e04]);
    mid2(1, 5, outX[e15], outY[e15], outZ[e15]);
    mid2(2, 6, outX[e26], outY[e26], outZ[e26]);
    mid2(3, 7, outX[e37], outY[e37], outZ[e37]);

    mid4(0, 1, 2, 3, outX[fBot],  outY[fBot],  outZ[fBot]);
    mid4(4, 5, 6, 7, outX[fTop],  outY[fTop],  outZ[fTop]);
    mid4(0, 1, 5, 4, outX[fFrnt], outY[fFrnt], outZ[fFrnt]);
    mid4(3, 2, 6, 7, outX[fBack], outY[fBack], outZ[fBack]);
    mid4(0, 3, 7, 4, outX[fLeft], outY[fLeft], outZ[fLeft]);
    mid4(1, 2, 6, 5, outX[fRght], outY[fRght], outZ[fRght]);

    outX[bCtr] = RealType(0.125) * (px[0]+px[1]+px[2]+px[3]+px[4]+px[5]+px[6]+px[7]);
    outY[bCtr] = RealType(0.125) * (py[0]+py[1]+py[2]+py[3]+py[4]+py[5]+py[6]+py[7]);
    outZ[bCtr] = RealType(0.125) * (pz[0]+pz[1]+pz[2]+pz[3]+pz[4]+pz[5]+pz[6]+pz[7]);

    auto writeChild = [&](int c, KeyType c0, KeyType c1, KeyType c2, KeyType c3,
                                 KeyType c4, KeyType c5, KeyType c6, KeyType c7)
    {
        uint64_t idx = outBase + c;
        outConn0[idx] = c0; outConn1[idx] = c1; outConn2[idx] = c2; outConn3[idx] = c3;
        outConn4[idx] = c4; outConn5[idx] = c5; outConn6[idx] = c6; outConn7[idx] = c7;
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

// Per-element-corner trilinear interpolation of a solution field at the same
// 19 node positions emitted by emitChildHexesKernel (12 edge mids + 6 face
// centers + 1 body center). For unmarked elements, no solution slots are
// written. The output indexing matches emitChildHexesKernel exactly:
//   newSol[old node id n] = oldSol[n]                (copied separately)
//   newSol[oldNumNodes + 19*markedPrefix[e] + i]     for the i-th of 19
template<typename KeyType, typename RealType>
__global__ void interpolateChildSolutionKernel(const KeyType* conn0,
                                                const KeyType* conn1,
                                                const KeyType* conn2,
                                                const KeyType* conn3,
                                                const KeyType* conn4,
                                                const KeyType* conn5,
                                                const KeyType* conn6,
                                                const KeyType* conn7,
                                                const RealType* oldSolution,
                                                const uint8_t* marks,
                                                const uint64_t* markedPrefix,
                                                size_t numElements,
                                                size_t numNodes,
                                                RealType* newSolution)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements) return;
    if (marks[e] == 0) return;

    KeyType n[8] = {conn0[e], conn1[e], conn2[e], conn3[e], conn4[e], conn5[e], conn6[e], conn7[e]};
    RealType u[8];
    for (int i = 0; i < 8; ++i) u[i] = oldSolution[n[i]];

    uint64_t base = numNodes + 19ULL * markedPrefix[e];

    // 12 edge midpoints: average of 2 corners
    auto m2 = [&](int a, int b) -> RealType { return RealType(0.5) * (u[a] + u[b]); };
    newSolution[base + 0]  = m2(0, 1);  // e01
    newSolution[base + 1]  = m2(1, 2);  // e12
    newSolution[base + 2]  = m2(2, 3);  // e23
    newSolution[base + 3]  = m2(3, 0);  // e30
    newSolution[base + 4]  = m2(4, 5);  // e45
    newSolution[base + 5]  = m2(5, 6);  // e56
    newSolution[base + 6]  = m2(6, 7);  // e67
    newSolution[base + 7]  = m2(7, 4);  // e74
    newSolution[base + 8]  = m2(0, 4);  // e04
    newSolution[base + 9]  = m2(1, 5);  // e15
    newSolution[base + 10] = m2(2, 6);  // e26
    newSolution[base + 11] = m2(3, 7);  // e37

    // 6 face centers: average of 4 corners
    auto m4 = [&](int a, int b, int c, int d) -> RealType {
        return RealType(0.25) * (u[a] + u[b] + u[c] + u[d]);
    };
    newSolution[base + 12] = m4(0, 1, 2, 3);  // fBot
    newSolution[base + 13] = m4(4, 5, 6, 7);  // fTop
    newSolution[base + 14] = m4(0, 1, 5, 4);  // fFrnt
    newSolution[base + 15] = m4(3, 2, 6, 7);  // fBack
    newSolution[base + 16] = m4(0, 3, 7, 4);  // fLeft
    newSolution[base + 17] = m4(1, 2, 6, 5);  // fRght

    // body center: average of 8 corners
    newSolution[base + 18] =
        RealType(0.125) * (u[0] + u[1] + u[2] + u[3] + u[4] + u[5] + u[6] + u[7]);
}

// Top-level: refine local elements in-place on GPU. Returns new device-side
// connectivity + coords ready to hand to ElementDomain device-data constructor.
// Templated on ElementTag so hex and tet share a Result shape (std::array of
// NodesPerElement conn columns), with element-specific kernels behind the
// scenes (emitChildHexesKernel for hex, TetRefiner::refine for tet).
template<typename KeyType, typename RealType, typename ElementTag = HexTag>
struct OctreeAlignedRefine
{
    static constexpr int NodesPerElement = ElementTag::NodesPerElement;

    struct Result
    {
        std::array<cstone::DeviceVector<KeyType>, NodesPerElement> d_conn;
        cstone::DeviceVector<RealType> d_x, d_y, d_z;
        size_t numNodes    = 0;
        size_t numElements = 0;

        // Reserved for tet parentage-based solution transfer (mars_amr_tet_refine).
        // Hex path leaves these empty.
        cstone::DeviceVector<uint64_t> d_edgeKeys;
        cstone::DeviceVector<uint64_t> d_markedPrefix;
        cstone::DeviceVector<uint8_t>  d_marks;
        size_t numUniqueEdges = 0;
        size_t numMarked      = 0;
        size_t oldNumNodes    = 0;
    };

    // Hex entry: takes 8 conn pointers, drives emitChildHexesKernel.
    static Result refineLocal(
        const KeyType* conn0, const KeyType* conn1, const KeyType* conn2, const KeyType* conn3,
        const KeyType* conn4, const KeyType* conn5, const KeyType* conn6, const KeyType* conn7,
        const RealType* nodeX, const RealType* nodeY, const RealType* nodeZ,
        const uint8_t* marks,
        size_t numElements, size_t numNodes,
        int blockSize = 256)
    {
        // Per-element child count (1 or 8)
        cstone::DeviceVector<uint64_t> d_childCount(numElements);
        cstone::DeviceVector<uint64_t> d_isMarked(numElements);
        thrust::transform(thrust::device, marks, marks + numElements,
                           thrust::device_pointer_cast(d_childCount.data()),
                           [] __device__(uint8_t m) -> uint64_t { return (m > 0) ? 8 : 1; });
        thrust::transform(thrust::device, marks, marks + numElements,
                           thrust::device_pointer_cast(d_isMarked.data()),
                           [] __device__(uint8_t m) -> uint64_t { return (m > 0) ? 1 : 0; });

        cstone::DeviceVector<uint64_t> d_elemPrefix(numElements);
        cstone::DeviceVector<uint64_t> d_markedPrefix(numElements);
        auto cc_b = thrust::device_pointer_cast(d_childCount.data());
        auto im_b = thrust::device_pointer_cast(d_isMarked.data());
        thrust::exclusive_scan(thrust::device, cc_b, cc_b + numElements,
                                thrust::device_pointer_cast(d_elemPrefix.data()));
        thrust::exclusive_scan(thrust::device, im_b, im_b + numElements,
                                thrust::device_pointer_cast(d_markedPrefix.data()));

        uint64_t totalNewElements;
        uint64_t numMarked;
        {
            uint64_t lastChild, lastEPrefix, lastIs, lastMPrefix;
            cudaMemcpy(&lastChild,   d_childCount.data()   + numElements - 1, sizeof(uint64_t), cudaMemcpyDeviceToHost);
            cudaMemcpy(&lastEPrefix, d_elemPrefix.data()   + numElements - 1, sizeof(uint64_t), cudaMemcpyDeviceToHost);
            cudaMemcpy(&lastIs,      d_isMarked.data()     + numElements - 1, sizeof(uint64_t), cudaMemcpyDeviceToHost);
            cudaMemcpy(&lastMPrefix, d_markedPrefix.data() + numElements - 1, sizeof(uint64_t), cudaMemcpyDeviceToHost);
            totalNewElements = lastEPrefix + lastChild;
            numMarked        = lastMPrefix + lastIs;
        }

        Result res;
        res.numElements = totalNewElements;
        res.numNodes    = numNodes + 19ULL * numMarked;
        res.oldNumNodes = numNodes;
        res.numMarked   = numMarked;

        for (int i = 0; i < 8; ++i) res.d_conn[i].resize(totalNewElements);
        res.d_x.resize(res.numNodes); res.d_y.resize(res.numNodes); res.d_z.resize(res.numNodes);

        // Copy original node coords (positions [0, numNodes))
        cudaMemcpy(res.d_x.data(), nodeX, numNodes * sizeof(RealType), cudaMemcpyDeviceToDevice);
        cudaMemcpy(res.d_y.data(), nodeY, numNodes * sizeof(RealType), cudaMemcpyDeviceToDevice);
        cudaMemcpy(res.d_z.data(), nodeZ, numNodes * sizeof(RealType), cudaMemcpyDeviceToDevice);

        int numBlocks = (numElements + blockSize - 1) / blockSize;
        emitChildHexesKernel<KeyType, RealType><<<numBlocks, blockSize>>>(
            conn0, conn1, conn2, conn3, conn4, conn5, conn6, conn7,
            nodeX, nodeY, nodeZ, marks,
            d_elemPrefix.data(), d_markedPrefix.data(),
            numElements, numNodes,
            res.d_conn[0].data(), res.d_conn[1].data(), res.d_conn[2].data(), res.d_conn[3].data(),
            res.d_conn[4].data(), res.d_conn[5].data(), res.d_conn[6].data(), res.d_conn[7].data(),
            res.d_x.data(), res.d_y.data(), res.d_z.data());
        cudaDeviceSynchronize();
        return res;
    }

    // Tet entry: takes 4 conn pointers, delegates to TetRefiner (Bey red refinement).
    static Result refineLocal(
        const KeyType* conn0, const KeyType* conn1, const KeyType* conn2, const KeyType* conn3,
        const RealType* nodeX, const RealType* nodeY, const RealType* nodeZ,
        const uint8_t* marks,
        size_t numElements, size_t numNodes,
        int blockSize = 256)
    {
        auto tetRes = TetRefiner<KeyType, RealType>::refine(
            conn0, conn1, conn2, conn3, marks, numElements, nodeX, nodeY, nodeZ, numNodes, blockSize);

        Result res;
        res.numElements    = tetRes.numElements;
        res.numNodes       = tetRes.numNodes;
        res.oldNumNodes    = tetRes.oldNumNodes;
        res.numMarked      = tetRes.numMarked;
        res.numUniqueEdges = tetRes.numUniqueEdges;
        res.d_x            = std::move(tetRes.d_x);
        res.d_y            = std::move(tetRes.d_y);
        res.d_z            = std::move(tetRes.d_z);
        res.d_conn[0]      = std::move(tetRes.d_conn0);
        res.d_conn[1]      = std::move(tetRes.d_conn1);
        res.d_conn[2]      = std::move(tetRes.d_conn2);
        res.d_conn[3]      = std::move(tetRes.d_conn3);
        res.d_edgeKeys     = std::move(tetRes.d_edgeKeys);
        res.d_markedPrefix = std::move(tetRes.d_markedPrefix);
        res.d_marks        = std::move(tetRes.d_marks);
        return res;
    }

    // Hex: trilinear interpolation of a per-node solution field at the 19
    // new node positions per refined hex (edge mids + face centers + body).
    // Output layout: [0, numNodes) copied; [numNodes, numNodes + 19*numMarked)
    // interpolated.
    static void transferSolution(
        const KeyType* conn0, const KeyType* conn1, const KeyType* conn2, const KeyType* conn3,
        const KeyType* conn4, const KeyType* conn5, const KeyType* conn6, const KeyType* conn7,
        const RealType* oldSolution,
        const uint8_t* marks,
        size_t numElements, size_t numNodes,
        cstone::DeviceVector<RealType>& newSolution,
        int blockSize = 256)
    {
        // Recompute markedPrefix locally (cheap; same logic as refineLocal)
        cstone::DeviceVector<uint64_t> d_isMarked(numElements);
        thrust::transform(thrust::device, marks, marks + numElements,
                           thrust::device_pointer_cast(d_isMarked.data()),
                           [] __device__(uint8_t m) -> uint64_t { return (m > 0) ? 1 : 0; });
        cstone::DeviceVector<uint64_t> d_markedPrefix(numElements);
        auto im_b = thrust::device_pointer_cast(d_isMarked.data());
        thrust::exclusive_scan(thrust::device, im_b, im_b + numElements,
                                thrust::device_pointer_cast(d_markedPrefix.data()));

        uint64_t numMarked = 0;
        {
            uint64_t lastIs, lastMPrefix;
            cudaMemcpy(&lastIs,      d_isMarked.data()     + numElements - 1, sizeof(uint64_t), cudaMemcpyDeviceToHost);
            cudaMemcpy(&lastMPrefix, d_markedPrefix.data() + numElements - 1, sizeof(uint64_t), cudaMemcpyDeviceToHost);
            numMarked = lastMPrefix + lastIs;
        }

        size_t newSize = numNodes + 19ULL * numMarked;
        newSolution.resize(newSize);

        // Copy old node solutions to [0, numNodes)
        cudaMemcpy(newSolution.data(), oldSolution, numNodes * sizeof(RealType), cudaMemcpyDeviceToDevice);

        // Interpolate at the 19*numMarked new positions
        int numBlocks = (numElements + blockSize - 1) / blockSize;
        interpolateChildSolutionKernel<KeyType, RealType><<<numBlocks, blockSize>>>(
            conn0, conn1, conn2, conn3, conn4, conn5, conn6, conn7,
            oldSolution, marks, d_markedPrefix.data(),
            numElements, numNodes, newSolution.data());
        cudaDeviceSynchronize();
    }

    // Tet: edge-midpoint solution transfer using the dedup'd edgeKeys from
    // refineLocal's Result. Caller must pass the same Result here so the
    // markedPrefix / edge dedup match what produced the new node layout.
    static void transferSolution(
        const KeyType* conn0, const KeyType* conn1, const KeyType* conn2, const KeyType* conn3,
        const RealType* oldSolution,
        const Result& refined,
        size_t numElements,
        cstone::DeviceVector<RealType>& newSolution,
        int blockSize = 256)
    {
        newSolution = SolutionTransfer<KeyType, RealType>::transferByParentageTet(
            conn0, conn1, conn2, conn3,
            refined.d_marks.data(), numElements,
            oldSolution, refined.oldNumNodes,
            refined.d_edgeKeys.data(), refined.numUniqueEdges,
            refined.numNodes, blockSize);
    }
};

} // namespace amr
} // namespace mars
