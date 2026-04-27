#pragma once

#include <thrust/device_vector.h>
#include <thrust/for_each.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/sort.h>
#include <thrust/binary_search.h>
#include <thrust/transform.h>
#include <cmath>

#include "cstone/cuda/cuda_utils.hpp"
#include "cstone/sfc/sfc.hpp"

namespace mars
{
namespace amr
{

template<typename KeyType, typename RealType>
__global__ void transferSfcKernel(const KeyType* oldSfcKeys,
                                  const RealType* oldX,
                                  const RealType* oldY,
                                  const RealType* oldZ,
                                  const RealType* oldNodeSolution,
                                  size_t oldNumNodes,
                                  const RealType* newX,
                                  const RealType* newY,
                                  const RealType* newZ,
                                  RealType* newSolution,
                                  size_t newNumNodes,
                                  RealType bxMin,
                                  RealType bxMax,
                                  RealType byMin,
                                  RealType byMax,
                                  RealType bzMin,
                                  RealType bzMax)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= newNumNodes) return;

    RealType px = newX[i];
    RealType py = newY[i];
    RealType pz = newZ[i];

    constexpr unsigned maxCoord = (1u << 10) - 1;
    unsigned ix = 0, iy = 0, iz = 0;
    RealType rangeX = bxMax - bxMin;
    RealType rangeY = byMax - byMin;
    RealType rangeZ = bzMax - bzMin;

    if (rangeX > RealType(0))
        ix = min(unsigned((px - bxMin) / rangeX * maxCoord), maxCoord);
    if (rangeY > RealType(0))
        iy = min(unsigned((py - byMin) / rangeY * maxCoord), maxCoord);
    if (rangeZ > RealType(0))
        iz = min(unsigned((pz - bzMin) / rangeZ * maxCoord), maxCoord);

    KeyType key = 0;
    for (int b = 0; b < 10; ++b)
    {
        key |= (KeyType((ix >> b) & 1) << (3 * b));
        key |= (KeyType((iy >> b) & 1) << (3 * b + 1));
        key |= (KeyType((iz >> b) & 1) << (3 * b + 2));
    }

    size_t lo = 0, hi = oldNumNodes;
    while (lo < hi)
    {
        size_t mid = (lo + hi) / 2;
        if (oldSfcKeys[mid] < key)
            lo = mid + 1;
        else
            hi = mid;
    }

    // 32 is enough because SFC preserves spatial locality: nearby nodes
    // in 3D map to nearby positions in the sorted key array
    constexpr int NEIGHBORHOOD = 32;
    size_t searchStart = (lo > NEIGHBORHOOD) ? lo - NEIGHBORHOOD : 0;
    size_t searchEnd   = min(lo + NEIGHBORHOOD, oldNumNodes);

    RealType bestDist = RealType(1e30);
    RealType bestVal  = RealType(0);

    for (size_t j = searchStart; j < searchEnd; ++j)
    {
        RealType dx   = px - oldX[j];
        RealType dy   = py - oldY[j];
        RealType dz   = pz - oldZ[j];
        RealType dist = dx * dx + dy * dy + dz * dz;

        if (dist < bestDist)
        {
            bestDist = dist;
            bestVal  = oldNodeSolution[j];
        }
    }

    newSolution[i] = bestVal;
}

// Interpolation weights from hex refinement parentage:
//   edge midpoint = avg of 2 endpoints
//   face center   = avg of 4 face corners
//   body center   = avg of 8 element corners
template<typename KeyType, typename RealType>
__global__ void transferByParentageKernel(const KeyType* oldConn0,
                                          const KeyType* oldConn1,
                                          const KeyType* oldConn2,
                                          const KeyType* oldConn3,
                                          const KeyType* oldConn4,
                                          const KeyType* oldConn5,
                                          const KeyType* oldConn6,
                                          const KeyType* oldConn7,
                                          const uint8_t* marks,
                                          const uint64_t* markedPrefix,
                                          const RealType* oldNodeSolution,
                                          RealType* newNodeSolution,
                                          KeyType edgeBaseNode,
                                          KeyType faceBaseNode,
                                          KeyType bodyBaseNode,
                                          const uint64_t* sortedEdgeKeys,
                                          size_t numUniqueEdges,
                                          const uint64_t* sortedFaceKeysLo,
                                          const uint64_t* sortedFaceKeysHi,
                                          size_t numUniqueFaces,
                                          size_t numElements,
                                          size_t oldNumNodes)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements || marks[e] <= 0) return;

    KeyType n[8] = {oldConn0[e], oldConn1[e], oldConn2[e], oldConn3[e],
                    oldConn4[e], oldConn5[e], oldConn6[e], oldConn7[e]};

    RealType u[8];
    for (int i = 0; i < 8; ++i)
        u[i] = (n[i] < oldNumNodes) ? oldNodeSolution[n[i]] : RealType(0);

    // Must match the key encoding used by HexRefiner
    auto findEdge = [&](KeyType a, KeyType b) -> KeyType
    {
        uint64_t lo_val = (a < b) ? a : b;
        uint64_t hi_val = (a < b) ? b : a;
        uint64_t key    = (hi_val << 32) | (lo_val & 0xFFFFFFFF);
        size_t lo = 0, hi = numUniqueEdges;
        while (lo < hi)
        {
            size_t mid = (lo + hi) / 2;
            if (sortedEdgeKeys[mid] < key) lo = mid + 1;
            else hi = mid;
        }
        return edgeBaseNode + lo;
    };

    auto findFace = [&](KeyType a, KeyType b, KeyType c, KeyType d) -> KeyType
    {
        uint64_t sa = a, sb = b, sc = c, sd = d;
        if (sa > sb) { uint64_t t = sa; sa = sb; sb = t; }
        if (sc > sd) { uint64_t t = sc; sc = sd; sd = t; }
        if (sa > sc) { uint64_t t = sa; sa = sc; sc = t; }
        if (sb > sd) { uint64_t t = sb; sb = sd; sd = t; }
        if (sb > sc) { uint64_t t = sb; sb = sc; sc = t; }
        uint64_t kLo = (sa << 32) | (sb & 0xFFFFFFFF);
        uint64_t kHi = (sc << 32) | (sd & 0xFFFFFFFF);
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

    struct EdgeDef { int a, b; };
    EdgeDef edges[12] = {{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7}};
    for (int ei = 0; ei < 12; ++ei)
    {
        KeyType midNode              = findEdge(n[edges[ei].a], n[edges[ei].b]);
        newNodeSolution[midNode] = RealType(0.5) * (u[edges[ei].a] + u[edges[ei].b]);
    }

    struct FaceDef { int a, b, c, d; };
    FaceDef faces[6] = {{0,1,2,3},{4,5,6,7},{0,1,5,4},{3,2,6,7},{0,3,7,4},{1,2,6,5}};
    for (int fi = 0; fi < 6; ++fi)
    {
        KeyType fNode             = findFace(n[faces[fi].a], n[faces[fi].b], n[faces[fi].c], n[faces[fi].d]);
        newNodeSolution[fNode] = RealType(0.25) * (u[faces[fi].a] + u[faces[fi].b] + u[faces[fi].c] + u[faces[fi].d]);
    }

    KeyType bNode              = bodyBaseNode + markedPrefix[e];
    newNodeSolution[bNode] = RealType(0.125) * (u[0] + u[1] + u[2] + u[3] + u[4] + u[5] + u[6] + u[7]);
}

// Tet4 parentage transfer: only edge midpoints (no face/body centers)
// Interpolation: edge midpoint = avg of 2 endpoints
template<typename KeyType, typename RealType>
__global__ void transferByParentageTetKernel(const KeyType* oldConn0,
                                             const KeyType* oldConn1,
                                             const KeyType* oldConn2,
                                             const KeyType* oldConn3,
                                             const uint8_t* marks,
                                             const RealType* oldNodeSolution,
                                             RealType* newNodeSolution,
                                             KeyType edgeBaseNode,
                                             const uint64_t* sortedEdgeKeys,
                                             size_t numUniqueEdges,
                                             size_t numElements,
                                             size_t oldNumNodes)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements || marks[e] <= 0) return;

    KeyType n[4] = {oldConn0[e], oldConn1[e], oldConn2[e], oldConn3[e]};
    RealType u[4];
    for (int i = 0; i < 4; ++i)
        u[i] = (n[i] < oldNumNodes) ? oldNodeSolution[n[i]] : RealType(0);

    auto findEdge = [&](KeyType a, KeyType b) -> KeyType
    {
        uint64_t lo_val = (a < b) ? a : b;
        uint64_t hi_val = (a < b) ? b : a;
        uint64_t key    = (hi_val << 32) | (lo_val & 0xFFFFFFFF);
        size_t lo = 0, hi = numUniqueEdges;
        while (lo < hi)
        {
            size_t mid = (lo + hi) / 2;
            if (sortedEdgeKeys[mid] < key) lo = mid + 1;
            else hi = mid;
        }
        return edgeBaseNode + lo;
    };

    // 6 edges of tet4
    struct EdgeDef { int a, b; };
    EdgeDef edges[6] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
    for (int ei = 0; ei < 6; ++ei)
    {
        KeyType midNode              = findEdge(n[edges[ei].a], n[edges[ei].b]);
        newNodeSolution[midNode] = RealType(0.5) * (u[edges[ei].a] + u[edges[ei].b]);
    }
}

template<typename RealType>
__global__ void copyOldNodeSolutionsKernel(const RealType* oldNodeSolution,
                                           RealType* newNodeSolution,
                                           size_t oldNumNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= oldNumNodes) return;
    newNodeSolution[i] = oldNodeSolution[i];
}

template<typename KeyType, typename RealType>
__global__ void transferOctreeKernel(const KeyType* treeLeaves,
                                     size_t numLeaves,
                                     const size_t* layout,
                                     const RealType* oldX,
                                     const RealType* oldY,
                                     const RealType* oldZ,
                                     const RealType* oldNodeSolution,
                                     size_t oldNumNodes,
                                     const RealType* newX,
                                     const RealType* newY,
                                     const RealType* newZ,
                                     RealType* newSolution,
                                     size_t newNumNodes,
                                     RealType bxMin,
                                     RealType bxMax,
                                     RealType byMin,
                                     RealType byMax,
                                     RealType bzMin,
                                     RealType bzMax)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= newNumNodes) return;

    RealType px = newX[i];
    RealType py = newY[i];
    RealType pz = newZ[i];

    constexpr unsigned maxCoord = (1u << 10) - 1;
    unsigned ix = 0, iy = 0, iz = 0;
    RealType rangeX = bxMax - bxMin;
    RealType rangeY = byMax - byMin;
    RealType rangeZ = bzMax - bzMin;

    if (rangeX > RealType(0))
        ix = min(unsigned((px - bxMin) / rangeX * maxCoord), maxCoord);
    if (rangeY > RealType(0))
        iy = min(unsigned((py - byMin) / rangeY * maxCoord), maxCoord);
    if (rangeZ > RealType(0))
        iz = min(unsigned((pz - bzMin) / rangeZ * maxCoord), maxCoord);

    KeyType key = 0;
    for (int b = 0; b < 10; ++b)
    {
        key |= (KeyType((ix >> b) & 1) << (3 * b));
        key |= (KeyType((iy >> b) & 1) << (3 * b + 1));
        key |= (KeyType((iz >> b) & 1) << (3 * b + 2));
    }

    size_t lo = 0, hi = numLeaves;
    while (lo < hi)
    {
        size_t mid = (lo + hi) / 2;
        if (treeLeaves[mid + 1] <= key)
            lo = mid + 1;
        else
            hi = mid;
    }

    // ±1 leaf: a point near a leaf boundary may be closest to a node in the neighbor leaf
    size_t leafStart = (lo > 0) ? lo - 1 : 0;
    size_t leafEnd   = min(lo + 2, numLeaves);

    size_t searchStart = layout[leafStart];
    size_t searchEnd   = layout[leafEnd];
    searchEnd          = min(searchEnd, oldNumNodes);

    RealType bestDist = RealType(1e30);
    RealType bestVal  = RealType(0);

    for (size_t j = searchStart; j < searchEnd; ++j)
    {
        RealType dx   = px - oldX[j];
        RealType dy   = py - oldY[j];
        RealType dz   = pz - oldZ[j];
        RealType dist = dx * dx + dy * dy + dz * dz;

        if (dist < bestDist)
        {
            bestDist = dist;
            bestVal  = oldNodeSolution[j];
        }
    }

    newSolution[i] = bestVal;
}

enum class TransferMethod
{
    SfcLookup,
    OctreeTraversal,
    Parentage
};

template<typename KeyType, typename RealType>
class SolutionTransfer
{
public:
    static cstone::DeviceVector<RealType> transferByParentage(
        const KeyType* oldConn0,
        const KeyType* oldConn1,
        const KeyType* oldConn2,
        const KeyType* oldConn3,
        const KeyType* oldConn4,
        const KeyType* oldConn5,
        const KeyType* oldConn6,
        const KeyType* oldConn7,
        const uint8_t* marks,
        const uint64_t* markedPrefix,
        size_t numElements,
        const RealType* oldNodeSolution,
        size_t oldNumNodes,
        const uint64_t* sortedEdgeKeys,
        size_t numUniqueEdges,
        const uint64_t* sortedFaceKeysLo,
        const uint64_t* sortedFaceKeysHi,
        size_t numUniqueFaces,
        size_t numMarked,
        size_t newNumNodes,
        int blockSize = 256)
    {
        KeyType edgeBaseNode = oldNumNodes;
        KeyType faceBaseNode = oldNumNodes + numUniqueEdges;
        KeyType bodyBaseNode = oldNumNodes + numUniqueEdges + numUniqueFaces;

        cstone::DeviceVector<RealType> d_newSolution(newNumNodes, RealType(0));

        int numBlocks1 = (oldNumNodes + blockSize - 1) / blockSize;
        copyOldNodeSolutionsKernel<<<numBlocks1, blockSize>>>(oldNodeSolution, d_newSolution.data(), oldNumNodes);

        int numBlocks2 = (numElements + blockSize - 1) / blockSize;
        transferByParentageKernel<<<numBlocks2, blockSize>>>(
            oldConn0, oldConn1, oldConn2, oldConn3, oldConn4, oldConn5, oldConn6, oldConn7, marks, markedPrefix,
            oldNodeSolution, d_newSolution.data(), edgeBaseNode, faceBaseNode, bodyBaseNode, sortedEdgeKeys,
            numUniqueEdges, sortedFaceKeysLo, sortedFaceKeysHi, numUniqueFaces, numElements, oldNumNodes);

        cudaDeviceSynchronize();
        return d_newSolution;
    }

    // Tet4 parentage transfer: edge midpoints only
    static cstone::DeviceVector<RealType> transferByParentageTet(
        const KeyType* oldConn0, const KeyType* oldConn1,
        const KeyType* oldConn2, const KeyType* oldConn3,
        const uint8_t* marks, size_t numElements,
        const RealType* oldNodeSolution, size_t oldNumNodes,
        const uint64_t* sortedEdgeKeys, size_t numUniqueEdges,
        size_t newNumNodes, int blockSize = 256)
    {
        cstone::DeviceVector<RealType> d_newSolution(newNumNodes, RealType(0));

        int numBlocks1 = (oldNumNodes + blockSize - 1) / blockSize;
        copyOldNodeSolutionsKernel<<<numBlocks1, blockSize>>>(oldNodeSolution, d_newSolution.data(), oldNumNodes);

        KeyType edgeBaseNode = oldNumNodes;
        int numBlocks2       = (numElements + blockSize - 1) / blockSize;
        transferByParentageTetKernel<<<numBlocks2, blockSize>>>(
            oldConn0, oldConn1, oldConn2, oldConn3, marks,
            oldNodeSolution, d_newSolution.data(), edgeBaseNode,
            sortedEdgeKeys, numUniqueEdges, numElements, oldNumNodes);

        cudaDeviceSynchronize();
        return d_newSolution;
    }

    static void transferSfc(const cstone::DeviceVector<KeyType>& d_oldSfcKeys,
                            const cstone::DeviceVector<RealType>& d_oldX,
                            const cstone::DeviceVector<RealType>& d_oldY,
                            const cstone::DeviceVector<RealType>& d_oldZ,
                            const RealType* d_oldNodeSolution,
                            size_t oldNumNodes,
                            const cstone::DeviceVector<RealType>& d_newX,
                            const cstone::DeviceVector<RealType>& d_newY,
                            const cstone::DeviceVector<RealType>& d_newZ,
                            cstone::DeviceVector<RealType>& d_newSolution,
                            const cstone::Box<RealType>& box,
                            int blockSize = 256)
    {
        size_t newNumNodes = d_newX.size();
        d_newSolution.resize(newNumNodes);

        int numBlocks = (newNumNodes + blockSize - 1) / blockSize;
        transferSfcKernel<KeyType, RealType><<<numBlocks, blockSize>>>(
            d_oldSfcKeys.data(), d_oldX.data(), d_oldY.data(), d_oldZ.data(), d_oldNodeSolution, oldNumNodes,
            d_newX.data(), d_newY.data(), d_newZ.data(), d_newSolution.data(), newNumNodes,
            box.xmin(), box.xmax(), box.ymin(), box.ymax(), box.zmin(), box.zmax());
        cudaDeviceSynchronize();
    }

    // More robust than SFC lookup for non-uniform meshes where
    // the fixed neighborhood heuristic may miss distant nodes
    template<typename AcceleratorTag>
    static void transferOctree(const ElementDomain<HexTag, RealType, KeyType, AcceleratorTag>& oldDomain,
                               const RealType* d_oldNodeSolution,
                               const cstone::DeviceVector<RealType>& d_newX,
                               const cstone::DeviceVector<RealType>& d_newY,
                               const cstone::DeviceVector<RealType>& d_newZ,
                               cstone::DeviceVector<RealType>& d_newSolution,
                               int blockSize = 256)
    {
        size_t oldNumNodes = oldDomain.getNodeCount();
        size_t newNumNodes = d_newX.size();
        d_newSolution.resize(newNumNodes);

        const auto& cstoneDomain = oldDomain.getDomain();
        auto treeLeaves          = cstoneDomain.focusTree().treeLeavesAcc();
        auto layout              = cstoneDomain.layout();
        size_t numLeaves         = treeLeaves.size() - 1;

        const auto& d_x = oldDomain.getNodeX();
        const auto& d_y = oldDomain.getNodeY();
        const auto& d_z = oldDomain.getNodeZ();
        const auto& box = oldDomain.getBoundingBox();

        cstone::DeviceVector<size_t> d_layout(layout.begin(), layout.end());

        int numBlocks = (newNumNodes + blockSize - 1) / blockSize;
        transferOctreeKernel<KeyType, RealType><<<numBlocks, blockSize>>>(
            treeLeaves.data(), numLeaves, d_layout.data(), d_x.data(), d_y.data(), d_z.data(), d_oldNodeSolution,
            oldNumNodes, d_newX.data(), d_newY.data(), d_newZ.data(), d_newSolution.data(), newNumNodes, box.xmin(),
            box.xmax(), box.ymin(), box.ymax(), box.zmin(), box.zmax());
        cudaDeviceSynchronize();
    }
};

} // namespace amr
} // namespace mars
