#pragma once

#include <vector>
#include <algorithm>
#include <limits>
#include <cstdint>
#include <cstdio>
#include <cuda_runtime.h>

namespace mars {
namespace fem {

// =============================================================================
// Element graph coloring for atomic-free CVFEM assembly.
//
// Two hex elements "conflict" if they share at least one node.
// A valid coloring assigns colors such that no two same-color elements conflict.
// Within a single-color kernel launch, elements never share nodes, so CSR
// scatter can use direct += instead of atomicAdd.
//
// For a structured hex mesh: exactly 8 colors suffice (parity of i,j,k indices).
// For unstructured hex meshes: typically 8-12 colors via greedy coloring.
// =============================================================================
struct CvfemColoringData {
    int  numColors  = 0;
    int* d_colorPerm    = nullptr; // device: element permutation grouped by color
    std::vector<int> h_colorOffsets;  // host: [0, c0_end, c1_end, ..., N]

    CvfemColoringData() = default;
    CvfemColoringData(const CvfemColoringData&) = delete;
    CvfemColoringData& operator=(const CvfemColoringData&) = delete;

    ~CvfemColoringData() {
        if (d_colorPerm) { cudaFree(d_colorPerm); d_colorPerm = nullptr; }
    }

    int colorSize(int c) const {
        return h_colorOffsets[c + 1] - h_colorOffsets[c];
    }
    const int* colorStart(int c) const { return d_colorPerm + h_colorOffsets[c]; }
};

// =============================================================================
// Host-side greedy graph coloring.
// Connectivity arrays must be host pointers, each of length numElements.
// numNodes = total node count (for sizing the reverse map).
// =============================================================================
template<typename KeyType>
void buildCvfemColoring(
    const KeyType* h_conn0, const KeyType* h_conn1,
    const KeyType* h_conn2, const KeyType* h_conn3,
    const KeyType* h_conn4, const KeyType* h_conn5,
    const KeyType* h_conn6, const KeyType* h_conn7,
    size_t numElements,
    size_t numNodes,
    CvfemColoringData& out)
{
    const KeyType* conn[8] = {
        h_conn0, h_conn1, h_conn2, h_conn3,
        h_conn4, h_conn5, h_conn6, h_conn7
    };

    // ------------------------------------------------------------------
    // Step 1: Build node→elements reverse map (CSR format)
    // ------------------------------------------------------------------
    std::vector<int> nodeElemCount(numNodes, 0);
    for (size_t e = 0; e < numElements; ++e)
        for (int n = 0; n < 8; ++n)
            ++nodeElemCount[static_cast<size_t>(conn[n][e])];

    std::vector<int> nodeElemOff(numNodes + 1, 0);
    for (size_t v = 0; v < numNodes; ++v)
        nodeElemOff[v + 1] = nodeElemOff[v] + nodeElemCount[v];

    std::vector<int> nodeElem(nodeElemOff[numNodes]);
    std::fill(nodeElemCount.begin(), nodeElemCount.end(), 0);
    for (size_t e = 0; e < numElements; ++e)
        for (int n = 0; n < 8; ++n) {
            size_t v = static_cast<size_t>(conn[n][e]);
            nodeElem[nodeElemOff[v] + nodeElemCount[v]++] = static_cast<int>(e);
        }

    // ------------------------------------------------------------------
    // Step 2: Greedy coloring using tick-based marking
    // Each element e gets the smallest color not used by any already-colored
    // element that shares a node with e.
    // ------------------------------------------------------------------
    std::vector<int> elemColor(numElements, -1);
    std::vector<int> usedBuf(32, 0);  // grows if needed (>8 colors unlikely)
    int numColors = 0;
    int tick = 1;  // never reset — 2^31 > any practical element count

    for (size_t e = 0; e < numElements; ++e) {
        int maxNbColor = -1;

        // Mark colors of all neighbor elements via shared nodes
        for (int n = 0; n < 8; ++n) {
            size_t v = static_cast<size_t>(conn[n][e]);
            for (int k = nodeElemOff[v]; k < nodeElemOff[v + 1]; ++k) {
                int nb = nodeElem[k];
                if (nb != static_cast<int>(e) && elemColor[nb] >= 0) {
                    int c = elemColor[nb];
                    if (c >= static_cast<int>(usedBuf.size()))
                        usedBuf.resize(c + 2, 0);
                    usedBuf[c] = tick;
                    if (c > maxNbColor) maxNbColor = c;
                }
            }
        }

        // Find smallest unused color
        int c = 0;
        while (c <= maxNbColor &&
               c < static_cast<int>(usedBuf.size()) &&
               usedBuf[c] == tick) ++c;

        elemColor[e] = c;
        if (c >= numColors) numColors = c + 1;
        if (c >= static_cast<int>(usedBuf.size()))
            usedBuf.resize(c + 2, 0);
        ++tick;
    }

    // ------------------------------------------------------------------
    // Step 3: Build color permutation and host offsets
    // ------------------------------------------------------------------
    std::vector<int> colorCount(numColors, 0);
    for (size_t e = 0; e < numElements; ++e)
        ++colorCount[elemColor[e]];

    std::vector<int> h_offsets(numColors + 1, 0);
    for (int c = 0; c < numColors; ++c)
        h_offsets[c + 1] = h_offsets[c] + colorCount[c];

    std::vector<int> h_perm(numElements);
    std::fill(colorCount.begin(), colorCount.end(), 0);
    for (size_t e = 0; e < numElements; ++e) {
        int c = elemColor[e];
        h_perm[h_offsets[c] + colorCount[c]++] = static_cast<int>(e);
    }

    // ------------------------------------------------------------------
    // Step 4: Upload permutation to device
    // ------------------------------------------------------------------
    if (out.d_colorPerm) { cudaFree(out.d_colorPerm); out.d_colorPerm = nullptr; }

    cudaMalloc(&out.d_colorPerm, numElements * sizeof(int));
    cudaMemcpy(out.d_colorPerm, h_perm.data(),
               numElements * sizeof(int), cudaMemcpyHostToDevice);

    out.numColors     = numColors;
    out.h_colorOffsets = std::move(h_offsets);

    std::printf("CVFEM coloring: %d colors for %zu elements\n", numColors, numElements);
    for (int c = 0; c < numColors; ++c)
        std::printf("  color %d: %d elements\n", c, out.colorSize(c));
}

} // namespace fem
} // namespace mars
