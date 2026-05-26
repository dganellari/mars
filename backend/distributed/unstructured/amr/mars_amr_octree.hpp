#pragma once

#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/reduce.h>
#include <thrust/transform.h>
#include <thrust/count.h>
#include <thrust/copy.h>
#include <mpi.h>

#include "cstone/cuda/cuda_utils.hpp"
#include "cstone/domain/domain.hpp"
#include "cstone/tree/csarray.hpp"

namespace mars
{
namespace amr
{

template<typename KeyType>
__global__ void assignLevelsFromTreeKernel(const KeyType* treeLeaves,
                                           size_t numLeaves,
                                           const KeyType* elemSfcKeys,
                                           int* elemLevels,
                                           size_t numElements)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements) return;

    KeyType key = elemSfcKeys[e];

    size_t lo = 0, hi = numLeaves;
    while (lo < hi)
    {
        size_t mid = (lo + hi) / 2;
        if (treeLeaves[mid + 1] <= key)
            lo = mid + 1;
        else
            hi = mid;
    }

    // SFC range of a leaf at level L is nodeRange(0) / 8^L
    if (lo < numLeaves)
    {
        KeyType range = treeLeaves[lo + 1] - treeLeaves[lo];
        int level     = 0;
        KeyType nr    = cstone::nodeRange<KeyType>(0);
        while (nr > range && level < 21)
        {
            nr >>= 3;
            level++;
        }
        elemLevels[e] = level;
    }
    else
    {
        elemLevels[e] = 0;
    }
}

template<typename RealType>
__global__ void markFromErrorKernel(const RealType* errorPerElement,
                                    const int* elemLevels,
                                    RealType refineThreshold,
                                    RealType coarsenThreshold,
                                    int maxLevel,
                                    uint8_t* marks,
                                    size_t numElements)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements) return;

    RealType err = errorPerElement[e];
    int level    = elemLevels[e];

    if (err > refineThreshold && level < maxLevel)
        marks[e] = 1;
    else if (err < coarsenThreshold && level > 0)
        marks[e] = 2;
    else
        marks[e] = 0;
}

// Enforce 1-irregular constraint: if a neighbor will be >1 level finer,
// propagate the refinement mark to prevent hanging nodes.
template<typename KeyType>
__global__ void enforce1IrregularKernel(const int* elemLevels,
                                        const uint8_t* marksIn,
                                        uint8_t* marksOut,
                                        const KeyType* conn0,
                                        const KeyType* conn1,
                                        const KeyType* conn2,
                                        const KeyType* conn3,
                                        const KeyType* conn4,
                                        const KeyType* conn5,
                                        const KeyType* conn6,
                                        const KeyType* conn7,
                                        const KeyType* nodeToElemOffsets,
                                        const KeyType* nodeToElemList,
                                        size_t numElements)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements) return;

    uint8_t myMark = marksIn[e];
    int myLevel   = elemLevels[e];
    int myNewLevel = myLevel + ((myMark > 0) ? 1 : 0);

    KeyType nodes[8] = {conn0[e], conn1[e], conn2[e], conn3[e], conn4[e], conn5[e], conn6[e], conn7[e]};

    uint8_t mark = myMark;

    for (int ni = 0; ni < 8; ++ni)
    {
        KeyType node  = nodes[ni];
        KeyType start = nodeToElemOffsets[node];
        KeyType end   = nodeToElemOffsets[node + 1];

        for (KeyType j = start; j < end; ++j)
        {
            KeyType neighbor = nodeToElemList[j];
            if (neighbor == e) continue;

            int neighborLevel    = elemLevels[neighbor];
            uint8_t neighborMark  = marksIn[neighbor];
            int neighborNewLevel = neighborLevel + ((neighborMark > 0) ? 1 : 0);

            if (neighborNewLevel > myNewLevel + 1)
            {
                mark = 1;
            }
        }
    }

    marksOut[e] = mark;
}

// Tet variant: 4 corner nodes per element (vs 8 for hex). Same algorithm.
template<typename KeyType>
__global__ void enforce1IrregularKernelTet(const int* elemLevels,
                                            const uint8_t* marksIn,
                                            uint8_t* marksOut,
                                            const KeyType* conn0,
                                            const KeyType* conn1,
                                            const KeyType* conn2,
                                            const KeyType* conn3,
                                            const KeyType* nodeToElemOffsets,
                                            const KeyType* nodeToElemList,
                                            size_t numElements)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements) return;

    uint8_t myMark = marksIn[e];
    int myLevel    = elemLevels[e];
    int myNewLevel = myLevel + ((myMark > 0) ? 1 : 0);

    KeyType nodes[4] = {conn0[e], conn1[e], conn2[e], conn3[e]};

    uint8_t mark = myMark;

    for (int ni = 0; ni < 4; ++ni)
    {
        KeyType node  = nodes[ni];
        KeyType start = nodeToElemOffsets[node];
        KeyType end   = nodeToElemOffsets[node + 1];

        for (KeyType j = start; j < end; ++j)
        {
            KeyType neighbor = nodeToElemList[j];
            if (neighbor == e) continue;

            int neighborLevel    = elemLevels[neighbor];
            uint8_t neighborMark = marksIn[neighbor];
            int neighborNewLevel = neighborLevel + ((neighborMark > 0) ? 1 : 0);

            if (neighborNewLevel > myNewLevel + 1) { mark = 1; }
        }
    }

    marksOut[e] = mark;
}

template<typename KeyType, typename RealType>
class AmrOctree
{
public:
    struct Config
    {
        RealType refineFraction  = 0.3;  // Doerfler: refine top 30% of error
        RealType coarsenFraction = 0.03; // Coarsen bottom 3%
        int maxLevel             = 5;    // Max refinement depth
        int blockSize            = 256;
        int irregularityIters    = 3;    // Iterations for 1-irregular propagation
    };

    template<typename ElementTag, typename AcceleratorTag>
    void initialize(const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain)
    {
        size_t numElements = domain.getElementCount();
        d_elemLevels_.resize(numElements);
        if (numElements == 0) return;

        // Derive per-element refinement level from cstone's focus-tree leaves.
        // Each element has an SFC key (from sync()); the leaf containing that
        // key has an SFC range nodeRange(0) / 8^L, which gives L.
        //
        // Before this fix, levels were zeroed on every initialize() and never
        // updated, so markFromErrorKernel always saw level=0. That broke
        // Doerfler refinement: every cell looked "refinable" forever, so each
        // adapt step refined the entire mesh uniformly (or stopped only when
        // cstone's bucketSize/tree-depth cap kicked in). After this fix, the
        // marker can correctly gate at maxLevel and reach a true adaptive
        // distribution of base + L1 + L2 cells.
        const auto& cstoneDom = domain.getDomain();
        const auto  leavesView = cstoneDom.focusTree().treeLeaves();
        size_t      numLeaves  = (leavesView.size() > 0) ? (leavesView.size() - 1) : 0;
        const KeyType* d_leaves = leavesView.data();

        const auto& d_elemKeys = domain.getElementSfcCodes();
        if (numLeaves == 0 || d_elemKeys.size() < numElements)
        {
            thrust::fill(thrust::device, d_elemLevels_.begin(), d_elemLevels_.end(), 0);
            return;
        }

        int numBlocks = (numElements + config_.blockSize - 1) / config_.blockSize;
        assignLevelsFromTreeKernel<KeyType><<<numBlocks, config_.blockSize>>>(
            d_leaves, numLeaves,
            thrust::raw_pointer_cast(d_elemKeys.data()),
            thrust::raw_pointer_cast(d_elemLevels_.data()),
            numElements);
        cudaDeviceSynchronize();
    }

    // Doerfler marking: refine top fraction of error, coarsen bottom fraction
    cstone::DeviceVector<uint8_t> markForRefinement(const RealType* d_errorPerElement, size_t numElements)
    {
        RealType maxError =
            thrust::reduce(thrust::device_pointer_cast(d_errorPerElement),
                           thrust::device_pointer_cast(d_errorPerElement + numElements), RealType(0),
                           thrust::maximum<RealType>());

        RealType refineThreshold  = config_.refineFraction * maxError;
        RealType coarsenThreshold = config_.coarsenFraction * maxError;

        cstone::DeviceVector<uint8_t> d_marks(numElements);
        int numBlocks = (numElements + config_.blockSize - 1) / config_.blockSize;

        markFromErrorKernel<<<numBlocks, config_.blockSize>>>(d_errorPerElement, d_elemLevels_.data(), refineThreshold,
                                                               coarsenThreshold, config_.maxLevel, d_marks.data(),
                                                               numElements);
        cudaDeviceSynchronize();

        return d_marks;
    }

    // Iteratively propagate marks so no two neighbors differ by >1 level.
    // Without this, hanging nodes break conformity of the mesh.
    template<typename ElementTag, typename AcceleratorTag>
    void enforce1Irregular(cstone::DeviceVector<uint8_t>& d_marks,
                           const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain)
    {
        size_t numElements = domain.getElementCount();
        int numBlocks      = (numElements + config_.blockSize - 1) / config_.blockSize;

        // Use local node IDs, not SFC keys, for indexing into nodeToElemOffsets
        const auto& d_conn = domain.getElementToNodeConnectivity();
        const auto& n2eOff = domain.getNodeToElementOffsets();
        const auto& n2eLst = domain.getNodeToElementList();

        cstone::DeviceVector<uint8_t> d_marksTemp(numElements);

        for (int iter = 0; iter < config_.irregularityIters; ++iter)
        {
            if constexpr (ElementTag::NodesPerElement == 8)
            {
                enforce1IrregularKernel<<<numBlocks, config_.blockSize>>>(
                    d_elemLevels_.data(), d_marks.data(), d_marksTemp.data(), std::get<0>(d_conn).data(),
                    std::get<1>(d_conn).data(), std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
                    std::get<4>(d_conn).data(), std::get<5>(d_conn).data(), std::get<6>(d_conn).data(),
                    std::get<7>(d_conn).data(), n2eOff.data(), n2eLst.data(), numElements);
            }
            else if constexpr (ElementTag::NodesPerElement == 4)
            {
                enforce1IrregularKernelTet<<<numBlocks, config_.blockSize>>>(
                    d_elemLevels_.data(), d_marks.data(), d_marksTemp.data(), std::get<0>(d_conn).data(),
                    std::get<1>(d_conn).data(), std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
                    n2eOff.data(), n2eLst.data(), numElements);
            }
            else
            {
                static_assert(ElementTag::NodesPerElement == 8 || ElementTag::NodesPerElement == 4,
                              "enforce1Irregular: only hex8 and tet4 are supported");
            }
            cudaDeviceSynchronize();

            auto m_b = thrust::device_pointer_cast(d_marks.data());
            auto m_e = thrust::device_pointer_cast(d_marks.data() + numElements);
            auto t_b = thrust::device_pointer_cast(d_marksTemp.data());
            auto t_e = thrust::device_pointer_cast(d_marksTemp.data() + numElements);
            bool changed = !thrust::equal(thrust::device, m_b, m_e, t_b);
            thrust::copy(thrust::device, t_b, t_e, m_b);

            if (!changed) break;
        }
    }

    const cstone::DeviceVector<int>& elementLevels() const { return d_elemLevels_; }

    Config& config() { return config_; }
    const Config& config() const { return config_; }

private:
    Config config_;
    cstone::DeviceVector<int> d_elemLevels_;
};

} // namespace amr
} // namespace mars
