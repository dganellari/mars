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
#include <type_traits>
#include <utility>

#include "cstone/cuda/cuda_utils.hpp"
#include "cstone/domain/domain.hpp"
#include "cstone/tree/csarray.hpp"
#include "mars_amr_tree_layout.hpp"
#include "../utils/mars_gpu_reduce.hpp"

namespace mars
{
namespace amr
{

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

template<typename KeyType, int N>
struct AmrCorners
{
    const KeyType* c[N];
};

// Shortest edge of each element: hexahedra in the refiner's corner numbering, tetrahedra all six edges
template<int N, typename KeyType, typename RealType>
__global__ void shortestEdgeKernel(AmrCorners<KeyType, N> conn, const RealType* x, const RealType* y,
                                   const RealType* z, size_t numElements, RealType* size)
{
    static_assert(N == 8 || N == 4, "hexahedra and tetrahedra");
    constexpr int numEdges     = N == 8 ? 12 : 6;
    constexpr int edges[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
                                  {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
    constexpr int tetEdges[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
    size_t e = blockIdx.x * size_t(blockDim.x) + threadIdx.x;
    if (e >= numElements) return;
    RealType shortest2 = 0;
    for (int q = 0; q < numEdges; ++q)
    {
        const KeyType a = conn.c[N == 8 ? edges[q][0] : tetEdges[q][0]][e];
        const KeyType b = conn.c[N == 8 ? edges[q][1] : tetEdges[q][1]][e];
        const RealType dx = x[b] - x[a], dy = y[b] - y[a], dz = z[b] - z[a];
        const RealType l2 = dx * dx + dy * dy + dz * dz;
        shortest2         = (q == 0 || l2 < shortest2) ? l2 : shortest2;
    }
    size[e] = sqrt(shortest2);
}

// An element is marked too when, after this refinement, an element sharing one of its nodes would be more than
// three times smaller: one level apart is a factor 2, two levels a factor 4. Sizes rather than levels: the levels
// of cornerstone's tree follow its bucket size, and a coarse element and its refined neighbour can share a leaf.
template<int N, typename KeyType, typename RealType>
__global__ void enforce1IrregularKernel(const RealType* size, const uint8_t* marksIn, uint8_t* marksOut,
                                        AmrCorners<KeyType, N> conn, const KeyType* nodeToElemOffsets,
                                        const KeyType* nodeToElemList, size_t numElements)
{
    size_t e = blockIdx.x * size_t(blockDim.x) + threadIdx.x;
    if (e >= numElements) return;

    const RealType mySize = size[e] * (marksIn[e] > 0 ? RealType(0.5) : RealType(1));
    uint8_t mark          = marksIn[e];
    for (int ni = 0; ni < N; ++ni)
    {
        const KeyType node = conn.c[ni][e];
        for (KeyType j = nodeToElemOffsets[node]; j < nodeToElemOffsets[node + 1]; ++j)
        {
            const KeyType neighbor = nodeToElemList[j];
            if (neighbor == e) continue;
            const RealType neighborSize = size[neighbor] * (marksIn[neighbor] > 0 ? RealType(0.5) : RealType(1));
            if (mySize > RealType(3) * neighborSize) mark = 1;
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

        // Acquire the current view after each rebuild; synchronization can replace its buffers.
        const auto tree = domain.getDomain().octreeProperties();
        if (tree.numLeafNodes == 0)
        {
            thrust::fill(thrust::device, d_elemLevels_.begin(), d_elemLevels_.end(), 0);
            return;
        }
        assign_leaf_levels<<<tree.numLeafNodes, config_.blockSize>>>(tree, d_elemLevels_.data());
        cudaCheckError();
    }

    // Doerfler marking: refine top fraction of error, coarsen bottom fraction
    cstone::DeviceVector<uint8_t> markForRefinement(const RealType* d_errorPerElement, size_t numElements)
    {
        RealType maxError = numElements ?
            mars::reduce(thrust::device_pointer_cast(d_errorPerElement),
                           thrust::device_pointer_cast(d_errorPerElement + numElements), RealType(0),
                           thrust::maximum<RealType>()) : RealType(0);
        // the same threshold on every rank, so the marked set does not depend on the partition
        MPI_Allreduce(MPI_IN_PLACE, &maxError, 1, std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT, MPI_MAX,
                      MPI_COMM_WORLD);
        // Empty ranks still take part in the global threshold reduction.
        if (numElements == 0) return {};

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

    // Propagates marks so no two elements sharing a node differ by more than one level after refinement: then a
    // hanging node sits at the midpoint of a coarse edge or the centre of a coarse face (HangingNodes). Sweeps this
    // rank's elements up to irregularityIters times, stopping when nothing changes; returns whether a mark of an
    // owned element changed. The marks of halo
    // elements come from their owners, so the caller repeats exchange and sweep until no rank changes a mark.
    template<typename ElementTag, typename AcceleratorTag>
    bool enforce1Irregular(cstone::DeviceVector<uint8_t>& d_marks,
                           const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain)
    {
        size_t numElements = domain.getElementCount();
        int numBlocks      = (numElements + config_.blockSize - 1) / config_.blockSize;
        bool ownChanged    = false;
        const size_t first = domain.startIndex(), last = domain.endIndex();

        // Use local node IDs, not SFC keys, for indexing into nodeToElemOffsets
        constexpr int N = ElementTag::NodesPerElement;
        static_assert(N == 8 || N == 4, "enforce1Irregular: only hex8 and tet4 are supported");
        const auto& d_conn = domain.getElementToNodeConnectivity();
        const auto& n2eOff = domain.getNodeToElementOffsets();
        const auto& n2eLst = domain.getNodeToElementList();
        AmrCorners<KeyType, N> corners;
        [&]<std::size_t... I>(std::index_sequence<I...>) { ((corners.c[I] = std::get<I>(d_conn).data()), ...); }
        (std::make_index_sequence<N>{});

        cstone::DeviceVector<uint8_t> d_marksTemp(numElements);
        cstone::DeviceVector<RealType> d_size(numElements);
        if (numElements > 0)
            shortestEdgeKernel<N, KeyType, RealType><<<numBlocks, config_.blockSize>>>(
                corners, domain.getNodeX().data(), domain.getNodeY().data(), domain.getNodeZ().data(), numElements,
                d_size.data());

        for (int iter = 0; iter < config_.irregularityIters && numElements > 0; ++iter)
        {
            enforce1IrregularKernel<N, KeyType, RealType><<<numBlocks, config_.blockSize>>>(
                d_size.data(), d_marks.data(), d_marksTemp.data(), corners, n2eOff.data(), n2eLst.data(), numElements);
            cudaDeviceSynchronize();

            auto m_b = thrust::device_pointer_cast(d_marks.data());
            auto m_e = thrust::device_pointer_cast(d_marks.data() + numElements);
            auto t_b = thrust::device_pointer_cast(d_marksTemp.data());
            auto t_e = thrust::device_pointer_cast(d_marksTemp.data() + numElements);
            bool changed = !thrust::equal(thrust::device, m_b, m_e, t_b);
            ownChanged |= changed && !thrust::equal(thrust::device, m_b + first, m_b + last, t_b + first);
            thrust::copy(thrust::device, t_b, t_e, m_b);

            if (!changed) break;
        }
        return ownChanged;
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
