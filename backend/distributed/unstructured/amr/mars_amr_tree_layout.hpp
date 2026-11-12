#pragma once

#include "cstone/tree/octree.hpp"
#include "cstone/util/array.hpp"

namespace mars::amr
{
template<class KeyType>
HOST_DEVICE_FUN util::array<size_t, 2> owned_leaf_range(const KeyType* leaves,
    cstone::TreeNodeIndex leaf, const KeyType* keys, size_t first, size_t last)
{
    if (first == last) return {first, last};
    const auto begin = stl::lower_bound(keys + first, keys + last, leaves[leaf]);
    const auto end = stl::lower_bound(begin, keys + last, leaves[leaf + 1]);
    return {size_t(begin - keys), size_t(end - keys)};
}

template<class RealType>
HOST_DEVICE_FUN unsigned weighted_leaf_count(RealType error, RealType maximum, unsigned bucket)
{
    const RealType scale = maximum > RealType(0) ? RealType(1) + error / maximum * RealType(bucket) : RealType(1);
    const auto weighted = static_cast<unsigned>(bucket * scale);
    return weighted > bucket ? weighted : bucket;
}

// The layout indexes elements, including halos; it must come from the latest sync.
template<class KeyType, class RealType>
HOST_DEVICE_FUN void assign_leaf_level(cstone::OctreeNsView<RealType, KeyType> tree,
                                      cstone::TreeNodeIndex leaf, unsigned lane, unsigned stride, int* levels)
{
    const int level = cstone::treeLevel(tree.leaves[leaf + 1] - tree.leaves[leaf]);
    for (size_t e = size_t(tree.layout[leaf]) + lane; e < tree.layout[leaf + 1]; e += stride)
        levels[e] = level;
}

#if defined(__CUDACC__) || defined(__HIPCC__)
template<class KeyType, class RealType>
__global__ void assign_leaf_levels(cstone::OctreeNsView<RealType, KeyType> tree, int* levels)
{
    const auto leaf = cstone::TreeNodeIndex(blockIdx.x);
    if (leaf < tree.numLeafNodes) assign_leaf_level(tree, leaf, threadIdx.x, blockDim.x, levels);
}
// Work once per leaf, using Cornerstone's sorted element keys instead of reconstructing node keys.
template<class KeyType, class RealType>
__global__ void weight_owned_leaves(const KeyType* leaves, cstone::TreeNodeIndex num_leaves,
                                   const KeyType* keys, size_t first, size_t last, const RealType* error,
                                   RealType global_max, unsigned bucket, unsigned* counts)
{
    const auto leaf = cstone::TreeNodeIndex(blockIdx.x);
    if (leaf >= num_leaves) return;
    __shared__ RealType partial[256];
    __shared__ size_t range[2];
    if (threadIdx.x == 0)
    {
        const auto bounds = owned_leaf_range(leaves, leaf, keys, first, last);
        range[0] = bounds[0];
        range[1] = bounds[1];
    }
    __syncthreads();
    RealType maximum = 0;
    for (size_t e = range[0] + threadIdx.x; e < range[1]; e += blockDim.x)
        maximum = error[e] > maximum ? error[e] : maximum;
    partial[threadIdx.x] = maximum;
    __syncthreads();
    for (unsigned width = blockDim.x / 2; width; width /= 2)
    {
        if (threadIdx.x < width && partial[threadIdx.x + width] > partial[threadIdx.x])
            partial[threadIdx.x] = partial[threadIdx.x + width];
        __syncthreads();
    }
    if (threadIdx.x == 0) counts[leaf] = weighted_leaf_count(partial[0], global_max, bucket);
}

template<class KeyType>
__global__ void mark_owned_leaves(const KeyType* leaves, cstone::TreeNodeIndex num_leaves,
                                 const KeyType* keys, size_t first, size_t last,
                                 const cstone::TreeNodeIndex* node_ops, uint8_t* marks)
{
    const auto leaf = cstone::TreeNodeIndex(blockIdx.x);
    if (leaf >= num_leaves) return;
    __shared__ size_t range[2];
    if (threadIdx.x == 0)
    {
        const auto bounds = owned_leaf_range(leaves, leaf, keys, first, last);
        range[0] = bounds[0];
        range[1] = bounds[1];
    }
    __syncthreads();
    const auto operation = node_ops[leaf + 1] - node_ops[leaf];
    const uint8_t mark = operation >= 8 ? 1 : operation == 0 ? 2 : 0;
    for (size_t e = range[0] + threadIdx.x; e < range[1]; e += blockDim.x) marks[e] = mark;
}

#endif
} // namespace mars::amr
