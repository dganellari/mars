#pragma once

#include <thrust/device_ptr.h>
#include <thrust/functional.h>
#include <mpi.h>
#include <type_traits>

#include "cstone/domain/domain.hpp"
#include "cstone/tree/csarray.hpp"
#include "mars_amr_tree_layout.hpp"
#include "../utils/mars_gpu_reduce.hpp"

namespace mars::amr
{
// Cornerstone supplies spatial split decisions; mesh closure remains a separate step.
template<class KeyType, class RealType>
class OctreeNativeAmr
{
public:
    struct Config
    {
        unsigned bucketSize = 64;
        RealType refineFraction = RealType(0.3);
        int blockSize = 256;
    };

    template<class ElementTag, class AcceleratorTag>
    cstone::DeviceVector<uint8_t> markFromOctree(
        const RealType* error, const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain,
        const Config& config = Config{})
    {
        // This kernel's shared reduction requires a power-of-two block of at most 256 threads.
        if (!(config.refineFraction >= 0 && config.refineFraction <= 1) || config.bucketSize == 0 ||
            config.bucketSize > 65535 || config.blockSize < 1 || config.blockSize > 256 ||
            (config.blockSize & (config.blockSize - 1)))
            mars::abortAllRanks(EXIT_FAILURE);
        const size_t first = domain.startIndex(), last = domain.endIndex();
        RealType maximum = 0;
        if (first != last)
            maximum = mars::reduce(thrust::device_pointer_cast(error + first),
                                   thrust::device_pointer_cast(error + last), RealType(0),
                                   thrust::maximum<RealType>());
        MPI_Allreduce(MPI_IN_PLACE, &maximum, 1, std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT,
                      MPI_MAX, MPI_COMM_WORLD);

        const auto tree = domain.getDomain().globalTree();
        const auto* keys = domain.getElementSfcCodes().data();
        d_weighted_counts_.resize(tree.numLeafNodes);
        d_node_ops_.resize(tree.numLeafNodes + 1);
        cstone::DeviceVector<uint8_t> marks(domain.getElementCount(), 0);
        if (tree.numLeafNodes == 0) return marks;
        weight_owned_leaves<<<tree.numLeafNodes, config.blockSize>>>(tree.leaves, tree.numLeafNodes, keys, first, last,
                                                                     error, maximum, config.bucketSize,
                                                                     config.refineFraction, d_weighted_counts_.data());
        cudaCheckError();

        // Assignments end on global-leaf boundaries, so every populated leaf has one owner.
        // Counts stay >= bucket: sibling merging is impossible. Remote leaves may therefore
        // keep the floor count without changing the decisions for this rank's owned leaves.
        cstone::computeNodeOpsGpu(cstone::execution::gpuDefaultStream, tree.leaves, tree.numLeafNodes,
                                  d_weighted_counts_.data(), config.bucketSize, d_node_ops_.data());
        mark_owned_leaves<<<tree.numLeafNodes, config.blockSize>>>(tree.leaves, tree.numLeafNodes, keys, first, last,
                                                                   d_node_ops_.data(), d_weighted_counts_.data(),
                                                                   config.bucketSize, marks.data());
        cudaCheckError();
        return marks;
    }

private:
    // These match Cornerstone's existing global tree; reuse their capacity between adapt steps.
    cstone::DeviceVector<unsigned> d_weighted_counts_;
    cstone::DeviceVector<cstone::TreeNodeIndex> d_node_ops_;
};
} // namespace mars::amr
