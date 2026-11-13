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

template<typename KeyType, typename RealType>
__global__ void markFromErrorKernel(const RealType* errorPerElement,
                                    const KeyType* elemLevels,
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

template<typename KeyType, typename RealType>
class AmrOctree
{
public:
    struct Config
    {
        RealType refineFraction  = 0.3;  // Maximum-relative threshold (historical Doerfler API name)
        RealType coarsenFraction = 0.03; // Coarsen bottom 3%
        int maxLevel             = 5;    // Max refinement depth
        int blockSize            = 256;
    };

    template<typename ElementTag, typename AcceleratorTag>
    void initialize(const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain)
    {
        const auto& generations = domain.refinement_generations();
        d_elemLevels_           = cstone::DeviceVector<KeyType>(domain.getElementCount(), KeyType(0));
        if (!generations.empty())
            thrust::copy(thrust::device, thrust::device_pointer_cast(generations.data()),
                         thrust::device_pointer_cast(generations.data() + generations.size()),
                         thrust::device_pointer_cast(d_elemLevels_.data()));
    }

    // Maximum-relative marking; this is not bulk Doerfler marking.
    cstone::DeviceVector<uint8_t>
    markForRefinement(const RealType* d_errorPerElement, size_t numElements, size_t first = 0, size_t last = SIZE_MAX)
    {
        if (last == SIZE_MAX) last = numElements;
        RealType maxError = first != last ? mars::reduce(thrust::device_pointer_cast(d_errorPerElement + first),
                                                         thrust::device_pointer_cast(d_errorPerElement + last),
                                                         RealType(0), thrust::maximum<RealType>())
                                          : RealType(0);
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
        cudaCheckError();

        return d_marks;
    }

    const cstone::DeviceVector<KeyType>& elementLevels() const { return d_elemLevels_; }

    Config& config() { return config_; }
    const Config& config() const { return config_; }

private:
    Config config_;
    cstone::DeviceVector<KeyType> d_elemLevels_;
};

} // namespace amr
} // namespace mars
