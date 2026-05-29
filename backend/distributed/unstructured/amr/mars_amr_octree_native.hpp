#pragma once

#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/transform.h>
#include <thrust/scan.h>
#include <thrust/fill.h>
#include <mpi.h>

#include "cstone/cuda/cuda_utils.hpp"
#include "cstone/domain/domain.hpp"
#include "cstone/tree/csarray.hpp"

namespace mars
{
namespace amr
{

// Compute element SFC key as min of its corner-node SFC keys (representative node convention)
template<typename KeyType>
__device__ KeyType computeElemSfc(const KeyType* nodeSfc,
                                  const KeyType* conn0, const KeyType* conn1,
                                  const KeyType* conn2, const KeyType* conn3,
                                  const KeyType* conn4, const KeyType* conn5,
                                  const KeyType* conn6, const KeyType* conn7,
                                  size_t e)
{
    KeyType k = nodeSfc[conn0[e]];
    KeyType k1 = nodeSfc[conn1[e]]; if (k1 < k) k = k1;
    KeyType k2 = nodeSfc[conn2[e]]; if (k2 < k) k = k2;
    KeyType k3 = nodeSfc[conn3[e]]; if (k3 < k) k = k3;
    KeyType k4 = nodeSfc[conn4[e]]; if (k4 < k) k = k4;
    KeyType k5 = nodeSfc[conn5[e]]; if (k5 < k) k = k5;
    KeyType k6 = nodeSfc[conn6[e]]; if (k6 < k) k = k6;
    KeyType k7 = nodeSfc[conn7[e]]; if (k7 < k) k = k7;
    return k;
}

// 4-corner overload for tet.
template<typename KeyType>
__device__ KeyType computeElemSfcTet(const KeyType* nodeSfc,
                                      const KeyType* conn0, const KeyType* conn1,
                                      const KeyType* conn2, const KeyType* conn3,
                                      size_t e)
{
    KeyType k = nodeSfc[conn0[e]];
    KeyType k1 = nodeSfc[conn1[e]]; if (k1 < k) k = k1;
    KeyType k2 = nodeSfc[conn2[e]]; if (k2 < k) k = k2;
    KeyType k3 = nodeSfc[conn3[e]]; if (k3 < k) k = k3;
    return k;
}

// Scatter per-element error to per-leaf max using atomic CAS
// (atomicMax not available for double)
template<typename KeyType, typename RealType>
__global__ void scatterErrorToLeavesKernel(const RealType* errorPerElement,
                                           const KeyType* nodeSfc,
                                           const KeyType* conn0, const KeyType* conn1,
                                           const KeyType* conn2, const KeyType* conn3,
                                           const KeyType* conn4, const KeyType* conn5,
                                           const KeyType* conn6, const KeyType* conn7,
                                           const KeyType* treeLeaves,
                                           size_t numLeaves,
                                           RealType* leafError,
                                           size_t startIdx,
                                           size_t endIdx)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x + startIdx;
    if (e >= endIdx) return;

    KeyType key  = computeElemSfc(nodeSfc, conn0, conn1, conn2, conn3, conn4, conn5, conn6, conn7, e);
    RealType err = errorPerElement[e];

    size_t lo = 0, hi = numLeaves;
    while (lo < hi)
    {
        size_t mid = (lo + hi) / 2;
        if (treeLeaves[mid + 1] <= key) lo = mid + 1;
        else hi = mid;
    }

    if (lo < numLeaves)
    {
        unsigned long long* addr = reinterpret_cast<unsigned long long*>(&leafError[lo]);
        unsigned long long assumed, old_val;
        old_val = *addr;
        do
        {
            assumed          = old_val;
            RealType old_err = __longlong_as_double(assumed);
            if (err <= old_err) break;
            old_val = atomicCAS(addr, assumed, __double_as_longlong(err));
        } while (assumed != old_val);
    }
}

// Tet variant of scatterErrorToLeavesKernel (4 corner nodes per element).
template<typename KeyType, typename RealType>
__global__ void scatterErrorToLeavesKernelTet(const RealType* errorPerElement,
                                               const KeyType* nodeSfc,
                                               const KeyType* conn0, const KeyType* conn1,
                                               const KeyType* conn2, const KeyType* conn3,
                                               const KeyType* treeLeaves,
                                               size_t numLeaves,
                                               RealType* leafError,
                                               size_t startIdx,
                                               size_t endIdx)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x + startIdx;
    if (e >= endIdx) return;

    KeyType key  = computeElemSfcTet(nodeSfc, conn0, conn1, conn2, conn3, e);
    RealType err = errorPerElement[e];

    size_t lo = 0, hi = numLeaves;
    while (lo < hi)
    {
        size_t mid = (lo + hi) / 2;
        if (treeLeaves[mid + 1] <= key) lo = mid + 1;
        else hi = mid;
    }

    if (lo < numLeaves)
    {
        unsigned long long* addr = reinterpret_cast<unsigned long long*>(&leafError[lo]);
        unsigned long long assumed, old_val;
        old_val = *addr;
        do
        {
            assumed          = old_val;
            RealType old_err = __longlong_as_double(assumed);
            if (err <= old_err) break;
            old_val = atomicCAS(addr, assumed, __double_as_longlong(err));
        } while (assumed != old_val);
    }
}

// Convert per-leaf error into pseudo-counts for cstone's rebalanceDecision.
// High-error leaves get inflated counts so cstone splits them.
// bucketSize controls the split threshold.
template<typename RealType>
__global__ void errorToWeightedCountsKernel(const RealType* leafError,
                                             const unsigned* realCounts,
                                             unsigned* weightedCounts,
                                             size_t numLeaves,
                                             RealType maxError,
                                             unsigned bucketSize)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numLeaves) return;

    unsigned count = realCounts[i];
    RealType err   = leafError[i];

    // Scale factor: leaves at maxError get count multiplied by bucketSize+1 (guaranteed split).
    // Leaves at zero error keep their real count.
    RealType scale = RealType(1);
    if (maxError > RealType(0))
        scale = RealType(1) + err / maxError * RealType(bucketSize);

    unsigned weighted = static_cast<unsigned>(count * scale);
    // Ensure at least the real count (never artificially merge)
    weightedCounts[i] = (weighted > count) ? weighted : count;
}

// Convert cstone nodeOps (split=8, merge=0, keep=1) to per-element marks
template<typename KeyType>
__global__ void leafOpsToElementMarksKernel(const KeyType* nodeSfc,
                                             const KeyType* conn0, const KeyType* conn1,
                                             const KeyType* conn2, const KeyType* conn3,
                                             const KeyType* conn4, const KeyType* conn5,
                                             const KeyType* conn6, const KeyType* conn7,
                                             const KeyType* treeLeaves,
                                             const cstone::TreeNodeIndex* nodeOps,
                                             size_t numLeaves,
                                             uint8_t* elementMarks,
                                             size_t startIdx,
                                             size_t endIdx)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x + startIdx;
    if (e >= endIdx) return;

    KeyType key = computeElemSfc(nodeSfc, conn0, conn1, conn2, conn3, conn4, conn5, conn6, conn7, e);

    size_t lo = 0, hi = numLeaves;
    while (lo < hi)
    {
        size_t mid = (lo + hi) / 2;
        if (treeLeaves[mid + 1] <= key) lo = mid + 1;
        else hi = mid;
    }

    if (lo < numLeaves)
    {
        // nodeOps is after exclusive_scan, so it holds output indices, not raw decisions.
        // We need the raw decisions. Recover: if nodeOps[lo+1] - nodeOps[lo] >= 8 -> split
        cstone::TreeNodeIndex diff = nodeOps[lo + 1] - nodeOps[lo];
        if (diff >= 8)
            elementMarks[e] = 1;
        else if (diff == 0)
            elementMarks[e] = 2;
        else
            elementMarks[e] = 0;
    }
    else
    {
        elementMarks[e] = 0;
    }
}

// Tet variant of leafOpsToElementMarksKernel (4 corner nodes per element).
template<typename KeyType>
__global__ void leafOpsToElementMarksKernelTet(const KeyType* nodeSfc,
                                                const KeyType* conn0, const KeyType* conn1,
                                                const KeyType* conn2, const KeyType* conn3,
                                                const KeyType* treeLeaves,
                                                const cstone::TreeNodeIndex* nodeOps,
                                                size_t numLeaves,
                                                uint8_t* elementMarks,
                                                size_t startIdx,
                                                size_t endIdx)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x + startIdx;
    if (e >= endIdx) return;

    KeyType key = computeElemSfcTet(nodeSfc, conn0, conn1, conn2, conn3, e);

    size_t lo = 0, hi = numLeaves;
    while (lo < hi)
    {
        size_t mid = (lo + hi) / 2;
        if (treeLeaves[mid + 1] <= key) lo = mid + 1;
        else hi = mid;
    }

    if (lo < numLeaves)
    {
        cstone::TreeNodeIndex diff = nodeOps[lo + 1] - nodeOps[lo];
        if (diff >= 8)      elementMarks[e] = 1;
        else if (diff == 0) elementMarks[e] = 2;
        else                elementMarks[e] = 0;
    }
    else
    {
        elementMarks[e] = 0;
    }
}

// Octree-native AMR: the octree IS the AMR data structure.
//
// Flow:
//   1. Scatter per-element error to cstone octree leaves (GPU atomic max)
//   2. Convert leaf errors to weighted counts (high error = heavy leaf)
//   3. Feed weighted counts into cstone's computeNodeOpsGpu (rebalanceDecision)
//   4. cstone decides split/merge per leaf based on its own bucket-size logic
//   5. Map leaf decisions back to per-element marks
//
// This means cstone's octree resolution directly controls mesh resolution.
// DD and AMR are unified: the same tree that partitions the domain also
// decides where to refine.
template<typename KeyType, typename RealType>
class OctreeNativeAmr
{
public:
    struct Config
    {
        unsigned bucketSize = 64;
        int blockSize       = 256;
    };

    template<typename ElementTag, typename AcceleratorTag>
    cstone::DeviceVector<uint8_t> markFromOctree(
        const RealType* d_errorPerElement,
        const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain,
        const Config& config = Config{})
    {
        int rank, nranks;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nranks);
        std::cerr << "OCTNATIVE r" << rank << " entry" << std::endl; std::cerr.flush();

        size_t numElements = domain.getElementCount();
        size_t startIdx    = domain.startIndex();
        size_t endIdx      = domain.endIndex();
        size_t localCount  = endIdx - startIdx;
        std::cerr << "OCTNATIVE r" << rank << " got domain bounds" << std::endl; std::cerr.flush();

        const auto& cstoneDomain = domain.getDomain();
        std::cerr << "OCTNATIVE r" << rank << " got cstone domain" << std::endl; std::cerr.flush();
        // Use the GLOBAL tree (identical on all ranks). cstone's focused tree
        // differs per rank outside the assigned slice, so it can't be used
        // directly for an MPI allreduce.
        auto globalTreeView      = cstoneDomain.globalTree();
        std::cerr << "OCTNATIVE r" << rank << " got globalTree" << std::endl; std::cerr.flush();
        size_t numLeaves         = globalTreeView.numLeafNodes;
        const KeyType* treeLeaves = globalTreeView.leaves;

        const auto& d_sfcMap = domain.getLocalToGlobalSfcMap();
        std::cerr << "OCTNATIVE r" << rank << " got sfcMap" << std::endl; std::cerr.flush();
        std::cerr << "OCTNATIVE r" << rank << " numLeaves=" << numLeaves
                  << " local=" << localCount << std::endl; std::cerr.flush();

        // Step 1: scatter LOCAL element errors into the global-leaf array.
        //         Only local elements [startIdx, endIdx) contribute here;
        //         halo elements would be double-counted by their owners.
        cstone::DeviceVector<RealType> d_leafError(numLeaves, RealType(0));

        // Get element-to-node connectivity (local node IDs) for SFC computation
        const auto& d_conn = domain.getElementToNodeConnectivity();

        int numBlocks = (localCount + config.blockSize - 1) / config.blockSize;
        if (numBlocks > 0)
        {
            if constexpr (ElementTag::NodesPerElement == 8)
            {
                scatterErrorToLeavesKernel<KeyType, RealType><<<numBlocks, config.blockSize>>>(
                    d_errorPerElement, d_sfcMap.data(),
                    std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
                    std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
                    std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
                    std::get<6>(d_conn).data(), std::get<7>(d_conn).data(),
                    treeLeaves, numLeaves,
                    d_leafError.data(), startIdx, endIdx);
            }
            else if constexpr (ElementTag::NodesPerElement == 4)
            {
                scatterErrorToLeavesKernelTet<KeyType, RealType><<<numBlocks, config.blockSize>>>(
                    d_errorPerElement, d_sfcMap.data(),
                    std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
                    std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
                    treeLeaves, numLeaves,
                    d_leafError.data(), startIdx, endIdx);
            }
            cudaDeviceSynchronize();
        }

        std::cerr << "OCTNATIVE r" << rank << " step1 done" << std::endl; std::cerr.flush();

        // Step 2: Allreduce(MAX) over the full global-leaf error array so all
        //         ranks see consistent per-leaf error values. Sizes match
        //         because numLeaves is the same global tree on every rank.
        if (nranks > 1)
        {
            std::vector<RealType> h_local(numLeaves), h_global(numLeaves);
            cudaMemcpy(h_local.data(), d_leafError.data(), numLeaves * sizeof(RealType),
                       cudaMemcpyDeviceToHost);

            auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
            std::cerr << "OCTNATIVE r" << rank << " calling allreduce, n=" << numLeaves << std::endl; std::cerr.flush();
            MPI_Allreduce(h_local.data(), h_global.data(), numLeaves, mpiType, MPI_MAX, MPI_COMM_WORLD);
            std::cerr << "OCTNATIVE r" << rank << " allreduce done" << std::endl; std::cerr.flush();

            cudaMemcpy(d_leafError.data(), h_global.data(), numLeaves * sizeof(RealType),
                       cudaMemcpyHostToDevice);
        }

        // Step 3: global max error
        auto le_b             = thrust::device_pointer_cast(d_leafError.data());
        auto le_e             = thrust::device_pointer_cast(d_leafError.data() + numLeaves);
        RealType maxLeafError = thrust::reduce(thrust::device, le_b, le_e,
                                                RealType(0), thrust::maximum<RealType>());

        // Step 4: weighted counts. cstone's `OctreeView` (returned by
        // `globalTree()`) does not expose per-leaf particle counts publicly,
        // so we use a uniform placeholder = bucketSize. The error term in
        // `errorToWeightedCountsKernel` drives the split decision — leaves
        // with high error get inflated above bucketSize and split.
        cstone::DeviceVector<unsigned> d_realCounts(numLeaves, config.bucketSize);

        cstone::DeviceVector<unsigned> d_weightedCounts(numLeaves);
        int leafBlocks = (numLeaves + config.blockSize - 1) / config.blockSize;
        errorToWeightedCountsKernel<<<leafBlocks, config.blockSize>>>(
            d_leafError.data(), d_realCounts.data(), d_weightedCounts.data(),
            numLeaves, maxLeafError, config.bucketSize);
        cudaDeviceSynchronize();

        // Step 5: cstone's rebalanceDecision (same function it uses for octree update).
        // Identical inputs across ranks => identical decisions => no comm needed.
        cstone::DeviceVector<cstone::TreeNodeIndex> d_nodeOps(numLeaves + 1);
        cstone::computeNodeOpsGpu(treeLeaves, numLeaves, d_weightedCounts.data(),
                                   config.bucketSize, d_nodeOps.data());

        // Step 6: map leaf decisions to element marks for LOCAL elements only.
        cstone::DeviceVector<uint8_t> d_marks(numElements, 0);
        if (numBlocks > 0)
        {
            if constexpr (ElementTag::NodesPerElement == 8)
            {
                leafOpsToElementMarksKernel<<<numBlocks, config.blockSize>>>(
                    d_sfcMap.data(),
                    std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
                    std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
                    std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
                    std::get<6>(d_conn).data(), std::get<7>(d_conn).data(),
                    treeLeaves, d_nodeOps.data(),
                    numLeaves, d_marks.data(), startIdx, endIdx);
            }
            else if constexpr (ElementTag::NodesPerElement == 4)
            {
                leafOpsToElementMarksKernelTet<<<numBlocks, config.blockSize>>>(
                    d_sfcMap.data(),
                    std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
                    std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
                    treeLeaves, d_nodeOps.data(),
                    numLeaves, d_marks.data(), startIdx, endIdx);
            }
            cudaDeviceSynchronize();
        }

        return d_marks;
    }
};

} // namespace amr
} // namespace mars
