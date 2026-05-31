#pragma once

#include <memory>
#include <chrono>
#include <vector>
#include <mpi.h>

#include "backend/distributed/unstructured/domain.hpp"
#include "mars_amr_octree.hpp"
#include "mars_amr_octree_native.hpp"
#include "mars_amr_hex_refine.hpp"
#include "mars_amr_tet_refine.hpp"
#include "mars_amr_octree_refine.hpp"
#include "mars_amr_solution_transfer.hpp"
#include "mars_amr_error_indicator.hpp"

namespace mars
{
namespace amr
{

struct AmrStats
{
    size_t elementsBeforeLocal  = 0;
    size_t elementsAfterLocal   = 0;
    size_t elementsBeforeGlobal = 0;
    size_t elementsAfterGlobal  = 0;
    size_t nodesBeforeGlobal    = 0;
    size_t nodesAfterGlobal     = 0;
    size_t elementsRefined      = 0;
    int level                   = 0;
    double errorNorm            = 0;
    float totalTimeMs           = 0;

    // Sub-phase timings (ms). Each phase is wall-clock with cudaDeviceSync +
    // MPI_Barrier so it reflects all-rank work, not just rank-0.
    float markTimeMs            = 0;   // marking (Doerfler / OctreeNative)
    float refineTimeMs          = 0;   // OctreeAlignedRefine::refineLocal
    float transferTimeMs        = 0;   // solution transfer (interp + presync coord save)
    float rebuildTimeMs         = 0;   // rebuildDomainFromDevice (cstone sync)
    float reorderTimeMs         = 0;   // post-sync key re-encode + binary-search lookup
    float haloFillTimeMs        = 0;   // exchangeNodeHalo on transferred fields
};

// How to decide which elements to refine
enum class MarkingStrategy
{
    Doerfler,     // Error-driven Doerfler marking + 1-irregular enforcement
    OctreeNative  // Cstone octree leaf-level split/merge decisions drive element marks
};

// Refiner traits: selects the right refiner + transfer for each element type
template<typename ElementTag, typename KeyType, typename RealType>
struct RefinerTraits;

template<typename KeyType, typename RealType>
struct RefinerTraits<HexTag, KeyType, RealType>
{
    using Refiner = HexRefiner<KeyType, RealType>;
    using Result  = typename Refiner::Result;

    static Result refine(const ElementDomain<HexTag, RealType, KeyType, cstone::GpuTag>& domain,
                         const uint8_t* marks, size_t numElements,
                         const RealType* nodeX, const RealType* nodeY, const RealType* nodeZ,
                         size_t numNodes, int blockSize)
    {
        const auto& d_conn = domain.getElementToNodeConnectivity();
        return Refiner::refine(
            std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
            std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
            std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
            std::get<6>(d_conn).data(), std::get<7>(d_conn).data(),
            marks, numElements, nodeX, nodeY, nodeZ, numNodes, blockSize);
    }

    // Refine only the local element range [startIdx, endIdx). Node coords
    // remain the full per-rank node array (refinement may reference nodes
    // outside the local element range when local elements share corners
    // with halo elements). The output connectivity references node IDs in
    // [0, numNodes) plus newly created nodes, suitable for handing back to
    // ElementDomain/cstone for a fresh sync.
    static Result refineLocal(const ElementDomain<HexTag, RealType, KeyType, cstone::GpuTag>& domain,
                              const uint8_t* marks, size_t startIdx, size_t endIdx,
                              const RealType* nodeX, const RealType* nodeY, const RealType* nodeZ,
                              size_t numNodes, int blockSize)
    {
        const auto& d_conn = domain.getElementToNodeConnectivity();
        size_t local       = endIdx - startIdx;
        return Refiner::refine(
            std::get<0>(d_conn).data() + startIdx, std::get<1>(d_conn).data() + startIdx,
            std::get<2>(d_conn).data() + startIdx, std::get<3>(d_conn).data() + startIdx,
            std::get<4>(d_conn).data() + startIdx, std::get<5>(d_conn).data() + startIdx,
            std::get<6>(d_conn).data() + startIdx, std::get<7>(d_conn).data() + startIdx,
            marks + startIdx, local, nodeX, nodeY, nodeZ, numNodes, blockSize);
    }

    static cstone::DeviceVector<RealType> transferField(
        const ElementDomain<HexTag, RealType, KeyType, cstone::GpuTag>& domain,
        const Result& refined, const RealType* oldField, size_t numElements, int blockSize)
    {
        const auto& d_conn = domain.getElementToNodeConnectivity();
        return SolutionTransfer<KeyType, RealType>::transferByParentage(
            std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
            std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
            std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
            std::get<6>(d_conn).data(), std::get<7>(d_conn).data(),
            refined.d_marks.data(), refined.d_markedPrefix.data(), numElements,
            oldField, refined.oldNumNodes,
            refined.d_edgeKeys.data(), refined.numUniqueEdges,
            refined.d_faceKeysLo.data(), refined.d_faceKeysHi.data(), refined.numUniqueFaces,
            refined.numMarked, refined.numNodes, blockSize);
    }
};

template<typename KeyType, typename RealType>
struct RefinerTraits<TetTag, KeyType, RealType>
{
    using Refiner = TetRefiner<KeyType, RealType>;
    using Result  = typename Refiner::Result;

    static Result refine(const ElementDomain<TetTag, RealType, KeyType, cstone::GpuTag>& domain,
                         const uint8_t* marks, size_t numElements,
                         const RealType* nodeX, const RealType* nodeY, const RealType* nodeZ,
                         size_t numNodes, int blockSize)
    {
        const auto& d_conn = domain.getElementToNodeConnectivity();
        return Refiner::refine(
            std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
            std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
            marks, numElements, nodeX, nodeY, nodeZ, numNodes, blockSize);
    }

    static Result refineLocal(const ElementDomain<TetTag, RealType, KeyType, cstone::GpuTag>& domain,
                              const uint8_t* marks, size_t startIdx, size_t endIdx,
                              const RealType* nodeX, const RealType* nodeY, const RealType* nodeZ,
                              size_t numNodes, int blockSize)
    {
        const auto& d_conn = domain.getElementToNodeConnectivity();
        size_t local       = endIdx - startIdx;
        return Refiner::refine(
            std::get<0>(d_conn).data() + startIdx, std::get<1>(d_conn).data() + startIdx,
            std::get<2>(d_conn).data() + startIdx, std::get<3>(d_conn).data() + startIdx,
            marks + startIdx, local, nodeX, nodeY, nodeZ, numNodes, blockSize);
    }

    static cstone::DeviceVector<RealType> transferField(
        const ElementDomain<TetTag, RealType, KeyType, cstone::GpuTag>& domain,
        const Result& refined, const RealType* oldField, size_t numElements, int blockSize)
    {
        const auto& d_conn = domain.getElementToNodeConnectivity();
        return SolutionTransfer<KeyType, RealType>::transferByParentageTet(
            std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
            std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
            refined.d_marks.data(), numElements,
            oldField, refined.oldNumNodes,
            refined.d_edgeKeys.data(), refined.numUniqueEdges,
            refined.numNodes, blockSize);
    }
};

// Compute bounding box from device coordinate arrays via thrust + MPI
template<typename RealType>
cstone::Box<RealType> computeBoundingBoxGpu(cstone::DeviceVector<RealType>& d_x,
                                             cstone::DeviceVector<RealType>& d_y,
                                             cstone::DeviceVector<RealType>& d_z,
                                             RealType padding = 0.05)
{
    RealType lxmin = thrust::reduce(d_x.begin(), d_x.end(), RealType(1e30), thrust::minimum<RealType>());
    RealType lxmax = thrust::reduce(d_x.begin(), d_x.end(), RealType(-1e30), thrust::maximum<RealType>());
    RealType lymin = thrust::reduce(d_y.begin(), d_y.end(), RealType(1e30), thrust::minimum<RealType>());
    RealType lymax = thrust::reduce(d_y.begin(), d_y.end(), RealType(-1e30), thrust::maximum<RealType>());
    RealType lzmin = thrust::reduce(d_z.begin(), d_z.end(), RealType(1e30), thrust::minimum<RealType>());
    RealType lzmax = thrust::reduce(d_z.begin(), d_z.end(), RealType(-1e30), thrust::maximum<RealType>());

    auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    RealType gxmin, gxmax, gymin, gymax, gzmin, gzmax;
    MPI_Allreduce(&lxmin, &gxmin, 1, mpiType, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&lxmax, &gxmax, 1, mpiType, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&lymin, &gymin, 1, mpiType, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&lymax, &gymax, 1, mpiType, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&lzmin, &gzmin, 1, mpiType, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&lzmax, &gzmax, 1, mpiType, MPI_MAX, MPI_COMM_WORLD);

    RealType rx = gxmax - gxmin, ry = gymax - gymin, rz = gzmax - gzmin;
    return cstone::Box<RealType>(gxmin - padding * rx, gxmax + padding * rx,
                                  gymin - padding * ry, gymax + padding * ry,
                                  gzmin - padding * rz, gzmax + padding * rz);
}

// Generic AMR manager: works with hex, tet, or any supported element type.
// Supports two marking strategies: Doerfler (error-driven) and OctreeNative (cstone-driven).
template<typename ElementTag, typename KeyType, typename RealType>
class AmrManager
{
public:
    using Domain    = ElementDomain<ElementTag, RealType, KeyType, cstone::GpuTag>;
    using DomainPtr = std::unique_ptr<Domain>;
    using DevVector = cstone::DeviceVector<RealType>;
    using Traits    = RefinerTraits<ElementTag, KeyType, RealType>;

    struct Config
    {
        int maxLevels            = 5;
        RealType refineFraction  = 0.3;
        RealType coarsenFraction = 0.03;
        RealType errorTolerance  = 1e-6;
        int bucketSize           = 64;
        int blockSize            = 256;
        int irregularityIters    = 3;
        MarkingStrategy strategy = MarkingStrategy::OctreeNative;
    };

    AmrManager(const Config& config = Config{}) : config_(config), currentLevel_(0)
    {
        typename AmrOctree<KeyType, RealType>::Config octConfig;
        octConfig.refineFraction    = config_.refineFraction;
        octConfig.coarsenFraction   = config_.coarsenFraction;
        octConfig.maxLevel          = config_.maxLevels;
        octConfig.blockSize         = config_.blockSize;
        octConfig.irregularityIters = config_.irregularityIters;
        amrOctree_.config()         = octConfig;
    }

    // Initialization timings populated by initialize(). All times in ms,
    // wall-clock with cudaDeviceSync + MPI_Barrier so they reflect all-rank
    // work, not just rank-0.
    //
    // domainSyncTimeMs file read + bbox + GPU upload + cstone sync (SFC sort
    //                  + redistribute + element halo build). Currently lumped
    //                  because ElementDomain's file ctor inlines these steps.
    //                  TODO: split file I/O into its own bucket.
    // haloTopoTimeMs   lazy build of HaloData + NodeHaloTopology (Allgatherv
    //                  + Alltoallv).
    // adjacencyTimeMs  element-to-node and node-to-element CSR build.
    // coordCacheTimeMs SFC-decoded node coord caching.
    // octreeTimeMs     AmrOctree state initialization (refinement-level
    //                  tracking).
    struct InitTimings {
        float domainSyncTimeMs = 0;
        float haloTopoTimeMs   = 0;
        float adjacencyTimeMs  = 0;
        float coordCacheTimeMs = 0;
        float octreeTimeMs     = 0;
        float totalMs          = 0;
    };
    const InitTimings& initTimings() const { return initTimings_; }

    void initialize(const std::string& meshFile, int rank, int numRanks,
                    int periodicAxesMask = 0,
                    RealType periodicBoxLo = RealType(0),
                    RealType periodicBoxHi = RealType(0))
    {
        rank_     = rank;
        numRanks_ = numRanks;
        auto t0 = std::chrono::high_resolution_clock::now();
        auto lap = [&t0]() -> float {
            cudaDeviceSynchronize();
            MPI_Barrier(MPI_COMM_WORLD);
            auto t1 = std::chrono::high_resolution_clock::now();
            float ms = std::chrono::duration<float, std::milli>(t1 - t0).count();
            t0 = t1;
            return ms;
        };

        // File ctor lumps file I/O with bbox + sync; ElementDomain prints an
        // intermediate "Created bounding box" message that splits them in the
        // log, but we can't separate the times from outside without restructuring.
        // For now: keep them combined in domainSyncTimeMs; fileReadTimeMs stays 0.
        // periodicAxesMask / periodicBoxLo / periodicBoxHi are threaded through
        // to ElementDomain for MFEM-style mesh-level periodic identification.
        domain_ = std::make_unique<Domain>(meshFile, rank, numRanks, true,
                                            config_.bucketSize,
                                            /*bucketSizeFocus*/ 8,
                                            periodicAxesMask,
                                            periodicBoxLo, periodicBoxHi);
        initTimings_.domainSyncTimeMs = lap();

        (void)domain_->getNodeOwnershipMap();   // builds HaloData + NodeHaloTopology
        initTimings_.haloTopoTimeMs = lap();

        (void)domain_->getElementToNodeConnectivity();  // builds adjacency CSRs
        initTimings_.adjacencyTimeMs = lap();

        domain_->cacheNodeCoordinates();
        initTimings_.coordCacheTimeMs = lap();

        amrOctree_.initialize(*domain_);
        initTimings_.octreeTimeMs = lap();

        initTimings_.totalMs = initTimings_.domainSyncTimeMs +
                                initTimings_.haloTopoTimeMs +
                                initTimings_.adjacencyTimeMs +
                                initTimings_.coordCacheTimeMs +
                                initTimings_.octreeTimeMs;
        currentLevel_ = 0;
    }

    void initialize(DomainPtr domain, int rank, int numRanks)
    {
        rank_     = rank;
        numRanks_ = numRanks;
        domain_   = std::move(domain);
        // Match working CVFEM init order: ownership first, then conn, coords, AMR state
        std::cerr << "AMR r" << rank << " getNodeOwnershipMap..." << std::endl; std::cerr.flush();
        (void)domain_->getNodeOwnershipMap();
        cudaDeviceSynchronize();
        std::cerr << "AMR r" << rank << " ownership done" << std::endl; std::cerr.flush();
        (void)domain_->getElementToNodeConnectivity();
        std::cerr << "AMR r" << rank << " adjacency done" << std::endl; std::cerr.flush();
        domain_->cacheNodeCoordinates();
        std::cerr << "AMR r" << rank << " coords done" << std::endl; std::cerr.flush();
        amrOctree_.initialize(*domain_);
        std::cerr << "AMR r" << rank << " octree done" << std::endl; std::cerr.flush();
        currentLevel_ = 0;
    }

    Domain& domain() { return *domain_; }
    const Domain& domain() const { return *domain_; }
    AmrOctree<KeyType, RealType>& octree() { return amrOctree_; }

    AmrStats adaptMesh(const RealType* d_errorPerElement,
                       const RealType* d_nodeSolution,
                       DevVector& d_newNodeSolution)
    {
        std::vector<const RealType*> oldFields = {d_nodeSolution};
        std::vector<DevVector*> newFields      = {&d_newNodeSolution};
        return adaptMeshMultiField(d_errorPerElement, oldFields, newFields);
    }

    AmrStats adaptMeshMultiField(const RealType* d_errorPerElement,
                                  const std::vector<const RealType*>& oldFields,
                                  const std::vector<DevVector*>& newFields)
    {
        AmrStats stats;
        stats.level = currentLevel_;
        auto totalStart = std::chrono::high_resolution_clock::now();

        size_t numElements = domain_->getElementCount();
        size_t numNodes    = domain_->getNodeCount();
        stats.elementsBeforeLocal = domain_->localElementCount();
        stats.nodesBeforeGlobal   = numNodes;

        std::cerr << "ADAPT r" << rank_ << " before barrier" << std::endl; std::cerr.flush();
        MPI_Barrier(MPI_COMM_WORLD);
        std::cerr << "ADAPT r" << rank_ << " after barrier, before allreduce" << std::endl; std::cerr.flush();
        MPI_Allreduce(&stats.elementsBeforeLocal, &stats.elementsBeforeGlobal, 1,
                       MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
        std::cerr << "ADAPT r" << rank_ << " after allreduce" << std::endl; std::cerr.flush();

        size_t startIdx = domain_->startIndex();
        size_t endIdx   = domain_->endIndex();

        // Phase timer: cudaDeviceSync + MPI_Barrier on each lap so the wall
        // time reflects all-rank work, not just rank-0.
        auto phaseT0 = std::chrono::high_resolution_clock::now();
        auto phaseLap = [&phaseT0]() -> float {
            cudaDeviceSynchronize();
            MPI_Barrier(MPI_COMM_WORLD);
            auto t1 = std::chrono::high_resolution_clock::now();
            float ms = std::chrono::duration<float, std::milli>(t1 - phaseT0).count();
            phaseT0 = t1;
            return ms;
        };

        cstone::DeviceVector<uint8_t> d_marks;
        if (config_.strategy == MarkingStrategy::OctreeNative)
        {
            std::cerr << "ADAPT r" << rank_ << " entering markFromOctree" << std::endl; std::cerr.flush();
            // OctreeNative path: cstone's octree is the AMR data structure.
            // Cross-rank consistency comes from the octree itself - all ranks
            // see the same leaf decisions because they reduce leaf errors via
            // MPI_Allreduce on the assigned slice, then run cstone's
            // computeNodeOpsGpu identically. Per-leaf split conformity is the
            // octree's invariant; no extra mark exchange or 1-irregular needed.
            d_marks = octreeNativeAmr_.markFromOctree(d_errorPerElement, *domain_);
            std::cerr << "ADAPT r" << rank_ << " markFromOctree returned" << std::endl; std::cerr.flush();
        }
        else
        {
            // Doerfler path: error threshold + 1-irregular enforcement at
            // mesh level. Multi-rank consistency requires halo exchange of
            // marks both before and after 1-irregular propagation.
            // cstone instantiates haloExchangeGpu for `int` (not `uint8_t`),
            // so we promote marks to int for the exchange step.
            d_marks = amrOctree_.markForRefinement(d_errorPerElement, numElements);

            if (startIdx > 0)
                cudaMemset(d_marks.data(), 0, startIdx * sizeof(uint8_t));
            if (endIdx < numElements)
                cudaMemset(d_marks.data() + endIdx, 0, (numElements - endIdx) * sizeof(uint8_t));

            auto exchangeMarks = [&]()
            {
                if (numRanks_ <= 1) return;
                // Promote to int for the exchange
                cstone::DeviceVector<int> d_marksInt(numElements);
                thrust::transform(thrust::device,
                                   thrust::device_pointer_cast(d_marks.data()),
                                   thrust::device_pointer_cast(d_marks.data() + numElements),
                                   thrust::device_pointer_cast(d_marksInt.data()),
                                   [] __device__(uint8_t m) -> int { return int(m); });
                cstone::DeviceVector<int> sendBuf, recvBuf;
                domain_->exchangeHalos(std::tie(d_marksInt), sendBuf, recvBuf);
                thrust::transform(thrust::device,
                                   thrust::device_pointer_cast(d_marksInt.data()),
                                   thrust::device_pointer_cast(d_marksInt.data() + numElements),
                                   thrust::device_pointer_cast(d_marks.data()),
                                   [] __device__(int m) -> uint8_t { return uint8_t(m); });
            };

            exchangeMarks();
            amrOctree_.enforce1Irregular(d_marks, *domain_);

            if (startIdx > 0)
                cudaMemset(d_marks.data(), 0, startIdx * sizeof(uint8_t));
            if (endIdx < numElements)
                cudaMemset(d_marks.data() + endIdx, 0, (numElements - endIdx) * sizeof(uint8_t));

            exchangeMarks();
        }
        stats.markTimeMs = phaseLap();

        auto marks_b = thrust::device_pointer_cast(d_marks.data());
        auto marks_e = thrust::device_pointer_cast(d_marks.data() + d_marks.size());
        size_t localRefined =
            thrust::count_if(thrust::device, marks_b, marks_e, [] __device__(uint8_t m) { return m > 0; });
        MPI_Allreduce(&localRefined, &stats.elementsRefined, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

        // Copy error to a DevVector so globalErrorNorm can iterate it
        DevVector d_errCopy(numElements);
        cudaMemcpy(d_errCopy.data(), d_errorPerElement, numElements * sizeof(RealType), cudaMemcpyDeviceToDevice);
        stats.errorNorm = ErrorIndicator<KeyType, RealType>::globalErrorNorm(d_errCopy);

        if (stats.elementsRefined == 0)
        {
            if (rank_ == 0) std::cout << "  AMR: No elements marked. Mesh unchanged.\n";
            for (size_t f = 0; f < oldFields.size(); ++f)
            {
                newFields[f]->resize(numNodes);
                cudaMemcpy(newFields[f]->data(), oldFields[f], numNodes * sizeof(RealType), cudaMemcpyDeviceToDevice);
            }
            auto totalEnd     = std::chrono::high_resolution_clock::now();
            stats.totalTimeMs = std::chrono::duration<float, std::milli>(totalEnd - totalStart).count();
            return stats;
        }

        const auto& d_x = domain_->getNodeX();
        const auto& d_y = domain_->getNodeY();
        const auto& d_z = domain_->getNodeZ();
        const auto& d_conn = domain_->getElementToNodeConnectivity();

        // Octree-aligned refinement: each rank emits child elements only for
        // its locally-owned, marked elements. Coords + connectivity flow to
        // cstone via the device-data ElementDomain constructor below; cstone
        // re-derives SFC keys, redistributes by SFC, and rebuilds halos.
        // No host round-trip, no per-rank "global mesh" reconstruction.
        // Element-type dispatch is in doRefineLocal / doTransferSolution so
        // the hex 8-conn body never gets parsed in the TetTag instantiation
        // (and vice versa).
        using Refine = OctreeAlignedRefine<KeyType, RealType, ElementTag>;
        size_t localCount = endIdx - startIdx;
        typename Refine::Result refined =
            doRefineLocal(d_conn, d_x, d_y, d_z, d_marks, startIdx, localCount, numNodes);
        stats.refineTimeMs = phaseLap();

        // Solution transfer onto the pre-cstone-sync node layout. Hex uses
        // trilinear at 19 child positions per refined parent; tet uses linear
        // at edge midpoints (numUniqueEdges new positions). After
        // rebuildDomainFromDevice, cstone redistributes nodes by SFC and dedupes
        // coincident-coord emissions; both element types are bit-exact-symmetric
        // across rank boundaries because all ranks compute the same interpolant
        // at a shared parent edge / face / vertex.
        std::vector<cstone::DeviceVector<RealType>> transferred(oldFields.size());
        for (size_t f = 0; f < oldFields.size(); ++f)
        {
            doTransferSolution(d_conn, oldFields[f], d_marks, refined,
                                startIdx, localCount, numNodes, transferred[f]);
        }

        // Save pre-sync coords before they're moved into cstone by rebuild.
        // We need them after rebuild to re-encode SFC keys with the NEW box that
        // cstone derives from the refined geometry. Encoding keys after rebuild
        // (rather than pinning the box) keeps cstone fully autonomous.
        size_t preNumNodes = refined.numNodes;
        cstone::DeviceVector<RealType> d_preX(preNumNodes), d_preY(preNumNodes), d_preZ(preNumNodes);
        cudaMemcpy(d_preX.data(), refined.d_x.data(), preNumNodes * sizeof(RealType), cudaMemcpyDeviceToDevice);
        cudaMemcpy(d_preY.data(), refined.d_y.data(), preNumNodes * sizeof(RealType), cudaMemcpyDeviceToDevice);
        cudaMemcpy(d_preZ.data(), refined.d_z.data(), preNumNodes * sizeof(RealType), cudaMemcpyDeviceToDevice);
        stats.transferTimeMs = phaseLap();

        rebuildDomainFromDevice(refined);
        amrOctree_.initialize(*domain_);
        stats.rebuildTimeMs = phaseLap();

        // Now encode pre-sync keys in the NEW box (the one cstone just derived).
        // Same physical coords -> same keys cstone computed during sync, so the
        // post-sync local nodes will find their pre-sync match by key lookup.
        const cstone::Box<RealType>& newBox = domain_->getBoundingBox();
        cstone::DeviceVector<KeyType> d_preKeys(preNumNodes);
        generateSfcKeys<KeyType, RealType>(
            d_preX.data(), d_preY.data(), d_preZ.data(),
            d_preKeys.data(), preNumNodes, newBox);

        // Sort (preKey, value) pairs by key for binary search post-sync.
        // Sort each field's value array against a fresh copy of the keys, since
        // sort_by_key permutes both. For the keys, capture once.
        cstone::DeviceVector<KeyType> d_sortedKeys;
        std::vector<cstone::DeviceVector<RealType>> d_sortedVals(oldFields.size());
        for (size_t f = 0; f < oldFields.size(); ++f)
        {
            cstone::DeviceVector<KeyType> kCopy = d_preKeys;
            d_sortedVals[f] = transferred[f];
            auto kb = thrust::device_pointer_cast(kCopy.data());
            auto ke = kb + preNumNodes;
            auto vb = thrust::device_pointer_cast(d_sortedVals[f].data());
            thrust::sort_by_key(thrust::device, kb, ke, vb);
            if (f == 0) d_sortedKeys = std::move(kCopy);
        }

        // Post-sync reordering: for each new local node, look up its SFC key
        // in the sorted pre-sync (key, value) table. Owner nodes find their
        // own emission; ghost nodes either find a peer's emission (boundary
        // nodes are emitted by all ranks touching them) or get 0 and are filled
        // by exchangeNodeHalo below.
        size_t newNumNodes = domain_->getNodeCount();
        const auto& d_newKeys = domain_->getLocalToGlobalSfcMap();
        for (size_t f = 0; f < oldFields.size(); ++f)
        {
            newFields[f]->resize(newNumNodes);

            const KeyType*  sortedKeys = thrust::raw_pointer_cast(d_sortedKeys.data());
            const RealType* sortedVals = thrust::raw_pointer_cast(d_sortedVals[f].data());
            const KeyType*  newKeys    = thrust::raw_pointer_cast(d_newKeys.data());
            RealType*       outVals    = thrust::raw_pointer_cast(newFields[f]->data());
            size_t          nSorted    = preNumNodes;

            thrust::for_each(thrust::device,
                              thrust::counting_iterator<size_t>(0),
                              thrust::counting_iterator<size_t>(newNumNodes),
                              [sortedKeys, sortedVals, newKeys, outVals, nSorted] __device__ (size_t i) {
                                  KeyType target = newKeys[i];
                                  size_t lo = 0, hi = nSorted;
                                  while (lo < hi) {
                                      size_t mid = (lo + hi) >> 1;
                                      if (sortedKeys[mid] < target) lo = mid + 1;
                                      else hi = mid;
                                  }
                                  outVals[i] = (lo < nSorted && sortedKeys[lo] == target)
                                                ? sortedVals[lo]
                                                : RealType(0);
                              });
            cudaDeviceSynchronize();

            // Fill any node not reached by the lookup (ghosts whose owner is
            // a peer that didn't emit them locally) via the node halo.
            domain_->exchangeNodeHalo(*newFields[f]);
        }
        stats.reorderTimeMs = phaseLap();   // includes halo fill (single loop)

        stats.elementsAfterLocal = domain_->localElementCount();
        stats.nodesAfterGlobal   = domain_->getNodeCount();
        MPI_Allreduce(&stats.elementsAfterLocal, &stats.elementsAfterGlobal, 1,
                       MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

        currentLevel_++;
        auto totalEnd     = std::chrono::high_resolution_clock::now();
        stats.totalTimeMs = std::chrono::duration<float, std::milli>(totalEnd - totalStart).count();
        return stats;
    }

    bool shouldContinue(const AmrStats& stats) const
    {
        return currentLevel_ < config_.maxLevels && stats.errorNorm > config_.errorTolerance &&
               stats.elementsRefined > 0;
    }

    int currentLevel() const { return currentLevel_; }
    Config& config() { return config_; }
    const Config& config() const { return config_; }

    static void printStats(const AmrStats& stats, int rank)
    {
        if (rank != 0) return;
        std::cout << "  AMR Level " << stats.level << ":\n";
        std::cout << "    Elements: " << stats.elementsBeforeGlobal << " -> " << stats.elementsAfterGlobal << "\n";
        std::cout << "    Nodes:    " << stats.nodesBeforeGlobal << " -> " << stats.nodesAfterGlobal << "\n";
        std::cout << "    Refined:  " << stats.elementsRefined << " elements\n";
        std::cout << "    Error:    " << std::scientific << stats.errorNorm << "\n";
        std::cout << "    Time:     " << std::fixed << stats.totalTimeMs << " ms\n";
        std::cout << "      mark:        " << std::fixed << stats.markTimeMs     << " ms\n";
        std::cout << "      refine:      " << std::fixed << stats.refineTimeMs   << " ms\n";
        std::cout << "      transfer:    " << std::fixed << stats.transferTimeMs << " ms\n";
        std::cout << "      rebuild:     " << std::fixed << stats.rebuildTimeMs  << " ms (cstone sync)\n";
        std::cout << "      reorder:     " << std::fixed << stats.reorderTimeMs  << " ms (key re-encode + lookup + halo fill)\n";
        std::cout << std::defaultfloat;
    }

private:
    // Rebuild ElementDomain from refined device data — no host round-trip.
    // Hands refined coords + connectivity straight to cstone via the
    // device-data ElementDomain constructor. cstone's sync redistributes
    // elements globally based on their SFC keys, establishing halos and
    // restoring partition consistency.
    //
    // Two paths, runtime-selected via env MARS_AMR_REUSE_DOMAIN:
    //
    //   default (env unset / 0):
    //     Reconstructs the ElementDomain from scratch. The new cstone::Domain
    //     has firstCall_=true and triggers focusTree_.converge(). At cube256/16
    //     this dominates AMR rebuild time (~80% of cstone sync = ~3 s at 977M).
    //
    //   MARS_AMR_REUSE_DOMAIN=1:
    //     EXPERIMENTAL / NOT BIT-EXACT. Reuses the existing cstone::Domain
    //     via resyncFromDevice. firstCall_ stays false. Drift observed in
    //     three iterations (5.9% to 78% L2 norm). cstone's internal state
    //     model is SPH-centric and cannot be retrofitted via downstream patches;
    //     proper fix requires upstream Domain::amrSync(). Kept here as a known
    //     follow-up; see docs/reference/10_cstone_domain_reuse_status.md.
    template<typename ResultType>
    void rebuildDomainFromDevice(ResultType& refined)
    {
        const char* reuseEnv = std::getenv("MARS_AMR_REUSE_DOMAIN");
        const bool reuseDomain = (reuseEnv != nullptr && std::string(reuseEnv) != "0"
                                  && domain_); // can only reuse if a Domain already exists

        typename Domain::DeviceCoordsTuple coords = std::make_tuple(
            std::move(refined.d_x), std::move(refined.d_y), std::move(refined.d_z));

        // Refined.d_conn is std::array<DeviceVector<KeyType>, NodesPerElement>.
        // Hex: pack-expand the 8 conn columns into the 8-tuple constructor arg.
        // Tet: same shape, just 4 columns.
        if constexpr (std::is_same_v<ElementTag, HexTag>)
        {
            typename Domain::DeviceConnectivityTuple conn = std::make_tuple(
                std::move(refined.d_conn[0]), std::move(refined.d_conn[1]),
                std::move(refined.d_conn[2]), std::move(refined.d_conn[3]),
                std::move(refined.d_conn[4]), std::move(refined.d_conn[5]),
                std::move(refined.d_conn[6]), std::move(refined.d_conn[7]));
            if (reuseDomain)
            {
                domain_->resyncFromDevice(std::move(coords), std::move(conn));
            }
            else
            {
                domain_ = std::make_unique<Domain>(std::move(coords), std::move(conn),
                                                    rank_, numRanks_, config_.bucketSize);
            }
        }
        else if constexpr (std::is_same_v<ElementTag, TetTag>)
        {
            typename Domain::DeviceConnectivityTuple conn = std::make_tuple(
                std::move(refined.d_conn[0]), std::move(refined.d_conn[1]),
                std::move(refined.d_conn[2]), std::move(refined.d_conn[3]));
            if (reuseDomain)
            {
                domain_->resyncFromDevice(std::move(coords), std::move(conn));
            }
            else
            {
                domain_ = std::make_unique<Domain>(std::move(coords), std::move(conn),
                                                    rank_, numRanks_, config_.bucketSize);
            }
        }

        // Lazy-init order matching mars_cvfem_graph: ownership first, then conn,
        // coords, then AMR octree state. resyncFromDevice already cleared lazy
        // state, so first access triggers re-init.
        (void)domain_->getNodeOwnershipMap();
        cudaDeviceSynchronize();
        (void)domain_->getElementToNodeConnectivity();
        domain_->cacheNodeCoordinates();
    }

    // Element-type-dispatched refinement helper. Templated on the conn tuple
    // so the 8-conn body is only parsed when `ConnTupleT` is a hex 8-tuple
    // (and vice versa for tet). Without this layer the if constexpr branches
    // inside adaptMeshMultiField are non-dependent at the point of parsing
    // and nvcc complains about std::get<7> on a 4-tuple even for TetTag.
    template<typename ConnTupleT, typename CoordVecT, typename MarkVecT>
    typename OctreeAlignedRefine<KeyType, RealType, ElementTag>::Result
    doRefineLocal(const ConnTupleT& d_conn,
                  const CoordVecT& d_x, const CoordVecT& d_y, const CoordVecT& d_z,
                  const MarkVecT& d_marks,
                  size_t startIdx, size_t localCount, size_t numNodes)
    {
        using Refine = OctreeAlignedRefine<KeyType, RealType, ElementTag>;
        if constexpr (std::tuple_size_v<ConnTupleT> == 8)
        {
            return Refine::refineLocal(
                std::get<0>(d_conn).data() + startIdx, std::get<1>(d_conn).data() + startIdx,
                std::get<2>(d_conn).data() + startIdx, std::get<3>(d_conn).data() + startIdx,
                std::get<4>(d_conn).data() + startIdx, std::get<5>(d_conn).data() + startIdx,
                std::get<6>(d_conn).data() + startIdx, std::get<7>(d_conn).data() + startIdx,
                d_x.data(), d_y.data(), d_z.data(),
                d_marks.data() + startIdx, localCount, numNodes, config_.blockSize);
        }
        else if constexpr (std::tuple_size_v<ConnTupleT> == 4)
        {
            return Refine::refineLocal(
                std::get<0>(d_conn).data() + startIdx, std::get<1>(d_conn).data() + startIdx,
                std::get<2>(d_conn).data() + startIdx, std::get<3>(d_conn).data() + startIdx,
                d_x.data(), d_y.data(), d_z.data(),
                d_marks.data() + startIdx, localCount, numNodes, config_.blockSize);
        }
    }

    template<typename ConnTupleT, typename MarkVecT, typename RefinedT, typename OutVecT>
    void doTransferSolution(const ConnTupleT& d_conn,
                             const RealType* oldField,
                             const MarkVecT& d_marks,
                             const RefinedT& refined,
                             size_t startIdx, size_t localCount, size_t numNodes,
                             OutVecT& transferred)
    {
        using Refine = OctreeAlignedRefine<KeyType, RealType, ElementTag>;
        if constexpr (std::tuple_size_v<ConnTupleT> == 8)
        {
            Refine::transferSolution(
                std::get<0>(d_conn).data() + startIdx, std::get<1>(d_conn).data() + startIdx,
                std::get<2>(d_conn).data() + startIdx, std::get<3>(d_conn).data() + startIdx,
                std::get<4>(d_conn).data() + startIdx, std::get<5>(d_conn).data() + startIdx,
                std::get<6>(d_conn).data() + startIdx, std::get<7>(d_conn).data() + startIdx,
                oldField, d_marks.data() + startIdx,
                localCount, numNodes, transferred, config_.blockSize);
        }
        else if constexpr (std::tuple_size_v<ConnTupleT> == 4)
        {
            Refine::transferSolution(
                std::get<0>(d_conn).data() + startIdx, std::get<1>(d_conn).data() + startIdx,
                std::get<2>(d_conn).data() + startIdx, std::get<3>(d_conn).data() + startIdx,
                oldField, refined, localCount, transferred, config_.blockSize);
        }
    }

    Config config_;
    int currentLevel_ = 0;
    int rank_         = 0;
    int numRanks_     = 1;
    InitTimings initTimings_;

    DomainPtr domain_;
    AmrOctree<KeyType, RealType> amrOctree_;
    OctreeNativeAmr<KeyType, RealType> octreeNativeAmr_;
};

} // namespace amr
} // namespace mars
