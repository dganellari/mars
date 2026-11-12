#pragma once

#include <memory>
#include <optional>
#include <chrono>
#include <vector>
#include <mpi.h>

#include "backend/distributed/unstructured/domain.hpp"
#include "mars_amr_octree.hpp"
#include "mars_amr_keyed_position.hpp"
#include "mars_amr_tet_audit.hpp"
#include "mars_amr_octree_native.hpp"
#include "mars_amr_hex_refine.hpp"
#include "mars_amr_tet_refine.hpp"
#include "mars_amr_octree_refine.hpp"
#include "mars_amr_tet_bisection.hpp"
#include "mars_amr_field_migration.hpp"
#include "mars_amr_solution_transfer.hpp"
#include "mars_amr_error_indicator.hpp"
#include "backend/distributed/unstructured/utils/mars_gpu_reduce.hpp"

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

// The keyed point of every local node; counts the nodes whose point does not encode back to their key
template<typename KeyType, typename RealType>
__global__ void keyedPositionsKernel(const KeyType* keys, size_t n, KeyedPosition<KeyType, RealType> pos, RealType* x,
                                     RealType* y, RealType* z, unsigned long long* mismatches)
{
    size_t i = blockIdx.x * size_t(blockDim.x) + threadIdx.x;
    if (i >= n) return;
    pos(keys[i], x[i], y[i], z[i]);
    if (pos.encode(x[i], y[i], z[i]) != keys[i]) atomicAdd(mismatches, 1ull);
}

// Every emitted node to the keyed point of the cell it lies in
template<typename KeyType, typename RealType>
__global__ void snapToKeysKernel(size_t n, KeyedPosition<KeyType, RealType> pos, RealType* x, RealType* y, RealType* z)
{
    size_t i = blockIdx.x * size_t(blockDim.x) + threadIdx.x;
    if (i >= n) return;
    pos(pos.encode(x[i], y[i], z[i]), x[i], y[i], z[i]);
}

enum class MarkingStrategy
{
    Doerfler,     // Error-driven Doerfler marking + 1-irregular enforcement
    OctreeNative  // Cstone octree leaf-level split/merge decisions drive element marks
};

enum class TetRefinement
{
    Red,
    ConformingBisection
};

// Refiner traits: selects the right refiner + transfer for each element type
template<typename ElementTag, typename KeyType, typename RealType>
struct RefinerTraits;

template<typename KeyType, typename RealType>
struct RefinerTraits<HexTag, KeyType, RealType>
{
    using Refiner = HexRefiner<KeyType, RealType>;
    using Result  = typename Refiner::Result;

    static Result refine(const ElementDomain<HexTag, RealType, KeyType, cstone::execution::Gpu>& domain,
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
    static Result refineLocal(const ElementDomain<HexTag, RealType, KeyType, cstone::execution::Gpu>& domain,
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
        const ElementDomain<HexTag, RealType, KeyType, cstone::execution::Gpu>& domain,
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

    static Result refine(const ElementDomain<TetTag, RealType, KeyType, cstone::execution::Gpu>& domain,
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

    static Result refineLocal(const ElementDomain<TetTag, RealType, KeyType, cstone::execution::Gpu>& domain,
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
        const ElementDomain<TetTag, RealType, KeyType, cstone::execution::Gpu>& domain,
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
cstone::Box<RealType> computeBoundingBox(cstone::DeviceVector<RealType>& d_x,
                                             cstone::DeviceVector<RealType>& d_y,
                                             cstone::DeviceVector<RealType>& d_z,
                                             RealType padding = 0.05)
{
    RealType lxmin = mars::reduce(d_x.begin(), d_x.end(), RealType(1e30), thrust::minimum<RealType>());
    RealType lxmax = mars::reduce(d_x.begin(), d_x.end(), RealType(-1e30), thrust::maximum<RealType>());
    RealType lymin = mars::reduce(d_y.begin(), d_y.end(), RealType(1e30), thrust::minimum<RealType>());
    RealType lymax = mars::reduce(d_y.begin(), d_y.end(), RealType(-1e30), thrust::maximum<RealType>());
    RealType lzmin = mars::reduce(d_z.begin(), d_z.end(), RealType(1e30), thrust::minimum<RealType>());
    RealType lzmax = mars::reduce(d_z.begin(), d_z.end(), RealType(-1e30), thrust::maximum<RealType>());

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
    using Domain    = ElementDomain<ElementTag, RealType, KeyType, cstone::execution::Gpu>;
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
        TetRefinement tet_refinement   = TetRefinement::Red;
        bool stage_refinement_messages = false;
        // Opt-in: keep the Exodus side sets on-device in the Domain so FEM/FSI/viz
        // can query local boundary node lists without re-reading the mesh file.
        bool storeSideSets       = false;
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
                                            bucketSizeFocus,
                                            periodicAxesMask,
                                            periodicBoxLo, periodicBoxHi,
                                            config_.storeSideSets);
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
        if (config_.tet_refinement == TetRefinement::ConformingBisection)
        {
            const char* reuse = std::getenv("MARS_AMR_REUSE_DOMAIN");
            int invalid       = !std::is_same_v<ElementTag, TetTag> || domain_->numBlocks() != 1 ||
                          domain_->hasBoundaryInfo() || domain_->getPeriodicAxesMask() != 0 ||
                          (reuse && std::string(reuse) != "0") || (numRanks_ > 1 && !domain_->sfcOwnership());
            MPI_Allreduce(MPI_IN_PLACE, &invalid, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            if (invalid)
                throw std::runtime_error("Conforming bisection requires single-block, nonperiodic tetrahedra with face "
                                         "labels, SFC ownership and fresh domain rebuilds.");
        }
        uint64_t fields = oldFields.size(), min_fields = fields, max_fields = fields;
        MPI_Allreduce(MPI_IN_PLACE, &min_fields, 1, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &max_fields, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
        int invalid = min_fields != max_fields || oldFields.size() != newFields.size() ||
                      (!d_errorPerElement && domain_->getElementCount()) ||
                      (!oldFields.empty() && numRanks_ > 1 && !domain_->sfcOwnership());
        for (const RealType* field : oldFields)
            if (!field && domain_->getNodeCount()) invalid = 1;
        for (const DevVector* field : newFields)
            if (!field) invalid = 1;
        MPI_Allreduce(MPI_IN_PLACE, &invalid, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (invalid)
            throw std::runtime_error(
                "AMR needs valid, equally sized field lists on every rank and SFC node ownership.");
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
            // Each global leaf is owned by one rank; only the maximum error is reduced globally.
            // Mesh balance is enforced below, independently of the spatial marking policy.
            d_marks = octreeNativeAmr_.markFromOctree(d_errorPerElement, *domain_);
            std::cerr << "ADAPT r" << rank_ << " markFromOctree returned" << std::endl; std::cerr.flush();
        }
        else
        {
            d_marks = amrOctree_.markForRefinement(d_errorPerElement, numElements);
        }

        // Hanging-node constraints require mesh balance for either marking policy.
        // Coarsening is not implemented: an element marked for it (2) stays as it is. The refiners split every
        // element with a non-zero mark.
        thrust::transform(thrust::device, thrust::device_pointer_cast(d_marks.data()),
                          thrust::device_pointer_cast(d_marks.data() + numElements),
                          thrust::device_pointer_cast(d_marks.data()),
                          [] __host__ __device__(uint8_t m) -> uint8_t { return m == 1 ? 1 : 0; });

        if (startIdx > 0)
            cudaMemset(d_marks.data(), 0, startIdx * sizeof(uint8_t));
        if (endIdx < numElements)
            cudaMemset(d_marks.data() + endIdx, 0, (numElements - endIdx) * sizeof(uint8_t));

        auto exchangeMarks = [&]()
        {
            if (numRanks_ <= 1) return;
            // Cornerstone provides halo exchange for int, not uint8_t.
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

        if (config_.tet_refinement != TetRefinement::ConformingBisection)
        {
            // A mark set near a rank boundary changes what the neighbouring rank must mark, so sweep with
            // fresh halo marks until no rank changes a mark: the refined mesh is then 2:1 balanced.
            constexpr int maxSweeps = 64;
            for (int sweep = 0;; ++sweep)
            {
                exchangeMarks();
                int changed = amrOctree_.enforce1Irregular(d_marks, *domain_) ? 1 : 0;
                MPI_Allreduce(MPI_IN_PLACE, &changed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                if (!changed) break;
                if (sweep + 1 == maxSweeps)
                {
                    if (rank_ == 0)
                        std::fprintf(stderr, "AMR: 2:1 balance did not settle in %d sweeps\n", maxSweeps);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
            }
        }

        if (startIdx > 0)
            cudaMemset(d_marks.data(), 0, startIdx * sizeof(uint8_t));
        if (endIdx < numElements)
            cudaMemset(d_marks.data() + endIdx, 0, (numElements - endIdx) * sizeof(uint8_t));

        if (config_.tet_refinement != TetRefinement::ConformingBisection) exchangeMarks();
        stats.markTimeMs = phaseLap();

        // Ghost marks are copies; only owned parents contribute to the refinement count.
        auto marks_b = thrust::device_pointer_cast(d_marks.data() + startIdx);
        auto marks_e = thrust::device_pointer_cast(d_marks.data() + endIdx);
        bool conforming = config_.tet_refinement == TetRefinement::ConformingBisection;
        uint64_t local_refined = mars::countIf(marks_b, marks_e,
            [conforming] __device__(uint8_t m) { return conforming ? m == 1 : m > 0; });
        uint64_t global_refined = 0;
        MPI_Allreduce(&local_refined, &global_refined, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
        stats.elementsRefined = global_refined;

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

        // Non-conforming refinement starts from the keyed points of the nodes (KeyedPosition) and keeps the box.
        const cstone::Box<RealType> keptBox = domain_->getBoundingBox();
        DevVector keyedX, keyedY, keyedZ;
        if (!conforming) keyedPositions(keptBox, keyedX, keyedY, keyedZ);
        const auto& d_x = conforming ? domain_->getNodeX() : keyedX;
        const auto& d_y = conforming ? domain_->getNodeY() : keyedY;
        const auto& d_z = conforming ? domain_->getNodeZ() : keyedZ;
        const auto& d_conn = domain_->getElementToNodeConnectivity();

        std::optional<TetAuditSnapshot<KeyType>> audit_before, audit_emitted;
        if constexpr (std::is_same_v<ElementTag, TetTag>)
            if (!conforming && tet_audit_enabled())
            {
                audit_before.emplace(audit_tets<KeyType, RealType>("refine-input",
                    {std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
                     std::get<2>(d_conn).data(), std::get<3>(d_conn).data()}, startIdx, endIdx, keptBox,
                    d_x.data(), d_y.data(), d_z.data(), numNodes, false));
                audit_require(!audit_before->totals.invalid, "invalid parent geometry", *audit_before);
            }

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
        if constexpr (std::is_same_v<ElementTag, TetTag>)
            if (audit_before)
            {
                auto raw = audit_tets<KeyType, RealType>("refine-raw",
                    {refined.d_conn[0].data(), refined.d_conn[1].data(), refined.d_conn[2].data(),
                     refined.d_conn[3].data()}, 0, refined.numElements, keptBox,
                    refined.d_x.data(), refined.d_y.data(), refined.d_z.data(), refined.numNodes, false);
                const uint64_t expected = audit_before->totals.count + 7 * global_refined;
                if (rank_ == 0)
                    std::printf("[amr-mesh-audit] level=%d old=%llu refined=%llu expected=%llu emitted=%llu\n",
                        currentLevel_, (unsigned long long)audit_before->totals.count,
                        (unsigned long long)global_refined, (unsigned long long)expected,
                        (unsigned long long)raw.totals.count);
                audit_require(refined.numMarked == local_refined && raw.totals.count == expected,
                              "owned refinement count mismatch", raw, &*audit_before);
                audit_red_children(*audit_before, raw, d_marks.data() + startIdx,
                                   refined.d_markedPrefix.empty() ? nullptr : refined.d_markedPrefix.data());
                double tolerance = audit_volume_tolerance(audit_before->totals, raw.totals);
                audit_require(!raw.totals.invalid &&
                    std::abs(raw.totals.volume - audit_before->totals.volume) <= tolerance &&
                    std::abs(raw.totals.absolute_volume - audit_before->totals.absolute_volume) <= tolerance,
                    "raw refinement changed volume or produced invalid cells", raw, &*audit_before);
                audit_emitted.emplace(std::move(raw));
            }
        if (config_.tet_refinement == TetRefinement::ConformingBisection)
            MPI_Allreduce(&refined.numMarked, &stats.elementsRefined, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
        else if (refined.numNodes > 0)
        {
            snapToKeysKernel<KeyType, RealType><<<unsigned((refined.numNodes + 255) / 256), 256>>>(
                refined.numNodes, KeyedPosition<KeyType, RealType>{keptBox}, refined.d_x.data(), refined.d_y.data(),
                refined.d_z.data());
            cudaCheckError();
        }
        if constexpr (std::is_same_v<ElementTag, TetTag>)
            if (audit_emitted)
            {
                auto snapped = audit_tets<KeyType, RealType>("refine-snapped",
                    {refined.d_conn[0].data(), refined.d_conn[1].data(), refined.d_conn[2].data(),
                     refined.d_conn[3].data()}, 0, refined.numElements, keptBox,
                    refined.d_x.data(), refined.d_y.data(), refined.d_z.data(), refined.numNodes);
                audit_require(audit_same_cells(audit_emitted->totals, snapped.totals),
                              "snapping changed node identity or collapsed a tet", snapped, &*audit_emitted);
                if (rank_ == 0)
                    std::printf("[amr-mesh-audit] snapping_volume_delta=%.17g\n",
                                snapped.totals.volume - audit_emitted->totals.volume);
                audit_emitted = std::move(snapped);
                audit_before.reset();
            }
        stats.refineTimeMs = phaseLap();

        // Interpolate from synchronized parent values, then migrate to the rebuilt
        // ownership. New copies of existing hanging nodes must not replace old values.
        std::vector<cstone::DeviceVector<RealType>> transferred(oldFields.size());
        for (size_t f = 0; f < oldFields.size(); ++f)
        {
            DevVector synchronized_field(numNodes);
            if (numNodes)
                cudaCheckCall(cudaMemcpy(synchronized_field.data(), oldFields[f], numNodes * sizeof(RealType),
                                         cudaMemcpyDeviceToDevice));
            domain_->exchangeNodeHalo(synchronized_field);
            doTransferSolution(d_conn, synchronized_field.data(), d_marks, refined,
                               startIdx, localCount, numNodes, transferred[f]);
        }

        // Save pre-sync coords before they're moved into cstone by rebuild.
        // We need them after rebuild to re-encode SFC keys with the NEW box that
        // cstone derives from the refined geometry. LEPP checks that encoding
        // before rebuild, using the constructor's exact box policy.
        size_t preNumNodes = refined.numNodes;
        std::vector<std::string> boundary_names;
        LeppDeviceExecution::Vector<LeppFaceRecord<KeyType>> saved_faces;
        if constexpr (std::is_same_v<ElementTag, TetTag>)
            if (conforming)
            {
                boundary_names = domain_->sideSetNames();
                saved_faces    = save_lepp_boundary<KeyType>(refined);
            }

        cstone::DeviceVector<uint8_t> migration_sources;
        if (!oldFields.empty())
        {
            // Deep halo nodes may have stale fields. Their old owner supplies
            // the value; every locally emitted midpoint has current endpoints.
            migration_sources = cstone::DeviceVector<uint8_t>(preNumNodes, uint8_t(1));
            if (numNodes)
                cudaCheckCall(cudaMemcpy(migration_sources.data(), domain_->getNodeOwnershipMap().data(),
                                         numNodes * sizeof(uint8_t), cudaMemcpyDeviceToDevice));
        }
        cstone::DeviceVector<RealType> d_preX(preNumNodes), d_preY(preNumNodes), d_preZ(preNumNodes);
        if (preNumNodes)
            cudaCheckCall(cudaMemcpy(d_preX.data(), refined.d_x.data(), preNumNodes * sizeof(RealType),
                                     cudaMemcpyDeviceToDevice));
        if (preNumNodes)
            cudaCheckCall(cudaMemcpy(d_preY.data(), refined.d_y.data(), preNumNodes * sizeof(RealType),
                                     cudaMemcpyDeviceToDevice));
        if (preNumNodes)
            cudaCheckCall(cudaMemcpy(d_preZ.data(), refined.d_z.data(), preNumNodes * sizeof(RealType),
                                     cudaMemcpyDeviceToDevice));
        stats.transferTimeMs = phaseLap();

        cstone::DeviceVector<KeyType> d_preKeys(preNumNodes);
        if (conforming)
        {
            auto rebuilt_box =
                device_mesh_box(d_preX.data(), d_preY.data(), d_preZ.data(), preNumNodes, MPI_COMM_WORLD);
            generateSfcKeys<KeyType, RealType>(d_preX.data(), d_preY.data(), d_preZ.data(), d_preKeys.data(),
                                               preNumNodes, rebuilt_box);
            auto node_owner    = domain_->sfcNodeOwner();
            node_owner.meshBox = rebuilt_box;
            auto route = AmrNodeOwner<KeyType, decltype(node_owner)>{node_owner, node_owner.numRanks};
            validate_lepp_node_keys<LeppDeviceExecution>(d_preKeys.data(), d_preX.data(), d_preY.data(), d_preZ.data(),
                                                         preNumNodes, route, MPI_COMM_WORLD,
                                                         config_.stage_refinement_messages);
        }
        rebuildDomainFromDevice(refined, conforming ? nullptr : &keptBox);
        if constexpr (std::is_same_v<ElementTag, TetTag>)
            if (audit_emitted)
            {
                const auto& conn = domain_->getElementToNodeConnectivity();
                auto rebuilt = audit_tets<KeyType, RealType>("connectivity-rebuilt",
                    {std::get<0>(conn).data(), std::get<1>(conn).data(), std::get<2>(conn).data(),
                     std::get<3>(conn).data()}, domain_->startIndex(), domain_->endIndex(), keptBox,
                    domain_->getNodeX().data(), domain_->getNodeY().data(), domain_->getNodeZ().data(),
                    domain_->getNodeCount(), true, domain_->getLocalToGlobalSfcMap().data());
                audit_compare(*audit_emitted, rebuilt);
                audit_emitted.reset();
            }
        amrOctree_.initialize(*domain_);
        stats.rebuildTimeMs = phaseLap();

        // Now encode pre-sync keys in the NEW box (the one cstone just derived).
        // Same physical coords -> same keys cstone computed during sync, so the
        // post-sync local nodes will find their pre-sync match by key lookup.
        const cstone::Box<RealType>& newBox = domain_->getBoundingBox();
        if (!conforming)
            generateSfcKeys<KeyType, RealType>(d_preX.data(), d_preY.data(), d_preZ.data(), d_preKeys.data(),
                                               preNumNodes, newBox);
        if constexpr (std::is_same_v<ElementTag, TetTag>)
            if (conforming)
                restore_lepp_boundary<KeyType>(*domain_, boundary_names, saved_faces, d_preKeys.data(),
                                               config_.stage_refinement_messages);

        // Repartitioning can move a node beyond its old owner's halo. Route
        // authoritative values to the new owner before filling ghost copies.
        const size_t newNumNodes = domain_->getNodeCount();
        const auto& d_newKeys = domain_->getLocalToGlobalSfcMap();
        std::unique_ptr<AmrFieldMigration<LeppDeviceExecution, KeyType, RealType>> migration;
        if (!oldFields.empty())
        {
            require_collectively(numRanks_ == 1 || domain_->sfcOwnership(),
                                 "AMR migration requires SFC ownership", MPI_COMM_WORLD);
            auto owner = domain_->sfcNodeOwner();
            migration = std::make_unique<AmrFieldMigration<LeppDeviceExecution, KeyType, RealType>>(
                d_preKeys.data(), preNumNodes, migration_sources.data(), d_newKeys.data(),
                domain_->getNodeOwnershipMap().data(), newNumNodes,
                AmrNodeOwner<KeyType, decltype(owner)>{owner, numRanks_}, MPI_COMM_WORLD,
                config_.stage_refinement_messages, conforming ? SIZE_MAX : numNodes);
        }
        for (size_t f = 0; f < oldFields.size(); ++f)
        {
            migration->apply(transferred[f].data(), transferred[f].size(), *newFields[f]);
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
    // cornerstone focus-tree bucket size of every domain the manager builds
    static constexpr unsigned bucketSizeFocus = 8;
    // The keyed point (KeyedPosition) of every local node in box, the box its keys were encoded with
    void keyedPositions(const cstone::Box<RealType>& box, DevVector& x, DevVector& y, DevVector& z) const
    {
        const size_t n     = domain_->getNodeCount();
        const auto& keys   = domain_->getLocalToGlobalSfcMap();
        unsigned long long mismatches = 0;
        x.resize(n);
        y.resize(n);
        z.resize(n);
        thrust::device_vector<unsigned long long> d_mismatches(1, 0ull);
        if (n > 0)
            keyedPositionsKernel<KeyType, RealType><<<unsigned((n + 255) / 256), 256>>>(
                keys.data(), n, KeyedPosition<KeyType, RealType>{box}, x.data(), y.data(), z.data(),
                thrust::raw_pointer_cast(d_mismatches.data()));
        cudaCheckError();
        cudaCheckCall(cudaMemcpy(&mismatches, thrust::raw_pointer_cast(d_mismatches.data()), sizeof(mismatches), cudaMemcpyDeviceToHost));
        MPI_Allreduce(MPI_IN_PLACE, &mismatches, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        if (mismatches > 0)
        {
            if (rank_ == 0)
                std::fprintf(stderr, "AMR: %llu keyed node positions do not encode back to their keys\n", mismatches);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // Rebuild from device children through Cornerstone's partition, focused-tree and halo machinery.
    // The opt-in reuse experiment remains unvalidated; fixed-box refinement always constructs a new domain.
    // Refresh all tree views after this call because sync can replace their backing buffers.
    template<typename ResultType>
    void rebuildDomainFromDevice(ResultType& refined, const cstone::Box<RealType>* fixedBox)
    {
        const char* reuseEnv = std::getenv("MARS_AMR_REUSE_DOMAIN");
        const bool reuseDomain = (reuseEnv != nullptr && std::string(reuseEnv) != "0"
                                  && domain_ && !fixedBox); // can only reuse if a Domain already exists

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
                domain_ = std::make_unique<Domain>(std::move(coords), std::move(conn), rank_, numRanks_,
                                                   config_.bucketSize, bucketSizeFocus, fixedBox);
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
                domain_ = std::make_unique<Domain>(std::move(coords), std::move(conn), rank_, numRanks_,
                                                   config_.bucketSize, bucketSizeFocus, fixedBox);
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
            if (config_.tet_refinement == TetRefinement::ConformingBisection)
                return TetBisectionRefiner<KeyType, RealType>::refine(*domain_, d_marks.data(),
                                                                      config_.stage_refinement_messages);
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
