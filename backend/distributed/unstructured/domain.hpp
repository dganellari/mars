#pragma once

#include <tuple>
#include <utility>

namespace cstone
{
using std::get;
}

// Cornerstone includes
#include "cstone/domain/domain.hpp"
#include "cstone/cuda/cuda_utils.hpp"
#include "cstone/sfc/sfc.hpp"
#include "cstone/domain/assignment.hpp"

#include "domain_cuda_impl.hpp"

// stl includes
// #include <adios2.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/extrema.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/copy.h>
#include <thrust/sequence.h>
#include <thrust/reduce.h>
#include <thrust/scatter.h>
#include <thrust/scan.h>
#include <thrust/fill.h>
#include <thrust/iterator/constant_iterator.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <limits>
#include <mpi.h>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

// Mars includes
#include "mars.hpp"
#include "mars_globals.hpp"
#include "mars_domain_utils.hpp"
// #include "mars_read_mesh_adios2.hpp"
#include "mars_read_mesh_binary.hpp"
#include "mars_read_mfem_mesh.hpp"
#include "mars_read_exodus_mesh.hpp"

namespace mars
{

// Element type tags
struct TetTag
{
    static constexpr int NodesPerElement = 4;
    static constexpr const char* Name    = "Tetrahedron";
};

struct HexTag
{
    static constexpr int NodesPerElement = 8;
    static constexpr const char* Name    = "Hexahedron";
};

struct TriTag
{
    static constexpr int NodesPerElement = 3;
    static constexpr const char* Name    = "Triangle";
};

struct QuadTag
{
    static constexpr int NodesPerElement = 4;
    static constexpr const char* Name    = "Quadrilateral";
};

// Forward declarations of CUDA kernels
template<typename ElementTag, typename KeyType, typename RealType>
__global__ void findRepresentativeNodesKernel(const KeyType* indices0,
                                              const KeyType* indices1,
                                              const KeyType* indices2,
                                              const KeyType* indices3,
                                              const KeyType* indices4,
                                              const KeyType* indices5,
                                              const KeyType* indices6,
                                              const KeyType* indices7,
                                              const KeyType* sfcCodes,
                                              KeyType* elemToNodeMap,
                                              int numElements);

template<typename ElementTag, typename KeyType, typename RealType>
__global__ void extractRepCoordinatesKernel(const RealType* x,
                                            const RealType* y,
                                            const RealType* z,
                                            const RealType* h,
                                            const KeyType* elemToNodeMap,
                                            RealType* elemX,
                                            RealType* elemY,
                                            RealType* elemZ,
                                            RealType* elemH,
                                            int numElements);

template<typename ElementTag, typename KeyType, typename RealType>
__global__ void computeCharacteristicSizesKernel(const RealType* x,
                                                 const RealType* y,
                                                 const RealType* z,
                                                 const KeyType* indices0,
                                                 const KeyType* indices1,
                                                 const KeyType* indices2,
                                                 const KeyType* indices3,
                                                 const KeyType* indices4,
                                                 const KeyType* indices5,
                                                 const KeyType* indices6,
                                                 const KeyType* indices7,
                                                 int* nodeTetCount,
                                                 RealType* h,
                                                 int numElements);

template<typename KeyType, typename RealType>
__global__ void finalizeCharacteristicSizesKernel(RealType* h, int* nodeTetCount, int numNodes);

// Forward declarations of CUDA kernels
template<typename RealType>
__global__ void
transformCharacteristicSizesKernel(RealType* d_h, size_t size, RealType meshFactor, RealType minH, RealType maxH);

template<typename RealType>
__global__ void fillCharacteristicSizesKernel(RealType* d_h, size_t size, RealType value);

// After other forward declarations
template<typename KeyType, typename RealType>
void generateSfcKeys(const RealType* x,
                     const RealType* y,
                     const RealType* z,
                     KeyType* keys,
                     size_t numKeys,
                     const cstone::Box<RealType>& box);

// GPU function to rebuild connectivity with SFC keys
// This will be implemented in domain.cu
// template<typename KeyType, typename SfcConnTple, typename ConnTuple>
// void rebuildElementConnectivity(SfcConnTuple& d_conn_keys_, ConnTuple& d_conn_, size_t newElementCount);

template<typename ElementTag, typename KeyType, typename RealType>
__global__ void buildSfcConnectivity(const KeyType* indices0,
                                     const KeyType* indices1,
                                     const KeyType* indices2,
                                     const KeyType* indices3,
                                     const KeyType* indices4,
                                     const KeyType* indices5,
                                     const KeyType* indices6,
                                     const KeyType* indices7,
                                     const KeyType* nodeSfcCodes,
                                     KeyType* conn_key0,
                                     KeyType* conn_key1,
                                     KeyType* conn_key2,
                                     KeyType* conn_key3,
                                     KeyType* conn_key4,
                                     KeyType* conn_key5,
                                     KeyType* conn_key6,
                                     KeyType* conn_key7,
                                     int numElements);

// Forward declarations for kernels that handle original coordinates
template<typename ElementTag, typename KeyType, typename RealType>
__global__ void extractElementNodeCoordsKernel(const RealType* x,
                                               const RealType* y,
                                               const RealType* z,
                                               const KeyType* idx0,
                                               const KeyType* idx1,
                                               const KeyType* idx2,
                                               const KeyType* idx3,
                                               const KeyType* idx4,
                                               const KeyType* idx5,
                                               const KeyType* idx6,
                                               const KeyType* idx7,
                                               RealType* elemOrigX0, RealType* elemOrigY0, RealType* elemOrigZ0,
                                               RealType* elemOrigX1, RealType* elemOrigY1, RealType* elemOrigZ1,
                                               RealType* elemOrigX2, RealType* elemOrigY2, RealType* elemOrigZ2,
                                               RealType* elemOrigX3, RealType* elemOrigY3, RealType* elemOrigZ3,
                                               RealType* elemOrigX4, RealType* elemOrigY4, RealType* elemOrigZ4,
                                               RealType* elemOrigX5, RealType* elemOrigY5, RealType* elemOrigZ5,
                                               RealType* elemOrigX6, RealType* elemOrigY6, RealType* elemOrigZ6,
                                               RealType* elemOrigX7, RealType* elemOrigY7, RealType* elemOrigZ7,
                                               int numElements);

template<typename ElementTag, typename KeyType, typename RealType>
__global__ void rebuildNodeCoordsFromElementsKernel(
    const KeyType* connKey0, const KeyType* connKey1, const KeyType* connKey2, const KeyType* connKey3,
    const KeyType* connKey4, const KeyType* connKey5, const KeyType* connKey6, const KeyType* connKey7,
    const RealType* elemOrigX0, const RealType* elemOrigY0, const RealType* elemOrigZ0,
    const RealType* elemOrigX1, const RealType* elemOrigY1, const RealType* elemOrigZ1,
    const RealType* elemOrigX2, const RealType* elemOrigY2, const RealType* elemOrigZ2,
    const RealType* elemOrigX3, const RealType* elemOrigY3, const RealType* elemOrigZ3,
    const RealType* elemOrigX4, const RealType* elemOrigY4, const RealType* elemOrigZ4,
    const RealType* elemOrigX5, const RealType* elemOrigY5, const RealType* elemOrigZ5,
    const RealType* elemOrigX6, const RealType* elemOrigY6, const RealType* elemOrigZ6,
    const RealType* elemOrigX7, const RealType* elemOrigY7, const RealType* elemOrigZ7,
    RealType* nodeX, RealType* nodeY, RealType* nodeZ,
    int numElements);

template<typename KeyType>
__global__ void decodeSfcToIntegersKernel(const KeyType* keys, unsigned* x, unsigned* y, unsigned* z, size_t numKeys);

template<typename RealType, typename KeyType>
__global__ void integerToPhysicalKernel(const unsigned* ix,
                                        const unsigned* iy,
                                        const unsigned* iz,
                                        RealType* x,
                                        RealType* y,
                                        RealType* z,
                                        size_t numKeys,
                                        cstone::Box<RealType> box);

template<typename KeyType>
__global__ void convertSfcToNodeIndicesKernel(
    const KeyType* sfcIndices, KeyType* nodeIndices, const KeyType* particleKeys, size_t numElements, size_t numNodes);

// Forward declaration for the shared SFC decoding function
template<typename KeyType, typename RealType>
MARS_HOST_DEVICE std::tuple<RealType, RealType, RealType> decodeSfcToPhysical(KeyType sfcKey,
                                                                              const cstone::Box<RealType>& box);

// Better: use std::array for compile-time indexing
template<typename KeyType, size_t NodesPerElement>
struct ConnPtrs
{
    const KeyType* ptrs[NodesPerElement];
};

template<typename KeyType, size_t NodesPerElement>
__global__ void
flattenConnectivityKernel(ConnPtrs<KeyType, NodesPerElement> conn, KeyType* flat_keys, size_t numElements);

template<typename KeyType>
__global__ void mapSfcToLocalIdKernel(
    const KeyType* sfc_conn, KeyType* local_conn, const KeyType* sorted_sfc, size_t num_elements, size_t num_nodes);

__global__ void initNodeOwnershipKernel(uint8_t* nodeOwnership, size_t nodeCount);

template<typename KeyType, typename RealType>
__global__ void decodeAllNodesKernel(
    const KeyType* sfcKeys, RealType* x, RealType* y, RealType* z, size_t numNodes, cstone::Box<RealType> box);

// Template struct to select the correct vector type based on accelerator tag
template<typename T, typename AcceleratorTag>
struct VectorSelector
{
    // Default to host vector for any tag
    using type = std::vector<T>;
};

// Specialization for GPU tag - use cornerstone device vector
template<typename T>
struct VectorSelector<T, cstone::GpuTag>
{
    using type = cstone::DeviceVector<T>;
};

// Specialization for host tag - use cornerstone host vector
template<typename T>
struct VectorSelector<T, cstone::CpuTag>
{
    using type = std::vector<T>;
};

// Helper for creating element connectivity type for device memory
template<typename ElementTag, typename AcceleratorTag, typename T>
struct ConnectivityTupleHelper;

// Specialization for tetrahedra (device)
template<typename AcceleratorTag, typename T>
struct ConnectivityTupleHelper<TetTag, AcceleratorTag, T>
{
    using DeviceVector = typename VectorSelector<T, AcceleratorTag>::type;

    using type = std::tuple<DeviceVector, DeviceVector, DeviceVector, DeviceVector>;
};

// Specialization for hexahedra (device)
template<typename AcceleratorTag, typename T>
struct ConnectivityTupleHelper<HexTag, AcceleratorTag, T>
{
    using DeviceVector = typename VectorSelector<T, AcceleratorTag>::type;

    using type = std::tuple<DeviceVector,
                            DeviceVector,
                            DeviceVector,
                            DeviceVector,
                            DeviceVector,
                            DeviceVector,
                            DeviceVector,
                            DeviceVector>;
};

// Specializations for triangle (device)
template<typename AcceleratorTag, typename T>
struct ConnectivityTupleHelper<TriTag, AcceleratorTag, T>
{
    using DeviceVector = typename VectorSelector<T, AcceleratorTag>::type;

    using type = std::tuple<DeviceVector, DeviceVector, DeviceVector>;
};

// Specializations for quad (device)
template<typename AcceleratorTag, typename T>
struct ConnectivityTupleHelper<QuadTag, AcceleratorTag, T>
{
    using DeviceVector = typename VectorSelector<T, AcceleratorTag>::type;

    using type = std::tuple<DeviceVector, DeviceVector, DeviceVector, DeviceVector>;
};

// Connectivity tuple helper for host vectors - always std::vector regardless of AcceleratorTag
template<typename ElementTag, typename T>
struct HostConnectivityTupleHelper;

// Specialization for tetrahedra (host)
template<typename T>
struct HostConnectivityTupleHelper<TetTag, T>
{
    using HostVector = typename VectorSelector<T, cstone::CpuTag>::type;
    using type       = std::tuple<HostVector, HostVector, HostVector, HostVector>;
};

// Specialization for hexahedra (host)
template<typename T>
struct HostConnectivityTupleHelper<HexTag, T>
{
    using HostVector = typename VectorSelector<T, cstone::CpuTag>::type;
    using type =
        std::tuple<HostVector, HostVector, HostVector, HostVector, HostVector, HostVector, HostVector, HostVector>;
};

// Specialization for triangles (host)
template<typename T>
struct HostConnectivityTupleHelper<TriTag, T>
{
    using HostVector = typename VectorSelector<T, cstone::CpuTag>::type;
    using type       = std::tuple<HostVector, HostVector, HostVector>;
};

// Specialization for quads (host)
template<typename T>
struct HostConnectivityTupleHelper<QuadTag, T>
{
    using HostVector = typename VectorSelector<T, cstone::CpuTag>::type;
    using type       = std::tuple<HostVector, HostVector, HostVector, HostVector>;
};

// Forward declaration for helper structs
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
class ElementDomain;

// ============================================================================
// AdjacencyData: Manage node-element and element-node mappings
// Lazily initialized for memory efficiency
// ============================================================================
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
struct AdjacencyData
{
    static constexpr int NodesPerElement = ElementTag::NodesPerElement;

    template<typename T>
    using DeviceVector            = typename VectorSelector<T, AcceleratorTag>::type;
    using DeviceConnectivityTuple = typename ConnectivityTupleHelper<ElementTag, AcceleratorTag, KeyType>::type;

    // Node-to-Element inverse mapping (CSR format)
    DeviceVector<KeyType> d_nodeToElementOffsets_;
    DeviceVector<KeyType> d_nodeToElementList_;

    // Element-to-Node connectivity with dense local IDs [0...N-1]
    DeviceConnectivityTuple d_conn_local_ids_;

    // Build adjacency structures from parent domain
    AdjacencyData(const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain);

    void buildNodeToElementMap(const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain);
    void createElementToNodeLocalIdMap(const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain);
};

// ============================================================================
// HaloData: Manages halo (ghost) elements and node ownership for MPI ranks
// Only needed in multi-rank simulations
// ============================================================================
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
struct HaloData
{
    template<typename T>
    using DeviceVector = typename VectorSelector<T, AcceleratorTag>::type;

    // Node ownership: 0=pure ghost (halo only), 1=pure owned (local only), 2=shared (local+halo boundary)
    DeviceVector<uint8_t> d_nodeOwnership_;

    // Cached list of halo element indices for faster iteration
    DeviceVector<KeyType> d_haloElementIndices_;

    // Build halo structures from parent domain
    HaloData(const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain);

    void buildNodeOwnership(const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain);
    void buildHaloElementIndices(const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain);
};

// ============================================================================
// NodeHaloTopology: per-neighbor send/recv lists of LOCAL node IDs for direct
// CUDA-aware MPI exchange of node-DOF arrays (replaces cstone-element-halo
// transport for the per-iteration solve communication path).
//
// Layout:
//   peers_                  list of peer ranks (those we send to or receive from)
//   sendOffsets_[p+1]       CSR offsets into sendNodeIds_ for peer index p
//   sendNodeIds_            flat list of local owned-node IDs to send
//   recvOffsets_[p+1]       CSR offsets into recvNodeIds_ for peer index p
//   recvNodeIds_            flat list of local ghost-node IDs to receive
//
// Exchange: pack via thrust::gather; MPI Isend/Irecv on device pointers
// (CUDA-aware MPI required); scatter via thrust::scatter on receive.
// ============================================================================
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
struct NodeHaloTopology
{
    template<typename T>
    using DeviceVector = typename VectorSelector<T, AcceleratorTag>::type;

    std::vector<int>      peers_;
    std::vector<int>      sendOffsets_;        // size = peers.size() + 1, host
    std::vector<int>      recvOffsets_;
    DeviceVector<int>     sendNodeIds_;        // size = sendOffsets_.back()
    DeviceVector<int>     recvNodeIds_;        // size = recvOffsets_.back()
    // Persistent send/recv staging buffers. Pre-sized in NodeHaloTopology ctor;
    // never resized on hot path. mutable so exchangeNodeHalo() can be const,
    // matching cstone Halos::exchangeHalos pattern; safe because pre-sizing
    // removes the race that bit v1.
    mutable DeviceVector<RealType> sendBuf_;
    mutable DeviceVector<RealType> recvBuf_;
    // Epoch counter for unique MPI tags across multiple exchangeNodeHalo() calls
    // between MPI collectives (matches cstone Halos::haloEpoch_ convention).
    mutable int epoch_{0};

    NodeHaloTopology() = default;
    NodeHaloTopology(const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain);

    // Gate 1 helper: build via the device-resident path consuming cstone's
    // incomingHaloIndices/outgoingHaloIndices instead of host MPI_Allgatherv.
    // Selected at runtime via env MARS_NODEHALO_V2; both paths populate the
    // same fields so callers don't change.
    void buildFromCstoneHalos(const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain);

    // Gate 1 validation: build with both paths and diff per-peer node-id lists.
    // Prints first mismatch to stderr; returns true if identical (after sort).
    static bool validateAgainstHost(const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain);
};

// ============================================================================
// CoordinateCache: Pre-decoded node coordinates from SFC keys
// Optional caching to avoid repeated SFC decoding (trading memory for speed)
// ============================================================================
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
struct CoordinateCache
{
    template<typename T>
    using DeviceVector = typename VectorSelector<T, AcceleratorTag>::type;

    // Cached physical coordinates
    DeviceVector<RealType> d_node_x_;
    DeviceVector<RealType> d_node_y_;
    DeviceVector<RealType> d_node_z_;

    // Build coordinate cache from parent domain
    CoordinateCache(const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain);
};

// ============================================================================
// OriginalCoordinates: Exact original node coordinates from mesh file
// Alternative to SFC-decoded coordinates for high-precision applications
// ============================================================================
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
struct OriginalCoordinates
{
    template<typename T>
    using DeviceVector = typename VectorSelector<T, AcceleratorTag>::type;

    // Original coordinates from mesh file (exact precision)
    DeviceVector<RealType> d_node_x_;
    DeviceVector<RealType> d_node_y_;
    DeviceVector<RealType> d_node_z_;

    // Constructor: allocate storage for nodeCount nodes
    OriginalCoordinates(size_t nodeCount);
};

// Main domain class templated on element type, real type, key type, and accelerator type
template<typename ElementTag     = TetTag,
         typename RealType       = float,
         typename KeyType        = unsigned,
         typename AcceleratorTag = cstone::GpuTag>
class ElementDomain
{
public:
    static constexpr int NodesPerElement = ElementTag::NodesPerElement;

    using DomainType = cstone::Domain<KeyType, RealType, AcceleratorTag>;

    // Template alias for appropriate vector types
    template<typename T>
    using DeviceVector = typename VectorSelector<T, AcceleratorTag>::type;

    // Host vector types are always std::vector
    template<typename T>
    using HostVector = typename VectorSelector<T, cstone::CpuTag>::type;

    // SoA data structures using tuples - device versions
    using DeviceCoordsTuple       = std::tuple<DeviceVector<RealType>, DeviceVector<RealType>, DeviceVector<RealType>>;
    using DevicePropsTuple        = std::tuple<DeviceVector<RealType>>;
    using DeviceConnectivityTuple = typename ConnectivityTupleHelper<ElementTag, AcceleratorTag, KeyType>::type;
    using DeviceBoundaryTuple     = std::tuple<DeviceVector<uint8_t>>; // (isBoundaryNode)

    // SoA data structures using tuples - host versions
    using HostCoordsTuple       = std::tuple<HostVector<RealType>, HostVector<RealType>, HostVector<RealType>>;
    using HostPropsTuple        = std::tuple<HostVector<RealType>>;
    using HostConnectivityTuple = typename HostConnectivityTupleHelper<ElementTag, KeyType>::type;
    using HostBoundaryTuple     = std::tuple<HostVector<uint8_t>>; // (isBoundaryNode)

    // Constructor from file
    ElementDomain(const std::string& meshFile, int rank, int numRanks, bool storeOriginalCoords = false,
                  int bucketSize = 64, unsigned bucketSizeFocus = 8);

    // Constructor from mesh data (for MFEM or other formats) - automatically computes bounding box
    ElementDomain(const HostCoordsTuple& h_coords, const HostConnectivityTuple& h_conn,
                  int rank, int numRanks, int bucketSize = 64, bool storeOriginalCoords = false,
                  unsigned bucketSizeFocus = 8);

    // Constructor from mesh data with boundary info - automatically computes bounding box
    ElementDomain(const HostCoordsTuple& h_coords,
                  const HostConnectivityTuple& h_conn,
                  const HostBoundaryTuple& h_boundary,
                  int rank,
                  int numRanks,
                  int bucketSize = 64,
                  unsigned bucketSizeFocus = 8);

    // Constructor from mesh data with explicit bounding box (for backward compatibility)
    ElementDomain(const HostCoordsTuple& h_coords,
                  const HostConnectivityTuple& h_conn,
                  const cstone::Box<RealType>& box,
                  int rank,
                  int numRanks,
                  int bucketSize = 64,
                  unsigned bucketSizeFocus = 8);

    // Device-data constructor (no host round-trip).
    // The caller hands in already-on-device coords + connectivity. Bbox is
    // computed on GPU + MPI. Used for AMR mesh rebuild after refinement.
    // The input vectors are consumed (moved into the domain's owned buffers
    // for sync scratch space; coords are passed straight to sync()).
    ElementDomain(DeviceCoordsTuple&& d_coords,
                  DeviceConnectivityTuple&& d_conn,
                  int rank,
                  int numRanks,
                  int bucketSize = 64,
                  unsigned bucketSizeFocus = 8);

    // Per-phase timings populated by the ctor (in ms, rank-local).
    // All cudaDeviceSynchronize-fenced so they reflect actual elapsed time.
    // Sub-phases sum to ~ctor time minus any unaccounted overhead (logging, etc).
    float readMeshTimeMs   = 0.0f;
    float bboxTimeMs       = 0.0f;
    float h2dTimeMs        = 0.0f;
    float charSizeTimeMs   = 0.0f;
    float syncTimeMs       = 0.0f;
    float nodeHaloTimeMs   = 0.0f;

    // GPU-accelerated calculation of characteristic sizes
    void calculateCharacteristicSizes(const DeviceConnectivityTuple& d_conn_, const DeviceCoordsTuple& d_coords_);

    // Transfer data to GPU for computations
    void transferDataToGPU(const HostCoordsTuple& h_coords_,
                           const HostConnectivityTuple& h_conn_,
                           DeviceCoordsTuple& d_coords_,
                           DeviceConnectivityTuple& d_conn_);

    // Transfer data to GPU for computations (including boundary info)
    void transferDataToGPU(const HostCoordsTuple& h_coords_,
                           const HostConnectivityTuple& h_conn_,
                           const HostBoundaryTuple& h_boundary_,
                           DeviceCoordsTuple& d_coords_,
                           DeviceConnectivityTuple& d_conn_,
                           DeviceBoundaryTuple& d_boundary_);

    // Domain synchronization (following cornerstone API pattern)
    void sync(const DeviceConnectivityTuple& d_conn_, const DeviceCoordsTuple& d_coords_);

    // Resync against a refined mesh WITHOUT reconstructing the cstone Domain.
    // Avoids firing focusTree_.converge() (which is firstCall_-only and dominates
    // AMR rebuild time at scale: ~80% of cstone sync at 977M elements/16 ranks).
    //
    // Contract:
    //   - Replaces the per-instance mesh state (coords, conn, props, sfc keys)
    //   - Recomputes bounding box via MPI_Allreduce (cheap, < 1 ms)
    //   - Resets all lazy-cached sub-objects (adjacency_, halo_, nodeHaloTopo_,
    //     coordCache_, originalCoords_) so they re-init against the new mesh
    //   - Reuses the existing cstone Domain object: subsequent sync() observes
    //     firstCall_=false and skips focusTree_.converge()
    //   - cstone's do-while retry loop in sync() handles the element-count
    //     change via its updateTree convergence check
    //
    // Required for Gordon-Bell scale: at 10^5+ ranks, firstCall_=true is also
    // gated on a host-replicated global octree Allgather, which is the actual
    // scaling wall. Domain reuse eliminates that as well.
    void resyncFromDevice(DeviceCoordsTuple&& d_coords_new,
                          DeviceConnectivityTuple&& d_conn_new);

    // Access to connectivity - template parameter to select index array
    template<int I>
    const DeviceVector<KeyType>& indices() const
    {
        return std::get<I>(d_conn_keys_);
    }

    template<int I>
    DeviceVector<KeyType>& getConnectivity()
    {
        return std::get<I>(d_conn_keys_);
    }

    // User-friendly coordinate access (uses cached coords if available, otherwise decodes on-the-fly)
    auto getElementNodeCoordinates(size_t elemIdx)
    {
        std::array<RealType, NodesPerElement> x, y, z;

        // Check if we have cached coordinates
        bool hasCachedCoords = (coordCache_ != nullptr);

        if (hasCachedCoords)
        {
            // Use cached coordinates via local IDs
            const auto& conn_local = getElementToNodeConnectivity();
            for (int i = 0; i < NodesPerElement; ++i)
            {
                KeyType localId = [&]()
                {
                    if (i == 0) return std::get<0>(conn_local)[elemIdx];
                    if (i == 1) return std::get<1>(conn_local)[elemIdx];
                    if (i == 2) return std::get<2>(conn_local)[elemIdx];
                    if (i == 3 && NodesPerElement > 3) return std::get<3>(conn_local)[elemIdx];
                    return KeyType(0);
                }();
                x[i] = coordCache_->d_node_x_[localId];
                y[i] = coordCache_->d_node_y_[localId];
                z[i] = coordCache_->d_node_z_[localId];
            }
        }
        else
        {
            // Decode from SFC keys on-the-fly
            for (int i = 0; i < NodesPerElement; ++i)
            {
                KeyType sfcKey = [&]()
                {
                    if (i == 0) return std::get<0>(d_conn_keys_)[elemIdx];
                    if (i == 1) return std::get<1>(d_conn_keys_)[elemIdx];
                    if (i == 2) return std::get<2>(d_conn_keys_)[elemIdx];
                    if (i == 3 && NodesPerElement > 3) return std::get<3>(d_conn_keys_)[elemIdx];
                    return KeyType(0);
                }();
                auto [xi, yi, zi] = sfcToPhysicalCoordinate(sfcKey);
                x[i]              = xi;
                y[i]              = yi;
                z[i]              = zi;
            }
        }

        return std::make_tuple(x, y, z);
    }

    // Helper to get local node ID from connectivity (requires adjacency to be built)
    template<int I>
    KeyType getLocalNodeId(size_t elemIdx) const
    {
        static_assert(I >= 0 && I < NodesPerElement, "Index out of range");
        if (elemIdx >= elementCount_) return KeyType(0);
        ensureAdjacency();
        return std::get<I>(adjacency_->d_conn_local_ids_)[elemIdx];
    }

    // Get element/node counts
    std::size_t getNodeCount() const { return nodeCount_; }
    std::size_t getElementCount() const { return elementCount_; }

    // Access to cornerstone domain
    DomainType& getDomain() { return *domain_; }
    const DomainType& getDomain() const { return *domain_; }

    MARS_HOST_DEVICE
    const cstone::Box<RealType>& getBoundingBox() const { return box_; }

    // Start and end indices for local work assignment
    std::size_t startIndex() const { return domain_->startIndex(); }
    std::size_t endIndex() const { return domain_->endIndex(); }
    std::size_t localElementCount() const { return endIndex() - startIndex(); }
    std::size_t haloElementCount() const { return getElementCount() - localElementCount(); }

    // Getter for halo count
    size_t getHaloElementCount() const
    {
        ensureHalo();
        return halo_->d_haloElementIndices_.size();
    }

    // Get ranges for local elements (can modify in FEM assembly)
    std::pair<size_t, size_t> localElementRange() const { return {startIndex(), endIndex()}; }

    // Get all halo ranges (read-only for FEM)
    std::vector<std::pair<size_t, size_t>> haloElementRanges() const
    {
        std::vector<std::pair<size_t, size_t>> ranges;
        if (startIndex() > 0)
        {
            ranges.emplace_back(0, startIndex()); // Halos before local
        }
        if (endIndex() < getElementCount())
        {
            ranges.emplace_back(endIndex(), getElementCount()); // Halos after local
        }
        return ranges;
    }

    // Check if an element is local (for FEM: can assemble contributions)
    bool isLocalElement(size_t elementIndex) const { return elementIndex >= startIndex() && elementIndex < endIndex(); }

    // Check if an element is a halo (for FEM: read-only, used for stencils)
    bool isHaloElement(size_t elementIndex) const { return !isLocalElement(elementIndex); }

    // Helper methods
    void readMeshDataSoA(const std::string& meshFile, HostCoordsTuple& h_coords_, HostConnectivityTuple& h_conn_);

    // Real coordinate conversion from SFC
    MARS_HOST_DEVICE std::tuple<RealType, RealType, RealType> sfcToPhysicalCoordinate(KeyType sfcKey) const;
    MARS_HOST_DEVICE RealType sfcToPhysicalCoordinateX(KeyType sfcKey) const;
    MARS_HOST_DEVICE RealType sfcToPhysicalCoordinateY(KeyType sfcKey) const;
    MARS_HOST_DEVICE RealType sfcToPhysicalCoordinateZ(KeyType sfcKey) const;

    // Integer coordinate conversion from SFC (for spatial operations/hashing)
    MARS_HOST_DEVICE std::tuple<unsigned, unsigned, unsigned> sfcToSpatialCoordinate(KeyType sfcKey) const;
    MARS_HOST_DEVICE unsigned sfcToSpatialCoordinateX(KeyType sfcKey) const;
    MARS_HOST_DEVICE unsigned sfcToSpatialCoordinateY(KeyType sfcKey) const;
    MARS_HOST_DEVICE unsigned sfcToSpatialCoordinateZ(KeyType sfcKey) const;

    // On-demand connectivity access
    template<int I>
    DeviceVector<KeyType> getConnectivity() const;
    DeviceConnectivityTuple getConnectivity() const;

    // Single connectivity access
    template<int I>
    MARS_HOST_DEVICE KeyType getConnectivity(size_t elementIndex) const;

    std::tuple<KeyType, KeyType, KeyType, KeyType> getConnectivity(size_t elementIndex) const;

    // getter for ownership map (lazy)
    const DeviceVector<uint8_t>& getNodeOwnershipMap() const
    {
        ensureHalo();
        return halo_->d_nodeOwnership_;
    }

    // getter for halo element indices (lazy)
    const DeviceVector<KeyType>& getHaloElementIndices() const
    {
        ensureHalo();
        return halo_->d_haloElementIndices_;
    }

    // build adjacency after sync (explicit call, no longer in sync())
    void buildAdjacency() { ensureAdjacency(); }

    // optional: explicitly request coordinate caching
    void cacheNodeCoordinates()
    {
        // If using original coordinates, they're already stored - no need to cache
        if (originalCoords_) return;
        ensureCoordinateCache();
    }

    // Exchange halo data for DOF vectors (MPI communication to sum shared node contributions)
    template<class... Vectors, class SendBuffer, class ReceiveBuffer>
    void exchangeHalos(std::tuple<Vectors&...> arrays, SendBuffer& sendBuffer, ReceiveBuffer& receiveBuffer) const
    {
        domain_->exchangeHalos(arrays, sendBuffer, receiveBuffer);
    }

    // Per-NODE halo exchange built on top of cstone's per-ELEMENT halo.
    //
    // Inputs:
    //   nodeArray  size = nodeCount, holds per-node DOF values. Owned-node slots
    //              must contain authoritative local values; ghost-node slots are
    //              overwritten with their owners' values on return.
    //   nodeToDof  size = nodeCount; node -> DOF index in nodeArray. Pass nullptr
    //              if nodeArray is indexed directly by local node id.
    //
    // Algorithm: pack per-element corner arrays gated by node ownership (sentinel
    // NaN where THIS rank doesn't own the corner node). cstone halo-exchanges the
    // 8/4 corner arrays element-keyed. Unpack writes into ghost-node slots only,
    // and only when the received value is non-sentinel. This guarantees ghost-DOF
    // values come from each ghost node's actual owner rank (not from an
    // intermediate rank that holds the same node as a stale ghost).
    //
    // Correctness invariant: for every ghost node N on this rank, the cstone
    // halo must include at least one element owned by N's owner. Holds for
    // typical FEM partitions where halo width covers full neighbor element layers.
    template<class VectorType>
    void exchangeNodeHalo(VectorType& nodeArray, const int* nodeToDof = nullptr) const
    {
        if (numRanks_ == 1) return;

        using T = typename VectorType::value_type;
        static_assert(std::is_same_v<T, RealType>,
                      "exchangeNodeHalo: vector element type must match domain RealType");

        ensureNodeHaloTopo();
        const auto& topo = *nodeHaloTopo_;

        // Buffers pre-sized in NodeHaloTopology ctor; no hot-path resize.
        size_t sendTotal = topo.sendOffsets_.empty() ? 0 : size_t(topo.sendOffsets_.back());
        size_t recvTotal = topo.recvOffsets_.empty() ? 0 : size_t(topo.recvOffsets_.back());

        // Pack: index-driven gather. No sentinel: the recv list on the peer
        // side determines which slots get overwritten on unpack, so any value
        // we pack here that doesn't correspond to a valid DOF will simply
        // never be read by anyone.
        T*           arr  = thrust::raw_pointer_cast(nodeArray.data());
        const int*   snds = thrust::raw_pointer_cast(topo.sendNodeIds_.data());
        T*           sbuf = thrust::raw_pointer_cast(topo.sendBuf_.data());
        if (sendTotal > 0)
        {
            thrust::for_each(thrust::device,
                              thrust::counting_iterator<size_t>(0),
                              thrust::counting_iterator<size_t>(sendTotal),
                              [arr, snds, sbuf, nodeToDof] __device__ (size_t i) {
                                  int n = snds[i];
                                  int idx = (nodeToDof != nullptr) ? nodeToDof[n] : n;
                                  // idx<0 means node is not a DOF on this rank.
                                  // Pack arr[0] as a placeholder; the corresponding
                                  // recv-side slot (if any) will be filtered the
                                  // same way on unpack so the value is discarded.
                                  sbuf[i] = (idx >= 0) ? arr[idx] : arr[0];
                              });
            cudaDeviceSynchronize();
        }

        // Post recvs and sends. CUDA-aware MPI: pass device pointers directly.
        // Tag = nodeHaloTagBase + epoch_, matching cstone's haloEpoch_ pattern
        // so multiple exchanges between MPI collectives don't tag-collide.
        T*           rbuf = thrust::raw_pointer_cast(topo.recvBuf_.data());
        auto mpiType = std::is_same_v<T, double> ? MPI_DOUBLE : MPI_FLOAT;
        constexpr int nodeHaloTagBase = 0x4d52;  // "MR" — disambiguate from cstone tags
        const int tag = nodeHaloTagBase + topo.epoch_;
        std::vector<MPI_Request> reqs;
        reqs.reserve(2 * topo.peers_.size());

        for (size_t p = 0; p < topo.peers_.size(); ++p)
        {
            int peer = topo.peers_[p];
            int rcnt = topo.recvOffsets_[p+1] - topo.recvOffsets_[p];
            if (rcnt > 0)
            {
                MPI_Request r;
                MPI_Irecv(rbuf + topo.recvOffsets_[p], rcnt, mpiType, peer,
                          tag, MPI_COMM_WORLD, &r);
                reqs.push_back(r);
            }
        }
        for (size_t p = 0; p < topo.peers_.size(); ++p)
        {
            int peer = topo.peers_[p];
            int scnt = topo.sendOffsets_[p+1] - topo.sendOffsets_[p];
            if (scnt > 0)
            {
                MPI_Request r;
                MPI_Isend(sbuf + topo.sendOffsets_[p], scnt, mpiType, peer,
                          tag, MPI_COMM_WORLD, &r);
                reqs.push_back(r);
            }
        }
        if (!reqs.empty()) MPI_Waitall(int(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
        ++topo.epoch_;

        // Scatter recvBuf into nodeArray[nodeToDof[recvNodeIds]].
        // idx<0 case skipped: the slot was packed by the peer with a placeholder
        // and we don't have a real DOF to write it to.
        const int* rnds = thrust::raw_pointer_cast(topo.recvNodeIds_.data());
        if (recvTotal > 0)
        {
            thrust::for_each(thrust::device,
                              thrust::counting_iterator<size_t>(0),
                              thrust::counting_iterator<size_t>(recvTotal),
                              [arr, rnds, rbuf, nodeToDof] __device__ (size_t i) {
                                  int n = rnds[i];
                                  int idx = (nodeToDof != nullptr) ? nodeToDof[n] : n;
                                  if (idx >= 0) arr[idx] = rbuf[i];
                              });
            cudaDeviceSynchronize();
        }
        return;
    }

    // Device connectivity functions
    template<int I>
    __device__ KeyType getConnectivityDevice(size_t elementIndex) const;

    int numRanks() const { return numRanks_; }
    int rank() const { return rank_; }

    //! Returns the map from a dense local node ID to its sparse global SFC Key.
    const DeviceVector<KeyType>& getLocalToGlobalSfcMap() const
    {
        ensureSfcMap();
        return d_localToGlobalSfcMap_;
    }

    const DeviceVector<KeyType>& getLocalToGlobalNodeMap() const { return d_localToGlobalNodeMap_; }

    //! Returns the Element->Node connectivity table using dense local node IDs (lazy).
    const DeviceConnectivityTuple& getElementToNodeConnectivity() const
    {
        ensureAdjacency();
        return adjacency_->d_conn_local_ids_;
    }

    //! Returns the offsets for the Node->Element CSR map (lazy).
    const DeviceVector<KeyType>& getNodeToElementOffsets() const
    {
        ensureAdjacency();
        return adjacency_->d_nodeToElementOffsets_;
    }

    //! Returns the list for the Node->Element CSR map (lazy).
    const DeviceVector<KeyType>& getNodeToElementList() const
    {
        ensureAdjacency();
        return adjacency_->d_nodeToElementList_;
    }

    //! Returns cached node coordinates (lazy, initializes if not already cached).
    const DeviceVector<RealType>& getNodeX() const
    {
        if (originalCoords_) {
            return originalCoords_->d_node_x_;
        }
        ensureCoordinateCache();
        return coordCache_->d_node_x_;
    }

    const DeviceVector<RealType>& getNodeY() const
    {
        if (originalCoords_) {
            return originalCoords_->d_node_y_;
        }
        ensureCoordinateCache();
        return coordCache_->d_node_y_;
    }

    const DeviceVector<RealType>& getNodeZ() const
    {
        if (originalCoords_) {
            return originalCoords_->d_node_z_;
        }
        ensureCoordinateCache();
        return coordCache_->d_node_z_;
    }

    //! Returns boundary node flags (true if node is on boundary).
    const DeviceVector<uint8_t>& getBoundaryNodes() const { return std::get<0>(d_boundary_); }

    //! Returns true if boundary information is available.
    bool hasBoundaryInfo() const { return std::get<0>(d_boundary_).size() > 0; }

private:
    int rank_;
    int numRanks_;
    std::unique_ptr<DomainType> domain_;
    std::size_t elementCount_ = 0;
    std::size_t nodeCount_    = 0;

    // Flag to store original coordinates instead of SFC-decoded coords
    bool storeOriginalCoords_ = false;

    // Core mesh data (always allocated)
    // element unique identifiers (SFC codes)
    DeviceVector<KeyType> d_elemSfcCodes_;

    // Device data in SoA format
    DevicePropsTuple d_props_; // (h)
    DeviceConnectivityTuple
        d_conn_keys_; // (sfc_i0, sfc_i1, sfc_i2, ...) depends on element type; store sfc instead of coordinates
    DeviceBoundaryTuple d_boundary_; // (isBoundaryNode) - optional boundary node flags
    // TODO: check if domain_->box() is a device function in cstone to avoid storing box_ here
    cstone::Box<RealType> box_; // bounding box for the domain

    // useful mapping from element to node indices
    DeviceVector<KeyType> d_elemToNodeMap_;

    // Lazily-initialized components (allocated on demand)
    // Maps a dense local node ID [0...N-1] to its sparse global SFC Key
    DeviceVector<KeyType> d_localToGlobalSfcMap_;

    // Maps a dense local node ID [0...N-1] to its global node index
    DeviceVector<KeyType> d_localToGlobalNodeMap_;

    // Friend declarations to allow helper structs access to private members
    friend struct AdjacencyData<ElementTag, RealType, KeyType, AcceleratorTag>;
    friend struct HaloData<ElementTag, RealType, KeyType, AcceleratorTag>;
    friend struct NodeHaloTopology<ElementTag, RealType, KeyType, AcceleratorTag>;
    friend struct CoordinateCache<ElementTag, RealType, KeyType, AcceleratorTag>;
    friend struct OriginalCoordinates<ElementTag, RealType, KeyType, AcceleratorTag>;

    // Free-function builder for NodeHaloTopology's host path needs access to halo_ and friends.
    template<typename ET, typename RT, typename KT, typename AT>
    friend void buildNodeHaloTopologyHostPath(NodeHaloTopology<ET, RT, KT, AT>& topo,
                                              const ElementDomain<ET, RT, KT, AT>& domain);

    mutable std::unique_ptr<AdjacencyData<ElementTag, RealType, KeyType, AcceleratorTag>> adjacency_;
    mutable std::unique_ptr<HaloData<ElementTag, RealType, KeyType, AcceleratorTag>> halo_;
    mutable std::unique_ptr<NodeHaloTopology<ElementTag, RealType, KeyType, AcceleratorTag>> nodeHaloTopo_;
    mutable std::unique_ptr<CoordinateCache<ElementTag, RealType, KeyType, AcceleratorTag>> coordCache_;
    mutable std::unique_ptr<OriginalCoordinates<ElementTag, RealType, KeyType, AcceleratorTag>> originalCoords_;

    void initializeConnectivityKeys();

    // Helper method for creating SFC map
    void createLocalToGlobalSfcMap();

    // Lazy initialization methods (called by public getters)
    void ensureSfcMap() const;
    void ensureAdjacency() const;
    void ensureHalo() const;
    void ensureNodeHaloTopo() const;
    void ensureCoordinateCache() const;


    void syncImpl(DeviceVector<KeyType>& d_nodeSfcCodes,
                  const DeviceConnectivityTuple& d_conn_,
                  int blockSize,
                  int numBlocks)
    {
        if constexpr (std::is_same_v<ElementTag, TetTag>)
        {
            auto& d_i0 = std::get<0>(d_conn_);
            auto& d_i1 = std::get<1>(d_conn_);
            auto& d_i2 = std::get<2>(d_conn_);
            auto& d_i3 = std::get<3>(d_conn_);

            // Find representative nodes
            findRepresentativeNodesKernel<TetTag, KeyType, RealType><<<numBlocks, blockSize>>>(
                thrust::raw_pointer_cast(d_i0.data()), thrust::raw_pointer_cast(d_i1.data()),
                thrust::raw_pointer_cast(d_i2.data()), thrust::raw_pointer_cast(d_i3.data()), nullptr, nullptr, nullptr,
                nullptr, thrust::raw_pointer_cast(d_nodeSfcCodes.data()),
                thrust::raw_pointer_cast(d_elemToNodeMap_.data()), elementCount_);
            cudaCheckError();

            // Build SFC connectivity
            buildSfcConnectivity<TetTag, KeyType, RealType><<<numBlocks, blockSize>>>(
                thrust::raw_pointer_cast(d_i0.data()), thrust::raw_pointer_cast(d_i1.data()),
                thrust::raw_pointer_cast(d_i2.data()), thrust::raw_pointer_cast(d_i3.data()), nullptr, nullptr, nullptr,
                nullptr, thrust::raw_pointer_cast(d_nodeSfcCodes.data()),
                thrust::raw_pointer_cast(std::get<0>(d_conn_keys_).data()),
                thrust::raw_pointer_cast(std::get<1>(d_conn_keys_).data()),
                thrust::raw_pointer_cast(std::get<2>(d_conn_keys_).data()),
                thrust::raw_pointer_cast(std::get<3>(d_conn_keys_).data()), nullptr, nullptr, nullptr, nullptr,
                elementCount_);
            cudaCheckError();
        }
        else if constexpr (std::is_same_v<ElementTag, HexTag>)
        {
            auto& d_i0 = std::get<0>(d_conn_);
            auto& d_i1 = std::get<1>(d_conn_);
            auto& d_i2 = std::get<2>(d_conn_);
            auto& d_i3 = std::get<3>(d_conn_);
            auto& d_i4 = std::get<4>(d_conn_);
            auto& d_i5 = std::get<5>(d_conn_);
            auto& d_i6 = std::get<6>(d_conn_);
            auto& d_i7 = std::get<7>(d_conn_);

            // Find representative nodes
            findRepresentativeNodesKernel<HexTag, KeyType, RealType>
                <<<numBlocks, blockSize>>>(thrust::raw_pointer_cast(d_i0.data()), thrust::raw_pointer_cast(d_i1.data()),
                                           thrust::raw_pointer_cast(d_i2.data()), thrust::raw_pointer_cast(d_i3.data()),
                                           thrust::raw_pointer_cast(d_i4.data()), thrust::raw_pointer_cast(d_i5.data()),
                                           thrust::raw_pointer_cast(d_i6.data()), thrust::raw_pointer_cast(d_i7.data()),
                                           thrust::raw_pointer_cast(d_nodeSfcCodes.data()),
                                           thrust::raw_pointer_cast(d_elemToNodeMap_.data()), elementCount_);
            cudaCheckError();

            // Build SFC connectivity for hex
            buildSfcConnectivity<HexTag, KeyType, RealType>
                <<<numBlocks, blockSize>>>(thrust::raw_pointer_cast(d_i0.data()), thrust::raw_pointer_cast(d_i1.data()),
                                           thrust::raw_pointer_cast(d_i2.data()), thrust::raw_pointer_cast(d_i3.data()),
                                           thrust::raw_pointer_cast(d_i4.data()), thrust::raw_pointer_cast(d_i5.data()),
                                           thrust::raw_pointer_cast(d_i6.data()), thrust::raw_pointer_cast(d_i7.data()),
                                           thrust::raw_pointer_cast(d_nodeSfcCodes.data()),
                                           thrust::raw_pointer_cast(std::get<0>(d_conn_keys_).data()),
                                           thrust::raw_pointer_cast(std::get<1>(d_conn_keys_).data()),
                                           thrust::raw_pointer_cast(std::get<2>(d_conn_keys_).data()),
                                           thrust::raw_pointer_cast(std::get<3>(d_conn_keys_).data()),
                                           thrust::raw_pointer_cast(std::get<4>(d_conn_keys_).data()),
                                           thrust::raw_pointer_cast(std::get<5>(d_conn_keys_).data()),
                                           thrust::raw_pointer_cast(std::get<6>(d_conn_keys_).data()),
                                           thrust::raw_pointer_cast(std::get<7>(d_conn_keys_).data()), elementCount_);
            cudaCheckError();
        }
        else if constexpr (std::is_same_v<ElementTag, TriTag>)
        {
            auto& d_i0 = std::get<0>(d_conn_);
            auto& d_i1 = std::get<1>(d_conn_);
            auto& d_i2 = std::get<2>(d_conn_);

            // Find representative nodes for triangles
            findRepresentativeNodesKernel<TriTag, KeyType, RealType><<<numBlocks, blockSize>>>(
                thrust::raw_pointer_cast(d_i0.data()), thrust::raw_pointer_cast(d_i1.data()),
                thrust::raw_pointer_cast(d_i2.data()), nullptr, // This one stays nullptr since tri only has 3 nodes
                thrust::raw_pointer_cast(d_nodeSfcCodes.data()), thrust::raw_pointer_cast(d_elemToNodeMap_.data()),
                elementCount_);
            cudaCheckError();

            // Build SFC connectivity for triangles (would need 3-node version)
            // buildSfcConnectivity<TriTag, KeyType, RealType><<<...>>>(...);
        }
        else if constexpr (std::is_same_v<ElementTag, QuadTag>)
        {
            auto& d_i0 = std::get<0>(d_conn_);
            auto& d_i1 = std::get<1>(d_conn_);
            auto& d_i2 = std::get<2>(d_conn_);
            auto& d_i3 = std::get<3>(d_conn_);

            // Find representative nodes for quads
            findRepresentativeNodesKernel<QuadTag, KeyType, RealType>
                <<<numBlocks, blockSize>>>(thrust::raw_pointer_cast(d_i0.data()), thrust::raw_pointer_cast(d_i1.data()),
                                           thrust::raw_pointer_cast(d_i2.data()), thrust::raw_pointer_cast(d_i3.data()),
                                           thrust::raw_pointer_cast(d_nodeSfcCodes.data()),
                                           thrust::raw_pointer_cast(d_elemToNodeMap_.data()), elementCount_);
            cudaCheckError();

            // Build SFC connectivity for quads
            buildSfcConnectivity<QuadTag, KeyType, RealType>
                <<<numBlocks, blockSize>>>(thrust::raw_pointer_cast(d_i0.data()), thrust::raw_pointer_cast(d_i1.data()),
                                           thrust::raw_pointer_cast(d_i2.data()), thrust::raw_pointer_cast(d_i3.data()),
                                           thrust::raw_pointer_cast(d_nodeSfcCodes.data()),
                                           thrust::raw_pointer_cast(std::get<0>(d_conn_keys_).data()),
                                           thrust::raw_pointer_cast(std::get<1>(d_conn_keys_).data()),
                                           thrust::raw_pointer_cast(std::get<2>(d_conn_keys_).data()),
                                           thrust::raw_pointer_cast(std::get<3>(d_conn_keys_).data()), elementCount_);
            cudaCheckError();
        }
    }
};

template<typename RealType>
cstone::Box<RealType> computeGlobalBoundingBox(const std::string& meshDir)
{
    // Use the SAME extension logic as readMeshWithElementPartitioning to ensure consistency
    std::string ext;
    if constexpr (std::is_same_v<RealType, float>) { ext = "float32"; }
    else if constexpr (std::is_same_v<RealType, double>) { ext = "double"; }
    else { throw std::runtime_error("Unsupported RealType for coordinate reading"); }

    std::ifstream x_file(meshDir + "/x." + ext, std::ios::binary);
    std::ifstream y_file(meshDir + "/y." + ext, std::ios::binary);
    std::ifstream z_file(meshDir + "/z." + ext, std::ios::binary);

    if (!x_file || !y_file || !z_file)
    {
        throw std::runtime_error("Failed to open coordinate files for global bounding box computation: " + meshDir +
                                 " (looking for ." + ext + " files)");
    }

    // Get total node count
    x_file.seekg(0, std::ios::end);
    size_t fileSize = x_file.tellg();
    x_file.seekg(0, std::ios::beg);

    size_t totalNodes;
    if (ext == "float32") { totalNodes = fileSize / sizeof(float); }
    else { totalNodes = fileSize / sizeof(double); }

    // Read ALL data at once (like createBoundingBox approach)
    if (ext == "float32")
    {
        std::vector<float> x_data(totalNodes), y_data(totalNodes), z_data(totalNodes);

        x_file.read(reinterpret_cast<char*>(x_data.data()), totalNodes * sizeof(float));
        y_file.read(reinterpret_cast<char*>(y_data.data()), totalNodes * sizeof(float));
        z_file.read(reinterpret_cast<char*>(z_data.data()), totalNodes * sizeof(float));

        // Use STL algorithms (like createBoundingBox)
        auto [xmin_it, xmax_it] = std::minmax_element(x_data.begin(), x_data.end());
        auto [ymin_it, ymax_it] = std::minmax_element(y_data.begin(), y_data.end());
        auto [zmin_it, zmax_it] = std::minmax_element(z_data.begin(), z_data.end());

        // Convert to RealType if needed
        RealType xmin = static_cast<RealType>(*xmin_it);
        RealType xmax = static_cast<RealType>(*xmax_it);
        RealType ymin = static_cast<RealType>(*ymin_it);
        RealType ymax = static_cast<RealType>(*ymax_it);
        RealType zmin = static_cast<RealType>(*zmin_it);
        RealType zmax = static_cast<RealType>(*zmax_it);

        // Add padding
        RealType padding = 0.05;
        RealType xRange  = xmax - xmin;
        RealType yRange  = ymax - ymin;
        RealType zRange  = zmax - zmin;

        return cstone::Box<RealType>(xmin - padding * xRange, xmax + padding * xRange, ymin - padding * yRange,
                                     ymax + padding * yRange, zmin - padding * zRange, zmax + padding * zRange);
    }
    else
    {
        // Similar for double files
        std::vector<double> x_data(totalNodes), y_data(totalNodes), z_data(totalNodes);

        x_file.read(reinterpret_cast<char*>(x_data.data()), totalNodes * sizeof(double));
        y_file.read(reinterpret_cast<char*>(y_data.data()), totalNodes * sizeof(double));
        z_file.read(reinterpret_cast<char*>(z_data.data()), totalNodes * sizeof(double));

        auto [xmin_it, xmax_it] = std::minmax_element(x_data.begin(), x_data.end());
        auto [ymin_it, ymax_it] = std::minmax_element(y_data.begin(), y_data.end());
        auto [zmin_it, zmax_it] = std::minmax_element(z_data.begin(), z_data.end());

        RealType xmin = static_cast<RealType>(*xmin_it);
        RealType xmax = static_cast<RealType>(*xmax_it);
        RealType ymin = static_cast<RealType>(*ymin_it);
        RealType ymax = static_cast<RealType>(*ymax_it);
        RealType zmin = static_cast<RealType>(*zmin_it);
        RealType zmax = static_cast<RealType>(*zmax_it);

        RealType padding = 0.05;
        RealType xRange  = xmax - xmin;
        RealType yRange  = ymax - ymin;
        RealType zRange  = zmax - zmin;

        return cstone::Box<RealType>(xmin - padding * xRange, xmax + padding * xRange, ymin - padding * yRange,
                                     ymax + padding * yRange, zmin - padding * zRange, zmax + padding * zRange);
    }
}

// Warning: This would create a local bounding box if called on a part on local coords.
template<typename RealType, typename CoordTuple>
cstone::Box<RealType> createBoundingBox(const CoordTuple& coords, RealType padding = 0.05)
{
    std::cout << "WARNING: Creating LOCAL bounding box!" << std::endl;
    std::cout << "    This may cause coordinate displacement in multi-rank SFC applications!" << std::endl;
    std::cout << "    Consider using computeGlobalBoundingBox() instead." << std::endl;

    // Convert device vectors to host vectors first if needed
    std::vector<RealType> h_x, h_y, h_z;

    // Get the vectors either directly or convert from device to host
    if constexpr (std::is_same_v<std::decay_t<decltype(std::get<0>(coords))>, std::vector<RealType>>)
    {
        // Host vectors - direct reference
        const auto& x = std::get<0>(coords);
        const auto& y = std::get<1>(coords);
        const auto& z = std::get<2>(coords);

        // Find min/max coordinates
        auto minmax_x = std::minmax_element(x.begin(), x.end());
        auto minmax_y = std::minmax_element(y.begin(), y.end());
        auto minmax_z = std::minmax_element(z.begin(), z.end());

        // Create box with padding
        RealType x_min = *minmax_x.first - padding;
        RealType x_max = *minmax_x.second + padding;
        RealType y_min = *minmax_y.first - padding;
        RealType y_max = *minmax_y.second + padding;
        RealType z_min = *minmax_z.first - padding;
        RealType z_max = *minmax_z.second + padding;

        return cstone::Box<RealType>(x_min, x_max, y_min, y_max, z_min, z_max);
    }
    else
    {
        // Device vectors - convert to host
        h_x = toHost(std::get<0>(coords));
        h_y = toHost(std::get<1>(coords));
        h_z = toHost(std::get<2>(coords));

        // Find min/max coordinates on host vectors
        auto minmax_x = std::minmax_element(h_x.begin(), h_x.end());
        auto minmax_y = std::minmax_element(h_y.begin(), h_y.end());
        auto minmax_z = std::minmax_element(h_z.begin(), h_z.end());

        // Create box with padding
        RealType x_min = *minmax_x.first - padding;
        RealType x_max = *minmax_x.second + padding;
        RealType y_min = *minmax_y.first - padding;
        RealType y_max = *minmax_y.second + padding;
        RealType z_min = *minmax_z.first - padding;
        RealType z_max = *minmax_z.second + padding;

        return cstone::Box<RealType>(x_min, x_max, y_min, y_max, z_min, z_max);
    }
}

// Compute global bounding box from coordinate data (host vectors)
template<typename RealType>
cstone::Box<RealType> computeGlobalBoundingBoxFromCoords(const std::vector<RealType>& x_coords,
                                                         const std::vector<RealType>& y_coords,
                                                         const std::vector<RealType>& z_coords,
                                                         RealType padding = 0.05)
{
    if (x_coords.empty() || y_coords.empty() || z_coords.empty())
    {
        throw std::runtime_error("Cannot compute bounding box from empty coordinate data");
    }

    // Find LOCAL min/max coordinates
    auto [xmin_it, xmax_it] = std::minmax_element(x_coords.begin(), x_coords.end());
    auto [ymin_it, ymax_it] = std::minmax_element(y_coords.begin(), y_coords.end());
    auto [zmin_it, zmax_it] = std::minmax_element(z_coords.begin(), z_coords.end());

    RealType local_xmin = *xmin_it;
    RealType local_xmax = *xmax_it;
    RealType local_ymin = *ymin_it;
    RealType local_ymax = *ymax_it;
    RealType local_zmin = *zmin_it;
    RealType local_zmax = *zmax_it;

    // Compute GLOBAL min/max across all ranks using MPI
    RealType global_xmin, global_xmax, global_ymin, global_ymax, global_zmin, global_zmax;
    MPI_Allreduce(&local_xmin, &global_xmin, 1, std::is_same_v<RealType, float> ? MPI_FLOAT : MPI_DOUBLE, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&local_xmax, &global_xmax, 1, std::is_same_v<RealType, float> ? MPI_FLOAT : MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&local_ymin, &global_ymin, 1, std::is_same_v<RealType, float> ? MPI_FLOAT : MPI_DOUBLE, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&local_ymax, &global_ymax, 1, std::is_same_v<RealType, float> ? MPI_FLOAT : MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&local_zmin, &global_zmin, 1, std::is_same_v<RealType, float> ? MPI_FLOAT : MPI_DOUBLE, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&local_zmax, &global_zmax, 1, std::is_same_v<RealType, float> ? MPI_FLOAT : MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);

    // Add padding based on GLOBAL ranges
    RealType xRange = global_xmax - global_xmin;
    RealType yRange = global_ymax - global_ymin;
    RealType zRange = global_zmax - global_zmin;

    return cstone::Box<RealType>(global_xmin - padding * xRange, global_xmax + padding * xRange,
                                 global_ymin - padding * yRange, global_ymax + padding * yRange,
                                 global_zmin - padding * zRange, global_zmax + padding * zRange);
}

template<typename KeyType, typename RealType>
void testSfcPrecision(const cstone::Box<RealType>& box, int rank = 0)
{
#ifndef NDEBUG // Only run in debug builds (when NDEBUG is NOT defined)
    RealType domainSizeX   = box.xmax() - box.xmin();
    RealType domainSizeY   = box.ymax() - box.ymin();
    RealType domainSizeZ   = box.zmax() - box.zmin();
    RealType maxDomainSize = std::max({domainSizeX, domainSizeY, domainSizeZ});

    constexpr int totalBits  = sizeof(KeyType) * 8;
    constexpr int bitsPerDim = totalBits / 3;
    RealType sfcPrecision    = maxDomainSize / (1ULL << bitsPerDim);

    if (rank == 0)
    {
        std::cout << "=== SFC Precision Analysis (DEBUG) ===" << std::endl;
        std::cout << "KeyType: " << (sizeof(KeyType) == 4 ? "uint32_t" : "uint64_t") << " (" << totalBits << " bits)"
                  << std::endl;
        std::cout << "Domain size: [" << domainSizeX << " x " << domainSizeY << " x " << domainSizeZ << "]"
                  << std::endl;
        std::cout << "SFC precision: " << sfcPrecision << " units" << std::endl;

        if (sizeof(KeyType) == 4 && sfcPrecision > 0.001)
        {
            std::cout << "⚠️  WARNING: Consider using uint64_t KeyType!" << std::endl;
        }
        else { std::cout << "✅ SFC precision should be sufficient." << std::endl; }
        std::cout << "=======================================" << std::endl;
    }
#endif // !NDEBUG
}

// Log free/total GPU memory at every rank around the cstone sync. Lets us
// pinpoint which rank hit the OOM ceiling and how close the survivors got.
inline void logGpuMemAroundSync(int rank, const char* tag)
{
    size_t freeBytes = 0, totalBytes = 0;
    cudaError_t err = cudaMemGetInfo(&freeBytes, &totalBytes);
    if (err != cudaSuccess) return;
    double freeGiB  = static_cast<double>(freeBytes)  / (1024.0 * 1024.0 * 1024.0);
    double totalGiB = static_cast<double>(totalBytes) / (1024.0 * 1024.0 * 1024.0);
    std::cout << "[GPU mem r" << rank << " " << tag << "] free=" << freeGiB
              << " GiB / " << totalGiB << " GiB" << std::endl;
}

// Element domain constructor
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::ElementDomain(const std::string& meshFile,
                                                                            int rank,
                                                                            int numRanks,
                                                                            bool storeOriginalCoords,
                                                                            int bucketSize,
                                                                            unsigned bucketSizeFocus)
    : rank_(rank)
    , numRanks_(numRanks)
    , box_(0, 1)
    , storeOriginalCoords_(storeOriginalCoords)
{
    if (rank == 0) {
        std::cout << "ElementDomain: Using bucketSize=" << bucketSize
                  << ", bucketSizeFocus=" << bucketSizeFocus << " (passed as parameters)" << std::endl;
    }

    // Host data in SoA format
    HostCoordsTuple h_coords;     // (x, y, z)
    HostConnectivityTuple h_conn; // (i0, i1, i2, ...) depends on element type

    // Device data (will be filled by transferDataToGPU)
    DeviceConnectivityTuple d_conn_;
    DeviceCoordsTuple d_coords_;

    using clk = std::chrono::high_resolution_clock;
    auto t0 = clk::now();

    // Read the mesh in SoA format (host-side: disk + GPU dedup + remap)
    readMeshDataSoA(meshFile, h_coords, h_conn);
    auto t1 = clk::now();
    readMeshTimeMs = std::chrono::duration<float, std::milli>(t1 - t0).count();

    RealType theta = 0.5;

    box_ = computeGlobalBoundingBoxFromCoords<RealType>(std::get<0>(h_coords), std::get<1>(h_coords),
                                                        std::get<2>(h_coords));
    auto t2 = clk::now();
    bboxTimeMs = std::chrono::duration<float, std::milli>(t2 - t1).count();

    std::cout << "Rank " << rank_ << ": Created bounding box: [" << box_.xmin() << "," << box_.xmax() << "] ["
              << box_.ymin() << "," << box_.ymax() << "] [" << box_.zmin() << "," << box_.zmax() << "]" << std::endl;

    // Test SFC precision for this domain (debug only)
    testSfcPrecision<KeyType, RealType>(box_, rank_);

    domain_ = std::make_unique<DomainType>(rank, numRanks, bucketSize, bucketSizeFocus, theta, box_);

    // Transfer data to GPU before sync
    auto t3 = clk::now();
    transferDataToGPU(h_coords, h_conn, d_coords_, d_conn_);
    cudaDeviceSynchronize();
    auto t4 = clk::now();
    h2dTimeMs = std::chrono::duration<float, std::milli>(t4 - t3).count();

    // Calculate characteristic sizes
    calculateCharacteristicSizes(d_conn_, d_coords_);
    cudaDeviceSynchronize();
    auto t5 = clk::now();
    charSizeTimeMs = std::chrono::duration<float, std::milli>(t5 - t4).count();

    std::cout << "Rank " << rank_ << ": Before sync - bucketSize=" << bucketSize
              << ", bucketSizeFocus=" << bucketSizeFocus << ", theta=" << theta << std::endl;

    logGpuMemAroundSync(rank_, "pre-sync");
    sync(d_conn_, d_coords_);
    cudaDeviceSynchronize();
    logGpuMemAroundSync(rank_, "post-sync");
    auto t6 = clk::now();
    syncTimeMs = std::chrono::duration<float, std::milli>(t6 - t5).count();

    std::cout << "Rank " << rank_ << " initialized " << ElementTag::Name << " domain with " << nodeCount_
              << " nodes and " << elementCount_ << " elements." << std::endl;
}

// Constructor from mesh data (for MFEM or other formats) - automatically computes bounding box
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::ElementDomain(const HostCoordsTuple& h_coords,
                                                                            const HostConnectivityTuple& h_conn,
                                                                            int rank,
                                                                            int numRanks,
                                                                            int bucketSize,
                                                                            bool storeOriginalCoords,
                                                                            unsigned bucketSizeFocus)
    : rank_(rank)
    , numRanks_(numRanks)
    , box_(0, 1)
    , storeOriginalCoords_(storeOriginalCoords)
{
    DeviceConnectivityTuple d_conn_;
    DeviceCoordsTuple d_coords_;

    // Set counts from input data
    nodeCount_    = std::get<0>(h_coords).size();
    elementCount_ = std::get<0>(h_conn).size();

    // Compute global bounding box from the provided coordinate data
    box_ = computeGlobalBoundingBoxFromCoords<RealType>(std::get<0>(h_coords), std::get<1>(h_coords),
                                                        std::get<2>(h_coords));

    RealType theta = 0.5;

    std::cout << "Rank " << rank_ << ": Computed bounding box from coordinates: [" << box_.xmin() << "," << box_.xmax()
              << "] [" << box_.ymin() << "," << box_.ymax() << "] [" << box_.zmin() << "," << box_.zmax() << "]"
              << std::endl;

    // Test SFC precision for this domain (debug only)
    testSfcPrecision<KeyType, RealType>(box_, rank_);

    domain_ = std::make_unique<DomainType>(rank, numRanks, bucketSize, bucketSizeFocus, theta, box_);

    // Transfer data to GPU before sync
    transferDataToGPU(h_coords, h_conn, d_coords_, d_conn_);

    // Calculate characteristic sizes
    calculateCharacteristicSizes(d_conn_, d_coords_);

    std::cout << "Rank " << rank_ << ": Before sync - bucketSize=" << bucketSize
              << ", bucketSizeFocus=" << bucketSizeFocus << ", theta=" << theta << ", input nodes=" << nodeCount_
              << ", input elements=" << elementCount_ << std::endl;

    logGpuMemAroundSync(rank_, "pre-sync");
    sync(d_conn_, d_coords_);
    logGpuMemAroundSync(rank_, "post-sync");

    std::cout << "Rank " << rank_ << " initialized " << ElementTag::Name << " domain with " << nodeCount_
              << " nodes and " << elementCount_ << " elements." << std::endl;
}

// Constructor from mesh data with boundary info - automatically computes bounding box
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::ElementDomain(const HostCoordsTuple& h_coords,
                                                                            const HostConnectivityTuple& h_conn,
                                                                            const HostBoundaryTuple& h_boundary,
                                                                            int rank,
                                                                            int numRanks,
                                                                            int bucketSize,
                                                                            unsigned bucketSizeFocus)
    : rank_(rank)
    , numRanks_(numRanks)
    , box_(0, 1)
{
    DeviceConnectivityTuple d_conn_;
    DeviceCoordsTuple d_coords_;

    // Set counts from input data
    nodeCount_    = std::get<0>(h_coords).size();
    elementCount_ = std::get<0>(h_conn).size();

    // Compute global bounding box from the provided coordinate data
    box_ = computeGlobalBoundingBoxFromCoords<RealType>(std::get<0>(h_coords), std::get<1>(h_coords),
                                                        std::get<2>(h_coords));

    RealType theta = 0.5;

    std::cout << "Rank " << rank_ << ": Computed bounding box from coordinates: [" << box_.xmin() << "," << box_.xmax()
              << "] [" << box_.ymin() << "," << box_.ymax() << "] [" << box_.zmin() << "," << box_.zmax() << "]"
              << std::endl;

    // Test SFC precision for this domain (debug only)
    testSfcPrecision<KeyType, RealType>(box_, rank_);

    domain_ = std::make_unique<DomainType>(rank, numRanks, bucketSize, bucketSizeFocus, theta, box_);

    // Transfer data to GPU before sync (including boundary data)
    transferDataToGPU(h_coords, h_conn, h_boundary, d_coords_, d_conn_, d_boundary_);

    // Calculate characteristic sizes
    calculateCharacteristicSizes(d_conn_, d_coords_);

    std::cout << "Rank " << rank_ << ": Before sync - bucketSize=" << bucketSize
              << ", bucketSizeFocus=" << bucketSizeFocus << ", theta=" << theta << ", input nodes=" << nodeCount_
              << ", input elements=" << elementCount_ << std::endl;

    logGpuMemAroundSync(rank_, "pre-sync");
    sync(d_conn_, d_coords_);
    logGpuMemAroundSync(rank_, "post-sync");

    std::cout << "Rank " << rank_ << " initialized " << ElementTag::Name << " domain with " << nodeCount_
              << " nodes and " << elementCount_ << " elements." << std::endl;
}

// Device-data constructor: no host round-trip
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::ElementDomain(DeviceCoordsTuple&& d_coords,
                                                                              DeviceConnectivityTuple&& d_conn,
                                                                              int rank,
                                                                              int numRanks,
                                                                              int bucketSize,
                                                                              unsigned bucketSizeFocus)
    : rank_(rank)
    , numRanks_(numRanks)
    , box_(0, 1)
{
    if constexpr (!std::is_same_v<AcceleratorTag, cstone::GpuTag>)
    {
        throw std::runtime_error("ElementDomain device-data constructor requires GpuTag");
    }

    auto& d_x = std::get<0>(d_coords);
    auto& d_y = std::get<1>(d_coords);
    auto& d_z = std::get<2>(d_coords);

    nodeCount_    = d_x.size();
    elementCount_ = std::get<0>(d_conn).size();

    // Bounding box on GPU + MPI
    RealType lxmin, lxmax, lymin, lymax, lzmin, lzmax;
    {
        auto xb = thrust::device_pointer_cast(d_x.data());
        auto yb = thrust::device_pointer_cast(d_y.data());
        auto zb = thrust::device_pointer_cast(d_z.data());
        lxmin = thrust::reduce(xb, xb + nodeCount_, std::numeric_limits<RealType>::max(),    thrust::minimum<RealType>());
        lxmax = thrust::reduce(xb, xb + nodeCount_, std::numeric_limits<RealType>::lowest(), thrust::maximum<RealType>());
        lymin = thrust::reduce(yb, yb + nodeCount_, std::numeric_limits<RealType>::max(),    thrust::minimum<RealType>());
        lymax = thrust::reduce(yb, yb + nodeCount_, std::numeric_limits<RealType>::lowest(), thrust::maximum<RealType>());
        lzmin = thrust::reduce(zb, zb + nodeCount_, std::numeric_limits<RealType>::max(),    thrust::minimum<RealType>());
        lzmax = thrust::reduce(zb, zb + nodeCount_, std::numeric_limits<RealType>::lowest(), thrust::maximum<RealType>());
    }

    auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    RealType gxmin, gxmax, gymin, gymax, gzmin, gzmax;
    MPI_Allreduce(&lxmin, &gxmin, 1, mpiType, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&lxmax, &gxmax, 1, mpiType, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&lymin, &gymin, 1, mpiType, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&lymax, &gymax, 1, mpiType, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&lzmin, &gzmin, 1, mpiType, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&lzmax, &gzmax, 1, mpiType, MPI_MAX, MPI_COMM_WORLD);

    RealType pad   = RealType(0.05);
    RealType rx    = gxmax - gxmin, ry = gymax - gymin, rz = gzmax - gzmin;
    box_ = cstone::Box<RealType>(gxmin - pad * rx, gxmax + pad * rx,
                                  gymin - pad * ry, gymax + pad * ry,
                                  gzmin - pad * rz, gzmax + pad * rz);

    RealType theta = 0.5;

    domain_ = std::make_unique<DomainType>(rank, numRanks, bucketSize, bucketSizeFocus, theta, box_);

    // d_props_ is initialized in calculateCharacteristicSizes; no transferDataToGPU needed
    DeviceConnectivityTuple d_conn_local = std::move(d_conn);
    DeviceCoordsTuple       d_coords_local = std::move(d_coords);
    calculateCharacteristicSizes(d_conn_local, d_coords_local);

    logGpuMemAroundSync(rank_, "pre-sync");
    sync(d_conn_local, d_coords_local);
    logGpuMemAroundSync(rank_, "post-sync");

    std::cout << "Rank " << rank_ << " initialized " << ElementTag::Name << " domain (device ctor) with "
              << nodeCount_ << " nodes and " << elementCount_ << " elements." << std::endl;
}

// Read mesh data in SoA format; uses element-based partitioning for better data locality
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::readMeshDataSoA(const std::string& meshFile,
                                                                                   HostCoordsTuple& h_coords_,
                                                                                   HostConnectivityTuple& h_conn_)
{
    int rank = rank_;
    int numRanks = numRanks_;
    std::cout << "Reading mesh data from " << meshFile << " on rank " << rank << std::endl;
    try
    {
        // Helper to handle coordinate conversion based on type
        // Now uses RealType directly since mesh readers use RealType for precision
        auto processCoordinates =
            [&](std::vector<RealType>& x_data, std::vector<RealType>& y_data, std::vector<RealType>& z_data)
        {
            auto& h_x = std::get<0>(h_coords_);
            auto& h_y = std::get<1>(h_coords_);
            auto& h_z = std::get<2>(h_coords_);

            // Direct move since mesh readers now use RealType
            h_x = std::move(x_data);
            h_y = std::move(y_data);
            h_z = std::move(z_data);
        };

        // Detect mesh file format based on extension
        // - .mesh = MFEM format
        // - .exo  = Exodus II format (netCDF)
        // - directory = binary mesh format
        namespace fs      = std::filesystem;
        bool isMFEMFile   = !fs::is_directory(meshFile) && meshFile.find(".mesh") != std::string::npos;
        bool isExodusFile = !fs::is_directory(meshFile) &&
                            (meshFile.find(".exo") != std::string::npos ||
                             meshFile.find(".e") != std::string::npos);

        // Use compile-time branching to select element type
        if constexpr (std::is_same_v<ElementTag, TetTag>)
        {
            if (isExodusFile)
            {
                // Use Exodus mesh reader for tet elements
                auto [readNodeCount, readElementCount, x_data, y_data, z_data, conn_tuple, localToGlobal,
                      boundaryNodes] =
                    mars::readExodusMeshWithElementPartitioning<4, RealType, KeyType>(meshFile, rank, numRanks);

                nodeCount_    = readNodeCount;
                elementCount_ = readElementCount;
                processCoordinates(x_data, y_data, z_data);
                h_conn_ = std::move(conn_tuple);

                // Store the global node indices
                d_localToGlobalNodeMap_.resize(localToGlobal.size());
                thrust::copy(localToGlobal.begin(), localToGlobal.end(),
                             thrust::device_pointer_cast(d_localToGlobalNodeMap_.data()));
            }
            else if (isMFEMFile)
            {
                // Use MFEM mesh reader
                auto [readNodeCount, readElementCount, x_data, y_data, z_data, conn_tuple, localToGlobal,
                      boundaryNodes] =
                    mars::readMFEMMeshWithElementPartitioning<4, RealType, KeyType>(meshFile, rank, numRanks);

                nodeCount_    = readNodeCount;
                elementCount_ = readElementCount;
                processCoordinates(x_data, y_data, z_data);
                h_conn_ = std::move(conn_tuple);

                // Store the global node indices
                d_localToGlobalNodeMap_.resize(localToGlobal.size());
                thrust::copy(localToGlobal.begin(), localToGlobal.end(),
                             thrust::device_pointer_cast(d_localToGlobalNodeMap_.data()));
            }
            else
            {
                // Use binary mesh reader
                auto [readNodeCount, readElementCount, x_data, y_data, z_data, conn_tuple, localToGlobal] =
                    mars::readMeshWithElementPartitioning<4, RealType, KeyType>(meshFile, rank, numRanks);

                nodeCount_    = readNodeCount;
                elementCount_ = readElementCount;
                processCoordinates(x_data, y_data, z_data);
                h_conn_ = std::move(conn_tuple);

                // Store the global node indices
                d_localToGlobalNodeMap_.resize(localToGlobal.size());
                thrust::copy(localToGlobal.begin(), localToGlobal.end(),
                             thrust::device_pointer_cast(d_localToGlobalNodeMap_.data()));
            }
        }
        else if constexpr (std::is_same_v<ElementTag, HexTag>)
        {
            if (isExodusFile)
            {
                // Use Exodus mesh reader for hex elements
                auto [readNodeCount, readElementCount, x_data, y_data, z_data, conn_tuple, localToGlobal,
                      boundaryNodes] =
                    mars::readExodusMeshWithElementPartitioning<8, RealType, KeyType>(meshFile, rank, numRanks);

                nodeCount_    = readNodeCount;
                elementCount_ = readElementCount;
                processCoordinates(x_data, y_data, z_data);
                h_conn_ = std::move(conn_tuple);

                // Store the global node indices
                d_localToGlobalNodeMap_.resize(localToGlobal.size());
                thrust::copy(localToGlobal.begin(), localToGlobal.end(),
                             thrust::device_pointer_cast(d_localToGlobalNodeMap_.data()));
            }
            else if (isMFEMFile)
            {
                // Use MFEM mesh reader for hex elements
                auto [readNodeCount, readElementCount, x_data, y_data, z_data, conn_tuple, localToGlobal,
                      boundaryNodes] =
                    mars::readMFEMMeshWithElementPartitioning<8, RealType, KeyType>(meshFile, rank, numRanks);

                nodeCount_    = readNodeCount;
                elementCount_ = readElementCount;
                processCoordinates(x_data, y_data, z_data);
                h_conn_ = std::move(conn_tuple);

                // Store the global node indices
                d_localToGlobalNodeMap_.resize(localToGlobal.size());
                thrust::copy(localToGlobal.begin(), localToGlobal.end(),
                             thrust::device_pointer_cast(d_localToGlobalNodeMap_.data()));
            }
            else
            {
                // Use binary mesh reader
                auto [readNodeCount, readElementCount, x_data, y_data, z_data, conn_tuple, localToGlobal] =
                    mars::readMeshWithElementPartitioning<8, RealType, KeyType>(meshFile, rank, numRanks);

                nodeCount_    = readNodeCount;
                elementCount_ = readElementCount;
                processCoordinates(x_data, y_data, z_data);
                h_conn_ = std::move(conn_tuple);

                // Store the global node indices
                d_localToGlobalNodeMap_.resize(localToGlobal.size());
                thrust::copy(localToGlobal.begin(), localToGlobal.end(),
                             thrust::device_pointer_cast(d_localToGlobalNodeMap_.data()));
            }
        }
        else if constexpr (std::is_same_v<ElementTag, TriTag>)
        {
            if (isMFEMFile)
            {
                // Use MFEM mesh reader for triangle elements
                // NOTE: Use RealType for coordinate precision (not float) to preserve exact values
                auto [readNodeCount, readElementCount, x_data, y_data, z_data, conn_tuple, localToGlobal,
                      boundaryNodes] =
                    mars::readMFEMMeshWithElementPartitioning<3, RealType, KeyType>(meshFile, rank, numRanks);

                nodeCount_    = readNodeCount;
                elementCount_ = readElementCount;
                processCoordinates(x_data, y_data, z_data);
                h_conn_ = std::move(conn_tuple);

                // Store the global node indices
                d_localToGlobalNodeMap_.resize(localToGlobal.size());
                thrust::copy(localToGlobal.begin(), localToGlobal.end(),
                             thrust::device_pointer_cast(d_localToGlobalNodeMap_.data()));
            }
            else
            {
                // Use binary mesh reader
                // NOTE: Use RealType for coordinate precision (not float) to preserve exact values
                auto [readNodeCount, readElementCount, x_data, y_data, z_data, conn_tuple, localToGlobal] =
                    mars::readMeshWithElementPartitioning<3, RealType, KeyType>(meshFile, rank, numRanks);

                nodeCount_    = readNodeCount;
                elementCount_ = readElementCount;
                processCoordinates(x_data, y_data, z_data);
                h_conn_ = std::move(conn_tuple);

                // Store the global node indices
                d_localToGlobalNodeMap_.resize(localToGlobal.size());
                thrust::copy(localToGlobal.begin(), localToGlobal.end(),
                             thrust::device_pointer_cast(d_localToGlobalNodeMap_.data()));
            }
        }
        else if constexpr (std::is_same_v<ElementTag, QuadTag>)
        {
            if (isMFEMFile)
            {
                // Use MFEM mesh reader for quad elements
                // NOTE: Use RealType for coordinate precision (not float) to preserve exact values
                auto [readNodeCount, readElementCount, x_data, y_data, z_data, conn_tuple, localToGlobal,
                      boundaryNodes] =
                    mars::readMFEMMeshWithElementPartitioning<4, RealType, KeyType>(meshFile, rank, numRanks);

                nodeCount_    = readNodeCount;
                elementCount_ = readElementCount;
                processCoordinates(x_data, y_data, z_data);
                h_conn_ = std::move(conn_tuple);

                // Store the global node indices
                d_localToGlobalNodeMap_.resize(localToGlobal.size());
                thrust::copy(localToGlobal.begin(), localToGlobal.end(),
                             thrust::device_pointer_cast(d_localToGlobalNodeMap_.data()));
            }
            else
            {
                // Use binary mesh reader
                // NOTE: Use RealType for coordinate precision (not float) to preserve exact values
                auto [readNodeCount, readElementCount, x_data, y_data, z_data, conn_tuple, localToGlobal] =
                    mars::readMeshWithElementPartitioning<4, RealType, KeyType>(meshFile, rank, numRanks);

                nodeCount_    = readNodeCount;
                elementCount_ = readElementCount;
                processCoordinates(x_data, y_data, z_data);
                h_conn_ = std::move(conn_tuple);

                // Store the global node indices
                d_localToGlobalNodeMap_.resize(localToGlobal.size());
                thrust::copy(localToGlobal.begin(), localToGlobal.end(),
                             thrust::device_pointer_cast(d_localToGlobalNodeMap_.data()));
            }
        }
        else
        {
            constexpr bool dependent_false = sizeof(ElementTag) != sizeof(ElementTag);
            static_assert(dependent_false, "Unsupported element type in readMeshDataSoA");
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error reading mesh from " << meshFile << ": " << e.what() << std::endl;
        throw;
    }
}

// Calculate characteristic sizes for elements - GPU-only version
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::calculateCharacteristicSizes(
    const DeviceConnectivityTuple& d_conn_, const DeviceCoordsTuple& d_coords_)
{
    // Work directly with device data
    auto& d_h = std::get<0>(d_props_);
    d_h.resize(nodeCount_);

    // Compute per-element h from actual edge lengths for Tet and Hex elements.
    // Cornerstone uses h as an SPH smoothing length (interaction radius = 2*h),
    // so h should reflect the local mesh size to avoid over-fetching ghost elements.
    if constexpr (std::is_same_v<AcceleratorTag, cstone::GpuTag> &&
                  (std::is_same_v<ElementTag, TetTag> || std::is_same_v<ElementTag, HexTag>))
    {
        DeviceVector<int> d_nodeTetCount(nodeCount_, 0);

        auto& d_x  = std::get<0>(d_coords_);
        auto& d_y  = std::get<1>(d_coords_);
        auto& d_z  = std::get<2>(d_coords_);
        auto& d_i0 = std::get<0>(d_conn_);
        auto& d_i1 = std::get<1>(d_conn_);
        auto& d_i2 = std::get<2>(d_conn_);
        auto& d_i3 = std::get<3>(d_conn_);

        constexpr int NodesPerElem = ElementTag::NodesPerElement;

        const KeyType* d_i4_ptr;
        const KeyType* d_i5_ptr;
        const KeyType* d_i6_ptr;
        const KeyType* d_i7_ptr;

        if constexpr (NodesPerElem > 4) { d_i4_ptr = thrust::raw_pointer_cast(std::get<4>(d_conn_).data()); }
        else { d_i4_ptr = nullptr; }

        if constexpr (NodesPerElem > 5) { d_i5_ptr = thrust::raw_pointer_cast(std::get<5>(d_conn_).data()); }
        else { d_i5_ptr = nullptr; }

        if constexpr (NodesPerElem > 6) { d_i6_ptr = thrust::raw_pointer_cast(std::get<6>(d_conn_).data()); }
        else { d_i6_ptr = nullptr; }

        if constexpr (NodesPerElem > 7) { d_i7_ptr = thrust::raw_pointer_cast(std::get<7>(d_conn_).data()); }
        else { d_i7_ptr = nullptr; }

        // accumulate edge lengths per node
        int blockSize = 256;
        int numBlocks = (elementCount_ + blockSize - 1) / blockSize;

        computeCharacteristicSizesKernel<ElementTag, KeyType, RealType><<<numBlocks, blockSize>>>(
            thrust::raw_pointer_cast(d_x.data()), thrust::raw_pointer_cast(d_y.data()),
            thrust::raw_pointer_cast(d_z.data()), thrust::raw_pointer_cast(d_i0.data()),
            thrust::raw_pointer_cast(d_i1.data()), thrust::raw_pointer_cast(d_i2.data()),
            thrust::raw_pointer_cast(d_i3.data()), d_i4_ptr, d_i5_ptr, d_i6_ptr, d_i7_ptr,
            thrust::raw_pointer_cast(d_nodeTetCount.data()), thrust::raw_pointer_cast(d_h.data()), elementCount_);

        cudaCheckError();

        // normalize by number of contributions
        numBlocks = (nodeCount_ + blockSize - 1) / blockSize;
        finalizeCharacteristicSizesKernel<KeyType, RealType><<<numBlocks, blockSize>>>(
            thrust::raw_pointer_cast(d_h.data()), thrust::raw_pointer_cast(d_nodeTetCount.data()), nodeCount_);

        cudaCheckError();

        constexpr RealType meshFactor = 1.0;
        constexpr RealType minH       = 1.0e-6;
        constexpr RealType maxH       = 1.0;

        transformCharacteristicSizesKernel<RealType>
            <<<numBlocks, blockSize>>>(thrust::raw_pointer_cast(d_h.data()), nodeCount_, meshFactor, minH, maxH);
        cudaCheckError();
    }
    else
    {
        // Fallback for Tri/Quad elements or CPU path: estimate h from mesh density.
        // Use cbrt(volume/count) as characteristic element size; factor 0.5 so 2h ~ element size.
        RealType domainVolume = (domain_->box().xmax() - domain_->box().xmin()) *
                                (domain_->box().ymax() - domain_->box().ymin()) *
                                (domain_->box().zmax() - domain_->box().zmin());
        RealType defaultH = RealType(0.5) * std::cbrt(domainVolume / static_cast<RealType>(elementCount_));

        int blockSize = 256;
        int numBlocks = (nodeCount_ + blockSize - 1) / blockSize;
        fillCharacteristicSizesKernel<RealType>
            <<<numBlocks, blockSize>>>(thrust::raw_pointer_cast(d_h.data()), nodeCount_, defaultH);
        cudaCheckError();
    }
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::sync(const DeviceConnectivityTuple& d_conn_,
                                                                        const DeviceCoordsTuple& d_coords_)
{
    if constexpr (std::is_same_v<AcceleratorTag, cstone::GpuTag>)
    {
        // Extract coordinates for clarity
        auto& d_x = std::get<0>(d_coords_);
        auto& d_y = std::get<1>(d_coords_);
        auto& d_z = std::get<2>(d_coords_);
        auto& d_h = std::get<0>(d_props_);

        // Calculate SFC codes for nodes
        DeviceVector<KeyType> d_nodeSfcCodes(nodeCount_);
        generateSfcKeys<KeyType, RealType>(thrust::raw_pointer_cast(d_x.data()), thrust::raw_pointer_cast(d_y.data()),
                                           thrust::raw_pointer_cast(d_z.data()),
                                           thrust::raw_pointer_cast(d_nodeSfcCodes.data()), nodeCount_,
                                           getBoundingBox());

        // Find representative nodes for each element
        d_elemToNodeMap_.resize(elementCount_);

        // initialize connectivity keys storage
        initializeConnectivityKeys();

        int blockSize = 256;
        int numBlocks = (elementCount_ + blockSize - 1) / blockSize;

        // Call element-specific implementation
        syncImpl(d_nodeSfcCodes, d_conn_, blockSize, numBlocks);

        // Create element arrays
        DeviceVector<RealType> d_elemX(elementCount_);
        DeviceVector<RealType> d_elemY(elementCount_);
        DeviceVector<RealType> d_elemZ(elementCount_);
        DeviceVector<RealType> d_elemH(elementCount_);
        d_elemSfcCodes_.resize(elementCount_);

        // Extract element coordinates
        extractRepCoordinatesKernel<ElementTag, KeyType, RealType><<<numBlocks, blockSize>>>(
            thrust::raw_pointer_cast(d_x.data()), thrust::raw_pointer_cast(d_y.data()),
            thrust::raw_pointer_cast(d_z.data()), thrust::raw_pointer_cast(d_h.data()),
            thrust::raw_pointer_cast(d_elemToNodeMap_.data()), thrust::raw_pointer_cast(d_elemX.data()),
            thrust::raw_pointer_cast(d_elemY.data()), thrust::raw_pointer_cast(d_elemZ.data()),
            thrust::raw_pointer_cast(d_elemH.data()), elementCount_);
        cudaCheckError();

        std::cout << "Rank " << rank_ << " syncing " << ElementTag::Name << " domain with " << elementCount_
                  << " elements." << std::endl;

        // Check if we need to sync with original coordinates (compile-time check for hex8)
        if constexpr (std::is_same_v<ElementTag, HexTag>) {
            if (storeOriginalCoords_) {
                // Extract original coordinates for all 8 nodes of each element (24 arrays total)
                DeviceVector<RealType> d_elemOrigX0(elementCount_), d_elemOrigY0(elementCount_), d_elemOrigZ0(elementCount_);
                DeviceVector<RealType> d_elemOrigX1(elementCount_), d_elemOrigY1(elementCount_), d_elemOrigZ1(elementCount_);
                DeviceVector<RealType> d_elemOrigX2(elementCount_), d_elemOrigY2(elementCount_), d_elemOrigZ2(elementCount_);
                DeviceVector<RealType> d_elemOrigX3(elementCount_), d_elemOrigY3(elementCount_), d_elemOrigZ3(elementCount_);
                DeviceVector<RealType> d_elemOrigX4(elementCount_), d_elemOrigY4(elementCount_), d_elemOrigZ4(elementCount_);
                DeviceVector<RealType> d_elemOrigX5(elementCount_), d_elemOrigY5(elementCount_), d_elemOrigZ5(elementCount_);
                DeviceVector<RealType> d_elemOrigX6(elementCount_), d_elemOrigY6(elementCount_), d_elemOrigZ6(elementCount_);
                DeviceVector<RealType> d_elemOrigX7(elementCount_), d_elemOrigY7(elementCount_), d_elemOrigZ7(elementCount_);

                // Extract original coordinates from input arrays
                auto& d_i0 = std::get<0>(d_conn_);
                auto& d_i1 = std::get<1>(d_conn_);
                auto& d_i2 = std::get<2>(d_conn_);
                auto& d_i3 = std::get<3>(d_conn_);
                auto& d_i4 = std::get<4>(d_conn_);
                auto& d_i5 = std::get<5>(d_conn_);
                auto& d_i6 = std::get<6>(d_conn_);
                auto& d_i7 = std::get<7>(d_conn_);


            extractElementNodeCoordsKernel<ElementTag, KeyType, RealType><<<numBlocks, blockSize>>>(
                thrust::raw_pointer_cast(d_x.data()), thrust::raw_pointer_cast(d_y.data()), thrust::raw_pointer_cast(d_z.data()),
                thrust::raw_pointer_cast(d_i0.data()), thrust::raw_pointer_cast(d_i1.data()),
                thrust::raw_pointer_cast(d_i2.data()), thrust::raw_pointer_cast(d_i3.data()),
                thrust::raw_pointer_cast(d_i4.data()), thrust::raw_pointer_cast(d_i5.data()),
                thrust::raw_pointer_cast(d_i6.data()), thrust::raw_pointer_cast(d_i7.data()),
                thrust::raw_pointer_cast(d_elemOrigX0.data()), thrust::raw_pointer_cast(d_elemOrigY0.data()), thrust::raw_pointer_cast(d_elemOrigZ0.data()),
                thrust::raw_pointer_cast(d_elemOrigX1.data()), thrust::raw_pointer_cast(d_elemOrigY1.data()), thrust::raw_pointer_cast(d_elemOrigZ1.data()),
                thrust::raw_pointer_cast(d_elemOrigX2.data()), thrust::raw_pointer_cast(d_elemOrigY2.data()), thrust::raw_pointer_cast(d_elemOrigZ2.data()),
                thrust::raw_pointer_cast(d_elemOrigX3.data()), thrust::raw_pointer_cast(d_elemOrigY3.data()), thrust::raw_pointer_cast(d_elemOrigZ3.data()),
                thrust::raw_pointer_cast(d_elemOrigX4.data()), thrust::raw_pointer_cast(d_elemOrigY4.data()), thrust::raw_pointer_cast(d_elemOrigZ4.data()),
                thrust::raw_pointer_cast(d_elemOrigX5.data()), thrust::raw_pointer_cast(d_elemOrigY5.data()), thrust::raw_pointer_cast(d_elemOrigZ5.data()),
                thrust::raw_pointer_cast(d_elemOrigX6.data()), thrust::raw_pointer_cast(d_elemOrigY6.data()), thrust::raw_pointer_cast(d_elemOrigZ6.data()),
                thrust::raw_pointer_cast(d_elemOrigX7.data()), thrust::raw_pointer_cast(d_elemOrigY7.data()), thrust::raw_pointer_cast(d_elemOrigZ7.data()),
                elementCount_);
            cudaCheckError();

            // Package into tuple for sync (direct tie, not std::ref)
            auto d_orig_coords = std::tie(
                d_elemOrigX0, d_elemOrigY0, d_elemOrigZ0,
                d_elemOrigX1, d_elemOrigY1, d_elemOrigZ1,
                d_elemOrigX2, d_elemOrigY2, d_elemOrigZ2,
                d_elemOrigX3, d_elemOrigY3, d_elemOrigZ3,
                d_elemOrigX4, d_elemOrigY4, d_elemOrigZ4,
                d_elemOrigX5, d_elemOrigY5, d_elemOrigZ5,
                d_elemOrigX6, d_elemOrigY6, d_elemOrigZ6,
                d_elemOrigX7, d_elemOrigY7, d_elemOrigZ7
            );

            // Sync with original coordinates
            syncDomainImplWithOrigCoords(domain_.get(), d_elemSfcCodes_, d_elemX, d_elemY, d_elemZ, d_elemH,
                                        elementCount_, d_conn_keys_, d_orig_coords);

            // Build adjacency to get local IDs (d_conn_keys_ -> d_conn_local_ids_)
            // This is needed because d_conn_keys_ contains SFC keys, not local node indices
            // IMPORTANT: Must call BEFORE allocating originalCoords_ because ensureAdjacency()
            // updates nodeCount_ via createLocalToGlobalSfcMap()
            ensureAdjacency();

            // Allocate node coordinate arrays using updated nodeCount_ (after ensureAdjacency)
            originalCoords_ = std::make_unique<OriginalCoordinates<ElementTag, RealType, KeyType, AcceleratorTag>>(nodeCount_);

            // Rebuild node arrays using post-sync LOCAL connectivity (d_conn_local_ids_)
            int newNumBlocks = (elementCount_ + blockSize - 1) / blockSize;
            rebuildNodeCoordsFromElementsKernel<ElementTag, KeyType, RealType><<<newNumBlocks, blockSize>>>(
                thrust::raw_pointer_cast(std::get<0>(adjacency_->d_conn_local_ids_).data()),
                thrust::raw_pointer_cast(std::get<1>(adjacency_->d_conn_local_ids_).data()),
                thrust::raw_pointer_cast(std::get<2>(adjacency_->d_conn_local_ids_).data()),
                thrust::raw_pointer_cast(std::get<3>(adjacency_->d_conn_local_ids_).data()),
                thrust::raw_pointer_cast(std::get<4>(adjacency_->d_conn_local_ids_).data()),
                thrust::raw_pointer_cast(std::get<5>(adjacency_->d_conn_local_ids_).data()),
                thrust::raw_pointer_cast(std::get<6>(adjacency_->d_conn_local_ids_).data()),
                thrust::raw_pointer_cast(std::get<7>(adjacency_->d_conn_local_ids_).data()),
                thrust::raw_pointer_cast(d_elemOrigX0.data()), thrust::raw_pointer_cast(d_elemOrigY0.data()), thrust::raw_pointer_cast(d_elemOrigZ0.data()),
                thrust::raw_pointer_cast(d_elemOrigX1.data()), thrust::raw_pointer_cast(d_elemOrigY1.data()), thrust::raw_pointer_cast(d_elemOrigZ1.data()),
                thrust::raw_pointer_cast(d_elemOrigX2.data()), thrust::raw_pointer_cast(d_elemOrigY2.data()), thrust::raw_pointer_cast(d_elemOrigZ2.data()),
                thrust::raw_pointer_cast(d_elemOrigX3.data()), thrust::raw_pointer_cast(d_elemOrigY3.data()), thrust::raw_pointer_cast(d_elemOrigZ3.data()),
                thrust::raw_pointer_cast(d_elemOrigX4.data()), thrust::raw_pointer_cast(d_elemOrigY4.data()), thrust::raw_pointer_cast(d_elemOrigZ4.data()),
                thrust::raw_pointer_cast(d_elemOrigX5.data()), thrust::raw_pointer_cast(d_elemOrigY5.data()), thrust::raw_pointer_cast(d_elemOrigZ5.data()),
                thrust::raw_pointer_cast(d_elemOrigX6.data()), thrust::raw_pointer_cast(d_elemOrigY6.data()), thrust::raw_pointer_cast(d_elemOrigZ6.data()),
                thrust::raw_pointer_cast(d_elemOrigX7.data()), thrust::raw_pointer_cast(d_elemOrigY7.data()), thrust::raw_pointer_cast(d_elemOrigZ7.data()),
                thrust::raw_pointer_cast(originalCoords_->d_node_x_.data()),
                thrust::raw_pointer_cast(originalCoords_->d_node_y_.data()),
                thrust::raw_pointer_cast(originalCoords_->d_node_z_.data()),
                elementCount_);
            cudaCheckError();
            } else {
                // Standard sync with SFC-decoded coordinates
                syncDomainImpl(domain_.get(), d_elemSfcCodes_, d_elemX, d_elemY, d_elemZ, d_elemH, elementCount_, d_conn_keys_);
            }
        } else {
            // For non-hex elements, always use standard sync
            syncDomainImpl(domain_.get(), d_elemSfcCodes_, d_elemX, d_elemY, d_elemZ, d_elemH, elementCount_, d_conn_keys_);
        }
    }
}

// Transfer data to GPU for computations
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::transferDataToGPU(
    const HostCoordsTuple& h_coords_,
    const HostConnectivityTuple& h_conn_,
    DeviceCoordsTuple& d_coords_,
    DeviceConnectivityTuple& d_conn_)
{
    if constexpr (std::is_same_v<AcceleratorTag, cstone::GpuTag>)
    {
        // Copy coordinate data from host to device
        copyTupleElements(d_coords_, h_coords_);

        // Copy connectivity data from host to device
        copyTupleElements(d_conn_, h_conn_);

        // Initialize properties (h values)
        d_props_ = std::tuple<DeviceVector<RealType>>(DeviceVector<RealType>(nodeCount_));
    }
}

// Transfer data to GPU for computations (including boundary info)
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::transferDataToGPU(
    const HostCoordsTuple& h_coords_,
    const HostConnectivityTuple& h_conn_,
    const HostBoundaryTuple& h_boundary_,
    DeviceCoordsTuple& d_coords_,
    DeviceConnectivityTuple& d_conn_,
    DeviceBoundaryTuple& d_boundary_)
{
    if constexpr (std::is_same_v<AcceleratorTag, cstone::GpuTag>)
    {
        // Copy coordinate data from host to device
        copyTupleElements(d_coords_, h_coords_);

        // Copy connectivity data from host to device
        copyTupleElements(d_conn_, h_conn_);

        // Copy boundary data from host to device
        std::get<0>(d_boundary_).resize(std::get<0>(h_boundary_).size());
        cudaMemcpy(std::get<0>(d_boundary_).data(), std::get<0>(h_boundary_).data(),
                   std::get<0>(h_boundary_).size() * sizeof(uint8_t), cudaMemcpyHostToDevice);

        // Initialize properties (h values)
        d_props_ = std::tuple<DeviceVector<RealType>>(DeviceVector<RealType>(nodeCount_));
    }
}

// Implementation of SFC-to-coordinate functions
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
MARS_HOST_DEVICE std::tuple<RealType, RealType, RealType>
ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::sfcToPhysicalCoordinate(KeyType sfcKey) const
{
    return decodeSfcToPhysical(sfcKey, getBoundingBox());
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
MARS_HOST_DEVICE RealType
ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::sfcToPhysicalCoordinateX(KeyType sfcKey) const
{
    // Convert raw key to SfcKind strong type
    auto sfcKindKey   = cstone::SfcKind<KeyType>(sfcKey);
    auto [ix, iy, iz] = cstone::decodeSfc(sfcKindKey);

    // Use SfcKind for maxTreeLevel
    constexpr unsigned maxCoord = (1u << cstone::maxTreeLevel<cstone::SfcKind<KeyType>>{}) - 1;
    RealType invMaxCoord        = RealType(1.0) / maxCoord;
    auto box                    = getBoundingBox();

    return box.xmin() + ix * invMaxCoord * (box.xmax() - box.xmin());
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
MARS_HOST_DEVICE RealType
ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::sfcToPhysicalCoordinateY(KeyType sfcKey) const
{
    // Convert raw key to SfcKind strong type
    auto sfcKindKey   = cstone::SfcKind<KeyType>(sfcKey);
    auto [ix, iy, iz] = cstone::decodeSfc(sfcKindKey);

    // Use SfcKind for maxTreeLevel
    constexpr unsigned maxCoord = (1u << cstone::maxTreeLevel<cstone::SfcKind<KeyType>>{}) - 1;
    RealType invMaxCoord        = RealType(1.0) / maxCoord;
    auto box                    = getBoundingBox();

    return box.ymin() + iy * invMaxCoord * (box.ymax() - box.ymin());
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
MARS_HOST_DEVICE RealType
ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::sfcToPhysicalCoordinateZ(KeyType sfcKey) const
{
    // Convert raw key to SfcKind strong type
    auto sfcKindKey   = cstone::SfcKind<KeyType>(sfcKey);
    auto [ix, iy, iz] = cstone::decodeSfc(sfcKindKey);

    // Use SfcKind for maxTreeLevel
    constexpr unsigned maxCoord = (1u << cstone::maxTreeLevel<cstone::SfcKind<KeyType>>{}) - 1;
    RealType invMaxCoord        = RealType(1.0) / maxCoord;
    auto box                    = getBoundingBox();

    return box.zmin() + iz * invMaxCoord * (box.zmax() - box.zmin());
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
MARS_HOST_DEVICE std::tuple<unsigned, unsigned, unsigned>
ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::sfcToSpatialCoordinate(KeyType sfcKey) const
{
    // Convert raw key to SfcKind strong type
    auto sfcKindKey   = cstone::SfcKind<KeyType>(sfcKey);
    auto [ix, iy, iz] = cstone::decodeSfc(sfcKindKey);
    return std::make_tuple(ix, iy, iz);
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
MARS_HOST_DEVICE unsigned
ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::sfcToSpatialCoordinateX(KeyType sfcKey) const
{
    // Convert raw key to SfcKind strong type
    auto sfcKindKey   = cstone::SfcKind<KeyType>(sfcKey);
    auto [ix, iy, iz] = cstone::decodeSfc(sfcKindKey);
    return ix;
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
MARS_HOST_DEVICE unsigned
ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::sfcToSpatialCoordinateY(KeyType sfcKey) const
{
    // Convert raw key to SfcKind strong type
    auto sfcKindKey   = cstone::SfcKind<KeyType>(sfcKey);
    auto [ix, iy, iz] = cstone::decodeSfc(sfcKindKey);
    return iy;
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
MARS_HOST_DEVICE unsigned
ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::sfcToSpatialCoordinateZ(KeyType sfcKey) const
{
    // Convert raw key to SfcKind strong type
    auto sfcKindKey   = cstone::SfcKind<KeyType>(sfcKey);
    auto [ix, iy, iz] = cstone::decodeSfc(sfcKindKey);
    return iz;
}

// host function to get connectivity for a specific element index
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
template<int I>
KeyType ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::getConnectivity(size_t elementIndex) const
{
    static_assert(I >= 0 && I < NodesPerElement, "Index out of range for element connectivity");

    if (elementIndex >= elementCount_) return KeyType(0);

    // Direct access to SFC key - no host transfer
    return std::get<I>(d_conn_keys_)[elementIndex];
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
template<int I>
typename ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::template DeviceVector<KeyType>
ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::getConnectivity() const
{
    static_assert(I >= 0 && I < NodesPerElement, "Index out of range for element connectivity");
    return std::get<I>(d_conn_keys_);
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
typename ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::DeviceConnectivityTuple
ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::getConnectivity() const
{
    return d_conn_keys_;
}

// Initialize d_conn_keys_ in constructor - missing resize calls
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::initializeConnectivityKeys()
{
    // Initialize all connectivity key vectors based on element type
    if constexpr (std::is_same_v<ElementTag, TetTag>)
    {
        std::get<0>(d_conn_keys_).resize(elementCount_);
        std::get<1>(d_conn_keys_).resize(elementCount_);
        std::get<2>(d_conn_keys_).resize(elementCount_);
        std::get<3>(d_conn_keys_).resize(elementCount_);
    }
    else if constexpr (std::is_same_v<ElementTag, HexTag>)
    {
        std::get<0>(d_conn_keys_).resize(elementCount_);
        std::get<1>(d_conn_keys_).resize(elementCount_);
        std::get<2>(d_conn_keys_).resize(elementCount_);
        std::get<3>(d_conn_keys_).resize(elementCount_);
        std::get<4>(d_conn_keys_).resize(elementCount_);
        std::get<5>(d_conn_keys_).resize(elementCount_);
        std::get<6>(d_conn_keys_).resize(elementCount_);
        std::get<7>(d_conn_keys_).resize(elementCount_);
    }
    else if constexpr (std::is_same_v<ElementTag, TriTag>)
    {
        std::get<0>(d_conn_keys_).resize(elementCount_);
        std::get<1>(d_conn_keys_).resize(elementCount_);
        std::get<2>(d_conn_keys_).resize(elementCount_);
    }
    else if constexpr (std::is_same_v<ElementTag, QuadTag>)
    {
        std::get<0>(d_conn_keys_).resize(elementCount_);
        std::get<1>(d_conn_keys_).resize(elementCount_);
        std::get<2>(d_conn_keys_).resize(elementCount_);
        std::get<3>(d_conn_keys_).resize(elementCount_);
    }
}

// ============================================================================
// Lazy initialization methods
// ============================================================================
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::ensureSfcMap() const
{
    if (d_localToGlobalSfcMap_.empty())
    {
        // Need to call non-const version, so cast away const
        const_cast<ElementDomain*>(this)->createLocalToGlobalSfcMap();
    }
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::ensureAdjacency() const
{
    if (!adjacency_)
    {
        // Adjacency needs the SFC map to be built first
        ensureSfcMap();
        adjacency_ = std::make_unique<AdjacencyData<ElementTag, RealType, KeyType, AcceleratorTag>>(*this);
    }
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::ensureHalo() const
{
    if (!halo_)
    {
        // Halo data needs the SFC map to be built first (for correct nodeCount)
        ensureSfcMap();
        // HaloData now works directly with SFC keys - no adjacency needed!
        halo_ = std::make_unique<HaloData<ElementTag, RealType, KeyType, AcceleratorTag>>(*this);
        // Build the per-node halo topology immediately. Its constructor also
        // fixes duplicate ownership (lowest-rank-wins via Allgatherv) which
        // must happen BEFORE any caller reads d_nodeOwnership_ to size DOF
        // mappings or build sparsity.
        if (numRanks_ > 1)
        {
            nodeHaloTopo_ = std::make_unique<NodeHaloTopology<ElementTag, RealType, KeyType, AcceleratorTag>>(*this);
        }
    }
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::ensureCoordinateCache() const
{
    if (!coordCache_)
    {
        coordCache_ = std::make_unique<CoordinateCache<ElementTag, RealType, KeyType, AcceleratorTag>>(*this);
        std::cout << "Rank " << rank_ << " cached coordinates lazily." << std::endl;
    }
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::ensureNodeHaloTopo() const
{
    if (!nodeHaloTopo_)
    {
        // ensureHalo() builds both halo_ and nodeHaloTopo_ in one shot.
        ensureHalo();
    }
}

// ============================================================================
// AdjacencyData implementation
// ============================================================================
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
AdjacencyData<ElementTag, RealType, KeyType, AcceleratorTag>::AdjacencyData(
    const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain)
{
    createElementToNodeLocalIdMap(domain);
    buildNodeToElementMap(domain);
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void AdjacencyData<ElementTag, RealType, KeyType, AcceleratorTag>::createElementToNodeLocalIdMap(
    const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain)
{
    if constexpr (std::is_same_v<AcceleratorTag, cstone::GpuTag>)
    {
        size_t elementCount = domain.getElementCount();

        // Resize local ID connectivity to match element count
        if constexpr (NodesPerElement >= 1) std::get<0>(d_conn_local_ids_).resize(elementCount);
        if constexpr (NodesPerElement >= 2) std::get<1>(d_conn_local_ids_).resize(elementCount);
        if constexpr (NodesPerElement >= 3) std::get<2>(d_conn_local_ids_).resize(elementCount);
        if constexpr (NodesPerElement >= 4) std::get<3>(d_conn_local_ids_).resize(elementCount);
        if constexpr (NodesPerElement >= 5) std::get<4>(d_conn_local_ids_).resize(elementCount);
        if constexpr (NodesPerElement >= 6) std::get<5>(d_conn_local_ids_).resize(elementCount);
        if constexpr (NodesPerElement >= 7) std::get<6>(d_conn_local_ids_).resize(elementCount);
        if constexpr (NodesPerElement >= 8) std::get<7>(d_conn_local_ids_).resize(elementCount);

        int blockSize = 256;
        int numBlocks = (elementCount + blockSize - 1) / blockSize;

        // Map SFC keys to local IDs for each connectivity array
        auto mapKernel = [&](auto I)
        {
            mapSfcToLocalIdKernel<<<numBlocks, blockSize>>>(
                thrust::raw_pointer_cast(std::get<I>(domain.d_conn_keys_).data()),
                thrust::raw_pointer_cast(std::get<I>(d_conn_local_ids_).data()),
                thrust::raw_pointer_cast(domain.d_localToGlobalSfcMap_.data()), elementCount, domain.getNodeCount());
            cudaCheckError();
        };

        // Call for each connectivity index
        [&]<size_t... Is>(std::index_sequence<Is...>)
        { (mapKernel(std::integral_constant<size_t, Is>{}), ...); }(std::make_index_sequence<NodesPerElement>{});
    }
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void AdjacencyData<ElementTag, RealType, KeyType, AcceleratorTag>::buildNodeToElementMap(
    const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain)
{
    if constexpr (std::is_same_v<AcceleratorTag, cstone::GpuTag>)
    {
        size_t elementCount           = domain.getElementCount();
        size_t nodeCount              = domain.getNodeCount();
        size_t numConnectivityEntries = elementCount * NodesPerElement;

        if (numConnectivityEntries == 0)
        {
            d_nodeToElementOffsets_.resize(nodeCount + 1);
            thrust::fill(
                thrust::device, thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(d_nodeToElementOffsets_.data())),
                thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(d_nodeToElementOffsets_.data()) + nodeCount + 1),
                KeyType(0));
            d_nodeToElementList_.resize(0);
            return;
        }

        // Step 1: Create flat list of (local_node_id, element_id) pairs
        DeviceVector<KeyType> flatNodeIds(numConnectivityEntries);
        DeviceVector<KeyType> flatElemIds(numConnectivityEntries);

        // Copy connectivity to flat arrays
        auto copyFlat = [&]<size_t I>(std::integral_constant<size_t, I>)
        {
            size_t offset = I * elementCount;
            thrust::copy(thrust::device,
                         thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(std::get<I>(d_conn_local_ids_).data())),
                         thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(std::get<I>(d_conn_local_ids_).data()) +
                                                     elementCount),
                         thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatNodeIds.data()) + offset));
            thrust::sequence(
                thrust::device, thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatElemIds.data()) + offset),
                thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatElemIds.data()) + offset + elementCount));
        };

        [&]<size_t... Is>(std::index_sequence<Is...>)
        { (copyFlat(std::integral_constant<size_t, Is>{}), ...); }(std::make_index_sequence<NodesPerElement>{});

        // Step 2: Sort by node ID to group all elements per node
        thrust::sort_by_key(
            thrust::device, thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatNodeIds.data())),
            thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatNodeIds.data()) + flatNodeIds.size()),
            thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatElemIds.data())));

        // Step 3: Build CSR structure using reduce_by_key
        DeviceVector<KeyType> uniqueNodeIds(nodeCount);
        DeviceVector<KeyType> countsPerNode(nodeCount);

        auto new_end = thrust::reduce_by_key(
            thrust::device, thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatNodeIds.data())),
            thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatNodeIds.data()) + flatNodeIds.size()),
            thrust::constant_iterator<KeyType>(1),
            thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(uniqueNodeIds.data())),
            thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(countsPerNode.data())));

        uniqueNodeIds.resize(new_end.first -
                             thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(uniqueNodeIds.data())));
        countsPerNode.resize(new_end.second -
                             thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(countsPerNode.data())));

        // Step 4: Build offset array (prefix sum of counts)
        d_nodeToElementOffsets_.resize(nodeCount + 1);
        thrust::fill(thrust::device,
                     thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(d_nodeToElementOffsets_.data())),
                     thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(d_nodeToElementOffsets_.data()) +
                                                 d_nodeToElementOffsets_.size()),
                     KeyType(0));

        thrust::scatter(
            thrust::device, thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(countsPerNode.data())),
            thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(countsPerNode.data()) + countsPerNode.size()),
            thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(uniqueNodeIds.data())),
            thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(d_nodeToElementOffsets_.data()) + 1));

        thrust::inclusive_scan(thrust::device,
                               thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(d_nodeToElementOffsets_.data())),
                               thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(d_nodeToElementOffsets_.data()) +
                                                           d_nodeToElementOffsets_.size()),
                               thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(d_nodeToElementOffsets_.data())));

        // Step 5: Copy the element list
        d_nodeToElementList_ = flatElemIds;
    }
}

// ============================================================================
// HaloData implementation
// ============================================================================
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
HaloData<ElementTag, RealType, KeyType, AcceleratorTag>::HaloData(
    const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain)
{
    // Build halo elements first, then ownership (so Phase 3 can detect shared nodes)
    buildHaloElementIndices(domain);
    buildNodeOwnership(domain);
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void HaloData<ElementTag, RealType, KeyType, AcceleratorTag>::buildHaloElementIndices(
    const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain)
{
    if constexpr (std::is_same_v<AcceleratorTag, cstone::GpuTag>)
    {
        auto haloRanges = domain.haloElementRanges();

        // Calculate total halo count
        size_t totalHalos = 0;
        for (auto [start, end] : haloRanges)
        {
            totalHalos += (end - start);
        }

        if (totalHalos == 0)
        {
            d_haloElementIndices_.resize(0);
            return;
        }

        d_haloElementIndices_.resize(totalHalos);

        // Flatten the ranges into a single contiguous vector
        size_t offset = 0;
        for (auto [start, end] : haloRanges)
        {
            size_t rangeSize = end - start;
            thrust::sequence(
                thrust::device,
                thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(d_haloElementIndices_.data()) + offset),
                thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(d_haloElementIndices_.data()) + offset +
                                            rangeSize),
                start);
            offset += rangeSize;
        }
    }
}

// ============================================================================
// CoordinateCache implementation
// ============================================================================
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
CoordinateCache<ElementTag, RealType, KeyType, AcceleratorTag>::CoordinateCache(
    const ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>& domain)
{
    if constexpr (std::is_same_v<AcceleratorTag, cstone::GpuTag>)
    {
        size_t nodeCount = domain.getNodeCount();
        d_node_x_.resize(nodeCount);
        d_node_y_.resize(nodeCount);
        d_node_z_.resize(nodeCount);

        // Decode all SFC keys to coordinates (parallel kernel)
        int blockSize = 256;
        int numBlocks = (nodeCount + blockSize - 1) / blockSize;
        decodeAllNodesKernel<<<numBlocks, blockSize>>>(
            thrust::raw_pointer_cast(domain.getLocalToGlobalSfcMap().data()),
            thrust::raw_pointer_cast(d_node_x_.data()), thrust::raw_pointer_cast(d_node_y_.data()),
            thrust::raw_pointer_cast(d_node_z_.data()), nodeCount, domain.getBoundingBox());
        cudaCheckError();

        std::cout << "Rank " << domain.rank() << " cached " << nodeCount << " node coordinates." << std::endl;
    }
}

// ============================================================================
// OriginalCoordinates implementation
// ============================================================================
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
OriginalCoordinates<ElementTag, RealType, KeyType, AcceleratorTag>::OriginalCoordinates(size_t nodeCount)
{
    // Allocate storage for exact original coordinates
    // Actual values will be filled during sync() after node reordering
    d_node_x_.resize(nodeCount);
    d_node_y_.resize(nodeCount);
    d_node_z_.resize(nodeCount);
}

// Helper to populate the ConnPtrs struct from tuple
template<typename KeyType, size_t... Is, typename Tuple>
auto makeConnPtrs(const Tuple& tuple, std::index_sequence<Is...>)
{
    ConnPtrs<KeyType, sizeof...(Is)> result;
    ((result.ptrs[Is] = thrust::raw_pointer_cast(std::get<Is>(tuple).data())), ...);
    return result;
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::createLocalToGlobalSfcMap()
{
    if constexpr (std::is_same_v<AcceleratorTag, cstone::GpuTag>)
    {
        // 1. Gather all node SFC keys from the element connectivity lists into a single flat vector
        size_t numConnectivityEntries = elementCount_ * NodesPerElement;
        if (numConnectivityEntries == 0)
        {
            nodeCount_ = 0;
            return;
        }

        DeviceVector<KeyType> allNodeKeys(numConnectivityEntries);

        auto conn_ptrs = makeConnPtrs<KeyType>(d_conn_keys_, std::make_index_sequence<NodesPerElement>{});
        // 2. Launch the optimized kernel to flatten the tuple-of-vectors
        int blockSize = 256;
        int numBlocks = (elementCount_ + blockSize - 1) / blockSize;
        flattenConnectivityKernel<<<numBlocks, blockSize>>>(conn_ptrs, thrust::raw_pointer_cast(allNodeKeys.data()),
                                                            elementCount_);
        cudaCheckError();

        // 3. Sort the keys using Thrust's highly optimized parallel sort
        auto start = thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(allNodeKeys.data()));
        auto end   = start + allNodeKeys.size();
        thrust::sort(thrust::device, start, end);

        // 4. Find the unique keys - moves duplicates to the end
        auto newEnd = thrust::unique(thrust::device, start, end);

        // 5. Resize to discard duplicates
        allNodeKeys.resize(newEnd - start);

        // 6. Move the result to the member variable
        d_localToGlobalSfcMap_ = std::move(allNodeKeys);

        // 7. Update the final node count
        nodeCount_ = d_localToGlobalSfcMap_.size();
    }
}

// ============================================================================
// resyncFromDevice: AMR rebuild path that reuses the cstone Domain object
// ============================================================================
// Diagnosis (cube256/16, AMR Level 1, 977M elements):
//   focusTree.converge:  2942 ms  (firstCall_-only; 80% of cstone sync time)
//   focusTree.updateTree: 421 ms
//   setupHalos:            78 ms
//   ... others, summing to ~650 ms outside converge
//
// Saving converge alone is ~2.5 s per AMR rebuild at this scale. At Gordon-Bell
// scale (10^5+ ranks), firstCall_=true also fires the global-octree Allgather
// inside distribute(), which scales as O(numRanks * bucketSize) host-side work
// and becomes the actual scaling wall. Reusing the Domain bypasses both.
//
// Mechanism:
//   - cstone::Domain::sync() runs `if (firstCall_) focusTree_.converge();`
//   - On second and subsequent calls, firstCall_=false, converge is skipped
//   - The do-while retry loop in cstone sync() handles element-count change
//     via updateTree's incremental rebalance (cstone's design intent: "particles
//     moved between syncs" — AMR is a more aggressive case of that)
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::resyncFromDevice(
    DeviceCoordsTuple&& d_coords_new,
    DeviceConnectivityTuple&& d_conn_new)
{
    if constexpr (!std::is_same_v<AcceleratorTag, cstone::GpuTag>)
    {
        throw std::runtime_error("resyncFromDevice requires GpuTag");
    }

    if (!domain_)
    {
        throw std::runtime_error("resyncFromDevice: cstone Domain has not been constructed yet");
    }

    auto& d_x = std::get<0>(d_coords_new);
    auto& d_y = std::get<1>(d_coords_new);
    auto& d_z = std::get<2>(d_coords_new);

    nodeCount_    = d_x.size();
    elementCount_ = std::get<0>(d_conn_new).size();

    // Recompute bounding box. AMR refinement preserves geometric extents in
    // theory, but recomputing is cheap (< 1 ms) and protects against any
    // floating-point drift in refinement coords.
    RealType lxmin, lxmax, lymin, lymax, lzmin, lzmax;
    {
        auto xb = thrust::device_pointer_cast(d_x.data());
        auto yb = thrust::device_pointer_cast(d_y.data());
        auto zb = thrust::device_pointer_cast(d_z.data());
        lxmin = thrust::reduce(xb, xb + nodeCount_, std::numeric_limits<RealType>::max(),    thrust::minimum<RealType>());
        lxmax = thrust::reduce(xb, xb + nodeCount_, std::numeric_limits<RealType>::lowest(), thrust::maximum<RealType>());
        lymin = thrust::reduce(yb, yb + nodeCount_, std::numeric_limits<RealType>::max(),    thrust::minimum<RealType>());
        lymax = thrust::reduce(yb, yb + nodeCount_, std::numeric_limits<RealType>::lowest(), thrust::maximum<RealType>());
        lzmin = thrust::reduce(zb, zb + nodeCount_, std::numeric_limits<RealType>::max(),    thrust::minimum<RealType>());
        lzmax = thrust::reduce(zb, zb + nodeCount_, std::numeric_limits<RealType>::lowest(), thrust::maximum<RealType>());
    }

    auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    RealType gxmin, gxmax, gymin, gymax, gzmin, gzmax;
    MPI_Allreduce(&lxmin, &gxmin, 1, mpiType, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&lxmax, &gxmax, 1, mpiType, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&lymin, &gymin, 1, mpiType, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&lymax, &gymax, 1, mpiType, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&lzmin, &gzmin, 1, mpiType, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&lzmax, &gzmax, 1, mpiType, MPI_MAX, MPI_COMM_WORLD);

    RealType pad = RealType(0.05);
    RealType rx  = gxmax - gxmin, ry = gymax - gymin, rz = gzmax - gzmin;
    box_ = cstone::Box<RealType>(gxmin - pad * rx, gxmax + pad * rx,
                                  gymin - pad * ry, gymax + pad * ry,
                                  gzmin - pad * rz, gzmax + pad * rz);
    // NOTE: cstone Domain captured the OLD box at construction; cstone's
    // Box is mostly used by sync() at call time (passed via box() accessor),
    // so the new bbox is picked up on the next sync(). If cstone caches box
    // internally between syncs, we may need a setter. Verified by inspection
    // that cstone Domain stores box in global_ and reads it via global_.box()
    // at sync time, so this is safe.

    // Reset all lazy-cached state. They'll re-init on next access against the
    // new mesh. Storage of the prior versions is freed here (RAII via unique_ptr).
    adjacency_.reset();
    halo_.reset();
    nodeHaloTopo_.reset();
    coordCache_.reset();
    originalCoords_.reset();

    // Clear per-instance device vectors that depend on element/node counts.
    // The actual storage will be re-allocated by sync() / lazy init.
    // (cstone::DeviceVector exposes resize() but not clear().)
    d_elemSfcCodes_.resize(0);
    d_elemToNodeMap_.resize(0);
    d_localToGlobalSfcMap_.resize(0);
    d_localToGlobalNodeMap_.resize(0);
    // d_conn_keys_ and d_props_ are rebuilt during sync(); explicit clear()
    // not strictly needed, but reduces peak memory footprint during the
    // transition.

    // Allocate per-element h-property storage at the new size. (sync() expects
    // d_props_ to exist; calculateCharacteristicSizes populates it.)
    std::get<0>(d_props_).resize(nodeCount_);

    // Move the new mesh into local variables (sync() reads them by const ref).
    DeviceConnectivityTuple d_conn_local   = std::move(d_conn_new);
    DeviceCoordsTuple       d_coords_local = std::move(d_coords_new);
    calculateCharacteristicSizes(d_conn_local, d_coords_local);

    // CRITICAL: cstone Domain's bufDesc_/layout_/focusTree state was
    // initialized on firstCall_ to the OLD particle distribution and does not
    // auto-reset on count change. setEndIndex() alone is insufficient (only
    // updates bufDesc_.end; .start/.size and focus-tree rebalanceStatus_ stay
    // stale → segfault on next updateTree precondition).
    //
    // The cstone patch (scripts/patch_cstone_amr_reset.py) adds
    // resetForAMRRefinement(newCount) which resets every persistent field that
    // makes "AMR-style sync after refinement" work:
    //   - bufDesc_.{start,end,size} → {0, newCount, newCount}
    //   - prevBufDesc_              → match new bufDesc_
    //   - layout_                   → {0, newCount}
    //   - focusTree_.{prevFocusStart,prevFocusEnd,rebalanceStatus_}
    //
    // After this, cstone's distribute() rebuilds global_ from scratch (it
    // does so regardless of firstCall_), updateMinMac re-arms rebalanceStatus_
    // back to `valid`, and updateTree runs incrementally on the new
    // distribution. firstCall_ stays `false` so the expensive converge() loop
    // is skipped — that's the win we're after.
    domain_->resetForAMRRefinement(static_cast<cstone::LocalIndex>(elementCount_));

    logGpuMemAroundSync(rank_, "pre-resync");
    // Calls cstone Domain::sync() under the hood. Since domain_ already had
    // a successful first sync, firstCall_=false on entry → converge skipped.
    sync(d_conn_local, d_coords_local);
    logGpuMemAroundSync(rank_, "post-resync");

    if (rank_ == 0)
    {
        std::cout << "Rank " << rank_ << " resynced " << ElementTag::Name
                  << " domain (Domain reuse, no firstCall_ converge) with "
                  << nodeCount_ << " nodes and " << elementCount_ << " elements." << std::endl;
    }
}

} // namespace mars
