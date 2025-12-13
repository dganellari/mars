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
#include <filesystem>
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

// Forward declaration - implementation in mars_face_topology.hpp
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
class FaceTopology;

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
    ElementDomain(const std::string& meshFile, int rank, int numRanks);

    // Constructor from mesh data (for MFEM or other formats) - automatically computes bounding box
    ElementDomain(const HostCoordsTuple& h_coords, const HostConnectivityTuple& h_conn, int rank, int numRanks);

    // Constructor from mesh data with boundary info - automatically computes bounding box
    ElementDomain(const HostCoordsTuple& h_coords,
                  const HostConnectivityTuple& h_conn,
                  const HostBoundaryTuple& h_boundary,
                  int rank,
                  int numRanks);

    // Constructor from mesh data with explicit bounding box (for backward compatibility)
    ElementDomain(const HostCoordsTuple& h_coords,
                  const HostConnectivityTuple& h_conn,
                  const cstone::Box<RealType>& box,
                  int rank,
                  int numRanks);

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
    void cacheNodeCoordinates() { ensureCoordinateCache(); }

    // Exchange halo data for DOF vectors (MPI communication to sum shared node contributions)
    template<class... Vectors, class SendBuffer, class ReceiveBuffer>
    void exchangeHalos(std::tuple<Vectors&...> arrays, SendBuffer& sendBuffer, ReceiveBuffer& receiveBuffer) const
    {
        domain_->exchangeHalos(arrays, sendBuffer, receiveBuffer);
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
        ensureCoordinateCache();
        return coordCache_->d_node_x_;
    }

    const DeviceVector<RealType>& getNodeY() const
    {
        ensureCoordinateCache();
        return coordCache_->d_node_y_;
    }

    const DeviceVector<RealType>& getNodeZ() const
    {
        ensureCoordinateCache();
        return coordCache_->d_node_z_;
    }

    //! Returns boundary node flags (true if node is on boundary).
    const DeviceVector<uint8_t>& getBoundaryNodes() const { return std::get<0>(d_boundary_); }

    //! Returns true if boundary information is available.
    bool hasBoundaryInfo() const { return std::get<0>(d_boundary_).size() > 0; }

    //! Face topology access (lazy initialization for CVFEM)
    void buildFaceTopology() { ensureFaceTopology(); }
    
    KeyType getNumFaces() const { 
        ensureFaceTopology(); 
        return faceTopology_->numFaces_; 
    }
    
    KeyType getNumBoundaryFaces() const { 
        ensureFaceTopology(); 
        return faceTopology_->numBoundaryFaces_; 
    }
    
    KeyType getNumInteriorFaces() const { 
        ensureFaceTopology(); 
        return faceTopology_->numInteriorFaces_; 
    }

    const DeviceVector<KeyType>& getElementToFaces() const {
        ensureFaceTopology();
        return faceTopology_->d_elementToFaces_;
    }

    const DeviceVector<KeyType>& getFaceToElementOffsets() const {
        ensureFaceTopology();
        return faceTopology_->d_faceToElementOffsets_;
    }

    const DeviceVector<KeyType>& getFaceToElementList() const {
        ensureFaceTopology();
        return faceTopology_->d_faceToElementList_;
    }

    const DeviceVector<KeyType>& getFaceNodes() const {
        ensureFaceTopology();
        return faceTopology_->d_faceNodes_;
    }

    const DeviceVector<RealType>& getFaceNormals(int component) const {
        ensureFaceTopology();
        if (component == 0) return faceTopology_->d_faceNormalX_;
        if (component == 1) return faceTopology_->d_faceNormalY_;
        return faceTopology_->d_faceNormalZ_;
    }

    const DeviceVector<RealType>& getFaceAreas() const {
        ensureFaceTopology();
        return faceTopology_->d_faceArea_;
    }

    const DeviceVector<RealType>& getFaceCentroids(int component) const {
        ensureFaceTopology();
        if (component == 0) return faceTopology_->d_faceCentroidX_;
        if (component == 1) return faceTopology_->d_faceCentroidY_;
        return faceTopology_->d_faceCentroidZ_;
    }

    const DeviceVector<uint8_t>& getIsBoundaryFace() const {
        ensureFaceTopology();
        return faceTopology_->d_isBoundaryFace_;
    }

private:
    int rank_;
    int numRanks_;
    std::unique_ptr<DomainType> domain_;
    std::size_t elementCount_ = 0;
    std::size_t nodeCount_    = 0;

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
    friend struct CoordinateCache<ElementTag, RealType, KeyType, AcceleratorTag>;
    friend struct FaceTopology<ElementTag, RealType, KeyType, AcceleratorTag>;

    mutable std::unique_ptr<AdjacencyData<ElementTag, RealType, KeyType, AcceleratorTag>> adjacency_;
    mutable std::unique_ptr<HaloData<ElementTag, RealType, KeyType, AcceleratorTag>> halo_;
    mutable std::unique_ptr<CoordinateCache<ElementTag, RealType, KeyType, AcceleratorTag>> coordCache_;
    mutable std::unique_ptr<FaceTopology<ElementTag, RealType, KeyType, AcceleratorTag>> faceTopology_;

    void initializeConnectivityKeys();

    // Helper method for creating SFC map
    void createLocalToGlobalSfcMap();

    // Lazy initialization methods (called by public getters)
    void ensureSfcMap() const;
    void ensureAdjacency() const;
    void ensureHalo() const;
    void ensureCoordinateCache() const;
    void ensureFaceTopology() const;

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

// Element domain constructor
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::ElementDomain(const std::string& meshFile,
                                                                            int rank,
                                                                            int numRanks)
    : rank_(rank)
    , numRanks_(numRanks)
    , box_(0, 1)
{
    // Host data in SoA format
    HostCoordsTuple h_coords;     // (x, y, z)
    HostConnectivityTuple h_conn; // (i0, i1, i2, ...) depends on element type

    // Device data (will be filled by transferDataToGPU)
    DeviceConnectivityTuple d_conn_;
    DeviceCoordsTuple d_coords_;

    // Read the mesh in SoA format
    readMeshDataSoA(meshFile, h_coords, h_conn);

    // Initialize cornerstone domain
    int bucketSize           = 64;
    unsigned bucketSizeFocus = 8;
    RealType theta           = 0.5;

    box_ = computeGlobalBoundingBoxFromCoords<RealType>(std::get<0>(h_coords), std::get<1>(h_coords),
                                                        std::get<2>(h_coords));

    std::cout << "Rank " << rank_ << ": Created bounding box: [" << box_.xmin() << "," << box_.xmax() << "] ["
              << box_.ymin() << "," << box_.ymax() << "] [" << box_.zmin() << "," << box_.zmax() << "]" << std::endl;

    // Test SFC precision for this domain (debug only)
    testSfcPrecision<KeyType, RealType>(box_, rank_);

    domain_ = std::make_unique<DomainType>(rank, numRanks, bucketSize, bucketSizeFocus, theta, box_);

    // Transfer data to GPU before sync
    transferDataToGPU(h_coords, h_conn, d_coords_, d_conn_);

    // Calculate characteristic sizes
    calculateCharacteristicSizes(d_conn_, d_coords_);

    std::cout << "Rank " << rank_ << ": Before sync - bucketSize=" << bucketSize
              << ", bucketSizeFocus=" << bucketSizeFocus << ", theta=" << theta << std::endl;

    // Perform sync
    sync(d_conn_, d_coords_);

    std::cout << "Rank " << rank_ << " initialized " << ElementTag::Name << " domain with " << nodeCount_
              << " nodes and " << elementCount_ << " elements." << std::endl;
}

// Constructor from mesh data (for MFEM or other formats) - automatically computes bounding box
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::ElementDomain(const HostCoordsTuple& h_coords,
                                                                            const HostConnectivityTuple& h_conn,
                                                                            int rank,
                                                                            int numRanks)
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

    // Initialize cornerstone domain
    int bucketSize           = 64;
    unsigned bucketSizeFocus = 8;
    RealType theta           = 0.5;

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

    // Perform sync
    sync(d_conn_, d_coords_);

    std::cout << "Rank " << rank_ << " initialized " << ElementTag::Name << " domain with " << nodeCount_
              << " nodes and " << elementCount_ << " elements." << std::endl;
}

// Constructor from mesh data with boundary info - automatically computes bounding box
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::ElementDomain(const HostCoordsTuple& h_coords,
                                                                            const HostConnectivityTuple& h_conn,
                                                                            const HostBoundaryTuple& h_boundary,
                                                                            int rank,
                                                                            int numRanks)
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

    // Initialize cornerstone domain
    int bucketSize           = 64;
    unsigned bucketSizeFocus = 8;
    RealType theta           = 0.5;

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

    // Perform sync
    sync(d_conn_, d_coords_);

    std::cout << "Rank " << rank_ << " initialized " << ElementTag::Name << " domain with " << nodeCount_
              << " nodes and " << elementCount_ << " elements." << std::endl;
}

// Read mesh data in SoA format; uses element-based partitioning for better data locality
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::readMeshDataSoA(const std::string& meshFile,
                                                                                   HostCoordsTuple& h_coords_,
                                                                                   HostConnectivityTuple& h_conn_)
{
    std::cout << "Reading mesh data from " << meshFile << " on rank " << rank_ << std::endl;
    try
    {
        // Helper to handle coordinate conversion based on type
        auto processCoordinates =
            [&](const std::vector<float>& x_data, const std::vector<float>& y_data, const std::vector<float>& z_data)
        {
            auto& h_x = std::get<0>(h_coords_);
            auto& h_y = std::get<1>(h_coords_);
            auto& h_z = std::get<2>(h_coords_);

            if constexpr (std::is_same_v<RealType, float>)
            {
                // If RealType is float, we can move directly
                h_x = std::move(x_data);
                h_y = std::move(y_data);
                h_z = std::move(z_data);
            }
            else
            {
                // Otherwise convert
                h_x.resize(nodeCount_);
                h_y.resize(nodeCount_);
                h_z.resize(nodeCount_);

                for (size_t i = 0; i < nodeCount_; ++i)
                {
                    h_x[i] = static_cast<RealType>(x_data[i]);
                    h_y[i] = static_cast<RealType>(y_data[i]);
                    h_z[i] = static_cast<RealType>(z_data[i]);
                }
            }
        };

        // Detect if meshFile is an MFEM file or binary mesh directory
        // MFEM meshes are files with .mesh extension
        // Binary meshes are directories containing coordinate files
        namespace fs    = std::filesystem;
        bool isMFEMFile = !fs::is_directory(meshFile) && meshFile.find(".mesh") != std::string::npos;

        // Use compile-time branching to select element type
        if constexpr (std::is_same_v<ElementTag, TetTag>)
        {
            if (isMFEMFile)
            {
                // Use MFEM mesh reader
                auto [readNodeCount, readElementCount, x_data, y_data, z_data, conn_tuple, localToGlobal,
                      boundaryNodes] =
                    mars::readMFEMMeshWithElementPartitioning<4, float, KeyType>(meshFile, rank_, numRanks_);

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
                    mars::readMeshWithElementPartitioning<4, float, KeyType>(meshFile, rank_, numRanks_);

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
            if (isMFEMFile)
            {
                // Use MFEM mesh reader for hex elements
                auto [readNodeCount, readElementCount, x_data, y_data, z_data, conn_tuple, localToGlobal,
                      boundaryNodes] =
                    mars::readMFEMMeshWithElementPartitioning<8, float, KeyType>(meshFile, rank_, numRanks_);

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
                    mars::readMeshWithElementPartitioning<8, float, KeyType>(meshFile, rank_, numRanks_);

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
                auto [readNodeCount, readElementCount, x_data, y_data, z_data, conn_tuple, localToGlobal,
                      boundaryNodes] =
                    mars::readMFEMMeshWithElementPartitioning<3, float, KeyType>(meshFile, rank_, numRanks_);

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
                    mars::readMeshWithElementPartitioning<3, float, KeyType>(meshFile, rank_, numRanks_);

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
                auto [readNodeCount, readElementCount, x_data, y_data, z_data, conn_tuple, localToGlobal,
                      boundaryNodes] =
                    mars::readMFEMMeshWithElementPartitioning<4, float, KeyType>(meshFile, rank_, numRanks_);

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
                    mars::readMeshWithElementPartitioning<4, float, KeyType>(meshFile, rank_, numRanks_);

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

    // For tetrahedra, compute based on average edge lengths
    if constexpr (std::is_same_v<ElementTag, TetTag> && std::is_same_v<AcceleratorTag, cstone::GpuTag>)
    {
        // Fix: Correct the device vector declaration
        DeviceVector<int> d_nodeTetCount(nodeCount_, 0);

        // Extract raw pointers for kernel
        auto& d_x  = std::get<0>(d_coords_);
        auto& d_y  = std::get<1>(d_coords_);
        auto& d_z  = std::get<2>(d_coords_);
        auto& d_i0 = std::get<0>(d_conn_);
        auto& d_i1 = std::get<1>(d_conn_);
        auto& d_i2 = std::get<2>(d_conn_);
        auto& d_i3 = std::get<3>(d_conn_);

        constexpr int NodesPerElem = ElementTag::NodesPerElement;

        // Extract additional connectivity pointers conditionally at compile time
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

        // First accumulate edge lengths per node
        int blockSize = 256;
        int numBlocks = (elementCount_ + blockSize - 1) / blockSize;

        computeCharacteristicSizesKernel<ElementTag, KeyType, RealType><<<numBlocks, blockSize>>>(
            thrust::raw_pointer_cast(d_x.data()), thrust::raw_pointer_cast(d_y.data()),
            thrust::raw_pointer_cast(d_z.data()), thrust::raw_pointer_cast(d_i0.data()),
            thrust::raw_pointer_cast(d_i1.data()), thrust::raw_pointer_cast(d_i2.data()),
            thrust::raw_pointer_cast(d_i3.data()), d_i4_ptr, d_i5_ptr, d_i6_ptr, d_i7_ptr,
            thrust::raw_pointer_cast(d_nodeTetCount.data()), thrust::raw_pointer_cast(d_h.data()), elementCount_);

        cudaCheckError();

        // Then normalize by number of contributions
        numBlocks = (nodeCount_ + blockSize - 1) / blockSize;
        finalizeCharacteristicSizesKernel<KeyType, RealType><<<numBlocks, blockSize>>>(
            thrust::raw_pointer_cast(d_h.data()), thrust::raw_pointer_cast(d_nodeTetCount.data()), nodeCount_);

        cudaCheckError();

        // For mesh-based methods (FEM/FDM)
        constexpr RealType meshFactor = 1.0;    // No reduction for FEM (adjust based on element order)
        constexpr RealType minH       = 1.0e-6; // Prevent extremely small values that cause instability
        constexpr RealType maxH       = 1.0;    // Upper bound based on problem domain

        // Use the proper kernel
        transformCharacteristicSizesKernel<RealType>
            <<<numBlocks, blockSize>>>(thrust::raw_pointer_cast(d_h.data()), nodeCount_, meshFactor, minH, maxH);
        cudaCheckError();
    }
    else
    {
        // For non-GPU or non-tetrahedral cases, use a default value
        RealType domainDiagonal = std::sqrt(std::pow(domain_->box().xmax() - domain_->box().xmin(), 2) +
                                            std::pow(domain_->box().ymax() - domain_->box().ymin(), 2) +
                                            std::pow(domain_->box().zmax() - domain_->box().zmin(), 2));
        RealType defaultH       = domainDiagonal * 0.01; // 1% of domain diagonal

        // Use the proper kernel
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

        std::cout << "Rank " << rank_ << " syncing" << ElementTag::Name << " domain with " << elementCount_
                  << " elements." << std::endl;
        syncDomainImpl(domain_.get(), d_elemSfcCodes_, d_elemX, d_elemY, d_elemZ, d_elemH, elementCount_, d_conn_keys_);
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
void ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::ensureFaceTopology() const
{
    if (!faceTopology_)
    {
        // Face topology needs coordinate cache and adjacency
        ensureCoordinateCache();
        ensureAdjacency();
        faceTopology_ = std::make_unique<FaceTopology<ElementTag, RealType, KeyType, AcceleratorTag>>(*this);
        std::cout << "Rank " << rank_ << " built face topology: " << faceTopology_->numFaces_ 
                  << " faces (" << faceTopology_->numBoundaryFaces_ << " boundary, " 
                  << faceTopology_->numInteriorFaces_ << " interior)" << std::endl;
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

} // namespace mars

// Include FaceTopology implementation (GPU version)
#include "mars_face_topology_gpu.hpp"

