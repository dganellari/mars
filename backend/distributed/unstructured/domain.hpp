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
#include <iostream>
#include <memory>
#include <vector>

// Mars includes
#include "mars.hpp"
#include "mars_globals.hpp"
#include "mars_domain_utils.hpp"
// #include "mars_read_mesh_adios2.hpp"
#include "mars_read_mesh_binary.hpp"

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
                                    const KeyType* nodeSfcCodes,
                                    KeyType* conn_key0,
                                    KeyType* conn_key1,
                                    KeyType* conn_key2,
                                    KeyType* conn_key3,
                                    int numElements);

template<typename KeyType>
__global__ void decodeSfcToIntegersKernel(const KeyType* keys,
                                          unsigned* x,
                                          unsigned* y, 
                                          unsigned* z,
                                          size_t numKeys);

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
__global__ void convertSfcToNodeIndicesKernel(const KeyType* sfcIndices,
                                              KeyType* nodeIndices,
                                              const KeyType* particleKeys,
                                              size_t numElements,
                                              size_t numNodes);
                                    
// Better: use std::array for compile-time indexing
template<typename KeyType, size_t NodesPerElement>
struct ConnPtrs {
    const KeyType* ptrs[NodesPerElement];
};

template<typename KeyType, size_t NodesPerElement>
__global__ void flattenConnectivityKernel(
    ConnPtrs<KeyType, NodesPerElement> conn,
    KeyType* flat_keys, 
    size_t numElements);


template<typename KeyType>
__global__ void mapSfcToLocalIdKernel(const KeyType* sfc_conn, 
                                      KeyType* local_conn, 
                                      const KeyType* sorted_sfc, 
                                      size_t num_elements, 
                                      size_t num_nodes);
 
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
    using type = std::tuple<HostVector, 
                            HostVector,
                            HostVector,
                            HostVector>;
};

// Specialization for hexahedra (host)
template<typename T>
struct HostConnectivityTupleHelper<HexTag, T>
{
    using HostVector = typename VectorSelector<T, cstone::CpuTag>::type;
    using type = std::tuple<HostVector, 
                            HostVector,
                            HostVector,
                            HostVector,
                            HostVector,
                            HostVector,
                            HostVector,
                            HostVector>;
};

// Specialization for triangles (host)
template<typename T>
struct HostConnectivityTupleHelper<TriTag, T>
{
    using HostVector = typename VectorSelector<T, cstone::CpuTag>::type;
    using type = std::tuple<HostVector, HostVector, HostVector>;
};

// Specialization for quads (host)
template<typename T>
struct HostConnectivityTupleHelper<QuadTag, T>
{
    using HostVector = typename VectorSelector<T, cstone::CpuTag>::type; 
    using type = std::tuple<HostVector, HostVector, HostVector, HostVector>;
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

    // SoA data structures using tuples - host versions
    using HostCoordsTuple       = std::tuple<HostVector<RealType>, HostVector<RealType>, HostVector<RealType>>;
    using HostPropsTuple        = std::tuple<HostVector<RealType>>;
    using HostConnectivityTuple = typename HostConnectivityTupleHelper<ElementTag, KeyType>::type;

    ElementDomain(const std::string& meshFile, int rank, int numRanks);

    // GPU-accelerated calculation of characteristic sizes
    void calculateCharacteristicSizes(const DeviceConnectivityTuple& d_conn_,
                                      const DeviceCoordsTuple& d_coords_);

    // Transfer data to GPU for computations
    void transferDataToGPU(const HostCoordsTuple& h_coords_, const HostConnectivityTuple& h_conn_, DeviceCoordsTuple& d_coords_,
                           DeviceConnectivityTuple& d_conn_);

    // Domain synchronization (following cornerstone API pattern)
    void sync(const DeviceConnectivityTuple& d_conn_,
              const DeviceCoordsTuple& d_coords_);

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
    std::size_t haloElementCount() const
    {
        return getElementCount() - localElementCount();
    }

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

    //build adjacency after sync
    void buildAdjacency();

    // Device connectivity functions
    template<int I>
    __device__ KeyType getConnectivityDevice(size_t elementIndex) const;

    int numRanks() const { return numRanks_; }
    int rank() const { return rank_; }

    //! Returns the map from a dense local node ID to its sparse global SFC Key.
    const DeviceVector<KeyType>& getLocalToGlobalSfcMap() const { return d_localToGlobalSfcMap_; }

    //! Returns the Element->Node connectivity table using dense local node IDs.
    const DeviceConnectivityTuple& getElementToNodeConnectivity() const { return d_conn_local_ids_; }

    //! Returns the offsets for the Node->Element CSR map.
    const DeviceVector<KeyType>& getNodeToElementOffsets() const { return d_nodeToElementOffsets_; }

    //! Returns the list for the Node->Element CSR map.
    const DeviceVector<KeyType>& getNodeToElementList() const { return d_nodeToElementList_; }
   
private:
    int rank_;
    int numRanks_;
    std::unique_ptr<DomainType> domain_;
    std::size_t elementCount_ = 0;
    std::size_t nodeCount_    = 0;

    // element unique identifiers (SFC codes)
    DeviceVector<KeyType> d_elemSfcCodes_;

    // Device data in SoA format
    DevicePropsTuple d_props_;       // (h)
    DeviceConnectivityTuple d_conn_keys_; // (sfc_i0, sfc_i1, sfc_i2, ...) depends on element type; store sfc instead of coordinates
    //TODO: check if domain_->box() is a device function in cstone to avoid storing box_ here
    cstone::Box<RealType> box_; // bounding box for the domain

    //useful mapping from element to node indices
    DeviceVector<KeyType> d_elemToNodeMap_;

    // Adjacency and mapping data structures
    // Maps a dense local node ID [0...N-1] to its sparse global SFC Key
    DeviceVector<KeyType> d_localToGlobalSfcMap_;
    // Element->Node connectivity using dense local IDs [0...N-1]
    DeviceConnectivityTuple d_conn_local_ids_;
    // Node-to-Element connectivity in CSR format
    DeviceVector<KeyType> d_nodeToElementOffsets_;
    DeviceVector<KeyType> d_nodeToElementList_;

    void initializeConnectivityKeys();

    // Helper methods for building maps
    void createLocalToGlobalSfcMap();
    void createElementToNodeLocalIdMap();
    void buildNodeToElementMap();

    void syncImpl(DeviceVector<KeyType>& d_nodeSfcCodes, const DeviceConnectivityTuple& d_conn_, int blockSize, int numBlocks)
    {
        if constexpr (std::is_same_v<ElementTag, TetTag>)
        {
            auto& d_i0 = std::get<0>(d_conn_);
            auto& d_i1 = std::get<1>(d_conn_);
            auto& d_i2 = std::get<2>(d_conn_);
            auto& d_i3 = std::get<3>(d_conn_);

            // Find representative nodes
            findRepresentativeNodesKernel<TetTag, KeyType, RealType>
                <<<numBlocks, blockSize>>>(thrust::raw_pointer_cast(d_i0.data()), thrust::raw_pointer_cast(d_i1.data()),
                                           thrust::raw_pointer_cast(d_i2.data()), thrust::raw_pointer_cast(d_i3.data()),
                                           thrust::raw_pointer_cast(d_nodeSfcCodes.data()),
                                           thrust::raw_pointer_cast(d_elemToNodeMap_.data()), elementCount_);
            cudaCheckError();

            // Build SFC connectivity
            buildSfcConnectivity<TetTag, KeyType, RealType><<<numBlocks, blockSize>>>(
                thrust::raw_pointer_cast(d_i0.data()), thrust::raw_pointer_cast(d_i1.data()),
                thrust::raw_pointer_cast(d_i2.data()), thrust::raw_pointer_cast(d_i3.data()),
                thrust::raw_pointer_cast(d_nodeSfcCodes.data()), 
                thrust::raw_pointer_cast(std::get<0>(d_conn_keys_).data()),
                thrust::raw_pointer_cast(std::get<1>(d_conn_keys_).data()),
                thrust::raw_pointer_cast(std::get<2>(d_conn_keys_).data()),
                thrust::raw_pointer_cast(std::get<3>(d_conn_keys_).data()), elementCount_);
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

            // Find representative nodes (would need 8-node version)
            findRepresentativeNodesKernel<HexTag, KeyType, RealType>
                <<<numBlocks, blockSize>>>(thrust::raw_pointer_cast(d_i0.data()), thrust::raw_pointer_cast(d_i1.data()),
                                           thrust::raw_pointer_cast(d_i2.data()), thrust::raw_pointer_cast(d_i3.data()),
                                           thrust::raw_pointer_cast(d_i4.data()), thrust::raw_pointer_cast(d_i5.data()),
                                           thrust::raw_pointer_cast(d_i6.data()), thrust::raw_pointer_cast(d_i7.data()),
                                           thrust::raw_pointer_cast(d_nodeSfcCodes.data()),
                                           thrust::raw_pointer_cast(d_elemToNodeMap_.data()), elementCount_);
            cudaCheckError();

            // Build SFC connectivity for hex (would need 8-node version)
            // buildSfcConnectivity<HexTag, KeyType, RealType><<<...>>>(...);
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
                thrust::raw_pointer_cast(d_nodeSfcCodes.data()),
                thrust::raw_pointer_cast(d_elemToNodeMap_.data()), elementCount_);
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
            buildSfcConnectivity<QuadTag, KeyType, RealType><<<numBlocks, blockSize>>>(
                thrust::raw_pointer_cast(d_i0.data()), thrust::raw_pointer_cast(d_i1.data()),
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
    // Try float32 first (since that's what most meshes use), then double
    std::vector<std::string> extensions;
    if constexpr (std::is_same_v<RealType, float>) {
        extensions = {"float32", "double"};
    } else {
        extensions = {"double", "float32"};
    }
    
    std::string ext;
    std::ifstream x_file, y_file, z_file;
    
    // Try different file extensions
    for (const auto& test_ext : extensions) {
        x_file.open(meshDir + "/x." + test_ext, std::ios::binary);
        y_file.open(meshDir + "/y." + test_ext, std::ios::binary);
        z_file.open(meshDir + "/z." + test_ext, std::ios::binary);
        
        if (x_file && y_file && z_file) {
            ext = test_ext;
            break;
        }
        
        x_file.close(); y_file.close(); z_file.close();
    }
    
    if (!x_file || !y_file || !z_file) {
        throw std::runtime_error("Failed to open coordinate files for global bounding box computation: " + meshDir);
    }

    // Get total node count
    x_file.seekg(0, std::ios::end);
    size_t fileSize = x_file.tellg();
    x_file.seekg(0, std::ios::beg);
    
    size_t totalNodes;
    if (ext == "float32") {
        totalNodes = fileSize / sizeof(float);
    } else {
        totalNodes = fileSize / sizeof(double);
    }

    // Read ALL data at once (like createBoundingBox approach)
    if (ext == "float32") {
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
        RealType xRange = xmax - xmin;
        RealType yRange = ymax - ymin;
        RealType zRange = zmax - zmin;

        return cstone::Box<RealType>(
            xmin - padding * xRange, xmax + padding * xRange,
            ymin - padding * yRange, ymax + padding * yRange,
            zmin - padding * zRange, zmax + padding * zRange
        );
    } else {
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
        RealType xRange = xmax - xmin;
        RealType yRange = ymax - ymin;
        RealType zRange = zmax - zmin;

        return cstone::Box<RealType>(
            xmin - padding * xRange, xmax + padding * xRange,
            ymin - padding * yRange, ymax + padding * yRange,
            zmin - padding * zRange, zmax + padding * zRange
        );
    }
}

// Warning: This would create a local bounding box if called on a part on local coords.
template<typename RealType, typename CoordTuple>
cstone::Box<RealType> createBoundingBox(const CoordTuple& coords, RealType padding = 0.05)
{
    std::cout << "⚠️  WARNING: Creating LOCAL bounding box!" << std::endl;
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

template<typename KeyType, typename RealType>
void testSfcPrecision(const cstone::Box<RealType>& box, int rank = 0)
{
#ifndef NDEBUG  // Only run in debug builds (when NDEBUG is NOT defined)
    RealType domainSizeX = box.xmax() - box.xmin();
    RealType domainSizeY = box.ymax() - box.ymin(); 
    RealType domainSizeZ = box.zmax() - box.zmin();
    RealType maxDomainSize = std::max({domainSizeX, domainSizeY, domainSizeZ});
    
    constexpr int totalBits = sizeof(KeyType) * 8;
    constexpr int bitsPerDim = totalBits / 3;
    RealType sfcPrecision = maxDomainSize / (1ULL << bitsPerDim);
    
    if (rank == 0) {
        std::cout << "=== SFC Precision Analysis (DEBUG) ===" << std::endl;
        std::cout << "KeyType: " << (sizeof(KeyType) == 4 ? "uint32_t" : "uint64_t") 
                  << " (" << totalBits << " bits)" << std::endl;
        std::cout << "Domain size: [" << domainSizeX << " x " << domainSizeY 
                  << " x " << domainSizeZ << "]" << std::endl;
        std::cout << "SFC precision: " << sfcPrecision << " units" << std::endl;
        
        if (sizeof(KeyType) == 4 && sfcPrecision > 0.001) {
            std::cout << "⚠️  WARNING: Consider using uint64_t KeyType!" << std::endl;
        } else {
            std::cout << "✅ SFC precision should be sufficient." << std::endl;
        }
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

    DeviceConnectivityTuple d_conn_; // (i0, i1, i2, ...) depends on element type
    DeviceCoordsTuple d_coords_; // (x, y, z)

    // Read the mesh in SoA format
    readMeshDataSoA(meshFile, h_coords, h_conn);

    // Initialize cornerstone domain
    int bucketSize           = 128;
    unsigned bucketSizeFocus = 64;
    RealType theta           = 0.5;

    // Create a bounding box with some padding
    box_ = computeGlobalBoundingBox<RealType>(meshFile);

    std::cout << "Rank " << rank_ << ": Created bounding box: [" 
          << box_.xmin() << "," << box_.xmax() << "] [" 
          << box_.ymin() << "," << box_.ymax() << "] [" 
          << box_.zmin() << "," << box_.zmax() << "]" << std::endl;

    // Test SFC precision for this domain (debug only)
    testSfcPrecision<KeyType, RealType>(box_, rank_);

    domain_ = std::make_unique<DomainType>(rank, numRanks, bucketSize, bucketSizeFocus, theta, box_);

    // Transfer data to GPU before sync
    transferDataToGPU(h_coords, h_conn, d_coords_, d_conn_);

    // Calculate characteristic sizes
    calculateCharacteristicSizes(d_conn_, d_coords_);

    // Perform sync
    sync(d_conn_, d_coords_);

    // build adjacency maps now that the domain is settled
    buildAdjacency();

    std::cout << "Rank " << rank_ << " initialized " << ElementTag::Name << " domain with " << nodeCount_
              << " nodes and " << elementCount_ << " elements and " << localElementCount() << " local elements." << std::endl;
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

        // Use compile-time branching to select element type
        if constexpr (std::is_same_v<ElementTag, TetTag>)
        {
            auto [readNodeCount, readElementCount, x_data, y_data, z_data, conn_tuple, localToGlobal] =
                mars::readMeshWithElementPartitioning<4, float, KeyType>(meshFile, rank_, numRanks_);

            nodeCount_    = readNodeCount;
            elementCount_ = readElementCount;
            processCoordinates(x_data, y_data, z_data);
            h_conn_ = std::move(conn_tuple);
        }
        else if constexpr (std::is_same_v<ElementTag, HexTag>)
        {
            auto [readNodeCount, readElementCount, x_data, y_data, z_data, conn_tuple, localToGlobal] =
                mars::readMeshWithElementPartitioning<8, float, KeyType>(meshFile, rank_, numRanks_);

            nodeCount_    = readNodeCount;
            elementCount_ = readElementCount;
            processCoordinates(x_data, y_data, z_data);
            h_conn_ = std::move(conn_tuple);
        }
        else if constexpr (std::is_same_v<ElementTag, TriTag>)
        {
            auto [readNodeCount, readElementCount, x_data, y_data, z_data, conn_tuple, localToGlobal] =
                mars::readMeshWithElementPartitioning<3, float, KeyType>(meshFile, rank_, numRanks_);

            nodeCount_    = readNodeCount;
            elementCount_ = readElementCount;
            processCoordinates(x_data, y_data, z_data);
            h_conn_ = std::move(conn_tuple);
        }
        else if constexpr (std::is_same_v<ElementTag, QuadTag>)
        {
            auto [readNodeCount, readElementCount, x_data, y_data, z_data, conn_tuple, localToGlobal] =
                mars::readMeshWithElementPartitioning<4, float, KeyType>(meshFile, rank_, numRanks_);

            nodeCount_    = readNodeCount;
            elementCount_ = readElementCount;
            processCoordinates(x_data, y_data, z_data);
            h_conn_ = std::move(conn_tuple);
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
void ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::calculateCharacteristicSizes(const DeviceConnectivityTuple& d_conn_, const DeviceCoordsTuple& d_coords_)
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

        // First accumulate edge lengths per node
        int blockSize = 256;
        int numBlocks = (elementCount_ + blockSize - 1) / blockSize;

        computeCharacteristicSizesKernel<TetTag, KeyType, RealType><<<numBlocks, blockSize>>>(
            thrust::raw_pointer_cast(d_x.data()), thrust::raw_pointer_cast(d_y.data()),
            thrust::raw_pointer_cast(d_z.data()), thrust::raw_pointer_cast(d_i0.data()),
            thrust::raw_pointer_cast(d_i1.data()), thrust::raw_pointer_cast(d_i2.data()),
            thrust::raw_pointer_cast(d_i3.data()), thrust::raw_pointer_cast(d_nodeTetCount.data()),
            thrust::raw_pointer_cast(d_h.data()), elementCount_);

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
                                           thrust::raw_pointer_cast(d_nodeSfcCodes.data()), nodeCount_, getBoundingBox());

        // Find representative nodes for each element
        d_elemToNodeMap_.resize(elementCount_);

        //initialize connectivity keys storage
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

        std::cout << "Rank " << rank_ << " syncing" << ElementTag::Name << " domain with "
                  << elementCount_ << " elements." << std::endl;
        syncDomainImpl(domain_.get(), d_elemSfcCodes_, d_elemX, d_elemY, d_elemZ, d_elemH, elementCount_, d_conn_keys_);
    }
}

// Transfer data to GPU for computations
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::transferDataToGPU(
    const HostCoordsTuple& h_coords_, const HostConnectivityTuple& h_conn_, DeviceCoordsTuple& d_coords_,
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

// Implementation of SFC-to-coordinate functions
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
MARS_HOST_DEVICE std::tuple<RealType, RealType, RealType> 
ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::sfcToPhysicalCoordinate(KeyType sfcKey) const
{
    // Convert raw key to SfcKind strong type
    auto sfcKindKey = cstone::SfcKind<KeyType>(sfcKey);
    auto [ix, iy, iz] = cstone::decodeSfc(sfcKindKey);
    
    // Use SfcKind for maxTreeLevel
    constexpr unsigned maxCoord = (1u << cstone::maxTreeLevel<cstone::SfcKind<KeyType>>{}) - 1;
    RealType invMaxCoord = RealType(1.0) / maxCoord;
    auto box = getBoundingBox();
    
    RealType x = box.xmin() + ix * invMaxCoord * (box.xmax() - box.xmin());
    RealType y = box.ymin() + iy * invMaxCoord * (box.ymax() - box.ymin());
    RealType z = box.zmin() + iz * invMaxCoord * (box.zmax() - box.zmin());
    
    return std::make_tuple(x, y, z);
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
MARS_HOST_DEVICE RealType 
ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::sfcToPhysicalCoordinateX(KeyType sfcKey) const
{
    // Convert raw key to SfcKind strong type
    auto sfcKindKey = cstone::SfcKind<KeyType>(sfcKey);
    auto [ix, iy, iz] = cstone::decodeSfc(sfcKindKey);
    
    // Use SfcKind for maxTreeLevel
    constexpr unsigned maxCoord = (1u << cstone::maxTreeLevel<cstone::SfcKind<KeyType>>{}) - 1;
    RealType invMaxCoord = RealType(1.0) / maxCoord;
    auto box = getBoundingBox();
    
    return box.xmin() + ix * invMaxCoord * (box.xmax() - box.xmin());
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
MARS_HOST_DEVICE RealType 
ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::sfcToPhysicalCoordinateY(KeyType sfcKey) const
{
    // Convert raw key to SfcKind strong type
    auto sfcKindKey = cstone::SfcKind<KeyType>(sfcKey);
    auto [ix, iy, iz] = cstone::decodeSfc(sfcKindKey);
    
    // Use SfcKind for maxTreeLevel
    constexpr unsigned maxCoord = (1u << cstone::maxTreeLevel<cstone::SfcKind<KeyType>>{}) - 1;
    RealType invMaxCoord = RealType(1.0) / maxCoord;
    auto box = getBoundingBox();
    
    return box.ymin() + iy * invMaxCoord * (box.ymax() - box.ymin());
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
MARS_HOST_DEVICE RealType 
ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::sfcToPhysicalCoordinateZ(KeyType sfcKey) const
{
    // Convert raw key to SfcKind strong type
    auto sfcKindKey = cstone::SfcKind<KeyType>(sfcKey);
    auto [ix, iy, iz] = cstone::decodeSfc(sfcKindKey);
    
    // Use SfcKind for maxTreeLevel
    constexpr unsigned maxCoord = (1u << cstone::maxTreeLevel<cstone::SfcKind<KeyType>>{}) - 1;
    RealType invMaxCoord = RealType(1.0) / maxCoord;
    auto box = getBoundingBox();
    
    return box.zmin() + iz * invMaxCoord * (box.zmax() - box.zmin());
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
MARS_HOST_DEVICE std::tuple<unsigned, unsigned, unsigned> 
ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::sfcToSpatialCoordinate(KeyType sfcKey) const
{
    // Convert raw key to SfcKind strong type
    auto sfcKindKey = cstone::SfcKind<KeyType>(sfcKey);
    auto [ix, iy, iz] = cstone::decodeSfc(sfcKindKey);
    return std::make_tuple(ix, iy, iz);
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
MARS_HOST_DEVICE unsigned 
ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::sfcToSpatialCoordinateX(KeyType sfcKey) const
{
    // Convert raw key to SfcKind strong type
    auto sfcKindKey = cstone::SfcKind<KeyType>(sfcKey);
    auto [ix, iy, iz] = cstone::decodeSfc(sfcKindKey);
    return ix;
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
MARS_HOST_DEVICE unsigned 
ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::sfcToSpatialCoordinateY(KeyType sfcKey) const
{
    // Convert raw key to SfcKind strong type
    auto sfcKindKey = cstone::SfcKind<KeyType>(sfcKey);
    auto [ix, iy, iz] = cstone::decodeSfc(sfcKindKey);
    return iy;
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
MARS_HOST_DEVICE unsigned 
ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::sfcToSpatialCoordinateZ(KeyType sfcKey) const
{
    // Convert raw key to SfcKind strong type
    auto sfcKindKey = cstone::SfcKind<KeyType>(sfcKey);
    auto [ix, iy, iz] = cstone::decodeSfc(sfcKindKey);
    return iz;
}

//host function to get connectivity for a specific element index
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

// Initialize d_conn_keys_ in constructor - missing resize calls
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::initializeConnectivityKeys()
{
    // Initialize all connectivity key vectors based on element type
    if constexpr (std::is_same_v<ElementTag, TetTag>) {
        std::get<0>(d_conn_keys_).resize(elementCount_);
        std::get<1>(d_conn_keys_).resize(elementCount_);
        std::get<2>(d_conn_keys_).resize(elementCount_);
        std::get<3>(d_conn_keys_).resize(elementCount_);
    }
    else if constexpr (std::is_same_v<ElementTag, HexTag>) {
        std::get<0>(d_conn_keys_).resize(elementCount_);
        std::get<1>(d_conn_keys_).resize(elementCount_);
        std::get<2>(d_conn_keys_).resize(elementCount_);
        std::get<3>(d_conn_keys_).resize(elementCount_);
        std::get<4>(d_conn_keys_).resize(elementCount_);
        std::get<5>(d_conn_keys_).resize(elementCount_);
        std::get<6>(d_conn_keys_).resize(elementCount_);
        std::get<7>(d_conn_keys_).resize(elementCount_);
    }
    else if constexpr (std::is_same_v<ElementTag, TriTag>) {
        std::get<0>(d_conn_keys_).resize(elementCount_);
        std::get<1>(d_conn_keys_).resize(elementCount_);
        std::get<2>(d_conn_keys_).resize(elementCount_);
    }
    else if constexpr (std::is_same_v<ElementTag, QuadTag>) {
        std::get<0>(d_conn_keys_).resize(elementCount_);
        std::get<1>(d_conn_keys_).resize(elementCount_);
        std::get<2>(d_conn_keys_).resize(elementCount_);
        std::get<3>(d_conn_keys_).resize(elementCount_);
    }
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::buildAdjacency()
{
    if constexpr (std::is_same_v<AcceleratorTag, cstone::GpuTag>)
    {
        // Step 1: Create the foundational mapping from sparse SFC keys to dense local IDs
        createLocalToGlobalSfcMap();

        // Step 2: Convert the SFC-based element connectivity to a local-ID-based one
        createElementToNodeLocalIdMap();

        // Step 3: Build the inverse mapping (Node -> Element) from the local-ID connectivity
        buildNodeToElementMap();
        
        std::cout << "Rank " << rank_ << " built adjacency maps. Final node count: " << nodeCount_ << std::endl;
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
        if (numConnectivityEntries == 0) {
            nodeCount_ = 0;
            return;
        }
        
        DeviceVector<KeyType> allNodeKeys(numConnectivityEntries);

        auto conn_ptrs = makeConnPtrs<KeyType>(d_conn_keys_, std::make_index_sequence<NodesPerElement>{});
        // 2. Launch the optimized kernel to flatten the tuple-of-vectors
        int blockSize = 256;
        int numBlocks = (elementCount_ + blockSize - 1) / blockSize;
        flattenConnectivityKernel<<<numBlocks, blockSize>>>(
            conn_ptrs,
            thrust::raw_pointer_cast(allNodeKeys.data()), 
            elementCount_
        );
        cudaCheckError();

        // 3. Sort the keys using Thrust's highly optimized parallel sort
        auto start = thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(allNodeKeys.data()));
        auto end = start + allNodeKeys.size();
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

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::buildNodeToElementMap()
{
    if constexpr (std::is_same_v<AcceleratorTag, cstone::GpuTag>)
    {
        size_t numConnectivityEntries = elementCount_ * NodesPerElement;
        if (numConnectivityEntries == 0) {
            d_nodeToElementOffsets_.resize(nodeCount_ + 1);
            thrust::fill(thrust::device, 
                        thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(d_nodeToElementOffsets_.data())),
                        thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(d_nodeToElementOffsets_.data()) + nodeCount_ + 1),
                        KeyType(0));
            d_nodeToElementList_.resize(0);
            return;
        }

        // Step 1: Create flat list of (local_node_id, element_id) pairs
        DeviceVector<KeyType> flatNodeIds(numConnectivityEntries);
        DeviceVector<KeyType> flatElemIds(numConnectivityEntries);

        // Copy all connectivity data to flat arrays
        if constexpr (NodesPerElement >= 1) {
            thrust::copy(thrust::device, 
                        thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(std::get<0>(d_conn_local_ids_).data())), 
                        thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(std::get<0>(d_conn_local_ids_).data()) + elementCount_),
                        thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatNodeIds.data())));
            thrust::sequence(thrust::device, 
                           thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatElemIds.data())), 
                           thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatElemIds.data()) + elementCount_));
        }
        
        if constexpr (NodesPerElement >= 2) {
            size_t offset = 1 * elementCount_;
            thrust::copy(thrust::device, 
                       thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(std::get<1>(d_conn_local_ids_).data())), 
                       thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(std::get<1>(d_conn_local_ids_).data()) + elementCount_),
                       thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatNodeIds.data()) + offset));
            thrust::sequence(thrust::device, 
                           thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatElemIds.data()) + offset), 
                           thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatElemIds.data()) + offset + elementCount_));
        }
        
        if constexpr (NodesPerElement >= 3) {
            size_t offset = 2 * elementCount_;
            thrust::copy(thrust::device, 
                       thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(std::get<2>(d_conn_local_ids_).data())), 
                       thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(std::get<2>(d_conn_local_ids_).data()) + elementCount_),
                       thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatNodeIds.data()) + offset));
            thrust::sequence(thrust::device, 
                           thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatElemIds.data()) + offset), 
                           thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatElemIds.data()) + offset + elementCount_));
        }
        
        if constexpr (NodesPerElement >= 4) {
            size_t offset = 3 * elementCount_;
            thrust::copy(thrust::device, 
                       thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(std::get<3>(d_conn_local_ids_).data())), 
                       thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(std::get<3>(d_conn_local_ids_).data()) + elementCount_),
                       thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatNodeIds.data()) + offset));
            thrust::sequence(thrust::device, 
                           thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatElemIds.data()) + offset), 
                           thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatElemIds.data()) + offset + elementCount_));
        }
        
        if constexpr (NodesPerElement >= 5) {
            size_t offset = 4 * elementCount_;
            thrust::copy(thrust::device, 
                       thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(std::get<4>(d_conn_local_ids_).data())), 
                       thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(std::get<4>(d_conn_local_ids_).data()) + elementCount_),
                       thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatNodeIds.data()) + offset));
            thrust::sequence(thrust::device, 
                           thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatElemIds.data()) + offset), 
                           thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatElemIds.data()) + offset + elementCount_));
        }
        
        if constexpr (NodesPerElement >= 6) {
            size_t offset = 5 * elementCount_;
            thrust::copy(thrust::device, 
                       thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(std::get<5>(d_conn_local_ids_).data())), 
                       thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(std::get<5>(d_conn_local_ids_).data()) + elementCount_),
                       thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatNodeIds.data()) + offset));
            thrust::sequence(thrust::device, 
                           thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatElemIds.data()) + offset), 
                           thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatElemIds.data()) + offset + elementCount_));
        }
        
        if constexpr (NodesPerElement >= 7) {
            size_t offset = 6 * elementCount_;
            thrust::copy(thrust::device, 
                       thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(std::get<6>(d_conn_local_ids_).data())), 
                       thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(std::get<6>(d_conn_local_ids_).data()) + elementCount_),
                       thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatNodeIds.data()) + offset));
            thrust::sequence(thrust::device, 
                           thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatElemIds.data()) + offset), 
                           thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatElemIds.data()) + offset + elementCount_));
        }
        
        if constexpr (NodesPerElement >= 8) {
            size_t offset = 7 * elementCount_;
            thrust::copy(thrust::device, 
                       thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(std::get<7>(d_conn_local_ids_).data())), 
                       thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(std::get<7>(d_conn_local_ids_).data()) + elementCount_),
                       thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatNodeIds.data()) + offset));
            thrust::sequence(thrust::device, 
                           thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatElemIds.data()) + offset), 
                           thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatElemIds.data()) + offset + elementCount_));
        }

        // Step 2: Sort by node ID to group all elements per node
        thrust::sort_by_key(thrust::device, 
                          thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatNodeIds.data())), 
                          thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatNodeIds.data()) + flatNodeIds.size()),
                          thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatElemIds.data())));

        // Step 3: Build CSR structure using reduce_by_key to count entries per node
        DeviceVector<KeyType> uniqueNodeIds(nodeCount_);
        DeviceVector<KeyType> countsPerNode(nodeCount_);
        
        // Count how many times each node appears
        auto new_end = thrust::reduce_by_key(
            thrust::device,
            thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatNodeIds.data())), 
            thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(flatNodeIds.data()) + flatNodeIds.size()),
            thrust::constant_iterator<KeyType>(1),
            thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(uniqueNodeIds.data())),
            thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(countsPerNode.data()))
        );
        
        // Resize to actual size
        uniqueNodeIds.resize(new_end.first - thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(uniqueNodeIds.data())));
        countsPerNode.resize(new_end.second - thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(countsPerNode.data())));

        // Step 4: Build offset array (prefix sum of counts)
        d_nodeToElementOffsets_.resize(nodeCount_ + 1);
        thrust::fill(thrust::device, 
                    thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(d_nodeToElementOffsets_.data())), 
                    thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(d_nodeToElementOffsets_.data()) + d_nodeToElementOffsets_.size()), 
                    KeyType(0));

        // Scatter counts into the offset array
        thrust::scatter(thrust::device,
                       thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(countsPerNode.data())), 
                       thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(countsPerNode.data()) + countsPerNode.size()),
                       thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(uniqueNodeIds.data())),
                       thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(d_nodeToElementOffsets_.data()) + 1));
        
        // Compute prefix sum to get offsets
        thrust::inclusive_scan(thrust::device,
                             thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(d_nodeToElementOffsets_.data())),
                             thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(d_nodeToElementOffsets_.data()) + d_nodeToElementOffsets_.size()),
                             thrust::device_ptr<KeyType>(thrust::raw_pointer_cast(d_nodeToElementOffsets_.data())));

        // Step 5: Copy the element list
        d_nodeToElementList_ = flatElemIds;
    }
}

template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::createElementToNodeLocalIdMap()
{
    if constexpr (std::is_same_v<AcceleratorTag, cstone::GpuTag>)
    {
        // Resize local ID connectivity to match element count
        if constexpr (NodesPerElement >= 1) std::get<0>(d_conn_local_ids_).resize(elementCount_);
        if constexpr (NodesPerElement >= 2) std::get<1>(d_conn_local_ids_).resize(elementCount_);
        if constexpr (NodesPerElement >= 3) std::get<2>(d_conn_local_ids_).resize(elementCount_);
        if constexpr (NodesPerElement >= 4) std::get<3>(d_conn_local_ids_).resize(elementCount_);
        if constexpr (NodesPerElement >= 5) std::get<4>(d_conn_local_ids_).resize(elementCount_);
        if constexpr (NodesPerElement >= 6) std::get<5>(d_conn_local_ids_).resize(elementCount_);
        if constexpr (NodesPerElement >= 7) std::get<6>(d_conn_local_ids_).resize(elementCount_);
        if constexpr (NodesPerElement >= 8) std::get<7>(d_conn_local_ids_).resize(elementCount_);

        int blockSize = 256;
        int numBlocks = (elementCount_ + blockSize - 1) / blockSize;

        // Map each connectivity array from SFC keys to local IDs
        if constexpr (NodesPerElement >= 1) {
            mapSfcToLocalIdKernel<<<numBlocks, blockSize>>>(
                thrust::raw_pointer_cast(std::get<0>(d_conn_keys_).data()),
                thrust::raw_pointer_cast(std::get<0>(d_conn_local_ids_).data()),
                thrust::raw_pointer_cast(d_localToGlobalSfcMap_.data()),
                elementCount_, nodeCount_);
            cudaCheckError();
        }
        if constexpr (NodesPerElement >= 2) {
            mapSfcToLocalIdKernel<<<numBlocks, blockSize>>>(
                thrust::raw_pointer_cast(std::get<1>(d_conn_keys_).data()),
                thrust::raw_pointer_cast(std::get<1>(d_conn_local_ids_).data()),
                thrust::raw_pointer_cast(d_localToGlobalSfcMap_.data()),
                elementCount_, nodeCount_);
            cudaCheckError();
        }
        if constexpr (NodesPerElement >= 3) {
            mapSfcToLocalIdKernel<<<numBlocks, blockSize>>>(
                thrust::raw_pointer_cast(std::get<2>(d_conn_keys_).data()),
                thrust::raw_pointer_cast(std::get<2>(d_conn_local_ids_).data()),
                thrust::raw_pointer_cast(d_localToGlobalSfcMap_.data()),
                elementCount_, nodeCount_);
            cudaCheckError();
        }
        if constexpr (NodesPerElement >= 4) {
            mapSfcToLocalIdKernel<<<numBlocks, blockSize>>>(
                thrust::raw_pointer_cast(std::get<3>(d_conn_keys_).data()),
                thrust::raw_pointer_cast(std::get<3>(d_conn_local_ids_).data()),
                thrust::raw_pointer_cast(d_localToGlobalSfcMap_.data()),
                elementCount_, nodeCount_);
            cudaCheckError();
        }
        if constexpr (NodesPerElement >= 5) {
            mapSfcToLocalIdKernel<<<numBlocks, blockSize>>>(
                thrust::raw_pointer_cast(std::get<4>(d_conn_keys_).data()),
                thrust::raw_pointer_cast(std::get<4>(d_conn_local_ids_).data()),
                thrust::raw_pointer_cast(d_localToGlobalSfcMap_.data()),
                elementCount_, nodeCount_);
            cudaCheckError();
        }
        if constexpr (NodesPerElement >= 6) {
            mapSfcToLocalIdKernel<<<numBlocks, blockSize>>>(
                thrust::raw_pointer_cast(std::get<5>(d_conn_keys_).data()),
                thrust::raw_pointer_cast(std::get<5>(d_conn_local_ids_).data()),
                thrust::raw_pointer_cast(d_localToGlobalSfcMap_.data()),
                elementCount_, nodeCount_);
            cudaCheckError();
        }
        if constexpr (NodesPerElement >= 7) {
            mapSfcToLocalIdKernel<<<numBlocks, blockSize>>>(
                thrust::raw_pointer_cast(std::get<6>(d_conn_keys_).data()),
                thrust::raw_pointer_cast(std::get<6>(d_conn_local_ids_).data()),
                thrust::raw_pointer_cast(d_localToGlobalSfcMap_.data()),
                elementCount_, nodeCount_);
            cudaCheckError();
        }
        if constexpr (NodesPerElement >= 8) {
            mapSfcToLocalIdKernel<<<numBlocks, blockSize>>>(
                thrust::raw_pointer_cast(std::get<7>(d_conn_keys_).data()),
                thrust::raw_pointer_cast(std::get<7>(d_conn_local_ids_).data()),
                thrust::raw_pointer_cast(d_localToGlobalSfcMap_.data()),
                elementCount_, nodeCount_);
            cudaCheckError();
        }
    }
}

} // namespace mars