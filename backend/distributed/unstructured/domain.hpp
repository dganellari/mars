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
#include "cstone/domain/assignment_gpu.cuh"

#include "domain_cuda_impl.hpp"

// stl includes
// #include <adios2.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/extrema.h>
#include <thrust/host_vector.h>
#include <algorithm>
#include <array>
#include <iostream>
#include <memory>
#include <vector>

// Mars includes
#include "mars.hpp"
#include "mars_globals.hpp"
#include "mars_cuda_utils.hpp"
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

// Forward declarations of CUDA kernels - add these with the other declarations
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
template<typename KeyType, typename SfcConnTuple, typename ConnTuple>
void rebuildElementConnectivity(SfcConnTuple& d_conn_keys_, ConnTuple& d_conn_, size_t newElementCount);

// Add this with your other forward declarations around line 120
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
    void calculateCharacteristicSizes();

    // Transfer data to GPU for computations
    void transferDataToGPU(const HostCoordsTuple& h_coords_, const HostConnectivityTuple& h_conn_);

    // Domain synchronization (following cornerstone API pattern)
    void sync();

    // CUDA kernel launcher specialization for GPU tag
    // template<typename KernelFunc>
    // void computeOnElements(KernelFunc kernel)
    // {
    //     int blockSize = 256;
    //     int numBlocks = (elementCount_ + blockSize - 1) / blockSize;

    //     // Extract raw pointers from tuple elements for kernel call
    //     auto& d_x = std::get<0>(d_coords_);
    //     auto& d_y = std::get<1>(d_coords_);
    //     auto& d_z = std::get<2>(d_coords_);

    //     // Get connectivity pointers based on element type
    //     if constexpr (std::is_same_v<ElementTag, TetTag>)
    //     {
    //         auto& d_i0 = std::get<0>(d_conn_);
    //         auto& d_i1 = std::get<1>(d_conn_);
    //         auto& d_i2 = std::get<2>(d_conn_);
    //         auto& d_i3 = std::get<3>(d_conn_);

    //         kernel<<<numBlocks, blockSize>>>(
    //             thrust::raw_pointer_cast(d_x.data()), thrust::raw_pointer_cast(d_y.data()),
    //             thrust::raw_pointer_cast(d_z.data()), thrust::raw_pointer_cast(d_i0.data()),
    //             thrust::raw_pointer_cast(d_i1.data()), thrust::raw_pointer_cast(d_i2.data()),
    //             thrust::raw_pointer_cast(d_i3.data()), elementCount_);
    //         cudaCheckError();
    //     }
    //     else if constexpr (std::is_same_v<ElementTag, HexTag>)
    //     {
    //         auto& d_i0 = std::get<0>(d_conn_);
    //         auto& d_i1 = std::get<1>(d_conn_);
    //         auto& d_i2 = std::get<2>(d_conn_);
    //         auto& d_i3 = std::get<3>(d_conn_);
    //         auto& d_i4 = std::get<4>(d_conn_);
    //         auto& d_i5 = std::get<5>(d_conn_);
    //         auto& d_i6 = std::get<6>(d_conn_);
    //         auto& d_i7 = std::get<7>(d_conn_);

    //         kernel<<<numBlocks, blockSize>>>(
    //             thrust::raw_pointer_cast(d_x.data()), thrust::raw_pointer_cast(d_y.data()),
    //             thrust::raw_pointer_cast(d_z.data()), thrust::raw_pointer_cast(d_i0.data()),
    //             thrust::raw_pointer_cast(d_i1.data()), thrust::raw_pointer_cast(d_i2.data()),
    //             thrust::raw_pointer_cast(d_i3.data()), thrust::raw_pointer_cast(d_i4.data()),
    //             thrust::raw_pointer_cast(d_i5.data()), thrust::raw_pointer_cast(d_i6.data()),
    //             thrust::raw_pointer_cast(d_i7.data()), elementCount_);
    //         cudaCheckError();
    //     }
    //     // Add cases for other element types
    // }

    // Access methods for coordinate data
    const DeviceVector<RealType>& x() const { return std::get<0>(d_coords_); }
    const DeviceVector<RealType>& y() const { return std::get<1>(d_coords_); }
    const DeviceVector<RealType>& z() const { return std::get<2>(d_coords_); }

    // Access to connectivity - template parameter to select index array
    template<int I>
    const DeviceVector<KeyType>& indices() const
    {
        return std::get<I>(d_conn_);
    }

    // Get element/node counts
    std::size_t getNodeCount() const { return nodeCount_; }
    std::size_t getElementCount() const { return elementCount_; }

    // Access to cornerstone domain
    DomainType& getDomain() { return *domain_; }

    // Start and end indices for local work assignment
    std::size_t startIndex() const { return domain_->startIndex(); }
    std::size_t endIndex() const { return domain_->endIndex(); }

    // Helper methods
    void readMeshDataSoA(const std::string& meshFile, HostCoordsTuple& h_coords_, HostConnectivityTuple& h_conn_);

private:
    int rank_;
    int numRanks_;
    std::unique_ptr<DomainType> domain_;
    std::size_t elementCount_ = 0;
    std::size_t nodeCount_    = 0;

    // Device data in SoA format
    DeviceCoordsTuple d_coords_;     // (x, y, z)
    DevicePropsTuple d_props_;       // (h)
    DeviceConnectivityTuple d_conn_; // (i0, i1, i2, ...) depends on element type
    DeviceConnectivityTuple
        d_conn_keys_; // (sfc_i0, sfc_i1, sfc_i2, ...) depends on element type; store sfc instead of coordinates
    DeviceVector<KeyType> d_elemToNodeMap_;
    DeviceVector<KeyType> d_elemSfcCodes_;

    void syncImpl(DeviceVector<KeyType>& d_nodeSfcCodes, int blockSize, int numBlocks)
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

// Create a bounding box from a coordinate tuple
template<typename RealType, typename CoordTuple>
cstone::Box<RealType> createBoundingBox(const CoordTuple& coords, RealType padding = 0.05)
{
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

// Element domain constructor
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::ElementDomain(const std::string& meshFile,
                                                                            int rank,
                                                                            int numRanks)
    : rank_(rank)
    , numRanks_(numRanks)
{
    // Host data in SoA format
    HostCoordsTuple h_coords;     // (x, y, z)
    HostConnectivityTuple h_conn; // (i0, i1, i2, ...) depends on element type

    // Read the mesh in SoA format
    readMeshDataSoA(meshFile, h_coords, h_conn);

    // Initialize cornerstone domain
    int bucketSize           = 64;
    unsigned bucketSizeFocus = 8;
    RealType theta           = 0.5;

    // Create a bounding box with some padding
    auto box = createBoundingBox<RealType>(h_coords);

    domain_ = std::make_unique<DomainType>(rank, numRanks, bucketSize, bucketSizeFocus, theta, box);

    // Transfer data to GPU before sync
    transferDataToGPU(h_coords, h_conn);

    // Calculate characteristic sizes
    calculateCharacteristicSizes();

    // Perform sync
    sync();

    std::cout << "Rank " << rank_ << " initialized " << ElementTag::Name << " domain with " << nodeCount_
              << " nodes and " << elementCount_ << " elements." << std::endl;
}

// Read mesh data in SoA format; uses element-based partitioning for better data locality
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::readMeshDataSoA(const std::string& meshFile,
                                                                                   HostCoordsTuple& h_coords_,
                                                                                   HostConnectivityTuple& h_conn_)
{
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
void ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::calculateCharacteristicSizes()
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
void ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::sync()
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
                                           thrust::raw_pointer_cast(d_nodeSfcCodes.data()), nodeCount_, domain_->box());

        // Find representative nodes for each element
        d_elemToNodeMap_.resize(elementCount_);
        int blockSize = 256;
        int numBlocks = (elementCount_ + blockSize - 1) / blockSize;

        // Call element-specific implementation
        syncImpl(d_nodeSfcCodes, blockSize, numBlocks);

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

        syncDomainImpl(domain_.get(), d_elemSfcCodes_, d_elemX, d_elemY, d_elemZ, d_elemH, elementCount_, d_conn_keys_);

        size_t newElementCount = domain_->endIndex() - domain_->startIndex();
        elementCount_          = newElementCount;

        if (newElementCount == 0) {
            return;
        }
        
        rebuildElementConnectivity<KeyType, DeviceConnectivityTuple, DeviceConnectivityTuple>(d_conn_keys_, d_conn_, newElementCount);
    }
}

// Transfer data to GPU for computations
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
void ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>::transferDataToGPU(
    const HostCoordsTuple& h_coords_, const HostConnectivityTuple& h_conn_)
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
} // namespace mars
