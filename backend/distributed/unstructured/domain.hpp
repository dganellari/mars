#pragma once

#include <tuple>
#include <utility>

// Use std::get for tuple access
using std::get;

// Cornerstone includes
#include "cstone/domain/domain.hpp"
#include <adios2.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/extrema.h>
#include <thrust/host_vector.h>
#include <algorithm>
#include <array>
#include <iostream>
#include <memory>
#include <vector>
#include "mars.hpp"
#include "mars_globals.hpp"
#include "mars_read_mesh_adios2.hpp"
#include "mars_read_mesh_binary.hpp"

using Real    = float;
using KeyType = unsigned;

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
template<typename ElementTag>
__global__ void findRepresentativeNodesKernel(const int* indices0,
                                              const int* indices1,
                                              const int* indices2,
                                              const int* indices3,
                                              const unsigned* sfcCodes,
                                              int* elemToNodeMap,
                                              int numElements);

template<typename ElementTag>
__global__ void extractRepCoordinatesKernel(const Real* x,
                                            const Real* y,
                                            const Real* z,
                                            const Real* h,
                                            const int* elemToNodeMap,
                                            Real* elemX,
                                            Real* elemY,
                                            Real* elemZ,
                                            Real* elemH,
                                            int numElements);

// Add these with the other forward declarations
template<typename ElementTag>
__global__ void computeCharacteristicSizesKernel(const Real* x,
                                                 const Real* y,
                                                 const Real* z,
                                                 const int* indices0,
                                                 const int* indices1,
                                                 const int* indices2,
                                                 const int* indices3,
                                                 int* nodeTetCount,
                                                 Real* h,
                                                 int numElements);

__global__ void finalizeCharacteristicSizesKernel(Real* h, int* nodeTetCount, int numNodes);

// After other forward declarations
template<typename KeyType>
void computeSfcKeysGpu(
    const Real* x, const Real* y, const Real* z, KeyType* keys, size_t numKeys, const cstone::Box<Real>& box);

// Template struct to select the correct vector type based on accelerator tag
template<typename T, typename AcceleratorTag>
struct VectorSelector
{
    // Default to host vector for any tag
    using type = std::vector<T>;
};

// Specialization for GPU tag - use device vector
template<typename T>
struct VectorSelector<T, cstone::GpuTag>
{
    using type = thrust::device_vector<T>;
};

// Helper for creating element connectivity type
template<typename ElementTag, typename AcceleratorTag>
struct ConnectivityTupleHelper;

// Specialization for tetrahedra
template<typename AcceleratorTag>
struct ConnectivityTupleHelper<TetTag, AcceleratorTag>
{
    template<typename T>
    using DeviceVector = typename VectorSelector<T, AcceleratorTag>::type;

    using type = std::tuple<DeviceVector<int>, DeviceVector<int>, DeviceVector<int>, DeviceVector<int>>;
};

// Specialization for hexahedra
template<typename AcceleratorTag>
struct ConnectivityTupleHelper<HexTag, AcceleratorTag>
{
    template<typename T>
    using DeviceVector = typename VectorSelector<T, AcceleratorTag>::type;

    using type = std::tuple<DeviceVector<int>,
                            DeviceVector<int>,
                            DeviceVector<int>,
                            DeviceVector<int>,
                            DeviceVector<int>,
                            DeviceVector<int>,
                            DeviceVector<int>,
                            DeviceVector<int>>;
};

// Specializations for triangle and quad would be similar
template<typename AcceleratorTag>
struct ConnectivityTupleHelper<TriTag, AcceleratorTag>
{
    template<typename T>
    using DeviceVector = typename VectorSelector<T, AcceleratorTag>::type;

    using type = std::tuple<DeviceVector<int>, DeviceVector<int>, DeviceVector<int>>;
};

template<typename AcceleratorTag>
struct ConnectivityTupleHelper<QuadTag, AcceleratorTag>
{
    template<typename T>
    using DeviceVector = typename VectorSelector<T, AcceleratorTag>::type;

    using type = std::tuple<DeviceVector<int>, DeviceVector<int>, DeviceVector<int>, DeviceVector<int>>;
};

// Main domain class templated on element type
template<typename ElementTag = TetTag, typename AcceleratorTag = cstone::GpuTag>
class ElementDomain
{
public:
    static constexpr int NodesPerElement = ElementTag::NodesPerElement;

    using DomainType = cstone::Domain<KeyType, Real, AcceleratorTag>;

    // Template alias for appropriate vector type
    template<typename T>
    using DeviceVector = typename VectorSelector<T, AcceleratorTag>::type;

    // SoA data structures using tuples
    using CoordsTuple       = std::tuple<DeviceVector<Real>, DeviceVector<Real>, DeviceVector<Real>>;
    using PropsTuple        = std::tuple<DeviceVector<Real>>;
    using ConnectivityTuple = typename ConnectivityTupleHelper<ElementTag, AcceleratorTag>::type;

    ElementDomain(const std::string& meshFile, int rank, int numRanks);

    // GPU-accelerated calculation of characteristic sizes
    void calculateCharacteristicSizes();

    // Map elements to representative nodes using SFC
    void mapElementsToNodes();

    // Transfer data to GPU for computations
    void transferDataToGPU(CoordsTuple& h_coords_, ConnectivityTuple& h_conn_);

    // Domain synchronization (following cornerstone API pattern)
    void sync();

    // Get all elements in a given octree node
    // std::vector<size_t> getElementsInOctreeNode(int octreeNodeIndex);

    // Launch a CUDA kernel to compute something on elements
    template<typename KernelFunc>
    void computeOnElements(KernelFunc kernel);

    // Access methods for coordinate data
    const DeviceVector<Real>& x() const { return std::get<0>(d_coords_); }
    const DeviceVector<Real>& y() const { return std::get<1>(d_coords_); }
    const DeviceVector<Real>& z() const { return std::get<2>(d_coords_); }

    // Access to connectivity - template parameter to select index array
    template<int I>
    const DeviceVector<int>& indices() const
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

private:
    int rank_;
    int numRanks_;
    std::unique_ptr<DomainType> domain_;
    std::size_t elementCount_ = 0;
    std::size_t nodeCount_    = 0;

    // Device data in SoA format
    CoordsTuple d_coords_;     // (x, y, z)
    PropsTuple d_props_;       // (h)
    ConnectivityTuple d_conn_; // (i0, i1, i2, ...) depends on element type
    DeviceVector<int> d_elemToNodeMap_;
    DeviceVector<KeyType> d_elemSfcCodes_;

    // Helper methods
    void readMeshDataSoA(const std::string& meshFile, CoordsTuple& h_coords_, ConnectivityTuple& h_conn_);
};

// CUDA kernel launcher specialization for GPU tag
template<typename ElementTag>
template<typename KernelFunc>
void ElementDomain<ElementTag, cstone::GpuTag>::computeOnElements(KernelFunc kernel)
{
    int blockSize = 256;
    int numBlocks = (elementCount_ + blockSize - 1) / blockSize;

    // Extract raw pointers from tuple elements for kernel call
    auto& d_x = std::get<0>(d_coords_);
    auto& d_y = std::get<1>(d_coords_);
    auto& d_z = std::get<2>(d_coords_);

    // Get connectivity pointers based on element type
    if constexpr (std::is_same_v<ElementTag, TetTag>)
    {
        auto& d_i0 = std::get<0>(d_conn_);
        auto& d_i1 = std::get<1>(d_conn_);
        auto& d_i2 = std::get<2>(d_conn_);
        auto& d_i3 = std::get<3>(d_conn_);

        kernel<<<numBlocks, blockSize>>>(thrust::raw_pointer_cast(d_x.data()), thrust::raw_pointer_cast(d_y.data()),
                                         thrust::raw_pointer_cast(d_z.data()), thrust::raw_pointer_cast(d_i0.data()),
                                         thrust::raw_pointer_cast(d_i1.data()), thrust::raw_pointer_cast(d_i2.data()),
                                         thrust::raw_pointer_cast(d_i3.data()), elementCount_);
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

        kernel<<<numBlocks, blockSize>>>(thrust::raw_pointer_cast(d_x.data()), thrust::raw_pointer_cast(d_y.data()),
                                         thrust::raw_pointer_cast(d_z.data()), thrust::raw_pointer_cast(d_i0.data()),
                                         thrust::raw_pointer_cast(d_i1.data()), thrust::raw_pointer_cast(d_i2.data()),
                                         thrust::raw_pointer_cast(d_i3.data()), thrust::raw_pointer_cast(d_i4.data()),
                                         thrust::raw_pointer_cast(d_i5.data()), thrust::raw_pointer_cast(d_i6.data()),
                                         thrust::raw_pointer_cast(d_i7.data()), elementCount_);
        cudaCheckError();
    }
    // Add cases for other element types
}

// CPU version that doesn't use CUDA syntax
template<typename ElementTag, typename AcceleratorTag>
template<typename KernelFunc>
void ElementDomain<ElementTag, AcceleratorTag>::computeOnElements(KernelFunc kernel)
{
    // CPU implementation - this would need to be implemented differently
    // Perhaps a serial loop or OpenMP parallelization
    std::cerr << "Warning: computeOnElements not implemented for CPU mode\n";
}

// Create a bounding box from a coordinate tuple
template<typename CoordTuple>
cstone::Box<Real> createBoundingBox(const CoordTuple& coords, Real padding = 0.05)
{
    // Get references to coordinate vectors
    const auto& x = std::get<0>(coords);
    const auto& y = std::get<1>(coords);
    const auto& z = std::get<2>(coords);

    // Find min/max coordinates
    auto minmax_x = std::minmax_element(x.begin(), x.end());
    auto minmax_y = std::minmax_element(y.begin(), y.end());
    auto minmax_z = std::minmax_element(z.begin(), z.end());

    // Create box with padding
    Real x_min = *minmax_x.first - padding;
    Real x_max = *minmax_x.second + padding;
    Real y_min = *minmax_y.first - padding;
    Real y_max = *minmax_y.second + padding;
    Real z_min = *minmax_z.first - padding;
    Real z_max = *minmax_z.second + padding;

    return cstone::Box<Real>(x_min, x_max, y_min, y_max, z_min, z_max);
}

// Element domain constructor
template<typename ElementTag, typename AcceleratorTag>
ElementDomain<ElementTag, AcceleratorTag>::ElementDomain(const std::string& meshFile, int rank, int numRanks)
    : rank_(rank)
    , numRanks_(numRanks)
{
    // Host data in SoA format
    CoordsTuple h_coords_;     // (x, y, z)
    ConnectivityTuple h_conn_; // (i0, i1, i2, ...) depends on element type
    // Read the mesh in SoA format
    readMeshDataSoA(meshFile, h_coords_, h_conn_);

    // Initialize cornerstone domain
    int bucketSize           = 64;
    unsigned bucketSizeFocus = 8;
    Real theta               = 0.5f;

    // Create a bounding box with some padding
    auto box = createBoundingBox(h_coords_);

    domain_ = std::make_unique<DomainType>(rank, numRanks, bucketSize, bucketSizeFocus, theta, box);

    // Transfer data to GPU before sync
    transferDataToGPU(h_coords_, h_conn_);

    // Calculate characteristic sizes
    calculateCharacteristicSizes();

    // Map elements to their representative nodes
    mapElementsToNodes();

    // Perform sync
    sync();

    std::cout << "Rank " << rank_ << " initialized " << ElementTag::Name << " domain with " << nodeCount_
              << " nodes and " << elementCount_ << " elements." << std::endl;
}

// Read mesh data in SoA format; assumes float32 for coordinates and int32 for indices
/* template <typename ElementTag, typename AcceleratorTag>
void ElementDomain<ElementTag, AcceleratorTag>::readMeshDataSoA(const std::string& meshFile, CoordsTuple& h_coords_,
                                                               ConnectivityTuple& h_conn_) {
    // Get references to coordinate vectors for clarity
    auto& h_x = std::get<0>(h_coords_);
    auto& h_y = std::get<1>(h_coords_);
    auto& h_z = std::get<2>(h_coords_);

    try {
        // Open the mesh directory and read coordinate files
        std::ifstream x_file(meshFile + "/x.float32", std::ios::binary);
        std::ifstream y_file(meshFile + "/y.float32", std::ios::binary);
        std::ifstream z_file(meshFile + "/z.float32", std::ios::binary);

        if (!x_file || !y_file || !z_file) {
            throw std::runtime_error("Failed to open coordinate files");
        }

        // Get file size to determine node count
        x_file.seekg(0, std::ios::end);
        size_t x_size = x_file.tellg() / sizeof(float);
        x_file.seekg(0, std::ios::beg);

        // Calculate this rank's portion
        size_t nodePerRank = x_size / numRanks_;
        size_t nodeStartIdx = rank_ * nodePerRank;
        size_t nodeEndIdx = (rank_ == numRanks_ - 1) ? x_size : nodeStartIdx + nodePerRank;
        nodeCount_ = nodeEndIdx - nodeStartIdx;

        // Read coordinate data
        std::vector<float> x_temp(nodeCount_), y_temp(nodeCount_), z_temp(nodeCount_);

        x_file.seekg(nodeStartIdx * sizeof(float));
        y_file.seekg(nodeStartIdx * sizeof(float));
        z_file.seekg(nodeStartIdx * sizeof(float));

        x_file.read(reinterpret_cast<char*>(x_temp.data()), nodeCount_ * sizeof(float));
        y_file.read(reinterpret_cast<char*>(y_temp.data()), nodeCount_ * sizeof(float));
        z_file.read(reinterpret_cast<char*>(z_temp.data()), nodeCount_ * sizeof(float));

        // Convert float to Real if needed
        h_x.resize(nodeCount_);
        h_y.resize(nodeCount_);
        h_z.resize(nodeCount_);

        for (size_t i = 0; i < nodeCount_; ++i) {
            h_x[i] = static_cast<Real>(x_temp[i]);
            h_y[i] = static_cast<Real>(y_temp[i]);
            h_z[i] = static_cast<Real>(z_temp[i]);
        }

        // Read connectivity files - one for each node index in the element
        // Open all index files based on element type
        std::vector<std::ifstream> index_files;
        for (int i = 0; i < NodesPerElement; ++i) {
            index_files.emplace_back(meshFile + "/i" + std::to_string(i) + ".int32", std::ios::binary);
            if (!index_files.back()) {
                throw std::runtime_error("Failed to open index file i" + std::to_string(i));
            }
        }

        // Get element count
        index_files[0].seekg(0, std::ios::end);
        size_t elem_size = index_files[0].tellg() / sizeof(int32_t);
        index_files[0].seekg(0, std::ios::beg);

        // Calculate this rank's element portion
        size_t elemPerRank = elem_size / numRanks_;
        size_t elemStartIdx = rank_ * elemPerRank;
        size_t elemEndIdx = (rank_ == numRanks_ - 1) ? elem_size : elemStartIdx + elemPerRank;
        elementCount_ = elemEndIdx - elemStartIdx;

        // Helper function to expand tuple
        auto resizeConnectivity = [this](auto& tuple, size_t size) {
            std::apply([size](auto&... vecs) { ((vecs.resize(size)), ...); }, tuple);
        };

        // Resize connectivity arrays
        resizeConnectivity(h_conn_, elementCount_);

        // Read element connectivity
        std::vector<int32_t> temp_indices(elementCount_);

        for (int i = 0; i < NodesPerElement; ++i) {
            index_files[i].seekg(elemStartIdx * sizeof(int32_t));
            index_files[i].read(reinterpret_cast<char*>(temp_indices.data()), elementCount_ * sizeof(int32_t));

            // Adjust indices for local numbering
            for (size_t j = 0; j < elementCount_; ++j) {
                temp_indices[j] -= nodeStartIdx;
                // Store in appropriate tuple element
                std::get<i>(h_conn_)[j] = temp_indices[j];
            }
        }

    } catch (const std::exception& e) {
        std::cerr << "Error reading mesh: " << e.what() << std::endl;
        throw;
    }
} */

// Read mesh data in SoA format; assumes float32 for coordinates and int32 for indices
template<typename ElementTag, typename AcceleratorTag>
void ElementDomain<ElementTag, AcceleratorTag>::readMeshDataSoA(const std::string& meshFile,
                                                                CoordsTuple& h_coords_,
                                                                ConnectivityTuple& h_conn_)
{
    // Get references to coordinate vectors for clarity
    auto& h_x = std::get<0>(h_coords_);
    auto& h_y = std::get<1>(h_coords_);
    auto& h_z = std::get<2>(h_coords_);

    try
    {
        // Read coordinates using utility function
        auto [readNodeCount, nodeStartIdx, x_data, y_data, z_data] =
            mars::readMeshCoordinatesBinary(meshFile, rank_, numRanks_);

        nodeCount_ = readNodeCount;

        // If Real is float, we can move directly
        if constexpr (std::is_same_v<Real, float>)
        {
            std::get<0>(h_coords_) = std::move(x_data);
            std::get<1>(h_coords_) = std::move(y_data);
            std::get<2>(h_coords_) = std::move(z_data);
        }
        else
        {
            // Otherwise convert
            h_x.resize(nodeCount_);
            h_y.resize(nodeCount_);
            h_z.resize(nodeCount_);

            for (size_t i = 0; i < nodeCount_; ++i)
            {
                h_x[i] = static_cast<Real>(x_data[i]);
                h_y[i] = static_cast<Real>(y_data[i]);
                h_z[i] = static_cast<Real>(z_data[i]);
            }
        }

        // Read connectivity using tuple-specific template function
        if constexpr (std::is_same_v<ElementTag, TetTag>)
        {
            auto [readElementCount, conn_tuple] =
                mars::readMeshConnectivityBinaryTuple<4>(meshFile, nodeStartIdx, rank_, numRanks_);
            elementCount_ = readElementCount;
            h_conn_       = std::move(conn_tuple); // Direct move - no copying needed
        }
        else if constexpr (std::is_same_v<ElementTag, HexTag>)
        {
            auto [readElementCount, conn_tuple] =
                mars::readMeshConnectivityBinaryTuple<8>(meshFile, nodeStartIdx, rank_, numRanks_);
            elementCount_ = readElementCount;
            h_conn_       = std::move(conn_tuple); // Direct move - no copying needed
        }
        else if constexpr (std::is_same_v<ElementTag, TriTag>)
        {
            auto [readElementCount, conn_tuple] =
                mars::readMeshConnectivityBinaryTuple<3>(meshFile, nodeStartIdx, rank_, numRanks_);
            elementCount_ = readElementCount;
            h_conn_       = std::move(conn_tuple);
        }
        else if constexpr (std::is_same_v<ElementTag, QuadTag>)
        {
            auto [readElementCount, conn_tuple] =
                mars::readMeshConnectivityBinaryTuple<4>(meshFile, nodeStartIdx, rank_, numRanks_);
            elementCount_ = readElementCount;
            h_conn_       = std::move(conn_tuple);
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error reading mesh: " << e.what() << std::endl;
        throw;
    }
}

// Calculate characteristic sizes for elements - GPU-only version
template<typename ElementTag, typename AcceleratorTag>
void ElementDomain<ElementTag, AcceleratorTag>::calculateCharacteristicSizes()
{
    // Work directly with device data
    auto& d_h = std::get<0>(d_props_);
    d_h.resize(nodeCount_, 0.0);

    // For tetrahedra, compute based on average edge lengths
    if constexpr (std::is_same_v<ElementTag, TetTag> && std::is_same_v<AcceleratorTag, cstone::GpuTag>)
    {
        // Device vectors for calculation
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

        computeCharacteristicSizesKernel<TetTag><<<numBlocks, blockSize>>>(
            thrust::raw_pointer_cast(d_x.data()), thrust::raw_pointer_cast(d_y.data()),
            thrust::raw_pointer_cast(d_z.data()), thrust::raw_pointer_cast(d_i0.data()),
            thrust::raw_pointer_cast(d_i1.data()), thrust::raw_pointer_cast(d_i2.data()),
            thrust::raw_pointer_cast(d_i3.data()), thrust::raw_pointer_cast(d_nodeTetCount.data()),
            thrust::raw_pointer_cast(d_h.data()), elementCount_);

        cudaCheckError();

        // Then normalize by number of contributions
        numBlocks = (nodeCount_ + blockSize - 1) / blockSize;
        finalizeCharacteristicSizesKernel<<<numBlocks, blockSize>>>(
            thrust::raw_pointer_cast(d_h.data()), thrust::raw_pointer_cast(d_nodeTetCount.data()), nodeCount_);

        cudaCheckError();

        // For mesh-based methods (FEM/FDM)
        constexpr Real meshFactor = 1.0;    // No reduction for FEM (adjust based on element order)
        constexpr Real minH       = 1.0e-6; // Prevent extremely small values that cause instability
        constexpr Real maxH       = 1.0;    // Upper bound based on problem domain

        thrust::transform(d_h.begin(), d_h.end(), d_h.begin(),
                          [=] __device__(Real val)
                          {
                              // Apply mesh factor (might be 1.0 for linear elements)
                              Real result = val * meshFactor;
                              // Apply min/max constraints for numerical stability
                              return max(minH, min(maxH, result));
                          });
    }
    else
    {
        // For non-GPU or non-tetrahedral cases, use a default value
        Real domainDiagonal = std::sqrt(std::pow(domain_->box().xmax() - domain_->box().xmin(), 2) +
                                        std::pow(domain_->box().ymax() - domain_->box().ymin(), 2) +
                                        std::pow(domain_->box().zmax() - domain_->box().zmin(), 2));
        Real defaultH       = domainDiagonal * 0.01; // 1% of domain diagonal
        thrust::fill(d_h.begin(), d_h.end(), defaultH);
    }
}

// Map elements to representative nodes
template<typename ElementTag, typename AcceleratorTag>
void ElementDomain<ElementTag, AcceleratorTag>::mapElementsToNodes()
{
    if constexpr (std::is_same_v<AcceleratorTag, cstone::GpuTag>)
    {
        // Get SFC codes from cornerstone domain or compute them
        DeviceVector<KeyType> d_nodeSfcCodes(nodeCount_);

        auto& d_x = std::get<0>(d_coords_);
        auto& d_y = std::get<1>(d_coords_);
        auto& d_z = std::get<2>(d_coords_);

        // Compute SFC codes on GPU
        computeSfcKeysGpu(thrust::raw_pointer_cast(d_x.data()), thrust::raw_pointer_cast(d_y.data()),
                          thrust::raw_pointer_cast(d_z.data()), thrust::raw_pointer_cast(d_nodeSfcCodes.data()),
                          nodeCount_, domain_->box());

        // Allocate element-to-node mapping
        d_elemToNodeMap_.resize(elementCount_);

        // Launch kernel to find representative nodes
        int blockSize = 256;
        int numBlocks = (elementCount_ + blockSize - 1) / blockSize;

        // Launch specialized kernel based on element type
        if constexpr (std::is_same_v<ElementTag, TetTag>)
        {
            auto& d_i0 = std::get<0>(d_conn_);
            auto& d_i1 = std::get<1>(d_conn_);
            auto& d_i2 = std::get<2>(d_conn_);
            auto& d_i3 = std::get<3>(d_conn_);

            findRepresentativeNodesKernel<TetTag>
                <<<numBlocks, blockSize>>>(thrust::raw_pointer_cast(d_i0.data()), thrust::raw_pointer_cast(d_i1.data()),
                                           thrust::raw_pointer_cast(d_i2.data()), thrust::raw_pointer_cast(d_i3.data()),
                                           thrust::raw_pointer_cast(d_nodeSfcCodes.data()),
                                           thrust::raw_pointer_cast(d_elemToNodeMap_.data()), elementCount_);
            cudaCheckError();
        }
        else if constexpr (std::is_same_v<ElementTag, HexTag>)
        {
            // Similar implementation for hex elements
        }
    }
}

// Transfer data to GPU for computations
template<typename ElementTag, typename AcceleratorTag>
void ElementDomain<ElementTag, AcceleratorTag>::transferDataToGPU(CoordsTuple& h_coords_, ConnectivityTuple& h_conn_)
{
    if constexpr (std::is_same_v<AcceleratorTag, cstone::GpuTag>)
    {
        // Copy coordinates
        std::get<0>(d_coords_) = std::get<0>(h_coords_); // x
        std::get<1>(d_coords_) = std::get<1>(h_coords_); // y
        std::get<2>(d_coords_) = std::get<2>(h_coords_); // z
        // Copy connectivity
        for (int i = 0; i < NodesPerElement; ++i)
        {
            std::get<i>(d_conn_) = std::get<i>(h_conn_);
        }
    }
}

// Domain synchronization
template<typename ElementTag, typename AcceleratorTag>
void ElementDomain<ElementTag, AcceleratorTag>::sync()
{
    if constexpr (std::is_same_v<AcceleratorTag, cstone::GpuTag>)
    {
        // Extract coordinates for clarity
        auto& d_x = std::get<0>(d_coords_);
        auto& d_y = std::get<1>(d_coords_);
        auto& d_z = std::get<2>(d_coords_);
        auto& d_h = std::get<0>(d_props_);

        // Extract element connectivity
        auto getConnIndices = [this](int i) -> DeviceVector<int>& { return std::get<i>(d_conn_); };

        // Calculate SFC codes for nodes
        DeviceVector<KeyType> d_nodeSfcCodes(nodeCount_);
        computeSfcKeysGpu(thrust::raw_pointer_cast(d_x.data()), thrust::raw_pointer_cast(d_y.data()),
                          thrust::raw_pointer_cast(d_z.data()), thrust::raw_pointer_cast(d_nodeSfcCodes.data()),
                          nodeCount_, domain_->box());

        // Find representative nodes for each element
        d_elemToNodeMap_.resize(elementCount_);
        int blockSize = 256;
        int numBlocks = (elementCount_ + blockSize - 1) / blockSize;

        if constexpr (std::is_same_v<ElementTag, TetTag>)
        {
            findRepresentativeNodesKernel<TetTag><<<numBlocks, blockSize>>>(
                thrust::raw_pointer_cast(getConnIndices(0).data()), thrust::raw_pointer_cast(getConnIndices(1).data()),
                thrust::raw_pointer_cast(getConnIndices(2).data()), thrust::raw_pointer_cast(getConnIndices(3).data()),
                thrust::raw_pointer_cast(d_nodeSfcCodes.data()), thrust::raw_pointer_cast(d_elemToNodeMap_.data()),
                elementCount_);
            cudaCheckError();
        }
        else if constexpr (std::is_same_v<ElementTag, HexTag>)
        {
            // Similar implementation for hexes
        }

        // Extract coordinates of representative nodes
        DeviceVector<Real> d_elemX(elementCount_);
        DeviceVector<Real> d_elemY(elementCount_);
        DeviceVector<Real> d_elemZ(elementCount_);
        DeviceVector<Real> d_elemH(elementCount_);
        d_elemSfcCodes_.resize(elementCount_);

        extractRepCoordinatesKernel<ElementTag><<<numBlocks, blockSize>>>(
            thrust::raw_pointer_cast(d_x.data()), thrust::raw_pointer_cast(d_y.data()),
            thrust::raw_pointer_cast(d_z.data()), thrust::raw_pointer_cast(d_h.data()),
            thrust::raw_pointer_cast(d_elemToNodeMap_.data()), thrust::raw_pointer_cast(d_elemX.data()),
            thrust::raw_pointer_cast(d_elemY.data()), thrust::raw_pointer_cast(d_elemZ.data()),
            thrust::raw_pointer_cast(d_elemH.data()), elementCount_);
        cudaCheckError();
        // Sync domain with element representatives
        DeviceVector<Real> scratch1, scratch2, scratch3;
        domain_->sync(d_elemSfcCodes_, d_elemX, d_elemY, d_elemZ, d_elemH, std::tie(d_elemToNodeMap_),
                      std::tie(scratch1, scratch2, scratch3));

        size_t localElements = domain_->endIndex() - domain_->startIndex();
        size_t totalElements = domain_->nParticles();
    }
    else
    {
        // CPU implementation - this would need to be implemented differently
        // Perhaps a serial loop or OpenMP parallelization
        std::cerr << "Warning: sync not implemented for CPU mode\n";

        // Would implement similar logic but with host vectors
    }
}
