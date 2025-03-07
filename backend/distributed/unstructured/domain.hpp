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

using Real = double;
using KeyType = unsigned;

// Forward declarations of CUDA kernels
__global__ void computeCharacteristicSizesKernel(double* nodes,
                                                 int* tets,
                                                 int* nodeTetCount,
                                                 double* h,
                                                 int localTetCount);

__global__ void finalizeCharacteristicSizesKernel(double* h, int* nodeTetCount, int numNodes);

__global__ void findRepresentativeNodesKernel(int* tets, unsigned* sfcCodes, int* tetToNodeMap, int numTets);

__global__ void computeTetrahedralVolumesKernel(double* x, double* y, double* z, int* connectivity, int numTets);

// Template struct to select the correct vector type based on accelerator tag
template<typename T, typename AcceleratorTag>
struct VectorSelector {
    // Default to host vector for any tag
    using type = std::vector<T>;
};

// Specialization for GPU tag - use device vector
template<typename T>
struct VectorSelector<T, cstone::GpuTag> {
    using type = thrust::device_vector<T>;
};

// TetrahedralDomain class with template support for CPU/GPU
template <typename AcceleratorTag = cstone::GpuTag>
class TetrahedralDomain {
public:
    using DomainType = cstone::Domain<KeyType, Real, AcceleratorTag>;

    // Template alias for appropriate vector type
    template<typename T>
    using DeviceVector = typename VectorSelector<T, AcceleratorTag>::type;

    TetrahedralDomain(const std::string& meshFile, int rank, int numRanks);

    // GPU-accelerated calculation of characteristic sizes
    void calculateCharacteristicSizes(const std::vector<double>& hostNodes,
                                      const std::vector<int>& hostTets,
                                      std::size_t localNodeCount,
                                      std::size_t localTetCount);

    // Map tetrahedra to representative nodes using SFC
    void mapTetrahedraToNodes();

    // Transfer data to GPU for computations
    void transferDataToGPU();

    // Domain synchronization (following cornerstone API pattern)
    void sync();

    /* // Exchange halos for field data
    template <typename T, class Vector = cstone::DeviceVector<T>>
    void exchangeHalos(Vector<T>& field) {
        Vector<T> sendBuffer, receiveBuffer;
        // Create a tuple of references as required by cornerstone
        domain_->exchangeHalos(std::tie(field), sendBuffer, receiveBuffer);
    } */

    // Get all tetrahedra in a given octree node
    std::vector<size_t> getTetrahedraInOctreeNode(int octreeNodeIndex);

    // Launch a CUDA kernel to compute something on tetrahedra
    template <typename KernelFunc>
    void computeOnTetrahedra(KernelFunc kernel);

    // Access methods for host data
    const std::vector<Real>& getX() const { return h_x_; }
    const std::vector<Real>& getY() const { return h_y_; }
    const std::vector<Real>& getZ() const { return h_z_; }
    std::size_t getNodeCount() const { return h_x_.size(); }
    std::size_t getTetCount() const { return h_tetCount_; }

    // Access to device data (for tests)
    DeviceVector<Real>& getDeviceX() { return d_x_; }
    DeviceVector<Real>& getDeviceY() { return d_y_; }
    DeviceVector<Real>& getDeviceZ() { return d_z_; }
    DeviceVector<int>& getDeviceConnectivity() { return d_connectivity_; }

    // Access to cornerstone domain
    DomainType& getDomain() { return *domain_; }

    // Start and end indices for local work assignment
    std::size_t startIndex() const { return domain_->startIndex(); }
    std::size_t endIndex() const { return domain_->endIndex(); }

private:
    int rank_;
    int numRanks_;
    std::unique_ptr<DomainType> domain_;

    // Host data
    std::vector<Real> h_x_, h_y_, h_z_, h_h_;
    std::vector<int> h_connectivity_;
    std::size_t h_tetCount_;

    // Mapping between tetrahedra and nodes
    std::vector<int> h_tetToNodeMap_;
    std::vector<std::vector<int>> h_nodeToTetMap_;

    // Device data - now properly templated
    DeviceVector<Real> d_x_, d_y_, d_z_, d_h_;
    DeviceVector<int> d_connectivity_;
    DeviceVector<int> d_tetToNodeMap_;
};

// CUDA kernel launcher specialization for GPU tag
template <>
template <typename KernelFunc>
void TetrahedralDomain<cstone::GpuTag>::computeOnTetrahedra(KernelFunc kernel) {
    int blockSize = 256;
    int numBlocks = (h_tetCount_ + blockSize - 1) / blockSize;

    kernel<<<numBlocks, blockSize>>>(thrust::raw_pointer_cast(d_x_.data()),
                                     thrust::raw_pointer_cast(d_y_.data()),
                                     thrust::raw_pointer_cast(d_z_.data()),
                                     thrust::raw_pointer_cast(d_connectivity_.data()),
                                     h_tetCount_);
}

// CPU version that doesn't use CUDA syntax
template <typename AcceleratorTag>
template <typename KernelFunc>
void TetrahedralDomain<AcceleratorTag>::computeOnTetrahedra(KernelFunc kernel) {
    // CPU implementation - this would need to be implemented differently
    // Perhaps a serial loop or OpenMP parallelization
    std::cerr << "Warning: computeOnTetrahedra not implemented for CPU mode\n";
}

// Implement TetrahedralDomain methods
template <typename AcceleratorTag>
TetrahedralDomain<AcceleratorTag>::TetrahedralDomain(const std::string& meshFile, int rank, int numRanks)
    : rank_(rank), numRanks_(numRanks) {
    // Read the tetrahedral mesh
    std::vector<double> hostNodes;
    std::vector<int> hostTets;
    std::size_t localNodeCount = 0;
    std::size_t localTetCount = 0;

    readMeshData(meshFile, rank, numRanks, hostNodes, hostTets, localNodeCount, localTetCount);

    // Extract x, y, z coordinates for cornerstone domain
    h_x_.resize(localNodeCount);
    h_y_.resize(localNodeCount);
    h_z_.resize(localNodeCount);

    for (size_t i = 0; i < localNodeCount; ++i) {
        h_x_[i] = hostNodes[i * 3];
        h_y_[i] = hostNodes[i * 3 + 1];
        h_z_[i] = hostNodes[i * 3 + 2];
    }

    // Calculate characteristic sizes
    calculateCharacteristicSizes(hostNodes, hostTets, localNodeCount, localTetCount);

    // Store connectivity
    h_connectivity_ = std::move(hostTets);
    h_tetCount_ = localTetCount;

    // Initialize cornerstone domain with node coordinates
    int bucketSize = 10;
    domain_ = std::make_unique<DomainType>(rank, numRanks, bucketSize);

    // Perform initial sync
    domain_->sync(h_x_, h_y_, h_z_, h_h_);

    // Map tetrahedra to their representative nodes
    mapTetrahedraToNodes();

    // Transfer data to GPU
    transferDataToGPU();

    std::cout << "Rank " << rank_ << " initialized domain with " << localNodeCount << " nodes and " << localTetCount
              << " tetrahedra." << std::endl;
}

// Calculate characteristic sizes for nodes using GPU
template <typename AcceleratorTag>
void TetrahedralDomain<AcceleratorTag>::calculateCharacteristicSizes(const std::vector<double>& hostNodes,
                                                     const std::vector<int>& hostTets,
                                                     std::size_t localNodeCount,
                                                     std::size_t localTetCount) {
    // Allocate device memory
    DeviceVector<double> d_nodes = hostNodes;
    DeviceVector<int> d_tets = hostTets;
    DeviceVector<double> d_h(localNodeCount, 0.0);
    DeviceVector<int> d_nodeTetCount(localNodeCount, 0);

    // Launch kernel to compute characteristic sizes
    int blockSize = 256;
    int numBlocks = (localTetCount + blockSize - 1) / blockSize;

    computeCharacteristicSizesKernel<<<numBlocks, blockSize>>>(
        thrust::raw_pointer_cast(d_nodes.data()),
        thrust::raw_pointer_cast(d_tets.data()),
        thrust::raw_pointer_cast(d_nodeTetCount.data()),
        thrust::raw_pointer_cast(d_h.data()),
        localTetCount);

    // Launch kernel to finalize sizes
    numBlocks = (localNodeCount + blockSize - 1) / blockSize;
    finalizeCharacteristicSizesKernel<<<numBlocks, blockSize>>>(
        thrust::raw_pointer_cast(d_h.data()), 
        thrust::raw_pointer_cast(d_nodeTetCount.data()), 
        localNodeCount);

    // Copy results back to host
    h_h_.resize(localNodeCount);
    thrust::copy(d_h.begin(), d_h.end(), h_h_.begin());
}

// Map tetrahedra to representative nodes using SFC
template <typename AcceleratorTag>
void TetrahedralDomain<AcceleratorTag>::mapTetrahedraToNodes() {
    // Get SFC codes from cornerstone domain
    const auto& sfcCodes = domain_->sfcCodes();
    h_tetToNodeMap_.resize(h_tetCount_);

    // Allocate device vectors
    DeviceVector<int> d_tets = h_connectivity_;
    DeviceVector<unsigned> d_sfcCodes(sfcCodes.begin(), sfcCodes.end());
    DeviceVector<int> d_tetToNodeMap(h_tetCount_);

    // Launch kernel to find representative nodes
    int blockSize = 256;
    int numBlocks = (h_tetCount_ + blockSize - 1) / blockSize;

    findRepresentativeNodesKernel<<<numBlocks, blockSize>>>(
        thrust::raw_pointer_cast(d_tets.data()),
        thrust::raw_pointer_cast(d_sfcCodes.data()),
        thrust::raw_pointer_cast(d_tetToNodeMap.data()),
        h_tetCount_);

    // Copy results back to host
    thrust::copy(d_tetToNodeMap.begin(), d_tetToNodeMap.end(), h_tetToNodeMap_.begin());

    // Build node to tetrahedra mapping
    h_nodeToTetMap_.resize(h_x_.size());
    for (size_t t = 0; t < h_tetCount_; t++) {
        h_nodeToTetMap_[h_tetToNodeMap_[t]].push_back(t);
    }
}

// Transfer data to GPU for computations
template <typename AcceleratorTag>
void TetrahedralDomain<AcceleratorTag>::transferDataToGPU() {
    d_x_ = h_x_;
    d_y_ = h_y_;
    d_z_ = h_z_;
    d_h_ = h_h_;
    d_connectivity_ = h_connectivity_;
    d_tetToNodeMap_ = h_tetToNodeMap_;
}

// Handle domain redistribution
template <typename AcceleratorTag>
void TetrahedralDomain<AcceleratorTag>::sync() {
    // Sync the cornerstone domain
    domain_->sync(h_x_, h_y_, h_z_, h_h_);

    // Re-map tetrahedra to nodes based on updated SFC
    mapTetrahedraToNodes();

    // Update GPU data
    transferDataToGPU();
}

// Get all tetrahedra in a given octree node
template <typename AcceleratorTag>
std::vector<size_t> TetrahedralDomain<AcceleratorTag>::getTetrahedraInOctreeNode(int octreeNodeIndex) {
    // Get nodes in this octree node from cornerstone
    auto nodeIndices = domain_->nodesInOctant(octreeNodeIndex);

    // Collect all tetrahedra represented by these nodes
    std::vector<size_t> result;
    for (int nodeIdx : nodeIndices) {
        for (int tetIdx : h_nodeToTetMap_[nodeIdx]) {
            result.push_back(tetIdx);
        }
    }
    return result;
}