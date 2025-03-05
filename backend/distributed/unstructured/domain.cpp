#include "cstone/domain/domain.hpp"
#include <adios2.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/extrema.h>
#include <thrust/host_vector.h>
#include <algorithm>
#include <array>
#include <iostream>
#include <vector>
#include "mars.hpp"
#include "mars_globals.hpp"
#include "mars_read_mesh_adios2.hpp"

// CUDA kernel to compute characteristic size (h) for each node
__global__ void computeCharacteristicSizesKernel(double* nodes,
                                                 int* tets,
                                                 int* nodeTetCount,
                                                 double* h,
                                                 int localTetCount) {
    int tetIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (tetIdx < localTetCount) {
        // For each tetrahedron, process all 6 edges
        for (int i = 0; i < 4; ++i) {
            int nodeIdx1 = tets[tetIdx * 4 + i];

            for (int j = i + 1; j < 4; ++j) {
                int nodeIdx2 = tets[tetIdx * 4 + j];

                // Calculate edge length
                double dx = nodes[nodeIdx1 * 3] - nodes[nodeIdx2 * 3];
                double dy = nodes[nodeIdx1 * 3 + 1] - nodes[nodeIdx2 * 3 + 1];
                double dz = nodes[nodeIdx1 * 3 + 2] - nodes[nodeIdx2 * 3 + 2];
                double edgeLength = sqrt(dx * dx + dy * dy + dz * dz);

                // Add to each node's total (using atomic to avoid race conditions)
                atomicAdd(&h[nodeIdx1], edgeLength);
                atomicAdd(&h[nodeIdx2], edgeLength);
                atomicAdd(&nodeTetCount[nodeIdx1], 1);
                atomicAdd(&nodeTetCount[nodeIdx2], 1);
            }
        }
    }
}

// CUDA kernel to finalize characteristic sizes
__global__ void finalizeCharacteristicSizesKernel(double* h, int* nodeTetCount, int numNodes) {
    int nodeIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (nodeIdx < numNodes) {
        if (nodeTetCount[nodeIdx] > 0) {
            h[nodeIdx] /= nodeTetCount[nodeIdx];
        } else {
            h[nodeIdx] = 0.01;  // Default for isolated nodes
        }
    }
}

// CUDA kernel to find the representative node for each tetrahedron
__global__ void findRepresentativeNodesKernel(int* tets, unsigned* sfcCodes, int* tetToNodeMap, int numTets) {
    int tetIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (tetIdx < numTets) {
        // Get the four nodes of this tetrahedron
        int node0 = tets[tetIdx * 4];
        int node1 = tets[tetIdx * 4 + 1];
        int node2 = tets[tetIdx * 4 + 2];
        int node3 = tets[tetIdx * 4 + 3];

        // Get SFC codes
        unsigned sfc0 = sfcCodes[node0];
        unsigned sfc1 = sfcCodes[node1];
        unsigned sfc2 = sfcCodes[node2];
        unsigned sfc3 = sfcCodes[node3];

        // Find minimum SFC code
        int repNode = node0;
        unsigned minSfc = sfc0;

        if (sfc1 < minSfc) {
            minSfc = sfc1;
            repNode = node1;
        }

        if (sfc2 < minSfc) {
            minSfc = sfc2;
            repNode = node2;
        }

        if (sfc3 < minSfc) {
            minSfc = sfc3;
            repNode = node3;
        }

        // Store the representative node
        tetToNodeMap[tetIdx] = repNode;
    }
}

// TetrahedralDomain class with CUDA support
using Real = double;
using KeyType = unsigned;

class TetrahedralDomain {
public:
    TetrahedralDomain(const std::string& meshFile, int rank, int numRanks)
        : rank_(rank), numRanks_(numRanks), domain_(rank, numRanks, 10) {
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

        // Initialize cornerstone domain with node coordinates
        domain_.sync(h_x_, h_y_, h_z_, h_h_);

        // Store connectivity
        h_connectivity_ = std::move(hostTets);
        h_tetCount_ = localTetCount;

        // Map tetrahedra to their representative nodes
        mapTetrahedraToNodes();

        // Transfer data to GPU
        transferDataToGPU();

        std::cout << "Rank " << rank_ << " initialized domain with " << localNodeCount << " nodes and " << localTetCount
                  << " tetrahedra." << std::endl;
    }

    // Calculate characteristic sizes for nodes using GPU
    void calculateCharacteristicSizes(const std::vector<double>& hostNodes,
                                      const std::vector<int>& hostTets,
                                      std::size_t localNodeCount,
                                      std::size_t localTetCount) {
        // Allocate device memory
        thrust::device_vector<double> d_nodes = hostNodes;
        thrust::device_vector<int> d_tets = hostTets;
        thrust::device_vector<double> d_h(localNodeCount, 0.0);
        thrust::device_vector<int> d_nodeTetCount(localNodeCount, 0);

        // Launch kernel to compute characteristic sizes
        int blockSize = 256;
        int numBlocks = (localTetCount + blockSize - 1) / blockSize;

        computeCharacteristicSizesKernel<<<numBlocks, blockSize>>>(thrust::raw_pointer_cast(d_nodes.data()),
                                                                   thrust::raw_pointer_cast(d_tets.data()),
                                                                   thrust::raw_pointer_cast(d_nodeTetCount.data()),
                                                                   thrust::raw_pointer_cast(d_h.data()),
                                                                   localTetCount);

        // Launch kernel to finalize sizes
        numBlocks = (localNodeCount + blockSize - 1) / blockSize;
        finalizeCharacteristicSizesKernel<<<numBlocks, blockSize>>>(
            thrust::raw_pointer_cast(d_h.data()), thrust::raw_pointer_cast(d_nodeTetCount.data()), localNodeCount);

        // Copy results back to host
        h_h_.resize(localNodeCount);
        thrust::copy(d_h.begin(), d_h.end(), h_h_.begin());
    }

    // Map tetrahedra to representative nodes using SFC
    void mapTetrahedraToNodes() {
        // Get SFC codes from cornerstone domain
        const auto& sfcCodes = domain_.sfcCodes();
        h_tetToNodeMap_.resize(h_tetCount_);

        // Allocate device vectors
        thrust::device_vector<int> d_tets = h_connectivity_;
        thrust::device_vector<unsigned> d_sfcCodes(sfcCodes.begin(), sfcCodes.end());
        thrust::device_vector<int> d_tetToNodeMap(h_tetCount_);

        // Launch kernel to find representative nodes
        int blockSize = 256;
        int numBlocks = (h_tetCount_ + blockSize - 1) / blockSize;

        findRepresentativeNodesKernel<<<numBlocks, blockSize>>>(thrust::raw_pointer_cast(d_tets.data()),
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
    void transferDataToGPU() {
        d_x_ = h_x_;
        d_y_ = h_y_;
        d_z_ = h_z_;
        d_h_ = h_h_;
        d_connectivity_ = h_connectivity_;
        d_tetToNodeMap_ = h_tetToNodeMap_;
    }

    // Handle domain redistribution
    void sync() {
        // Sync the cornerstone domain
        domain_.sync(h_x_, h_y_, h_z_, h_h_);

        // Re-map tetrahedra to nodes based on updated SFC
        mapTetrahedraToNodes();

        // Update GPU data
        transferDataToGPU();
    }

    // Get all tetrahedra in a given octree node
    std::vector<size_t> getTetrahedraInOctreeNode(int octreeNodeIndex) {
        // Get nodes in this octree node from cornerstone
        auto nodeIndices = domain_.nodesInOctant(octreeNodeIndex);

        // Collect all tetrahedra represented by these nodes
        std::vector<size_t> result;
        for (int nodeIdx : nodeIndices) {
            for (int tetIdx : h_nodeToTetMap_[nodeIdx]) {
                result.push_back(tetIdx);
            }
        }
        return result;
    }

    // Launch a CUDA kernel to compute something on tetrahedra
    template <typename KernelFunc>
    void computeOnTetrahedra(KernelFunc kernel) {
        int blockSize = 256;
        int numBlocks = (h_tetCount_ + blockSize - 1) / blockSize;

        kernel<<<numBlocks, blockSize>>>(thrust::raw_pointer_cast(d_x_.data()),
                                         thrust::raw_pointer_cast(d_y_.data()),
                                         thrust::raw_pointer_cast(d_z_.data()),
                                         thrust::raw_pointer_cast(d_connectivity_.data()),
                                         h_tetCount_);
    }

    // Access methods for host data
    const std::vector<Real>& getX() const { return h_x_; }
    const std::vector<Real>& getY() const { return h_y_; }
    const std::vector<Real>& getZ() const { return h_z_; }
    std::size_t getNodeCount() const { return h_x_.size(); }
    std::size_t getTetCount() const { return h_tetCount_; }

    // Access to cornerstone domain
    cstone::Domain<KeyType, Real>& getDomain() { return domain_; }

    // Start and end indices for local work assignment
    std::size_t startIndex() const { return domain_.startIndex(); }
    std::size_t endIndex() const { return domain_.endIndex(); }

private:
    int rank_;
    int numRanks_;

    // Cornerstone domain
    cstone::Domain<KeyType, Real> domain_;

    // Host data
    std::vector<Real> h_x_, h_y_, h_z_, h_h_;
    std::vector<int> h_connectivity_;
    std::size_t h_tetCount_;

    // Mapping between tetrahedra and nodes
    std::vector<int> h_tetToNodeMap_;
    std::vector<std::vector<int>> h_nodeToTetMap_;

    // Device data
    thrust::device_vector<Real> d_x_, d_y_, d_z_, d_h_;
    thrust::device_vector<int> d_connectivity_;
    thrust::device_vector<int> d_tetToNodeMap_;
};

// Example kernel to compute something on tetrahedra (e.g., volume)
__global__ void computeTetrahedralVolumesKernel(double* x, double* y, double* z, int* connectivity, int numTets) {
    int tetIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (tetIdx < numTets) {
        // Get the four nodes of this tetrahedron
        int n0 = connectivity[tetIdx * 4];
        int n1 = connectivity[tetIdx * 4 + 1];
        int n2 = connectivity[tetIdx * 4 + 2];
        int n3 = connectivity[tetIdx * 4 + 3];

        // Get node coordinates
        double x0 = x[n0], y0 = y[n0], z0 = z[n0];
        double x1 = x[n1], y1 = y[n1], z1 = z[n1];
        double x2 = x[n2], y2 = y[n2], z2 = z[n2];
        double x3 = x[n3], y3 = y[n3], z3 = z[n3];

        // Compute vectors for volume calculation
        double v1x = x1 - x0, v1y = y1 - y0, v1z = z1 - z0;
        double v2x = x2 - x0, v2y = y2 - y0, v2z = z2 - z0;
        double v3x = x3 - x0, v3y = y3 - y0, v3z = z3 - z0;

        // Compute volume using the scalar triple product
        double volume =
            fabs(v1x * (v2y * v3z - v2z * v3y) + v1y * (v2z * v3x - v2x * v3z) + v1z * (v2x * v3y - v2y * v3x)) / 6.0;

        // Volume calculation result could be stored or used here
    }
}

int main(int argc, char** argv) {
    int rank = 0, numRanks = 1;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    try {
        // Create the tetrahedral domain
        TetrahedralDomain domain("tetrahedral_mesh.bp", rank, numRanks);

        // Compute volume of all tetrahedra
        domain.computeOnTetrahedra(computeTetrahedralVolumesKernel);

        // Example iteration over octree nodes
        for (int octreeNode = 0; octreeNode < domain.getDomain().numOctants(); octreeNode++) {
            // Get tetrahedra in this octree node
            auto tetsInNode = domain.getTetrahedraInOctreeNode(octreeNode);

            if (rank == 0) {
                std::cout << "Octree node " << octreeNode << " contains " << tetsInNode.size() << " tetrahedra\n";
            }
        }

        // Sync to handle any changes
        domain.sync();

    } catch (const std::exception& e) {
        std::cerr << "Error on rank " << rank << ": " << e.what() << std::endl;
    }

    // Finalize MPI
    MPI_Finalize();
    return 0;
}
}
