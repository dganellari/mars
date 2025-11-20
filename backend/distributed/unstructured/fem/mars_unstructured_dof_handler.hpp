#pragma once

#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include <mpi.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/gather.h>
#include <thrust/scatter.h>
#include <thrust/partition.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/functional.h>
#include <map>
#include <set>

namespace mars {
namespace fem {

/**
 * @brief Simplified DOF handler for unstructured FEM distributed assembly
 * 
 * Manages ghost/boundary node communication patterns for MPI-parallel FEM.
 * Uses point-to-point MPI communication to exchange DOF values at shared nodes.
 */
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
class UnstructuredDofHandler {
public:
    using Domain = ElementDomain<ElementTag, RealType, KeyType, AcceleratorTag>;
    template<typename T>
    using DeviceVector = typename VectorSelector<T, AcceleratorTag>::type;
    template<typename T>
    using HostVector = std::vector<T>;
    
    UnstructuredDofHandler(Domain& domain, int rank, int numRanks)
        : domain_(domain)
        , rank_(rank)
        , numRanks_(numRanks)
        , initialized_(false)
    {}
    
    /**
     * @brief Initialize DOF handler - identifies boundary and ghost nodes
     * Must be called before scatter_add operations
     */
    void initialize() {
        if (initialized_) return;
        if (numRanks_ == 1) {
            initialized_ = true;
            return; // No communication needed for single rank
        }
        
        // Build boundary/ghost node lists from existing HaloData
        // Note: Halo should already be built by domain initialization
        buildNodeCommunicationPatterns();
        
        initialized_ = true;
    }
    
    /**
     * @brief Scatter-add operation: sum DOF values at shared nodes across ranks
     * 
     * After local assembly, shared nodes at partition boundaries have contributions
     * from multiple ranks. This function:
     * 1. Sends boundary node values to neighboring ranks
     * 2. Receives ghost node contributions from neighbors  
     * 3. Atomically adds received values to local DOFs
     * 
     * Uses GPU-aware MPI to avoid device-host-device copies.
     * 
     * @param dofValues Device vector of DOF values (size = number of nodes)
     */
    template<typename VectorType>
    void scatterAddGhostData(VectorType& dofValues) {
        if (numRanks_ == 1) return; // Nothing to do
        if (!initialized_) {
            throw std::runtime_error("UnstructuredDofHandler::scatterAddGhostData: must call initialize() first");
        }
        
        size_t nodeCount = domain_.getNodeCount();
        if (dofValues.size() != nodeCount) {
            throw std::runtime_error("UnstructuredDofHandler::scatterAddGhostData: vector size mismatch");
        }
        
        // Exchange data directly on device with GPU-aware MPI
        exchangeNodeData(dofValues);
    }
    
    size_t numBoundaryNodes() const { return boundaryNodes_.size(); }
    size_t numGhostNodes() const { return ghostNodes_.size(); }
    bool isInitialized() const { return initialized_; }
    
    /**
     * @brief Exchange node data via GPU-aware MPI point-to-point communication
     * 
     * Pattern:
     * 1. Pack boundary node values into device send buffer
     * 2. MPI_Isend device buffer to all other ranks (GPU-aware MPI)
     * 3. MPI_Irecv into device receive buffer from all other ranks
     * 4. Wait and atomically add received ghost contributions on GPU
     * 
     * Note: Currently broadcasts to all ranks. Optimization would identify
     * actual neighbor ranks from halo element structure.
     */
    template<typename VectorType>
    void exchangeNodeData(VectorType& dofValues) {
        using T = typename VectorType::value_type;
        
        size_t numBoundary = boundaryNodes_.size();
        size_t numGhost = ghostNodes_.size();
        
        // Allocate device buffers for send/recv
        VectorType sendBuf(numBoundary);
        VectorType recvBuf(numGhost);
        
        // Pack boundary values using thrust::gather on device
        DeviceVector<KeyType> d_boundaryNodes(boundaryNodes_.data(), boundaryNodes_.data() + boundaryNodes_.size());
        thrust::gather(thrust::device,
                      d_boundaryNodes.data(), d_boundaryNodes.data() + d_boundaryNodes.size(),
                      thrust::device_pointer_cast(dofValues.data()),
                      sendBuf.data());
        
        // Setup MPI requests for non-blocking send/recv
        std::vector<MPI_Request> requests;
        requests.reserve(2 * (numRanks_ - 1));
        
        // Get raw device pointers for GPU-aware MPI
        T* d_sendPtr = thrust::raw_pointer_cast(sendBuf.data());
        T* d_recvPtr = thrust::raw_pointer_cast(recvBuf.data());
        
        // Post non-blocking receives from all other ranks (GPU-aware)
        for (int r = 0; r < numRanks_; ++r) {
            if (r == rank_) continue;
            
            MPI_Request req;
            MPI_Irecv(d_recvPtr, numGhost, getMPIType<T>(),
                     r, 0, MPI_COMM_WORLD, &req);
            requests.push_back(req);
        }
        
        // Post non-blocking sends to all other ranks (GPU-aware)
        for (int r = 0; r < numRanks_; ++r) {
            if (r == rank_) continue;
            
            MPI_Request req;
            MPI_Isend(d_sendPtr, numBoundary, getMPIType<T>(),
                     r, 0, MPI_COMM_WORLD, &req);
            requests.push_back(req);
        }
        
        // Wait for all communications to complete
        MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
        
        // Add received ghost contributions to owned nodes on device using thrust::scatter
        DeviceVector<KeyType> d_ghostNodes(ghostNodes_.data(), ghostNodes_.data() + ghostNodes_.size());
        
        // Gather current values, add received values, scatter back
        VectorType currentVals(numGhost);
        thrust::gather(thrust::device,
                      d_ghostNodes.data(), d_ghostNodes.data() + d_ghostNodes.size(),
                      thrust::device_pointer_cast(dofValues.data()),
                      currentVals.data());
        
        thrust::transform(thrust::device,
                         currentVals.data(), currentVals.data() + currentVals.size(),
                         recvBuf.data(),
                         currentVals.data(),
                         thrust::plus<T>());
        
        thrust::scatter(thrust::device,
                       currentVals.data(), currentVals.data() + currentVals.size(),
                       d_ghostNodes.data(),
                       thrust::device_pointer_cast(dofValues.data()));
    }
    
    /**
     * @brief Build node communication patterns from element halos
     * 
     * Strategy:
     * 1. Identify owned nodes that appear in halo elements (boundary nodes to send)
     * 2. Identify ghost nodes that appear in halo elements (ghost nodes to receive)
     * 3. Build per-rank send/receive lists
     * 
     * All operations performed on GPU using thrust.
     */
    void buildNodeCommunicationPatterns() {
        const auto& nodeOwnership = domain_.getNodeOwnershipMap();
        const auto& haloElementIndices = domain_.getHaloElementIndices();
        const auto& conn_tuple = domain_.getElementToNodeConnectivity();
        
        size_t nodeCount = domain_.getNodeCount();
        size_t haloCount = haloElementIndices.size();
        
        if (haloCount == 0) {
            if (rank_ == 0) {
                std::cout << "   UnstructuredDofHandler: No halo elements, no communication needed\n";
            }
            return;
        }
        
        // Extract all 4 nodes from each halo element on device
        DeviceVector<KeyType> allHaloNodes(haloCount * 4);
        
        auto conn0_ptr = thrust::device_pointer_cast(std::get<0>(conn_tuple).data());
        auto conn1_ptr = thrust::device_pointer_cast(std::get<1>(conn_tuple).data());
        auto conn2_ptr = thrust::device_pointer_cast(std::get<2>(conn_tuple).data());
        auto conn3_ptr = thrust::device_pointer_cast(std::get<3>(conn_tuple).data());
        auto halo_ptr = thrust::device_pointer_cast(haloElementIndices.data());
        auto allHalo_ptr = thrust::raw_pointer_cast(allHaloNodes.data());
        
        // Gather nodes from halo elements
        thrust::for_each(thrust::device,
                        thrust::make_counting_iterator<size_t>(0),
                        thrust::make_counting_iterator<size_t>(haloCount),
                        [=] __device__ (size_t i) {
                            KeyType elemIdx = halo_ptr[i];
                            allHalo_ptr[i * 4 + 0] = conn0_ptr[elemIdx];
                            allHalo_ptr[i * 4 + 1] = conn1_ptr[elemIdx];
                            allHalo_ptr[i * 4 + 2] = conn2_ptr[elemIdx];
                            allHalo_ptr[i * 4 + 3] = conn3_ptr[elemIdx];
                        });
        
        // Sort and unique to get unique nodes
        thrust::sort(thrust::device, allHaloNodes.data(), allHaloNodes.data() + allHaloNodes.size());
        auto new_end = thrust::unique(thrust::device, allHaloNodes.data(), allHaloNodes.data() + allHaloNodes.size());
        size_t new_size = new_end - allHaloNodes.data();
        allHaloNodes.resize(new_size);
        
        // Partition into boundary (owned) and ghost (not owned) nodes
        auto ownership_ptr = thrust::device_pointer_cast(nodeOwnership.data());
        
        auto partition_point = thrust::partition(thrust::device,
                                                 allHaloNodes.data(),
                                                 allHaloNodes.data() + allHaloNodes.size(),
                                                 [=] __device__ (KeyType nodeIdx) {
                                                     return ownership_ptr[nodeIdx] == 1;
                                                 });
        
        // Copy boundary and ghost nodes to host vectors
        size_t numBoundary = partition_point - allHaloNodes.data();
        size_t numGhost = allHaloNodes.size() - numBoundary;
        
        boundaryNodes_.resize(numBoundary);
        ghostNodes_.resize(numGhost);
        
        thrust::copy(thrust::device, allHaloNodes.data(), allHaloNodes.data() + numBoundary, boundaryNodes_.begin());
        thrust::copy(thrust::device, allHaloNodes.data() + numBoundary, allHaloNodes.data() + allHaloNodes.size(), ghostNodes_.begin());
        
        if (rank_ == 0) {
            std::cout << "   UnstructuredDofHandler: " << boundaryNodes_.size()
                     << " boundary nodes, " << ghostNodes_.size() << " ghost nodes\n";
        }
    }
    
private:
    // Helper to get MPI datatype
    template<typename T>
    MPI_Datatype getMPIType() const {
        if constexpr (std::is_same_v<T, float>) return MPI_FLOAT;
        else if constexpr (std::is_same_v<T, double>) return MPI_DOUBLE;
        else if constexpr (std::is_same_v<T, int>) return MPI_INT;
        else static_assert(std::is_same_v<T, float>, "Unsupported MPI type");
        return MPI_FLOAT;
    }
    
    Domain& domain_;
    int rank_;
    int numRanks_;
    bool initialized_;
    
    // Communication patterns
    HostVector<KeyType> boundaryNodes_;  // Owned nodes to send
    HostVector<KeyType> ghostNodes_;      // Ghost nodes to receive
};

} // namespace fem
} // namespace mars
