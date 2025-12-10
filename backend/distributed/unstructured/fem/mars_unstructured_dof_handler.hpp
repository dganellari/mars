#pragma once

#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include <mpi.h>
#include <unordered_map>
#ifdef MARS_ENABLE_HYPRE
#include <HYPRE.h>
#endif
#include <functional>
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
#include <limits>
#include <algorithm>

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
        , numLocalDofs_(0)
        , numLocalGhostDofs_(0)
        , numGlobalDofs_(0)
        , ownedDofStart_(0)
        , ownedDofEnd_(0)
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
     * @brief Enumerate global DOFs for distributed assembly
     * Creates consistent global DOF numbering across ranks
     */
    void enumerate_dofs() {
        // Standard parallel FEM: Node owned by rank that owns a LOCAL element containing it
        // Use Cornerstone's element partition directly

        const auto& nodeSfcKeys = domain_.getLocalToGlobalSfcMap();
        size_t numTotalNodes = domain_.getNodeCount();
        size_t localElementCount = domain_.localElementCount();

        if (numTotalNodes == 0) {
            std::cerr << "Rank " << rank_ << ": Error - domain has 0 nodes\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Copy SFC keys to host
        thrust::host_vector<KeyType> h_sfcKeys(numTotalNodes);
        thrust::copy(thrust::device_pointer_cast(nodeSfcKeys.data()),
                    thrust::device_pointer_cast(nodeSfcKeys.data() + numTotalNodes),
                    h_sfcKeys.begin());

        // Get element connectivity (only LOCAL elements [0, localElementCount))
        const auto& conn_tuple = domain_.getElementToNodeConnectivity();
        const auto& conn0 = std::get<0>(conn_tuple);
        const auto& conn1 = std::get<1>(conn_tuple);
        const auto& conn2 = std::get<2>(conn_tuple);
        const auto& conn3 = std::get<3>(conn_tuple);

        // Copy LOCAL element connectivity to host
        thrust::host_vector<KeyType> h_conn0(localElementCount);
        thrust::host_vector<KeyType> h_conn1(localElementCount);
        thrust::host_vector<KeyType> h_conn2(localElementCount);
        thrust::host_vector<KeyType> h_conn3(localElementCount);
        thrust::copy_n(thrust::device_pointer_cast(conn0.data()), localElementCount, h_conn0.begin());
        thrust::copy_n(thrust::device_pointer_cast(conn1.data()), localElementCount, h_conn1.begin());
        thrust::copy_n(thrust::device_pointer_cast(conn2.data()), localElementCount, h_conn2.begin());
        thrust::copy_n(thrust::device_pointer_cast(conn3.data()), localElementCount, h_conn3.begin());

        // Find node SFC keys that appear in my LOCAL elements
        // Connectivity uses local node indices, so convert to SFC keys
        std::set<KeyType> myLocalNodeSfcs;
        for (size_t e = 0; e < localElementCount; ++e) {
            // Each conn value is a local node index
            KeyType idx0 = h_conn0[e];
            KeyType idx1 = h_conn1[e];
            KeyType idx2 = h_conn2[e];
            KeyType idx3 = h_conn3[e];

            // Convert to SFC keys
            if (idx0 < numTotalNodes) myLocalNodeSfcs.insert(h_sfcKeys[idx0]);
            if (idx1 < numTotalNodes) myLocalNodeSfcs.insert(h_sfcKeys[idx1]);
            if (idx2 < numTotalNodes) myLocalNodeSfcs.insert(h_sfcKeys[idx2]);
            if (idx3 < numTotalNodes) myLocalNodeSfcs.insert(h_sfcKeys[idx3]);
        }

        std::vector<KeyType> myOwnedNodeSfcs(myLocalNodeSfcs.begin(), myLocalNodeSfcs.end());

        // MPI: Gather all ranks' owned node SFC keys
        std::vector<int> allCounts(numRanks_);
        int myCount = static_cast<int>(myOwnedNodeSfcs.size());
        MPI_Allgather(&myCount, 1, MPI_INT, allCounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

        std::vector<int> displacements(numRanks_);
        size_t totalCandidates = 0;
        for (int r = 0; r < numRanks_; ++r) {
            displacements[r] = static_cast<int>(totalCandidates);
            totalCandidates += allCounts[r];
        }

        std::vector<KeyType> allOwnedNodeSfcs(totalCandidates);
        MPI_Datatype mpiKeyType = (sizeof(KeyType) == 8) ? MPI_UNSIGNED_LONG : MPI_UNSIGNED;
        MPI_Allgatherv(myOwnedNodeSfcs.data(), myCount, mpiKeyType,
                      allOwnedNodeSfcs.data(), allCounts.data(), displacements.data(),
                      mpiKeyType, MPI_COMM_WORLD);

        // Build map: node SFC → owner rank (balanced tie-breaking)
        // First, collect all ranks that have each node in local elements
        std::unordered_map<KeyType, std::vector<int>> nodeSfcToCandidates;
        for (int r = 0; r < numRanks_; ++r) {
            int start = displacements[r];
            int end = start + allCounts[r];
            for (int i = start; i < end; ++i) {
                KeyType sfc = allOwnedNodeSfcs[i];
                nodeSfcToCandidates[sfc].push_back(r);
            }
        }

        // Resolve ownership: for shared nodes, use SFC-based hash for balanced distribution
        std::unordered_map<KeyType, int> nodeSfcToOwner;
        for (const auto& [sfc, candidates] : nodeSfcToCandidates) {
            if (candidates.size() == 1) {
                // Unshared node - owned by the only candidate
                nodeSfcToOwner[sfc] = candidates[0];
            } else {
                // Shared node - use SFC hash modulo number of candidates
                size_t winner_idx = sfc % candidates.size();
                nodeSfcToOwner[sfc] = candidates[winner_idx];
            }
        }

        // Assign resolved ownership based on this mapping
        resolvedOwnership_.resize(numTotalNodes);
        size_t numTrueOwned = 0;
        for (size_t i = 0; i < numTotalNodes; ++i) {
            KeyType sfc = h_sfcKeys[i];
            auto it = nodeSfcToOwner.find(sfc);
            if (it != nodeSfcToOwner.end() && it->second == rank_) {
                resolvedOwnership_[i] = 1;
                numTrueOwned++;
            } else {
                resolvedOwnership_[i] = 0;
            }
        }

        // Verify ownership totals via MPI reduction
        size_t globalOwnedCount = 0;
        MPI_Reduce(&numTrueOwned, &globalOwnedCount, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

        // Also check total local nodes across all ranks
        size_t totalLocalNodes = 0;
        size_t myLocalNodes = numTotalNodes;
        MPI_Reduce(&myLocalNodes, &totalLocalNodes, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

        if (rank_ == 0) {
            std::cout << "\n=== Ownership Resolution Summary ===\n";
            std::cout << "Total unique nodes in domain: " << nodeSfcToOwner.size() << "\n";
            std::cout << "Total owned nodes across all ranks: " << globalOwnedCount << "\n";
            std::cout << "Total local nodes (with overlaps): " << totalLocalNodes << "\n";
            std::cout << "Average nodes per rank: " << (totalLocalNodes / numRanks_) << "\n";
            std::cout << "Missing nodes (if < original mesh): " << (5265 - nodeSfcToOwner.size()) << " nodes\n";
            if (globalOwnedCount != nodeSfcToOwner.size()) {
                std::cout << "ERROR: Ownership mismatch! Expected " << nodeSfcToOwner.size()
                          << " but got " << globalOwnedCount << "\n";
            }
        }

        std::cout << "Rank " << rank_ << ": After MPI consensus - "
                  << numTrueOwned << " owned nodes out of " << numTotalNodes
                  << " total nodes (from " << myOwnedNodeSfcs.size() << " candidates in local elements)\n";

        // Build unique global DOF list (sorted SFC keys)
        std::vector<KeyType> uniqueSfcKeys;
        uniqueSfcKeys.reserve(nodeSfcToOwner.size());
        for (const auto& [sfc, rank] : nodeSfcToOwner) {
            uniqueSfcKeys.push_back(sfc);
        }
        std::sort(uniqueSfcKeys.begin(), uniqueSfcKeys.end());

        numGlobalDofs_ = uniqueSfcKeys.size();

        // Build SFC-to-global-DOF map
        std::unordered_map<KeyType, KeyType> sfcToGlobalDof;
        for (size_t i = 0; i < uniqueSfcKeys.size(); ++i) {
            sfcToGlobalDof[uniqueSfcKeys[i]] = static_cast<KeyType>(i);
        }

        // Build localToGlobalDof_ for truly owned nodes (in SFC order)
        std::vector<std::pair<KeyType, size_t>> ownedNodesSfc;
        for (size_t i = 0; i < numTotalNodes; ++i) {
            if (resolvedOwnership_[i] == 1) {
                ownedNodesSfc.push_back({h_sfcKeys[i], i});
            }
        }
        std::sort(ownedNodesSfc.begin(), ownedNodesSfc.end());

        localToGlobalDof_.resize(numTrueOwned);
        for (size_t i = 0; i < ownedNodesSfc.size(); ++i) {
            localToGlobalDof_[i] = sfcToGlobalDof[ownedNodesSfc[i].first];
        }

        // Store owned DOF range
        if (numTrueOwned > 0) {
            ownedDofStart_ = localToGlobalDof_[0];
            ownedDofEnd_ = localToGlobalDof_[numTrueOwned - 1] + 1;
        } else {
            ownedDofStart_ = 0;
            ownedDofEnd_ = 0;
        }

        if (rank_ == 0) {
            std::cout << "Global DOFs: " << numGlobalDofs_ << " (resolved from "
                      << totalCandidates << " candidates)\n";
        }

        // Store maps for use in buildNodeToDofMappings
        sfcToGlobalDof_ = std::move(sfcToGlobalDof);

        // Build node-to-DOF mappings for assembly
        buildNodeToDofMappings();
    }
    
    /**
     * @brief Build node-to-DOF mappings for assembly
     * Creates mappings from node indices to local and global DOF indices
     * Uses resolved ownership (after MPI consensus) for non-overlapping DOF assignment
     */
    void buildNodeToDofMappings() {
        size_t numTotalNodes = domain_.getNodeCount();

        // Get SFC keys from domain
        const auto& nodeSfcKeys = domain_.getLocalToGlobalSfcMap();
        thrust::host_vector<KeyType> h_sfcKeys(numTotalNodes);
        thrust::copy(thrust::device_pointer_cast(nodeSfcKeys.data()),
                    thrust::device_pointer_cast(nodeSfcKeys.data() + numTotalNodes),
                    h_sfcKeys.begin());

        // Initialize mappings
        nodeToLocalDof_.assign(numTotalNodes, static_cast<KeyType>(-1));
        nodeToGlobalDof_.assign(numTotalNodes, static_cast<KeyType>(-1));

        // PASS 1: Build mappings for truly owned nodes (using resolvedOwnership_)
        std::vector<std::pair<KeyType, size_t>> ownedNodesSfc;  // (sfcKey, nodeIdx)
        for (size_t nodeIdx = 0; nodeIdx < numTotalNodes; ++nodeIdx) {
            if (resolvedOwnership_[nodeIdx] == 1) {
                ownedNodesSfc.push_back({h_sfcKeys[nodeIdx], nodeIdx});
            }
        }

        // Sort owned nodes by SFC key for spatial locality
        std::sort(ownedNodesSfc.begin(), ownedNodesSfc.end());

        // Assign DOF indices
        size_t localDofCounter = 0;
        for (const auto& [sfcKey, nodeIdx] : ownedNodesSfc) {
            nodeToLocalDof_[nodeIdx] = localDofCounter;
            auto it = sfcToGlobalDof_.find(sfcKey);
            if (it != sfcToGlobalDof_.end()) {
                nodeToGlobalDof_[nodeIdx] = it->second;
            } else {
                nodeToGlobalDof_[nodeIdx] = static_cast<KeyType>(-1);
            }
            localDofCounter++;
        }

        numLocalDofs_ = localDofCounter;

        // PASS 2: Assign local indices to ghost nodes (after owned DOFs)
        // Ghost = resolvedOwnership_ == 0 (includes nodes lost in tie-breaking)
        size_t ghostDofCounter = 0;
        ghostDofToNode_.clear();

        // Count ghosts
        size_t numGhosts = 0;
        for (size_t nodeIdx = 0; nodeIdx < numTotalNodes; ++nodeIdx) {
            if (resolvedOwnership_[nodeIdx] == 0) numGhosts++;
        }
        ghostDofToNode_.resize(numGhosts);

        for (size_t nodeIdx = 0; nodeIdx < numTotalNodes; ++nodeIdx) {
            if (resolvedOwnership_[nodeIdx] == 0) {  // Ghost
                nodeToLocalDof_[nodeIdx] = numLocalDofs_ + ghostDofCounter;
                // Ghost nodes also get global DOF from SFC map
                KeyType sfcKey = h_sfcKeys[nodeIdx];
                auto it = sfcToGlobalDof_.find(sfcKey);
                if (it != sfcToGlobalDof_.end()) {
                    nodeToGlobalDof_[nodeIdx] = it->second;
                } else {
                    nodeToGlobalDof_[nodeIdx] = static_cast<KeyType>(-1);
                }
                ghostDofToNode_[ghostDofCounter] = nodeIdx;
                ghostDofCounter++;
            }
        }

        numLocalGhostDofs_ = ghostDofCounter;

        if (rank_ == 0) {
            std::cout << "Built DOF mappings: " << numLocalDofs_ << " owned + "
                      << numLocalGhostDofs_ << " ghost DOFs (rank 0)\n";
        }
        
        // Exchange ghost DOF global indices with owning ranks
        if (numRanks_ > 1 && numLocalGhostDofs_ > 0) {
            exchangeGhostDofIndices();
        }
    }
    
    /**
     * @brief Exchange global DOF indices for ghost nodes with owning ranks
     * Uses node SFC keys to identify which rank owns each ghost node
     */
    void exchangeGhostDofIndices() {
        size_t numTotalNodes = domain_.getNodeCount();
        
        // Get node ownership and SFC keys
        const auto& nodeOwnership = domain_.getNodeOwnershipMap();
        const auto& nodeSfcKeys = domain_.getLocalToGlobalSfcMap();
        
        thrust::host_vector<uint8_t> h_ownership(numTotalNodes);
        thrust::host_vector<KeyType> h_sfcKeys(numTotalNodes);
        
        thrust::copy(thrust::device_pointer_cast(nodeOwnership.data()),
                    thrust::device_pointer_cast(nodeOwnership.data() + numTotalNodes),
                    h_ownership.begin());
        thrust::copy(thrust::device_pointer_cast(nodeSfcKeys.data()),
                    thrust::device_pointer_cast(nodeSfcKeys.data() + numTotalNodes),
                    h_sfcKeys.begin());
        
        // Build mapping from owned node SFC keys to global DOF indices for fast lookup
        std::unordered_map<KeyType, KeyType> ownedSfcToGlobalDof;
        for (size_t nodeIdx = 0; nodeIdx < numTotalNodes; ++nodeIdx) {
            if (h_ownership[nodeIdx] == 1 || h_ownership[nodeIdx] == 2) {  // Owned (1) or shared (2)
                ownedSfcToGlobalDof[h_sfcKeys[nodeIdx]] = nodeToGlobalDof_[nodeIdx];
            }
        }
        
        // Gather rank boundaries from all ranks to enable findRank locally
        std::vector<KeyType> rankBoundaries(numRanks_ + 1);
        KeyType myAssignmentStart = domain_.getDomain().assignmentStart();
        
        // Determine MPI datatype based on KeyType size
        MPI_Datatype mpiKeyType = (sizeof(KeyType) == 8) ? MPI_UNSIGNED_LONG : MPI_UNSIGNED;
        
        // Gather assignment start keys from all ranks
        MPI_Allgather(&myAssignmentStart, 1, mpiKeyType,
                     rankBoundaries.data(), 1, mpiKeyType,
                     MPI_COMM_WORLD);
        
        // Set the final boundary (maximum possible SFC key)
        rankBoundaries[numRanks_] = std::numeric_limits<KeyType>::max();
        
        // Lambda to find which rank owns a given SFC key (binary search)
        auto findRank = [&](KeyType key) -> int {
            auto it = std::upper_bound(rankBoundaries.begin(), rankBoundaries.end(), key);
            return static_cast<int>(it - rankBoundaries.begin()) - 1;
        };
        
        // Build communication pattern: which ghost nodes need data from which ranks?
        std::map<int, std::vector<size_t>> ghostNodesByRank;  // rank -> list of local ghost node indices
        std::map<int, std::vector<KeyType>> sfcKeysByRank;    // rank -> list of SFC keys to query
        
        for (size_t nodeIdx = 0; nodeIdx < numTotalNodes; ++nodeIdx) {
            if (h_ownership[nodeIdx] == 0) {  // Pure ghost node (state 0, not shared)
                KeyType sfcKey = h_sfcKeys[nodeIdx];
                
                // Find which rank owns this SFC key
                int ownerRank = findRank(sfcKey);
                
                if (ownerRank >= 0 && ownerRank < numRanks_ && ownerRank != rank_) {
                    ghostNodesByRank[ownerRank].push_back(nodeIdx);
                    sfcKeysByRank[ownerRank].push_back(sfcKey);
                }
            }
        }
        
        // Count neighbors for MPI exchange
        std::set<int> neighborRanks;
        for (auto& [rank, _] : ghostNodesByRank) {
            neighborRanks.insert(rank);
        }
        
        // Perform non-blocking MPI exchange
        std::vector<MPI_Request> sendRequests;
        std::vector<MPI_Request> recvRequests;
        std::map<int, std::vector<KeyType>> sendBuffers;
        
        // Ghost DOF exchange with neighbors (silent)
        
        // Use MPI_Alltoall to exchange counts reliably
        std::vector<int> sendCountsVec(numRanks_, 0);
        std::vector<int> recvCountsVec(numRanks_, 0);
        
        for (int neighborRank : neighborRanks) {
            sendCountsVec[neighborRank] = sfcKeysByRank[neighborRank].size();
        }
        
        MPI_Alltoall(sendCountsVec.data(), 1, MPI_INT,
                    recvCountsVec.data(), 1, MPI_INT,
                    MPI_COMM_WORLD);
        
        // Prepare receive buffers based on counts
        std::map<int, std::vector<KeyType>> recvBuffers;
        for (int r = 0; r < numRanks_; ++r) {
            if (r != rank_ && recvCountsVec[r] > 0) {
                recvBuffers[r].resize(recvCountsVec[r]);
            }
        }
        
        MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE);
        sendRequests.clear();
        
        
        // Now exchange SFC keys
        for (auto& [neighborRank, sfcKeys] : sfcKeysByRank) {
            sendBuffers[neighborRank] = sfcKeys;
            MPI_Request req;
            MPI_Isend(sendBuffers[neighborRank].data(), sendBuffers[neighborRank].size(), 
                     mpiKeyType, neighborRank, /*tag=*/101, MPI_COMM_WORLD, &req);
            sendRequests.push_back(req);
        }
        
        for (auto& [neighborRank, buffer] : recvBuffers) {
            MPI_Request req;
            MPI_Irecv(buffer.data(), buffer.size(), mpiKeyType, neighborRank, 
                     /*tag=*/101, MPI_COMM_WORLD, &req);
            recvRequests.push_back(req);
        }
        
        MPI_Waitall(recvRequests.size(), recvRequests.data(), MPI_STATUSES_IGNORE);
        recvRequests.clear();
        
        
        // Lookup global DOF indices for received SFC keys and send back
        std::map<int, std::vector<KeyType>> globalDofResponses;
        for (auto& [neighborRank, queriedSfcKeys] : recvBuffers) {
            std::vector<KeyType> responseDofs(queriedSfcKeys.size());
            for (size_t i = 0; i < queriedSfcKeys.size(); ++i) {
                auto it = ownedSfcToGlobalDof.find(queriedSfcKeys[i]);
                responseDofs[i] = (it != ownedSfcToGlobalDof.end()) ? it->second : static_cast<KeyType>(-1);
            }
            globalDofResponses[neighborRank] = responseDofs;
            
            MPI_Request req;
            MPI_Isend(globalDofResponses[neighborRank].data(), globalDofResponses[neighborRank].size(),
                     mpiKeyType, neighborRank, /*tag=*/102, MPI_COMM_WORLD, &req);
            sendRequests.push_back(req);
        }
        
        // Receive global DOF indices back
        std::map<int, std::vector<KeyType>> receivedGlobalDofs;
        for (auto& [neighborRank, ghostNodeList] : ghostNodesByRank) {
            receivedGlobalDofs[neighborRank].resize(ghostNodeList.size());
            MPI_Request req;
            MPI_Irecv(receivedGlobalDofs[neighborRank].data(), receivedGlobalDofs[neighborRank].size(),
                     mpiKeyType, neighborRank, /*tag=*/102, MPI_COMM_WORLD, &req);
            recvRequests.push_back(req);
        }
        
        MPI_Waitall(recvRequests.size(), recvRequests.data(), MPI_STATUSES_IGNORE);
        MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE);
        
        // Update global DOF indices for ghost nodes
        for (auto& [neighborRank, ghostNodeList] : ghostNodesByRank) {
            const auto& dofIndices = receivedGlobalDofs[neighborRank];
            for (size_t i = 0; i < ghostNodeList.size(); ++i) {
                size_t nodeIdx = ghostNodeList[i];
                nodeToGlobalDof_[nodeIdx] = dofIndices[i];
            }
        }
        
        if (rank_ == 0 || rank_ == 1) {
            size_t numResolvedGhosts = 0;
            for (size_t nodeIdx = 0; nodeIdx < numTotalNodes; ++nodeIdx) {
                if (h_ownership[nodeIdx] == 0 && nodeToGlobalDof_[nodeIdx] != static_cast<KeyType>(-1)) {
                    numResolvedGhosts++;
                }
            }
            // Ghost DOF resolution complete (silent)
        }
        
        // Build neighbor communication lists for halo exchange
        buildNeighborCommunication();
    }
    
    /**
     * @brief Build neighbor communication lists for halo exchange
     * Identifies neighboring ranks and prepares send/recv buffers
     */
    void buildNeighborCommunication() {
        size_t numTotalNodes = domain_.getNodeCount();
        
        // Get node ownership
        const auto& nodeOwnership = domain_.getNodeOwnershipMap();
        thrust::host_vector<uint8_t> h_ownership(numTotalNodes);
        thrust::copy(thrust::device_pointer_cast(nodeOwnership.data()),
                    thrust::device_pointer_cast(nodeOwnership.data() + numTotalNodes),
                    h_ownership.begin());
        
        // Gather rank boundaries from all ranks to enable findRank locally
        std::vector<KeyType> rankBoundaries(numRanks_ + 1);
        KeyType myAssignmentStart = domain_.getDomain().assignmentStart();
        
        // Determine MPI datatype based on KeyType size
        MPI_Datatype mpiKeyType = (sizeof(KeyType) == 8) ? MPI_UNSIGNED_LONG : MPI_UNSIGNED;
        
        // Gather assignment start keys from all ranks
        MPI_Allgather(&myAssignmentStart, 1, mpiKeyType,
                     rankBoundaries.data(), 1, mpiKeyType,
                     MPI_COMM_WORLD);
        
        // Set the final boundary (maximum possible SFC key)
        rankBoundaries[numRanks_] = std::numeric_limits<KeyType>::max();
        
        // Lambda to find which rank owns a given SFC key (binary search)
        auto findRank = [&](KeyType key) -> int {
            auto it = std::upper_bound(rankBoundaries.begin(), rankBoundaries.end(), key);
            return static_cast<int>(it - rankBoundaries.begin()) - 1;
        };
        
        const auto& nodeSfcKeys = domain_.getLocalToGlobalSfcMap();
        thrust::host_vector<KeyType> h_sfcKeys(numTotalNodes);
        thrust::copy(thrust::device_pointer_cast(nodeSfcKeys.data()),
                    thrust::device_pointer_cast(nodeSfcKeys.data() + numTotalNodes),
                    h_sfcKeys.begin());
        
        // Build lists of ghost nodes grouped by owning rank
        std::map<int, std::vector<size_t>> ghostNodesByRank;
        for (size_t nodeIdx = 0; nodeIdx < numTotalNodes; ++nodeIdx) {
            if (h_ownership[nodeIdx] == 0) {  // Pure ghost node (state 0)
                KeyType sfcKey = h_sfcKeys[nodeIdx];
                int ownerRank = findRank(sfcKey);
                
                if (ownerRank >= 0 && ownerRank < numRanks_ && ownerRank != rank_) {
                    ghostNodesByRank[ownerRank].push_back(nodeIdx);
                }
            }
        }
        
        // Store neighbor ranks and communication lists
        neighborRanks_.clear();
        recvGhostNodes_.clear();
        
        for (auto& [neighborRank, ghostNodes] : ghostNodesByRank) {
            neighborRanks_.push_back(neighborRank);
            recvGhostNodes_[neighborRank] = ghostNodes;
        }
        
        // Now determine which of our owned DOFs need to be sent to neighbors
        // This requires knowing which neighbors will request data from us
        // Use MPI_Alltoall to exchange neighbor counts
        std::vector<int> sendCounts(numRanks_, 0);
        std::vector<int> recvCounts(numRanks_, 0);
        
        for (int neighborRank : neighborRanks_) {
            recvCounts[neighborRank] = recvGhostNodes_[neighborRank].size();
        }
        
        MPI_Alltoall(recvCounts.data(), 1, MPI_INT, sendCounts.data(), 1, MPI_INT, MPI_COMM_WORLD);
        
        // Now we know how many DOFs each neighbor needs from us
        // Receive the global DOF indices they need
        sendOwnedNodes_.clear();
        
        std::vector<MPI_Request> sendRequests, recvRequests;
        std::map<int, std::vector<KeyType>> sendGlobalDofs, recvGlobalDofs;
        
        // Post receives for global DOF indices from neighbors
        for (int neighborRank = 0; neighborRank < numRanks_; ++neighborRank) {
            if (sendCounts[neighborRank] > 0) {
                recvGlobalDofs[neighborRank].resize(sendCounts[neighborRank]);
                MPI_Request req;
                MPI_Irecv(recvGlobalDofs[neighborRank].data(), sendCounts[neighborRank],
                         mpiKeyType, neighborRank, /*tag=*/200, MPI_COMM_WORLD, &req);
                recvRequests.push_back(req);
            }
        }
        
        // Send our ghost global DOF indices to their owners
        for (int neighborRank : neighborRanks_) {
            if (recvCounts[neighborRank] > 0) {
                sendGlobalDofs[neighborRank].resize(recvCounts[neighborRank]);
                for (size_t i = 0; i < recvGhostNodes_[neighborRank].size(); ++i) {
                    size_t nodeIdx = recvGhostNodes_[neighborRank][i];
                    sendGlobalDofs[neighborRank][i] = nodeToGlobalDof_[nodeIdx];
                }
                MPI_Request req;
                MPI_Isend(sendGlobalDofs[neighborRank].data(), recvCounts[neighborRank],
                         mpiKeyType, neighborRank, /*tag=*/200, MPI_COMM_WORLD, &req);
                sendRequests.push_back(req);
            }
        }
        
        MPI_Waitall(recvRequests.size(), recvRequests.data(), MPI_STATUSES_IGNORE);
        MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE);
        
        // Map received global DOF indices to our owned local DOF indices
        for (auto& [neighborRank, globalDofs] : recvGlobalDofs) {
            std::vector<size_t> ownedLocalDofs;
            for (KeyType globalDof : globalDofs) {
                // Find local DOF index for this global DOF
                for (size_t localDof = 0; localDof < numLocalDofs_; ++localDof) {
                    if (localToGlobalDof_[localDof] == globalDof) {
                        ownedLocalDofs.push_back(localDof);
                        break;
                    }
                }
            }
            sendOwnedNodes_[neighborRank] = ownedLocalDofs;
        }
        
        if (rank_ == 0) {
            std::cout << "Ghost exchange: " << neighborRanks_.size() << " neighbors (rank 0)\n";
        }

        if (false && rank_ == 0) {
            for (int neighborRank : neighborRanks_) {
                std::cout << "  -> Neighbor " << neighborRank
                          << ": send " << sendOwnedNodes_[neighborRank].size()
                          << ", recv " << recvGhostNodes_[neighborRank].size() << "\n";
            }
        }
    }
    
    /**
     * @brief Update ghost DOF values via halo exchange
     * Call this before each SpMV to ensure ghost values are current
     * 
     * @param dofVector Vector of DOF values (size: numLocalDofs + numLocalGhostDofs)
     */
    template<typename VectorType>
    void updateGhostDofValues(VectorType& dofVector) {
        if (numRanks_ == 1 || neighborRanks_.empty()) return;  // No communication needed
        
        // Copy vector to host for MPI communication
        std::vector<RealType> h_dofVector(dofVector.size());
        thrust::copy(thrust::device_pointer_cast(dofVector.data()),
                    thrust::device_pointer_cast(dofVector.data() + dofVector.size()),
                    h_dofVector.begin());
        
        // Determine MPI datatype based on RealType
        MPI_Datatype mpiRealType = (sizeof(RealType) == 8) ? MPI_DOUBLE : MPI_FLOAT;
        
        std::vector<MPI_Request> sendRequests, recvRequests;
        std::map<int, std::vector<RealType>> sendBuffers, recvBuffers;
        
        // Post receives for ghost values
        for (int neighborRank : neighborRanks_) {
            size_t numRecv = recvGhostNodes_[neighborRank].size();
            if (numRecv > 0) {
                recvBuffers[neighborRank].resize(numRecv);
                MPI_Request req;
                MPI_Irecv(recvBuffers[neighborRank].data(), numRecv, mpiRealType,
                         neighborRank, /*tag=*/300, MPI_COMM_WORLD, &req);
                recvRequests.push_back(req);
            }
        }
        
        // Send owned values to neighbors
        for (auto& [neighborRank, ownedLocalDofs] : sendOwnedNodes_) {
            size_t numSend = ownedLocalDofs.size();
            if (numSend > 0) {
                sendBuffers[neighborRank].resize(numSend);
                for (size_t i = 0; i < numSend; ++i) {
                    sendBuffers[neighborRank][i] = h_dofVector[ownedLocalDofs[i]];
                }
                MPI_Request req;
                MPI_Isend(sendBuffers[neighborRank].data(), numSend, mpiRealType,
                         neighborRank, /*tag=*/300, MPI_COMM_WORLD, &req);
                sendRequests.push_back(req);
            }
        }
        
        // Wait for receives and update ghost values
        MPI_Waitall(recvRequests.size(), recvRequests.data(), MPI_STATUSES_IGNORE);
        
        for (auto& [neighborRank, ghostNodes] : recvGhostNodes_) {
            const auto& recvData = recvBuffers[neighborRank];
            for (size_t i = 0; i < ghostNodes.size(); ++i) {
                size_t nodeIdx = ghostNodes[i];
                size_t localGhostDof = nodeToLocalDof_[nodeIdx];  // Local index in extended vector
                h_dofVector[localGhostDof] = recvData[i];
            }
        }
        
        MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE);
        
        // Copy updated vector back to device
        thrust::copy(h_dofVector.begin(), h_dofVector.end(),
                    thrust::device_pointer_cast(dofVector.data()));
    }
    
    /**
     * @brief Set boundary node data for GPU-based boundary detection
     */
    void set_boundary_data(const std::vector<uint8_t>& h_boundary_nodes) {
        d_boundary_nodes_.resize(h_boundary_nodes.size());
        thrust::copy(h_boundary_nodes.begin(), h_boundary_nodes.end(), 
                    thrust::device_pointer_cast(d_boundary_nodes_.data()));
    }
    
    /**
     * @brief Iterate over owned DOFs on boundary (GPU-based computation)
     * Similar to structured version's boundary_owned_dof_iterate
     */
    void boundary_owned_dof_iterate(std::function<void(size_t)> func) {
        // Try GPU-based detection first
        if (!d_boundary_nodes_.empty()) {
            // GPU-based boundary detection using stored boundary data
            size_t numTotalNodes = domain_.getNodeCount();
            if (d_boundary_nodes_.size() != numTotalNodes) {
                std::cerr << "ERROR: Boundary data size mismatch: " << d_boundary_nodes_.size() 
                          << " vs node count " << numTotalNodes << std::endl;
                return;
            }
            
            // Get ownership map on device
            const auto& nodeOwnership = domain_.getNodeOwnershipMap();
            
            // Create device vectors for computation
            DeviceVector<size_t> d_node_indices(numTotalNodes);
            thrust::sequence(thrust::device, 
                           thrust::device_pointer_cast(d_node_indices.data()),
                           thrust::device_pointer_cast(d_node_indices.data() + numTotalNodes));
            
            // Find nodes that are owned (ownership == 1) AND on boundary (boundary == 1)
            DeviceVector<size_t> d_boundary_dof_indices(numTotalNodes);
            auto end = thrust::copy_if(thrust::device,
                                      thrust::device_pointer_cast(d_node_indices.data()),
                                      thrust::device_pointer_cast(d_node_indices.data() + numTotalNodes),
                                      thrust::device_pointer_cast(d_boundary_dof_indices.data()),
                                      [nodeOwnership_ptr = thrust::raw_pointer_cast(nodeOwnership.data()),
                                       boundary_ptr = thrust::raw_pointer_cast(d_boundary_nodes_.data())]
                                      __device__ (size_t nodeIdx) {
                                          return nodeOwnership_ptr[nodeIdx] == 1 && boundary_ptr[nodeIdx] == 1;
                                      });
            
            // Calculate number of boundary DOFs copied
            size_t num_boundary_dofs = thrust::distance(thrust::device_pointer_cast(d_boundary_dof_indices.data()), end);
            
            // Copy boundary DOF node indices to host
            HostVector<size_t> h_boundary_node_indices(num_boundary_dofs);
            thrust::copy(thrust::device_pointer_cast(d_boundary_dof_indices.data()), 
                        thrust::device_pointer_cast(d_boundary_dof_indices.data() + num_boundary_dofs),
                        h_boundary_node_indices.begin());
            
            // Convert node indices to local DOF indices and call func
            for (size_t boundary_node_idx : h_boundary_node_indices) {
                // Find the local DOF index for this boundary node
                KeyType local_dof_idx = nodeToLocalDof_[boundary_node_idx];
                if (local_dof_idx >= 0 && local_dof_idx < numLocalDofs_) {
                    func(static_cast<size_t>(local_dof_idx));
                }
            }
        } else {
            // Fallback: Use domain boundary info (host-based)
            // Assume domain has boundary info since it was reported
            const auto& boundaryNodes = domain_.getBoundaryNodes();
            if (boundaryNodes.size() != domain_.getNodeCount()) {
                std::cerr << "ERROR: Boundary nodes size mismatch: " << boundaryNodes.size() 
                          << " vs node count " << domain_.getNodeCount() << std::endl;
                return;
            }
            
            thrust::host_vector<uint8_t> h_boundary_nodes(domain_.getNodeCount());
            thrust::copy(thrust::device_pointer_cast(boundaryNodes.data()),
                         thrust::device_pointer_cast(boundaryNodes.data() + domain_.getNodeCount()),
                         h_boundary_nodes.begin());
            
            // Copy ownership to host
            const auto& nodeOwnership = domain_.getNodeOwnershipMap();
            thrust::host_vector<uint8_t> h_ownership(domain_.getNodeCount());
            thrust::copy(thrust::device_pointer_cast(nodeOwnership.data()),
                         thrust::device_pointer_cast(nodeOwnership.data() + domain_.getNodeCount()),
                         h_ownership.begin());
            
            // Iterate over owned nodes only and check if they are on boundary
            size_t ownedNodeCounter = 0;
            for (size_t localNodeIdx = 0; localNodeIdx < domain_.getNodeCount(); ++localNodeIdx) {
                // Skip ghost nodes (not owned by this rank)
                if (h_ownership[localNodeIdx] != 1) continue;
                
                // Check if node is on boundary (topological detection like MFEM)
                if (h_boundary_nodes[localNodeIdx] == 1) {
                    // This is a boundary DOF - call the func with LOCAL DOF index
                    size_t localDofIdx = ownedNodeCounter;  // Local DOF index (0 to numLocalDofs-1)
                    func(localDofIdx);
                }
                ownedNodeCounter++;
            }
        }
    }
    
    /**
     * @brief Get node-to-local-DOF mapping
     * Maps node index to local DOF index (0 to numLocalDofs-1)
     * Returns -1 for ghost nodes
     */
    const std::vector<KeyType>& get_node_to_local_dof() const {
        return nodeToLocalDof_;
    }
    
    /**
     * @brief Get node-to-global-DOF mapping
     * Maps node index to global DOF index
     * Returns -1 for ghost nodes (unless ghost exchange is implemented)
     */
    const std::vector<KeyType>& get_node_to_global_dof() const {
        return nodeToGlobalDof_;
    }

    /**
     * @brief Get resolved ownership (after MPI consensus)
     * 1 = truly owned by this rank, 0 = ghost (including nodes lost in tie-breaking)
     */
    const std::vector<uint8_t>& get_resolved_ownership() const {
        return resolvedOwnership_;
    }

    /**
     * @brief Get coordinates for a DOF using SFC key
     * For unstructured meshes, SFC keys correspond to node indices
     */
    void get_dof_coordinates_from_sfc(KeyType sfc, double* point) const {
        // For unstructured meshes, SFC keys are node indices
        size_t nodeIdx = static_cast<size_t>(sfc);
        
        // Get coordinates from domain
        const auto& d_x = domain_.getNodeX();
        const auto& d_y = domain_.getNodeY();
        const auto& d_z = domain_.getNodeZ();
        
        // Copy to host for access (inefficient but functional)
        thrust::host_vector<float> h_x(1), h_y(1), h_z(1);
        thrust::copy(thrust::device_pointer_cast(d_x.data() + nodeIdx),
                     thrust::device_pointer_cast(d_x.data() + nodeIdx + 1),
                     h_x.begin());
        thrust::copy(thrust::device_pointer_cast(d_y.data() + nodeIdx),
                     thrust::device_pointer_cast(d_y.data() + nodeIdx + 1),
                     h_y.begin());
        thrust::copy(thrust::device_pointer_cast(d_z.data() + nodeIdx),
                     thrust::device_pointer_cast(d_z.data() + nodeIdx + 1),
                     h_z.begin());
        
        point[0] = h_x[0];
        point[1] = h_y[0];
        point[2] = h_z[0];
    }
    
    /**
     * @brief Convert local DOF index to global DOF index
     */
    size_t local_to_global(size_t local_dof) const {
        return localToGlobalDof_[local_dof];
    }
    
    /**
     * @brief Get total number of global DOFs
     */
    size_t get_num_global_dofs() const {
        return numGlobalDofs_;
    }
    
    /**
     * @brief Get number of local (owned) DOFs
     */
    size_t get_num_local_dofs() const {
        return numLocalDofs_;
    }
    
    /**
     * @brief Get total number of local DOFs including ghosts (for local matrix size)
     */
    size_t get_num_local_dofs_with_ghosts() const {
        return get_num_local_dofs() + numLocalGhostDofs_;
    }
    
    /**
     * @brief Map global DOF index to local ghost index (returns -1 if not a ghost)
     */
    KeyType global_to_local_ghost(KeyType globalDof) const {
        auto it = globalToLocalGhostDof_.find(globalDof);
        return (it != globalToLocalGhostDof_.end()) ? it->second : static_cast<KeyType>(-1);
    }
    
    /**
     * @brief Get global DOF index for a local ghost DOF index
     */
    KeyType get_ghost_global_dof(KeyType localGhostIdx) const {
        if (localGhostIdx >= ghostDofToNode_.size()) {
            return static_cast<KeyType>(-1);
        }
        KeyType nodeIdx = ghostDofToNode_[localGhostIdx];
        return nodeToGlobalDof_[nodeIdx];
    }
    
    /**
     * @brief Get owned DOF range
     */
    std::pair<size_t, size_t> get_owned_dof_range() const {
        return {ownedDofStart_, ownedDofEnd_};
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
        size_t expectedSize = get_num_local_dofs_with_ghosts();
        if (dofValues.size() != expectedSize) {
            std::cerr << "Rank " << rank_ << ": Error - vector size " << dofValues.size() 
                      << " != expected " << expectedSize << " (owned+ghost)\n";
            throw std::runtime_error("UnstructuredDofHandler::scatterAddGhostData: vector size mismatch");
        }
        
        // Use robust global exchange (Allgather) instead of broken point-to-point
        exchangeNodeDataGlobal(dofValues);
    }
    
    size_t numBoundaryNodes() const { return boundaryNodes_.size(); }
    size_t numGhostNodes() const { return ghostNodes_.size(); }
    bool isInitialized() const { return initialized_; }

    /**
     * @brief Robust global exchange of ghost contributions
     * 
     * 1. Collects all local ghost contributions (GlobalDOF, Value)
     * 2. Allgathers them to all ranks
     * 3. Each rank adds contributions for DOFs it owns
     */
    template<typename VectorType>
    void exchangeNodeDataGlobal(VectorType& dofValues) {
        using T = typename VectorType::value_type;
        
        // 1. Extract ghost values from device
        size_t numGhosts = numLocalGhostDofs_;
        if (numGhosts == 0) {
            // Even if we have no ghosts, we must participate in Allgather
            // But we can send empty data
        }

        // Copy ghost values to host
        // Ghost DOFs are stored after owned DOFs: [0..numLocalDofs_-1] owned, [numLocalDofs_..end] ghost
        std::vector<T> h_ghostValues(numGhosts);
        if (numGhosts > 0) {
            thrust::copy(thrust::device_pointer_cast(dofValues.data()) + numLocalDofs_, 
                         thrust::device_pointer_cast(dofValues.data()) + numLocalDofs_ + numGhosts, 
                         h_ghostValues.begin());
        }

        // Prepare send buffers: (GlobalDOF, Value) pairs
        // We split into two arrays for simpler MPI types
        std::vector<KeyType> sendGlobalDofs;
        std::vector<T> sendValues;
        sendGlobalDofs.reserve(numGhosts);
        sendValues.reserve(numGhosts);

        int skipped_invalid = 0;
        int skipped_zero = 0;
        for (size_t i = 0; i < numGhosts; ++i) {
            T val = h_ghostValues[i];
            if (std::abs(val) < 1e-20) {
                skipped_zero++;
                continue; // Optimization: only send non-zero contributions
            }
            KeyType globalDof = get_ghost_global_dof(i);
            if (globalDof == static_cast<KeyType>(-1)) {
                skipped_invalid++;
                if (rank_ == 0 && skipped_invalid <= 5) {
                    std::cerr << "Rank " << rank_ << ": Ghost DOF " << i << " has invalid global DOF (-1), value=" << val << std::endl;
                }
                continue;
            }
            sendGlobalDofs.push_back(globalDof);
            sendValues.push_back(val);
        }
        
        if (rank_ == 0 || rank_ == 1) {
            std::cout << "Rank " << rank_ << ": Ghost exchange - " << numGhosts << " ghosts, "
                      << skipped_zero << " zero, " << skipped_invalid << " invalid, "
                      << sendGlobalDofs.size() << " to send" << std::endl;
        }

        // 2. Allgatherv
        // First exchange counts
        int myCount = static_cast<int>(sendGlobalDofs.size());
        std::vector<int> allCounts(numRanks_);
        MPI_Allgather(&myCount, 1, MPI_INT, allCounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

        // Calculate displacements
        std::vector<int> displs(numRanks_);
        int totalCount = 0;
        for (int i = 0; i < numRanks_; ++i) {
            displs[i] = totalCount;
            totalCount += allCounts[i];
        }

        // Exchange Keys (GlobalDOFs)
        std::vector<KeyType> recvGlobalDofs(totalCount);
        // Use MPI_UINT64_T for KeyType (assuming uint64_t)
        MPI_Allgatherv(sendGlobalDofs.data(), myCount, MPI_UINT64_T,
                       recvGlobalDofs.data(), allCounts.data(), displs.data(), MPI_UINT64_T, MPI_COMM_WORLD);

        // Exchange Values
        std::vector<T> recvValues(totalCount);
        MPI_Allgatherv(sendValues.data(), myCount, getMPIType<T>(),
                       recvValues.data(), allCounts.data(), displs.data(), getMPIType<T>(), MPI_COMM_WORLD);

        // 3. Accumulate to owned DOFs
        std::vector<T> h_ownedUpdates(numLocalDofs_, 0.0);
        bool anyUpdates = false;
        int myContributions = 0;

        for (int i = 0; i < totalCount; ++i) {
            KeyType globalDof = recvGlobalDofs[i];
            T val = recvValues[i];

            // Check if I own this DOF
            if (globalDof >= ownedDofStart_ && globalDof < ownedDofEnd_) {
                size_t localIdx = globalDof - ownedDofStart_;
                h_ownedUpdates[localIdx] += val;
                anyUpdates = true;
                myContributions++;
            }
        }
        
        if (rank_ == 0 || rank_ == 1) {
            std::cout << "Rank " << rank_ << ": Received " << totalCount << " contributions, "
                      << myContributions << " for my DOFs [" << ownedDofStart_ << ", " << ownedDofEnd_ << ")" << std::endl;
        }

        // 4. Update device vector (owned part)
        if (anyUpdates) {
            // We need to add h_ownedUpdates to the existing values on device
            // Copy updates to device
            DeviceVector<T> d_updates = h_ownedUpdates;
            
            // Add to owned part of dofValues
            thrust::transform(thrust::device_pointer_cast(dofValues.data()), 
                              thrust::device_pointer_cast(dofValues.data()) + numLocalDofs_,
                              thrust::device_pointer_cast(d_updates.data()),
                              thrust::device_pointer_cast(dofValues.data()),
                              thrust::plus<T>());
        }
    }
    
    /**
     * @brief Exchange node data via GPU-aware MPI point-to-point communication
     * 
     * TODO: CRITICAL BUG - This implementation has buffer size mismatches!
     * - Rank A sends numBoundary_A values to Rank B
     * - Rank B expects numGhost_B values from Rank A
     * - If numBoundary_A != numGhost_B, MPI error/hang occurs
     * 
     * Correct approach:
     * 1. Build rank-to-rank neighbor graph from halo structure
     * 2. Exchange buffer sizes between neighbors
     * 3. Use matched send/recv pairs or MPI_Alltoallv with displacements
     * 
     * Pattern (when fixed):
     * 1. Pack boundary node values into device send buffer
     * 2. MPI_Isend device buffer to neighbor ranks (GPU-aware MPI)
     * 3. MPI_Irecv into device receive buffer from neighbors
     * 4. Wait and atomically add received ghost contributions on GPU
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
        
        // Validate inputs
        if (nodeOwnership.size() != nodeCount) {
            std::cerr << "Rank " << rank_ << ": Error in buildNodeCommunicationPatterns - nodeOwnership size mismatch\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        if (haloCount > 0) {
            // Check if connectivity data is available
            if (std::get<0>(conn_tuple).size() == 0) {
                std::cerr << "Rank " << rank_ << ": Error - connectivity data not available for halo processing\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            
            // Validate halo indices are within bounds
            size_t maxElementIndex = domain_.getElementCount();
            thrust::host_vector<KeyType> h_haloIndices(haloCount);
            thrust::copy(thrust::device_pointer_cast(haloElementIndices.data()),
                        thrust::device_pointer_cast(haloElementIndices.data() + haloCount),
                        h_haloIndices.begin());
            
            for (size_t i = 0; i < haloCount; ++i) {
                if (h_haloIndices[i] >= maxElementIndex) {
                    std::cerr << "Rank " << rank_ << ": Error - halo element index " << h_haloIndices[i] 
                              << " >= max elements " << maxElementIndex << "\n";
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
            }
        }
        
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
                                                     uint8_t own = ownership_ptr[nodeIdx];
                                                     return own == 1 || own == 2;  // Owned (1) or shared (2) nodes
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
    
    // Distributed DOF management
    HostVector<size_t> localToGlobalDof_;  // Local to global DOF mapping
    size_t numGlobalDofs_;                 // Total global DOFs
    size_t ownedDofStart_;                 // Start of owned DOF range
    size_t ownedDofEnd_;                   // End of owned DOF range
    
    // Node-to-DOF mappings for assembly
    std::vector<KeyType> nodeToLocalDof_;   // Node index -> local DOF index (row)
    std::vector<KeyType> nodeToGlobalDof_;  // Node index -> global DOF index (column)
    std::vector<KeyType> ghostDofToNode_;   // Local ghost DOF index -> Node index
    std::unordered_map<KeyType, KeyType> sfcToGlobalDof_;  // SFC key -> global DOF index
    std::vector<uint8_t> resolvedOwnership_;  // Resolved ownership after MPI consensus
    
    // Ghost DOF management for local-local matrix assembly
    std::map<KeyType, KeyType> globalToLocalGhostDof_;  // Global DOF -> local ghost DOF index
    size_t numLocalDofs_;        // Number of owned DOFs on this rank
    size_t numLocalGhostDofs_;   // Number of ghost DOFs (numLocalDofs + numGhostDofs = total local storage)
    
    // Neighbor communication for halo exchange
    std::vector<int> neighborRanks_;                              // List of neighboring MPI ranks
    std::map<int, std::vector<size_t>> recvGhostNodes_;          // rank -> local ghost node indices to receive
    std::map<int, std::vector<size_t>> sendOwnedNodes_;          // rank -> owned local DOF indices to send
    
    // Boundary detection data
    DeviceVector<uint8_t> d_boundary_nodes_;  // Boundary flags on device for GPU computation
};

} // namespace fem
} // namespace mars
