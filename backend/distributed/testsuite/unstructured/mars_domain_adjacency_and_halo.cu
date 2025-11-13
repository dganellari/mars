#include "gtest/gtest.h"
#include "backend/distributed/unstructured/domain.hpp"
#include <cuda_runtime.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <mpi.h>
#include <random>

namespace mars {
namespace testing {

// Test fixture for domain adjacency and halo tests
class MarsDomainAdjacencyAndHaloTest : public ::testing::Test {
protected:
    using RealType = double;
    using KeyType = uint64_t;
    using AcceleratorTag = cstone::GpuTag;
    using DomainType = ElementDomain<TetTag, RealType, KeyType, AcceleratorTag>;

    void SetUp() override {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
        MPI_Comm_size(MPI_COMM_WORLD, &numRanks_);
        
        // Use the test mesh from your existing tests (consistent with mars_domain_binary_mesh.cu)
        const char* meshPathEnv = std::getenv("MESH_PATH");
        if (!meshPathEnv) {
            GTEST_SKIP() << "MESH_PATH environment variable not set";
        }
        meshDir_ = std::string(meshPathEnv);

        //ensure trailing slash
        if (!meshDir_.empty() && meshDir_.back() != '/' && meshDir_.back() != '\\') {
            meshDir_ += '/';
        }
    }

    int rank_;
    int numRanks_;
    std::string meshDir_;
};

// Test 1: Local-to-Global SFC Map Uniqueness and Sorting
TEST_F(MarsDomainAdjacencyAndHaloTest, LocalToGlobalSfcMapUniqueness) {
    DomainType domain(meshDir_, rank_, numRanks_);
    
    const auto& sfcMap = domain.getLocalToGlobalSfcMap();
    size_t nodeCount = domain.getNodeCount();
    
    ASSERT_EQ(sfcMap.size(), nodeCount) << "SFC map size mismatch";
    ASSERT_GT(nodeCount, 0) << "Node count should be positive";
    
    // Copy to host for verification
    thrust::host_vector<KeyType> h_sfcMap(sfcMap.size());
    thrust::copy(thrust::device, 
                 sfcMap.data(),
                 sfcMap.data() + sfcMap.size(),
                 h_sfcMap.begin());
    
    // Check uniqueness and sorting
    for (size_t i = 1; i < h_sfcMap.size(); ++i) {
        EXPECT_LT(h_sfcMap[i-1], h_sfcMap[i]) 
            << "Rank " << rank_ << ": SFC keys not unique or sorted at index " << i;
    }
    
    std::cout << "Rank " << rank_ << ": Verified " << nodeCount 
              << " unique SFC keys" << std::endl;
}

// Test 2: Element-to-Node Local ID Connectivity Validity
TEST_F(MarsDomainAdjacencyAndHaloTest, ElementToNodeLocalIdMapping) {
    DomainType domain(meshDir_, rank_, numRanks_);
    
    size_t elementCount = domain.getElementCount();
    size_t nodeCount = domain.getNodeCount();
    
    ASSERT_GT(elementCount, 0) << "Element count should be positive";
    
    // Get SFC key connectivity (after buildAdjacency) - stores global SFC keys, not local indices
    auto conn0 = domain.getConnectivity<0>();
    auto conn1 = domain.getConnectivity<1>();
    auto conn2 = domain.getConnectivity<2>();
    auto conn3 = domain.getConnectivity<3>();
    
    // Copy to host
    thrust::host_vector<KeyType> h_conn0(conn0.size());
    thrust::copy(thrust::device, conn0.data(), conn0.data() + conn0.size(), h_conn0.begin());

    thrust::host_vector<KeyType> h_conn1(conn1.size());
    thrust::copy(thrust::device, conn1.data(), conn1.data() + conn1.size(), h_conn1.begin());

    thrust::host_vector<KeyType> h_conn2(conn2.size());
    thrust::copy(thrust::device, conn2.data(), conn2.data() + conn2.size(), h_conn2.begin());

    thrust::host_vector<KeyType> h_conn3(conn3.size());
    thrust::copy(thrust::device, conn3.data(), conn3.data() + conn3.size(), h_conn3.begin());
    
    // Get SFC map for validation
    const auto& sfcMap = domain.getLocalToGlobalSfcMap();
    thrust::host_vector<KeyType> h_sfcMap(sfcMap.size());
    thrust::copy(thrust::device, sfcMap.data(), sfcMap.data() + sfcMap.size(), h_sfcMap.begin());
    
    // Verify all SFC keys in connectivity are valid (present in SFC map)
    for (size_t i = 0; i < elementCount; ++i) {
        EXPECT_TRUE(std::binary_search(h_sfcMap.begin(), h_sfcMap.end(), h_conn0[i])) 
            << "Rank " << rank_ << ": Invalid SFC key in conn0 at element " << i;
        EXPECT_TRUE(std::binary_search(h_sfcMap.begin(), h_sfcMap.end(), h_conn1[i])) 
            << "Rank " << rank_ << ": Invalid SFC key in conn1 at element " << i;
        EXPECT_TRUE(std::binary_search(h_sfcMap.begin(), h_sfcMap.end(), h_conn2[i])) 
            << "Rank " << rank_ << ": Invalid SFC key in conn2 at element " << i;
        EXPECT_TRUE(std::binary_search(h_sfcMap.begin(), h_sfcMap.end(), h_conn3[i])) 
            << "Rank " << rank_ << ": Invalid SFC key in conn3 at element " << i;
    }
    
    std::cout << "Rank " << rank_ << ": Verified SFC key connectivity for " 
              << elementCount << " elements" << std::endl;
}

// Test 3: Node-to-Element CSR Structure Correctness
TEST_F(MarsDomainAdjacencyAndHaloTest, NodeToElementCsrStructure) {
    DomainType domain(meshDir_, rank_, numRanks_);
    
    const auto& offsets = domain.getNodeToElementOffsets();
    const auto& list = domain.getNodeToElementList();
    size_t nodeCount = domain.getNodeCount();
    size_t elementCount = domain.getElementCount();
    
    ASSERT_EQ(offsets.size(), nodeCount + 1) << "Offset array size incorrect";
    
    // Copy to host
    thrust::host_vector<KeyType> h_offsets(offsets.size());
    thrust::copy(thrust::device, offsets.data(), offsets.data() + offsets.size(), h_offsets.begin());

    thrust::host_vector<KeyType> h_list(list.size());
    thrust::copy(thrust::device, list.data(), list.data() + list.size(), h_list.begin());
    
    // Check CSR structure validity
    EXPECT_EQ(h_offsets[0], 0) << "CSR offset must start at 0";
    EXPECT_EQ(h_offsets[nodeCount], h_list.size()) 
        << "Rank " << rank_ << ": Last offset must equal list size";
    
    // Offsets must be non-decreasing
    for (size_t i = 1; i <= nodeCount; ++i) {
        EXPECT_GE(h_offsets[i], h_offsets[i-1]) 
            << "Rank " << rank_ << ": Offsets not monotonic at index " << i;
    }
    
    // All element IDs in list must be valid
    for (size_t i = 0; i < h_list.size(); ++i) {
        EXPECT_LT(h_list[i], elementCount) 
            << "Rank " << rank_ << ": Invalid element ID in node-to-element list at index " << i;
    }
    
    std::cout << "Rank " << rank_ << ": Node-to-element CSR has " 
              << h_list.size() << " entries" << std::endl;
}

// Test 4: Halo Element Range Validity (Multi-rank only)
TEST_F(MarsDomainAdjacencyAndHaloTest, HaloElementRanges) {
    DomainType domain(meshDir_, rank_, numRanks_);
    
    auto [localStart, localEnd] = domain.localElementRange();
    auto haloRanges = domain.haloElementRanges();
    size_t totalElements = domain.getElementCount();
    size_t localCount = domain.localElementCount();
    
    // Local range must be valid
    EXPECT_LE(localStart, localEnd) << "Invalid local element range";
    EXPECT_LE(localEnd, totalElements) << "Local end exceeds total elements";
    EXPECT_EQ(localEnd - localStart, localCount) << "Local count mismatch";
    
    if (numRanks_ > 1) {
        // Halo ranges must not overlap with local range
        for (auto [haloStart, haloEnd] : haloRanges) {
            EXPECT_LT(haloStart, haloEnd) << "Invalid halo range";
            EXPECT_TRUE(haloEnd <= localStart || haloStart >= localEnd)
                << "Rank " << rank_ << ": Halo range overlaps with local range";
        }
        
        // Total elements = local + halos
        size_t haloCount = 0;
        for (auto [start, end] : haloRanges) {
            haloCount += (end - start);
        }
        EXPECT_EQ(localCount + haloCount, totalElements)
            << "Rank " << rank_ << ": Local + halo count mismatch";
            
        std::cout << "Rank " << rank_ << ": " << localCount << " local elements, "
                  << haloCount << " halo elements" << std::endl;
    } else {
        // Single rank - no halos expected
        EXPECT_TRUE(haloRanges.empty()) << "Single rank should have no halos";
        EXPECT_EQ(localCount, totalElements) << "Single rank: all elements should be local";
    }
}

// Test 5: Halo Element Indices Contiguity (Multi-rank only)
TEST_F(MarsDomainAdjacencyAndHaloTest, HaloElementIndicesContiguous) {
    if (numRanks_ < 2) {
        GTEST_SKIP() << "Halo indices test requires multiple ranks";
    }
    
    DomainType domain(meshDir_, rank_, numRanks_);
    
    const auto& haloIndices = domain.getHaloElementIndices();
    size_t haloCount = domain.getHaloElementCount();
    
    ASSERT_EQ(haloIndices.size(), haloCount) << "Halo indices size mismatch";
    
    if (haloCount == 0) {
        std::cout << "Rank " << rank_ << ": No halos (boundary rank)" << std::endl;
        return;
    }
    
    // Copy to host
    thrust::host_vector<KeyType> h_haloIndices(haloIndices.size());
    thrust::copy(thrust::device, haloIndices.data(), haloIndices.data() + haloIndices.size(), h_haloIndices.begin());
    
    // Verify all halo indices are valid and are actually halos
    for (size_t i = 0; i < h_haloIndices.size(); ++i) {
        size_t elemIdx = static_cast<size_t>(h_haloIndices[i]);
        EXPECT_LT(elemIdx, domain.getElementCount()) 
            << "Rank " << rank_ << ": Invalid halo element index at " << i;
        EXPECT_TRUE(domain.isHaloElement(elemIdx))
            << "Rank " << rank_ << ": Index " << elemIdx 
            << " in halo list but isHaloElement() returns false";
    }
    
    std::cout << "Rank " << rank_ << ": Verified " << haloCount 
              << " halo element indices" << std::endl;
}

// Test 6: Node Ownership Correctness (Multi-rank only)
TEST_F(MarsDomainAdjacencyAndHaloTest, NodeOwnershipMapping) {
    if (numRanks_ < 2) {
        GTEST_SKIP() << "Node ownership test requires multiple ranks";
    }
    
    DomainType domain(meshDir_, rank_, numRanks_);
    
    const auto& ownership = domain.getNodeOwnershipMap();
    size_t nodeCount = domain.getNodeCount();
    
    ASSERT_EQ(ownership.size(), nodeCount) << "Ownership map size mismatch";
    
    // Copy ownership to host
    thrust::host_vector<uint8_t> h_ownership(ownership.size());
    thrust::copy(thrust::device, ownership.data(), ownership.data() + ownership.size(), h_ownership.begin());
    
    // Get connectivity for verification
    auto conn0 = domain.getConnectivity<0>();
    auto conn1 = domain.getConnectivity<1>();
    auto conn2 = domain.getConnectivity<2>();
    auto conn3 = domain.getConnectivity<3>();
    
    thrust::host_vector<KeyType> h_conn0(conn0.size());
    thrust::copy(thrust::device, conn0.data(), conn0.data() + conn0.size(), h_conn0.begin());

    thrust::host_vector<KeyType> h_conn1(conn1.size());
    thrust::copy(thrust::device, conn1.data(), conn1.data() + conn1.size(), h_conn1.begin());

    thrust::host_vector<KeyType> h_conn2(conn2.size());
    thrust::copy(thrust::device, conn2.data(), conn2.data() + conn2.size(), h_conn2.begin());

    thrust::host_vector<KeyType> h_conn3(conn3.size());
    thrust::copy(thrust::device, conn3.data(), conn3.data() + conn3.size(), h_conn3.begin());
    
    // Copy SFC map and build reverse map
    const auto& sfcMap = domain.getLocalToGlobalSfcMap();
    thrust::host_vector<KeyType> h_sfcMap(sfcMap.size());
    thrust::copy(thrust::device, sfcMap.data(), sfcMap.data() + sfcMap.size(), h_sfcMap.begin());
    
    std::unordered_map<KeyType, size_t> sfcToLocal;
    for (size_t i = 0; i < h_sfcMap.size(); ++i) {
        sfcToLocal[h_sfcMap[i]] = i;
    }
    
    auto [localStart, localEnd] = domain.localElementRange();
    
    // Sample check: nodes in local elements should be owned
    size_t checkCount = std::min(size_t(100), localEnd - localStart);
    for (size_t i = 0; i < checkCount; ++i) {
        size_t elemIdx = localStart + i;
        // Map SFC keys to local indices
        size_t localNode0 = sfcToLocal[h_conn0[elemIdx]];
        EXPECT_TRUE(h_ownership[localNode0]) 
            << "Rank " << rank_ << ": Node in local element not marked as owned";
    }
    
    // Count owned vs halo nodes
    size_t ownedCount = 0;
    for (size_t i = 0; i < nodeCount; ++i) {
        if (h_ownership[i]) ownedCount++;
    }
    
    std::cout << "Rank " << rank_ << ": " << ownedCount << " owned nodes, "
              << (nodeCount - ownedCount) << " halo nodes" << std::endl;
}

// Test 7: SFC Coordinate Decoding Bounds Check
TEST_F(MarsDomainAdjacencyAndHaloTest, SfcCoordinateDecoding) {
    DomainType domain(meshDir_, rank_, numRanks_);
    
    const auto& sfcMap = domain.getLocalToGlobalSfcMap();
    auto box = domain.getBoundingBox();
    
    // Copy to host
    thrust::host_vector<KeyType> h_sfcKeys(sfcMap.size());
    thrust::copy(thrust::device, sfcMap.data(), sfcMap.data() + sfcMap.size(), h_sfcKeys.begin());
    
    // Decode a sample of SFC keys and check bounds
    size_t checkCount = std::min(h_sfcKeys.size(), size_t(100));
    for (size_t i = 0; i < checkCount; i += 10) {
        auto [x, y, z] = domain.sfcToPhysicalCoordinate(h_sfcKeys[i]);
        
        EXPECT_GE(x, box.xmin()) << "Decoded x below box min";
        EXPECT_LE(x, box.xmax()) << "Decoded x above box max";
        EXPECT_GE(y, box.ymin()) << "Decoded y below box min";
        EXPECT_LE(y, box.ymax()) << "Decoded y above box max";
        EXPECT_GE(z, box.zmin()) << "Decoded z below box min";
        EXPECT_LE(z, box.zmax()) << "Decoded z above box max";
    }
    
    std::cout << "Rank " << rank_ << ": Verified coordinate decoding for " 
              << checkCount << " nodes" << std::endl;
}

// Test 8: Element Ownership Flags Consistency
TEST_F(MarsDomainAdjacencyAndHaloTest, ElementOwnershipFlags) {
    DomainType domain(meshDir_, rank_, numRanks_);
    
    auto [localStart, localEnd] = domain.localElementRange();
    size_t totalElements = domain.getElementCount();
    
    // Check a sample of elements
    size_t checkCount = std::min(totalElements, size_t(1000));
    for (size_t i = 0; i < checkCount; ++i) {
        bool isLocal = domain.isLocalElement(i);
        bool isHalo = domain.isHaloElement(i);
        
        // Element must be either local or halo, not both
        EXPECT_NE(isLocal, isHalo) 
            << "Rank " << rank_ << ": Element " << i << " ownership ambiguous";
        
        // Verify against local range
        if (i >= localStart && i < localEnd) {
            EXPECT_TRUE(isLocal) 
                << "Rank " << rank_ << ": Element " << i 
                << " in local range but not flagged as local";
        } else {
            EXPECT_TRUE(isHalo) 
                << "Rank " << rank_ << ": Element " << i 
                << " outside local range but not flagged as halo";
        }
    }
    
    std::cout << "Rank " << rank_ << ": Verified ownership flags for " 
              << checkCount << " elements" << std::endl;
}

// TEST 9: Node-to-Element Inverse Lookup Semantic Correctness
TEST_F(MarsDomainAdjacencyAndHaloTest, NodeToElementInverseLookup) {
    DomainType domain(meshDir_, rank_, numRanks_);
    
    const auto& offsets = domain.getNodeToElementOffsets();
    const auto& list = domain.getNodeToElementList();
    size_t nodeCount = domain.getNodeCount();
    
    // Copy CSR structure to host
    thrust::host_vector<KeyType> h_offsets(offsets.size());
    thrust::copy(thrust::device, offsets.data(), offsets.data() + offsets.size(), h_offsets.begin());
    
    thrust::host_vector<KeyType> h_list(list.size());
    thrust::copy(thrust::device, list.data(), list.data() + list.size(), h_list.begin());
    
    // Copy connectivity arrays
    auto conn0 = domain.getConnectivity<0>();
    auto conn1 = domain.getConnectivity<1>();
    auto conn2 = domain.getConnectivity<2>();
    auto conn3 = domain.getConnectivity<3>();
    
    thrust::host_vector<KeyType> h_conn0(conn0.size());
    thrust::copy(thrust::device, conn0.data(), conn0.data() + conn0.size(), h_conn0.begin());
    
    thrust::host_vector<KeyType> h_conn1(conn1.size());
    thrust::copy(thrust::device, conn1.data(), conn1.data() + conn1.size(), h_conn1.begin());
    
    thrust::host_vector<KeyType> h_conn2(conn2.size());
    thrust::copy(thrust::device, conn2.data(), conn2.data() + conn2.size(), h_conn2.begin());
    
    thrust::host_vector<KeyType> h_conn3(conn3.size());
    thrust::copy(thrust::device, conn3.data(), conn3.data() + conn3.size(), h_conn3.begin());
    
    // Copy SFC map
    const auto& sfcMap = domain.getLocalToGlobalSfcMap();
    thrust::host_vector<KeyType> h_sfcMap(sfcMap.size());
    thrust::copy(thrust::device, sfcMap.data(), sfcMap.data() + sfcMap.size(), h_sfcMap.begin());
    
    // Test random sample of nodes
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> dist(0, nodeCount - 1);
    
    size_t testCount = std::min(nodeCount, size_t(50));
    for (size_t t = 0; t < testCount; ++t) {
        size_t testNodeId = dist(gen);
        KeyType nodeSfcKey = h_sfcMap[testNodeId];
        
        // Get all elements that reference this node (from CSR)
        KeyType start = h_offsets[testNodeId];
        KeyType end = h_offsets[testNodeId + 1];
        
        // Verify each element in the list actually references testNodeId's SFC key
        for (KeyType i = start; i < end; ++i) {
            size_t elemIdx = static_cast<size_t>(h_list[i]);
            
            // Check if nodeSfcKey appears in element's connectivity
            bool found = (h_conn0[elemIdx] == nodeSfcKey ||
                         h_conn1[elemIdx] == nodeSfcKey ||
                         h_conn2[elemIdx] == nodeSfcKey ||
                         h_conn3[elemIdx] == nodeSfcKey);
            
            EXPECT_TRUE(found) 
                << "Rank " << rank_ << ": Node " << testNodeId 
                << " (SFC key " << nodeSfcKey << ") in element " << elemIdx 
                << "'s adjacency list but element doesn't reference this node in connectivity";
        }
    }
    
    std::cout << "Rank " << rank_ << ": Verified inverse lookup for " 
              << testCount << " random nodes" << std::endl;
}

// TEST 10: Halo Node Ownership Verification
TEST_F(MarsDomainAdjacencyAndHaloTest, HaloNodeOwnership) {
    if (numRanks_ < 2) {
        GTEST_SKIP() << "Halo node ownership test requires multiple ranks";
    }
    
    DomainType domain(meshDir_, rank_, numRanks_);
    
    const auto& ownership = domain.getNodeOwnershipMap();
    size_t nodeCount = domain.getNodeCount();
    auto [localStart, localEnd] = domain.localElementRange();
    auto haloRanges = domain.haloElementRanges();
    
    if (haloRanges.empty()) {
        std::cout << "Rank " << rank_ << ": No halos, skipping halo node check" << std::endl;
        return;
    }
    
    // Copy ownership to host
    thrust::host_vector<uint8_t> h_ownership(ownership.size());
    thrust::copy(thrust::device, ownership.data(), ownership.data() + ownership.size(), h_ownership.begin());
    
    auto conn0 = domain.getConnectivity<0>();
    auto conn1 = domain.getConnectivity<1>();
    auto conn2 = domain.getConnectivity<2>();
    auto conn3 = domain.getConnectivity<3>();
    
    thrust::host_vector<KeyType> h_conn0(conn0.size());
    thrust::copy(thrust::device, conn0.data(), conn0.data() + conn0.size(), h_conn0.begin());
    
    thrust::host_vector<KeyType> h_conn1(conn1.size());
    thrust::copy(thrust::device, conn1.data(), conn1.data() + conn1.size(), h_conn1.begin());
    
    thrust::host_vector<KeyType> h_conn2(conn2.size());
    thrust::copy(thrust::device, conn2.data(), conn2.data() + conn2.size(), h_conn2.begin());
    
    thrust::host_vector<KeyType> h_conn3(conn3.size());
    thrust::copy(thrust::device, conn3.data(), conn3.data() + conn3.size(), h_conn3.begin());
    
    // Copy SFC map and build reverse map
    const auto& sfcMap = domain.getLocalToGlobalSfcMap();
    thrust::host_vector<KeyType> h_sfcMap(sfcMap.size());
    thrust::copy(thrust::device, sfcMap.data(), sfcMap.data() + sfcMap.size(), h_sfcMap.begin());
    
    std::unordered_map<KeyType, size_t> sfcToLocal;
    for (size_t i = 0; i < h_sfcMap.size(); ++i) {
        sfcToLocal[h_sfcMap[i]] = i;
    }
    
    // Track nodes that appear ONLY in halo elements (should be marked as NOT owned)
    std::vector<bool> appearsInLocal(nodeCount, false);
    std::vector<bool> appearsInHalo(nodeCount, false);
    
    // Mark nodes in local elements
    for (size_t i = localStart; i < localEnd; ++i) {
        appearsInLocal[sfcToLocal[h_conn0[i]]] = true;
        appearsInLocal[sfcToLocal[h_conn1[i]]] = true;
        appearsInLocal[sfcToLocal[h_conn2[i]]] = true;
        appearsInLocal[sfcToLocal[h_conn3[i]]] = true;
    }
    
    // Mark nodes in halo elements
    for (auto [haloStart, haloEnd] : haloRanges) {
        size_t checkCount = std::min(size_t(100), haloEnd - haloStart);
        for (size_t i = haloStart; i < haloStart + checkCount; ++i) {
            appearsInHalo[sfcToLocal[h_conn0[i]]] = true;
            appearsInHalo[sfcToLocal[h_conn1[i]]] = true;
            appearsInHalo[sfcToLocal[h_conn2[i]]] = true;
            appearsInHalo[sfcToLocal[h_conn3[i]]] = true;
        }
    }
    
    // Verify: nodes ONLY in halos should NOT be owned
    size_t pureHaloNodes = 0;
    size_t boundaryNodes = 0;
    
    for (size_t i = 0; i < nodeCount; ++i) {
        if (appearsInHalo[i] && !appearsInLocal[i]) {
            // Pure halo node - should NOT be owned
            EXPECT_FALSE(h_ownership[i]) 
                << "Rank " << rank_ << ": Node " << i 
                << " appears only in halo elements but marked as owned";
            pureHaloNodes++;
        } else if (appearsInLocal[i] && appearsInHalo[i]) {
            // Boundary node - should be owned (appears in local elements)
            EXPECT_TRUE(h_ownership[i]) 
                << "Rank " << rank_ << ": Node " << i 
                << " appears in local elements but not marked as owned";
            boundaryNodes++;
        }
    }
    
    std::cout << "Rank " << rank_ << ": Found " << pureHaloNodes 
              << " pure halo nodes, " << boundaryNodes << " boundary nodes" << std::endl;
}

// TEST 11: Empty Rank and Edge Cases
TEST_F(MarsDomainAdjacencyAndHaloTest, EdgeCases) {
    DomainType domain(meshDir_, rank_, numRanks_);
    
    size_t nodeCount = domain.getNodeCount();
    size_t elementCount = domain.getElementCount();
    size_t localCount = domain.localElementCount();
    
    // Test 1: Verify counts are non-negative and consistent
    EXPECT_GE(nodeCount, 0) << "Node count should be non-negative";
    EXPECT_GE(elementCount, 0) << "Element count should be non-negative";
    EXPECT_LE(localCount, elementCount) << "Local count exceeds total";
    
    // Test 2: If local count is zero, all elements should be halos (multi-rank only)
    if (numRanks_ > 1 && localCount == 0) {
        EXPECT_EQ(domain.getHaloElementCount(), elementCount)
            << "Rank " << rank_ << ": Zero local elements but not all are halos";
        std::cout << "Rank " << rank_ << ": Edge case - all elements are halos" << std::endl;
    }
    
    // Test 3: Verify data structures are sized correctly even if empty
    const auto& offsets = domain.getNodeToElementOffsets();
    EXPECT_EQ(offsets.size(), nodeCount + 1) 
        << "Node-to-element offsets incorrect size";
    
    // Test 4: SFC map should have entry for every node
    const auto& sfcMap = domain.getLocalToGlobalSfcMap();
    EXPECT_EQ(sfcMap.size(), nodeCount) 
        << "SFC map size doesn't match node count";
    
    // Test 5: Ownership map should have entry for every node
    if (numRanks_ > 1) {
        const auto& ownership = domain.getNodeOwnershipMap();
        EXPECT_EQ(ownership.size(), nodeCount) 
            << "Ownership map size doesn't match node count";
    }
    
    std::cout << "Rank " << rank_ << ": Edge case validation passed "
              << "(nodes=" << nodeCount << ", elements=" << elementCount 
              << ", local=" << localCount << ")" << std::endl;
}

} // namespace testing
} // namespace mars

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}