#include "gtest/gtest.h"
#include "backend/distributed/unstructured/domain.hpp"
#include <cuda_runtime.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <mpi.h>

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
        
        // Use the test mesh from your existing tests
        meshDir_ = std::string(MARS_TEST_MESH_DIR) + "/cube_10k/";
    }

    void TearDown() override {
        cudaDeviceReset();
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
    thrust::host_vector<KeyType> h_sfcMap = sfcMap;
    
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
    
    // Get local ID connectivity (after buildAdjacency)
    auto conn0 = domain.getConnectivity<0>();
    auto conn1 = domain.getConnectivity<1>();
    auto conn2 = domain.getConnectivity<2>();
    auto conn3 = domain.getConnectivity<3>();
    
    // Copy to host
    thrust::host_vector<KeyType> h_conn0(conn0.begin(), conn0.end());
    thrust::host_vector<KeyType> h_conn1(conn1.begin(), conn1.end());
    thrust::host_vector<KeyType> h_conn2(conn2.begin(), conn2.end());
    thrust::host_vector<KeyType> h_conn3(conn3.begin(), conn3.end());
    
    // Verify all local IDs are in valid range [0, nodeCount)
    for (size_t i = 0; i < elementCount; ++i) {
        EXPECT_LT(h_conn0[i], nodeCount) 
            << "Rank " << rank_ << ": Invalid local ID in conn0 at element " << i;
        EXPECT_LT(h_conn1[i], nodeCount) 
            << "Rank " << rank_ << ": Invalid local ID in conn1 at element " << i;
        EXPECT_LT(h_conn2[i], nodeCount) 
            << "Rank " << rank_ << ": Invalid local ID in conn2 at element " << i;
        EXPECT_LT(h_conn3[i], nodeCount) 
            << "Rank " << rank_ << ": Invalid local ID in conn3 at element " << i;
    }
    
    std::cout << "Rank " << rank_ << ": Verified connectivity for " 
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
    thrust::host_vector<KeyType> h_offsets(offsets.begin(), offsets.end());
    thrust::host_vector<KeyType> h_list(list.begin(), list.end());
    
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
    thrust::host_vector<size_t> h_haloIndices(haloIndices.begin(), haloIndices.end());
    
    // Verify all halo indices are valid and are actually halos
    for (size_t i = 0; i < h_haloIndices.size(); ++i) {
        size_t elemIdx = h_haloIndices[i];
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
    
    // Copy to host
    thrust::host_vector<bool> h_ownership(ownership.begin(), ownership.end());
    
    // Get connectivity for verification
    auto conn0 = domain.getConnectivity<0>();
    auto conn1 = domain.getConnectivity<1>();
    auto conn2 = domain.getConnectivity<2>();
    auto conn3 = domain.getConnectivity<3>();
    
    thrust::host_vector<KeyType> h_conn0(conn0.begin(), conn0.end());
    thrust::host_vector<KeyType> h_conn1(conn1.begin(), conn1.end());
    thrust::host_vector<KeyType> h_conn2(conn2.begin(), conn2.end());
    thrust::host_vector<KeyType> h_conn3(conn3.begin(), conn3.end());
    
    auto [localStart, localEnd] = domain.localElementRange();
    
    // Sample check: nodes in local elements should be owned
    size_t checkCount = std::min(size_t(100), localEnd - localStart);
    for (size_t i = 0; i < checkCount; ++i) {
        size_t elemIdx = localStart + i;
        EXPECT_TRUE(h_ownership[h_conn0[elemIdx]]) 
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
    thrust::host_vector<KeyType> h_sfcKeys(sfcMap.begin(), sfcMap.end());
    
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

} // namespace testing
} // namespace mars

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}