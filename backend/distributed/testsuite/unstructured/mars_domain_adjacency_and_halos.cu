#include "gtest/gtest.h"
#include "backend/distributed/unstructured/domain.hpp"
#include <cuda_runtime.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

namespace mars {
namespace testing {

// Test fixture for domain adjacency and halo tests
class DomainAdjacencyTest : public ::testing::Test {
protected:
    using RealType = float;
    using KeyType = unsigned int;
    using AcceleratorTag = cstone::GpuTag;
    using DomainType = ElementDomain<TetTag, RealType, KeyType, AcceleratorTag>;

    void SetUp() override {
        // Initialize MPI if not already done
        int initialized;
        MPI_Initialized(&initialized);
        if (!initialized) {
            int argc = 0;
            char** argv = nullptr;
            MPI_Init(&argc, &argv);
        }
        
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
        MPI_Comm_size(MPI_COMM_WORLD, &numRanks_);
    }

    void TearDown() override {
        // Clean up CUDA
        cudaDeviceReset();
    }

    int rank_;
    int numRanks_;
};

// Test 1: Local-to-Global SFC Map Correctness
TEST_F(DomainAdjacencyTest, LocalToGlobalSfcMapUniqueness) {
    // Create domain from test mesh
    std::string meshDir = "path/to/test/mesh"; // Update with actual test mesh
    DomainType domain(meshDir, rank_, numRanks_);
    
    const auto& sfcMap = domain.getLocalToGlobalSfcMap();
    size_t nodeCount = domain.getNodeCount();
    
    ASSERT_EQ(sfcMap.size(), nodeCount) << "SFC map size mismatch";
    
    // Copy to host for verification
    thrust::host_vector<KeyType> h_sfcMap = sfcMap;
    
    // Check uniqueness (should be sorted and unique after buildAdjacency)
    for (size_t i = 1; i < h_sfcMap.size(); ++i) {
        EXPECT_LT(h_sfcMap[i-1], h_sfcMap[i]) 
            << "SFC keys not unique or sorted at index " << i;
    }
}

// Test 2: Element-to-Node Local ID Connectivity
TEST_F(DomainAdjacencyTest, ElementToNodeLocalIdMapping) {
    std::string meshDir = "path/to/test/mesh";
    DomainType domain(meshDir, rank_, numRanks_);
    
    const auto& connLocalIds = domain.getElementToNodeConnectivity();
    size_t elementCount = domain.getElementCount();
    size_t nodeCount = domain.getNodeCount();
    
    // Copy connectivity to host
    thrust::host_vector<KeyType> h_conn0 = std::get<0>(connLocalIds);
    thrust::host_vector<KeyType> h_conn1 = std::get<1>(connLocalIds);
    thrust::host_vector<KeyType> h_conn2 = std::get<2>(connLocalIds);
    thrust::host_vector<KeyType> h_conn3 = std::get<3>(connLocalIds);
    
    // Verify all local IDs are in valid range [0, nodeCount)
    for (size_t i = 0; i < elementCount; ++i) {
        EXPECT_LT(h_conn0[i], nodeCount) << "Invalid local ID in conn0 at element " << i;
        EXPECT_LT(h_conn1[i], nodeCount) << "Invalid local ID in conn1 at element " << i;
        EXPECT_LT(h_conn2[i], nodeCount) << "Invalid local ID in conn2 at element " << i;
        EXPECT_LT(h_conn3[i], nodeCount) << "Invalid local ID in conn3 at element " << i;
    }
}

// Test 3: Node-to-Element CSR Map Correctness
TEST_F(DomainAdjacencyTest, NodeToElementCsrStructure) {
    std::string meshDir = "path/to/test/mesh";
    DomainType domain(meshDir, rank_, numRanks_);
    
    const auto& offsets = domain.getNodeToElementOffsets();
    const auto& list = domain.getNodeToElementList();
    size_t nodeCount = domain.getNodeCount();
    size_t elementCount = domain.getElementCount();
    
    ASSERT_EQ(offsets.size(), nodeCount + 1) << "Offset array size incorrect";
    
    // Copy to host
    thrust::host_vector<KeyType> h_offsets = offsets;
    thrust::host_vector<KeyType> h_list = list;
    
    // Check CSR structure validity
    EXPECT_EQ(h_offsets[0], 0) << "CSR offset must start at 0";
    EXPECT_EQ(h_offsets[nodeCount], h_list.size()) 
        << "Last offset must equal list size";
    
    // Offsets must be non-decreasing
    for (size_t i = 1; i <= nodeCount; ++i) {
        EXPECT_GE(h_offsets[i], h_offsets[i-1]) 
            << "Offsets not monotonic at index " << i;
    }
    
    // All element IDs in list must be valid
    for (size_t i = 0; i < h_list.size(); ++i) {
        EXPECT_LT(h_list[i], elementCount) 
            << "Invalid element ID in node-to-element list at index " << i;
    }
    
    // Verify connectivity consistency: if node N is in element E,
    // then E should appear in node N's element list
    const auto& connLocalIds = domain.getElementToNodeConnectivity();
    thrust::host_vector<KeyType> h_conn0 = std::get<0>(connLocalIds);
    
    // Check a few random elements
    for (size_t elemIdx = 0; elemIdx < std::min(elementCount, size_t(100)); elemIdx += 10) {
        KeyType nodeId = h_conn0[elemIdx];
        KeyType start = h_offsets[nodeId];
        KeyType end = h_offsets[nodeId + 1];
        
        // Element should appear in this node's list
        bool found = false;
        for (KeyType j = start; j < end; ++j) {
            if (h_list[j] == elemIdx) {
                found = true;
                break;
            }
        }
        EXPECT_TRUE(found) << "Element " << elemIdx << " not in node " << nodeId << "'s element list";
    }
}

// Test 4: Halo Element Range Validity
TEST_F(DomainAdjacencyTest, HaloElementRanges) {
    if (numRanks_ < 2) {
        GTEST_SKIP() << "Halo test requires multiple ranks";
    }
    
    std::string meshDir = "path/to/test/mesh";
    DomainType domain(meshDir, rank_, numRanks_);
    
    auto [localStart, localEnd] = domain.localElementRange();
    auto haloRanges = domain.haloElementRanges();
    size_t totalElements = domain.getElementCount();
    
    // Local range must be valid
    EXPECT_LT(localStart, localEnd) << "Invalid local element range";
    EXPECT_LE(localEnd, totalElements) << "Local end exceeds total elements";
    
    // Halo ranges must not overlap with local range
    for (auto [haloStart, haloEnd] : haloRanges) {
        EXPECT_LT(haloStart, haloEnd) << "Invalid halo range";
        EXPECT_TRUE(haloEnd <= localStart || haloStart >= localEnd)
            << "Halo range overlaps with local range";
    }
    
    // Total elements = local + halos
    size_t haloCount = 0;
    for (auto [start, end] : haloRanges) {
        haloCount += (end - start);
    }
    EXPECT_EQ(localEnd - localStart + haloCount, totalElements)
        << "Local + halo count mismatch";
}

// Test 5: Halo Element Indices Contiguity
TEST_F(DomainAdjacencyTest, HaloElementIndicesContiguous) {
    if (numRanks_ < 2) {
        GTEST_SKIP() << "Halo test requires multiple ranks";
    }
    
    std::string meshDir = "path/to/test/mesh";
    DomainType domain(meshDir, rank_, numRanks_);
    
    const auto& haloIndices = domain.getHaloElementIndices();
    size_t haloCount = domain.getHaloElementCount();
    
    ASSERT_EQ(haloIndices.size(), haloCount) << "Halo indices size mismatch";
    
    // Copy to host
    thrust::host_vector<size_t> h_haloIndices = haloIndices;
    
    // Verify all halo indices are valid and are actually halos
    for (size_t i = 0; i < h_haloIndices.size(); ++i) {
        size_t elemIdx = h_haloIndices[i];
        EXPECT_LT(elemIdx, domain.getElementCount()) 
            << "Invalid halo element index at " << i;
        EXPECT_TRUE(domain.isHaloElement(elemIdx))
            << "Index " << elemIdx << " marked as halo but isHaloElement() returns false";
    }
}

// Test 6: Node Ownership Correctness
TEST_F(DomainAdjacencyTest, NodeOwnershipMapping) {
    if (numRanks_ < 2) {
        GTEST_SKIP() << "Node ownership test requires multiple ranks";
    }
    
    std::string meshDir = "path/to/test/mesh";
    DomainType domain(meshDir, rank_, numRanks_);
    
    const auto& ownership = domain.getNodeOwnershipMap();
    size_t nodeCount = domain.getNodeCount();
    
    ASSERT_EQ(ownership.size(), nodeCount) << "Ownership map size mismatch";
    
    // Copy to host
    thrust::host_vector<bool> h_ownership = ownership;
    
    // Verify ownership consistency with element ownership
    const auto& connLocalIds = domain.getElementToNodeConnectivity();
    thrust::host_vector<KeyType> h_conn0 = std::get<0>(connLocalIds);
    thrust::host_vector<KeyType> h_conn1 = std::get<1>(connLocalIds);
    thrust::host_vector<KeyType> h_conn2 = std::get<2>(connLocalIds);
    thrust::host_vector<KeyType> h_conn3 = std::get<3>(connLocalIds);
    
    auto [localStart, localEnd] = domain.localElementRange();
    
    // All nodes in local elements must be marked as owned
    for (size_t elemIdx = localStart; elemIdx < localEnd; ++elemIdx) {
        EXPECT_TRUE(h_ownership[h_conn0[elemIdx]]) 
            << "Node " << h_conn0[elemIdx] << " in local element not marked as owned";
        EXPECT_TRUE(h_ownership[h_conn1[elemIdx]]) 
            << "Node " << h_conn1[elemIdx] << " in local element not marked as owned";
        EXPECT_TRUE(h_ownership[h_conn2[elemIdx]]) 
            << "Node " << h_conn2[elemIdx] << " in local element not marked as owned";
        EXPECT_TRUE(h_ownership[h_conn3[elemIdx]]) 
            << "Node " << h_conn3[elemIdx] << " in local element not marked as owned";
    }
}

// Test 7: SFC Coordinate Decoding Accuracy
TEST_F(DomainAdjacencyTest, SfcCoordinateDecoding) {
    std::string meshDir = "path/to/test/mesh";
    DomainType domain(meshDir, rank_, numRanks_);
    
    const auto& sfcMap = domain.getLocalToGlobalSfcMap();
    auto box = domain.getBoundingBox();
    
    // Copy to host
    thrust::host_vector<KeyType> h_sfcKeys = sfcMap;
    
    // Decode a few SFC keys and check bounds
    for (size_t i = 0; i < std::min(h_sfcKeys.size(), size_t(100)); i += 10) {
        auto [x, y, z] = domain.sfcToPhysicalCoordinate(h_sfcKeys[i]);
        
        EXPECT_GE(x, box.xmin()) << "Decoded x below box min";
        EXPECT_LE(x, box.xmax()) << "Decoded x above box max";
        EXPECT_GE(y, box.ymin()) << "Decoded y below box min";
        EXPECT_LE(y, box.ymax()) << "Decoded y above box max";
        EXPECT_GE(z, box.zmin()) << "Decoded z below box min";
        EXPECT_LE(z, box.zmax()) << "Decoded z above box max";
    }
}

// Test 8: Element Ownership Flags
TEST_F(DomainAdjacencyTest, ElementOwnershipFlags) {
    std::string meshDir = "path/to/test/mesh";
    DomainType domain(meshDir, rank_, numRanks_);
    
    auto [localStart, localEnd] = domain.localElementRange();
    size_t totalElements = domain.getElementCount();
    
    // Check ownership flags for all elements
    for (size_t i = 0; i < totalElements; ++i) {
        bool isLocal = domain.isLocalElement(i);
        bool isHalo = domain.isHaloElement(i);
        
        // Element must be either local or halo, not both
        EXPECT_NE(isLocal, isHalo) << "Element " << i << " ownership ambiguous";
        
        // Verify against local range
        if (i >= localStart && i < localEnd) {
            EXPECT_TRUE(isLocal) << "Element " << i << " in local range but not flagged as local";
        } else {
            EXPECT_TRUE(isHalo) << "Element " << i << " outside local range but not flagged as halo";
        }
    }
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