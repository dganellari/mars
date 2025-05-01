#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <unordered_set>
#include "mars_read_mesh_binary.hpp"

namespace fs = std::filesystem;
using namespace mars;

// Test fixture that creates test mesh files
class MeshReadBinaryMPITest : public ::testing::Test {
protected:
    // Test directory for mesh files
    fs::path testDir;
    
    // Create test files with known data
    void SetUp() override {
        // Get MPI rank for unique directory naming
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        
        // Create a temporary directory with rank-specific name
        testDir = fs::temp_directory_path() / ("mars_mesh_binary_test_rank_" + std::to_string(rank));
        fs::create_directories(testDir);
        
        // Create coordinate files
        createCoordinateFile("x.float32", {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f});
        createCoordinateFile("y.float32", {0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f});
        createCoordinateFile("z.float32", {10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f});
        
        // Create connectivity files for tetrahedra
        createConnectivityFile("i0.int32", {0, 2, 4, 6});
        createConnectivityFile("i1.int32", {1, 3, 5, 7});
        createConnectivityFile("i2.int32", {2, 4, 6, 0});
        createConnectivityFile("i3.int32", {3, 5, 7, 1});
    }
    
    void TearDown() override {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        
        // Allow all ranks to synchronize before cleanup
        MPI_Barrier(MPI_COMM_WORLD);
        
        // Only rank 0 cleans up
        if (rank == 0) {
            fs::remove_all(testDir);
        }
    }
    
    // Helper to create a binary file with float data
    void createCoordinateFile(const std::string& filename, const std::vector<float>& data) {
        std::ofstream file((testDir / filename).string(), std::ios::binary);
        file.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(float));
        file.close();
    }
    
    // Helper to create a binary file with int data
    void createConnectivityFile(const std::string& filename, const std::vector<int>& data) {
        std::ofstream file((testDir / filename).string(), std::ios::binary);
        file.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(int));
        file.close();
    }
};

// Test multi-rank distribution with element-based partitioning
TEST_F(MeshReadBinaryMPITest, MultiRankDistribution) {
    // Create larger coordinate files for testing rank distribution
    std::vector<float> x_coords(100), y_coords(100), z_coords(100);
    for (int i = 0; i < 100; i++) {
        x_coords[i] = static_cast<float>(i);
        y_coords[i] = static_cast<float>(i) * 0.1f;
        z_coords[i] = static_cast<float>(i) * 10.0f;
    }
    
    createCoordinateFile("x.float32", x_coords);
    createCoordinateFile("y.float32", y_coords);
    createCoordinateFile("z.float32", z_coords);
    
    // Create connectivity files with specific patterns to test node sharing
    // Each element i uses nodes [2i, 2i+1, 2i+2, 2i+3]
    std::vector<int> i0(50), i1(50), i2(50), i3(50);
    for (int i = 0; i < 50; i++) {
        i0[i] = i * 2;       // 0, 2, 4, ...
        i1[i] = i * 2 + 1;   // 1, 3, 5, ...
        i2[i] = i * 2 + 2;   // 2, 4, 6, ...
        i3[i] = i * 2 + 3;   // 3, 5, 7, ...
    }
    
    createConnectivityFile("i0.int32", i0);
    createConnectivityFile("i1.int32", i1);
    createConnectivityFile("i2.int32", i2);
    createConnectivityFile("i3.int32", i3);
    
    // Test with various numbers of ranks
    std::vector<int> rankCounts = {2, 3, 4, 8};
    
    for (int numRanks : rankCounts) {
        std::cout << "\nTesting with " << numRanks << " ranks:" << std::endl;
        
        size_t totalElements = 0;
        size_t totalNodesAcrossRanks = 0;
        std::unordered_set<int> uniqueNodesAcrossRanks;
        
        // For each rank, read its portion of the mesh
        for (int rank = 0; rank < numRanks; rank++) {
            auto [nodeCount, elementCount, x, y, z, conn, localToGlobal] =
                mars::readMeshWithElementPartitioning<4, float>(testDir.string(), rank, numRanks);
            
            // Track total elements
            totalElements += elementCount;
            
            // Track unique nodes
            totalNodesAcrossRanks += nodeCount;
            
            // Expected elements per rank (roughly elementCount / numRanks)
            size_t expectedElements = 50 / numRanks;
            if (rank == numRanks - 1) {
                expectedElements += 50 % numRanks; // Last rank gets remainder
            }
            
            EXPECT_EQ(elementCount, expectedElements)
                << "Rank " << rank << " of " << numRanks << " has incorrect element count";
            
            // Verify connectivity is valid (indices within bounds)
            const auto& i0_conn = std::get<0>(conn);
            const auto& i1_conn = std::get<1>(conn);
            const auto& i2_conn = std::get<2>(conn);
            const auto& i3_conn = std::get<3>(conn);
            
            for (size_t i = 0; i < elementCount; i++) {
                EXPECT_LT(i0_conn[i], nodeCount);
                EXPECT_LT(i1_conn[i], nodeCount);
                EXPECT_LT(i2_conn[i], nodeCount);
                EXPECT_LT(i3_conn[i], nodeCount);
            }
            
            // Verify the first element's connectivity has been correctly mapped
            if (elementCount > 0) {
                size_t firstGlobalElementIdx = rank * (50 / numRanks);
                
                // The global nodes for the first element are [2*firstGlobalElementIdx, 2*firstGlobalElementIdx+1, ...]
                // But they should be mapped to local indices starting from 0
                EXPECT_EQ(i0_conn[0], 0);
                EXPECT_EQ(i1_conn[0], 1);
                
                // Verify coordinate mapping
                int globalFirstNode = i0[firstGlobalElementIdx];
                EXPECT_FLOAT_EQ(x[0], static_cast<float>(globalFirstNode));
            }
        }
        
        // All elements should be distributed
        EXPECT_EQ(totalElements, 50) 
            << "With " << numRanks << " ranks, total elements don't match expected";
    }
}

// Test explicitly for element-based partitioning with actual MPI ranks
TEST_F(MeshReadBinaryMPITest, ElementBasedPartitioning) {
    // Create a special test case where elements need nodes from other partitions
    // Element 0 uses nodes [0,1,2,3]
    // Element 1 uses nodes [2,3,4,5] - shares nodes with Element 0
    // Element 2 uses nodes [4,5,6,7] - shares nodes with Element 1
    // Element 3 uses nodes [6,7,0,1] - shares nodes with Element 2 and Element 0
    
    createConnectivityFile("i0.int32", {0, 2, 4, 6});
    createConnectivityFile("i1.int32", {1, 3, 5, 7});
    createConnectivityFile("i2.int32", {2, 4, 6, 0});
    createConnectivityFile("i3.int32", {3, 5, 7, 1});
    
    // Get actual MPI rank and size from Mars environment
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    std::cout << "Running on rank " << rank << " of " << size << std::endl;
    
    // Read this rank's portion of the mesh
    auto [nodeCount, elementCount, x, y, z, connectivity, localToGlobal] =
        mars::readMeshWithElementPartitioning<4, float>(testDir.string(), rank, size);
    
    // Gather element and node counts from all ranks for validation
    std::vector<size_t> elemCounts(size);
    std::vector<size_t> nodeCounts(size);
    
    MPI_Gather(&elementCount, 1, MPI_UNSIGNED_LONG, 
               elemCounts.data(), 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Gather(&nodeCount, 1, MPI_UNSIGNED_LONG, 
               nodeCounts.data(), 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    
    // Print distribution on rank 0
    if (rank == 0) {
        size_t totalElements = 0;
        for (int r = 0; r < size; r++) {
            std::cout << "Rank " << r << " got " << elemCounts[r] << " elements and "
                      << nodeCounts[r] << " nodes" << std::endl;
            totalElements += elemCounts[r];
        }
        
        // Verify total element count is correct
        EXPECT_EQ(totalElements, 4) << "Total number of elements should be 4";
    }
    
    // All ranks validate their own data
    const auto& i0 = std::get<0>(connectivity);
    const auto& i1 = std::get<1>(connectivity);
    const auto& i2 = std::get<2>(connectivity);
    const auto& i3 = std::get<3>(connectivity);
    
    // Most important validation: all indices are within bounds
    for (size_t i = 0; i < elementCount; i++) {
        EXPECT_GE(i0[i], 0) << "Rank " << rank << ", Element " << i << " has negative i0";
        EXPECT_LT(i0[i], nodeCount) << "Rank " << rank << ", Element " << i << " has i0 out of bounds";
        
        EXPECT_GE(i1[i], 0) << "Rank " << rank << ", Element " << i << " has negative i1";
        EXPECT_LT(i1[i], nodeCount) << "Rank " << rank << ", Element " << i << " has i1 out of bounds";
        
        EXPECT_GE(i2[i], 0) << "Rank " << rank << ", Element " << i << " has negative i2";
        EXPECT_LT(i2[i], nodeCount) << "Rank " << rank << ", Element " << i << " has i2 out of bounds";
        
        EXPECT_GE(i3[i], 0) << "Rank " << rank << ", Element " << i << " has negative i3";
        EXPECT_LT(i3[i], nodeCount) << "Rank " << rank << ", Element " << i << " has i3 out of bounds";
    }
    
    // Only check element-specific values if there are elements on this rank
    if (elementCount > 0) {
    // Create vectors of the original connectivity for this test
        std::vector<int> orig_i0 = {0, 2, 4, 6};
        std::vector<int> orig_i1 = {1, 3, 5, 7};
        std::vector<int> orig_i2 = {2, 4, 6, 0};
        std::vector<int> orig_i3 = {3, 5, 7, 1};
        
        // Verify that first node has a valid coordinate
        for (size_t i = 0; i < nodeCount; i++) {
            EXPECT_GE(x[i], 1.0f) << "Coordinate value should be valid";
            EXPECT_LE(x[i], 8.0f) << "Coordinate value should be valid";
        }
    }
    
    // Synchronize all processes before finishing the test
    MPI_Barrier(MPI_COMM_WORLD);
}

int main(int argc, char **argv) {
    mars::Env env(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}