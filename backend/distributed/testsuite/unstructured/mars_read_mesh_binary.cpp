#define USE_CUDA

#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>
#include "mars_read_mesh_binary.hpp"  // Include only this header

namespace fs = std::filesystem;

// Test fixture that creates test mesh files
class MeshReadBinaryTest : public ::testing::Test {
protected:
    // Test directory for mesh files
    fs::path testDir;
    int rank = 0;
    int numRanks = 1;
    
    // Create test files with known data
    void SetUp() override {
        // Create a temporary directory
        testDir = fs::temp_directory_path() / "mars_mesh_binary_test";
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
    
    // Clean up test files
    void TearDown() override {
        fs::remove_all(testDir);
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

// Test reading coordinates
TEST_F(MeshReadBinaryTest, ReadCoordinates) {
    // Test reading coordinates using the function directly
    auto [nodeCount, nodeStartIdx, x_data, y_data, z_data] = 
        mars::readMeshCoordinatesBinary(testDir.string(), rank, numRanks);
    
    // Verify node count
    EXPECT_EQ(nodeCount, 8);
    EXPECT_EQ(nodeStartIdx, 0);  // For rank 0
    
    // Verify coordinate values
    EXPECT_EQ(x_data.size(), 8);
    EXPECT_EQ(y_data.size(), 8);
    EXPECT_EQ(z_data.size(), 8);
    
    EXPECT_FLOAT_EQ(x_data[0], 1.0f);
    EXPECT_FLOAT_EQ(y_data[2], 0.3f);
    EXPECT_FLOAT_EQ(z_data[5], 60.0f);
}

// Test reading connectivity with tuple version
TEST_F(MeshReadBinaryTest, ReadConnectivityTuple) {
    // Test reading tetrahedral connectivity (4 nodes per element)
    auto [elementCount, conn_tuple] = 
        mars::readMeshConnectivityBinaryTuple<4>(testDir.string(), 0, rank, numRanks);
    
    // Verify element count
    EXPECT_EQ(elementCount, 4);
    
    // Verify tuple structure (should have 4 vectors)
    EXPECT_EQ(std::tuple_size_v<decltype(conn_tuple)>, 4);
    
    // Verify connectivity values
    EXPECT_EQ(std::get<0>(conn_tuple)[0], 0);
    EXPECT_EQ(std::get<1>(conn_tuple)[1], 3);
    EXPECT_EQ(std::get<2>(conn_tuple)[2], 6);
    EXPECT_EQ(std::get<3>(conn_tuple)[3], 1);
    
    // Test with triangles (3 nodes per element)
    // We'll use the same files but interpret first 3 indices as triangle connectivity
    auto [triCount, tri_tuple] = 
        mars::readMeshConnectivityBinaryTuple<3>(testDir.string(), 0, rank, numRanks);
    
    EXPECT_EQ(triCount, 4);
    EXPECT_EQ(std::tuple_size_v<decltype(tri_tuple)>, 3);
}

// Test reading connectivity with vector of vectors version
TEST_F(MeshReadBinaryTest, ReadConnectivityVectors) {
    // Test reading connectivity using the vector version
    auto [elementCount, indices] = 
        mars::readMeshConnectivityBinary(testDir.string(), 4, 0, rank, numRanks);
    
    // Verify element and node counts
    EXPECT_EQ(elementCount, 4);
    EXPECT_EQ(indices.size(), 4);  // 4 indices per tet
    EXPECT_EQ(indices[0].size(), 4);  // 4 elements
    
    // Verify connectivity values
    EXPECT_EQ(indices[0][0], 0);
    EXPECT_EQ(indices[1][1], 3);
    EXPECT_EQ(indices[2][2], 6);
    EXPECT_EQ(indices[3][3], 1);
}

// Test multi-rank distribution
TEST_F(MeshReadBinaryTest, MultiRankDistribution) {
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
    
    // Create larger connectivity files
    std::vector<int> i0(50), i1(50), i2(50), i3(50);
    for (int i = 0; i < 50; i++) {
        i0[i] = i;
        i1[i] = i + 1;
        i2[i] = i + 2;
        i3[i] = i + 3;
    }
    
    createConnectivityFile("i0.int32", i0);
    createConnectivityFile("i1.int32", i1);
    createConnectivityFile("i2.int32", i2);
    createConnectivityFile("i3.int32", i3);
    
    // Test rank 0 of 2
    {
        auto [nodeCount0, nodeStartIdx0, x0, y0, z0] = 
            mars::readMeshCoordinatesBinary(testDir.string(), 0, 2);
        
        EXPECT_EQ(nodeCount0, 50); // First half of nodes
        EXPECT_EQ(nodeStartIdx0, 0);
        EXPECT_FLOAT_EQ(x0[0], 0.0f);
        EXPECT_FLOAT_EQ(x0[49], 49.0f);
        
        auto [elemCount0, conn0] = 
            mars::readMeshConnectivityBinaryTuple<4>(testDir.string(), nodeStartIdx0, 0, 2);
        
        EXPECT_EQ(elemCount0, 25); // First half of elements
    }
    
    // Test rank 1 of 2
    {
        auto [nodeCount1, nodeStartIdx1, x1, y1, z1] = 
            mars::readMeshCoordinatesBinary(testDir.string(), 1, 2);
        
        EXPECT_EQ(nodeCount1, 50); // Second half of nodes
        EXPECT_EQ(nodeStartIdx1, 50);
        EXPECT_FLOAT_EQ(x1[0], 50.0f); // First node of rank 1
        
        auto [elemCount1, conn1] = 
            mars::readMeshConnectivityBinaryTuple<4>(testDir.string(), nodeStartIdx1, 1, 2);
        
        EXPECT_EQ(elemCount1, 25); // Second half of elements
        
        // Verify that indices are adjusted for local numbering
        // Element at rank 1 should have its indices adjusted by nodeStartIdx1
        EXPECT_EQ(std::get<0>(conn1)[0], i0[25] - nodeStartIdx1);
    }
}

// Test error handling for missing files
TEST_F(MeshReadBinaryTest, HandleMissingFiles) {
    fs::remove(testDir / "x.float32"); // Delete a coordinate file
    
    // Should throw an exception for missing file
    EXPECT_THROW({
        mars::readMeshCoordinatesBinary(testDir.string(), rank, numRanks);
    }, std::runtime_error);
    
    // Also test connectivity error
    fs::remove(testDir / "i0.int32");
    EXPECT_THROW({
        mars::readMeshConnectivityBinaryTuple<4>(testDir.string(), 0, rank, numRanks);
    }, std::runtime_error);
}

// Test to print mesh connectivity in a pretty format
TEST_F(MeshReadBinaryTest, PrintMeshConnectivity) {
    // Read coordinates and connectivity
    auto [nodeCount, nodeStartIdx, x_data, y_data, z_data] = 
        mars::readMeshCoordinatesBinary(testDir.string(), rank, numRanks);
    
    auto [elementCount, conn_tuple] = 
        mars::readMeshConnectivityBinaryTuple<4>(testDir.string(), 0, rank, numRanks);
    
    // Print mesh summary
    std::cout << "\n=== MESH SUMMARY ===" << std::endl;
    std::cout << "Nodes: " << nodeCount << ", Elements: " << elementCount << std::endl;
    
    // Print node coordinates
    std::cout << "\n=== NODE COORDINATES ===" << std::endl;
    std::cout << "Index  |     X     |     Y     |     Z     " << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    for (size_t i = 0; i < nodeCount; i++) {
        printf("%5zu | %9.4f | %9.4f | %9.4f\n", 
               i, x_data[i], y_data[i], z_data[i]);
    }
    
    // Extract connectivity arrays for easier access
    const auto& i0 = std::get<0>(conn_tuple);
    const auto& i1 = std::get<1>(conn_tuple);
    const auto& i2 = std::get<2>(conn_tuple);
    const auto& i3 = std::get<3>(conn_tuple);
    
    // Print element connectivity
    std::cout << "\n=== ELEMENT CONNECTIVITY ===" << std::endl;
    std::cout << "Element |  Node 0  |  Node 1  |  Node 2  |  Node 3  " << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    for (size_t e = 0; e < elementCount; e++) {
        printf("%7zu | %8d | %8d | %8d | %8d\n", 
               e, i0[e], i1[e], i2[e], i3[e]);
    }
    
    // Visualize the mesh structure (simplified ASCII representation)
    std::cout << "\n=== MESH VISUALIZATION ===" << std::endl;
    std::cout << "Each row represents an element with connections between its nodes:" << std::endl;
    
    for (size_t e = 0; e < elementCount; e++) {
        std::cout << "Element " << e << ": ";
        std::cout << i0[e] << " -- " << i1[e] << " -- " << i2[e] << " -- " << i3[e] << " -- " << i0[e];
        std::cout << " (tetrahedron)" << std::endl;
    }
    EXPECT_TRUE(elementCount > 0) << "Mesh connectivity printed for visual inspection"; 
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}