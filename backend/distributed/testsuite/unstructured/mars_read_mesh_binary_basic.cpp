#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>
#include <cmath>
#include "mars_read_mesh_binary.hpp"

namespace fs = std::filesystem;
using namespace mars;

// Test fixture that creates test mesh files
class MeshReadBinaryTest : public ::testing::Test {
protected:
    // Test directory for mesh files
    fs::path testDir;
    int rank = 0;
    int numRanks = 1;
    
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

// Test reading coordinates
TEST_F(MeshReadBinaryTest, ReadCoordinates) {
    // Test reading coordinates using the function directly
    auto [nodeCount, nodeStartIdx, x_data, y_data, z_data] = 
        mars::readMeshCoordinatesBinary<float>(testDir.string(), rank, numRanks);
    
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

// Test reading coordinates with double precision
TEST_F(MeshReadBinaryTest, ReadCoordinatesDouble) {
    // Convert test data to double precision
    std::vector<double> x_double = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    std::vector<double> y_double = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
    std::vector<double> z_double = {10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0};
    
    // Create double-precision binary files
    std::ofstream x_file((testDir / "x.double").string(), std::ios::binary);
    std::ofstream y_file((testDir / "y.double").string(), std::ios::binary);
    std::ofstream z_file((testDir / "z.double").string(), std::ios::binary);
    
    x_file.write(reinterpret_cast<const char*>(x_double.data()), x_double.size() * sizeof(double));
    y_file.write(reinterpret_cast<const char*>(y_double.data()), y_double.size() * sizeof(double));
    z_file.write(reinterpret_cast<const char*>(z_double.data()), z_double.size() * sizeof(double));

    x_file.close();
    y_file.close();
    z_file.close();
    
    // Test reading double-precision coordinates
    auto [nodeCount, nodeStartIdx, x_data, y_data, z_data] = 
        mars::readMeshCoordinatesBinary<double>(testDir.string(), rank, numRanks);
    
    // Verify node count
    EXPECT_EQ(nodeCount, 8);
    
    // Verify coordinate values with double precision
    EXPECT_DOUBLE_EQ(x_data[0], 1.0);
    EXPECT_DOUBLE_EQ(y_data[2], 0.3);
    EXPECT_DOUBLE_EQ(z_data[5], 60.0);
}

// Test error handling for missing files
TEST_F(MeshReadBinaryTest, HandleMissingFiles) {
    fs::remove(testDir / "x.float32"); // Delete a coordinate file
    
    // Should throw an exception for missing file
    EXPECT_THROW({
        mars::readMeshCoordinatesBinary<float>(testDir.string(), rank, numRanks);
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
        mars::readMeshCoordinatesBinary<float>(testDir.string(), rank, numRanks);
    
    auto [elementCount, conn_tuple] = 
        mars::readMeshConnectivityBinaryTuple<4>(testDir.string(), 0, rank, numRanks);
    
    // Print mesh summary
    std::cout << "\n=== MESH SUMMARY ===" << std::endl;
    std::cout << "Nodes: " << nodeCount << ", Elements: " << elementCount << std::endl;
    
    // Extract connectivity arrays for easier access
    const auto& i0 = std::get<0>(conn_tuple);
    const auto& i1 = std::get<1>(conn_tuple);
    const auto& i2 = std::get<2>(conn_tuple);
    const auto& i3 = std::get<3>(conn_tuple);
    
    // Verify something to make the test useful
    EXPECT_GT(nodeCount, 0) << "Node count should be positive";
    EXPECT_GT(elementCount, 0) << "Element count should be positive";
}

int main(int argc, char **argv) {
    mars::Env env(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}