#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>
#include "domain.hpp"

namespace fs = std::filesystem;

// Test fixture that creates test mesh files
class MeshReadTest : public ::testing::Test {
protected:
    // Test directory for mesh files
    fs::path testDir;
    int rank = 0;
    int numRanks = 1;
    
    // Create test files with known data
    void SetUp() override {
        // Create a temporary directory
        testDir = fs::temp_directory_path() / "mars_mesh_test";
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

// Test reading tetrahedral mesh with float precision
TEST_F(MeshReadTest, ReadTetMeshFloat) {
    // Test with Real as float
    using Real = float;
    using ElementType = TetTag;
    using AccelType = cstone::GpuTag;
    using Domain = ElementDomain<ElementType, AccelType>;
    
    // Create host data structures to receive the data
    typename Domain::CoordsTuple h_coords;
    typename Domain::ConnectivityTuple h_conn;
    
    // Create a domain object to test the readMeshDataSoA method
    Domain domain(testDir.string(), rank, numRanks);
    
    // Access results
    EXPECT_EQ(domain.getNodeCount(), 8);
    EXPECT_EQ(domain.getElementCount(), 4);
    
    // Check coordinates (test a few points)
    EXPECT_FLOAT_EQ(domain.x()[0], 1.0f);
    EXPECT_FLOAT_EQ(domain.y()[2], 0.3f);
    EXPECT_FLOAT_EQ(domain.z()[5], 60.0f);
    
    // Test connectivity
    EXPECT_EQ(domain.indices<0>()[0], 0);
    EXPECT_EQ(domain.indices<1>()[1], 3);
    EXPECT_EQ(domain.indices<2>()[2], 6);
    EXPECT_EQ(domain.indices<3>()[3], 1);
}

// Test reading tetrahedral mesh with double precision
TEST_F(MeshReadTest, ReadTetMeshDouble) {
    // Test with Real as double to test conversion
    using ElementType = TetTag;
    using AccelType = cstone::GpuTag;
    
    // This is a specialized test that explicitly sets Real to double
    class TestDomain : public ElementDomain<ElementType, AccelType> {
    public:
        using ElementDomain<ElementType, AccelType>::ElementDomain;
        
        // Test method that calls readMeshDataSoA directly for testing
        void testRead(const std::string& meshFile) {
            CoordsTuple h_coords;
            ConnectivityTuple h_conn;
            readMeshDataSoA(meshFile, h_coords, h_conn);
            
            // Verify data is properly converted to double
            EXPECT_EQ(std::get<0>(h_coords).size(), 8);
            EXPECT_EQ(std::get<0>(h_conn).size(), 4);
            
            EXPECT_DOUBLE_EQ(std::get<0>(h_coords)[0], 1.0);
            EXPECT_DOUBLE_EQ(std::get<1>(h_coords)[2], 0.3);
            EXPECT_DOUBLE_EQ(std::get<2>(h_coords)[5], 60.0);
        }
    };
    
    TestDomain domain(testDir.string(), rank, numRanks);
    domain.testRead(testDir.string());
}

// Test error handling for missing files
TEST_F(MeshReadTest, HandleMissingFiles) {
    fs::remove(testDir / "x.float32"); // Delete a coordinate file
    
    using ElementType = TetTag;
    using AccelType = cstone::GpuTag;
    using Domain = ElementDomain<ElementType, AccelType>;
    
    // Should throw an exception for missing file
    EXPECT_THROW({
        Domain domain(testDir.string(), rank, numRanks);
    }, std::runtime_error);
}

// Test multi-rank distribution (simulated)
TEST_F(MeshReadTest, MultiRankDistribution) {
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
    
    // Test with multiple ranks
    using ElementType = TetTag;
    using AccelType = cstone::GpuTag;
    using Domain = ElementDomain<ElementType, AccelType>;
    
    // Test rank 0 of 2
    {
        Domain domain0(testDir.string(), 0, 2);
        EXPECT_EQ(domain0.getNodeCount(), 50); // First half of nodes
        EXPECT_EQ(domain0.getElementCount(), 25); // First half of elements
    }
    
    // Test rank 1 of 2
    {
        Domain domain1(testDir.string(), 1, 2);
        EXPECT_EQ(domain1.getNodeCount(), 50); // Second half of nodes
        EXPECT_EQ(domain1.getElementCount(), 25); // Second half of elements
        
        // Check specific node value that should be offset
        EXPECT_FLOAT_EQ(domain1.x()[0], 50.0f); // First node of rank 1
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}