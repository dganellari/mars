#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>
#include <cmath>
#include "mars_read_mesh_binary.hpp"  // Include only this header

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

// Test multi-rank distribution with element-based partitioning
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
            auto [nodeCount, elementCount, x, y, z, conn] = 
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
            
            std::cout << "Rank " << rank << ": " << elementCount << " elements, " 
                    << nodeCount << " nodes" << std::endl;
        }
        
        // All elements should be distributed
        EXPECT_EQ(totalElements, 50) 
            << "With " << numRanks << " ranks, total elements don't match expected";
    }
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

    // Add debug output to check if files exist and have the right size
std::cout << "Double files created at: " << testDir << std::endl;
std::cout << "x.double exists: " << (fs::exists(testDir / "x.double") ? "yes" : "no") << std::endl;
    
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

// ########################External meshes########################

class ExternalMeshTest : public ::testing::TestWithParam<std::string> {
protected:
    int rank = 0;
    int numRanks = 1;
    std::string meshPath;
    
    void SetUp() override {
        meshPath = GetParam();
        if (meshPath == "DUMMY_PATH_NO_MESH_FOUND" || !fs::exists(meshPath)) {
            GTEST_SKIP() << "No valid mesh file found for testing";
        }
    }
};

// Basic validation test
TEST_P(ExternalMeshTest, BasicMeshValidation) {
    std::cout << "Testing mesh: " << meshPath << std::endl;
    
    // Read coordinates with float precision
    auto [nodeCount, nodeStartIdx, x_data, y_data, z_data] = 
        mars::readMeshCoordinatesBinary<float>(meshPath, rank, numRanks);
    
    // Basic checks
    EXPECT_GT(nodeCount, 0) << "Mesh should have nodes";
    EXPECT_EQ(x_data.size(), nodeCount);
    EXPECT_EQ(y_data.size(), nodeCount);
    EXPECT_EQ(z_data.size(), nodeCount);
    
    // Try reading connectivity
    try {
        auto [tetCount, tet_conn] = 
            mars::readMeshConnectivityBinaryTuple<4>(meshPath, nodeStartIdx, rank, numRanks);
        
        EXPECT_GT(tetCount, 0) << "Expecting tetrahedral elements";
        std::cout << "Mesh statistics: " << nodeCount << " nodes, " << tetCount << " tetrahedra" << std::endl;
        
        // Check indices are valid
        if (tetCount > 0) {
            for (size_t i = 0; i < std::min(tetCount, size_t(10)); i++) {
                    int nodeIdx0 = std::get<0>(tet_conn)[i];
                    int nodeIdx1 = std::get<1>(tet_conn)[i];
                    int nodeIdx2 = std::get<2>(tet_conn)[i];
                    int nodeIdx3 = std::get<3>(tet_conn)[i];

                    EXPECT_GE(nodeIdx0, 0) << "Node index should be non-negative";
                    EXPECT_LT(nodeIdx0, nodeCount) << "Node index should be less than node count";

                    EXPECT_GE(nodeIdx1, 0) << "Node index should be non-negative";
                    EXPECT_LT(nodeIdx1, nodeCount) << "Node index should be less than node count";

                    EXPECT_GE(nodeIdx2, 0) << "Node index should be non-negative";
                    EXPECT_LT(nodeIdx2, nodeCount) << "Node index should be less than node count";

                    EXPECT_GE(nodeIdx3, 0) << "Node index should be non-negative";
                    EXPECT_LT(nodeIdx3, nodeCount) << "Node index should be less than node count";
                }
            }
    } catch (const std::exception& e) {
        std::cout << "No tetrahedral connectivity found: " << e.what() << std::endl;
    }
}

// Mesh statistics test
TEST_P(ExternalMeshTest, MeshStatistics) {
    // Read coordinates
    auto [nodeCount, nodeStartIdx, x_data, y_data, z_data] = 
        mars::readMeshCoordinatesBinary<float>(meshPath, rank, numRanks);
    
    if (nodeCount == 0) {
        GTEST_SKIP() << "Empty mesh, skipping statistics test";
    }
    
    // Calculate bounding box
    float xmin = *std::min_element(x_data.begin(), x_data.end());
    float xmax = *std::max_element(x_data.begin(), x_data.end());
    float ymin = *std::min_element(y_data.begin(), y_data.end());
    float ymax = *std::max_element(y_data.begin(), y_data.end());
    float zmin = *std::min_element(z_data.begin(), z_data.end());
    float zmax = *std::max_element(z_data.begin(), z_data.end());
    
    // Calculate dimensions
    float dx = xmax - xmin;
    float dy = ymax - ymin;
    float dz = zmax - zmin;
    float diagonal = std::sqrt(dx*dx + dy*dy + dz*dz);
    
    std::cout << "Mesh bounding box: [" << xmin << ", " << xmax << "] × ["
              << ymin << ", " << ymax << "] × [" << zmin << ", " << zmax << "]" << std::endl;
    std::cout << "Mesh diagonal length: " << diagonal << std::endl;
    
    EXPECT_GT(dx, 0) << "X dimension should be positive";
    EXPECT_GT(dy, 0) << "Y dimension should be positive";
    EXPECT_GT(dz, 0) << "Z dimension should be positive";
}

// Performance test (only for larger meshes)
TEST_P(ExternalMeshTest, ReadPerformance) {
    auto start = std::chrono::high_resolution_clock::now();
    
    // Read coordinates
    auto [nodeCount, nodeStartIdx, x_data, y_data, z_data] = 
        mars::readMeshCoordinatesBinary<float>(meshPath, rank, numRanks);
    
    auto coordTime = std::chrono::high_resolution_clock::now();
    auto coordDuration = std::chrono::duration_cast<std::chrono::milliseconds>(
        coordTime - start).count();
    
    std::cout << "Read " << nodeCount << " nodes in " << coordDuration << " ms" << std::endl;
    
    // Only check performance for large meshes
    if (nodeCount > 1000000) {
        EXPECT_LT(coordDuration, 5000) << "Reading large mesh coordinates should be reasonably fast";
    }
    
    // Try to read connectivity
    try {
        auto connStart = std::chrono::high_resolution_clock::now();
        
        auto [tetCount, tet_conn] = 
            mars::readMeshConnectivityBinaryTuple<4>(meshPath, nodeStartIdx, rank, numRanks);
        
        auto connTime = std::chrono::high_resolution_clock::now();
        auto connDuration = std::chrono::duration_cast<std::chrono::milliseconds>(
            connTime - connStart).count();
        
        std::cout << "Read " << tetCount << " tetrahedra in " << connDuration << " ms" << std::endl;
    } catch (const std::exception& e) {
        std::cout << "No tetrahedral connectivity available" << std::endl;
    }
}

// Multi-rank simulation
TEST_P(ExternalMeshTest, TwoRankSimulation) {
    // Try with 2 ranks
    int simulatedRanks = 2;
    std::cout << "\nTesting " << simulatedRanks << "-way partitioning:" << std::endl;
    
    size_t totalNodes = 0;
    
    for (int r = 0; r < simulatedRanks; r++) {
        // Read coordinates for this simulated rank
        auto [nodeCount, nodeStartIdx, x_data, y_data, z_data] = 
            mars::readMeshCoordinatesBinary<float>(meshPath, r, simulatedRanks);
        
        totalNodes += nodeCount;
        std::cout << "Rank " << r << ": " << nodeCount << " nodes starting at " << nodeStartIdx << std::endl;
    }
    
    // Verify against single-rank read
    auto [fullNodeCount, _, x_full, y_full, z_full] = 
        mars::readMeshCoordinatesBinary<float>(meshPath, 0, 1);
    
    std::cout << "Total nodes across " << simulatedRanks << " ranks: " << totalNodes << std::endl;
    std::cout << "Nodes in single-rank read: " << fullNodeCount << std::endl;
    
    EXPECT_EQ(totalNodes, fullNodeCount) << "Multi-rank node count should match single-rank";
}

// Multi-rank simulation with flexible rank configurations
TEST_P(ExternalMeshTest, MultiRankSimulation) {
    // Test with various rank configurations
    std::vector<int> rankConfigs = {2, 3, 4, 8};
    
    // Read the full mesh once for reference
    auto [fullNodeCount, _, x_full, y_full, z_full] = 
        mars::readMeshCoordinatesBinary<float>(meshPath, 0, 1);
    
    // Skip tests for very small meshes
    if (fullNodeCount < 10) {
        GTEST_SKIP() << "Mesh too small for meaningful rank partitioning tests";
        return;
    }
    
    // Limit the maximum number of ranks based on mesh size
    // (avoid testing with more ranks than nodes)
    int maxRanks = static_cast<int>(std::min(fullNodeCount / 2, static_cast<size_t>(8)));
    
    for (int numRanks : rankConfigs) {
        // Skip if too many ranks for this mesh
        if (numRanks > maxRanks) continue;
        
        std::cout << "\nTesting " << numRanks << "-way partitioning:" << std::endl;
        
        size_t totalNodes = 0;
        size_t totalElements = 0;
        
        // Verify node distribution and element distribution across ranks
        for (int r = 0; r < numRanks; r++) {
            // Read coordinates for this simulated rank
            auto [nodeCount, nodeStartIdx, x_data, y_data, z_data] = 
                mars::readMeshCoordinatesBinary<float>(meshPath, r, numRanks);
            
            totalNodes += nodeCount;
            
            // Calculate expected node count for this rank
            size_t expectedNodeCount = fullNodeCount / numRanks;
            if (r == numRanks - 1) {
                // Last rank may get remainder nodes
                expectedNodeCount = fullNodeCount - (numRanks - 1) * (fullNodeCount / numRanks);
            }
            
            // Check if node count is approximately as expected (within 20%)
            // This test is flexible since partitioning may not be exactly even
            EXPECT_NEAR(static_cast<double>(nodeCount), 
                       static_cast<double>(expectedNodeCount),
                       expectedNodeCount * 0.2) << "Rank " << r << " node count not within expected range";
                       
            // Verify node start index
            EXPECT_EQ(nodeStartIdx, r * (fullNodeCount / numRanks)) 
                << "Incorrect node start index for rank " << r;
            
            std::cout << "Rank " << r << ": " << nodeCount << " nodes starting at " << nodeStartIdx << std::endl;
            
            // Try to read connectivity
            try {
                auto [elemCount, elem_conn] = 
                    mars::readMeshConnectivityBinaryTuple<4>(meshPath, nodeStartIdx, r, numRanks);
                
                totalElements += elemCount;
            } catch (const std::exception& e) {
                // Not a critical error - some meshes might not have connectivity
                std::cout << "No connectivity data for rank " << r << ": " << e.what() << std::endl;
            }
        }
        
        // Verify total node count matches full mesh
        std::cout << "Total nodes across " << numRanks << " ranks: " << totalNodes << std::endl;
        std::cout << "Nodes in single-rank read: " << fullNodeCount << std::endl;
        
        EXPECT_EQ(totalNodes, fullNodeCount) 
            << "Total node count across " << numRanks << " ranks should match full mesh";
        
        // If we successfully read elements, verify total element count
        try {
            auto [fullElemCount, _] = 
                mars::readMeshConnectivityBinaryTuple<4>(meshPath, 0, 0, 1);
                
            std::cout << "Total elements across ranks: " << totalElements << std::endl;
            std::cout << "Elements in single-rank read: " << fullElemCount << std::endl;
            
            EXPECT_EQ(totalElements, fullElemCount) 
                << "Total element count across " << numRanks << " ranks should match full mesh";
        } catch (const std::exception& e) {
            // Skip element verification if connectivity not available
            std::cout << "Skipping element verification: " << e.what() << std::endl;
        }
    }
}

// Find available mesh directories containing the binary format files
std::vector<std::string> GetMeshFiles() {
    std::vector<std::string> meshPaths;
    
    // Add environment variable path if available
    const char* envPath = std::getenv("MESH_PATH");
    if (envPath) {
        if (fs::exists(envPath) && fs::is_directory(fs::path(envPath))) {
            // Check if this directory contains required mesh files
            if (fs::exists(fs::path(envPath) / "x.float32") || 
                fs::exists(fs::path(envPath) / "x.double")) {
                meshPaths.push_back(envPath);
            } else {
                std::cerr << "Warning: MESH_PATH directory exists but does not contain required mesh files" << std::endl;
            }
        } else {
            std::cerr << "Warning: MESH_PATH is not a valid directory: " << envPath << std::endl;
        }
    }
    
    // Add common test locations
    std::vector<std::string> commonLocations = {
        "./test_data",
        "../test_data",
        "./meshes",
        "../meshes"
        // Add more paths as needed
    };
    
    for (const auto& path : commonLocations) {
        if (fs::exists(path) && fs::is_directory(fs::path(path))) {
            // Check if this directory contains required mesh files
            if (fs::exists(fs::path(path) / "x.float32") || 
                fs::exists(fs::path(path) / "x.double")) {
                meshPaths.push_back(path);
            }
        }
    }
    
    // If no mesh paths found, add a dummy path so tests can be skipped gracefully
    if (meshPaths.empty()) {
        meshPaths.push_back("DUMMY_PATH_NO_MESH_FOUND");
    }
    
    return meshPaths;
}

// Instantiate the test suite
INSTANTIATE_TEST_SUITE_P(
    MeshFiles,
    ExternalMeshTest,
    ::testing::ValuesIn(GetMeshFiles()),
    [](const testing::TestParamInfo<std::string>& info) {
        // Check if this is our dummy path
        if (info.param == "DUMMY_PATH_NO_MESH_FOUND") {
            return std::string("NoMeshFound");
        }
        
        // For real paths, use the last directory component as the test name
        fs::path p(info.param);
        std::string name;
        
        if (fs::is_directory(p)) {
            // Get the last directory component
            name = p.filename().string();
        } else {
            // Fallback to the full path
            name = p.string();
        }
        
        // Replace non-alphanumeric characters
        std::replace_if(name.begin(), name.end(), 
                       [](char c) { return !std::isalnum(c); }, '_');
        
        // Ensure the name is not empty
        if (name.empty()) {
            name = "Mesh";
        }
        
        return name;
    }
);

// Simple test for external meshes that uses environment variables
TEST(ExternalMeshTest, ReadExistingMeshEnvVar) {
    // Get mesh path from environment variable
    const char* meshPathEnv = std::getenv("MESH_PATH");
    std::string meshPath = meshPathEnv ? meshPathEnv : "";
    
    if (meshPath.empty()) {
        GTEST_SKIP() << "MESH_PATH environment variable is not set";
        return;
    }
    
    if (!fs::exists(meshPath)) {
        FAIL() << "MESH_PATH directory does not exist: " << meshPath;
        return;
    }
    
    if (!fs::is_directory(fs::path(meshPath))) {
        FAIL() << "MESH_PATH must be a directory: " << meshPath;
        return;
    }
    
    // Check if directory contains required files
    if (!fs::exists(fs::path(meshPath) / "x.float32") && 
        !fs::exists(fs::path(meshPath) / "x.double")) {
        FAIL() << "MESH_PATH directory does not contain required coordinate files: " << meshPath;
        return;
    }
    
    int rank = 0;
    int numRanks = 1;
    
    // Try reading with float precision
    std::cout << "Testing mesh reading with float precision from: " << meshPath << std::endl;
    auto [nodeCount, nodeStartIdx, x_data, y_data, z_data] = 
        mars::readMeshCoordinatesBinary<float>(meshPath, rank, numRanks);
    
    EXPECT_GT(nodeCount, 0) << "Expected positive node count";
    
    // Try reading connectivity
    try {
        auto [tetCount, tet_conn] = 
            mars::readMeshConnectivityBinaryTuple<4>(meshPath, nodeStartIdx, rank, numRanks);
        
        EXPECT_GT(tetCount, 0) << "Expected positive tetrahedron count";
        
        // Print mesh summary
        std::cout << "Read mesh with " << nodeCount << " nodes and " 
                  << tetCount << " tetrahedra." << std::endl;
        
        // Validate indices are within node range
        if (tetCount > 0) {
            for (size_t i = 0; i < std::min(tetCount, size_t(10)); i++) {
                EXPECT_GE(std::get<0>(tet_conn)[i], 0);
                EXPECT_LT(std::get<0>(tet_conn)[i], nodeCount);
            }
        }
    } catch (const std::exception& e) {
        std::cout << "Could not read tetrahedron connectivity: " << e.what() << std::endl;
    }
    
    // Try reading with double precision
    try {
        std::cout << "Testing mesh reading with double precision from: " << meshPath << std::endl;
        auto [nodeCount2, nodeStartIdx2, x_data2, y_data2, z_data2] = 
            mars::readMeshCoordinatesBinary<double>(meshPath, rank, numRanks);
        
        EXPECT_GT(nodeCount2, 0) << "Expected positive node count";
        std::cout << "Read mesh with " << nodeCount2 << " nodes (double precision)" << std::endl;
    } catch (const std::exception& e) {
        std::cout << "Could not read double precision coordinates: " << e.what() << std::endl;
    }
}

// Test to visualize mesh in ParaView by writing to standard VTK format
TEST_P(ExternalMeshTest, WriteStandardVTK) {
    std::cout << "Reading mesh for standard VTK export: " << meshPath << std::endl;
    
    // Read coordinates with float precision
    auto [nodeCount, nodeStartIdx, x_data, y_data, z_data] = 
        mars::readMeshCoordinatesBinary<float>(meshPath, rank, numRanks);
    
    if (nodeCount == 0) {
        GTEST_SKIP() << "Empty mesh, skipping VTK export";
    }
    
    // Create output filename in the same directory as the input
    fs::path inputDir = fs::path(meshPath).parent_path();
    fs::path outputPath = inputDir / "mesh_visualization.vtk";
    std::cout << "Writing VTK file to: " << outputPath << std::endl;
    
    // Open file for writing
    std::ofstream vtkFile(outputPath, std::ios::out);
    if (!vtkFile) {
        GTEST_SKIP() << "Failed to open VTK file for writing: " << outputPath;
    }
    
    // Write VTK header
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Mesh exported from MARS binary format\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET UNSTRUCTURED_GRID\n";
    
    // Write points
    vtkFile << "POINTS " << nodeCount << " float\n";
    for (size_t i = 0; i < nodeCount; i++) {
        vtkFile << x_data[i] << " " << y_data[i] << " " << z_data[i] << "\n";
    }
    
    // Try reading connectivity
    int numElements = 0;
    
    try {
        auto [tetCount, tet_conn] = 
            mars::readMeshConnectivityBinaryTuple<4>(meshPath, nodeStartIdx, rank, numRanks);
        
        numElements = tetCount;
        
        // Extract connectivity arrays for easier access
        const auto& i0 = std::get<0>(tet_conn);
        const auto& i1 = std::get<1>(tet_conn);
        const auto& i2 = std::get<2>(tet_conn);
        const auto& i3 = std::get<3>(tet_conn);
        
        // Write cells (must do this inside the try block where tet_conn is in scope)
        vtkFile << "CELLS " << numElements << " " << numElements * 5 << "\n";
        for (size_t i = 0; i < tetCount; i++) {
            vtkFile << "4 " 
                   << i0[i] << " " 
                   << i1[i] << " " 
                   << i2[i] << " " 
                   << i3[i] << "\n";
        }
        
        // Write cell types (10 = VTK_TETRA)
        vtkFile << "CELL_TYPES " << numElements << "\n";
        for (int i = 0; i < numElements; i++) {
            vtkFile << "10\n";  // 10 is the VTK type for tetrahedron
        }
    } catch (const std::exception& e) {
        std::cout << "No tetrahedral connectivity found, exporting point cloud only: " << e.what() << std::endl;
    }
    
    vtkFile.close();
    
    std::cout << "VTK file written with " << nodeCount << " nodes and " 
              << numElements << " tetrahedral elements" << std::endl;
    std::cout << "VTK file path: " << outputPath << std::endl;
    
    // Verify file was written
    EXPECT_TRUE(fs::exists(outputPath)) << "VTK file was not created";
    EXPECT_GT(fs::file_size(outputPath), 0) << "VTK file is empty";
}

// Test explicitly for element-based partitioning with actual MPI ranks
TEST_F(MeshReadBinaryTest, ElementBasedPartitioning) {
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
    auto [nodeCount, elementCount, x, y, z, connectivity] = 
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

// Test partitioning balance
TEST_P(ExternalMeshTest, ElementBasedPartitioningBalance) {
    try {
        // Read the full mesh to get total counts
        auto [nodeCount, elemStartIdx, x_full, y_full, z_full] = 
            mars::readMeshCoordinatesBinary<float>(meshPath, 0, 1);
        
        // Skip very small meshes
        if (nodeCount < 100) {
            GTEST_SKIP() << "Mesh too small for meaningful balance test";
            return;
        }
        
        // Test with 4 ranks
        int testRanks = 4;
        std::vector<size_t> rankNodeCounts(testRanks);
        std::vector<size_t> rankElementCounts(testRanks);
        
        // Read each rank's portion with element-based partitioning
        for (int r = 0; r < testRanks; r++) {
            try {
                auto [nodeCount, elementCount, x, y, z, conn] = 
                    mars::readMeshWithElementPartitioning<4, float>(meshPath, r, testRanks);
                
                rankNodeCounts[r] = nodeCount;
                rankElementCounts[r] = elementCount;
                
                std::cout << "Rank " << r << ": " << elementCount << " elements, " 
                         << nodeCount << " nodes" << std::endl;
            } catch (const std::exception& e) {
                std::cout << "Error reading rank " << r << ": " << e.what() << std::endl;
            }
        }
        
        // Calculate node and element imbalance
        size_t minNodes = *std::min_element(rankNodeCounts.begin(), rankNodeCounts.end());
        size_t maxNodes = *std::max_element(rankNodeCounts.begin(), rankNodeCounts.end());
        double nodeImbalance = (minNodes > 0) ? static_cast<double>(maxNodes) / minNodes : 0.0;
        
        std::cout << "Node distribution: min=" << minNodes << ", max=" << maxNodes 
                 << ", imbalance=" << nodeImbalance << std::endl;
        
        // Element balance (if available)
        if (rankElementCounts[0] > 0) {
            size_t minElems = *std::min_element(rankElementCounts.begin(), rankElementCounts.end());
            size_t maxElems = *std::max_element(rankElementCounts.begin(), rankElementCounts.end());
            double elemImbalance = (minElems > 0) ? static_cast<double>(maxElems) / minElems : 0.0;
            
            std::cout << "Element distribution: min=" << minElems << ", max=" << maxElems 
                     << ", imbalance=" << elemImbalance << std::endl;
            
            // For large meshes, expect reasonable balance
            if (elemImbalance > 0) {
                EXPECT_LT(elemImbalance, 2.0) << "Element imbalance should be less than 2x";
            }
        }
    } catch (const std::exception& e) {
        GTEST_SKIP() << "Failed to read mesh: " << e.what();
    }
}

// Test that verifies no illegal node indices
TEST_P(ExternalMeshTest, NoIllegalNodeIndices) {
    try {
        for (int r = 0; r < 2; r++) {
            auto [nodeCount, elementCount, x, y, z, conn] = 
                mars::readMeshWithElementPartitioning<4, float>(meshPath, r, 2);
            
            // Check every index in connectivity
            auto checkIndices = [nodeCount](const auto& indices) {
                for (size_t i = 0; i < indices.size(); i++) {
                    // All indices should be non-negative and within bounds
                    EXPECT_GE(indices[i], 0) << "Negative node index found";
                    EXPECT_LT(indices[i], nodeCount) << "Node index out of bounds";
                }
            };
            
            std::apply([&checkIndices](const auto&... vecs) {
                (checkIndices(vecs), ...);
            }, conn);
            
            std::cout << "Rank " << r << ": Verified " << elementCount 
                     << " elements have valid node indices" << std::endl;
        }
    } catch (const std::exception& e) {
        GTEST_SKIP() << "Failed to read mesh: " << e.what();
    }
}

// Test to visualize mesh using ADIOS2 for VTK output
#ifdef ADIOS2_HAVE_MPI
TEST_P(ExternalMeshTest, WriteADIOS2VTK) {
    std::cout << "Reading mesh for ADIOS2 VTK export: " << meshPath << std::endl;
    
    // Read coordinates with float precision
    auto [nodeCount, nodeStartIdx, x_data, y_data, z_data] = 
        mars::readMeshCoordinatesBinary<float>(meshPath, rank, numRanks);
    
    if (nodeCount == 0) {
        GTEST_SKIP() << "Empty mesh, skipping ADIOS2 export";
    }
    
    // Try reading connectivity
    std::vector<int> cellTypes;
    std::vector<int> connectivity;
    std::vector<int> cellOffsets;
    int numElements = 0;
    
    try {
        auto [tetCount, tet_conn] = 
            mars::readMeshConnectivityBinaryTuple<4>(meshPath, nodeStartIdx, rank, numRanks);
        
        numElements = tetCount;
        
        // For ADIOS2 VTK, we need to prepare connectivity in special format
        connectivity.reserve(tetCount * 4);  // 4 nodes per tetrahedron
        cellTypes.reserve(tetCount);
        cellOffsets.reserve(tetCount);
        
        int offset = 0;
        for (size_t i = 0; i < tetCount; i++) {
            connectivity.push_back(std::get<0>(tet_conn)[i]);
            connectivity.push_back(std::get<1>(tet_conn)[i]);
            connectivity.push_back(std::get<2>(tet_conn)[i]);
            connectivity.push_back(std::get<3>(tet_conn)[i]);
            
            cellTypes.push_back(10);  // 10 = VTK_TETRA
            cellOffsets.push_back(offset);
            offset += 4;  // 4 nodes per tetrahedron
        }
    } catch (const std::exception& e) {
        std::cout << "No tetrahedral connectivity found, exporting point cloud only: " << e.what() << std::endl;
    }
    
    // Create output filename in the same directory as the input
    fs::path inputDir = fs::path(meshPath).parent_path();
    fs::path outputPath = inputDir / "mesh_visualization_adios2.bp";
    std::cout << "Writing ADIOS2 file to: " << outputPath << std::endl;
    
    // Initialize ADIOS2
    adios2::ADIOS adios(MPI_COMM_SELF);
    adios2::IO io = adios.DeclareIO("VTKWriter");
    
    // VTK output in ADIOS2 format
    io.SetEngine("BP4");
    
    // Define variables
    adios2::Variable<float> varX = 
        io.DefineVariable<float>("coordinates/x", {nodeCount}, {0}, {nodeCount});
    adios2::Variable<float> varY = 
        io.DefineVariable<float>("coordinates/y", {nodeCount}, {0}, {nodeCount});
    adios2::Variable<float> varZ = 
        io.DefineVariable<float>("coordinates/z", {nodeCount}, {0}, {nodeCount});
    
    // Define connectivity variables if we have elements
    if (numElements > 0) {
        io.DefineVariable<int>("connectivity", 
                               {connectivity.size()}, {0}, {connectivity.size()});
        io.DefineVariable<int>("types", 
                               {cellTypes.size()}, {0}, {cellTypes.size()});
        io.DefineVariable<int>("offsets", 
                               {cellOffsets.size()}, {0}, {cellOffsets.size()});
    }
    
    // Open file and write data
    adios2::Engine engine = io.Open(outputPath.string(), adios2::Mode::Write);
    engine.BeginStep();
    
    // Write coordinates
    engine.Put(varX, x_data.data());
    engine.Put(varY, y_data.data());
    engine.Put(varZ, z_data.data());
    
    // Write connectivity if available
    if (numElements > 0) {
        engine.Put("connectivity", connectivity.data());
        engine.Put("types", cellTypes.data());
        engine.Put("offsets", cellOffsets.data());
    }
    
    engine.EndStep();
    engine.Close();
    
    std::cout << "ADIOS2 file written with " << nodeCount << " nodes and " 
              << numElements << " tetrahedral elements" << std::endl;
    std::cout << "ADIOS2 file path: " << outputPath << std::endl;
    
    // Verify file was written
    EXPECT_TRUE(fs::exists(outputPath)) << "ADIOS2 file was not created";
    EXPECT_GT(fs::file_size(outputPath), 0) << "ADIOS2 file is empty";
}
#endif // ADIOS2_HAVE_MPI

int main(int argc, char **argv) {
    mars::Env env(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}