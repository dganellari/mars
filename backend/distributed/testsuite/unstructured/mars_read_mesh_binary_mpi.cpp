#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <unistd.h>
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
        // One directory per process: ctest may run these tests in parallel, and each rank reads
        // its own copy of the files.
        testDir = fs::temp_directory_path() / ("mars_mesh_binary_mpi_test_" + std::to_string(getpid()));
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

// Each MPI rank reads its own element slice. Element e uses nodes 2e..2e+3, so neighbouring
// slices share two nodes. All collective calls come before the checks: a rank that stopped early
// would leave the others waiting.
TEST_F(MeshReadBinaryMPITest, MultiRankDistribution) {
    int rank, numRanks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    constexpr int numElements = 50;
    constexpr int numNodes    = 2 * numElements + 2;

    std::vector<float> x_coords(numNodes), y_coords(numNodes), z_coords(numNodes);
    for (int i = 0; i < numNodes; i++) {
        x_coords[i] = static_cast<float>(i);
        y_coords[i] = static_cast<float>(i) * 0.5f;
        z_coords[i] = static_cast<float>(i) * 10.0f;
    }
    createCoordinateFile("x.float32", x_coords);
    createCoordinateFile("y.float32", y_coords);
    createCoordinateFile("z.float32", z_coords);

    std::vector<int> i0(numElements), i1(numElements), i2(numElements), i3(numElements);
    for (int e = 0; e < numElements; e++) {
        i0[e] = 2 * e;
        i1[e] = 2 * e + 1;
        i2[e] = 2 * e + 2;
        i3[e] = 2 * e + 3;
    }
    createConnectivityFile("i0.int32", i0);
    createConnectivityFile("i1.int32", i1);
    createConnectivityFile("i2.int32", i2);
    createConnectivityFile("i3.int32", i3);

    auto [nodeCount, elementCount, x, y, z, conn, localToGlobal] =
        mars::readMeshWithElementPartitioning<4, float>(testDir.string(), rank, numRanks);

    unsigned long localElements = elementCount;
    unsigned long totalElements = 0;
    MPI_Allreduce(&localElements, &totalElements, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

    std::vector<int> nodeSeen(numNodes, 0);
    for (auto g : localToGlobal) {
        if (g < static_cast<unsigned>(numNodes)) { nodeSeen[g] = 1; }
    }
    MPI_Allreduce(MPI_IN_PLACE, nodeSeen.data(), numNodes, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    EXPECT_EQ(totalElements, static_cast<unsigned long>(numElements));
    for (int g = 0; g < numNodes; g++) {
        EXPECT_EQ(nodeSeen[g], 1) << "node " << g << " is on no rank";
    }

    // Contiguous slices; the last rank takes the remainder.
    size_t perRank   = numElements / numRanks;
    size_t firstElem = rank * perRank;
    size_t expected  = (rank == numRanks - 1) ? numElements - firstElem : perRank;
    EXPECT_EQ(elementCount, expected) << "rank " << rank << " of " << numRanks;
    EXPECT_EQ(nodeCount, elementCount > 0 ? 2 * elementCount + 2 : 0) << "rank " << rank;
    ASSERT_EQ(localToGlobal.size(), nodeCount) << "rank " << rank;

    const std::vector<unsigned>* corners[4] = {&std::get<0>(conn), &std::get<1>(conn), &std::get<2>(conn),
                                               &std::get<3>(conn)};
    for (size_t e = 0; e < elementCount; e++) {
        for (int k = 0; k < 4; k++) {
            unsigned local = (*corners[k])[e];
            ASSERT_LT(local, nodeCount) << "rank " << rank << ", element " << e << ", corner " << k;
            EXPECT_EQ(localToGlobal[local], 2 * (firstElem + e) + k)
                << "rank " << rank << ", element " << e << ", corner " << k;
        }
    }

    for (size_t l = 0; l < nodeCount; l++) {
        float g = static_cast<float>(localToGlobal[l]);
        EXPECT_FLOAT_EQ(x[l], g) << "rank " << rank << ", node " << localToGlobal[l];
        EXPECT_FLOAT_EQ(y[l], g * 0.5f) << "rank " << rank << ", node " << localToGlobal[l];
        EXPECT_FLOAT_EQ(z[l], g * 10.0f) << "rank " << rank << ", node " << localToGlobal[l];
    }
}

// One element on several ranks: all ranks but the last get no elements and must get no nodes.
TEST_F(MeshReadBinaryMPITest, RanksWithoutElements) {
    int rank, numRanks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    createConnectivityFile("i0.int32", {0});
    createConnectivityFile("i1.int32", {1});
    createConnectivityFile("i2.int32", {2});
    createConnectivityFile("i3.int32", {3});

    auto [nodeCount, elementCount, x, y, z, conn, localToGlobal] =
        mars::readMeshWithElementPartitioning<4, float>(testDir.string(), rank, numRanks);

    bool last = (rank == numRanks - 1);
    EXPECT_EQ(elementCount, last ? 1u : 0u) << "rank " << rank;
    EXPECT_EQ(nodeCount, last ? 4u : 0u) << "rank " << rank;
    EXPECT_EQ(x.size(), nodeCount) << "rank " << rank;
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