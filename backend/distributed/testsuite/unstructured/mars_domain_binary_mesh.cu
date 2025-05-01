#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <algorithm>
#include <vector>
#include "domain.hpp"

namespace fs = std::filesystem;
using namespace mars;

class ElementDomainTest : public ::testing::Test
{
protected:
    // Test directory for mesh files - rank-specific to avoid conflicts
    fs::path testDir;
    int rank;
    int numRanks;

    void SetUp() override
    {
        // Get actual MPI rank and size
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

        // Create a temporary directory with rank-specific name
        testDir = fs::temp_directory_path() / ("mars_domain_test_rank_" + std::to_string(rank));
        fs::create_directories(testDir);

        // Create basic coordinate files
        createCoordinateFile("x.float32", {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f});
        createCoordinateFile("y.float32", {0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f});
        createCoordinateFile("z.float32", {10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f});

        // Create connectivity files for tetrahedra
        createConnectivityFile("i0.int32", {0, 2, 4, 6});
        createConnectivityFile("i1.int32", {1, 3, 5, 7});
        createConnectivityFile("i2.int32", {2, 4, 6, 0});
        createConnectivityFile("i3.int32", {3, 5, 7, 1});

        // Ensure all ranks have created their files before proceeding
        MPI_Barrier(MPI_COMM_WORLD);
    }

    void TearDown() override
    {
        // Synchronize before cleanup
        MPI_Barrier(MPI_COMM_WORLD);
        fs::remove_all(testDir);
    }

    // Helper to create a binary file with float data
    void createCoordinateFile(const std::string& filename, const std::vector<float>& data)
    {
        std::ofstream file((testDir / filename).string(), std::ios::binary);
        file.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(float));
        file.close();
    }

    // Helper to create a binary file with int data
    void createConnectivityFile(const std::string& filename, const std::vector<int>& data)
    {
        std::ofstream file((testDir / filename).string(), std::ios::binary);
        file.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(int));
        file.close();
    }

    // Create a more complex mesh with node sharing for testing
    void createComplexMesh(int numNodes, int numElements)
    {
        // Create coordinate arrays
        std::vector<float> x_coords(numNodes);
        std::vector<float> y_coords(numNodes);
        std::vector<float> z_coords(numNodes);

        // Generate coordinates in a grid pattern
        for (int i = 0; i < numNodes; i++)
        {
            int x_idx = i % 10;
            int y_idx = (i / 10) % 10;
            int z_idx = i / 100;

            x_coords[i] = static_cast<float>(x_idx);
            y_coords[i] = static_cast<float>(y_idx);
            z_coords[i] = static_cast<float>(z_idx);
        }

        // Write coordinate files
        createCoordinateFile("x.float32", x_coords);
        createCoordinateFile("y.float32", y_coords);
        createCoordinateFile("z.float32", z_coords);

        // Create connectivity for tetrahedra - ensure node sharing between elements
        std::vector<int> i0(numElements);
        std::vector<int> i1(numElements);
        std::vector<int> i2(numElements);
        std::vector<int> i3(numElements);

        for (int i = 0; i < numElements; i++)
        {
            // Create tetrahedra that share vertices
            i0[i] = i % numNodes;                  // First vertex
            i1[i] = (i + 1) % numNodes;            // Second vertex
            i2[i] = (i + 10) % numNodes;           // Third vertex
            i3[i] = (i + numNodes - 1) % numNodes; // Fourth vertex
        }

        // Write connectivity files
        createConnectivityFile("i0.int32", i0);
        createConnectivityFile("i1.int32", i1);
        createConnectivityFile("i2.int32", i2);
        createConnectivityFile("i3.int32", i3);
    }
};

// Test basic domain creation and properties
TEST_F(ElementDomainTest, DomainCreation)
{
    if (rank != 0)
    {
        GTEST_SKIP() << "This test only runs on rank 0";
        return;
    }
    // Force single rank for this test to ensure consistent results
    using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
    Domain domain(testDir.string(), 0, 1);

    // Test basic properties
    EXPECT_EQ(domain.getNodeCount(), 8);
    EXPECT_EQ(domain.getElementCount(), 4);

    // Verify that coordinates were loaded correctly by checking values on host
    auto h_x = toHost(domain.x());
    auto h_y = toHost(domain.y());
    auto h_z = toHost(domain.z());

    EXPECT_FLOAT_EQ(h_x[0], 1.0f);
    EXPECT_FLOAT_EQ(h_y[2], 0.3f);
    EXPECT_FLOAT_EQ(h_z[5], 60.0f);

    // Verify connectivity
    auto h_i0 = toHost(domain.indices<0>());
    auto h_i1 = toHost(domain.indices<1>());
    auto h_i2 = toHost(domain.indices<2>());
    auto h_i3 = toHost(domain.indices<3>());

    // Check indices from first element
    EXPECT_EQ(h_i0[0], 0);
    EXPECT_EQ(h_i1[0], 1);
    EXPECT_EQ(h_i2[0], 2);
    EXPECT_EQ(h_i3[0], 3);

    // Check cornerstone domain exists
    EXPECT_NE(&domain.getDomain(), nullptr);

    MPI_Barrier(MPI_COMM_WORLD);
}

// Test basic domain creation and properties - works with any number of ranks
TEST_F(ElementDomainTest, DomainCreationMultiRank)
{
    // Create domain with actual rank/numRanks to test proper MPI support
    using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
    Domain domain(testDir.string(), rank, numRanks);

    // Test that we have a valid domain with some nodes and elements
    EXPECT_GT(domain.getNodeCount(), 0) << "Node count should be positive on rank " << rank;
    EXPECT_GT(domain.getElementCount(), 0) << "Element count should be positive on rank " << rank;

    // Gather total counts to verify all elements are accounted for
    size_t localElements  = domain.getElementCount();
    size_t globalElements = 0;
    MPI_Allreduce(&localElements, &globalElements, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

    // Total element count should match our test data regardless of distribution
    EXPECT_EQ(globalElements, 4) << "Total element count should be 4 across all ranks";

    // Verify that coordinates were loaded correctly by checking values are in expected ranges
    auto h_x = toHost(domain.x());
    auto h_y = toHost(domain.y());
    auto h_z = toHost(domain.z());

    // Instead of checking specific indices (which depend on distribution),
    // verify coordinate ranges and that arrays are properly sized
    EXPECT_EQ(h_x.size(), domain.getNodeCount());
    EXPECT_EQ(h_y.size(), domain.getNodeCount());
    EXPECT_EQ(h_z.size(), domain.getNodeCount());

    // Check all coordinates are in valid ranges
    for (size_t i = 0; i < domain.getNodeCount(); i++)
    {
        EXPECT_GE(h_x[i], 1.0f) << "X coordinate out of range on rank " << rank;
        EXPECT_LE(h_x[i], 8.0f) << "X coordinate out of range on rank " << rank;
        EXPECT_GE(h_y[i], 0.1f) << "Y coordinate out of range on rank " << rank;
        EXPECT_LE(h_y[i], 0.8f) << "Y coordinate out of range on rank " << rank;
        EXPECT_GE(h_z[i], 10.0f) << "Z coordinate out of range on rank " << rank;
        EXPECT_LE(h_z[i], 80.0f) << "Z coordinate out of range on rank " << rank;
    }

    // Verify connectivity - check that indices are valid
    auto h_i0 = toHost(domain.indices<0>());
    auto h_i1 = toHost(domain.indices<1>());
    auto h_i2 = toHost(domain.indices<2>());
    auto h_i3 = toHost(domain.indices<3>());

    // Check the connectivity arrays are properly sized
    EXPECT_EQ(h_i0.size(), domain.getElementCount());
    EXPECT_EQ(h_i1.size(), domain.getElementCount());
    EXPECT_EQ(h_i2.size(), domain.getElementCount());
    EXPECT_EQ(h_i3.size(), domain.getElementCount());

    // All node indices should be valid for this rank's portion
    for (size_t i = 0; i < domain.getElementCount(); i++)
    {
        EXPECT_LT(h_i0[i], domain.getNodeCount()) << "Node index out of bounds on rank " << rank;
        EXPECT_LT(h_i1[i], domain.getNodeCount()) << "Node index out of bounds on rank " << rank;
        EXPECT_LT(h_i2[i], domain.getNodeCount()) << "Node index out of bounds on rank " << rank;
        EXPECT_LT(h_i3[i], domain.getNodeCount()) << "Node index out of bounds on rank " << rank;
    }

    // Check cornerstone domain exists
    EXPECT_NE(&domain.getDomain(), nullptr);

    // Verify domain ranges make sense for local data
    EXPECT_GE(domain.startIndex(), 0) << "Start index invalid on rank " << rank;
    EXPECT_LE(domain.endIndex(), domain.getElementCount()) << "End index invalid on rank " << rank;

    // Make sure all tests complete before moving to next test
    MPI_Barrier(MPI_COMM_WORLD);
}

// Test SFC key generation
TEST_F(ElementDomainTest, SfcKeyGeneration)
{
    // Create a domain with a known bounding box
    using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
    Domain domain(testDir.string(), rank, numRanks);

    // Collect global statistics on the SFC keys
    // We need to extract element SFC codes - we'll use a custom GPU kernel to access them

    // First, get the domain dimensions
    auto box = domain.getDomain().box();

    // Verify box dimensions match our test data
    EXPECT_GE(box.xmin(), 0.0f);
    EXPECT_LE(box.xmax(), 10.0f);
    EXPECT_GE(box.ymin(), 0.0f);
    EXPECT_LE(box.ymax(), 1.0f);
    EXPECT_GE(box.zmin(), 0.0f);
    EXPECT_LE(box.zmax(), 100.0f);

    // For actual elements, just check that start/end indices are sane
    EXPECT_GE(domain.startIndex(), 0);
    EXPECT_LE(domain.endIndex(), domain.getElementCount());
}

// Test characteristic sizes computation
TEST_F(ElementDomainTest, CharacteristicSizes)
{
    // Create two domains with different precisions to test type handling
    using FloatDomain  = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
    using DoubleDomain = ElementDomain<TetTag, double, unsigned, cstone::GpuTag>;

    FloatDomain floatDomain(testDir.string(), rank, numRanks);
    DoubleDomain doubleDomain(testDir.string(), rank, numRanks);

    // Both should have computed characteristic sizes (h values)
    // We can't see them directly, but we can check the domain was created successfully
    EXPECT_EQ(floatDomain.getNodeCount(), 8);
    EXPECT_EQ(doubleDomain.getNodeCount(), 8);
}

// Test multi-rank domain partitioning with distribution
TEST_F(ElementDomainTest, MultiRankDistribution)
{
    // Skip if we're running with only 1 rank
    if (numRanks < 2)
    {
        GTEST_SKIP() << "This test requires at least 2 ranks";
        return;
    }

    // Create a larger, more complex mesh for meaningful partitioning
    createComplexMesh(100, 50);

    // Wait for all ranks to finish creating their meshes
    MPI_Barrier(MPI_COMM_WORLD);

    // Create domain with the current rank and number of ranks
    using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
    Domain domain(testDir.string(), rank, numRanks);

    // Each rank verifies its own portion of the domain
    EXPECT_GE(domain.getNodeCount(), 0);
    EXPECT_GE(domain.getElementCount(), 0);

    // Gather element and node counts to analyze distribution
    std::vector<size_t> allElementCounts(numRanks);
    std::vector<size_t> allNodeCounts(numRanks);

    size_t elemCount = domain.getElementCount();
    size_t nodeCount = domain.getNodeCount();

    MPI_Gather(&elemCount, 1, MPI_UNSIGNED_LONG, allElementCounts.data(), 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    MPI_Gather(&nodeCount, 1, MPI_UNSIGNED_LONG, allNodeCounts.data(), 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    // Rank 0 analyzes distribution
    if (rank == 0)
    {
        size_t totalElements = 0;
        size_t totalNodes    = 0;

        for (int r = 0; r < numRanks; r++)
        {
            std::cout << "Rank " << r << " has " << allElementCounts[r] << " elements and " << allNodeCounts[r]
                      << " nodes" << std::endl;
            totalElements += allElementCounts[r];
            totalNodes += allNodeCounts[r];
        }

        // Check distribution balance if multiple ranks have elements
        std::vector<size_t> nonZeroCounts;
        for (auto count : allElementCounts)
        {
            if (count > 0) nonZeroCounts.push_back(count);
        }

        if (nonZeroCounts.size() > 1)
        {
            size_t minElems      = *std::min_element(nonZeroCounts.begin(), nonZeroCounts.end());
            size_t maxElems      = *std::max_element(nonZeroCounts.begin(), nonZeroCounts.end());
            double elemImbalance = static_cast<double>(maxElems) / minElems;

            std::cout << "Element distribution: min=" << minElems << ", max=" << maxElems
                      << ", imbalance=" << elemImbalance << std::endl;

            // For a distributed mesh, expect reasonable balance
            EXPECT_LT(elemImbalance, 2.0) << "Element imbalance should be less than 2x";
        }

        // Verify total element count
        EXPECT_EQ(totalElements, 50) << "Total element count should match the mesh";
    }

    // Check validity of connectivity indices
    if (domain.getElementCount() > 0)
    {
        auto h_i0 = toHost(domain.indices<0>());
        auto h_i1 = toHost(domain.indices<1>());
        auto h_i2 = toHost(domain.indices<2>());
        auto h_i3 = toHost(domain.indices<3>());

        for (size_t i = 0; i < domain.getElementCount(); i++)
        {
            EXPECT_LT(h_i0[i], domain.getNodeCount()) << "Node index out of bounds on rank " << rank;
            EXPECT_LT(h_i1[i], domain.getNodeCount()) << "Node index out of bounds on rank " << rank;
            EXPECT_LT(h_i2[i], domain.getNodeCount()) << "Node index out of bounds on rank " << rank;
            EXPECT_LT(h_i3[i], domain.getNodeCount()) << "Node index out of bounds on rank " << rank;
        }
    }

    // Make sure all ranks synchronize before ending
    MPI_Barrier(MPI_COMM_WORLD);
}

// Test representative node mapping - this tests the core of the domain logic
// where each element is assigned to a representative node for tree building
TEST_F(ElementDomainTest, RepresentativeNodeMapping)
{
    // Create a small mesh where we know which nodes should be representatives
    std::vector<float> x = {0.0f, 1.0f, 0.0f, 1.0f};
    std::vector<float> y = {0.0f, 0.0f, 1.0f, 1.0f};
    std::vector<float> z = {0.0f, 0.0f, 0.0f, 1.0f};

    std::vector<int> i0 = {0}; // Single tet with nodes 0,1,2,3
    std::vector<int> i1 = {1};
    std::vector<int> i2 = {2};
    std::vector<int> i3 = {3};

    createCoordinateFile("x.float32", x);
    createCoordinateFile("y.float32", y);
    createCoordinateFile("z.float32", z);
    createConnectivityFile("i0.int32", i0);
    createConnectivityFile("i1.int32", i1);
    createConnectivityFile("i2.int32", i2);
    createConnectivityFile("i3.int32", i3);

    // Create domain
    using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
    Domain domain(testDir.string(), rank, numRanks);

    // Force sync to trigger representative node mapping
    domain.sync();

    // Verify the domain was created correctly
    EXPECT_EQ(domain.getNodeCount(), 4);
    EXPECT_EQ(domain.getElementCount(), 1);

    // We can't directly access elemToNodeMap_ since it's private,
    // but we can verify the domain functions correctly by checking
    // that we have sensible start/end indices
    EXPECT_EQ(domain.startIndex(), 0);
    EXPECT_EQ(domain.endIndex(), 1);
}

// Test domain synchronization with a changing mesh
TEST_F(ElementDomainTest, DomainSynchronization)
{
    // First, create a domain with a small mesh
    createCoordinateFile("x.float32", {0.0f, 1.0f, 0.0f, 1.0f});
    createCoordinateFile("y.float32", {0.0f, 0.0f, 1.0f, 1.0f});
    createCoordinateFile("z.float32", {0.0f, 0.0f, 0.0f, 1.0f});
    createConnectivityFile("i0.int32", {0});
    createConnectivityFile("i1.int32", {1});
    createConnectivityFile("i2.int32", {2});
    createConnectivityFile("i3.int32", {3});

    // Create domain
    using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
    Domain domain(testDir.string(), rank, numRanks);

    // Now modify the mesh to simulate a changing topology
    createCoordinateFile("x.float32", {0.0f, 2.0f, 0.0f, 2.0f, 1.0f, 3.0f});
    createCoordinateFile("y.float32", {0.0f, 0.0f, 2.0f, 2.0f, 1.0f, 1.0f});
    createCoordinateFile("z.float32", {0.0f, 0.0f, 0.0f, 2.0f, 1.0f, 1.0f});
    createConnectivityFile("i0.int32", {0, 1});
    createConnectivityFile("i1.int32", {1, 4});
    createConnectivityFile("i2.int32", {2, 5});
    createConnectivityFile("i3.int32", {3, 0});

    // The domain will handle this change when reading the new data
    Domain domain2(testDir.string(), rank, numRanks);

    // Verify the domain was updated correctly
    EXPECT_EQ(domain2.getNodeCount(), 6);
    EXPECT_EQ(domain2.getElementCount(), 2);

    // Check that the cornerstone domain is functional
    EXPECT_GE(domain2.startIndex(), 0);
    EXPECT_LE(domain2.endIndex(), domain2.getElementCount());
}

// Test with different key types (uint64_t)
TEST_F(ElementDomainTest, DifferentKeyTypes)
{
    // Create domains with different key types
    using DomainUint   = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
    using DomainUint64 = ElementDomain<TetTag, float, uint64_t, cstone::GpuTag>;

    DomainUint domainUint(testDir.string(), rank, numRanks);
    DomainUint64 domainUint64(testDir.string(), rank, numRanks);

    // Both should function correctly with the same mesh
    EXPECT_EQ(domainUint.getNodeCount(), domainUint64.getNodeCount());
    EXPECT_EQ(domainUint.getElementCount(), domainUint64.getElementCount());
}

// Main function that initializes MPI and runs the tests
int main(int argc, char** argv)
{
    // Initialize MPI through Mars environment
    mars::Env env(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
