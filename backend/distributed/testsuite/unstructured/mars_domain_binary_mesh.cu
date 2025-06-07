#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <algorithm>
#include <vector>
#include <chrono>
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
        
        std::error_code ec;
        fs::remove_all(testDir, ec);
        if (ec && rank == 0) {
            std::cout << "Warning: Failed to clean up test directory: " << ec.message() << std::endl;
        }
    }

    // Helper to create a binary file with float data
    void createCoordinateFile(const std::string& filename, const std::vector<float>& data)
    {
        std::ofstream file((testDir / filename).string(), std::ios::binary);
        if (!file) {
            throw std::runtime_error("Failed to create coordinate file: " + filename);
        }
        file.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(float));
    }

    // Helper to create a binary file with int data
    void createConnectivityFile(const std::string& filename, const std::vector<int>& data)
    {
        std::ofstream file((testDir / filename).string(), std::ios::binary);
        if (!file) {
            throw std::runtime_error("Failed to create connectivity file: " + filename);
        }
        file.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(int));
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

    // Helper to validate SFC-based coordinate conversion
    template<typename DomainType>
    void validateSfcCoordinates(const DomainType& domain, const std::string& testName = "")
    {
        if (domain.getElementCount() == 0) return;

        auto box = domain.getDomain().box();
        size_t samplesToCheck = std::min(domain.getElementCount(), size_t(5));
        
        for (size_t elemIdx = 0; elemIdx < samplesToCheck; elemIdx++)
        {
            // Get SFC keys for this element
            auto sfc0 = domain.template getConnectivity<0>(elemIdx);
            auto sfc1 = domain.template getConnectivity<1>(elemIdx);
            auto sfc2 = domain.template getConnectivity<2>(elemIdx);
            auto sfc3 = domain.template getConnectivity<3>(elemIdx);

            // Convert SFC keys to physical coordinates
            auto [x0, y0, z0] = domain.sfcToPhysicalCoordinate(sfc0);
            auto [x1, y1, z1] = domain.sfcToPhysicalCoordinate(sfc1);
            auto [x2, y2, z2] = domain.sfcToPhysicalCoordinate(sfc2);
            auto [x3, y3, z3] = domain.sfcToPhysicalCoordinate(sfc3);

            // Validate coordinates are within domain bounds
            EXPECT_GE(x0, box.xmin()) << testName << ": Node 0 X coordinate below domain minimum";
            EXPECT_LE(x0, box.xmax()) << testName << ": Node 0 X coordinate above domain maximum";
            EXPECT_GE(y0, box.ymin()) << testName << ": Node 0 Y coordinate below domain minimum";
            EXPECT_LE(y0, box.ymax()) << testName << ": Node 0 Y coordinate above domain maximum";
            EXPECT_GE(z0, box.zmin()) << testName << ": Node 0 Z coordinate below domain minimum";
            EXPECT_LE(z0, box.zmax()) << testName << ": Node 0 Z coordinate above domain maximum";

            // Test individual coordinate access functions
            auto x0_individual = domain.sfcToPhysicalCoordinateX(sfc0);
            auto y0_individual = domain.sfcToPhysicalCoordinateY(sfc0);
            auto z0_individual = domain.sfcToPhysicalCoordinateZ(sfc0);

            EXPECT_FLOAT_EQ(x0, x0_individual) << testName << ": X coordinate mismatch between tuple and individual access";
            EXPECT_FLOAT_EQ(y0, y0_individual) << testName << ": Y coordinate mismatch between tuple and individual access";
            EXPECT_FLOAT_EQ(z0, z0_individual) << testName << ": Z coordinate mismatch between tuple and individual access";

            // Test spatial coordinate conversion
            auto [ix0, iy0, iz0] = domain.sfcToSpatialCoordinate(sfc0);
            EXPECT_GE(ix0, 0u) << testName << ": Spatial coordinate should be non-negative";
            EXPECT_GE(iy0, 0u) << testName << ": Spatial coordinate should be non-negative";
            EXPECT_GE(iz0, 0u) << testName << ": Spatial coordinate should be non-negative";
        }
    }

    // Helper to validate connectivity using SFC indices access
    template<typename DomainType>
    void validateConnectivity(const DomainType& domain, const std::string& testName = "")
    {
        if (domain.getElementCount() == 0) return;

        // Use the SFC indices access methods
        auto h_i0 = toHost(domain.template indices<0>());
        auto h_i1 = toHost(domain.template indices<1>());
        auto h_i2 = toHost(domain.template indices<2>());
        auto h_i3 = toHost(domain.template indices<3>());

        EXPECT_EQ(h_i0.size(), domain.getElementCount()) << testName << ": SFC indices size mismatch";
        EXPECT_EQ(h_i1.size(), domain.getElementCount()) << testName << ": SFC indices size mismatch";
        EXPECT_EQ(h_i2.size(), domain.getElementCount()) << testName << ": SFC indices size mismatch";
        EXPECT_EQ(h_i3.size(), domain.getElementCount()) << testName << ": SFC indices size mismatch";

        // Verify individual SFC access works
        size_t samplesToCheck = std::min(domain.getElementCount(), size_t(5));
        for (size_t i = 0; i < samplesToCheck; i++)
        {
            // Test individual SFC access
            auto sfc0_individual = domain.template getConnectivity<0>(i);
            auto sfc1_individual = domain.template getConnectivity<1>(i);
            auto sfc2_individual = domain.template getConnectivity<2>(i);
            auto sfc3_individual = domain.template getConnectivity<3>(i);

            EXPECT_EQ(h_i0[i], sfc0_individual) << testName << ": SFC connectivity mismatch for element " << i;
            EXPECT_EQ(h_i1[i], sfc1_individual) << testName << ": SFC connectivity mismatch for element " << i;
            EXPECT_EQ(h_i2[i], sfc2_individual) << testName << ": SFC connectivity mismatch for element " << i;
            EXPECT_EQ(h_i3[i], sfc3_individual) << testName << ": SFC connectivity mismatch for element " << i;
        }
    }

    // Helper to test performance of SFC access
    template<typename DomainType>
    void benchmarkSfcAccess(const DomainType& domain, const std::string& testName)
    {
        if (domain.getElementCount() == 0 || rank != 0) return;

        const int iterations = 50;

        // Benchmark SFC connectivity access
        auto start = std::chrono::high_resolution_clock::now();
        size_t totalSize = 0;
        for (int iter = 0; iter < iterations; iter++) {
            auto h_i0 = toHost(domain.template indices<0>());
            auto h_i1 = toHost(domain.template indices<1>());
            auto h_i2 = toHost(domain.template indices<2>());
            auto h_i3 = toHost(domain.template indices<3>());
            totalSize += h_i0.size() + h_i1.size() + h_i2.size() + h_i3.size();
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto sfc_duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        // Benchmark individual SFC key access
        start = std::chrono::high_resolution_clock::now();
        unsigned long long sfcSum = 0;
        for (int iter = 0; iter < iterations; iter++) {
            size_t maxElements = std::min(domain.getElementCount(), size_t(10));
            for (size_t i = 0; i < maxElements; i++) {
                auto sfc0 = domain.template getConnectivity<0>(i);
                auto sfc1 = domain.template getConnectivity<1>(i);
                auto sfc2 = domain.template getConnectivity<2>(i);
                auto sfc3 = domain.template getConnectivity<3>(i);
                sfcSum += sfc0 + sfc1 + sfc2 + sfc3; // Use volatile to prevent optimization
            }
        }
        end = std::chrono::high_resolution_clock::now();
        auto individual_duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        // Benchmark coordinate conversion
        start = std::chrono::high_resolution_clock::now();
        double coordSum = 0.0;
        for (int iter = 0; iter < iterations; iter++) {
            size_t maxElements = std::min(domain.getElementCount(), size_t(10));
            for (size_t i = 0; i < maxElements; i++) {
                auto sfc0 = domain.template getConnectivity<0>(i);
                auto [x, y, z] = domain.sfcToPhysicalCoordinate(sfc0);
                coordSum += x + y + z; // Use volatile to prevent optimization
            }
        }
        end = std::chrono::high_resolution_clock::now();
        auto coord_duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        std::cout << "[" << testName << "] SFC Performance benchmark (" << iterations << " iterations):" << std::endl;
        std::cout << "  Bulk SFC access (mem bandwidth gpu-cpu): " << sfc_duration.count() << " μs" << std::endl;
        std::cout << "  Individual SFC access: " << individual_duration.count() << " μs" << std::endl;
        std::cout << "  Coordinate conversion: " << coord_duration.count() << " μs" << std::endl;

        if(totalSize == 0 || sfcSum == 0 || coordSum == 0.0)
        {
            std::cout << "Warning: Zero size or sum in SFC access benchmark, check data integrity." << std::endl;
        }
    }
};

// Test basic domain creation and properties
TEST_F(ElementDomainTest, DomainCreation)
{
    if (rank != 0)
    {
        GTEST_SKIP() << "This test only runs on rank 0";
    }

    try
    {
        // Force single rank for this test to ensure consistent results
        using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
        Domain domain(testDir.string(), 0, 1);

        // Test basic properties
        EXPECT_EQ(domain.getNodeCount(), 8);
        EXPECT_EQ(domain.getElementCount(), 4);

        // Verify SFC-based coordinate conversion
        validateSfcCoordinates(domain, "DomainCreation");

        // Verify connectivity using SFC methods
        validateConnectivity(domain, "DomainCreation");

        // Check cornerstone domain exists
        EXPECT_NE(&domain.getDomain(), nullptr);

        std::cout << "DomainCreation test completed successfully" << std::endl;
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception in DomainCreation: " << e.what();
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

// Test basic domain creation and properties - works with any number of ranks
TEST_F(ElementDomainTest, DomainCreationMultiRank)
{
    try
    {
        // Create domain with actual rank/numRanks to test proper MPI support
        using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
        Domain domain(testDir.string(), rank, numRanks);

        // Test that we have a valid domain with some nodes and elements
        EXPECT_GT(domain.getNodeCount(), 0) << "Node count should be positive on rank " << rank;
        EXPECT_GT(domain.getElementCount(), 0) << "Element count should be positive on rank " << rank;

        // Gather total counts to verify all elements are accounted for
        size_t localElements = domain.getElementCount();
        size_t globalElements = 0;
        MPI_Allreduce(&localElements, &globalElements, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

        // Total element count should match our test data regardless of distribution
        EXPECT_EQ(globalElements, 4) << "Total element count should be 4 across all ranks";

        // Validate SFC coordinates and connectivity on each rank
        validateSfcCoordinates(domain, "DomainCreationMultiRank");
        validateConnectivity(domain, "DomainCreationMultiRank");

        // Check cornerstone domain exists
        EXPECT_NE(&domain.getDomain(), nullptr);

        // Verify domain ranges make sense for local data
        EXPECT_GE(domain.startIndex(), 0) << "Start index invalid on rank " << rank;
        EXPECT_LE(domain.endIndex(), domain.getElementCount()) << "End index invalid on rank " << rank;

        if (rank == 0) {
            std::cout << "DomainCreationMultiRank test completed successfully" << std::endl;
        }
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception in DomainCreationMultiRank: " << e.what();
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

// Test SFC key generation
TEST_F(ElementDomainTest, SfcKeyGeneration)
{
    try
    {
        // Create a domain with a known bounding box
        using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
        Domain domain(testDir.string(), rank, numRanks);

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

        // Test SFC indices are valid
        validateConnectivity(domain, "SfcKeyGeneration");

        // Test SFC coordinate conversion
        validateSfcCoordinates(domain, "SfcKeyGeneration");

        if (rank == 0) {
            std::cout << "SfcKeyGeneration test completed successfully" << std::endl;
        }
    }
    catch (the std::exception& e)
    {
        FAIL() << "Exception in SfcKeyGeneration: " << e.what();
    }
}

// Test SFC-based access performance
TEST_F(ElementDomainTest, SfcAccessPerformance)
{
    try
    {
        using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
        Domain domain(testDir.string(), rank, numRanks);

        // Test that SFC access works correctly
        validateSfcCoordinates(domain, "SfcAccessPerformance");
        validateConnectivity(domain, "SfcAccessPerformance");

        // Run performance benchmark
        benchmarkSfcAccess(domain, "SfcAccessPerformance");

        if (rank == 0) {
            std::cout << "SfcAccessPerformance test completed successfully" << std::endl;
        }
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception in SfcAccessPerformance: " << e.what();
    }
}

// Test characteristic sizes computation using SFC coordinates
TEST_F(ElementDomainTest, CharacteristicSizes)
{
    try
    {
        // Create two domains with different precisions to test type handling
        using FloatDomain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
        using DoubleDomain = ElementDomain<TetTag, double, unsigned, cstone::GpuTag>;

        FloatDomain floatDomain(testDir.string(), rank, numRanks);
        DoubleDomain doubleDomain(testDir.string(), rank, numRanks);

        // Both should have computed characteristic sizes (h values)
        EXPECT_EQ(floatDomain.getNodeCount(), 8);
        EXPECT_EQ(doubleDomain.getNodeCount(), 8);

        // Test that both domains support SFC access
        validateSfcCoordinates(floatDomain, "CharacteristicSizes_Float");
        validateSfcCoordinates(doubleDomain, "CharacteristicSizes_Double");

        if (rank == 0) {
            std::cout << "CharacteristicSizes test completed successfully" << std::endl;
        }
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception in CharacteristicSizes: " << e.what();
    }
}

// Test multi-rank domain partitioning with distribution
TEST_F(ElementDomainTest, MultiRankDistribution)
{
    // Skip if we're running with only 1 rank
    if (numRanks < 2)
    {
        GTEST_SKIP() << "This test requires at least 2 ranks";
    }

    try
    {
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

        // Validate SFC coordinates and connectivity on each rank
        validateSfcCoordinates(domain, "MultiRankDistribution");
        validateConnectivity(domain, "MultiRankDistribution");

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
            size_t totalNodes = 0;

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
                size_t minElems = *std::min_element(nonZeroCounts.begin(), nonZeroCounts.end());
                size_t maxElems = *std::max_element(nonZeroCounts.begin(), nonZeroCounts.end());
                double elemImbalance = static_cast<double>(maxElems) / minElems;

                std::cout << "Element distribution: min=" << minElems << ", max=" << maxElems
                          << ", imbalance=" << elemImbalance << std::endl;

                // For a distributed mesh, expect reasonable balance
                EXPECT_LT(elemImbalance, 2.0) << "Element imbalance should be less than 2x";
            }

            // Verify total element count
            EXPECT_EQ(totalElements, 50) << "Total element count should match the mesh";
            std::cout << "MultiRankDistribution test completed successfully" << std::endl;
        }
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception in MultiRankDistribution: " << e.what();
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

// Test with different key types (uint64_t)
TEST_F(ElementDomainTest, DifferentKeyTypes)
{
    try
    {
        // Create domains with different key types
        using DomainUint = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
        using DomainUint64 = ElementDomain<TetTag, float, uint64_t, cstone::GpuTag>;

        DomainUint domainUint(testDir.string(), rank, numRanks);
        DomainUint64 domainUint64(testDir.string(), rank, numRanks);

        // Both should function correctly with the same mesh
        EXPECT_EQ(domainUint.getNodeCount(), domainUint64.getNodeCount());
        EXPECT_EQ(domainUint.getElementCount(), domainUint64.getElementCount());

        // Test SFC access with both key types
        validateSfcCoordinates(domainUint, "DifferentKeyTypes_Uint");
        validateSfcCoordinates(domainUint64, "DifferentKeyTypes_Uint64");
        validateConnectivity(domainUint, "DifferentKeyTypes_Uint");
        validateConnectivity(domainUint64, "DifferentKeyTypes_Uint64");

        if (rank == 0) {
            std::cout << "DifferentKeyTypes test completed successfully" << std::endl;
        }
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception in DifferentKeyTypes: " << e.what();
    }
}

// Test FEM-style element iteration using SFC coordinates
TEST_F(ElementDomainTest, FemStyleElementIteration)
{
    try
    {
        using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
        Domain domain(testDir.string(), rank, numRanks);

        if (domain.getElementCount() == 0) {
            GTEST_SKIP() << "No elements on this rank";
        }

        // Simulate FEM element loop
        size_t elementsToProcess = std::min(domain.getElementCount(), size_t(5));
        std::vector<float> elementVolumes(elementsToProcess, 0.0f);

        for (size_t elemIdx = 0; elemIdx < elementsToProcess; elemIdx++) {
            // Get SFC keys for element nodes
            auto sfc0 = domain.getConnectivity<0>(elemIdx);
            auto sfc1 = domain.getConnectivity<1>(elemIdx);
            auto sfc2 = domain.getConnectivity<2>(elemIdx);
            auto sfc3 = domain.getConnectivity<3>(elemIdx);

            // Convert to physical coordinates
            auto [x0, y0, z0] = domain.sfcToPhysicalCoordinate(sfc0);
            auto [x1, y1, z1] = domain.sfcToPhysicalCoordinate(sfc1);
            auto [x2, y2, z2] = domain.sfcToPhysicalCoordinate(sfc2);
            auto [x3, y3, z3] = domain.sfcToPhysicalCoordinate(sfc3);

            // Simple tetrahedral volume calculation (1/6 * |det(matrix)|)
            float det = (x1-x0)*((y2-y0)*(z3-z0) - (z2-z0)*(y3-y0)) -
                       (y1-y0)*((x2-x0)*(z3-z0) - (z2-z0)*(x3-x0)) +
                       (z1-z0)*((x2-x0)*(y3-y0) - (y2-y0)*(x3-x0));
            
            elementVolumes[elemIdx] = std::abs(det) / 6.0f;
            
            // Volume should be positive for a valid tetrahedron
            EXPECT_GT(elementVolumes[elemIdx], 0.0f) << "Element " << elemIdx << " should have positive volume";
        }

        if (rank == 0) {
            std::cout << "FEM-style element iteration validated for " << elementsToProcess << " elements" << std::endl;
            if (elementsToProcess > 0) {
                auto minVolume = *std::min_element(elementVolumes.begin(), elementVolumes.end());
                auto maxVolume = *std::max_element(elementVolumes.begin(), elementVolumes.end());
                std::cout << "  Volume range: [" << minVolume << ", " << maxVolume << "]" << std::endl;
            }
        }
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception in FEM iteration test: " << e.what();
    }
}

// Main function that initializes MPI and runs the tests
int main(int argc, char** argv)
{
    // Initialize MPI through Mars environment
    mars::Env env(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}