#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <algorithm>
#include <vector>
#include <chrono>
#include <mpi.h>
#include "domain.hpp"

namespace fs = std::filesystem;
using namespace mars;

class ExternalMeshDomainTest : public ::testing::Test
{
protected:
    int rank;
    int numRanks;
    std::string meshPath;
    int deviceCount = 0;

    void SetUp() override
    {
        // Get MPI rank and size
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
        
        // Find mesh path once for all tests
        meshPath = getMeshPath();
        
        // Check device availability once
        cudaGetDeviceCount(&deviceCount);
        
        // Print diagnostic message only once
        if (rank == 0) {
            std::cout << "Test setup: Found " << deviceCount << " CUDA devices" << std::endl;
            std::cout << "Using mesh at: " << (meshPath.empty() ? "none" : meshPath) << std::endl;
        }
    }

private:
    // Get mesh paths from environment variable
    std::string getMeshPath() const
    {
        const char* meshPathEnv = std::getenv("MESH_PATH");
        std::string path = meshPathEnv ? meshPathEnv : "";

        if (path.empty() || !fs::exists(path))
        {
            // Fallback to common test locations
            std::vector<std::string> commonLocations = {
                "./test_data", "../test_data", "./meshes", "../meshes",
                "../../test_data", "../../meshes"
            };

            for (const auto& loc : commonLocations)
            {
                if (fs::exists(loc) && fs::is_directory(fs::path(loc)))
                {
                    if (hasRequiredMeshFiles(loc))
                    {
                        return loc;
                    }
                }
            }
        }
        return path;
    }
    
    bool hasRequiredMeshFiles(const std::string& path) const
    {
        fs::path basePath(path);
        return fs::exists(basePath / "x.float32") || fs::exists(basePath / "x.double");
    }
    
protected:
    // Centralized check for test prerequisites
    void checkPrerequisites(const std::string& testName)
    {
        if (meshPath.empty() || !fs::exists(meshPath))
        {
            GTEST_SKIP() << testName << ": No valid mesh directory found";
        }

        if (!hasRequiredMeshFiles(meshPath))
        {
            GTEST_SKIP() << testName << ": Mesh directory does not contain required coordinate files";
        }

        if (deviceCount == 0)
        {
            GTEST_SKIP() << testName << ": No CUDA devices available or device initialization failed";
        }
    }

    // Helper to validate domain basic properties
    template<typename DomainType>
    void validateBasicDomainProperties(const DomainType& domain)
    {
        EXPECT_GT(domain.getNodeCount(), 0) << "Domain should have nodes";
        EXPECT_GT(domain.getElementCount(), 0) << "Domain should have elements";
        EXPECT_GE(domain.startIndex(), 0) << "Invalid start index";
        EXPECT_LE(domain.endIndex(), domain.getElementCount()) << "Invalid end index";
        EXPECT_NE(&domain.getDomain(), nullptr) << "Cornerstone domain should exist";
    }

    // Helper to validate bounding box
    template<typename DomainType>
    void validateBoundingBox(const DomainType& domain)
    {
        auto box = domain.getDomain().box();
        EXPECT_LT(box.xmin(), box.xmax()) << "Invalid X dimension in bounding box";
        EXPECT_LT(box.ymin(), box.ymax()) << "Invalid Y dimension in bounding box";
        EXPECT_LT(box.zmin(), box.zmax()) << "Invalid Z dimension in bounding box";
    }

    // Helper to validate SFC-to-coordinate conversion
    template<typename DomainType>
    void validateSfcCoordinateConversion(const DomainType& domain, size_t maxSamples = 10)
    {
        if (domain.getElementCount() == 0) return;

        auto box = domain.getDomain().box();
        size_t samplesToCheck = std::min(domain.getElementCount(), maxSamples);
        
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
            EXPECT_GE(x0, box.xmin()) << "Node 0 X coordinate below domain minimum";
            EXPECT_LE(x0, box.xmax()) << "Node 0 X coordinate above domain maximum";
            EXPECT_GE(y0, box.ymin()) << "Node 0 Y coordinate below domain minimum";
            EXPECT_LE(y0, box.ymax()) << "Node 0 Y coordinate above domain maximum";
            EXPECT_GE(z0, box.zmin()) << "Node 0 Z coordinate below domain minimum";
            EXPECT_LE(z0, box.zmax()) << "Node 0 Z coordinate above domain maximum";

            // Test individual coordinate access functions
            auto x0_individual = domain.sfcToPhysicalCoordinateX(sfc0);
            auto y0_individual = domain.sfcToPhysicalCoordinateY(sfc0);
            auto z0_individual = domain.sfcToPhysicalCoordinateZ(sfc0);

            EXPECT_FLOAT_EQ(x0, x0_individual) << "X coordinate mismatch between tuple and individual access";
            EXPECT_FLOAT_EQ(y0, y0_individual) << "Y coordinate mismatch between tuple and individual access";
            EXPECT_FLOAT_EQ(z0, z0_individual) << "Z coordinate mismatch between tuple and individual access";

            // Test spatial coordinate conversion
            auto [ix0, iy0, iz0] = domain.sfcToSpatialCoordinate(sfc0);
            EXPECT_GE(ix0, 0u) << "Spatial coordinate should be non-negative";
            EXPECT_GE(iy0, 0u) << "Spatial coordinate should be non-negative";
            EXPECT_GE(iz0, 0u) << "Spatial coordinate should be non-negative";
        }
    }

    // Helper to validate SFC connectivity storage
    template<typename DomainType>
    void validateSfcConnectivity(const DomainType& domain, size_t maxSamples = 5)
    {
        if (domain.getElementCount() == 0) return;

        // Direct SFC connectivity access (stored data)
        auto h_sfc_i0 = toHost(domain.template indices<0>());
        auto h_sfc_i1 = toHost(domain.template indices<1>());
        auto h_sfc_i2 = toHost(domain.template indices<2>());
        auto h_sfc_i3 = toHost(domain.template indices<3>());

        EXPECT_EQ(h_sfc_i0.size(), domain.getElementCount());
        EXPECT_EQ(h_sfc_i1.size(), domain.getElementCount());
        EXPECT_EQ(h_sfc_i2.size(), domain.getElementCount());
        EXPECT_EQ(h_sfc_i3.size(), domain.getElementCount());

        // Check first few elements for valid SFC codes
        size_t samplesToCheck = std::min(domain.getElementCount(), maxSamples);
        for (size_t i = 0; i < samplesToCheck; i++)
        {
            // SFC codes should be valid (non-negative)
            EXPECT_GE(h_sfc_i0[i], 0u) << "Invalid SFC code in element " << i << " node 0";
            EXPECT_GE(h_sfc_i1[i], 0u) << "Invalid SFC code in element " << i << " node 1";
            EXPECT_GE(h_sfc_i2[i], 0u) << "Invalid SFC code in element " << i << " node 2";
            EXPECT_GE(h_sfc_i3[i], 0u) << "Invalid SFC code in element " << i << " node 3";

            // Test individual SFC access
            auto sfc0_individual = domain.template getConnectivity<0>(i);
            auto sfc1_individual = domain.template getConnectivity<1>(i);
            auto sfc2_individual = domain.template getConnectivity<2>(i);
            auto sfc3_individual = domain.template getConnectivity<3>(i);

            EXPECT_EQ(h_sfc_i0[i], sfc0_individual) << "SFC connectivity mismatch for element " << i;
            EXPECT_EQ(h_sfc_i1[i], sfc1_individual) << "SFC connectivity mismatch for element " << i;
            EXPECT_EQ(h_sfc_i2[i], sfc2_individual) << "SFC connectivity mismatch for element " << i;
            EXPECT_EQ(h_sfc_i3[i], sfc3_individual) << "SFC connectivity mismatch for element " << i;
        }
    }

    // Helper to print domain statistics across ranks
    template<typename DomainType>
    void printDomainStatistics(const DomainType& domain, const std::string& testName)
    {
        size_t localElementCount = domain.getElementCount();
        size_t localNodeCount = domain.getNodeCount();
        size_t totalElementCount = 0;
        size_t totalNodeCount = 0;

        MPI_Allreduce(&localElementCount, &totalElementCount, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&localNodeCount, &totalNodeCount, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

        if (rank == 0)
        {
            std::cout << "[" << testName << "] Domain statistics:" << std::endl;
            std::cout << "  Total across all ranks: " << totalNodeCount << " nodes, " 
                      << totalElementCount << " elements" << std::endl;
            std::cout << "  Rank " << rank << ": " << localNodeCount << " nodes, " 
                      << localElementCount << " elements" << std::endl;
        }
    }
};

// Basic test for external mesh domain creation with new SFC-based architecture
TEST_F(ExternalMeshDomainTest, BasicSfcDomainCreation)
{
    checkPrerequisites("BasicSfcDomainCreation");

    try
    {
        // Create domain with tetrahedral elements
        using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
        Domain domain(meshPath, rank, numRanks);

        // Use helper methods for validation
        validateBasicDomainProperties(domain);
        validateBoundingBox(domain);
        validateSfcCoordinateConversion(domain);
        validateSfcConnectivity(domain);
        printDomainStatistics(domain, "BasicSfcDomainCreation");
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception creating SFC domain: " << e.what();
    }
}

// Test SFC-based coordinate conversion consistency
TEST_F(ExternalMeshDomainTest, SfcCoordinateConsistency)
{
    checkPrerequisites("SfcCoordinateConsistency");

    try
    {
        using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
        Domain domain(meshPath, rank, numRanks);

        if (domain.getElementCount() == 0) {
            GTEST_SKIP() << "No elements on this rank";
        }

        // Test first element
        size_t elemIdx = 0;
        auto sfc0 = domain.getConnectivity<0>(elemIdx);
        auto sfc1 = domain.getConnectivity<1>(elemIdx);
        
        // Test consistency across multiple calls
        auto [x1, y1, z1] = domain.sfcToPhysicalCoordinate(sfc0);
        auto [x2, y2, z2] = domain.sfcToPhysicalCoordinate(sfc0);
        
        EXPECT_FLOAT_EQ(x1, x2) << "SFC coordinate conversion should be deterministic";
        EXPECT_FLOAT_EQ(y1, y2) << "SFC coordinate conversion should be deterministic";
        EXPECT_FLOAT_EQ(z1, z2) << "SFC coordinate conversion should be deterministic";

        // Test individual vs tuple access consistency
        auto x_individual = domain.sfcToPhysicalCoordinateX(sfc0);
        auto y_individual = domain.sfcToPhysicalCoordinateY(sfc0);
        auto z_individual = domain.sfcToPhysicalCoordinateZ(sfc0);

        EXPECT_FLOAT_EQ(x1, x_individual) << "X coordinate should match between tuple and individual access";
        EXPECT_FLOAT_EQ(y1, y_individual) << "Y coordinate should match between tuple and individual access";
        EXPECT_FLOAT_EQ(z1, z_individual) << "Z coordinate should match between tuple and individual access";

        // Test spatial coordinate consistency
        auto [ix1, iy1, iz1] = domain.sfcToSpatialCoordinate(sfc0);
        auto [ix2, iy2, iz2] = domain.sfcToSpatialCoordinate(sfc0);
        
        EXPECT_EQ(ix1, ix2) << "Spatial coordinates should be deterministic";
        EXPECT_EQ(iy1, iy2) << "Spatial coordinates should be deterministic";
        EXPECT_EQ(iz1, iz2) << "Spatial coordinates should be deterministic";

        if (rank == 0) {
            std::cout << "SFC coordinate conversion consistency validated" << std::endl;
        }
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception in SFC consistency test: " << e.what();
    }
}

// Test different precisions with SFC architecture
TEST_F(ExternalMeshDomainTest, SfcDifferentPrecisions)
{
    checkPrerequisites("SfcDifferentPrecisions");

    try
    {
        // Test with float precision and uint32 SFC keys
        {
            using FloatDomain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
            FloatDomain domain(meshPath, rank, numRanks);
            
            validateBasicDomainProperties(domain);
            validateSfcCoordinateConversion(domain, 5);
            
            if (rank == 0) {
                std::cout << "Float precision with uint32 SFC keys validated" << std::endl;
            }
        }

        // Test with double precision and uint64 SFC keys if x.double exists
        if (fs::exists(fs::path(meshPath) / "x.double"))
        {
            using DoubleDomain = ElementDomain<TetTag, double, uint64_t, cstone::GpuTag>;
            DoubleDomain domain(meshPath, rank, numRanks);
            
            validateBasicDomainProperties(domain);
            validateSfcCoordinateConversion(domain, 5);
            
            if (rank == 0) {
                std::cout << "Double precision with uint64 SFC keys validated" << std::endl;
            }
        }
        else if (rank == 0)
        {
            std::cout << "Skipping double precision test - x.double not found" << std::endl;
        }
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception in SFC precision test: " << e.what();
    }
}

// Test SFC connectivity access patterns and validation
TEST_F(ExternalMeshDomainTest, SfcConnectivityValidation)
{
    checkPrerequisites("SfcConnectivityValidation");

    try
    {
        using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
        Domain domain(meshPath, rank, numRanks);

        if (domain.getElementCount() == 0) {
            GTEST_SKIP() << "No elements on this rank";
        }

        // Test 1: Validate that bulk and individual access return same values
        auto h_sfc_i0 = toHost(domain.indices<0>());
        auto h_sfc_i1 = toHost(domain.indices<1>());
        auto h_sfc_i2 = toHost(domain.indices<2>());
        auto h_sfc_i3 = toHost(domain.indices<3>());

        size_t samplesToCheck = std::min(domain.getElementCount(), size_t(10));
        for (size_t i = 0; i < samplesToCheck; i++) {
            auto sfc0_individual = domain.getConnectivity<0>(i);
            auto sfc1_individual = domain.getConnectivity<1>(i);
            auto sfc2_individual = domain.getConnectivity<2>(i);
            auto sfc3_individual = domain.getConnectivity<3>(i);

            EXPECT_EQ(h_sfc_i0[i], sfc0_individual) << "Bulk vs individual SFC access mismatch for element " << i << " node 0";
            EXPECT_EQ(h_sfc_i1[i], sfc1_individual) << "Bulk vs individual SFC access mismatch for element " << i << " node 1";
            EXPECT_EQ(h_sfc_i2[i], sfc2_individual) << "Bulk vs individual SFC access mismatch for element " << i << " node 2";
            EXPECT_EQ(h_sfc_i3[i], sfc3_individual) << "Bulk vs individual SFC access mismatch for element " << i << " node 3";
        }

        // Test 2: Validate SFC keys are reasonable (non-zero for real meshes)
        bool foundNonZeroSfc = false;
        for (size_t i = 0; i < samplesToCheck; i++) {
            if (h_sfc_i0[i] > 0 || h_sfc_i1[i] > 0 || h_sfc_i2[i] > 0 || h_sfc_i3[i] > 0) {
                foundNonZeroSfc = true;
                break;
            }
        }
        EXPECT_TRUE(foundNonZeroSfc) << "Expected at least some non-zero SFC keys in real mesh data";

        // Test 3: Validate coordinate conversion produces reasonable results
        for (size_t i = 0; i < samplesToCheck; i++) {
            auto sfc0 = domain.getConnectivity<0>(i);
            auto [x, y, z] = domain.sfcToPhysicalCoordinate(sfc0);
            
            auto box = domain.getDomain().box();
            EXPECT_GE(x, box.xmin()) << "Converted coordinate X should be within domain bounds";
            EXPECT_LE(x, box.xmax()) << "Converted coordinate X should be within domain bounds";
            EXPECT_GE(y, box.ymin()) << "Converted coordinate Y should be within domain bounds";
            EXPECT_LE(y, box.ymax()) << "Converted coordinate Y should be within domain bounds";
            EXPECT_GE(z, box.zmin()) << "Converted coordinate Z should be within domain bounds";
            EXPECT_LE(z, box.zmax()) << "Converted coordinate Z should be within domain bounds";
        }

        // Test 4: Check that different nodes in same element have different coordinates
        if (samplesToCheck > 0) {
            size_t elemIdx = 0;
            auto sfc0 = domain.getConnectivity<0>(elemIdx);
            auto sfc1 = domain.getConnectivity<1>(elemIdx);
            auto sfc2 = domain.getConnectivity<2>(elemIdx);
            auto sfc3 = domain.getConnectivity<3>(elemIdx);

            auto [x0, y0, z0] = domain.sfcToPhysicalCoordinate(sfc0);
            auto [x1, y1, z1] = domain.sfcToPhysicalCoordinate(sfc1);
            auto [x2, y2, z2] = domain.sfcToPhysicalCoordinate(sfc2);
            auto [x3, y3, z3] = domain.sfcToPhysicalCoordinate(sfc3);

            // For a valid tetrahedron, nodes should not all be at the same location
            bool nodesDiffer = (x0 != x1 || y0 != y1 || z0 != z1) ||
                              (x0 != x2 || y0 != y2 || z0 != z2) ||
                              (x0 != x3 || y0 != y3 || z0 != z3);
            
            EXPECT_TRUE(nodesDiffer) << "Element nodes should have different coordinates";
        }

        if (rank == 0) {
            std::cout << "SFC connectivity validation completed successfully" << std::endl;
            std::cout << "  Validated " << samplesToCheck << " elements" << std::endl;
            std::cout << "  Total elements on rank 0: " << domain.getElementCount() << std::endl;
        }

    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception in SFC connectivity validation: " << e.what();
    }
}

// Optional: Real GPU performance test (if you want actual performance measurement)
TEST_F(ExternalMeshDomainTest, GpuSfcCoordinateConversionPerformance)
{
    checkPrerequisites("GpuSfcCoordinateConversionPerformance");

    try
    {
        using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
        Domain domain(meshPath, rank, numRanks);

        if (domain.getElementCount() == 0) {
            GTEST_SKIP() << "No elements on this rank";
        }

        // Only run performance test on reasonable mesh sizes
        if (domain.getElementCount() < 1000) {
            GTEST_SKIP() << "Mesh too small for meaningful performance test";
        }

        // Test GPU-based coordinate conversion performance
        size_t numCoordinates = domain.getElementCount() * 4; // 4 nodes per tetrahedron
        
        // Create device vectors for the test
        using DeviceVector = typename Domain::template DeviceVector<unsigned>;
        using DeviceVectorReal = typename Domain::template DeviceVector<float>;
        
        DeviceVector d_sfcKeys(numCoordinates);
        DeviceVectorReal d_coordinates(numCoordinates * 3); // x, y, z for each coordinate

        // Fill SFC keys from domain connectivity (this is setup, not measured)
        std::vector<unsigned> h_sfcKeys(numCoordinates);
        for (size_t i = 0; i < domain.getElementCount(); ++i) {
            h_sfcKeys[i * 4 + 0] = domain.getConnectivity<0>(i);
            h_sfcKeys[i * 4 + 1] = domain.getConnectivity<1>(i);
            h_sfcKeys[i * 4 + 2] = domain.getConnectivity<2>(i);
            h_sfcKeys[i * 4 + 3] = domain.getConnectivity<3>(i);
        }
        copyToDevice(d_sfcKeys, h_sfcKeys);

        // GPU kernel for coordinate conversion
        auto coordinateKernel = [] __device__ (const unsigned* sfcKeys, float* coords, 
                                              size_t numKeys, cstone::Box<float> box) 
        {
            int tid = blockIdx.x * blockDim.x + threadIdx.x;
            if (tid >= numKeys) return;
            
            auto [ix, iy, iz] = cstone::decodeSfc(sfcKeys[tid]);
            constexpr unsigned maxCoord = (1u << cstone::maxTreeLevel<unsigned>{}) - 1;
            float invMaxCoord = 1.0f / maxCoord;
            
            coords[tid * 3 + 0] = box.xmin() + ix * invMaxCoord * (box.xmax() - box.xmin());
            coords[tid * 3 + 1] = box.ymin() + iy * invMaxCoord * (box.ymax() - box.ymin());
            coords[tid * 3 + 2] = box.zmin() + iz * invMaxCoord * (box.zmax() - box.zmin());
        };

        // Warm up GPU
        int blockSize = 256;
        int numBlocks = (numCoordinates + blockSize - 1) / blockSize;
        coordinateKernel<<<numBlocks, blockSize>>>(
            thrust::raw_pointer_cast(d_sfcKeys.data()),
            thrust::raw_pointer_cast(d_coordinates.data()),
            numCoordinates,
            domain.getDomain().box()
        );
        cudaDeviceSynchronize();

        // Measure performance
        const int iterations = 100;
        cudaEvent_t start, stop;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);
        
        cudaEventRecord(start);
        for (int iter = 0; iter < iterations; iter++) {
            coordinateKernel<<<numBlocks, blockSize>>>(
                thrust::raw_pointer_cast(d_sfcKeys.data()),
                thrust::raw_pointer_cast(d_coordinates.data()),
                numCoordinates,
                domain.getDomain().box()
            );
        }
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        
        float milliseconds = 0;
        cudaEventElapsedTime(&milliseconds, start, stop);
        
        // Validate results (copy back small sample)
        auto h_coords = toHost(d_coordinates);
        auto box = domain.getDomain().box();
        
        for (size_t i = 0; i < std::min(size_t(100), numCoordinates); i++) {
            float x = h_coords[i * 3 + 0];
            float y = h_coords[i * 3 + 1];
            float z = h_coords[i * 3 + 2];
            
            EXPECT_GE(x, box.xmin()) << "GPU converted coordinate should be within bounds";
            EXPECT_LE(x, box.xmax()) << "GPU converted coordinate should be within bounds";
            EXPECT_GE(y, box.ymin()) << "GPU converted coordinate should be within bounds";
            EXPECT_LE(y, box.ymax()) << "GPU converted coordinate should be within bounds";
            EXPECT_GE(z, box.zmin()) << "GPU converted coordinate should be within bounds";
            EXPECT_LE(z, box.zmax()) << "GPU converted coordinate should be within bounds";
        }

        if (rank == 0) {
            float avgTimePerIteration = milliseconds / iterations;
            float coordsPerSecond = (numCoordinates * 1000.0f) / avgTimePerIteration;
            
            std::cout << "GPU SFC coordinate conversion performance:" << std::endl;
            std::cout << "  Converted " << numCoordinates << " coordinates" << std::endl;
            std::cout << "  Average time per iteration: " << avgTimePerIteration << " ms" << std::endl;
            std::cout << "  Throughput: " << coordsPerSecond / 1e6 << " million coordinates/second" << std::endl;
        }

        cudaEventDestroy(start);
        cudaEventDestroy(stop);

    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception in GPU performance test: " << e.what();
    }
}

// Test edge cases and error handling for SFC operations
TEST_F(ExternalMeshDomainTest, SfcEdgeCasesAndErrorHandling)
{
    checkPrerequisites("SfcEdgeCasesAndErrorHandling");

    try
    {
        using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
        Domain domain(meshPath, rank, numRanks);

        if (domain.getElementCount() > 0) {
            // Test boundary connectivity access
            auto sfc0 = domain.getConnectivity<0>(domain.getElementCount() - 1);
            auto sfc1 = domain.getConnectivity<1>(domain.getElementCount() - 1);
            auto sfc2 = domain.getConnectivity<2>(domain.getElementCount() - 1);
            auto sfc3 = domain.getConnectivity<3>(domain.getElementCount() - 1);
            
            // All SFC codes should be valid
            EXPECT_GE(sfc0, 0u) << "Valid SFC code should be non-negative";
            EXPECT_GE(sfc1, 0u) << "Valid SFC code should be non-negative";
            EXPECT_GE(sfc2, 0u) << "Valid SFC code should be non-negative";
            EXPECT_GE(sfc3, 0u) << "Valid SFC code should be non-negative";

            // Test coordinate conversion for boundary element
            auto [x0, y0, z0] = domain.sfcToPhysicalCoordinate(sfc0);
            auto box = domain.getDomain().box();
            
            EXPECT_GE(x0, box.xmin()) << "Boundary coordinate should be within domain";
            EXPECT_LE(x0, box.xmax()) << "Boundary coordinate should be within domain";
            
            // Test out-of-bounds connectivity access (should return default values)
            auto oob_sfc = domain.getConnectivity<0>(domain.getElementCount() + 1000);
            EXPECT_EQ(oob_sfc, 0u) << "Out-of-bounds SFC access should return 0";

            // Test coordinate conversion with zero SFC key
            auto [zero_x, zero_y, zero_z] = domain.sfcToPhysicalCoordinate(0u);
            // Should return valid coordinates (bottom-left corner of domain)
            EXPECT_GE(zero_x, box.xmin()) << "Zero SFC should map to valid coordinate";
            EXPECT_GE(zero_y, box.ymin()) << "Zero SFC should map to valid coordinate";
            EXPECT_GE(zero_z, box.zmin()) << "Zero SFC should map to valid coordinate";
        }

        if (rank == 0) {
            std::cout << "SFC edge cases and error handling validated" << std::endl;
        }
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception in SFC edge case test: " << e.what();
    }
}

// Test FEM-style element iteration using SFC coordinates
TEST_F(ExternalMeshDomainTest, FemStyleElementIteration)
{
    checkPrerequisites("FemStyleElementIteration");

    try
    {
        using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
        Domain domain(meshPath, rank, numRanks);

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

// Main function
int main(int argc, char** argv)
{
    // Initialize MPI through Mars environment
    mars::Env env(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    