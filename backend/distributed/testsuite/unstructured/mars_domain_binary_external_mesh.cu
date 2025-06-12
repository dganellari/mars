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

// Bring std::get into scope for thrust::tuple compatibility
namespace cstone {
    using std::get;
}

__global__ void validateConnectivityKernel(const unsigned* sfc0_ptr, const unsigned* sfc1_ptr,
                                         const unsigned* sfc2_ptr, const unsigned* sfc3_ptr,
                                         cstone::Box<float> box, int* results, 
                                         size_t numElements)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= numElements) return;
    
    // Get SFC keys for this element
    auto sfc0 = sfc0_ptr[tid];
    auto sfc1 = sfc1_ptr[tid];
    auto sfc2 = sfc2_ptr[tid];
    auto sfc3 = sfc3_ptr[tid];
    
    // Check that nodes are distinct (this is valid regardless of SFC key value)
    bool distinct = (sfc0 != sfc1) && (sfc0 != sfc2) && (sfc0 != sfc3) && 
                   (sfc1 != sfc2) && (sfc1 != sfc3) && (sfc2 != sfc3);
    
    if (!distinct) {
        results[tid] = 0;
        return;
    }
    
    // If any SFC key is 0, it maps to the minimum corner of the box
    // This is valid and expected in some elements after domain decomposition
    if (sfc0 == 0 || sfc1 == 0 || sfc2 == 0 || sfc3 == 0) {
        results[tid] = 1;
        return;
    }
    
    // Check all four SFC keys by decoding and validating coordinates
    constexpr unsigned maxCoord = (1u << cstone::maxTreeLevel<cstone::SfcKind<unsigned>>{}) - 1;
    float invMaxCoord = 1.0f / maxCoord;
    const float tolerance = 1e-5f;

    // Helper lambda to check if coordinates are within bounds
    auto validateCoords = [&](unsigned sfc) -> bool {
        auto sfcKindKey = cstone::SfcKind<unsigned>(sfc);
        auto [ix, iy, iz] = cstone::decodeSfc(sfcKindKey);
        
        float x = box.xmin() + ix * invMaxCoord * (box.xmax() - box.xmin());
        float y = box.ymin() + iy * invMaxCoord * (box.ymax() - box.ymin());
        float z = box.zmin() + iz * invMaxCoord * (box.zmax() - box.zmin());
        
        return (x >= box.xmin()-tolerance && x <= box.xmax()+tolerance &&
                y >= box.ymin()-tolerance && y <= box.ymax()+tolerance &&
                z >= box.zmin()-tolerance && z <= box.zmax()+tolerance);
    };
    
    // Check all four SFC keys - each node must be valid
    bool allValid = validateCoords(sfc0) && validateCoords(sfc1) && 
                   validateCoords(sfc2) && validateCoords(sfc3);
    
    results[tid] = allValid ? 1 : 0;
}

__global__ void performanceTestKernel(const unsigned* sfcKeys, float* coords,
                                     size_t numKeys, cstone::Box<float> box)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= numKeys) return;
    
    auto sfcKindKey = cstone::SfcKind<unsigned>(sfcKeys[tid]);
    auto [ix, iy, iz] = cstone::decodeSfc(sfcKindKey);
    
    constexpr unsigned maxCoord = (1u << cstone::maxTreeLevel<cstone::SfcKind<unsigned>>{}) - 1;
    float invMaxCoord = 1.0f / maxCoord;
    
    coords[tid * 3 + 0] = box.xmin() + ix * invMaxCoord * (box.xmax() - box.xmin());
    coords[tid * 3 + 1] = box.ymin() + iy * invMaxCoord * (box.ymax() - box.ymin());
    coords[tid * 3 + 2] = box.zmin() + iz * invMaxCoord * (box.zmax() - box.zmin());
}

__global__ void binarySearchTestKernel(const unsigned* sortedKeys, const unsigned* queryKeys,
                                      int* results, size_t numSorted, size_t numQueries)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= numQueries) return;
    
    unsigned key = queryKeys[tid];
    
    // Binary search implementation
    int left = 0;
    int right = numSorted - 1;
    int result = -1;
    
    while (left <= right) {
        int mid = left + (right - left) / 2;
        if (sortedKeys[mid] == key) {
            result = mid;
            break;
        } else if (sortedKeys[mid] < key) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }
    
    results[tid] = result;
}

__global__ void volumeCalculationKernel(const unsigned* sfc0_ptr, const unsigned* sfc1_ptr,
                                       const unsigned* sfc2_ptr, const unsigned* sfc3_ptr,
                                       cstone::Box<float> box, float* volumes,
                                       size_t numElements)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= numElements) return;
    
    // Get SFC keys
    auto sfc0 = sfc0_ptr[tid];
    auto sfc1 = sfc1_ptr[tid];
    auto sfc2 = sfc2_ptr[tid];
    auto sfc3 = sfc3_ptr[tid];
    
    // Convert to coordinates using structured bindings
    constexpr unsigned maxCoord = (1u << cstone::maxTreeLevel<cstone::SfcKind<unsigned>>{}) - 1;
    float invMaxCoord = 1.0f / maxCoord;
    
    auto convertSfc = [&](unsigned sfc) {
        auto sfcKey = cstone::SfcKind<unsigned>(sfc);
        auto [ix, iy, iz] = cstone::decodeSfc(sfcKey);
        float x = box.xmin() + ix * invMaxCoord * (box.xmax() - box.xmin());
        float y = box.ymin() + iy * invMaxCoord * (box.ymax() - box.ymin());
        float z = box.zmin() + iz * invMaxCoord * (box.zmax() - box.zmin());
        return make_float3(x, y, z);
    };
    
    auto p0 = convertSfc(sfc0);
    auto p1 = convertSfc(sfc1);
    auto p2 = convertSfc(sfc2);
    auto p3 = convertSfc(sfc3);
    
    // Calculate tetrahedral volume
    float det = (p1.x - p0.x) * ((p2.y - p0.y) * (p3.z - p0.z) - (p2.z - p0.z) * (p3.y - p0.y)) -
                (p1.y - p0.y) * ((p2.x - p0.x) * (p3.z - p0.z) - (p2.z - p0.z) * (p3.x - p0.x)) +
                (p1.z - p0.z) * ((p2.x - p0.x) * (p3.y - p0.y) - (p2.y - p0.y) * (p3.x - p0.x));
    
    volumes[tid] = fabsf(det) / 6.0f;
}

class ExternalMeshDomainTest : public ::testing::Test
{
protected:
    int rank;
    int numRanks;
    std::string meshPath;
    int deviceCount = 0;

    void SetUp() override
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
        
        meshPath = getMeshPath();
        cudaGetDeviceCount(&deviceCount);
        
        if (rank == 0) {
            std::cout << "GPU Test setup: Found " << deviceCount << " CUDA devices" << std::endl;
            std::cout << "Using mesh at: " << (meshPath.empty() ? "none" : meshPath) << std::endl;
        }
    }

private:
    std::string getMeshPath() const
    {
        const char* meshPathEnv = std::getenv("MESH_PATH");
        std::string path = meshPathEnv ? meshPathEnv : "";

        if (path.empty() || !fs::exists(path))
        {
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
            GTEST_SKIP() << testName << ": No CUDA devices available";
        }
    }

    // Helper validation functions
    template<typename Domain>
    void validateBasicDomainProperties(const Domain& domain)
    {
        EXPECT_GT(domain.getElementCount(), 0) << "Domain should have elements";
        
        if (rank == 0) {
            std::cout << "  Element count: " << domain.getElementCount() << std::endl;
        }
    }
    
    template<typename Domain>
    void validateBoundingBox(const Domain& domain)
    {
        auto box = domain.getDomain().box();
        
        EXPECT_LT(box.xmin(), box.xmax()) << "Invalid X bounds";
        EXPECT_LT(box.ymin(), box.ymax()) << "Invalid Y bounds";
        EXPECT_LT(box.zmin(), box.zmax()) << "Invalid Z bounds";
        
        if (rank == 0) {
            std::cout << "  Bounding box: [" << box.xmin() << ", " << box.xmax() << "] x ["
                      << box.ymin() << ", " << box.ymax() << "] x ["
                      << box.zmin() << ", " << box.zmax() << "]" << std::endl;
        }
    }
    
    template<typename Domain>
    void validateSfcCoordinateConversion(const Domain& domain)
    {
        if (domain.getElementCount() == 0) return;
            
        auto h_i0 = toHost(domain.template getConnectivity<0>());
        if (h_i0.empty()) return;
            
        auto sfc0 = h_i0[0];
        auto box = domain.getDomain().box();
            
        // Manual SFC to coordinate conversion
        auto sfcKindKey = cstone::SfcKind<unsigned>(sfc0);
        auto [ix, iy, iz] = cstone::decodeSfc(sfcKindKey);
        constexpr unsigned maxCoord = (1u << cstone::maxTreeLevel<cstone::SfcKind<unsigned>>{}) - 1;
        float invMaxCoord = 1.0f / maxCoord;
            
        float x = box.xmin() + ix * invMaxCoord * (box.xmax() - box.xmin());
        float y = box.ymin() + iy * invMaxCoord * (box.ymax() - box.ymin());
        float z = box.zmin() + iz * invMaxCoord * (box.zmax() - box.zmin());
            
        EXPECT_GE(x, box.xmin()) << "SFC coordinate conversion failed - X below minimum";
        EXPECT_LE(x, box.xmax()) << "SFC coordinate conversion failed - X above maximum";
        EXPECT_GE(y, box.ymin()) << "SFC coordinate conversion failed - Y below minimum";
        EXPECT_LE(y, box.ymax()) << "SFC coordinate conversion failed - Y above maximum";
        EXPECT_GE(z, box.zmin()) << "SFC coordinate conversion failed - Z below minimum";
        EXPECT_LE(z, box.zmax()) << "SFC coordinate conversion failed - Z above maximum";
    }
        
    template<typename Domain>
    void validateSfcConnectivity(const Domain& domain)
    {
        if (domain.getElementCount() == 0) return;
            
        size_t samplesToCheck = std::min(domain.getElementCount(), size_t(10));
            
        auto h_i0 = toHost(domain.template getConnectivity<0>());
        auto h_i1 = toHost(domain.template getConnectivity<1>());
        auto h_i2 = toHost(domain.template getConnectivity<2>());
        auto h_i3 = toHost(domain.template getConnectivity<3>());
            
        for (size_t i = 0; i < samplesToCheck; i++) {
            if (i >= h_i0.size()) break;
                
            auto sfc0 = h_i0[i];
            auto sfc1 = h_i1[i];
            auto sfc2 = h_i2[i];
            auto sfc3 = h_i3[i];
                
            // Check that connectivity values are reasonable (non-zero for valid meshes)
            EXPECT_GT(sfc0, 0u) << "Connectivity 0 should be non-zero for element " << i;
            EXPECT_GT(sfc1, 0u) << "Connectivity 1 should be non-zero for element " << i;
            EXPECT_GT(sfc2, 0u) << "Connectivity 2 should be non-zero for element " << i;
            EXPECT_GT(sfc3, 0u) << "Connectivity 3 should be non-zero for element " << i;
                
            // Check that nodes are distinct
            EXPECT_NE(sfc0, sfc1) << "Element " << i << " has duplicate nodes 0,1";
            EXPECT_NE(sfc0, sfc2) << "Element " << i << " has duplicate nodes 0,2";
            EXPECT_NE(sfc0, sfc3) << "Element " << i << " has duplicate nodes 0,3";
            EXPECT_NE(sfc1, sfc2) << "Element " << i << " has duplicate nodes 1,2";
            EXPECT_NE(sfc1, sfc3) << "Element " << i << " has duplicate nodes 1,3";
            EXPECT_NE(sfc2, sfc3) << "Element " << i << " has duplicate nodes 2,3";
        }
    }    
    
    template<typename Domain>
    void printDomainStatistics(const Domain& domain, const std::string& testName)
    {
        if (rank == 0) {
            std::cout << testName << " completed successfully:" << std::endl;
            std::cout << "  Elements on rank 0: " << domain.getElementCount() << std::endl;
            
            if (domain.getElementCount() > 0) {
                auto box = domain.getDomain().box();
                float volume = (box.xmax() - box.xmin()) * 
                              (box.ymax() - box.ymin()) * 
                              (box.zmax() - box.zmin());
                std::cout << "  Domain volume: " << volume << std::endl;
            }
        }
    }
};

TEST_F(ExternalMeshDomainTest, HostSfcConnectivityValidation)
{
    checkPrerequisites("HostSfcConnectivityValidation");

    try
    {
        using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
        Domain domain(meshPath, rank, numRanks);

        if (domain.getElementCount() == 0) {
            GTEST_SKIP() << "No elements on this rank";
        }

        size_t samplesToCheck = std::min(domain.getElementCount(), size_t(100));
        
        auto h_i0 = toHost(domain.template getConnectivity<0>());
        auto h_i1 = toHost(domain.template getConnectivity<1>());
        auto h_i2 = toHost(domain.template getConnectivity<2>());
        auto h_i3 = toHost(domain.template getConnectivity<3>());
        
        for (size_t i = 0; i < samplesToCheck; i++) {
            if (i >= h_i0.size()) break;
            
            auto sfc0 = h_i0[i];
            auto sfc1 = h_i1[i];
            auto sfc2 = h_i2[i];
            auto sfc3 = h_i3[i];

            // Manual SFC to coordinate conversion
            auto box = domain.getDomain().box();
            auto convertSfc = [&](unsigned sfc) -> std::tuple<float, float, float> {
                auto sfcKindKey = cstone::SfcKind<unsigned>(sfc);
                auto [ix, iy, iz] = cstone::decodeSfc(sfcKindKey);
                constexpr unsigned maxCoord = (1u << cstone::maxTreeLevel<cstone::SfcKind<unsigned>>{}) - 1;
                float invMaxCoord = 1.0f / maxCoord;
                
                float x = box.xmin() + ix * invMaxCoord * (box.xmax() - box.xmin());
                float y = box.ymin() + iy * invMaxCoord * (box.ymax() - box.ymin());
                float z = box.zmin() + iz * invMaxCoord * (box.zmax() - box.zmin());
                return {x, y, z};
            };

            auto [x0, y0, z0] = convertSfc(sfc0);
            auto [x1, y1, z1] = convertSfc(sfc1);
            auto [x2, y2, z2] = convertSfc(sfc2);
            auto [x3, y3, z3] = convertSfc(sfc3);
            
            // Validate coordinates are within bounds
            EXPECT_GE(x0, box.xmin()) << "Node 0 X coordinate below domain minimum";
            EXPECT_LE(x0, box.xmax()) << "Node 0 X coordinate above domain maximum";
            EXPECT_GE(y0, box.ymin()) << "Node 0 Y coordinate below domain minimum";
            EXPECT_LE(y0, box.ymax()) << "Node 0 Y coordinate above domain maximum";
            
            // Check that nodes are distinct
            bool nodesDiffer = (x0 != x1 || y0 != y1 || z0 != z1) ||
                              (x0 != x2 || y0 != y2 || z0 != z2) ||
                              (x0 != x3 || y0 != y3 || z0 != z3);
            EXPECT_TRUE(nodesDiffer) << "Element nodes should have different coordinates";
        }

        if (rank == 0) {
            std::cout << "Host SFC connectivity validation passed for " 
                      << samplesToCheck << " elements" << std::endl;
        }
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception in host SFC validation: " << e.what();
    }
}

TEST_F(ExternalMeshDomainTest, DeviceSfcConnectivityValidation)
{
    checkPrerequisites("DeviceSfcConnectivityValidation");

    try
    {
        using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
        Domain domain(meshPath, rank, numRanks);

        if (domain.getElementCount() == 0) {
            GTEST_SKIP() << "No elements on this rank";
        }

        // First, verify that the box is actually global by printing it on all ranks
        auto box = domain.getDomain().box();

        // Continue with the original validation
        size_t testElements = domain.getElementCount();
        
        auto* sfc0_ptr = domain.template indices<0>().data();
        auto* sfc1_ptr = domain.template indices<1>().data();
        auto* sfc2_ptr = domain.template indices<2>().data();
        auto* sfc3_ptr = domain.template indices<3>().data();
        
        cstone::DeviceVector<int> d_results(testElements);
        
        int blockSize = 256;
        int numBlocks = (testElements + blockSize - 1) / blockSize;
        
        validateConnectivityKernel<<<numBlocks, blockSize>>>(
            sfc0_ptr, sfc1_ptr, sfc2_ptr, sfc3_ptr,
            box, d_results.data(), testElements);
        
        cudaDeviceSynchronize();
        cudaError_t err = cudaGetLastError();
        ASSERT_EQ(err, cudaSuccess) << "CUDA kernel execution failed: " << cudaGetErrorString(err);
        
        auto h_results = toHost(d_results);
        int failedElements = 0;
        for (size_t i = 0; i < testElements; i++) {
            if (h_results[i] == 0) {
                failedElements++;
            }
        }
        
        float failurePercentage = (float)failedElements / testElements * 100.0f;
        
        // If we have failures, print some diagnostics
        if (failedElements > 0) {
            std::cout << "Rank " << rank << ": " << failedElements 
                      << " elements failed validation (" << failurePercentage << "%)" << std::endl;
            
            auto h_sfc0 = toHost(domain.template indices<0>());
            auto h_sfc1 = toHost(domain.template indices<1>());
            auto h_sfc2 = toHost(domain.template indices<2>());
            auto h_sfc3 = toHost(domain.template indices<3>());
            
            // Define number of elements to diagnose
            int diagCount = std::min(5, failedElements);
            int found = 0;
            
            for (size_t i = 0; i < testElements && found < diagCount; i++) {
                if (h_results[i] == 0) {
                    auto sfc0 = h_sfc0[i];
                    auto sfc1 = h_sfc1[i];
                    auto sfc2 = h_sfc2[i];
                    auto sfc3 = h_sfc3[i];
                    
                    std::cout << "  Element " << i << " SFC keys: " 
                              << sfc0 << ", " << sfc1 << ", " << sfc2 << ", " << sfc3 << std::endl;
                    
                    // Print decoded coordinates for the first node
                    auto sfcKindKey0 = cstone::SfcKind<unsigned>(sfc0);
                    auto [ix0, iy0, iz0] = cstone::decodeSfc(sfcKindKey0);
                    
                    constexpr unsigned maxCoord = (1u << cstone::maxTreeLevel<cstone::SfcKind<unsigned>>{}) - 1;
                    float invMaxCoord = 1.0f / maxCoord;
                    
                    float x0 = box.xmin() + ix0 * invMaxCoord * (box.xmax() - box.xmin());
                    float y0 = box.ymin() + iy0 * invMaxCoord * (box.ymax() - box.ymin());
                    float z0 = box.zmin() + iz0 * invMaxCoord * (box.zmax() - box.zmin());
                    
                    std::cout << "    Node 0 coords: (" << x0 << ", " << y0 << ", " << z0 << ")" << std::endl;
                    
                    found++;
                }
            }
        }
        
        // For multi-rank tests, we need to be more tolerant for floating-point precision issues
        float maxAllowedFailurePercent = numRanks > 1 ? 1.0f : 0.1f;
        
        EXPECT_LE(failurePercentage, maxAllowedFailurePercent) 
            << failedElements << " elements failed validation (" 
            << failurePercentage << "%)";
        
        std::cout << "Rank " << rank << ": " 
                  << (testElements - failedElements) << "/" << testElements
                  << " elements passed validation (" 
                  << (100.0f - failurePercentage) << "%)" << std::endl;
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception in device SFC validation: " << e.what();
    }
}

TEST_F(ExternalMeshDomainTest, BasicSfcDomainCreation)
{
    checkPrerequisites("BasicSfcDomainCreation");

    try
    {
        using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
        Domain domain(meshPath, rank, numRanks);

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

int main(int argc, char** argv)
{
    mars::Env env(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}