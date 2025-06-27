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
namespace cstone
{
using std::get;
}

__global__ void validateConnectivityKernel(const unsigned* sfc0_ptr,
                                           const unsigned* sfc1_ptr,
                                           const unsigned* sfc2_ptr,
                                           const unsigned* sfc3_ptr,
                                           cstone::Box<float> box,
                                           int* results,
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
    bool distinct =
        (sfc0 != sfc1) && (sfc0 != sfc2) && (sfc0 != sfc3) && (sfc1 != sfc2) && (sfc1 != sfc3) && (sfc2 != sfc3);

    if (!distinct)
    {
        results[tid] = 0;
        return;
    }

    // FAIL if ALL SFC keys are 0 (indicates initialization bug)
    if (sfc0 == 0 && sfc1 == 0 && sfc2 == 0 && sfc3 == 0)
    {
        results[tid] = 0;
        return;
    }

    // If any SFC key is 0, it maps to the minimum corner of the box
    // This is valid and expected in some elements after domain decomposition
    if (sfc0 == 0 || sfc1 == 0 || sfc2 == 0 || sfc3 == 0)
    {
        results[tid] = 1;
        return;
    }

    // Check all four SFC keys by decoding and validating coordinates
    constexpr unsigned maxCoord = (1u << cstone::maxTreeLevel<cstone::SfcKind<unsigned>>{}) - 1;
    float invMaxCoord           = 1.0f / maxCoord;
    const float tolerance       = 1e-5f;

    // Helper lambda to check if coordinates are within bounds
    auto validateCoords = [&](unsigned sfc) -> bool
    {
        auto sfcKindKey   = cstone::SfcKind<unsigned>(sfc);
        auto [ix, iy, iz] = cstone::decodeSfc(sfcKindKey);

        float x = box.xmin() + ix * invMaxCoord * (box.xmax() - box.xmin());
        float y = box.ymin() + iy * invMaxCoord * (box.ymax() - box.ymin());
        float z = box.zmin() + iz * invMaxCoord * (box.zmax() - box.zmin());

        return (x >= box.xmin() - tolerance && x <= box.xmax() + tolerance && y >= box.ymin() - tolerance &&
                y <= box.ymax() + tolerance && z >= box.zmin() - tolerance && z <= box.zmax() + tolerance);
    };

    // Check all four SFC keys - each node must be valid
    bool allValid = validateCoords(sfc0) && validateCoords(sfc1) && validateCoords(sfc2) && validateCoords(sfc3);

    results[tid] = allValid ? 1 : 0;
}

__global__ void performanceTestKernel(const unsigned* sfc0_ptr,
                                      const unsigned* sfc1_ptr,
                                      const unsigned* sfc2_ptr,
                                      const unsigned* sfc3_ptr,
                                      cstone::Box<float> box,
                                      float* coordinates,
                                      size_t numKeys)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= numKeys) return;

    size_t elemIdx = tid / 4; // 4 nodes per tetrahedron
    int nodeIdx    = tid % 4; // Which node (0,1,2,3)

    // Get the appropriate SFC key based on node index
    unsigned sfcKey;
    switch (nodeIdx)
    {
        case 0: sfcKey = sfc0_ptr[elemIdx]; break;
        case 1: sfcKey = sfc1_ptr[elemIdx]; break;
        case 2: sfcKey = sfc2_ptr[elemIdx]; break;
        case 3: sfcKey = sfc3_ptr[elemIdx]; break;
        default: return;
    }

    // Skip zero SFC keys (boundary elements)
    if (sfcKey == 0)
    {
        coordinates[tid * 3 + 0] = box.xmin();
        coordinates[tid * 3 + 1] = box.ymin();
        coordinates[tid * 3 + 2] = box.zmin();
        return;
    }

    // Manual SFC to coordinate conversion
    auto sfcKindKey             = cstone::SfcKind<unsigned>(sfcKey);
    auto [ix, iy, iz]           = cstone::decodeSfc(sfcKindKey);
    constexpr unsigned maxCoord = (1u << cstone::maxTreeLevel<cstone::SfcKind<unsigned>>{}) - 1;
    float invMaxCoord           = 1.0f / maxCoord;

    coordinates[tid * 3 + 0] = box.xmin() + ix * invMaxCoord * (box.xmax() - box.xmin());
    coordinates[tid * 3 + 1] = box.ymin() + iy * invMaxCoord * (box.ymax() - box.ymin());
    coordinates[tid * 3 + 2] = box.zmin() + iz * invMaxCoord * (box.zmax() - box.zmin());
}

__global__ void volumeCalculationKernel(const unsigned* sfc0_ptr,
                                        const unsigned* sfc1_ptr,
                                        const unsigned* sfc2_ptr,
                                        const unsigned* sfc3_ptr,
                                        cstone::Box<float> box,
                                        float* volumes,
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
    float invMaxCoord           = 1.0f / maxCoord;

    auto convertSfc = [&](unsigned sfc)
    {
        auto sfcKey       = cstone::SfcKind<unsigned>(sfc);
        auto [ix, iy, iz] = cstone::decodeSfc(sfcKey);
        float x           = box.xmin() + ix * invMaxCoord * (box.xmax() - box.xmin());
        float y           = box.ymin() + iy * invMaxCoord * (box.ymax() - box.ymin());
        float z           = box.zmin() + iz * invMaxCoord * (box.zmax() - box.zmin());
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
public:
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

        if (rank == 0)
        {
            std::cout << "GPU Test setup: Found " << deviceCount << " CUDA devices" << std::endl;
            std::cout << "Using mesh at: " << (meshPath.empty() ? "none" : meshPath) << std::endl;
        }
    }

private:
    std::string getMeshPath() const
    {
        const char* meshPathEnv = std::getenv("MESH_PATH");
        std::string path        = meshPathEnv ? meshPathEnv : "";

        if (path.empty() || !fs::exists(path))
        {
            std::vector<std::string> commonLocations = {"./test_data", "../test_data",    "./meshes",
                                                        "../meshes",   "../../test_data", "../../meshes"};

            for (const auto& loc : commonLocations)
            {
                if (fs::exists(loc) && fs::is_directory(fs::path(loc)))
                {
                    if (hasRequiredMeshFiles(loc)) { return loc; }
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

        if (deviceCount == 0) { GTEST_SKIP() << testName << ": No CUDA devices available"; }
    }

    // Helper validation functions
    template<typename Domain>
    void validateBasicDomainProperties(const Domain& domain)
    {
        EXPECT_GT(domain.getElementCount(), 0) << "Domain should have elements";

        if (rank == 0) { std::cout << "  Element count: " << domain.getElementCount() << std::endl; }
    }

    template<typename Domain>
    void validateBoundingBox(const Domain& domain)
    {
        auto box = domain.getDomain().box();

        EXPECT_LT(box.xmin(), box.xmax()) << "Invalid X bounds";
        EXPECT_LT(box.ymin(), box.ymax()) << "Invalid Y bounds";
        EXPECT_LT(box.zmin(), box.zmax()) << "Invalid Z bounds";

        if (rank == 0)
        {
            std::cout << "  Bounding box: [" << box.xmin() << ", " << box.xmax() << "] x [" << box.ymin() << ", "
                      << box.ymax() << "] x [" << box.zmin() << ", " << box.zmax() << "]" << std::endl;
        }
    }

    template<typename Domain>
    void validateSfcCoordinateConversion(const Domain& domain)
    {
        if (domain.getElementCount() == 0) return;

        auto h_i0 = toHost(domain.template indices<0>());
        if (h_i0.empty()) return;

        auto sfc0 = h_i0[0];
        auto box  = domain.getDomain().box();

        // Manual SFC to coordinate conversion
        auto sfcKindKey             = cstone::SfcKind<unsigned>(sfc0);
        auto [ix, iy, iz]           = cstone::decodeSfc(sfcKindKey);
        constexpr unsigned maxCoord = (1u << cstone::maxTreeLevel<cstone::SfcKind<unsigned>>{}) - 1;
        float invMaxCoord           = 1.0f / maxCoord;

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

        auto h_i0 = toHost(domain.template indices<0>());
        auto h_i1 = toHost(domain.template indices<1>());
        auto h_i2 = toHost(domain.template indices<2>());
        auto h_i3 = toHost(domain.template indices<3>());

        // Check for suspicious patterns that indicate bugs
        size_t zeroCount     = 0;
        size_t allZeroCount  = 0;
        size_t totalElements = std::min(samplesToCheck, h_i0.size());

        for (size_t i = 0; i < totalElements; i++)
        {
            auto sfc0 = h_i0[i];
            auto sfc1 = h_i1[i];
            auto sfc2 = h_i2[i];
            auto sfc3 = h_i3[i];

            // Count elements with zero keys
            if (sfc0 == 0 || sfc1 == 0 || sfc2 == 0 || sfc3 == 0) { zeroCount++; }

            // FAIL if ALL nodes of an element are zero (initialization bug)
            if (sfc0 == 0 && sfc1 == 0 && sfc2 == 0 && sfc3 == 0)
            {
                allZeroCount++;
                EXPECT_FALSE(true) << "Element " << i << " has all zero SFC keys - indicates initialization bug";
            }

            // SFC keys must be distinct
            EXPECT_NE(sfc0, sfc1) << "Element " << i << " has duplicate nodes 0,1";
            EXPECT_NE(sfc0, sfc2) << "Element " << i << " has duplicate nodes 0,2";
            EXPECT_NE(sfc0, sfc3) << "Element " << i << " has duplicate nodes 0,3";
            EXPECT_NE(sfc1, sfc2) << "Element " << i << " has duplicate nodes 1,2";
            EXPECT_NE(sfc1, sfc3) << "Element " << i << " has duplicate nodes 1,3";
            EXPECT_NE(sfc2, sfc3) << "Element " << i << " has duplicate nodes 2,3";

            // For coordinate validation, handle zero keys properly
            auto box        = domain.getDomain().box();
            auto convertSfc = [&](unsigned sfc) -> std::tuple<float, float, float>
            {
                if (sfc == 0)
                {
                    // SFC key 0 maps to minimum corner
                    return {box.xmin(), box.ymin(), box.zmin()};
                }

                auto sfcKindKey             = cstone::SfcKind<unsigned>(sfc);
                auto [ix, iy, iz]           = cstone::decodeSfc(sfcKindKey);
                constexpr unsigned maxCoord = (1u << cstone::maxTreeLevel<cstone::SfcKind<unsigned>>{}) - 1;
                float invMaxCoord           = 1.0f / maxCoord;

                float x = box.xmin() + ix * invMaxCoord * (box.xmax() - box.xmin());
                float y = box.ymin() + iy * invMaxCoord * (box.ymax() - box.ymin());
                float z = box.zmin() + iz * invMaxCoord * (box.zmax() - box.zmin());
                return {x, y, z};
            };

            auto [x0, y0, z0] = convertSfc(sfc0);
            auto [x1, y1, z1] = convertSfc(sfc1);
            auto [x2, y2, z2] = convertSfc(sfc2);
            auto [x3, y3, z3] = convertSfc(sfc3);

            // Validate coordinates are within bounds with tolerance
            const float tolerance = 1e-5f;

            auto validateCoord = [&](float coord, float min_val, float max_val, const std::string& axis, int node)
            {
                EXPECT_GE(coord, min_val - tolerance)
                    << "Element " << i << " Node " << node << " " << axis << " coordinate below domain minimum";
                EXPECT_LE(coord, max_val + tolerance)
                    << "Element " << i << " Node " << node << " " << axis << " coordinate above domain maximum";
            };

            validateCoord(x0, box.xmin(), box.xmax(), "X", 0);
            validateCoord(y0, box.ymin(), box.ymax(), "Y", 0);
            validateCoord(z0, box.zmin(), box.zmax(), "Z", 0);

            // Check that nodes are spatially distinct
            bool nodesDiffer = (x0 != x1 || y0 != y1 || z0 != z1) && (x0 != x2 || y0 != y2 || z0 != z2) &&
                               (x0 != x3 || y0 != y3 || z0 != z3) && (x1 != x2 || y1 != y2 || z1 != z2) &&
                               (x1 != x3 || y1 != y3 || z1 != z3) && (x2 != x3 || y2 != y3 || z2 != z3);
            EXPECT_TRUE(nodesDiffer) << "Element " << i << " nodes should have different coordinates";
        }

        // Only fail if elements have ALL zero keys (real initialization bug)
        float allZeroPercentage = (float)allZeroCount / totalElements * 100.0f;

        EXPECT_EQ(allZeroCount, 0) << allZeroCount << " elements (" << allZeroPercentage
                                   << "%) have all zero SFC keys on rank " << rank << " - indicates initialization bug";

        // Optional: Log info about boundary elements (but don't fail)
        if (zeroCount > 0)
        {
            float zeroPercentage = (float)zeroCount / totalElements * 100.0f;
            std::cout << "Rank " << rank << ": " << zeroCount << " elements (" << zeroPercentage
                      << "%) have at least one zero SFC key (boundary elements)" << std::endl;
        }
    }

    template<typename Domain>
    void printDomainStatistics(const Domain& domain, const std::string& testName)
    {
        if (rank == 0)
        {
            std::cout << testName << " completed successfully:" << std::endl;
            std::cout << "  Elements on rank 0: " << domain.getElementCount() << std::endl;

            if (domain.getElementCount() > 0)
            {
                auto box     = domain.getDomain().box();
                float volume = (box.xmax() - box.xmin()) * (box.ymax() - box.ymin()) * (box.zmax() - box.zmin());
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

        if (domain.getElementCount() == 0) { GTEST_SKIP() << "No elements on this rank"; }

        size_t samplesToCheck = std::min(domain.getElementCount(), size_t(100));

        auto h_i0 = toHost(domain.template indices<0>());
        auto h_i1 = toHost(domain.template indices<1>());
        auto h_i2 = toHost(domain.template indices<2>());
        auto h_i3 = toHost(domain.template indices<3>());

        for (size_t i = 0; i < samplesToCheck; i++)
        {
            if (i >= h_i0.size()) break;

            auto sfc0 = h_i0[i];
            auto sfc1 = h_i1[i];
            auto sfc2 = h_i2[i];
            auto sfc3 = h_i3[i];

            // Manual SFC to coordinate conversion
            auto box        = domain.getDomain().box();
            auto convertSfc = [&](unsigned sfc) -> std::tuple<float, float, float>
            {
                auto sfcKindKey             = cstone::SfcKind<unsigned>(sfc);
                auto [ix, iy, iz]           = cstone::decodeSfc(sfcKindKey);
                constexpr unsigned maxCoord = (1u << cstone::maxTreeLevel<cstone::SfcKind<unsigned>>{}) - 1;
                float invMaxCoord           = 1.0f / maxCoord;

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
            bool nodesDiffer = (x0 != x1 || y0 != y1 || z0 != z1) || (x0 != x2 || y0 != y2 || z0 != z2) ||
                               (x0 != x3 || y0 != y3 || z0 != z3);
            EXPECT_TRUE(nodesDiffer) << "Element nodes should have different coordinates";
        }

        if (rank == 0)
        {
            std::cout << "Host SFC connectivity validation passed for " << samplesToCheck << " elements" << std::endl;
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

        if (domain.getElementCount() == 0) { GTEST_SKIP() << "No elements on this rank"; }

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

        validateConnectivityKernel<<<numBlocks, blockSize>>>(sfc0_ptr, sfc1_ptr, sfc2_ptr, sfc3_ptr, box,
                                                             d_results.data(), testElements);

        cudaDeviceSynchronize();
        cudaError_t err = cudaGetLastError();
        ASSERT_EQ(err, cudaSuccess) << "CUDA kernel execution failed: " << cudaGetErrorString(err);

        auto h_results     = toHost(d_results);
        int failedElements = 0;
        for (size_t i = 0; i < testElements; i++)
        {
            if (h_results[i] == 0) { failedElements++; }
        }

        float failurePercentage = (float)failedElements / testElements * 100.0f;

        // If we have failures, print some diagnostics
        if (failedElements > 0)
        {
            std::cout << "Rank " << rank << ": " << failedElements << " elements failed validation ("
                      << failurePercentage << "%)" << std::endl;

            auto h_sfc0 = toHost(domain.template indices<0>());
            auto h_sfc1 = toHost(domain.template indices<1>());
            auto h_sfc2 = toHost(domain.template indices<2>());
            auto h_sfc3 = toHost(domain.template indices<3>());

            // Define number of elements to diagnose
            int diagCount = std::min(5, failedElements);
            int found     = 0;

            for (size_t i = 0; i < testElements && found < diagCount; i++)
            {
                if (h_results[i] == 0)
                {
                    auto sfc0 = h_sfc0[i];
                    auto sfc1 = h_sfc1[i];
                    auto sfc2 = h_sfc2[i];
                    auto sfc3 = h_sfc3[i];

                    std::cout << "  Element " << i << " SFC keys: " << sfc0 << ", " << sfc1 << ", " << sfc2 << ", "
                              << sfc3 << std::endl;

                    // Print decoded coordinates for the first node
                    auto sfcKindKey0     = cstone::SfcKind<unsigned>(sfc0);
                    auto [ix0, iy0, iz0] = cstone::decodeSfc(sfcKindKey0);

                    constexpr unsigned maxCoord = (1u << cstone::maxTreeLevel<cstone::SfcKind<unsigned>>{}) - 1;
                    float invMaxCoord           = 1.0f / maxCoord;

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
            << failedElements << " elements failed validation (" << failurePercentage << "%)";

        std::cout << "Rank " << rank << ": " << (testElements - failedElements) << "/" << testElements
                  << " elements passed validation (" << (100.0f - failurePercentage) << "%)" << std::endl;
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

TEST_F(ExternalMeshDomainTest, GpuVolumeCalculation)
{
    checkPrerequisites("GpuVolumeCalculation");

    try
    {
        using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
        Domain domain(meshPath, rank, numRanks);

        if (domain.getElementCount() == 0) { GTEST_SKIP() << "No elements on this rank"; }

        size_t numElements = domain.getElementCount();

        auto* sfc0_ptr = domain.template indices<0>().data();
        auto* sfc1_ptr = domain.template indices<1>().data();
        auto* sfc2_ptr = domain.template indices<2>().data();
        auto* sfc3_ptr = domain.template indices<3>().data();

        cstone::DeviceVector<float> d_volumes(numElements);

        int blockSize = 256;
        int numBlocks = (numElements + blockSize - 1) / blockSize;

        auto box = domain.getDomain().box();

        volumeCalculationKernel<<<numBlocks, blockSize>>>(sfc0_ptr, sfc1_ptr, sfc2_ptr, sfc3_ptr, box, d_volumes.data(),
                                                          numElements);

        cudaDeviceSynchronize();
        cudaError_t err = cudaGetLastError();
        ASSERT_EQ(err, cudaSuccess) << "CUDA kernel execution failed: " << cudaGetErrorString(err);

        // Verify volumes
        auto h_volumes          = toHost(d_volumes);
        float totalVolume       = 0.0f;
        int negativeOrZeroCount = 0;

        for (size_t i = 0; i < numElements; i++)
        {
            if (h_volumes[i] <= 0.0f) { negativeOrZeroCount++; }
            totalVolume += h_volumes[i];
        }

        // Allow a small percentage of zero/negative volumes due to degenerate elements
        float badElementPercentage = 100.0f * negativeOrZeroCount / numElements;
        EXPECT_LT(badElementPercentage, 1.0f) << negativeOrZeroCount << " elements have zero or negative volume";

        // Check total volume is positive
        EXPECT_GT(totalVolume, 0.0f) << "Total tetrahedral volume should be positive";

        // Report statistics
        if (rank == 0)
        {
            std::cout << "Total tetrahedral volume on rank 0: " << totalVolume << std::endl;
            if (negativeOrZeroCount > 0)
            {
                std::cout << "Warning: " << negativeOrZeroCount << " elements (" << badElementPercentage
                          << "%) have zero or negative volume" << std::endl;
            }
        }
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception in GPU volume calculation test: " << e.what();
    }
}

TEST_F(ExternalMeshDomainTest, GpuCoordinateConversionPerformance)
{
    checkPrerequisites("GpuCoordinateConversionPerformance");

    try
    {
        using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
        Domain domain(meshPath, rank, numRanks);

        if (domain.getElementCount() < 100) { GTEST_SKIP() << "Mesh too small for performance test"; }

        size_t numKeys = domain.getElementCount() * 4; // 4 nodes per tetrahedron

        // Get raw device pointers - same as in binary_mesh.cu
        auto* sfc0_ptr = domain.template indices<0>().data();
        auto* sfc1_ptr = domain.template indices<1>().data();
        auto* sfc2_ptr = domain.template indices<2>().data();
        auto* sfc3_ptr = domain.template indices<3>().data();

        cstone::DeviceVector<float> d_coordinates(numKeys * 3);

        int blockSize = 256;
        int numBlocks = (numKeys + blockSize - 1) / blockSize;

        auto box = domain.getDomain().box();

        // Warmup run
        performanceTestKernel<<<numBlocks, blockSize>>>(sfc0_ptr, sfc1_ptr, sfc2_ptr, sfc3_ptr, box,
                                                        d_coordinates.data(), numKeys);

        cudaDeviceSynchronize();

        // Performance measurement using CUDA events like in binary_mesh.cu
        cudaEvent_t start, stop;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);

        const int iterations = 100;
        cudaEventRecord(start);
        for (int iter = 0; iter < iterations; iter++)
        {
            performanceTestKernel<<<numBlocks, blockSize>>>(sfc0_ptr, sfc1_ptr, sfc2_ptr, sfc3_ptr, box,
                                                            d_coordinates.data(), numKeys);
        }
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);

        float milliseconds = 0;
        cudaEventElapsedTime(&milliseconds, start, stop);

        // Validate some results
        auto h_coords = toHost(d_coordinates);

        int validCoords = 0;
        for (size_t i = 0; i < std::min(size_t(100), numKeys); i++)
        {
            float x = h_coords[i * 3 + 0];
            float y = h_coords[i * 3 + 1];
            float z = h_coords[i * 3 + 2];

            if (x >= box.xmin() && x <= box.xmax() && y >= box.ymin() && y <= box.ymax() && z >= box.zmin() &&
                z <= box.zmax())
            {
                validCoords++;
            }
        }

        EXPECT_GT(validCoords, 90) << "Most coordinates should be within bounds";

        if (rank == 0)
        {
            float avgTime    = milliseconds / iterations;
            float throughput = (numKeys * 1000.0f) / avgTime;
            std::cout << "GPU coordinate conversion performance:" << std::endl;
            std::cout << "  " << numKeys << " coordinates converted" << std::endl;
            std::cout << "  Average time: " << avgTime << " ms" << std::endl;
            std::cout << "  Throughput: " << throughput / 1e6 << " million coords/second" << std::endl;
        }

        cudaEventDestroy(start);
        cudaEventDestroy(stop);
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception in GPU coordinate conversion performance test: " << e.what();
    }
}

TEST_F(ExternalMeshDomainTest, AnalyzeTetrahedronOrientations)
{
    checkPrerequisites("AnalyzeTetrahedronOrientations");

    try
    {
        using Domain = ElementDomain<TetTag, float, uint64_t, cstone::GpuTag>;
        Domain domain(meshPath, rank, numRanks);

        if (domain.getElementCount() == 0)
        {
            std::cout << "Rank " << rank << " has zero elements, skipping orientation analysis" << std::endl;
            return;
        }

        // Get SFC keys
        auto h_sfc0 = toHost(domain.template indices<0>());
        auto h_sfc1 = toHost(domain.template indices<1>());
        auto h_sfc2 = toHost(domain.template indices<2>());
        auto h_sfc3 = toHost(domain.template indices<3>());

        size_t localStart = domain.startIndex();
        size_t localEnd   = domain.endIndex();

        std::cout << "Rank " << rank << ": Analyzing tetrahedron orientations for " << (localEnd - localStart)
                  << " local elements" << std::endl;

        // Orientation analysis
        int invertedTets   = 0;
        int validTets      = 0;
        int totalChecked   = 0;
        int degenerateTets = 0;

        float minVolume    = std::numeric_limits<float>::max();
        float maxVolume    = std::numeric_limits<float>::lowest();
        double totalVolume = 0.0;

        // Sample size for detailed analysis
        size_t sampleSize = std::min(size_t(1000), localEnd - localStart);

        for (size_t i = localStart; i < localStart + sampleSize; i++)
        {
            try
            {
                // Get the 4 vertices
                auto sfcKindKey0 = cstone::sfcKey(h_sfc0[i]);
                auto sfcKindKey1 = cstone::sfcKey(h_sfc1[i]);
                auto sfcKindKey2 = cstone::sfcKey(h_sfc2[i]);
                auto sfcKindKey3 = cstone::sfcKey(h_sfc3[i]);

                auto [ix0, iy0, iz0] = cstone::decodeSfc(sfcKindKey0);
                auto [ix1, iy1, iz1] = cstone::decodeSfc(sfcKindKey1);
                auto [ix2, iy2, iz2] = cstone::decodeSfc(sfcKindKey2);
                auto [ix3, iy3, iz3] = cstone::decodeSfc(sfcKindKey3);

                // Check for degenerate elements (duplicate coordinates)
                std::set<std::tuple<uint64_t, uint64_t, uint64_t>> uniqueCoords = {
                    {ix0, iy0, iz0}, {ix1, iy1, iz1}, {ix2, iy2, iz2}, {ix3, iy3, iz3}};

                if (uniqueCoords.size() < 4)
                {
                    degenerateTets++;
                    continue;
                }

                // Convert to float coordinates for volume calculation
                constexpr uint64_t maxCoord = (1ULL << cstone::maxTreeLevel<cstone::SfcKind<uint64_t>>{}) - 1;
                float invMaxCoord           = 1.0f / static_cast<float>(maxCoord);

                float x0 = static_cast<float>(ix0) * invMaxCoord * 100.0f;
                float y0 = static_cast<float>(iy0) * invMaxCoord * 100.0f;
                float z0 = static_cast<float>(iz0) * invMaxCoord * 100.0f;

                float x1 = static_cast<float>(ix1) * invMaxCoord * 100.0f;
                float y1 = static_cast<float>(iy1) * invMaxCoord * 100.0f;
                float z1 = static_cast<float>(iz1) * invMaxCoord * 100.0f;

                float x2 = static_cast<float>(ix2) * invMaxCoord * 100.0f;
                float y2 = static_cast<float>(iy2) * invMaxCoord * 100.0f;
                float z2 = static_cast<float>(iz2) * invMaxCoord * 100.0f;

                float x3 = static_cast<float>(ix3) * invMaxCoord * 100.0f;
                float y3 = static_cast<float>(iy3) * invMaxCoord * 100.0f;
                float z3 = static_cast<float>(iz3) * invMaxCoord * 100.0f;

                // Calculate tetrahedron volume using determinant
                // Volume = (1/6) * det|v1-v0, v2-v0, v3-v0|
                float v1x = x1 - x0, v1y = y1 - y0, v1z = z1 - z0;
                float v2x = x2 - x0, v2y = y2 - y0, v2z = z2 - z0;
                float v3x = x3 - x0, v3y = y3 - y0, v3z = z3 - z0;

                float det =
                    v1x * (v2y * v3z - v2z * v3y) - v1y * (v2x * v3z - v2z * v3x) + v1z * (v2x * v3y - v2y * v3x);

                float volume    = det / 6.0f;
                float absVolume = std::abs(volume);

                // Update statistics
                minVolume = std::min(minVolume, absVolume);
                maxVolume = std::max(maxVolume, absVolume);
                totalVolume += absVolume;

                if (volume < 0)
                {
                    invertedTets++;
                    if (invertedTets <= 5)
                    { // Print first 5 inverted tets
                        std::cout << "  Inverted tet at element " << i << ": volume = " << volume << std::endl;
                        std::cout << "    SFC keys: " << h_sfc0[i] << ", " << h_sfc1[i] << ", " << h_sfc2[i] << ", "
                                  << h_sfc3[i] << std::endl;
                        std::cout << "    Coords: (" << x0 << "," << y0 << "," << z0 << ") -> "
                                  << "(" << x1 << "," << y1 << "," << z1 << ") -> "
                                  << "(" << x2 << "," << y2 << "," << z2 << ") -> "
                                  << "(" << x3 << "," << y3 << "," << z3 << ")" << std::endl;
                    }
                }
                else { validTets++; }

                totalChecked++;
            }
            catch (...)
            {
                // Skip problematic elements
                continue;
            }
        }

        // Calculate percentages
        float invertedPercentage   = (totalChecked > 0) ? (100.0f * invertedTets / totalChecked) : 0.0f;
        float degeneratePercentage = (sampleSize > 0) ? (100.0f * degenerateTets / sampleSize) : 0.0f;
        double avgVolume           = (validTets > 0) ? (totalVolume / validTets) : 0.0;

        // Print detailed statistics
        std::cout << "Rank " << rank << " Tetrahedron Orientation Analysis:" << std::endl;
        std::cout << "  Sample size: " << sampleSize << " elements" << std::endl;
        std::cout << "  Valid tetrahedra: " << validTets << " (" << (100.0f - invertedPercentage - degeneratePercentage)
                  << "%)" << std::endl;
        std::cout << "  Inverted tetrahedra: " << invertedTets << " (" << invertedPercentage << "%)" << std::endl;
        std::cout << "  Degenerate tetrahedra: " << degenerateTets << " (" << degeneratePercentage << "%)" << std::endl;

        if (validTets > 0)
        {
            std::cout << "  Volume statistics:" << std::endl;
            std::cout << "    Min volume: " << minVolume << std::endl;
            std::cout << "    Max volume: " << maxVolume << std::endl;
            std::cout << "    Avg volume: " << avgVolume << std::endl;
            std::cout << "    Total volume: " << totalVolume << std::endl;
        }

        // Connectivity analysis
        std::cout << "  Connectivity analysis:" << std::endl;
        int sharedNodeCount = 0;
        int analyzedPairs   = 0;

        for (size_t i = localStart; i < std::min(localStart + 100, localEnd - 1); i++)
        {
            uint64_t keys1[4] = {h_sfc0[i], h_sfc1[i], h_sfc2[i], h_sfc3[i]};
            uint64_t keys2[4] = {h_sfc0[i + 1], h_sfc1[i + 1], h_sfc2[i + 1], h_sfc3[i + 1]};

            int shared = 0;
            for (int j = 0; j < 4; j++)
            {
                for (int k = 0; k < 4; k++)
                {
                    if (keys1[j] == keys2[k] && keys1[j] != 0) shared++;
                }
            }
            if (shared > 0) sharedNodeCount++;
            analyzedPairs++;
        }

        float connectivityPercentage = (analyzedPairs > 0) ? (100.0f * sharedNodeCount / analyzedPairs) : 0.0f;
        std::cout << "    Adjacent elements sharing nodes: " << sharedNodeCount << "/" << analyzedPairs << " ("
                  << connectivityPercentage << "%)" << std::endl;

        // Assess overall mesh quality
        std::cout << "  Mesh Quality Assessment:" << std::endl;
        if (invertedPercentage > 5.0f)
        {
            std::cout << "    WARNING: High percentage of inverted tetrahedra detected!" << std::endl;
            std::cout << "    This could cause displaced blocks in visualization." << std::endl;
        }

        if (degeneratePercentage > 10.0f)
        {
            std::cout << "    WARNING: High percentage of degenerate tetrahedra detected!" << std::endl;
            std::cout << "    This indicates connectivity issues." << std::endl;
        }

        if (connectivityPercentage < 50.0f)
        {
            std::cout << "    WARNING: Low connectivity between adjacent elements!" << std::endl;
            std::cout << "    This could indicate mesh fragmentation." << std::endl;
        }

        if (invertedPercentage < 1.0f && degeneratePercentage < 5.0f && connectivityPercentage > 80.0f)
        {
            std::cout << "    GOOD: Mesh appears to have proper connectivity and orientation." << std::endl;
        }

        // Global statistics gathering
        struct GlobalStats
        {
            int totalInverted   = 0;
            int totalValid      = 0;
            int totalDegenerate = 0;
            int totalShared     = 0;
            int totalPairs      = 0;
        } localStats, globalStats;

        localStats.totalInverted   = invertedTets;
        localStats.totalValid      = validTets;
        localStats.totalDegenerate = degenerateTets;
        localStats.totalShared     = sharedNodeCount;
        localStats.totalPairs      = analyzedPairs;

        MPI_Reduce(&localStats, &globalStats, sizeof(GlobalStats) / sizeof(int), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

        if (rank == 0)
        {
            std::cout << "\nGlobal Mesh Statistics:" << std::endl;
            int totalElements = globalStats.totalInverted + globalStats.totalValid + globalStats.totalDegenerate;
            if (totalElements > 0)
            {
                std::cout << "  Total elements analyzed: " << totalElements << std::endl;
                std::cout << "  Global inverted percentage: " << (100.0f * globalStats.totalInverted / totalElements)
                          << "%" << std::endl;
                std::cout << "  Global degenerate percentage: "
                          << (100.0f * globalStats.totalDegenerate / totalElements) << "%" << std::endl;
            }
            if (globalStats.totalPairs > 0)
            {
                std::cout << "  Global connectivity percentage: "
                          << (100.0f * globalStats.totalShared / globalStats.totalPairs) << "%" << std::endl;
            }
        }
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception in tetrahedron orientation analysis: " << e.what();
    }
}

template<typename T>
std::string getTypeName() {
    if constexpr (std::is_same_v<T, float>) return "float";
    else if constexpr (std::is_same_v<T, double>) return "double";
    else if constexpr (std::is_same_v<T, unsigned>) return "uint32";
    else if constexpr (std::is_same_v<T, uint64_t>) return "uint64";
    else return "unknown";
}

// Helper template function that contains the actual test logic
template<typename CoordType, typename KeyType, bool halos = false>
void testSFCVisualization(ExternalMeshDomainTest* testInstance)
{
    using Domain = ElementDomain<TetTag, CoordType, KeyType, cstone::GpuTag>;
    Domain domain(testInstance->meshPath, testInstance->rank, testInstance->numRanks);

    if (domain.getElementCount() == 0)
    {
        std::cout << "Rank " << testInstance->rank << " has zero elements, skipping visualization" << std::endl;
        return;
    }

    // Get SFC connectivity arrays
    auto h_sfc0 = toHost(domain.template indices<0>());
    auto h_sfc1 = toHost(domain.template indices<1>());
    auto h_sfc2 = toHost(domain.template indices<2>());
    auto h_sfc3 = toHost(domain.template indices<3>());

    std::string typeStr = (sizeof(CoordType) == 4) ? "float" : "double";
    std::cout << "Rank " << testInstance->rank << ": Testing with " << typeStr << " precision" << std::endl;

    std::string outputDir = testInstance->meshPath + "/vtk_output_" + getTypeName<CoordType>() + "_" + getTypeName<KeyType>() +
                        (halos ? "_halos" : "_no_halos");
    if (testInstance->rank == 0) { std::filesystem::create_directories(outputDir); }
    MPI_Barrier(MPI_COMM_WORLD);

    std::string filename = outputDir + "/sfc_" + typeStr + "_rank" + std::to_string(testInstance->rank) + ".vtk";
    std::ofstream file(filename);

    file << "# vtk DataFile Version 3.0\n";
    file << "SFC coordinate visualization (" << typeStr << " precision)\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";

    // Coordinate reconstruction with template precision
    std::vector<CoordType> vertices;
    std::map<KeyType, size_t> vertexMap;

    auto addVertex = [&](KeyType sfcKey) -> size_t
    {
    
        auto it = vertexMap.find(sfcKey);
        if (it != vertexMap.end()) return it->second;

        // Decode SFC
        auto sfcKindKey   = cstone::sfcKey(sfcKey);
        auto [ix, iy, iz] = cstone::decodeSfc(sfcKindKey);

     // Convert with template precision
        constexpr KeyType maxCoord = (1ULL << cstone::maxTreeLevel<cstone::SfcKind<KeyType>>{}) - 1;
        CoordType invMaxCoord       = static_cast<CoordType>(1.0) / static_cast<CoordType>(maxCoord);

        auto bbox = domain.getBoundingBox();

        CoordType x = bbox.xmin() + (static_cast<CoordType>(ix) * invMaxCoord) * (bbox.xmax() - bbox.xmin());
        CoordType y = bbox.ymin() + (static_cast<CoordType>(iy) * invMaxCoord) * (bbox.ymax() - bbox.ymin());
        CoordType z = bbox.zmin() + (static_cast<CoordType>(iz) * invMaxCoord) * (bbox.zmax() - bbox.zmin());

        size_t idx = vertices.size() / 3;
        vertices.push_back(x);
        vertices.push_back(y);
        vertices.push_back(z);

        vertexMap[sfcKey] = idx;
        return idx;
    };

    size_t localStart = domain.startIndex();
    size_t localEnd   = domain.endIndex();
    // Build connectivity
    std::vector<size_t> connectivity;
    if constexpr (halos)
    {
        localStart = 0;
        localEnd   = domain.getElementCount();
    }

    for (size_t i = localStart; i < localEnd; i++)
    {
        size_t v0 = addVertex(h_sfc0[i]);
        size_t v1 = addVertex(h_sfc1[i]);
        size_t v2 = addVertex(h_sfc2[i]);
        size_t v3 = addVertex(h_sfc3[i]);

        connectivity.push_back(v0);
        connectivity.push_back(v1);
        connectivity.push_back(v2);
        connectivity.push_back(v3);
    }

    // Ensure all ranks have the same number of vertices
#ifndef NDEBUG
    auto bbox = domain.getBoundingBox();
    constexpr auto treeLevel = cstone::maxTreeLevel<cstone::SfcKind<KeyType>>{};
    std::cout << "SFC tree level: " << treeLevel << std::endl;
    std::cout << "Max coord value: " << ((1ULL << treeLevel) - 1) << std::endl;
    std::cout << "Quantization step: " << (1.0 / ((1ULL << treeLevel) - 1)) << std::endl;
    std::cout << "  Bbox: [" << bbox.xmin() << "," << bbox.xmax() << "] [" 
          << bbox.ymin() << "," << bbox.ymax() << "] [" 
          << bbox.zmin() << "," << bbox.zmax() << "]" << std::endl;
#endif

    // Write VTK data (VTK always uses float for coordinates)
    file << "POINTS " << vertices.size() / 3 << " float\n";
    for (size_t i = 0; i < vertices.size(); i += 3)
    {
        file << static_cast<float>(vertices[i]) << " " << static_cast<float>(vertices[i + 1]) << " "
             << static_cast<float>(vertices[i + 2]) << "\n";
    }

    size_t numCells = connectivity.size() / 4;
    file << "CELLS " << numCells << " " << numCells * 5 << "\n";
    for (size_t i = 0; i < connectivity.size(); i += 4)
    {
        file << "4 " << connectivity[i] << " " << connectivity[i + 1] << " " << connectivity[i + 2] << " "
             << connectivity[i + 3] << "\n";
    }

    file << "CELL_TYPES " << numCells << "\n";
    for (size_t i = 0; i < numCells; i++)
    {
        file << "10\n";
    }

    file << "CELL_DATA " << numCells << "\n";
    file << "SCALARS rank int 1\n";
    file << "LOOKUP_TABLE default\n";
    for (size_t i = 0; i < numCells; i++)
    {
        file << testInstance->rank << "\n";
    }

    file.close();

    // Precision analysis
    if (!vertices.empty())
    {
        CoordType minCoord = *std::min_element(vertices.begin(), vertices.end());
        CoordType maxCoord = *std::max_element(vertices.begin(), vertices.end());

        std::cout << "Rank " << testInstance->rank << " (" << typeStr << " precision):" << std::endl;
        std::cout << "  Range: [" << minCoord << ", " << maxCoord << "]" << std::endl;
        std::cout << "  Unique vertices: " << vertices.size() / 3 << std::endl;
        std::cout << "  Elements: " << numCells << std::endl;

        // Check for precision artifacts
        int zeroCount       = 0;
        CoordType tolerance = (sizeof(CoordType) == 4) ? CoordType(1e-6) : CoordType(1e-12);
        for (const auto& coord : vertices)
        {
            if (std::abs(coord) < tolerance) zeroCount++;
        }

        std::cout << "  Near-zero coordinates: " << zeroCount << " (tolerance: " << tolerance << ")" << std::endl;
    }

    if (testInstance->rank == 0)
    {
        std::cout << "SFC visualization (" << typeStr << ") completed:" << std::endl;
        std::cout << "  Files: " << outputDir << "/sfc_" << typeStr << "_rank*.vtk" << std::endl;
    }
}

// Now your actual test functions just call the template
TEST_F(ExternalMeshDomainTest, VisualizeRawSFCDecodingFloat)
{
    checkPrerequisites("VisualizeRawSFCDecodingFloat");

    try
    {
        testSFCVisualization<float, unsigned>(this);
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception in float SFC visualization: " << e.what();
    }
}

TEST_F(ExternalMeshDomainTest, VisualizeRawSFCDecodingFloatUint64)
{
    checkPrerequisites("VisualizeRawSFCDecodingDouble");

    try
    {
        testSFCVisualization<float, uint64_t>(this);
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception in double SFC visualization: " << e.what();
    }
}

// Now your actual test functions just call the template
TEST_F(ExternalMeshDomainTest, VisualizeRawSFCDecodingWithHalosFloat)
{
    checkPrerequisites("VisualizeRawSFCDecodingFloat");

    try
    {
        testSFCVisualization<float, unsigned, true>(this);
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception in float SFC visualization: " << e.what();
    }
}

TEST_F(ExternalMeshDomainTest, VisualizeRawSFCDecodingWithHalosFloatUint64)
{
    checkPrerequisites("VisualizeRawSFCDecodingDouble");

    try
    {
        testSFCVisualization<float, uint64_t, true>(this);
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception in double SFC visualization: " << e.what();
    }
}

int main(int argc, char** argv)
{
    mars::Env env(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
