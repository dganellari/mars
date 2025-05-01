#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <algorithm>
#include <vector>
#include "domain.hpp"

namespace fs = std::filesystem;
using namespace mars;

class ExternalMeshDomainTest : public ::testing::Test
{
protected:
    int rank;
    int numRanks;

    void SetUp() override
    {
        // Get MPI rank and size
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
        printf("Rank %d of %d\n", rank, numRanks);
    }

    // Get mesh paths from environment variable
    std::string getMeshPath() const
    {
        const char* meshPathEnv = std::getenv("MESH_PATH");
        std::string meshPath    = meshPathEnv ? meshPathEnv : "";

        if (meshPath.empty() || !fs::exists(meshPath))
        {
            // Fallback to common test locations
            std::vector<std::string> commonLocations = {"./test_data", "../test_data", "./meshes", "../meshes"};

            for (const auto& path : commonLocations)
            {
                if (fs::exists(path) && fs::is_directory(fs::path(path)))
                {
                    // Check if directory contains required mesh files
                    if (fs::exists(fs::path(path) / "x.float32") || fs::exists(fs::path(path) / "x.double"))
                    {
                        return path;
                    }
                }
            }
        }

        return meshPath;
    }
};

// Basic test for external mesh domain creation
TEST_F(ExternalMeshDomainTest, BasicDomainCreation)
{
    // Get mesh path
    std::string meshPath = getMeshPath();
    if (meshPath.empty() || !fs::exists(meshPath))
    {
        GTEST_SKIP() << "No valid mesh directory found";
        return;
    }

    // Check if directory contains required files
    if (!fs::exists(fs::path(meshPath) / "x.float32") && !fs::exists(fs::path(meshPath) / "x.double"))
    {
        GTEST_SKIP() << "Mesh directory does not contain required coordinate files";
        return;
    }

    std::cout << "Creating domain from mesh at: " << meshPath << std::endl;

    // Create a catch block to handle exceptions
    try
    {
        // Set device - critical for multi-rank runs
        int deviceCount;
        cudaGetDeviceCount(&deviceCount);
        if (deviceCount == 0)
        {
            GTEST_SKIP() << "No CUDA devices available, skipping test";
            return;
        }

        // Create domain with tetrahedral elements
        using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
        Domain domain(meshPath, rank, numRanks);

        // Basic validation tests
        EXPECT_GT(domain.getNodeCount(), 0) << "Domain should have nodes";
        EXPECT_GT(domain.getElementCount(), 0) << "Domain should have elements";

        // Print basic information
        std::cout << "Rank " << rank << " domain: " << domain.getNodeCount() << " nodes, " << domain.getElementCount()
                  << " elements" << std::endl;

        // Gather total counts across ranks
        size_t localElementCount = domain.getElementCount();
        size_t localNodeCount    = domain.getNodeCount();
        size_t totalElementCount = 0;
        size_t totalNodeCount    = 0;

        MPI_Allreduce(&localElementCount, &totalElementCount, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&localNodeCount, &totalNodeCount, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

        if (rank == 0)
        {
            std::cout << "Total across all ranks: " << totalNodeCount << " nodes, " << totalElementCount << " elements"
                      << std::endl;
        }

        // Verify coordinate data
        auto h_x = toHost(domain.x());
        auto h_y = toHost(domain.y());
        auto h_z = toHost(domain.z());

        EXPECT_EQ(h_x.size(), domain.getNodeCount());
        EXPECT_EQ(h_y.size(), domain.getNodeCount());
        EXPECT_EQ(h_z.size(), domain.getNodeCount());

        // Check connectivity indices are valid
        if (domain.getElementCount() > 0)
        {
            auto h_i0 = toHost(domain.indices<0>());
            auto h_i1 = toHost(domain.indices<1>());
            auto h_i2 = toHost(domain.indices<2>());
            auto h_i3 = toHost(domain.indices<3>());

            EXPECT_EQ(h_i0.size(), domain.getElementCount());

            // Check first element's node indices are valid
            for (size_t i = 0; i < std::min(domain.getElementCount(), size_t(5)); i++)
            {
                EXPECT_LT(h_i0[i], domain.getNodeCount()) << "Invalid node index in element " << i;
                EXPECT_LT(h_i1[i], domain.getNodeCount()) << "Invalid node index in element " << i;
                EXPECT_LT(h_i2[i], domain.getNodeCount()) << "Invalid node index in element " << i;
                EXPECT_LT(h_i3[i], domain.getNodeCount()) << "Invalid node index in element " << i;
            }
        }

        // Check domain bounding box
        auto box = domain.getDomain().box();
        std::cout << "Domain box: [" << box.xmin() << ", " << box.xmax() << "] x [" << box.ymin() << ", " << box.ymax()
                  << "] x [" << box.zmin() << ", " << box.zmax() << "]" << std::endl;

        EXPECT_LT(box.xmin(), box.xmax()) << "Invalid X dimension in bounding box";
        EXPECT_LT(box.ymin(), box.ymax()) << "Invalid Y dimension in bounding box";
        EXPECT_LT(box.zmin(), box.zmax()) << "Invalid Z dimension in bounding box";
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception creating domain: " << e.what();
    }
}

// Test with different precisions - float and double
TEST_F(ExternalMeshDomainTest, DifferentPrecisions)
{
    // Get mesh path
    std::string meshPath = getMeshPath();
    if (meshPath.empty() || !fs::exists(meshPath))
    {
        GTEST_SKIP() << "No valid mesh directory found";
        return;
    }

    // Check for device availability
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount == 0)
    {
        GTEST_SKIP() << "No CUDA devices available, skipping test";
        return;
    }
    cudaSetDevice(rank % deviceCount);

    try
    {
        // Test with float precision
        {
            using FloatDomain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
            FloatDomain domain(meshPath, rank, numRanks);

            std::cout << "Float domain on rank " << rank << ": " << domain.getNodeCount() << " nodes, "
                      << domain.getElementCount() << " elements" << std::endl;

            EXPECT_GT(domain.getNodeCount(), 0);
        }

        // Test with double precision if x.double exists
        if (fs::exists(fs::path(meshPath) / "x.double"))
        {
            using DoubleDomain = ElementDomain<TetTag, double, uint64_t, cstone::GpuTag>;
            DoubleDomain domain(meshPath, rank, numRanks);

            std::cout << "Double domain on rank " << rank << ": " << domain.getNodeCount() << " nodes, "
                      << domain.getElementCount() << " elements" << std::endl;

            EXPECT_GT(domain.getNodeCount(), 0);
        }
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception in precision test: " << e.what();
    }
}

// Main function
int main(int argc, char** argv)
{
    // Initialize MPI through Mars environment
    mars::Env env(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
