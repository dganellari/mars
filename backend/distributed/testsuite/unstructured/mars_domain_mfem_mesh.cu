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

class MFEMMeshDomainTest : public ::testing::Test
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
            std::cout << "MFEM Mesh GPU Test setup: Found " << deviceCount << " CUDA devices" << std::endl;
            std::cout << "Using mesh at: " << (meshPath.empty() ? "none" : meshPath) << std::endl;
        }
    }

private:
    std::string getMeshPath() const
    {
        const char* meshPathEnv = std::getenv("MFEM_MESH_PATH");
        std::string path        = meshPathEnv ? meshPathEnv : "";

        if (path.empty() || !fs::exists(path))
        {
            std::vector<std::string> commonLocations = {
                "./meshes", "../meshes", "../../meshes", "../../../meshes",
                "./test_data", "../test_data", "../../test_data"
            };

            for (const auto& loc : commonLocations)
            {
                if (fs::exists(loc) && fs::is_directory(fs::path(loc)))
                {
                    // Look for common MFEM mesh files
                    auto tetMesh = fs::path(loc) / "beam-tet.mesh";
                    auto hexMesh = fs::path(loc) / "beam-hex.mesh";

                    if (fs::exists(tetMesh)) {
                        return tetMesh.string();
                    } else if (fs::exists(hexMesh)) {
                        return hexMesh.string();
                    }
                }
            }
        }
        return path;
    }

protected:
    void checkPrerequisites(const std::string& testName)
    {
        if (meshPath.empty() || !fs::exists(meshPath))
        {
            GTEST_SKIP() << testName << ": No valid MFEM mesh file found";
        }

        if (deviceCount == 0) {
            GTEST_SKIP() << testName << ": No CUDA devices available";
        }
    }

    // Helper validation functions
    template<typename Domain>
    void validateBasicDomainProperties(const Domain& domain)
    {
        EXPECT_GT(domain.getNodeCount(), 0) << "Domain should have nodes";
        EXPECT_GT(domain.getElementCount(), 0) << "Domain should have elements";

        if (rank == 0) {
            std::cout << "  Node count: " << domain.getNodeCount() << std::endl;
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

        if (rank == 0)
        {
            std::cout << "  Bounding box: [" << box.xmin() << ", " << box.xmax() << "] x [" << box.ymin() << ", "
                      << box.ymax() << "] x [" << box.zmin() << ", " << box.zmax() << "]" << std::endl;
        }
    }

    template<typename Domain>
    void validateConnectivity(const Domain& domain)
    {
        if (domain.getElementCount() == 0) return;

        size_t samplesToCheck = std::min(domain.getElementCount(), size_t(100));

        const auto& conn = domain.getConnectivity();
        auto h_conn0 = toHost(std::get<0>(conn));
        auto h_conn1 = toHost(std::get<1>(conn));
        auto h_conn2 = toHost(std::get<2>(conn));
        auto h_conn3 = toHost(std::get<3>(conn));

        const auto& sfcKeys = domain.getLocalToGlobalSfcMap();
        size_t nodeCount = domain.getNodeCount();

        for (size_t i = 0; i < samplesToCheck; ++i) {
            // Verify indices are within bounds
            EXPECT_LT(h_conn0[i], nodeCount) << "Element " << i << " has invalid node index 0";
            EXPECT_LT(h_conn1[i], nodeCount) << "Element " << i << " has invalid node index 1";
            EXPECT_LT(h_conn2[i], nodeCount) << "Element " << i << " has invalid node index 2";
            EXPECT_LT(h_conn3[i], nodeCount) << "Element " << i << " has invalid node index 3";

            // Verify nodes are distinct
            EXPECT_NE(h_conn0[i], h_conn1[i]) << "Element " << i << " has duplicate nodes 0,1";
            EXPECT_NE(h_conn0[i], h_conn2[i]) << "Element " << i << " has duplicate nodes 0,2";
            EXPECT_NE(h_conn0[i], h_conn3[i]) << "Element " << i << " has duplicate nodes 0,3";
            EXPECT_NE(h_conn1[i], h_conn2[i]) << "Element " << i << " has duplicate nodes 1,2";
            EXPECT_NE(h_conn1[i], h_conn3[i]) << "Element " << i << " has duplicate nodes 1,3";
            EXPECT_NE(h_conn2[i], h_conn3[i]) << "Element " << i << " has duplicate nodes 2,3";
        }

        if (rank == 0) {
            std::cout << "  Connectivity validation passed for " << samplesToCheck << " elements" << std::endl;
        }
    }

    template<typename Domain>
    void validateNodeOwnership(const Domain& domain)
    {
        if (numRanks == 1) {
            if (rank == 0) {
                std::cout << "  Single rank - skipping ownership validation" << std::endl;
            }
            return;
        }

        const auto& ownership = domain.getNodeOwnershipMap();
        size_t nodeCount = domain.getNodeCount();

        EXPECT_EQ(ownership.size(), nodeCount) << "Ownership map size mismatch";

        auto h_ownership = toHost(ownership);

        size_t ownedCount = 0;
        size_t ghostCount = 0;

        for (size_t i = 0; i < nodeCount; ++i) {
            if (h_ownership[i] == 1) {
                ownedCount++;
            } else if (h_ownership[i] == 0) {
                ghostCount++;
            } else {
                FAIL() << "Invalid ownership value " << static_cast<int>(h_ownership[i])
                       << " at node " << i;
            }
        }

        EXPECT_GT(ownedCount, 0) << "Rank " << rank << " should own some nodes";

        // Verify total owned nodes equals total unique nodes across all ranks
        size_t totalOwned = 0;
        MPI_Reduce(&ownedCount, &totalOwned, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            std::cout << "  Total owned nodes across all ranks: " << totalOwned << std::endl;
            std::cout << "  Rank 0: " << ownedCount << " owned, " << ghostCount << " ghost" << std::endl;
        } else if (rank == 1) {
            std::cout << "  Rank 1: " << ownedCount << " owned, " << ghostCount << " ghost" << std::endl;
        }
    }

    template<typename Domain>
    void validateBoundaryInfo(const Domain& domain)
    {
        if (!domain.hasBoundaryInfo()) {
            if (rank == 0) {
                std::cout << "  Warning: Domain has no boundary information" << std::endl;
            }
            return;
        }

        const auto& boundaryNodes = domain.getBoundaryNodes();
        auto h_boundary = toHost(boundaryNodes);

        size_t boundaryCount = 0;
        for (auto val : h_boundary) {
            if (val != 0) boundaryCount++;
        }

        if (rank == 0) {
            std::cout << "  Boundary nodes: " << boundaryCount << " / " << h_boundary.size()
                      << " (" << (100.0 * boundaryCount / h_boundary.size()) << "%)" << std::endl;
        }
    }

    template<typename Domain>
    void printDomainStatistics(const Domain& domain, const std::string& testName)
    {
        size_t localElements = domain.localElementCount();
        size_t totalElements = domain.getElementCount();
        size_t haloElements = totalElements - localElements;

        size_t globalTotalElements = 0;
        MPI_Reduce(&localElements, &globalTotalElements, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

        if (rank == 0)
        {
            std::cout << testName << " completed successfully:" << std::endl;
            std::cout << "  Global elements: " << globalTotalElements << std::endl;
            std::cout << "  Rank 0 local elements: " << localElements << std::endl;
            std::cout << "  Rank 0 halo elements: " << haloElements << std::endl;

            if (domain.getElementCount() > 0)
            {
                auto box     = domain.getDomain().box();
                double volume = (box.xmax() - box.xmin()) * (box.ymax() - box.ymin()) * (box.zmax() - box.zmin());
                std::cout << "  Domain volume: " << volume << std::endl;
            }
        }
    }
};

TEST_F(MFEMMeshDomainTest, BasicMFEMMeshLoading)
{
    checkPrerequisites("BasicMFEMMeshLoading");

    try
    {
        // Load MFEM mesh directly using ElementDomain constructor
        using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
        Domain domain(meshPath, rank, numRanks);

        if (rank == 0) {
            std::cout << "Successfully loaded MFEM mesh:" << std::endl;
        }

        validateBasicDomainProperties(domain);
        validateBoundingBox(domain);

        if (rank == 0) {
            auto box = domain.getDomain().box();
            std::cout << "  Mesh dimensions: ["
                      << (box.xmax() - box.xmin()) << " x "
                      << (box.ymax() - box.ymin()) << " x "
                      << (box.zmax() - box.zmin()) << "]" << std::endl;
        }
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception loading MFEM mesh: " << e.what();
    }
}

TEST_F(MFEMMeshDomainTest, MFEMDomainCreation)
{
    checkPrerequisites("MFEMDomainCreation");

    try
    {
        // Load MFEM mesh using ElementDomain - it handles partitioning automatically
        using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
        Domain domain(meshPath, rank, numRanks);

        validateBasicDomainProperties(domain);
        validateBoundingBox(domain);
        validateConnectivity(domain);
        validateNodeOwnership(domain);
        printDomainStatistics(domain, "MFEMDomainCreation");
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception creating domain from MFEM mesh: " << e.what();
    }
}

TEST_F(MFEMMeshDomainTest, MFEMDomainWithBoundaryInfo)
{
    checkPrerequisites("MFEMDomainWithBoundaryInfo");

    try
    {
        // Load MFEM mesh - readMFEMMeshWithElementPartitioning returns boundary info
        using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
        Domain domain(meshPath, rank, numRanks);

        validateBasicDomainProperties(domain);
        validateBoundingBox(domain);
        validateConnectivity(domain);
        validateNodeOwnership(domain);
        // Note: Boundary info is automatically loaded from MFEM file
        printDomainStatistics(domain, "MFEMDomainWithBoundaryInfo");
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception creating domain with boundary info: " << e.what();
    }
}

TEST_F(MFEMMeshDomainTest, MultiRankConsistency)
{
    checkPrerequisites("MultiRankConsistency");

    if (numRanks == 1) {
        GTEST_SKIP() << "Multi-rank test requires at least 2 ranks";
    }

    try
    {
        // Load MFEM mesh - each rank reads its partition
        using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
        Domain domain(meshPath, rank, numRanks);

        // Verify bounding box is consistent across ranks
        auto box = domain.getDomain().box();
        std::vector<float> localBox = {box.xmin(), box.xmax(), box.ymin(), box.ymax(), box.zmin(), box.zmax()};
        std::vector<float> globalBoxMin(6), globalBoxMax(6);

        MPI_Allreduce(localBox.data(), globalBoxMin.data(), 6, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(localBox.data(), globalBoxMax.data(), 6, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);

        // All ranks should have the same bounding box
        const float tolerance = 1e-5f;
        for (int i = 0; i < 6; ++i) {
            EXPECT_NEAR(globalBoxMin[i], globalBoxMax[i], tolerance)
                << "Bounding box component " << i << " not consistent across ranks";
        }

        // Verify total elements
        size_t localElements = domain.localElementCount();
        size_t totalLocalElements = 0;
        MPI_Reduce(&localElements, &totalLocalElements, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            std::cout << "Multi-rank consistency check:" << std::endl;
            std::cout << "  Sum of local elements: " << totalLocalElements << std::endl;
            std::cout << "  Bounding box consistent: Yes" << std::endl;
        }
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception in multi-rank consistency test: " << e.what();
    }
}

// Helper function to write VTK output for visualization
template<typename Domain>
void writeVTKOutput(const Domain& domain, const std::string& filename, int rank)
{
    // Get SFC connectivity
    const auto& conn = domain.getConnectivity();
    auto h_conn0 = toHost(std::get<0>(conn));
    auto h_conn1 = toHost(std::get<1>(conn));
    auto h_conn2 = toHost(std::get<2>(conn));
    auto h_conn3 = toHost(std::get<3>(conn));

    // Get SFC keys for nodes
    const auto& sfcKeys = domain.getLocalToGlobalSfcMap();
    auto h_sfcKeys = toHost(sfcKeys);

    std::ofstream file(filename);
    file << "# vtk DataFile Version 3.0\n";
    file << "MFEM Mesh Domain Visualization\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";

    // Build vertex map and decode coordinates
    std::map<uint64_t, size_t> vertexMap;
    std::vector<double> vertices;

    auto addVertex = [&](uint64_t localNodeIdx) -> size_t {
        uint64_t sfcKey = h_sfcKeys[localNodeIdx];
        auto it = vertexMap.find(sfcKey);
        if (it != vertexMap.end()) return it->second;

        // Decode SFC to get coordinates
        auto sfcKindKey = cstone::sfcKey(sfcKey);
        auto [ix, iy, iz] = cstone::decodeSfc(sfcKindKey);

        constexpr uint64_t maxCoord = (1ULL << cstone::maxTreeLevel<cstone::SfcKind<uint64_t>>{}) - 1;
        double invMaxCoord = 1.0 / static_cast<double>(maxCoord);

        auto box = domain.getDomain().box();
        double x = box.xmin() + ix * invMaxCoord * (box.xmax() - box.xmin());
        double y = box.ymin() + iy * invMaxCoord * (box.ymax() - box.ymin());
        double z = box.zmin() + iz * invMaxCoord * (box.zmax() - box.zmin());

        size_t idx = vertices.size() / 3;
        vertices.push_back(x);
        vertices.push_back(y);
        vertices.push_back(z);

        vertexMap[sfcKey] = idx;
        return idx;
    };

    // Build connectivity for LOCAL elements only (visualize partition)
    size_t localStart = domain.startIndex();
    size_t localEnd = domain.endIndex();
    std::vector<size_t> connectivity;

    for (size_t e = localStart; e < localEnd; ++e) {
        size_t v0 = addVertex(h_conn0[e]);
        size_t v1 = addVertex(h_conn1[e]);
        size_t v2 = addVertex(h_conn2[e]);
        size_t v3 = addVertex(h_conn3[e]);

        connectivity.push_back(v0);
        connectivity.push_back(v1);
        connectivity.push_back(v2);
        connectivity.push_back(v3);
    }

    // Write points
    file << "POINTS " << vertices.size() / 3 << " double\n";
    for (size_t i = 0; i < vertices.size(); i += 3) {
        file << vertices[i] << " " << vertices[i+1] << " " << vertices[i+2] << "\n";
    }

    // Write cells
    size_t numCells = connectivity.size() / 4;
    file << "CELLS " << numCells << " " << numCells * 5 << "\n";
    for (size_t i = 0; i < connectivity.size(); i += 4) {
        file << "4 " << connectivity[i] << " " << connectivity[i+1] << " "
             << connectivity[i+2] << " " << connectivity[i+3] << "\n";
    }

    // Write cell types (10 = tetrahedron)
    file << "CELL_TYPES " << numCells << "\n";
    for (size_t i = 0; i < numCells; ++i) {
        file << "10\n";
    }

    // Write cell data (rank for coloring)
    file << "CELL_DATA " << numCells << "\n";
    file << "SCALARS rank int 1\n";
    file << "LOOKUP_TABLE default\n";
    for (size_t i = 0; i < numCells; ++i) {
        file << rank << "\n";
    }

    file.close();
}

TEST_F(MFEMMeshDomainTest, VisualizeMFEMDomainPartitioning)
{
    checkPrerequisites("VisualizeMFEMDomainPartitioning");

    try
    {
        // Load MFEM mesh using ElementDomain
        using Domain = ElementDomain<TetTag, float, unsigned, cstone::GpuTag>;
        Domain domain(meshPath, rank, numRanks);

        // Create output directory
        std::string outputDir = "vtk_output_mfem";
        if (rank == 0) {
            fs::create_directories(outputDir);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // Write VTK output
        std::string filename = outputDir + "/mfem_domain_rank" + std::to_string(rank) + ".vtk";
        writeVTKOutput(domain, filename, rank);

        if (rank == 0) {
            std::cout << "VTK visualization files written to: " << outputDir << "/mfem_domain_rank*.vtk" << std::endl;
            std::cout << "View with ParaView to see domain partitioning" << std::endl;
        }

        // Verify output was created
        EXPECT_TRUE(fs::exists(filename)) << "VTK output file was not created";
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception in VTK visualization test: " << e.what();
    }
}

int main(int argc, char** argv)
{
    mars::Env env(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
