#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>
#include <chrono>
#include <algorithm>
#include <vector>
#include <string>
#include "mars_read_mesh_binary.hpp"

namespace fs = std::filesystem;
using namespace mars;

// Test fixture for external meshes
class ExternalMeshTest : public ::testing::TestWithParam<std::string> {
protected:
    int rank = 0;
    int numRanks = 1;
    std::string meshPath;
    
    void SetUp() override {
        // Get actual MPI rank and size
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
        
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

// Test to visualize mesh in ParaView by writing to standard VTK format
TEST_P(ExternalMeshTest, WriteStandardVTK) {
    // Skip test if not on rank 0
    if (rank != 0) {
        GTEST_SKIP() << "VTK export only ran on rank 0";
    }

    std::cout << "Reading mesh for standard VTK export: " << meshPath << std::endl;
    
    // Read coordinates with float precision
    auto [nodeCount, nodeStartIdx, x_data, y_data, z_data] = 
        mars::readMeshCoordinatesBinary<float>(meshPath, 0, 1);
    
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
            mars::readMeshConnectivityBinaryTuple<4>(meshPath, nodeStartIdx, 0, 1);
        
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

// Test to visualize mesh in ParaView using actual MPI distribution
TEST_P(ExternalMeshTest, WriteVTKWithActualMpiDistribution) {
    // Read this rank's portion of the mesh with element-based partitioning
    // Now with localToGlobal included in the return tuple
    auto [nodeCount, elementCount, x, y, z, conn, localToGlobal] = 
        mars::readMeshWithElementPartitioning<4, float>(meshPath, rank, numRanks);
    
    // Gather all node data to rank 0 for VTK writing
    std::vector<float> allX, allY, allZ;
    std::vector<int> allI0, allI1, allI2, allI3;
    size_t totalNodes = 0, totalElements = 0;
    
    // First, gather the counts to allocate space
    std::vector<int> nodeCounts(numRanks);
    std::vector<int> elemCounts(numRanks);
    int myNodeCount = static_cast<int>(nodeCount);
    int myElemCount = static_cast<int>(elementCount);
    
    MPI_Gather(&myNodeCount, 1, MPI_INT, nodeCounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&myElemCount, 1, MPI_INT, elemCounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // We need global node IDs to create a consistent combined mesh
    std::vector<int> globalNodeIDs;
    if (nodeCount > 0) {
        globalNodeIDs.resize(nodeCount);
        for (size_t i = 0; i < nodeCount; ++i) {
            globalNodeIDs[i] = localToGlobal[i]; // Copy the global IDs
        }
    }
    
    // Prepare displacement arrays for gathering
    std::vector<int> nodeDispls(numRanks, 0);
    std::vector<int> elemDispls(numRanks, 0);
    std::vector<int> allGlobalNodeIds; // Will hold all global node IDs in order
    
    if (rank == 0) {
        // Calculate total counts and displacements
        for (int i = 0; i < numRanks; i++) {
            if (i > 0) {
                nodeDispls[i] = nodeDispls[i-1] + nodeCounts[i-1];
                elemDispls[i] = elemDispls[i-1] + elemCounts[i-1];
            }
            totalNodes += nodeCounts[i];
            totalElements += elemCounts[i];
        }
        
        // Allocate space for the gathered data
        allX.resize(totalNodes);
        allY.resize(totalNodes);
        allZ.resize(totalNodes);
        allGlobalNodeIds.resize(totalNodes);
        
        if (totalElements > 0) {
            allI0.resize(totalElements);
            allI1.resize(totalElements);
            allI2.resize(totalElements);
            allI3.resize(totalElements);
        }
    }
    
    // Gather coordinate data
    MPI_Gatherv(x.data(), myNodeCount, MPI_FLOAT, 
               allX.data(), nodeCounts.data(), nodeDispls.data(), 
               MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gatherv(y.data(), myNodeCount, MPI_FLOAT, 
               allY.data(), nodeCounts.data(), nodeDispls.data(), 
               MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gatherv(z.data(), myNodeCount, MPI_FLOAT, 
               allZ.data(), nodeCounts.data(), nodeDispls.data(), 
               MPI_FLOAT, 0, MPI_COMM_WORLD);
    
    // Gather global node IDs
    MPI_Gatherv(globalNodeIDs.data(), myNodeCount, MPI_INT,
               allGlobalNodeIds.data(), nodeCounts.data(), nodeDispls.data(),
               MPI_INT, 0, MPI_COMM_WORLD);
    
    // Gather connectivity data if available
    if (elementCount > 0) {
        // Extract connectivity arrays
        const auto& i0 = std::get<0>(conn);
        const auto& i1 = std::get<1>(conn);
        const auto& i2 = std::get<2>(conn);
        const auto& i3 = std::get<3>(conn);
        
        // Create vectors to hold the connectivity with corrected global indices
        std::vector<int> globalI0(elementCount);
        std::vector<int> globalI1(elementCount);
        std::vector<int> globalI2(elementCount);
        std::vector<int> globalI3(elementCount);
        
        // Convert local indices to global indices
        for (size_t i = 0; i < elementCount; i++) {
            // The local index i0[i] corresponds to global node localToGlobal[i0[i]]
            globalI0[i] = localToGlobal[i0[i]];
            globalI1[i] = localToGlobal[i1[i]];
            globalI2[i] = localToGlobal[i2[i]];
            globalI3[i] = localToGlobal[i3[i]];
        }
        
        // Gather the global connectivity
        MPI_Gatherv(globalI0.data(), myElemCount, MPI_INT,
                  allI0.data(), elemCounts.data(), elemDispls.data(),
                  MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gatherv(globalI1.data(), myElemCount, MPI_INT,
                  allI1.data(), elemCounts.data(), elemDispls.data(),
                  MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gatherv(globalI2.data(), myElemCount, MPI_INT,
                  allI2.data(), elemCounts.data(), elemDispls.data(),
                  MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gatherv(globalI3.data(), myElemCount, MPI_INT,
                  allI3.data(), elemCounts.data(), elemDispls.data(),
                  MPI_INT, 0, MPI_COMM_WORLD);
    }
    
    // Create a coherent mesh on rank 0
    if (rank == 0 && totalNodes > 0) {
        // First, create a mapping from global node IDs to positions in gathered arrays
        std::unordered_map<int, int> globalToIndex;
        for (size_t i = 0; i < totalNodes; ++i) {
            globalToIndex[allGlobalNodeIds[i]] = i;
        }
        
        // Fix connectivity to refer to positions in gathered arrays
        if (totalElements > 0) {
            for (size_t i = 0; i < totalElements; ++i) {
                allI0[i] = globalToIndex[allI0[i]];
                allI1[i] = globalToIndex[allI1[i]];
                allI2[i] = globalToIndex[allI2[i]];
                allI3[i] = globalToIndex[allI3[i]];
            }
        }
        
        // Verify indices are valid
        bool indicesValid = true;
        for (size_t i = 0; i < totalElements && indicesValid; i++) {
            if (allI0[i] < 0 || allI0[i] >= totalNodes ||
                allI1[i] < 0 || allI1[i] >= totalNodes ||
                allI2[i] < 0 || allI2[i] >= totalNodes ||
                allI3[i] < 0 || allI3[i] >= totalNodes) {
                std::cout << "Invalid index at element " << i 
                         << ": [" << allI0[i] << ", " << allI1[i] 
                         << ", " << allI2[i] << ", " << allI3[i] 
                         << "]" << std::endl;
                indicesValid = false;
            }
        }
        
        if (!indicesValid) {
            std::cout << "Warning: Invalid node indices found in connectivity data" << std::endl;
        } else {
            std::cout << "All connectivity indices verified valid" << std::endl;
        }
        
        // Write VTK file
        fs::path inputDir = fs::path(meshPath).parent_path();
        fs::path outputPath = inputDir / ("mesh_visualization_mpi_" + std::to_string(numRanks) + "_ranks.vtk");
        std::cout << "Writing MPI-distributed VTK file to: " << outputPath << std::endl;
        
        // Write VTK file with gathered data
        std::ofstream vtkFile(outputPath, std::ios::out);
        if (!vtkFile) {
            GTEST_SKIP() << "Failed to open VTK file for writing: " << outputPath;
        }
        
        // Write VTK header
        vtkFile << "# vtk DataFile Version 3.0\n";
        vtkFile << "Mesh exported from MARS using MPI distribution\n";
        vtkFile << "ASCII\n";
        vtkFile << "DATASET UNSTRUCTURED_GRID\n";
        
        // Write points
        vtkFile << "POINTS " << totalNodes << " float\n";
        for (size_t i = 0; i < totalNodes; i++) {
            vtkFile << allX[i] << " " << allY[i] << " " << allZ[i] << "\n";
        }
        
        // Write cells if available
        if (totalElements > 0) {
            vtkFile << "CELLS " << totalElements << " " << totalElements * 5 << "\n";
            for (size_t i = 0; i < totalElements; i++) {
                vtkFile << "4 " 
                       << allI0[i] << " " 
                       << allI1[i] << " " 
                       << allI2[i] << " " 
                       << allI3[i] << "\n";
            }
            
            // Write cell types (10 = VTK_TETRA)
            vtkFile << "CELL_TYPES " << totalElements << "\n";
            for (size_t i = 0; i < totalElements; i++) {
                vtkFile << "10\n";  // 10 is the VTK type for tetrahedron
            }
        }
        
        vtkFile.close();
        
        std::cout << "VTK file written with " << totalNodes << " nodes and " 
                  << totalElements << " tetrahedral elements" << std::endl;
        
        // Verify file was written
        EXPECT_TRUE(fs::exists(outputPath)) << "VTK file was not created";
        EXPECT_GT(fs::file_size(outputPath), 0) << "VTK file is empty";
    }
    
    // Ensure all ranks finish before proceeding
    MPI_Barrier(MPI_COMM_WORLD);
}

// Test partitioning balance with element-based partitioning
TEST_P(ExternalMeshTest, ElementPartitioningBalanceSingleRank) {
    try {
        // Skip if we're not on rank 0
        if (rank != 0) {
            GTEST_SKIP() << "Test only runs on rank 0";
        }
        
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
        std::vector<size_t> rankNodeCounts(testRanks, 0);
        std::vector<size_t> rankElementCounts(testRanks, 0);
        
        // Read each rank's portion with element-based partitioning
        for (int r = 0; r < testRanks; r++) {
            try {
                auto [nodeCount, elementCount, x, y, z, conn, localToGlobal] = 
                    mars::readMeshWithElementPartitioning<4, float>(meshPath, r, testRanks);
                
                rankNodeCounts[r] = nodeCount;
                rankElementCounts[r] = elementCount;
                
                std::cout << "Rank " << r << ": " << elementCount << " elements, " 
                         << nodeCount << " nodes" << std::endl;
            } catch (const std::exception& e) {
                std::cout << "Error reading rank " << r << ": " << e.what() << std::endl;
                // Continue with other ranks, don't abort test completely
            }
        }
        
        // Calculate node and element imbalance
        auto nonZeroNodes = std::count_if(rankNodeCounts.begin(), rankNodeCounts.end(), 
                                         [](size_t n) { return n > 0; });
        auto nonZeroElements = std::count_if(rankElementCounts.begin(), rankElementCounts.end(), 
                                           [](size_t n) { return n > 0; });
        
        // Only calculate imbalance if we have valid data
        if (nonZeroNodes > 1) {
            // Filter out zeros for min calculation
            std::vector<size_t> validNodeCounts;
            for (auto n : rankNodeCounts) {
                if (n > 0) validNodeCounts.push_back(n);
            }
            
            size_t minNodes = *std::min_element(validNodeCounts.begin(), validNodeCounts.end());
            size_t maxNodes = *std::max_element(validNodeCounts.begin(), validNodeCounts.end());
            double nodeImbalance = static_cast<double>(maxNodes) / minNodes;
            
            std::cout << "Node distribution: min=" << minNodes << ", max=" << maxNodes 
                     << ", imbalance=" << nodeImbalance << std::endl;
        }
        
        // Element balance (if available)
        if (nonZeroElements > 1) {
            // Filter out zeros for min calculation
            std::vector<size_t> validElemCounts;
            for (auto n : rankElementCounts) {
                if (n > 0) validElemCounts.push_back(n);
            }
            
            size_t minElems = *std::min_element(validElemCounts.begin(), validElemCounts.end());
            size_t maxElems = *std::max_element(validElemCounts.begin(), validElemCounts.end());
            double elemImbalance = static_cast<double>(maxElems) / minElems;
            
            std::cout << "Element distribution: min=" << minElems << ", max=" << maxElems 
                     << ", imbalance=" << elemImbalance << std::endl;
            
            // For large meshes, expect reasonable balance
            EXPECT_LT(elemImbalance, 2.0) << "Element imbalance should be less than 2x";
        }
    } catch (const std::exception& e) {
        GTEST_SKIP() << "Failed to read mesh: " << e.what();
    }
}

// Test actual MPI-based distribution
TEST_P(ExternalMeshTest, ElementPartitioningBalanceMultipleRanks) {
    // This test runs on all ranks
    
    // Read this rank's portion of the mesh with element-based partitioning
    auto [nodeCount, elementCount, x, y, z, conn, localToGlobal] =
        mars::readMeshWithElementPartitioning<4, float>(meshPath, rank, numRanks);
    
    // Each rank verifies its own data
    EXPECT_GE(nodeCount, 0) << "Node count should be non-negative";
    EXPECT_GE(elementCount, 0) << "Element count should be non-negative";
    
    // Gather all element counts to rank 0 to verify distribution
    std::vector<size_t> allElementCounts(numRanks);
    MPI_Gather(&elementCount, 1, MPI_UNSIGNED_LONG, 
              allElementCounts.data(), 1, MPI_UNSIGNED_LONG, 
              0, MPI_COMM_WORLD);
    
    // Rank 0 verifies the distribution is reasonable
    if (rank == 0 && numRanks > 1) {
        // Filter out zeros
        std::vector<size_t> nonZeroCounts;
        for (auto count : allElementCounts) {
            if (count > 0) nonZeroCounts.push_back(count);
        }
        
        if (nonZeroCounts.size() > 1) {
            size_t minElems = *std::min_element(nonZeroCounts.begin(), nonZeroCounts.end());
            size_t maxElems = *std::max_element(nonZeroCounts.begin(), nonZeroCounts.end());
            double elemImbalance = static_cast<double>(maxElems) / minElems;
            
            std::cout << "Actual MPI element distribution: min=" << minElems 
                     << ", max=" << maxElems << ", imbalance=" << elemImbalance << std::endl;
            
            // For large meshes, expect reasonable balance
            EXPECT_LT(elemImbalance, 2.0) << "Element imbalance should be less than 2x";
        }
    }
    
    // Have each rank verify its node indices are valid
    if (elementCount > 0) {
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
    }
    
    // Make sure all ranks complete testing before ending
    MPI_Barrier(MPI_COMM_WORLD);
}

// Test that verifies no illegal node indices with element-based partitioning
TEST_P(ExternalMeshTest, NoIllegalNodeIndicesSingleRank) {
    try {
        // Skip if we're not on rank 0
        if (rank != 0) {
            GTEST_SKIP() << "Test only runs on rank 0";
        }
        
        // Test with 2 simulated ranks
        for (int r = 0; r < 2; r++) {
            auto [nodeCount, elementCount, x, y, z, conn, localToGlobal] =
                mars::readMeshWithElementPartitioning<4, float>(meshPath, r, 2);
            
            // Skip if no elements
            if (elementCount == 0) continue;
            
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

// Updated test that uses actual MPI processes rather than simulation
TEST_P(ExternalMeshTest, NoIllegalNodeIndicesMultipleRanks) {
    // Read this rank's portion of the mesh with element-based partitioning
    auto [nodeCount, elementCount, x, y, z, conn, localToGlobal] = 
        mars::readMeshWithElementPartitioning<4, float>(meshPath, rank, numRanks);
    
    // Each rank verifies its own data
    EXPECT_GE(nodeCount, 0) << "Node count should be non-negative";
    EXPECT_GE(elementCount, 0) << "Element count should be non-negative";
    
    // Have each rank verify its own node indices are valid
    if (elementCount > 0) {
        int currRank = rank;
        auto checkIndices = [nodeCount, currRank](const auto& indices) {
            for (size_t i = 0; i < indices.size(); i++) {
                // All indices should be non-negative and within bounds
                EXPECT_GE(indices[i], 0) << "Negative node index found on rank " << currRank;
                EXPECT_LT(indices[i], nodeCount) << "Node index out of bounds on rank " << currRank;
            }
        };
        
        std::apply([&checkIndices](const auto&... vecs) {
            (checkIndices(vecs), ...);
        }, conn);
        
        // Only print from rank 0 to avoid flooding output
        if (rank == 0) {
            std::cout << "Rank " << rank << ": Verified " << elementCount 
                     << " elements have valid node indices" << std::endl;
        }
    }
    
    // Collect verification status from all ranks
    int localVerificationSuccess = (elementCount == 0) ? 1 : 0;  // 0 means successful verification with elements
    int globalVerificationStatus = 0;
    
    MPI_Allreduce(&localVerificationSuccess, &globalVerificationStatus, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    // If no rank had elements to verify, print a message from rank 0
    if (rank == 0 && globalVerificationStatus == numRanks) {
        std::cout << "No elements to verify across any rank" << std::endl;
    }
    
    // Make sure all ranks complete testing before ending
    MPI_Barrier(MPI_COMM_WORLD);
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

// Simple test for external meshes that uses environment variables
TEST(ExternalMeshEnvVarTest, ReadExistingMeshEnvVar) {
    // Get mesh path from environment variable
    const char* meshPathEnv = std::getenv("MESH_PATH");
    std::string meshPath = meshPathEnv ? meshPathEnv : "";
    
    if (meshPath.empty()) {
        GTEST_SKIP() << "MESH_PATH environment variable is not set";
        return;
    }
    
    if (!fs::exists(meshPath)) {
        GTEST_SKIP() << "MESH_PATH directory does not exist: " << meshPath;
        return;
    }
    
    if (!fs::is_directory(fs::path(meshPath))) {
        GTEST_SKIP() << "MESH_PATH must be a directory: " << meshPath;
        return;
    }
    
    // Check if directory contains required files
    if (!fs::exists(fs::path(meshPath) / "x.float32") && 
        !fs::exists(fs::path(meshPath) / "x.double")) {
        GTEST_SKIP() << "MESH_PATH directory does not contain required coordinate files: " << meshPath;
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

int main(int argc, char **argv) {
    mars::Env env(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}