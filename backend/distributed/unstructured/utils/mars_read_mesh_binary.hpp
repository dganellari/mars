#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <tuple>
#include <filesystem>

namespace mars {

// Helper function to create a tuple of N vector<int>
template<int N>
auto createNVectors(size_t size) {
    if constexpr (N == 3) {
        return std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>(
            std::vector<int>(size), std::vector<int>(size), std::vector<int>(size));
    } else if constexpr (N == 4) {
        return std::tuple<std::vector<int>, std::vector<int>, std::vector<int>, std::vector<int>>(
            std::vector<int>(size), std::vector<int>(size), std::vector<int>(size), std::vector<int>(size));
    } else if constexpr (N == 8) {
        return std::tuple<std::vector<int>, std::vector<int>, std::vector<int>, std::vector<int>,
                          std::vector<int>, std::vector<int>, std::vector<int>, std::vector<int>>(
            std::vector<int>(size), std::vector<int>(size), std::vector<int>(size), std::vector<int>(size),
            std::vector<int>(size), std::vector<int>(size), std::vector<int>(size), std::vector<int>(size));
    }
}

/**
 * Read binary mesh coordinates from separate x, y, z files
 * @param meshDir Directory containing mesh files
 * @param rank Current MPI rank
 * @param numRanks Total MPI ranks
 * @return Tuple containing: node count, node start index, coordinates (x,y,z)
 */
template<typename RealType=float>
inline std::tuple<size_t, size_t, std::vector<RealType>, std::vector<RealType>, std::vector<RealType>> 
readMeshCoordinatesBinary(const std::string& meshDir, int rank, int numRanks) {
    // Determine the file extension based on RealType
    std::string ext;
    if constexpr (std::is_same_v<RealType, float>) {
        ext = "float32";
    } else if constexpr (std::is_same_v<RealType, double>) {
        ext = "double";
    } else {
        throw std::runtime_error("Unsupported RealType for coordinate reading");
    }
    
    // Open the mesh directory and read coordinate files
    std::string x_path = meshDir + "/x." + ext;
    std::string y_path = meshDir + "/y." + ext;
    std::string z_path = meshDir + "/z." + ext;
    
    std::ifstream x_file(x_path, std::ios::binary);
    std::ifstream y_file(y_path, std::ios::binary);
    std::ifstream z_file(z_path, std::ios::binary);

    if (!x_file || !y_file || !z_file) {
        throw std::runtime_error("Failed to open coordinate files: " + x_path);
    }

    // Get file size to determine node count
    x_file.seekg(0, std::ios::end);
    size_t total_nodes = x_file.tellg() / sizeof(RealType);
    x_file.seekg(0, std::ios::beg);

    // Calculate this rank's portion
    size_t nodePerRank = total_nodes / numRanks;
    size_t nodeStartIdx = rank * nodePerRank;
    size_t nodeEndIdx = (rank == numRanks - 1) ? total_nodes : nodeStartIdx + nodePerRank;
    size_t nodeCount = nodeEndIdx - nodeStartIdx;

    // Read coordinate data
    std::vector<RealType> x_data(nodeCount), y_data(nodeCount), z_data(nodeCount);

    x_file.seekg(nodeStartIdx * sizeof(RealType));
    y_file.seekg(nodeStartIdx * sizeof(RealType));
    z_file.seekg(nodeStartIdx * sizeof(RealType));

    x_file.read(reinterpret_cast<char*>(x_data.data()), nodeCount * sizeof(RealType));
    y_file.read(reinterpret_cast<char*>(y_data.data()), nodeCount * sizeof(RealType));
    z_file.read(reinterpret_cast<char*>(z_data.data()), nodeCount * sizeof(RealType));

    // Use std::move to avoid copying
    return {nodeCount, nodeStartIdx, std::move(x_data), std::move(y_data), std::move(z_data)};
}

/**
 * Read binary element connectivity from separate index files
 * @tparam N Number of nodes per element
 * @param meshDir Directory containing mesh files
 * @param nodeStartIdx Starting node index for this rank
 * @param rank Current MPI rank
 * @param numRanks Total MPI ranks
 * @return Tuple of connectivity vectors, one per node index
 */
template<int N>
auto readMeshConnectivityBinaryTuple(const std::string& meshDir, size_t nodeStartIdx, int rank, int numRanks) {
    // Calculate element range for this rank
    std::ifstream test_file(meshDir + "/i0.int32", std::ios::binary);
    if (!test_file) {
        throw std::runtime_error("Failed to open index file i0");
    }
    
    test_file.seekg(0, std::ios::end);
    size_t total_elements = test_file.tellg() / sizeof(int32_t);
    test_file.seekg(0, std::ios::beg);
    test_file.close();
    
    size_t elemPerRank = total_elements / numRanks;
    size_t elemStartIdx = rank * elemPerRank;
    size_t elemEndIdx = (rank == numRanks - 1) ? total_elements : elemStartIdx + elemPerRank;
    size_t elementCount = elemEndIdx - elemStartIdx;
    
    // Create the tuple of vectors
    auto result = createNVectors<N>(elementCount);
    
    // Read data into each vector in the tuple
    int i = 0;
    auto read_into_tuple = [&](auto& vec) {
        std::ifstream idx_file(meshDir + "/i" + std::to_string(i++) + ".int32", std::ios::binary);
        if (!idx_file) {
            throw std::runtime_error("Failed to open index file i" + std::to_string(i-1));
        }
        
        std::vector<int32_t> temp_indices(elementCount);
        idx_file.seekg(elemStartIdx * sizeof(int32_t));
        idx_file.read(reinterpret_cast<char*>(temp_indices.data()), elementCount * sizeof(int32_t));
        
        // Adjust indices for local numbering
        for (size_t j = 0; j < elementCount; ++j) {
            vec[j] = temp_indices[j] - nodeStartIdx;
        }
    };
    
    std::apply([&](auto&... vecs) { (read_into_tuple(vecs), ...); }, result);
    
    return std::make_tuple(elementCount, std::move(result));
}

/**
 * Read binary element connectivity from separate index files
 * @param meshDir Directory containing mesh files
 * @param nodesPerElement Number of nodes per element
 * @param nodeStartIdx Starting node index for this rank
 * @param rank Current MPI rank
 * @param numRanks Total MPI ranks
 * @return Tuple containing: element count, element indices arrays
 */
inline std::tuple<size_t, std::vector<std::vector<int>>> 
readMeshConnectivityBinary(const std::string& meshDir, int nodesPerElement, size_t nodeStartIdx, int rank, int numRanks) {
    // Open all index files based on element type
    std::vector<std::ifstream> index_files;
    for (int i = 0; i < nodesPerElement; ++i) {
        index_files.emplace_back(meshDir + "/i" + std::to_string(i) + ".int32", std::ios::binary);
        if (!index_files.back()) {
            throw std::runtime_error("Failed to open index file i" + std::to_string(i));
        }
    }

    // Get element count
    index_files[0].seekg(0, std::ios::end);
    size_t total_elements = index_files[0].tellg() / sizeof(int32_t);
    index_files[0].seekg(0, std::ios::beg);

    // Calculate this rank's element portion
    size_t elemPerRank = total_elements / numRanks;
    size_t elemStartIdx = rank * elemPerRank;
    size_t elemEndIdx = (rank == numRanks - 1) ? total_elements : elemStartIdx + elemPerRank;
    size_t elementCount = elemEndIdx - elemStartIdx;

    // Read element connectivity
    std::vector<std::vector<int>> indices(nodesPerElement, std::vector<int>(elementCount));
    std::vector<int32_t> temp_indices(elementCount);

    for (int i = 0; i < nodesPerElement; ++i) {
        index_files[i].seekg(elemStartIdx * sizeof(int32_t));
        index_files[i].read(reinterpret_cast<char*>(temp_indices.data()), elementCount * sizeof(int32_t));

        // Adjust indices for local numbering
        for (size_t j = 0; j < elementCount; ++j) {
            indices[i][j] = temp_indices[j] - nodeStartIdx;
        }
    }

    return {elementCount, indices};
}

} // namespace mars