#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <tuple>
#include <filesystem>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>

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
 * Read mesh using element-based partitioning to ensure all elements have their nodes
 * @tparam N Number of nodes per element
 * @tparam RealType Type for coordinate data (float or double)
 * @param meshDir Directory containing mesh files
 * @param rank Current MPI rank
 * @param numRanks Total MPI ranks
 * @return Tuple containing: node count, element count, coordinates (x,y,z), and element connectivity
 */
template<int N, typename RealType=float>
inline auto readMeshWithElementPartitioning(const std::string& meshDir, int rank, int numRanks) {
    // Determine the file extension based on RealType
    std::string ext;
    if constexpr (std::is_same_v<RealType, float>) {
        ext = "float32";
    } else if constexpr (std::is_same_v<RealType, double>) {
        ext = "double";
    } else {
        throw std::runtime_error("Unsupported RealType for coordinate reading");
    }

    // STEP 1: Read element connectivity for this rank's portion
    // Open the first index file to determine total element count
    std::ifstream test_file(meshDir + "/i0.int32", std::ios::binary);
    if (!test_file) {
        throw std::runtime_error("Failed to open index file i0");
    }
    
    test_file.seekg(0, std::ios::end);
    size_t total_elements = test_file.tellg() / sizeof(int32_t);
    test_file.seekg(0, std::ios::beg);
    test_file.close();
    
    // Calculate this rank's element portion
    size_t elemPerRank = total_elements / numRanks;
    size_t elemStartIdx = rank * elemPerRank;
    size_t elemEndIdx = (rank == numRanks - 1) ? total_elements : elemStartIdx + elemPerRank;
    size_t elementCount = elemEndIdx - elemStartIdx;
    
    // Read element connectivity
    auto rawConnectivity = createNVectors<N>(elementCount);
    
    int i = 0;
    auto read_connectivity = [&](auto& vec) {
        std::ifstream idx_file(meshDir + "/i" + std::to_string(i++) + ".int32", std::ios::binary);
        if (!idx_file) {
            throw std::runtime_error("Failed to open index file i" + std::to_string(i-1));
        }
        
        idx_file.seekg(elemStartIdx * sizeof(int32_t));
        idx_file.read(reinterpret_cast<char*>(vec.data()), elementCount * sizeof(int32_t));
    };
    
    std::apply([&](auto&... vecs) { (read_connectivity(vecs), ...); }, rawConnectivity);
    
    // STEP 2: Identify unique nodes needed by this rank
    std::unordered_set<int> uniqueNodes;
    
    auto collect_nodes = [&](const auto& vec) {
        for (auto nodeId : vec) {
            uniqueNodes.insert(nodeId);
        }
    };
    
    std::apply([&](const auto&... vecs) { (collect_nodes(vecs), ...); }, rawConnectivity);
    
    // Create a sorted vector of unique node IDs for efficient file reading
    std::vector<int> neededNodes(uniqueNodes.begin(), uniqueNodes.end());
    std::sort(neededNodes.begin(), neededNodes.end());
    size_t nodeCount = neededNodes.size();
    
    // STEP 3: Create global-to-local node ID mapping
    std::unordered_map<int, int> globalToLocal;
    for (size_t i = 0; i < neededNodes.size(); i++) {
        globalToLocal[neededNodes[i]] = i;
    }
    
    // STEP 4: Read only needed node coordinates
    std::vector<RealType> x_data(nodeCount), y_data(nodeCount), z_data(nodeCount);
    
    std::ifstream x_file(meshDir + "/x." + ext, std::ios::binary);
    std::ifstream y_file(meshDir + "/y." + ext, std::ios::binary);
    std::ifstream z_file(meshDir + "/z." + ext, std::ios::binary);

    if (!x_file || !y_file || !z_file) {
        throw std::runtime_error("Failed to open coordinate files");
    }
    
    // Read coordinates for each needed node
    for (size_t i = 0; i < neededNodes.size(); i++) {
        int globalNodeId = neededNodes[i];
        
        x_file.seekg(globalNodeId * sizeof(RealType));
        y_file.seekg(globalNodeId * sizeof(RealType));
        z_file.seekg(globalNodeId * sizeof(RealType));
        
        x_file.read(reinterpret_cast<char*>(&x_data[i]), sizeof(RealType));
        y_file.read(reinterpret_cast<char*>(&y_data[i]), sizeof(RealType));
        z_file.read(reinterpret_cast<char*>(&z_data[i]), sizeof(RealType));
    }
    
    // STEP 5: Adjust connectivity to use local indices
    auto localConnectivity = createNVectors<N>(elementCount);
    
    auto adjust_indices = [&](const auto& global_vec, auto& local_vec) {
        for (size_t j = 0; j < elementCount; j++) {
            local_vec[j] = globalToLocal[global_vec[j]];
        }
    };
    
    // Apply the index adjustment to each connectivity vector
    std::apply([&](const auto&... global_vecs) {
        std::apply([&](auto&... local_vecs) {
            (adjust_indices(global_vecs, local_vecs), ...);
        }, localConnectivity);
    }, rawConnectivity);
    
    return std::make_tuple(nodeCount, elementCount, 
                           std::move(x_data), std::move(y_data), std::move(z_data), 
                           std::move(localConnectivity));
}

} // namespace mars