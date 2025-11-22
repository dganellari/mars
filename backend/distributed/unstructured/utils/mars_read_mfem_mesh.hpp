#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <tuple>
#include <map>
#include <sstream>
#include <algorithm>
#include "mars.hpp"

namespace mars
{

/**
 * @brief Read MFEM mesh format and return data compatible with mars::ElementDomain
 * 
 * This follows the same pattern as readMeshWithElementPartitioning() but reads
 * MFEM mesh files instead of binary format. Returns coordinates and connectivity
 * in the same SoA format expected by ElementDomain.
 * 
 * @tparam N Number of nodes per element (4 for tetrahedra, 8 for hexahedra, etc.)
 * @tparam RealType Type for coordinate data (float or double)
 * @tparam KeyType Type for connectivity indices
 * @param meshFile Path to MFEM mesh file
 * @param rank Current MPI rank
 * @param numRanks Total MPI ranks
 * @return Tuple containing: node count, element count, x_coords, y_coords, z_coords, connectivity tuple, local_to_global_map
 */
template<int N, typename RealType = float, typename KeyType = unsigned>
inline auto readMFEMMeshWithElementPartitioning(const std::string& meshFile, int rank, int numRanks)
{
    std::ifstream file(meshFile);
    if (!file) {
        throw std::runtime_error("Failed to open MFEM mesh file: " + meshFile);
    }

    // Helper to create N-tuple of vectors
    auto createNVectors = [](size_t size) {
        if constexpr (N == 3) {
            return std::tuple<std::vector<KeyType>, std::vector<KeyType>, std::vector<KeyType>>(
                std::vector<KeyType>(size), std::vector<KeyType>(size), std::vector<KeyType>(size));
        } else if constexpr (N == 4) {
            return std::tuple<std::vector<KeyType>, std::vector<KeyType>, std::vector<KeyType>, std::vector<KeyType>>(
                std::vector<KeyType>(size), std::vector<KeyType>(size), std::vector<KeyType>(size), std::vector<KeyType>(size));
        } else if constexpr (N == 8) {
            return std::tuple<std::vector<KeyType>, std::vector<KeyType>, std::vector<KeyType>, std::vector<KeyType>,
                              std::vector<KeyType>, std::vector<KeyType>, std::vector<KeyType>, std::vector<KeyType>>(
                std::vector<KeyType>(size), std::vector<KeyType>(size), std::vector<KeyType>(size), std::vector<KeyType>(size),
                std::vector<KeyType>(size), std::vector<KeyType>(size), std::vector<KeyType>(size), std::vector<KeyType>(size));
        }
    };

    std::string line;
    int dimension = 0;
    size_t total_elements = 0;
    size_t total_vertices = 0;

    // Parse header to find dimension
    while (std::getline(file, line)) {
        if (line.find("dimension") != std::string::npos) {
            file >> dimension;
            std::getline(file, line); // consume rest of line
            break;
        }
    }

    if (dimension == 0) {
        throw std::runtime_error("Failed to read dimension from MFEM mesh");
    }

    // Parse elements section
    while (std::getline(file, line)) {
        if (line.find("elements") != std::string::npos) {
            file >> total_elements;
            std::getline(file, line); // consume rest of line
            break;
        }
    }

    if (total_elements == 0) {
        throw std::runtime_error("No elements found in MFEM mesh");
    }

    // Parse boundary section (skip for now, but consume the data)
    while (std::getline(file, line)) {
        if (line.find("boundary") != std::string::npos) {
            size_t num_boundary;
            file >> num_boundary;
            std::getline(file, line);

            // Skip boundary elements
            for (size_t i = 0; i < num_boundary; i++) {
                std::getline(file, line);
            }
            break;
        }
    }

    // Parse vertices section - MUST BE READ BEFORE processing elements
    while (std::getline(file, line)) {
        if (line.find("vertices") != std::string::npos) {
            file >> total_vertices;
            std::getline(file, line); // consume rest of line
            break;
        }
    }

    if (total_vertices == 0) {
        throw std::runtime_error("No vertices found in MFEM mesh");
    }

    // STEP 1: Calculate this rank's element portion and read elements
    size_t elemPerRank = total_elements / numRanks;
    size_t elemStartIdx = rank * elemPerRank;
    size_t elemEndIdx = (rank == numRanks - 1) ? total_elements : elemStartIdx + elemPerRank;
    size_t elementCount = elemEndIdx - elemStartIdx;

    // Go back and read elements
    file.clear();
    file.seekg(0, std::ios::beg);
    
    // Skip to elements section
    while (std::getline(file, line)) {
        if (line.find("elements") != std::string::npos) {
            size_t num_elems;
            file >> num_elems;
            std::getline(file, line);
            break;
        }
    }

    // Read all elements but only store this rank's portion
    auto rawConnectivity = createNVectors(elementCount);
    std::vector<int> element_attributes(elementCount);

    for (size_t i = 0; i < total_elements; i++) {
        int attr, geom_type;
        file >> attr >> geom_type;

        std::vector<KeyType> nodes(N);
        for (int j = 0; j < N; j++) {
            int node_id;
            file >> node_id;
            nodes[j] = static_cast<KeyType>(node_id);
        }

        // Only store if this element belongs to current rank
        if (i >= elemStartIdx && i < elemEndIdx) {
            size_t local_idx = i - elemStartIdx;
            element_attributes[local_idx] = attr;

            // Store in SoA format
            if constexpr (N == 4) {
                std::get<0>(rawConnectivity)[local_idx] = nodes[0];
                std::get<1>(rawConnectivity)[local_idx] = nodes[1];
                std::get<2>(rawConnectivity)[local_idx] = nodes[2];
                std::get<3>(rawConnectivity)[local_idx] = nodes[3];
            } else if constexpr (N == 8) {
                for (int j = 0; j < 8; j++) {
                    if (j == 0) std::get<0>(rawConnectivity)[local_idx] = nodes[j];
                    else if (j == 1) std::get<1>(rawConnectivity)[local_idx] = nodes[j];
                    else if (j == 2) std::get<2>(rawConnectivity)[local_idx] = nodes[j];
                    else if (j == 3) std::get<3>(rawConnectivity)[local_idx] = nodes[j];
                    else if (j == 4) std::get<4>(rawConnectivity)[local_idx] = nodes[j];
                    else if (j == 5) std::get<5>(rawConnectivity)[local_idx] = nodes[j];
                    else if (j == 6) std::get<6>(rawConnectivity)[local_idx] = nodes[j];
                    else if (j == 7) std::get<7>(rawConnectivity)[local_idx] = nodes[j];
                }
            } else if constexpr (N == 3) {
                std::get<0>(rawConnectivity)[local_idx] = nodes[0];
                std::get<1>(rawConnectivity)[local_idx] = nodes[1];
                std::get<2>(rawConnectivity)[local_idx] = nodes[2];
            }
        }
    }

    // STEP 2: Identify unique nodes needed by this rank
    std::vector<KeyType> uniqueNodes;
    {
        std::vector<bool> nodeUsed(total_vertices, false);

        // First pass: mark which nodes are used
        auto mark_nodes = [&](const auto& vec) {
            for (auto nodeId : vec) {
                if (nodeId < nodeUsed.size()) {
                    nodeUsed[nodeId] = true;
                }
            }
        };

        std::apply([&](const auto&... vecs) { (mark_nodes(vecs), ...); }, rawConnectivity);

        // Collect used nodes
        for (size_t i = 0; i < nodeUsed.size(); i++) {
            if (nodeUsed[i]) {
                uniqueNodes.push_back(static_cast<KeyType>(i));
            }
        }
    }

    std::sort(uniqueNodes.begin(), uniqueNodes.end());
    size_t nodeCount = uniqueNodes.size();

    // STEP 3: Create global-to-local node ID mapping
    std::map<KeyType, KeyType> globalToLocal;
    for (size_t i = 0; i < uniqueNodes.size(); i++) {
        globalToLocal[uniqueNodes[i]] = static_cast<KeyType>(i);
    }

    // STEP 4: Read only needed node coordinates
    // Go back to vertices section
    file.clear();
    file.seekg(0, std::ios::beg);
    
    while (std::getline(file, line)) {
        if (line.find("vertices") != std::string::npos) {
            file >> total_vertices;
            std::getline(file, line);
            break;
        }
    }

    // Read vertex dimension (should match mesh dimension)
    int vertex_dim;
    file >> vertex_dim;
    std::getline(file, line);

    if (vertex_dim != dimension) {
        throw std::runtime_error("Vertex dimension does not match mesh dimension");
    }

    // STEP 4: Read only needed node coordinates
    std::vector<RealType> x_data(nodeCount), y_data(nodeCount), z_data(nodeCount);
    
    // Create a map for quick lookup
    std::map<KeyType, size_t> nodeToLocalIdx;
    for (size_t i = 0; i < uniqueNodes.size(); i++) {
        nodeToLocalIdx[uniqueNodes[i]] = i;
    }

    // Read all vertices but only store needed ones
    for (size_t i = 0; i < total_vertices; i++) {
        double x, y, z = 0.0;
        file >> x >> y;
        if (dimension == 3) {
            file >> z;
        }

        auto it = nodeToLocalIdx.find(static_cast<KeyType>(i));
        if (it != nodeToLocalIdx.end()) {
            size_t local_idx = it->second;
            x_data[local_idx] = static_cast<RealType>(x);
            y_data[local_idx] = static_cast<RealType>(y);
            z_data[local_idx] = static_cast<RealType>(z);
        }
    }

    // STEP 5: Adjust connectivity to use local indices
    auto localConnectivity = createNVectors(elementCount);

    auto adjust_indices = [&](const auto& global_vec, auto& local_vec) {
        for (size_t j = 0; j < elementCount; j++) {
            local_vec[j] = globalToLocal[global_vec[j]];
        }
    };

    // Apply the index adjustment to each connectivity vector
    std::apply(
        [&](const auto&... global_vecs) {
            std::apply([&](auto&... local_vecs) { (adjust_indices(global_vecs, local_vecs), ...); }, localConnectivity);
        },
        rawConnectivity);

    // STEP 5: Read boundary elements to collect boundary nodes
    std::unordered_set<KeyType> boundaryNodes;
    
    // Go back to boundary section
    file.clear();
    file.seekg(0, std::ios::beg);
    
    while (std::getline(file, line)) {
        if (line.find("boundary") != std::string::npos) {
            size_t num_boundary;
            file >> num_boundary;
            std::getline(file, line);
            
            // Read all boundary elements and collect their nodes
            for (size_t i = 0; i < num_boundary; i++) {
                int attr, geom_type;
                file >> attr >> geom_type;
                
                // Read boundary element nodes (triangles for 3D tets, edges for 2D tris)
                if (dimension == 3) {
                    // Triangular boundary elements
                    KeyType n0, n1, n2;
                    file >> n0 >> n1 >> n2;
                    boundaryNodes.insert(n0);
                    boundaryNodes.insert(n1);
                    boundaryNodes.insert(n2);
                } else if (dimension == 2) {
                    // Edge boundary elements
                    KeyType n0, n1;
                    file >> n0 >> n1;
                    boundaryNodes.insert(n0);
                    boundaryNodes.insert(n1);
                }
            }
            break;
        }
    }

    file.close();

    return std::make_tuple(nodeCount, elementCount, std::move(x_data), std::move(y_data), std::move(z_data),
                           std::move(localConnectivity), std::move(uniqueNodes), std::move(boundaryNodes));
}

/**
 * @brief Apply uniform refinement to MFEM mesh data
 * 
 * Performs uniform refinement by splitting each element into smaller elements.
 * For tetrahedra: each tet splits into 8 child tets via edge midpoint subdivision.
 * 
 * @tparam N Number of nodes per element
 * @tparam RealType Type for coordinate data
 * @tparam KeyType Type for connectivity indices
 * @param x_data X coordinates (will be modified)
 * @param y_data Y coordinates (will be modified)
 * @param z_data Z coordinates (will be modified)
 * @param connectivity Element connectivity (will be modified)
 * @param nodeCount Current node count (will be updated)
 * @param elementCount Current element count (will be updated)
 */
// template<int N, typename RealType, typename KeyType>
// inline void uniformRefineMFEMMesh(std::vector<RealType>& x_data,
//                                   std::vector<RealType>& y_data,
//                                   std::vector<RealType>& z_data,
//                                   auto& connectivity,
//                                   size_t& nodeCount,
//                                   size_t& elementCount)
// {
//     if constexpr (N != 4) {
//         throw std::runtime_error("Uniform refinement currently only supported for tetrahedra (N=4)");
//     }

//     size_t old_elem_count = elementCount;
//     size_t old_node_count = nodeCount;

//     // Extract old connectivity
//     auto& i0 = std::get<0>(connectivity);
//     auto& i1 = std::get<1>(connectivity);
//     auto& i2 = std::get<2>(connectivity);
//     auto& i3 = std::get<3>(connectivity);

//     // Map to store edge midpoints: (min_node, max_node) -> new_node_id
//     std::map<std::pair<KeyType, KeyType>, KeyType> edge_to_midpoint;
//     KeyType next_node_id = old_node_count;

//     // Helper to get or create edge midpoint
//     auto get_edge_midpoint = [&](KeyType n0, KeyType n1) -> KeyType {
//         KeyType min_n = std::min(n0, n1);
//         KeyType max_n = std::max(n0, n1);
//         auto edge = std::make_pair(min_n, max_n);

//         auto it = edge_to_midpoint.find(edge);
//         if (it != edge_to_midpoint.end()) {
//             return it->second;
//         }

//         // Create new midpoint node
//         KeyType mid_id = next_node_id++;
//         edge_to_midpoint[edge] = mid_id;

//         // Add coordinates
//         x_data.push_back((x_data[n0] + x_data[n1]) / 2);
//         y_data.push_back((y_data[n0] + y_data[n1]) / 2);
//         z_data.push_back((z_data[n0] + z_data[n1]) / 2);

//         return mid_id;
//     };

//     // New connectivity storage (8 children per parent)
//     std::vector<KeyType> new_i0, new_i1, new_i2, new_i3;
//     new_i0.reserve(old_elem_count * 8);
//     new_i1.reserve(old_elem_count * 8);
//     new_i2.reserve(old_elem_count * 8);
//     new_i3.reserve(old_elem_count * 8);

//     // Refine each tetrahedron
//     for (size_t e = 0; e < old_elem_count; e++) {
//         KeyType v0 = i0[e], v1 = i1[e], v2 = i2[e], v3 = i3[e];

//         // Get edge midpoints (6 edges in a tet)
//         KeyType m01 = get_edge_midpoint(v0, v1);
//         KeyType m02 = get_edge_midpoint(v0, v2);
//         KeyType m03 = get_edge_midpoint(v0, v3);
//         KeyType m12 = get_edge_midpoint(v1, v2);
//         KeyType m13 = get_edge_midpoint(v1, v3);
//         KeyType m23 = get_edge_midpoint(v2, v3);

//         // Create 8 child tetrahedra (4 corner + 4 octahedral)
//         // Corner tets at original vertices
//         new_i0.push_back(v0);  new_i1.push_back(m01); new_i2.push_back(m02); new_i3.push_back(m03);
//         new_i0.push_back(m01); new_i1.push_back(v1);  new_i2.push_back(m12); new_i3.push_back(m13);
//         new_i0.push_back(m02); new_i1.push_back(m12); new_i2.push_back(v2);  new_i3.push_back(m23);
//         new_i0.push_back(m03); new_i1.push_back(m13); new_i2.push_back(m23); new_i3.push_back(v3);

//         // Octahedral split into 4 tets (interior)
//         new_i0.push_back(m01); new_i1.push_back(m02); new_i2.push_back(m03); new_i3.push_back(m12);
//         new_i0.push_back(m01); new_i1.push_back(m02); new_i2.push_back(m12); new_i3.push_back(m13);
//         new_i0.push_back(m02); new_i1.push_back(m03); new_i2.push_back(m12); new_i3.push_back(m23);
//         new_i0.push_back(m02); new_i1.push_back(m12); new_i2.push_back(m13); new_i3.push_back(m23);
//     }

//     // Update connectivity with refined mesh
//     i0 = std::move(new_i0);
//     i1 = std::move(new_i1);
//     i2 = std::move(new_i2);
//     i3 = std::move(new_i3);

//     // Update counts
//     nodeCount = x_data.size();
//     elementCount = i0.size();
// }

} // namespace mars
