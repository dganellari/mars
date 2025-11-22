#pragma once

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <iostream>
#include <map>

namespace mars {

/**
 * Simple MFEM mesh loader for tetrahedral meshes
 * 
 * Supports MFEM mesh format v1.0 with:
 * - dimension
 * - elements (attribute, geometry_type, vertex_indices)
 * - boundary (attribute, geometry_type, vertex_indices)
 * - vertices (count, dimension, coordinates)
 */
class MFEMMeshLoader {
public:
    std::vector<std::array<double, 3>> vertices;
    std::vector<std::array<int, 4>> elements;  // tetrahedra only
    std::vector<std::array<int, 3>> boundary_faces;
    std::vector<int> boundary_vertices;
    
    int dimension = 0;
    
    bool load(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Cannot open mesh file " << filename << std::endl;
            return false;
        }
        
        std::string line;
        
        // Read header
        std::getline(file, line);
        if (line.find("MFEM mesh") == std::string::npos) {
            std::cerr << "Error: Not a valid MFEM mesh file" << std::endl;
            return false;
        }
        
        // Read sections
        while (std::getline(file, line)) {
            // Skip comments and empty lines
            if (line.empty() || line[0] == '#') continue;
            
            if (line == "dimension") {
                std::getline(file, line);
                dimension = std::stoi(line);
                if (dimension != 3) {
                    std::cerr << "Error: Only 3D meshes supported" << std::endl;
                    return false;
                }
            }
            else if (line == "elements") {
                if (!read_elements(file)) return false;
            }
            else if (line == "boundary") {
                if (!read_boundary(file)) return false;
            }
            else if (line == "vertices") {
                if (!read_vertices(file)) return false;
            }
        }
        
        // Extract boundary vertices from boundary faces
        extract_boundary_vertices();
        
        // If no boundary faces found, detect geometrically
        if (boundary_faces.empty()) {
            detect_boundary_geometrically();
        }
        
        std::cout << "Loaded MFEM mesh:" << std::endl;
        std::cout << "  Vertices: " << vertices.size() << std::endl;
        std::cout << "  Elements: " << elements.size() << std::endl;
        std::cout << "  Boundary faces: " << boundary_faces.size() << std::endl;
        std::cout << "  Boundary vertices: " << boundary_vertices.size() << std::endl;
        std::cout << "  Interior vertices: " << (vertices.size() - boundary_vertices.size()) << std::endl;
        
        return true;
    }
    
    // Uniform mesh refinement - splits each tet into 8 child tets
    void uniform_refinement() {
        std::vector<std::array<double, 3>> new_vertices = vertices;
        std::vector<std::array<int, 4>> new_elements;
        
        // Build edge map for edge midpoints
        std::map<std::pair<int,int>, int> edge_to_midpoint;
        
        auto get_edge_midpoint = [&](int v0, int v1) -> int {
            if (v0 > v1) std::swap(v0, v1);
            auto edge = std::make_pair(v0, v1);
            
            if (edge_to_midpoint.find(edge) == edge_to_midpoint.end()) {
                // Create new vertex at edge midpoint
                int new_idx = new_vertices.size();
                std::array<double, 3> midpoint;
                for (int i = 0; i < 3; i++) {
                    midpoint[i] = 0.5 * (vertices[v0][i] + vertices[v1][i]);
                }
                new_vertices.push_back(midpoint);
                edge_to_midpoint[edge] = new_idx;
            }
            return edge_to_midpoint[edge];
        };
        
        // Refine each tetrahedron
        for (const auto& elem : elements) {
            int v0 = elem[0], v1 = elem[1], v2 = elem[2], v3 = elem[3];
            
            // Get edge midpoints (6 edges per tet)
            int m01 = get_edge_midpoint(v0, v1);
            int m02 = get_edge_midpoint(v0, v2);
            int m03 = get_edge_midpoint(v0, v3);
            int m12 = get_edge_midpoint(v1, v2);
            int m13 = get_edge_midpoint(v1, v3);
            int m23 = get_edge_midpoint(v2, v3);
            
            // 8 child tets (corner tets + 4 octahedral tets)
            new_elements.push_back({v0, m01, m02, m03});  // corner 0
            new_elements.push_back({m01, v1, m12, m13});  // corner 1
            new_elements.push_back({m02, m12, v2, m23});  // corner 2
            new_elements.push_back({m03, m13, m23, v3});  // corner 3
            
            // Central octahedron split into 4 tets
            new_elements.push_back({m01, m02, m03, m13});
            new_elements.push_back({m01, m02, m12, m13});
            new_elements.push_back({m02, m03, m13, m23});
            new_elements.push_back({m02, m12, m13, m23});
        }
        
        // Update mesh
        vertices = new_vertices;
        elements = new_elements;
        
        // Recompute boundary (simple approach: mark all original boundary vertices + midpoints on boundary edges)
        boundary_faces.clear();
        extract_boundary_vertices();
    }
    
private:
    bool read_elements(std::ifstream& file) {
        std::string line;
        std::getline(file, line);
        int num_elements = std::stoi(line);
        
        for (int i = 0; i < num_elements; i++) {
            std::getline(file, line);
            std::istringstream iss(line);
            
            int attribute, geom_type;
            iss >> attribute >> geom_type;
            
            // TETRAHEDRON = 4
            if (geom_type != 4) {
                std::cerr << "Warning: Skipping non-tetrahedral element" << std::endl;
                continue;
            }
            
            std::array<int, 4> elem;
            iss >> elem[0] >> elem[1] >> elem[2] >> elem[3];
            elements.push_back(elem);
        }
        
        return true;
    }
    
    bool read_boundary(std::ifstream& file) {
        std::string line;
        std::getline(file, line);
        int num_boundary = std::stoi(line);
        
        for (int i = 0; i < num_boundary; i++) {
            std::getline(file, line);
            std::istringstream iss(line);
            
            int attribute, geom_type;
            iss >> attribute >> geom_type;
            
            // TRIANGLE = 2
            if (geom_type == 2) {
                std::array<int, 3> face;
                iss >> face[0] >> face[1] >> face[2];
                boundary_faces.push_back(face);
            }
            // Ignore other boundary element types
        }
        
        return true;
    }
    
    bool read_vertices(std::ifstream& file) {
        std::string line;
        std::getline(file, line);
        int num_vertices = std::stoi(line);
        
        std::getline(file, line);
        int vertex_dim = std::stoi(line);
        
        if (vertex_dim != 3) {
            std::cerr << "Error: Only 3D vertices supported" << std::endl;
            return false;
        }
        
        for (int i = 0; i < num_vertices; i++) {
            std::getline(file, line);
            std::istringstream iss(line);
            
            std::array<double, 3> v;
            iss >> v[0] >> v[1] >> v[2];
            vertices.push_back(v);
        }
        
        return true;
    }
    
    void extract_boundary_vertices() {
        std::vector<bool> is_boundary(vertices.size(), false);
        
        // Mark all vertices that appear in boundary faces
        for (const auto& face : boundary_faces) {
            is_boundary[face[0]] = true;
            is_boundary[face[1]] = true;
            is_boundary[face[2]] = true;
        }
        
        // Collect boundary vertex indices
        for (int i = 0; i < vertices.size(); i++) {
            if (is_boundary[i]) {
                boundary_vertices.push_back(i);
            }
        }
        
        // Sort for consistency
        std::sort(boundary_vertices.begin(), boundary_vertices.end());
    }
    
    void detect_boundary_geometrically() {
        // For tetrahedral meshes, detect boundary faces geometrically
        // A face is boundary if it appears in exactly one tetrahedron
        
        // Face representation: sorted triple of node indices
        using Face = std::array<int, 3>;
        std::map<Face, int> face_count;
        
        // Count face occurrences
        for (const auto& elem : elements) {
            // Tet faces: (0,1,2), (0,1,3), (0,2,3), (1,2,3)
            std::vector<Face> faces = {
                {elem[0], elem[1], elem[2]},
                {elem[0], elem[1], elem[3]},
                {elem[0], elem[2], elem[3]},
                {elem[1], elem[2], elem[3]}
            };
            
            for (auto& face : faces) {
                std::sort(face.begin(), face.end());
                face_count[face]++;
            }
        }
        
        // Collect boundary faces (count == 1)
        for (const auto& [face, count] : face_count) {
            if (count == 1) {
                boundary_faces.push_back(face);
            }
        }
        
        // Extract boundary vertices from detected boundary faces
        extract_boundary_vertices();
    }
};

} // namespace mars
