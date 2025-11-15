# Connectivity Management

Connectivity management handles element-to-element relationships and topological information for unstructured meshes, providing the foundation for mesh traversal and numerical computations.

## Purpose

Connectivity information enables:

- Mesh traversal algorithms
- Finite element assembly
- Boundary condition application
- Mesh adaptation operations
- Parallel data dependencies

## Key Components

### ConnectivityManager Class
```cpp
class ConnectivityManager {
public:
    std::vector<std::vector<size_t>> element_to_nodes;
    std::vector<std::vector<size_t>> node_to_elements;
    std::vector<std::vector<size_t>> element_to_faces;
    std::vector<std::vector<size_t>> face_to_elements;
    
    void build_connectivity(const std::vector<Element>& elements);
    std::vector<size_t> get_adjacent_elements(size_t element_id) const;
    std::vector<size_t> get_elements_at_node(size_t node_id) const;
};
```

### Topological Relationships
- Element-to-node connectivity (primary)
- Node-to-element relationships (derived)
- Face-based element adjacency
- Boundary face identification

## Usage Example

```cpp
// Build connectivity information
domain.build_connectivity();

// Access connectivity data
auto connectivity = domain.get_connectivity_manager();

// Get nodes of an element
auto elem_nodes = connectivity.element_to_nodes[element_id];
std::cout << "Element " << element_id << " has nodes: ";
for (auto node : elem_nodes) {
    std::cout << node << " ";
}

// Get elements sharing a node
auto node_elements = connectivity.get_elements_at_node(node_id);
std::cout << "Node " << node_id << " belongs to elements: ";
for (auto elem : node_elements) {
    std::cout << elem << " ";
}

// Find adjacent elements
auto adjacent = connectivity.get_adjacent_elements(element_id);
std::cout << "Adjacent elements: ";
for (auto adj : adjacent) {
    std::cout << adj << " ";
}
```

## Building Connectivity

### Element-to-Node Mapping
```cpp
void build_element_to_nodes(std::vector<std::vector<size_t>>& elem_to_nodes,
                          const std::vector<Element>& elements) {
    elem_to_nodes.resize(elements.size());
    for (size_t i = 0; i < elements.size(); ++i) {
        elem_to_nodes[i] = elements[i].nodes;
    }
}
```

### Node-to-Element Mapping
```cpp
void build_node_to_elements(std::vector<std::vector<size_t>>& node_to_elem,
                          const std::vector<Element>& elements,
                          size_t num_nodes) {
    node_to_elem.resize(num_nodes);
    for (size_t elem_id = 0; elem_id < elements.size(); ++elem_id) {
        for (size_t node : elements[elem_id].nodes) {
            node_to_elem[node].push_back(elem_id);
        }
    }
}
```

### Face-Based Adjacency
```cpp
struct Face {
    std::vector<size_t> nodes;
    size_t element_id;
    int local_face_index;
};

std::vector<std::vector<size_t>> build_face_adjacency(
    const std::vector<Element>& elements) {
    
    std::map<std::set<size_t>, std::vector<size_t>> face_to_elements;
    
    // Extract faces from all elements
    for (size_t elem_id = 0; elem_id < elements.size(); ++elem_id) {
        auto faces = extract_faces(elements[elem_id]);
        for (const auto& face : faces) {
            std::set<size_t> face_key(face.nodes.begin(), face.nodes.end());
            face_to_elements[face_key].push_back(elem_id);
        }
    }
    
    // Build element-to-element adjacency
    std::vector<std::vector<size_t>> adjacency(elements.size());
    for (const auto& [face_key, elem_list] : face_to_elements) {
        if (elem_list.size() == 2) {
            // Interior face
            adjacency[elem_list[0]].push_back(elem_list[1]);
            adjacency[elem_list[1]].push_back(elem_list[0]);
        }
        // Boundary faces have only one element
    }
    
    return adjacency;
}
```

## Applications

### Finite Element Assembly
```cpp
// Assemble global matrix using connectivity
void assemble_matrix(SparseMatrix& global_matrix,
                    const ConnectivityManager& connectivity,
                    const std::vector<Element>& elements) {
    
    for (size_t elem_id = 0; elem_id < elements.size(); ++elem_id) {
        // Get local element matrix
        auto local_matrix = compute_element_matrix(elements[elem_id]);
        
        // Get global node indices
        auto global_nodes = connectivity.element_to_nodes[elem_id];
        
        // Add to global matrix
        for (size_t i = 0; i < global_nodes.size(); ++i) {
            for (size_t j = 0; j < global_nodes.size(); ++j) {
                global_matrix.add(global_nodes[i], global_nodes[j], 
                                local_matrix(i,j));
            }
        }
    }
}
```

### Mesh Traversal
```cpp
// Breadth-first traversal
void traverse_mesh(size_t start_element,
                  const ConnectivityManager& connectivity) {
    
    std::queue<size_t> queue;
    std::vector<bool> visited(connectivity.element_to_nodes.size(), false);
    
    queue.push(start_element);
    visited[start_element] = true;
    
    while (!queue.empty()) {
        size_t current = queue.front(); queue.pop();
        
        // Process current element
        process_element(current);
        
        // Add unvisited neighbors
        for (size_t neighbor : connectivity.get_adjacent_elements(current)) {
            if (!visited[neighbor]) {
                visited[neighbor] = true;
                queue.push(neighbor);
            }
        }
    }
}
```

### Boundary Detection
```cpp
std::vector<size_t> find_boundary_faces(
    const ConnectivityManager& connectivity) {
    
    std::vector<size_t> boundary_faces;
    
    for (size_t face_id = 0; face_id < connectivity.face_to_elements.size(); ++face_id) {
        if (connectivity.face_to_elements[face_id].size() == 1) {
            // Face belongs to only one element (boundary)
            boundary_faces.push_back(face_id);
        }
    }
    
    return boundary_faces;
}
```

## Performance Optimizations

### Memory Layout
- Compressed storage for sparse connectivity
- Cache-aligned data structures
- Minimal memory overhead

### Query Optimization
- Pre-computed adjacency lists
- Fast lookup tables
- Lazy evaluation for derived relationships

### Parallel Considerations
- Local connectivity within partitions
- Halo connectivity for inter-partition relationships
- Consistent global indexing

## Integration with Other Components

- **Adjacency Structures**: Builds upon basic connectivity
- **Halo Management**: Uses connectivity for ghost element identification
- **Characteristic Sizes**: Leverages connectivity for quality metrics

## Error Handling

- Validate connectivity consistency
- Check for orphaned nodes/elements
- Detect topological errors
- Provide debugging utilities

## See Also

- [ElementDomain Overview](ElementDomain-Overview.md)
- [Adjacency Structures](Adjacency-Structures.md)
- [Halo Management](Halo-Management.md)