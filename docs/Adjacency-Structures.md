# Adjacency Structures

Adjacency structures provide efficient neighbor finding and connectivity information for unstructured mesh elements, enabling fast traversal and local operations.

## Purpose

Adjacency information is crucial for:

- Finite element assembly
- Mesh traversal algorithms
- Boundary condition application
- Mesh quality assessment
- Parallel communication patterns

## Key Components

### AdjacencyData Class
```cpp
class AdjacencyData {
public:
    std::vector<std::vector<size_t>> element_to_element;
    std::vector<std::vector<size_t>> element_to_faces;
    std::vector<std::vector<size_t>> face_to_elements;
    
    void build_adjacency(const std::vector<Element>& elements);
    std::vector<size_t> get_neighbors(size_t element_id) const;
};
```

### Face-Based Adjacency
- Identifies shared faces between elements
- Builds element-to-element connectivity
- Supports arbitrary polyhedral elements

## Usage Example

```cpp
// Build adjacency for domain
domain.build_adjacency();

// Get adjacency data
auto adjacency = domain.get_adjacency_data();

// Find neighbors of element 0
auto neighbors = adjacency.get_neighbors(0);
std::cout << "Element 0 has " << neighbors.size() << " neighbors: ";
for (auto n : neighbors) {
    std::cout << n << " ";
}

// Access element-to-element connectivity
const auto& elem_to_elem = adjacency.element_to_element;
for (size_t i = 0; i < elem_to_elem.size(); ++i) {
    std::cout << "Element " << i << " neighbors: ";
    for (auto n : elem_to_elem[i]) {
        std::cout << n << " ";
    }
    std::cout << std::endl;
}
```

## Building Adjacency

### Algorithm Overview
1. **Face Extraction**: Get all faces for each element
2. **Face Matching**: Find elements sharing each face
3. **Connectivity Building**: Create element-to-element maps
4. **Boundary Detection**: Identify boundary faces

### Face Representation
```cpp
struct Face {
    std::vector<size_t> nodes;  // Node indices
    size_t element_id;          // Owning element
    size_t local_face_id;       // Face index within element
};
```

## Performance Optimizations

### Memory Layout
- Compressed sparse row (CSR) format for adjacency
- Minimal memory overhead (O(number of adjacencies))
- Cache-friendly data structures

### Lazy Construction
- Adjacency built only when first requested
- Incremental updates for mesh adaptation
- Parallel construction algorithms

## Applications

### Finite Element Assembly
```cpp
// Use adjacency for element loop
for (size_t elem_id = 0; elem_id < domain.size(); ++elem_id) {
    auto neighbors = adjacency.get_neighbors(elem_id);
    
    // Assemble local matrix with neighboring elements
    for (auto neighbor_id : neighbors) {
        // Compute element-element coupling
        assemble_coupling(elem_id, neighbor_id);
    }
}
```

### Mesh Traversal
```cpp
// Breadth-first search using adjacency
std::queue<size_t> queue;
std::vector<bool> visited(domain.size(), false);

queue.push(start_element);
visited[start_element] = true;

while (!queue.empty()) {
    size_t current = queue.front(); queue.pop();
    
    // Process current element
    process_element(current);
    
    // Add unvisited neighbors
    for (auto neighbor : adjacency.get_neighbors(current)) {
        if (!visited[neighbor]) {
            visited[neighbor] = true;
            queue.push(neighbor);
        }
    }
}
```

## Integration with Other Components

- **Halo Management**: Uses adjacency to identify ghost elements
- **Characteristic Sizes**: Leverages adjacency for local mesh quality
- **GPU Acceleration**: Adjacency structures optimized for GPU access

## See Also

- [ElementDomain Overview](ElementDomain-Overview.md)
- [Halo Management](Halo-Management.md)
- [Connectivity Management](Connectivity-Management.md)