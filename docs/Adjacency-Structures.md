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

### AdjacencyData Struct (GPU)
```cpp
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
struct AdjacencyData {
    using DeviceVector = typename VectorSelector<KeyType, AcceleratorTag>::type;
    
    DeviceVector d_nodeToElementOffsets_;  // CSR row offsets
    DeviceVector d_nodeToElementList_;     // CSR column indices
    DeviceVector d_conn_local_ids_;        // SFC keys mapped to local IDs
    
    void buildNodeToElementMap(const ElementDomain& domain);
    void createElementToNodeLocalIdMap(const ElementDomain& domain);
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

### Algorithm Overview (GPU-based)
1. **Flatten Connectivity**: Convert tuple-of-vectors to flat device array (GPU kernel)
2. **Build Node-to-Element CSR**: Use Thrust `sort_by_key`, `reduce_by_key`, `exclusive_scan` on device
3. **Map SFC to Local IDs**: Convert SFC keys to dense local indices (GPU kernel)
4. **Store in Device Memory**: All CSR structures remain on GPU (`DeviceVector`)

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