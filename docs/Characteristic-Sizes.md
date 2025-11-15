# Characteristic Sizes

Characteristic sizes provide mesh quality metrics and sizing functions for unstructured meshes, enabling adaptive refinement and error estimation.

## Purpose

Characteristic sizes are essential for:

- Mesh quality assessment
- Adaptive mesh refinement
- Error estimation in numerical methods
- Load balancing decisions
- Convergence monitoring

## Key Components

### Characteristic Size Computation (GPU)
```cpp
// GPU kernel computes sizes for all elements in parallel
template<typename ElementTag, typename KeyType, typename RealType>
__global__ void computeCharacteristicSizesKernel(
    const RealType* x, const RealType* y, const RealType* z,
    const KeyType* indices0, /* ... other connectivity pointers */
    RealType* h,              // Output: nodal sizes
    int* nodeTetCount,        // Atomic counter for averaging
    int numElements);

// Stored in device memory
DeviceVector<RealType> d_h_;  // Characteristic sizes per node
```

### Size Metrics
- **Element Size**: Volume, surface area, or characteristic length
- **Nodal Size**: Averaged sizes at vertices
- **Edge Length**: Local mesh resolution

## Usage Example

```cpp
// Compute characteristic sizes
domain.compute_characteristic_sizes();

// Get size data
auto sizes = domain.get_characteristic_sizes();

// Element sizes
std::cout << "Average element size: " << 
    std::accumulate(sizes.element_sizes.begin(), sizes.element_sizes.end(), 0.0) / 
    sizes.element_sizes.size() << std::endl;

// Nodal sizes for mesh quality
for (size_t i = 0; i < domain.get_coordinates().size(); ++i) {
    double size = sizes.get_nodal_size(i);
    if (size < min_size_threshold) {
        std::cout << "Small element at node " << i << ": " << size << std::endl;
    }
}
```

## Size Computation Methods

### Edge-Based Sizing (GPU Kernel)
For each element, characteristic size computed from edge lengths:
```cpp
// Inside computeCharacteristicSizesKernel
// Decode SFC keys to physical coordinates
auto [x0, y0, z0] = decodeSfcToPhysical(indices0[elemIdx], box);
auto [x1, y1, z1] = decodeSfcToPhysical(indices1[elemIdx], box);
// ... repeat for all nodes

// Compute edge lengths and characteristic size
RealType dx = x1 - x0, dy = y1 - y0, dz = z1 - z0;
RealType edgeLen = sqrt(dx*dx + dy*dy + dz*dz);
RealType elemSize = /* average of all edge lengths */;

// Atomically accumulate at nodes for averaging
atomicAdd(&h[node0], elemSize);
atomicAdd(&nodeTetCount[node0], 1);
```

### Edge-Based Sizing
```cpp
double average_edge_length(const Element& elem, 
                          const CoordinateCache& coords) {
    double total_length = 0.0;
    int edge_count = 0;
    
    // Sum lengths of all edges
    for (size_t i = 0; i < elem.nodes.size(); ++i) {
        for (size_t j = i + 1; j < elem.nodes.size(); ++j) {
            Point p1 = coords.get_coordinate(elem.nodes[i]);
            Point p2 = coords.get_coordinate(elem.nodes[j]);
            total_length += distance(p1, p2);
            edge_count++;
        }
    }
    
    return total_length / edge_count;
}
```

### Nodal Size Computation
```cpp
void compute_nodal_sizes(std::vector<double>& nodal_sizes,
                        const std::vector<Element>& elements,
                        const std::vector<double>& element_sizes) {
    std::vector<double> sum_sizes(nodal_sizes.size(), 0.0);
    std::vector<int> count(nodal_sizes.size(), 0);
    
    // Accumulate element sizes at nodes
    for (size_t i = 0; i < elements.size(); ++i) {
        for (size_t node : elements[i].nodes) {
            sum_sizes[node] += element_sizes[i];
            count[node]++;
        }
    }
    
    // Average at each node
    for (size_t i = 0; i < nodal_sizes.size(); ++i) {
        if (count[i] > 0) {
            nodal_sizes[i] = sum_sizes[i] / count[i];
        }
    }
}
```

## Applications

### Mesh Quality Assessment
```cpp
// Check mesh quality using size ratios
void assess_mesh_quality(const CharacteristicSize& sizes,
                        const AdjacencyData& adjacency) {
    for (size_t elem_id = 0; elem_id < sizes.element_sizes.size(); ++elem_id) {
        double elem_size = sizes.get_element_size(elem_id);
        auto neighbors = adjacency.get_neighbors(elem_id);
        
        double max_ratio = 1.0;
        for (size_t neighbor : neighbors) {
            double neighbor_size = sizes.get_element_size(neighbor);
            double ratio = std::max(elem_size, neighbor_size) / 
                          std::min(elem_size, neighbor_size);
            max_ratio = std::max(max_ratio, ratio);
        }
        
        if (max_ratio > quality_threshold) {
            std::cout << "Poor quality element: " << elem_id 
                      << " (ratio: " << max_ratio << ")" << std::endl;
        }
    }
}
```

### Adaptive Refinement
```cpp
// Identify elements for refinement
std::vector<size_t> elements_to_refine;
for (size_t i = 0; i < sizes.element_sizes.size(); ++i) {
    if (sizes.get_element_size(i) > max_desired_size) {
        elements_to_refine.push_back(i);
    }
}

// Identify elements for coarsening
std::vector<size_t> elements_to_coarsen;
for (size_t i = 0; i < sizes.element_sizes.size(); ++i) {
    if (sizes.get_element_size(i) < min_desired_size) {
        elements_to_coarsen.push_back(i);
    }
}
```

## Performance Considerations

### Computation Cost
- O(N) for element size computation
- O(N) for nodal size averaging
- Incremental updates for mesh adaptation

### Memory Usage
- Additional storage for size arrays
- Typically small overhead (2-3x element count)

### Parallel Computation
- Local computation with halo exchange
- Consistent sizing across partition boundaries
- Load-balanced size calculations

## Integration with Other Components

- **Adjacency Structures**: For neighbor-based quality checks
- **Halo Management**: For consistent sizing across partitions
- **GPU Acceleration**: Parallel size computations

## See Also

- [ElementDomain Overview](ElementDomain-Overview.md)
- [Adjacency Structures](Adjacency-Structures.md)
- [Coordinate Caching](Coordinate-Caching.md)