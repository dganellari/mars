# Coordinate Caching

Coordinate caching provides efficient storage and access to mesh vertex coordinates, optimizing memory usage and computational performance for unstructured meshes.

## Purpose

Coordinate caching addresses performance bottlenecks in:

- Repeated coordinate lookups during element processing
- Memory bandwidth limitations
- Cache locality optimization
- GPU memory management

## Key Components

### CoordinateCache Struct (GPU SoA)
```cpp
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
struct CoordinateCache {
    using DeviceVector = typename VectorSelector<RealType, AcceleratorTag>::type;
    
    DeviceVector d_node_x_;  // X coordinates in device memory
    DeviceVector d_node_y_;  // Y coordinates in device memory
    DeviceVector d_node_z_;  // Z coordinates in device memory
    
    void cacheNodeCoordinates(const ElementDomain& domain);
};
```

### Memory Layout (GPU SoA Only)
- **Structure of Arrays (SoA)**: Three separate `DeviceVector<RealType>` for X, Y, Z
- **Coalesced GPU Access**: Threads access contiguous memory within each array
- **No AoS Option**: SoA mandatory for GPU performance
- **Decoded from SFC**: Coordinates extracted from SFC keys via `decodeSfcToPhysical()`

## Usage Example

```cpp
// Initialize coordinate cache
CoordinateCache cache;
std::vector<Point> coordinates = domain.get_coordinates();
cache.set_coordinates(coordinates);

// Single coordinate access
Point p = cache.get_coordinate(node_id);
std::cout << "Point: (" << p.x << ", " << p.y << ", " << p.z << ")" << std::endl;

// Bulk coordinate access
std::vector<size_t> node_list = element.get_nodes();
std::vector<Point> element_coords;
cache.get_coordinates(node_list, element_coords);

// Process element coordinates
for (const auto& coord : element_coords) {
    // Compute element properties
    process_coordinate(coord);
}
```

## Caching Strategies

### Lazy Loading
- Coordinates loaded on first access
- Memory allocated only when needed
- Automatic cleanup when unused

### Lazy GPU Allocation
```cpp
// Coordinates decoded from SFC keys only when first accessed
const auto& getNodeCoordinates() {
    if (!coordinate_cache_) {
        coordinate_cache_ = std::make_unique<CoordinateCache>(...);
        coordinate_cache_->cacheNodeCoordinates(*this);
    }
    return coordinate_cache_;
}
```

### Memory Pool Allocation
- Pre-allocated memory pools
- Reduced allocation overhead
- Better memory locality

## Performance Optimizations

### Cache-Friendly Layout
```cpp
// Structure of Arrays for better cache performance
struct SOACoordinates {
    std::vector<float> x, y, z;
    
    Point get(size_t i) const {
        return {x[i], y[i], z[i]};
    }
};
```

### SIMD Operations
- Vectorized coordinate processing
- Aligned memory allocation
- Compiler auto-vectorization hints

### GPU Memory Management
```cpp
// CUDA memory allocation
float *d_x, *d_y, *d_z;
cudaMalloc(&d_x, N * sizeof(float));
cudaMalloc(&d_y, N * sizeof(float));
cudaMalloc(&d_z, N * sizeof(float));

// Copy coordinates to device
cudaMemcpy(d_x, h_x.data(), N * sizeof(float), cudaMemcpyHostToDevice);
// ... similar for y, z

// Use in CUDA kernel
__global__ void process_coordinates(float* x, float* y, float* z, size_t N) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < N) {
        // Process coordinate (x[i], y[i], z[i])
    }
}
```

## Memory Management

### Memory Usage Tracking
```cpp
size_t CoordinateCache::memory_usage() const {
    size_t total = 0;
    total += coordinates.capacity() * sizeof(Point);
    total += index_map.capacity() * sizeof(size_t);
    total += prefetch_buffer.capacity() * sizeof(Point);
    return total;
}
```

### Automatic Optimization
- Monitor access patterns
- Switch layouts based on usage
- Compress unused coordinates
- Memory defragmentation

## Integration with Other Components

### Element Processing
```cpp
// Efficient element coordinate access
void process_element(const Element& elem, const CoordinateCache& cache) {
    std::vector<Point> coords;
    cache.get_coordinates(elem.nodes, coords);
    
    // Compute element centroid
    Point centroid = {0, 0, 0};
    for (const auto& p : coords) {
        centroid.x += p.x;
        centroid.y += p.y;
        centroid.z += p.z;
    }
    centroid.x /= coords.size();
    centroid.y /= coords.size();
    centroid.z /= coords.size();
}
```

### Adjacency Computations
- Cached coordinates for distance calculations
- Spatial indexing for neighbor finding
- Coordinate transformations

## Error Handling

- Bounds checking for coordinate access
- Memory allocation failure handling
- Coordinate validation (NaN, infinity checks)
- Cache consistency verification

## See Also

- [ElementDomain Overview](ElementDomain-Overview.md)
- [GPU Acceleration](GPU-Acceleration.md)
- [Characteristic Sizes](Characteristic-Sizes.md)