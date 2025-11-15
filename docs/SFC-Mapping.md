# SFC Mapping

Space-Filling Curve (SFC) mapping provides load balancing and spatial locality for unstructured mesh partitioning in parallel computations.

## Purpose

SFC mapping addresses the challenge of distributing mesh elements across parallel ranks while:

- Minimizing inter-partition communication
- Maintaining spatial locality
- Achieving load balance
- Supporting dynamic repartitioning

## Supported Curves

### Hilbert Curve
- Excellent spatial locality preservation
- Good load balancing properties
- Higher computational cost for mapping

### Morton (Z-Order) Curve
- Fast computation
- Good clustering properties
- Simpler implementation

## Key Components

### SFCMapping Class
```cpp
class SFCMapping {
public:
    std::vector<size_t> map_elements_to_sfc(
        const std::vector<Element>& elements,
        const std::vector<Point>& coordinates);
    
    std::vector<size_t> partition_sfc_indices(
        const std::vector<size_t>& sfc_indices, 
        size_t num_partitions);
};
```

### Element-to-SFC Mapping
- Identifies lowest SFC corner node of each element (not centroid)
- Maps 3D corner coordinates to 1D SFC key via Cornerstone encoding
- GPU kernel (`findRepresentativeNodesKernel`) processes all elements in parallel
- Handles various element types (tets, hexes, triangles, quads)

## Usage Example

```cpp
// Create SFC mapper
SFCMapping mapper;

// Get elements and coordinates from domain
auto elements = domain.get_elements();
auto coords = domain.get_coordinates();

// Map elements to SFC indices
auto sfc_indices = mapper.map_elements_to_sfc(elements, coords);

// Partition based on SFC
auto partitions = mapper.partition_sfc_indices(sfc_indices, num_ranks);

// Assign elements to local rank
std::vector<Element> local_elements;
for (size_t i = 0; i < elements.size(); ++i) {
    if (partitions[i] == rank) {
        local_elements.push_back(elements[i]);
    }
}
```

## Algorithm Details

### SFC Encoding (GPU-based)
1. **Corner Identification**: Find lowest SFC corner per element (GPU kernel)
2. **Coordinate Decoding**: Extract physical coordinates from SFC keys
3. **Key Generation**: Cornerstone `sfcKey()` encoding on GPU
4. **Sorting**: Thrust `sort_by_key()` orders elements by SFC index on device

### Load Balancing
- **Equal Partitioning**: Divide SFC range equally among ranks
- **Weight Consideration**: Account for element computational cost
- **Imbalance Detection**: Monitor and correct load imbalances

## Performance Characteristics

- **Mapping Time**: O(N log N) for N elements
- **Memory Usage**: O(N) additional storage for indices
- **Communication**: Minimal inter-rank data transfer
- **Scalability**: Good parallel scaling up to thousands of ranks

## Integration with ElementDomain

SFC mapping is automatically used during domain construction:

```cpp
// SFC mapping happens during domain initialization
ElementDomain<TetTag, float, unsigned> domain(mesh_path, rank, num_ranks);

// The SFC map is available for queries
auto sfc_map = domain.get_sfc_mapping();
auto local_sfc_range = sfc_map.get_local_range(rank);
```

## Integration with Cornerstone

- **SFC Key Encoding**: Uses Cornerstone `sfcKey<KeyType>(x, y, z, box)` for encoding
- **Box Management**: Global bounding box from Cornerstone domain
- **Decoding**: `decodeSfcToPhysical()` converts keys back to coordinates (host/device compatible)
- **Local-to-Global Mapping**: Sparse mapping via device vectors for distributed ownership

## See Also

- [ElementDomain Overview](ElementDomain-Overview.md)
- [Mesh Reading and Partitioning](Mesh-Reading-and-Partitioning.md)
- [Multi-Rank Support](Multi-Rank-Support.md)