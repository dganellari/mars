# Halo Management

Halo management handles ghost cells and boundary data exchange in parallel unstructured mesh computations, ensuring correct handling of inter-partition boundaries.

## Purpose

In parallel computing, halo (ghost) elements are needed for:

- Finite element computations across partition boundaries
- Maintaining solution continuity
- Parallel data exchange
- Load balancing verification

## Key Components

### HaloData Struct (GPU)
```cpp
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
struct HaloData {
    using DeviceVector = typename VectorSelector<uint8_t, AcceleratorTag>::type;
    
    DeviceVector d_nodeOwnership_;          // Node ownership flags (0=local, 1=halo)
    DeviceVector d_haloElementIndices_;     // Indices of halo elements
    
    void buildNodeOwnership(const ElementDomain& domain, int rank, int numRanks);
    void buildHaloElementIndices(const ElementDomain& domain);
};
```

### Halo Exchange Protocol
- Identifies elements needed across partition boundaries
- Manages send/receive buffers
- Handles non-blocking MPI communications

## Usage Example

```cpp
// Build halo data after adjacency construction
domain.build_halo_data();

// Get halo information
auto halo = domain.get_halo_data();

// Check halo sizes
std::cout << "Local halo elements: " << halo.local_halo_elements.size() << std::endl;
std::cout << "Remote halo elements: " << halo.remote_halo_elements.size() << std::endl;

// Exchange solution data
std::vector<double> solution = domain.get_solution_data();
halo.exchange_data(solution);

// Now solution includes ghost values
std::cout << "Solution size with halo: " << solution.size() << std::endl;
```

## Halo Construction Algorithm

### Step 1: Identify Boundary Elements
```cpp
// Find elements adjacent to partition boundaries
for (size_t elem_id = 0; elem_id < domain.size(); ++elem_id) {
    auto neighbors = adjacency.get_neighbors(elem_id);
    for (auto neighbor : neighbors) {
        if (domain.get_owner_rank(neighbor) != rank) {
            // This element needs halo data
            boundary_elements.insert(elem_id);
            break;
        }
    }
}
```

### Step 2: Determine Halo Layers
- **Single Layer**: Only immediately adjacent elements
- **Multi-Layer**: Extended halo for higher-order methods
- **Adaptive**: Halo size based on stencil requirements

### Step 3: Communication Setup
- Build send/receive rank lists
- Allocate communication buffers
- Set up MPI data types

## Halo Element Identification

### GPU-based Node Ownership
```cpp
// Mark nodes as local (0) or halo (1) based on element ownership
// Uses Cornerstone domain to determine rank boundaries
__global__ void markLocalElementNodesKernel(
    uint8_t* nodeOwnership,
    const KeyType* sfc_conn,
    const KeyType* localToGlobalSfcMap,
    size_t localElementCount);

// Extract halo element indices from ownership flags
void buildHaloElementIndices(DeviceVector& d_haloIndices,
                            const DeviceVector& d_nodeOwnership,
                            size_t elementCount);
```

**Note**: Actual MPI data exchange is handled by Cornerstone domain, not directly by MARS.

### Asynchronous Exchange
- Overlapping computation with communication
- Non-blocking operations for performance
- Careful synchronization to avoid race conditions

## Performance Considerations

### Memory Overhead
- Halo elements typically 5-15% of total elements
- Additional storage for communication buffers
- Trade-off between memory and communication frequency

### Communication Optimization
- Minimize message count through aggregation
- Use derived MPI data types for complex data
- Overlap communication with computation

### Load Balancing Impact
- Halo size affects parallel efficiency
- Monitor communication vs. computation ratio
- Adjust partitioning to minimize halo overhead

## Integration with GPU Acceleration

For GPU-accelerated computations:

```cpp
// Transfer halo data to/from GPU
halo.exchange_gpu_data(device_solution, host_buffer);

// GPU computation can proceed with halo data
kernel_compute<<<blocks, threads>>>(device_solution, halo_size);
```

## Error Handling

- Detect communication failures
- Validate data consistency across partitions
- Handle load imbalances gracefully
- Provide debugging information for halo issues

## See Also

- [ElementDomain Overview](ElementDomain-Overview.md)
- [Adjacency Structures](Adjacency-Structures.md)
- [Multi-Rank Support](Multi-Rank-Support.md)