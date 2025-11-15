# Multi-Rank Support

Multi-rank support provides MPI-based distributed computing capabilities for unstructured mesh processing, enabling parallel execution across multiple compute nodes.

## Purpose

Multi-rank support enables:

- Large-scale parallel computations
- Distributed mesh processing
- Load-balanced parallel execution
- Inter-node communication optimization
- Scalable algorithms for massive meshes

## Key Components

### Cornerstone Domain Integration
```cpp
// MARS integrates with Cornerstone for MPI coordination
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
class ElementDomain {
private:
    cstone::Domain<KeyType, RealType, AcceleratorTag> cstone_domain_;
    int rank_, numRanks_;
    
    // Local-to-global SFC mapping for distributed elements
    DeviceVector d_localToGlobalSfcMap_;  // Sparse ownership
    
    // Sync with Cornerstone after mesh loading
    void sync() {
        // Build SFC keys and integrate with Cornerstone domain
        // All MPI coordination handled by Cornerstone
    }
};
```

### MPI Integration
- Communicator management
- Collective operations
- Point-to-point communication
- Custom data types for mesh structures

## Usage Example

```cpp
// Initialize MPI
MPI_Init(&argc, &argv);

// Create multi-rank manager
MultiRankManager mpi_mgr;

// Create domain with MPI support
ElementDomain<TetTag, float, unsigned> domain("mesh_dir", 
                                             mpi_mgr.get_rank(), 
                                             mpi_mgr.get_num_ranks());

// Domain automatically handles partitioning
std::cout << "Rank " << mpi_mgr.get_rank() 
          << " has " << domain.get_elements().size() << " elements" << std::endl;

// Perform parallel computation
domain.build_adjacency();
domain.compute_characteristic_sizes();

// Global reduction for verification
double local_sum = std::accumulate(domain.get_characteristic_sizes().begin(),
                                  domain.get_characteristic_sizes().end(), 0.0);
double global_sum;
mpi_mgr.allreduce(&local_sum, &global_sum, 1, MPI_SUM);

if (mpi_mgr.get_rank() == 0) {
    std::cout << "Global characteristic size sum: " << global_sum << std::endl;
}

MPI_Finalize();
```

## Domain Decomposition

### Per-Rank Binary File Reading
```cpp
// Each rank reads its own pre-partitioned binary file
auto meshData = readMeshDataSoA<float, unsigned>(
    meshPath, rank, numRanks,
    [](const std::vector<float>& coords_in) {
        // Convert float to RealType if needed
        std::vector<RealType> coords_out(coords_in.begin(), coords_in.end());
        return coords_out;
    });

// Transfer to GPU immediately
d_x_.assign(meshData.x.begin(), meshData.x.end());
d_y_.assign(meshData.y.begin(), meshData.y.end());
d_z_.assign(meshData.z.begin(), meshData.z.end());
// ... connectivity also copied to device

// Generate SFC keys on GPU
generateSfcKeys(d_x_.data(), d_y_.data(), d_z_.data(),
               d_localToGlobalSfcMap_.data(), elementCount, box);
```

### Load Balancing
```cpp
// Dynamic load balancing
void rebalance_load(const std::vector<double>& computational_costs,
                   std::vector<std::vector<size_t>>& rank_elements,
                   MultiRankManager& mpi_mgr) {
    
    // Compute current load imbalance
    double local_load = std::accumulate(computational_costs.begin(),
                                       computational_costs.end(), 0.0);
    double max_load, avg_load;
    mpi_mgr.allreduce(&local_load, &max_load, 1, MPI_MAX);
    mpi_mgr.allreduce(&local_load, &avg_load, 1, MPI_SUM);
    avg_load /= mpi_mgr.get_num_ranks();
    
    double imbalance = max_load / avg_load;
    
    if (imbalance > imbalance_threshold) {
        // Trigger load rebalancing
        redistribute_elements(rank_elements, computational_costs, mpi_mgr);
    }
}
```

## Communication Patterns

### GPU-Resident Halo Management
```cpp
// Halo elements identified on GPU, data never leaves device
struct HaloData {
    DeviceVector<uint8_t> d_nodeOwnership_;      // 0=local, 1=halo
    DeviceVector<size_t> d_haloElementIndices_;  // Halo element IDs
    
    void buildNodeOwnership(const ElementDomain& domain, int rank, int numRanks) {
        // Mark nodes based on Cornerstone domain boundaries
        // All operations on GPU
        markLocalElementNodesKernel<<<...>>>(
            d_nodeOwnership_.data(),
            domain.d_conn_sfc_keys_,
            domain.d_localToGlobalSfcMap_.data(),
            domain.localElementCount());
    }
};

// Note: Actual MPI exchange handled by Cornerstone, not exposed in MARS API
```

### Collective Operations
```cpp
// Global statistics computation
void compute_global_statistics(const ElementDomain& domain,
                              MultiRankManager& mpi_mgr) {
    
    // Local statistics
    size_t local_elements = domain.get_elements().size();
    double local_volume = compute_local_volume(domain);
    
    // Global sums
    size_t global_elements;
    double global_volume;
    
    mpi_mgr.allreduce(&local_elements, &global_elements, 1, MPI_SUM);
    mpi_mgr.allreduce(&local_volume, &global_volume, 1, MPI_SUM);
    
    // Global min/max
    double local_min_size = *std::min_element(domain.get_characteristic_sizes().begin(),
                                             domain.get_characteristic_sizes().end());
    double local_max_size = *std::max_element(domain.get_characteristic_sizes().begin(),
                                             domain.get_characteristic_sizes().end());
    
    double global_min_size, global_max_size;
    mpi_mgr.allreduce(&local_min_size, &global_min_size, 1, MPI_MIN);
    mpi_mgr.allreduce(&local_max_size, &global_max_size, 1, MPI_MAX);
    
    if (mpi_mgr.get_rank() == 0) {
        std::cout << "Global mesh statistics:" << std::endl;
        std::cout << "  Total elements: " << global_elements << std::endl;
        std::cout << "  Total volume: " << global_volume << std::endl;
        std::cout << "  Size range: [" << global_min_size << ", " << global_max_size << "]" << std::endl;
    }
}
```

## Performance Optimization

### Communication Overlap
```cpp
// Overlap communication with computation
void overlapped_processing(ElementDomain& domain, MultiRankManager& mpi_mgr) {
    // Start halo exchange asynchronously
    auto halo_future = std::async(std::launch::async, [&]() {
        domain.exchange_halo_data();
    });
    
    // Perform local computations while communicating
    domain.compute_local_properties();
    
    // Wait for communication to complete
    halo_future.wait();
    
    // Use halo data for boundary computations
    domain.compute_boundary_properties();
}
```

### Memory Management
```cpp
// Distributed memory allocation
class DistributedAllocator {
public:
    void* allocate(size_t size, MultiRankManager& mpi_mgr) {
        // Allocate memory across ranks
        size_t local_size = size / mpi_mgr.get_num_ranks();
        size_t remainder = size % mpi_mgr.get_num_ranks();
        
        if (mpi_mgr.get_rank() < remainder) {
            local_size++;
        }
        
        return malloc(local_size);
    }
};
```

## Scalability Features

### Weak Scaling
- Problem size increases with processor count
- Constant work per processor
- Communication scales appropriately

### Strong Scaling
- Fixed problem size, increasing processors
- Parallel efficiency measurement
- Communication overhead analysis

### Load Balancing Metrics
```cpp
void analyze_load_balance(const std::vector<double>& local_times,
                         MultiRankManager& mpi_mgr) {
    
    double local_time = local_times[mpi_mgr.get_rank()];
    double max_time, avg_time;
    
    mpi_mgr.allreduce(&local_time, &max_time, 1, MPI_MAX);
    mpi_mgr.allreduce(&local_time, &avg_time, 1, MPI_SUM);
    avg_time /= mpi_mgr.get_num_ranks();
    
    double efficiency = avg_time / max_time;
    double imbalance = (max_time - avg_time) / avg_time;
    
    if (mpi_mgr.get_rank() == 0) {
        std::cout << "Load balance analysis:" << std::endl;
        std::cout << "  Parallel efficiency: " << efficiency * 100 << "%" << std::endl;
        std::cout << "  Load imbalance: " << imbalance * 100 << "%" << std::endl;
    }
}
```

## Error Handling

### MPI Error Checking
```cpp
#define MPI_CHECK(call) \
    do { \
        int error = call; \
        if (error != MPI_SUCCESS) { \
            char error_string[MPI_MAX_ERROR_STRING]; \
            int length; \
            MPI_Error_string(error, error_string, &length); \
            std::cerr << "MPI error: " << error_string << std::endl; \
            MPI_Abort(MPI_COMM_WORLD, error); \
        } \
    } while(0)

MPI_CHECK(MPI_Allreduce(sendbuf, recvbuf, count, MPI_DOUBLE, MPI_SUM, comm));
```

### Deadlock Prevention
- Non-blocking communication patterns
- Proper request management
- Timeout mechanisms for debugging

## Integration with Other Components

- **Halo Management**: MPI-based data exchange
- **SFC Mapping**: Distributed load balancing
- **GPU Acceleration**: Multi-GPU, multi-node configurations

## See Also

- [ElementDomain Overview](ElementDomain-Overview.md)
- [Halo Management](Halo-Management.md)
- [SFC Mapping](SFC-Mapping.md)