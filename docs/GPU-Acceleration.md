# GPU Acceleration

GPU acceleration provides CUDA-based parallel processing for unstructured mesh operations, enabling high-performance computing on NVIDIA GPUs.

## Purpose

**GPU execution is mandatory** in MARS unstructured mesh system. All data structures and algorithms are GPU-native:

- Mesh storage in device memory (`DeviceVector` via `VectorSelector<T, GpuTag>`)
- SFC key generation and sorting (Thrust algorithms)
- Adjacency/halo building (GPU kernels + Thrust primitives)
- Characteristic size computation (parallel GPU kernels)
- No CPU fallback - `AcceleratorTag = GpuTag` required

## Key Components

### VectorSelector Pattern
```cpp
// Template-based device vector selection
template<typename T, typename AcceleratorTag>
struct VectorSelector;

// Specialization for GPU (Cornerstone DeviceVector)
template<typename T>
struct VectorSelector<T, cstone::GpuTag> {
    using type = cstone::DeviceVector<T>;
};

// Usage in ElementDomain
using DeviceVector = typename VectorSelector<KeyType, AcceleratorTag>::type;
DeviceVector d_localToGlobalSfcMap_;  // Lives in GPU memory
DeviceVector d_conn_local_ids_;       // All connectivity on GPU
```

### CUDA Kernel Organization
- Element-wise kernels for local computations
- Reduction kernels for global operations
- Memory coalescing optimizations
- Shared memory utilization

## Usage Example

```cpp
// Initialize GPU accelerator
GPUAccelerator gpu;
gpu.initialize(0);  // Use GPU device 0

// Transfer mesh data to GPU
auto elements = domain.get_elements();
auto coords = domain.get_coordinates();

Element* d_elements = gpu.allocate_device<Element>(elements.size());
float* d_coords = gpu.allocate_device<float>(coords.size() * 3);

gpu.copy_to_device(elements.data(), d_elements, elements.size());
gpu.copy_to_device(coords.data(), d_coords, coords.size() * 3);

// Launch element processing kernel
gpu.launch_element_processing(d_elements, elements.size());

// Transfer results back
std::vector<double> results(elements.size());
double* d_results = gpu.allocate_device<double>(results.size());
gpu.copy_to_host(d_results, results.data(), results.size());
```

## CUDA Kernel Examples

### SFC Representative Node Kernel
```cpp
template<typename ElementTag, typename KeyType, typename RealType>
__global__ void findRepresentativeNodesKernel(
    const KeyType* indices0,  // First node SFC keys
    const KeyType* indices1,  // Second node SFC keys
    // ... other node indices
    KeyType* repNodes,        // Output: lowest SFC key per element
    size_t numElements) {
    
    size_t elemIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (elemIdx >= numElements) return;
    
    // Find minimum SFC key among element nodes
    KeyType minKey = indices0[elemIdx];
    minKey = min(minKey, indices1[elemIdx]);
    // ... check all nodes
    
    repNodes[elemIdx] = minKey;  // Lowest corner, not centroid
}
```

### Adjacency Computation Kernel
```cpp
__global__ void build_adjacency_kernel(const Element* elements,
                                     size_t num_elements,
                                     int* adjacency_matrix,
                                     size_t pitch) {
    size_t elem_i = blockIdx.x * blockDim.x + threadIdx.x;
    size_t elem_j = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (elem_i < num_elements && elem_j < num_elements && elem_i != elem_j) {
        // Check if elements share a face
        bool adjacent = check_face_adjacency(elements[elem_i], elements[elem_j]);
        
        if (adjacent) {
            // Atomic update for thread safety
            atomicAdd(&adjacency_matrix[elem_i * pitch + elem_j / 32], 
                     1 << (elem_j % 32));
        }
    }
}
```

## Memory Management

### Unified Memory
```cpp
// CUDA Unified Memory for simplified data management
cudaMallocManaged(&unified_elements, num_elements * sizeof(Element));
cudaMallocManaged(&unified_coords, num_coords * sizeof(float) * 3);

// Data accessible from both CPU and GPU
initialize_mesh_data(unified_elements, unified_coords);

// GPU computation
process_elements_kernel<<<blocks, threads>>>(unified_elements, unified_coords);

// No explicit data transfer needed
```

### Memory Pool Allocation
```cpp
class GPUMemoryPool {
private:
    std::vector<void*> allocations;
    std::vector<size_t> sizes;
    
public:
    template<typename T>
    T* allocate(size_t count) {
        size_t bytes = count * sizeof(T);
        T* ptr;
        cudaMalloc(&ptr, bytes);
        allocations.push_back(ptr);
        sizes.push_back(bytes);
        return ptr;
    }
    
    void deallocate_all() {
        for (auto ptr : allocations) {
            cudaFree(ptr);
        }
        allocations.clear();
        sizes.clear();
    }
};
```

## Performance Optimizations

### Kernel Launch Configuration
```cpp
// Optimal thread block size
const int BLOCK_SIZE = 256;

// Calculate grid dimensions
dim3 blocks((num_elements + BLOCK_SIZE - 1) / BLOCK_SIZE);
dim3 threads(BLOCK_SIZE);

// Launch kernel with optimal configuration
process_kernel<<<blocks, threads>>>(d_data, num_elements);
```

### Memory Coalescing
```cpp
// Structure of Arrays (SoA) for better memory access
struct SOACoordinates {
    float* x;
    float* y; 
    float* z;
};

// Coalesced memory access in kernel
__global__ void process_coords_kernel(SOACoordinates coords, size_t N) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < N) {
        float val = coords.x[i] + coords.y[i] + coords.z[i];
        // Coalesced access pattern
    }
}
```

### Shared Memory Usage
```cpp
__global__ void element_kernel(const Element* elements, 
                              const float* coords,
                              float* results) {
    __shared__ float shared_coords[32 * 3];  // 32 elements * 3 coords
    
    size_t local_id = threadIdx.x;
    size_t elem_id = blockIdx.x * blockDim.x + local_id;
    
    // Load coordinates into shared memory
    if (elem_id < num_elements) {
        const Element& elem = elements[elem_id];
        for (int i = 0; i < elem.num_nodes; ++i) {
            size_t node_idx = elem.nodes[i];
            shared_coords[local_id * 3 + i * 3] = coords[node_idx * 3];
            shared_coords[local_id * 3 + i * 3 + 1] = coords[node_idx * 3 + 1];
            shared_coords[local_id * 3 + i * 3 + 2] = coords[node_idx * 3 + 2];
        }
    }
    
    __syncthreads();
    
    // Process using shared memory
    if (elem_id < num_elements) {
        results[elem_id] = process_shared_coords(&shared_coords[local_id * 3]);
    }
}
```

## Integration with CPU Components

### Hybrid CPU-GPU Processing
```cpp
void hybrid_mesh_processing(const ElementDomain& domain) {
    // CPU preprocessing
    auto elements = domain.get_elements();
    auto adjacency = domain.get_adjacency_data();
    
    // GPU acceleration for compute-intensive parts
    GPUAccelerator gpu;
    gpu.process_elements(elements);
    gpu.build_adjacency_gpu(elements);
    
    // CPU postprocessing
    domain.update_from_gpu_results(gpu.get_results());
}
```

### Asynchronous Operations
```cpp
// Overlap computation with data transfer
cudaStream_t stream;
cudaStreamCreate(&stream);

// Asynchronous data transfer
cudaMemcpyAsync(d_elements, h_elements, size, cudaMemcpyHostToDevice, stream);

// Asynchronous kernel execution
kernel<<<blocks, threads, 0, stream>>>(d_elements);

// Asynchronous result transfer
cudaMemcpyAsync(h_results, d_results, size, cudaMemcpyDeviceToHost, stream);

// CPU work while GPU is busy
do_cpu_work();

// Wait for GPU completion
cudaStreamSynchronize(stream);
```

## Error Handling and Debugging

### CUDA Error Checking
```cpp
#define CUDA_CHECK(call) \
    do { \
        cudaError_t error = call; \
        if (error != cudaSuccess) { \
            std::cerr << "CUDA error: " << cudaGetErrorString(error) \
                      << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
            exit(1); \
        } \
    } while(0)

CUDA_CHECK(cudaMalloc(&d_data, size));
CUDA_CHECK(cudaMemcpy(d_data, h_data, size, cudaMemcpyHostToDevice));
```

### GPU Memory Debugging
```cpp
// Enable CUDA memory checking
cudaSetDeviceFlags(cudaDeviceMapHost);

// Use cuda-memcheck for memory errors
// cuda-memcheck ./program

// Check for kernel launch errors
cudaError_t error = cudaGetLastError();
if (error != cudaSuccess) {
    std::cerr << "Kernel launch error: " << cudaGetErrorString(error) << std::endl;
}
```

## Performance Benchmarking

### Profiling Tools
- **Nsight Systems**: Timeline profiling
- **Nsight Compute**: Kernel analysis
- **NVProf**: Command-line profiling

### Performance Metrics
- Memory throughput (GB/s)
- Compute utilization (%)
- Memory utilization (%)
- PCIe transfer bandwidth

## See Also

- [ElementDomain Overview](ElementDomain-Overview.md)
- [Coordinate Caching](Coordinate-Caching.md)
- [Multi-Rank Support](Multi-Rank-Support.md)