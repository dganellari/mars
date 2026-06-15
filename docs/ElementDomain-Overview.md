# ElementDomain Overview

The `ElementDomain` class is the central component of MARS's unstructured mesh system, providing a unified interface for managing mesh elements, their properties, and relationships.

## Purpose

`ElementDomain` serves as the main container and manager for unstructured mesh data, handling:

- Element storage and access
- Coordinate management
- Partitioning information
- Lazy initialization of dependent components
- Type-safe element operations

## Template Parameters

```cpp
template<typename ElementTag     = TetTag,
         typename RealType       = float,
         typename KeyType        = unsigned,
         typename AcceleratorTag = cstone::GpuTag>
class ElementDomain;
```

- **ElementTag**: Element type tag (e.g., `TetTag` for tetrahedra, `HexTag` for hexahedra)
- **RealType**: Floating-point type for coordinates (e.g., `float`, `double`)
- **KeyType**: Unsigned integer type for SFC keys (e.g., `unsigned`, `uint64_t`)
- **AcceleratorTag**: Backend tag (`cstone::GpuTag` for CUDA/HIP, `cstone::CpuTag` for host)

## Key Features

### Element Management
- Stores mesh elements as SFC keys in device memory (`DeviceVector<KeyType>`)
- Elements represented by lowest SFC corner node, not centroids
- Connectivity stored as tuple of device vectors via `VectorSelector<T, GpuTag>`

### Coordinate Handling
- Vertex coordinates stored as Structure-of-Arrays (SoA) on GPU
- Lazy coordinate caching in device memory for performance
- No host mirroring - data stays on GPU after initialization

### Partitioning Support
- Domain decomposition via Cornerstone SFC partitioning
- Local-to-global SFC map for sparse element ownership
- MPI operations only for metadata sync, data stays on GPU

## Usage Example

```cpp
// Create domain for tetrahedral mesh (read + partition + cstone sync on construction)
ElementDomain<TetTag, double, uint64_t, cstone::GpuTag> domain("mesh_directory", rank, numRanks);

// Mesh sizes this rank sees (owned + halo)
std::cout << "Elements: " << domain.getElementCount() << std::endl;
std::cout << "Nodes:    " << domain.getNodeCount()    << std::endl;

// Cache decoded node coordinates (SoA, device-side), then access them
domain.cacheNodeCoordinates();
const auto& d_x = domain.getNodeX();   // DeviceVector<RealType>
const auto& d_y = domain.getNodeY();
const auto& d_z = domain.getNodeZ();

// Element-to-node connectivity (local node IDs per element), built lazily
const auto& d_conn = domain.getElementToNodeConnectivity();

// Per-element corner SFC keys, e.g. for element 0 (host accessor):
auto keys = domain.getConnectivity(0);  // tuple of NodesPerElement SFC keys
```

## Lazy Composition

The `ElementDomain` uses lazy initialization for GPU memory efficiency:

- GPU components (`AdjacencyData`, `HaloData`, `CoordinateCache`) allocated on first access
- Reduces VRAM footprint and startup time
- Thread-safe initialization with mutex protection
- Uses `std::unique_ptr` for optional GPU allocations

## Dependencies

- **Mesh Reading**: For initial mesh loading
- **SFC Mapping**: For load balancing
- **Adjacency Structures**: For neighbor relationships
- **Halo Management**: For parallel ghost cells

## Implementation Notes

- Uses shared pointers for component ownership
- Implements copy-on-write semantics where appropriate
- Provides const-correct access methods
- Supports both host and device (GPU) operations

## See Also

- [Mesh Reading and Partitioning](Mesh-Reading-and-Partitioning.md)
- [SFC Mapping](SFC-Mapping.md)
- [Adjacency Structures](Adjacency-Structures.md)