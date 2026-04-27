# Unstructured Meshes in MARS

Welcome to the MARS Unstructured Meshes Documentation. This section covers the unstructured mesh handling capabilities in the MARS (Multilevel Adaptive Refinement Solver) framework.

## Overview

MARS provides comprehensive support for unstructured mesh processing, including:

- **ElementDomain**: Core class for managing mesh elements and their properties
- **Mesh Reading**: Support for various mesh formats with automatic partitioning
- **SFC Mapping**: Space-filling curve mapping for load balancing
- **Adjacency Structures**: Efficient neighbor finding and connectivity
- **Halo Management**: Ghost cell handling for parallel computations (element halo + per-node halo)
- **Coordinate Caching**: Optimized coordinate storage and access
- **Characteristic Sizes**: Mesh quality metrics and sizing functions
- **Connectivity Management**: Element-to-element relationships
- **GPU Acceleration**: CUDA-based parallel processing
- **Multi-Rank Support**: MPI-based distributed computing
- **AMR**: GPU-native multi-rank adaptive mesh refinement with solution transfer

## Key Components

### Core Infrastructure
- [ElementDomain Overview](ElementDomain-Overview.md) - Main class for unstructured mesh management
- [Mesh Reading and Partitioning](Mesh-Reading-and-Partitioning.md) - Input handling and domain decomposition
- [SFC Mapping](SFC-Mapping.md) - Load balancing through space-filling curves

### Advanced Features
- [Adjacency Structures](Adjacency-Structures.md) - Neighbor relationships and connectivity
- [Halo Management](Halo-Management.md) - Ghost cell management for parallel processing
- [Node Halo Topology](Node-Halo-Topology.md) - Per-node halo via direct CUDA-aware MPI
- [Coordinate Caching](Coordinate-Caching.md) - Efficient coordinate storage
- [Characteristic Sizes](Characteristic-Sizes.md) - Mesh quality and sizing
- [Connectivity Management](Connectivity-Management.md) - Element relationships

### Adaptive Mesh Refinement
- [AMR Module](AMR-Module.md) - Refinement pipeline (mark → refine → rebuild → transfer)
- [Solution Transfer](Solution-Transfer.md) - Field interpolation across AMR levels (warm-start CG)

### Performance & Parallelism
- [GPU Acceleration](GPU-Acceleration.md) - CUDA implementation details
- [Multi-Rank Support](Multi-Rank-Support.md) - MPI parallel processing

## Quick Start

```cpp
#include <mars.hpp>

// Hex8 mesh, double precision, uint64_t SFC keys, GPU
using Domain = mars::ElementDomain<mars::HexTag, double, uint64_t, cstone::GpuTag>;

// Read + partition + build cstone domain
Domain domain(meshFile, rank, numRanks);

// Lazy-built data: triggers HaloData (ownership) + NodeHaloTopology
const auto& d_nodeOwnership = domain.getNodeOwnershipMap();   // size = nodeCount
const auto& d_conn          = domain.getElementToNodeConnectivity();  // local node IDs

// Cache decoded node coordinates
domain.cacheNodeCoordinates();
const auto& d_x = domain.getNodeX();
const auto& d_y = domain.getNodeY();
const auto& d_z = domain.getNodeZ();

// Distributed CG halo callback
auto haloExchange = [&domain, dofMap = d_node_to_dof.data()]
    (cstone::DeviceVector<double>& p) {
        domain.exchangeNodeHalo(p, dofMap);
    };
```

For an end-to-end distributed AMR + Poisson example, see
[`examples/distributed/unstructured/mars_amr_cvfem_graph.cu`](../examples/distributed/unstructured/mars_amr_cvfem_graph.cu)
and the [AMR Module](AMR-Module.md) walkthrough.

## Architecture

The unstructured mesh system in MARS is built on Cornerstone octree and uses a lazy composition pattern where GPU components are initialized on-demand to minimize VRAM usage and startup time. Key architectural decisions include:

- **GPU-Native Design**: All data structures live in device memory (Cornerstone `DeviceVector`)
- **Template-based Design**: Type-safe mesh element handling via `AcceleratorTag`
- **Lazy Initialization**: GPU components (adjacency, halo, coordinates) built only when needed
- **SFC-Centric**: Elements identified by space-filling curve keys, not integer indices
- **Thrust Algorithms**: CSR building, sorting, and reductions use Thrust primitives
- **MPI Partitioning**: Multi-rank support via Cornerstone domain decomposition

## Contributing

When contributing to the unstructured mesh system:

1. Follow the existing template patterns for new mesh element types
2. Implement lazy initialization for new components
3. Add comprehensive tests for parallel functionality
4. Update this documentation for new features

## See Also

- [MARS Core Documentation](../README.md)
- [Examples](../../examples/)
- [API Reference](../../core/mars.hpp)