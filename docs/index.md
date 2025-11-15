# Unstructured Meshes in MARS

Welcome to the MARS Unstructured Meshes Documentation. This section covers the unstructured mesh handling capabilities in the MARS (Multilevel Adaptive Refinement Solver) framework.

## Overview

MARS provides comprehensive support for unstructured mesh processing, including:

- **ElementDomain**: Core class for managing mesh elements and their properties
- **Mesh Reading**: Support for various mesh formats with automatic partitioning
- **SFC Mapping**: Space-filling curve mapping for load balancing
- **Adjacency Structures**: Efficient neighbor finding and connectivity
- **Halo Management**: Ghost cell handling for parallel computations
- **Coordinate Caching**: Optimized coordinate storage and access
- **Characteristic Sizes**: Mesh quality metrics and sizing functions
- **Connectivity Management**: Element-to-element relationships
- **GPU Acceleration**: CUDA-based parallel processing
- **Multi-Rank Support**: MPI-based distributed computing

## Key Components

### Core Infrastructure
- [ElementDomain Overview](ElementDomain-Overview.md) - Main class for unstructured mesh management
- [Mesh Reading and Partitioning](Mesh-Reading-and-Partitioning.md) - Input handling and domain decomposition
- [SFC Mapping](SFC-Mapping.md) - Load balancing through space-filling curves

### Advanced Features
- [Adjacency Structures](Adjacency-Structures.md) - Neighbor relationships and connectivity
- [Halo Management](Halo-Management.md) - Ghost cell management for parallel processing
- [Coordinate Caching](Coordinate-Caching.md) - Efficient coordinate storage
- [Characteristic Sizes](Characteristic-Sizes.md) - Mesh quality and sizing
- [Connectivity Management](Connectivity-Management.md) - Element relationships

### Performance & Parallelism
- [GPU Acceleration](GPU-Acceleration.md) - CUDA implementation details
- [Multi-Rank Support](Multi-Rank-Support.md) - MPI parallel processing

## Quick Start

```cpp
#include <mars.hpp>

// Create an ElementDomain for tetrahedral meshes
ElementDomain<TetTag, float, unsigned> domain("mesh_dir", rank, numRanks);

// Access mesh elements
auto elements = domain.get_elements();
auto coordinates = domain.get_coordinates();

// Build adjacency structures
domain.build_adjacency();

// Get characteristic sizes
auto sizes = domain.get_characteristic_sizes();
```

## Architecture

The unstructured mesh system in MARS uses a lazy composition pattern where components are initialized on-demand to minimize memory usage and startup time. Key architectural decisions include:

- **Template-based Design**: Type-safe mesh element handling
- **Lazy Initialization**: Components built only when needed
- **Memory Efficiency**: Shared ownership and caching
- **Parallel Support**: MPI and CUDA integration
- **Extensibility**: Plugin architecture for new mesh types

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