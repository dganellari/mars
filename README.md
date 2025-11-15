[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause) [![Build status](https://ci.appveyor.com/api/projects/status/a6kjacwk5e5pd4by/branch/development?svg=true)](https://ci.appveyor.com/project/zulianp/mars/branch/development) [![Documentation](https://readthedocs.org/projects/mesh-adaptive-refinement-for-supercomputing-mars/badge/?version=latest)](https://mesh-adaptive-refinement-for-supercomputing-mars.readthedocs.io/en/latest/?badge=latest)


# M.A.R.S #
## Mesh Adaptive Refinement for Supercomputing ##

**[Read the Full Documentation](https://mesh-adaptive-refinement-for-supercomputing-mars.readthedocs.io/en/latest/)**

MARS is an open-source mesh management library designed to handle N-dimensional elements (N <= 4). 
MARS is developed in C++ and makes use of template meta-programming to have compile time dimensions of elements and vectors, thus allowing for both compile time performance optimizations and concise and reusable code.

The main features of MARS consist of:

1. Parallel mesh generation

2. Adaptive mesh refinement using bisection algorithms

3. Conforming mesh data-structure

4. Mesh quality estimators to study the output of different mesh-refinement strategies

5. Performance portable algorithms and data-structures targetting different accelerators

6. Performance portable space filling curves algorithms for efficient mesh management.

7. Unstructured mesh support through cornerstone library, enabling octree-based adaptive mesh refinement for complex geometries.

MARS targets multi-core CPUs and GPUs using the C++ Kokkos programming model. The mesh is entirely constructed and stored on the device (GPUs). This enables libraries using MARS to perform further operations directly on the device, avoiding going through the host. 

Currently, MARS supports as its performance portable, parallel, adaptive refinement based algorithm the LEPP (Longest edge propagation path) from Rivara. Mesh generation is fully supported in parallel.

Performance portable forest of octrees and space filling curves algorithms for adaptive mesh refinement are being planned.

## Downloading MARS and its dependencies ##

Clone the repository and its submodules. MARS relies on googletest and google/benchmark.

`git clone --recurse-submodules https://bitbucket.org/zulianp/mars.git`
or for older git versions

`git clone https://bitbucket.org/zulianp/mars.git && cd mars && git submodule update --init --recursive`

Compiling M.A.R.S for serial usage:

	- cd mars/
	- mkdir build
	- cd build
	- cmake ..
	- make

## MARS Kokkos requirements ##

Mars depends on both Kokkos and Kokkos Kernels libraries. 

It will automatically find Kokkos if installed into your system. It can work with kokkos standalone or with kokkos from the Trilinos library. 

Mars looks for KOKKOS_DIR or TRILINOS_DIR into the environment variables. 
When using Trilinos it will find them from Trilinos in $TRILINOS_DIR otherwise it will look for kokkos and kokkos kernels installations at $KOKKOS_DIR.

Use -DMARS_ENABLE_KOKKOS=ON to use the feature. For more details check CMakeLists.txt.

The default when compiling MARS with Kokkos without specifing any other CMAKE flag is the Kokkos/OpenMP execution space. Kokkos should also be compiled with OpenMP support. Otherwise the default is the serial execution space.

To compile for CUDA the Cmake flag needs to be set: MARS_ENABLE_CUDA=ON. An example would be: 
```
cmake -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_BUILD_TYPE=Release -DMARS_ENABLE_KOKKOS=ON -DMARS_ENABLE_CUDA=ON ..
```

If compiled for CUDA then Kokkos should also be compiled with CUDA (Kokkos_ENABLE_CUDA=ON) and CUDA_LAMBDA (Kokkos_ENABLE_CUDA_LAMBDA=ON) support.

## Unstructured Mesh Support ##

MARS supports GPU-native unstructured meshes through integration with the Cornerstone library, enabling space-filling curve (SFC) based mesh management for complex geometries and distributed simulations.

### Key Features

- **GPU-Native Architecture**: All data structures live in device memory (`DeviceVector` via Cornerstone)
- **SFC-Based Partitioning**: Elements identified by space-filling curve keys for optimal load balancing
- **Lazy Composition**: Components (adjacency, halo, coordinates) allocated on-demand to minimize VRAM usage
- **Thrust Algorithms**: CSR building, sorting, and reductions use GPU-optimized Thrust primitives
- **MPI Integration**: Multi-rank support via Cornerstone domain decomposition
- **Element Support**: Tetrahedra, hexahedra, triangles, and quadrilaterals

### Quick Start

```cpp
#include "domain.hpp"

// Create GPU-native unstructured domain
ElementDomain<TetTag, float, unsigned> domain("mesh_dir", rank, numRanks);

// Components built lazily on first access
auto offsets = domain.getNodeToElementOffsets();  // Builds adjacency
auto coords = domain.getNodeCoordinates();        // Caches coordinates
auto sizes = domain.getCharacteristicSizes();     // Computes sizes
```

### Build Configuration

To enable unstructured support:
- Set `-DMARS_ENABLE_UNSTRUCTURED=ON` during CMake configuration
- Cornerstone is fetched automatically if not found on the system
- GPU support requires `-DMARS_ENABLE_CUDA=ON` or `-DMARS_ENABLE_HIP=ON`

Example CMake command for unstructured with CUDA:

```bash
cmake .. \
  -DMARS_ENABLE_KOKKOS=OFF \
  -DMARS_ENABLE_CUDA=ON \
  -DMARS_ENABLE_TESTS=ON \
  -DMARS_ENABLE_UNSTRUCTURED=ON \
  -DCMAKE_CUDA_ARCHITECTURES=90
```

### Documentation

For comprehensive guides and API references, see:

- **[Unstructured Meshes Documentation](https://mesh-adaptive-refinement-for-supercomputing-mars.readthedocs.io/en/latest/)**
- **[ElementDomain Overview](https://mesh-adaptive-refinement-for-supercomputing-mars.readthedocs.io/en/latest/ElementDomain-Overview/)** - Core mesh management class
- **[Mesh Reading & Partitioning](https://mesh-adaptive-refinement-for-supercomputing-mars.readthedocs.io/en/latest/Mesh-Reading-and-Partitioning/)** - Binary mesh format and loading
- **[SFC Mapping](https://mesh-adaptive-refinement-for-supercomputing-mars.readthedocs.io/en/latest/SFC-Mapping/)** - Space-filling curve based load balancing
- **[Adjacency Structures](https://mesh-adaptive-refinement-for-supercomputing-mars.readthedocs.io/en/latest/Adjacency-Structures/)** - CSR-based neighbor finding
- **[Halo Management](https://mesh-adaptive-refinement-for-supercomputing-mars.readthedocs.io/en/latest/Halo-Management/)** - Ghost element handling
- **[Coordinate Caching](https://mesh-adaptive-refinement-for-supercomputing-mars.readthedocs.io/en/latest/Coordinate-Caching/)** - GPU SoA coordinate storage
- **[Characteristic Sizes](https://mesh-adaptive-refinement-for-supercomputing-mars.readthedocs.io/en/latest/Characteristic-Sizes/)** - Mesh quality metrics
- **[GPU Acceleration](https://mesh-adaptive-refinement-for-supercomputing-mars.readthedocs.io/en/latest/GPU-Acceleration/)** - CUDA kernel implementation
- **[Multi-Rank Support](https://mesh-adaptive-refinement-for-supercomputing-mars.readthedocs.io/en/latest/Multi-Rank-Support/)** - MPI distributed computing

### Implementation Details

The unstructured backend uses:
- **Lazy initialization** for memory efficiency (adjacency, halo, coordinates built on-demand)
- **SFC keys as connectivity** for sparse global element identification
- **Thrust-based CSR** building via `sort_by_key`, `reduce_by_key`, `exclusive_scan`
- **Lowest SFC corner** representation (not centroids) for element identification
- **Friend access patterns** for zero-copy GPU operations between components

For more details, see the `backend/distributed/unstructured` directory and its testsuite.

# Contributors
Ganellari Daniel, Zulian Patrick and Ramelli Dylan.

# License
The software is realized with NO WARRANTY and it is licenzed under [BSD 3-Clause license](https://opensource.org/licenses/BSD-3-Clause)

# Copyright
Copyright (c) 2015 Institute of Computational Science - USI Università della Svizzera Italiana, ETH-Z Eidgenössische Technische Hochschule Zürich

## Cite MARS ##

If you use the MARS Serial backend please use the following bibliographic entry

```
#!bibtex

@misc{mars_serial,
    author = {Zulian, Patrick and Ganellari, Daniel and Rovi, Gabriele and Ramelli, Dylan},
    title = {{MARS} - {M}esh {A}daptive {R}efinement for {S}upercomputing. {G}it repository},
    url = {https://bitbucket.org/zulianp/mars},
    year = {2018}
}
```

If you use the MARS Distributed backends (Kokkos, AMR and Unstructured) please use the following bibliographic entry


```
#!bibtex

@misc{mars_distributed,
    author = {Ganellari, Daniel and Zulian, Patrick and Rovi, Gabriele and Ramelli, Dylan},
    title = {{MARS} - {M}esh {A}daptive {R}efinement for {S}upercomputing. {G}it repository},
    url = {https://bitbucket.org/zulianp/mars},
    year = {2018}
}
```




