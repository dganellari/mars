[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause) [![Documentation](https://readthedocs.org/projects/mesh-adaptive-refinement-for-supercomputing-mars/badge/?version=latest)](https://mesh-adaptive-refinement-for-supercomputing-mars.readthedocs.io/en/latest/?badge=latest)


# M.A.R.S #
## Mesh Adaptive Refinement for Supercomputing ##

**[Read the Full Documentation](https://mesh-adaptive-refinement-for-supercomputing-mars.readthedocs.io/en/latest/)**

MARS is an open-source, GPU-native mesh management library for N-dimensional elements
(N <= 4). It is developed in C++20 and makes heavy use of template meta-programming so
element dimensions, floating-point precision, and SFC key types are compile-time
parameters — giving both compile-time performance optimizations and concise, reusable
code.

The main features of MARS consist of:

1. GPU-native unstructured meshes — the mesh is built and stored entirely on the device
   (CUDA / HIP), with no host round-trips after load. Built on the cornerstone-octree
   library.

2. Space-filling-curve (SFC) domain decomposition — elements are identified by their
   lowest SFC corner key and load-balanced across ranks via cornerstone.

3. GPU-native finite-element and CVFEM assembly — element → DOF map → CSR sparsity →
   assembled matrix, all on the device, with multiple optimized assembly kernels
   (tensor-product, shared-memory, tensor-core variants).

4. Distributed multi-rank execution via MPI, including a per-node halo for solver
   communication (CUDA-aware MPI) on top of the cornerstone element halo.

5. GPU-native adaptive mesh refinement (mark → refine → rebuild → solution transfer);
   experimental, single-rank for now.

6. Lazy composition — adjacency, halo, and coordinate caches are built on first access
   to minimize VRAM and startup time.

MARS targets multi-core CPUs and GPUs (NVIDIA via CUDA, AMD via HIP); a Kokkos backend
covers the older structured-mesh path. Because the mesh stays on the device, libraries
built on MARS can run further operations directly on the GPU without going through the
host.

## Releases & status ##

Current release: **v0.1.0** — see [CHANGELOG.md](CHANGELOG.md). For the honest
stable / experimental / unsupported breakdown (what to rely on and what is still under
development), read [KNOWN_LIMITATIONS.md](KNOWN_LIMITATIONS.md). The major version is
`0`, so the public API may change between minor releases.

## Downloading MARS and its dependencies ##

Clone the repository. MARS has no git submodules — its dependencies
(cornerstone-octree, googletest, google/benchmark) are fetched automatically by
CMake at configure time, so a plain clone is all you need:

`git clone https://github.com/dganellari/mars.git`

A network connection is required at configure time for the dependency fetch.

GPU build (the main use case; for AMD use `-DMARS_ENABLE_HIP=ON` instead of the CUDA flags).
The unstructured backend is on by default in a GPU build:

```bash
cd mars
cmake -B build -DMARS_ENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=90
cmake --build build -j
```

CPU-only build. This builds the core library; the unstructured backend needs CUDA or HIP:

```bash
cd mars
cmake -B build
cmake --build build -j
```

### Checking your build

A CUDA build configured with `-DMARS_ENABLE_TESTS=ON -DMARS_ENABLE_FEM_EXAMPLES=ON` registers
release checks that run the documented drivers end to end on meshes generated at test time:

```bash
cd build
ctest -L release
```

They check that hex and tet assembly give the same matrix and RHS norms on 1 rank and on N
ranks (`-DMARS_RELEASE_TEST_RANKS=N`, default 4), that a CVFEM Poisson solve converges, and that
10-step lid-driven cavity and channel Navier–Stokes runs finish on 1 and N ranks without a failed
solve or NaN. They need python3 with numpy and an MPI launcher, and take a few minutes on one
GPU. ctest starts every GPU run through the MPI launcher CMake found (`mpiexec`, or `srun` on
Slurm), so on a Slurm cluster either run ctest inside an allocation:

```bash
salloc -A <account> -N 1 -t 00:30:00      # add the partition/GPU flags your site needs
ctest -L release -V
```

or give the launcher your site's flags once at configure time and run ctest from the login node:

```bash
cmake -B build -DMPIEXEC_EXECUTABLE=$(which srun) \
  "-DMPIEXEC_PREFLAGS=--account=<account>;--time=00:10:00;--nodes=1"
```

The drivers pick GPU `rank % deviceCount` themselves, so no GPU-binding wrapper is required. The
long (1-2 hour) Poiseuille validation against the analytic profile is opt-in: configure with
`-DMARS_ENABLE_VALIDATION_TESTS=ON`, then run `ctest -L validation`.

To use MARS from another CMake project, install it and point `CMAKE_PREFIX_PATH` at the
install prefix (see `examples/usage_from_external_cmake_project/`):

```bash
cmake --install build --prefix <prefix>
```

```cmake
find_package(Mars REQUIRED)
target_link_libraries(my_app PRIVATE Mars::mars)
```

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
- **Element Support**: Tetrahedra and hexahedra (triangle/quadrilateral tags exist but are not implemented yet)

### Quick Start

```cpp
#include "domain.hpp"

// Create GPU-native unstructured domain (read + partition + cstone sync).
// Template params: <ElementTag, RealType, KeyType, AcceleratorTag>.
ElementDomain<HexTag, double, uint64_t, cstone::GpuTag> domain("mesh_dir", rank, numRanks);

// Components built lazily on first access (all device-side):
const auto& offsets = domain.getNodeToElementOffsets();   // builds adjacency (CSR)
domain.cacheNodeCoordinates();                            // caches decoded coords
const auto& d_x = domain.getNodeX();                      // SoA node coordinates
const auto& d_conn = domain.getElementToNodeConnectivity(); // local node IDs per element
const auto& d_owner = domain.getNodeOwnershipMap();       // 0=ghost, 1=owned, 2=shared
```

### Build Configuration

Unstructured support:
- Needs `-DMARS_ENABLE_CUDA=ON` or `-DMARS_ENABLE_HIP=ON`; `MARS_ENABLE_UNSTRUCTURED`
  is then ON by default
- Cornerstone is fetched automatically if not found on the system

Example CMake command for unstructured with CUDA:

```bash
cmake .. \
  -DMARS_ENABLE_KOKKOS=OFF \
  -DMARS_ENABLE_CUDA=ON \
  -DMARS_ENABLE_TESTS=ON \
  -DMARS_ENABLE_UNSTRUCTURED=ON \
  -DMARS_ENABLE_FEM_EXAMPLES=ON \
  -DCMAKE_CUDA_ARCHITECTURES=90
make -j
```

`CMAKE_CUDA_ARCHITECTURES` is device-specific: 90 for GH200/H100, 80 for A100.

The block above already enables `MARS_ENABLE_FEM_EXAMPLES`, needed for the
CVFEM / FEM example drivers (Poisson, CVFEM assembly, the high-order
matrix-free gates). Other optional add-ons:

- `-DMARS_ENABLE_HYPRE=ON` — optional BoomerAMG-preconditioned solvers. The
  Navier–Stokes drivers default to CG and only need it for `--solver=hypre` (without
  it that option stops with an error); the segregated SIMPLE driver requires it.
- `-DMARS_ENABLE_ADIOS2=ON`, `-DMARS_ENABLE_VTK=ON` — extra I/O backends.

See `cmake/MarsOptions.cmake` and `cmake/MarsDependencies.cmake` for the full
option list.

### Documentation

For comprehensive guides and API references, see:

- **[Unstructured Meshes Documentation](https://mesh-adaptive-refinement-for-supercomputing-mars.readthedocs.io/en/latest/)**

**Getting started & tutorials**

- **[Quickstart](https://mesh-adaptive-refinement-for-supercomputing-mars.readthedocs.io/en/latest/Quickstart/)** - Clone, build, generate a mesh, run your first GPU assembly
- **[Poiseuille Channel Flow — a From-Scratch CFD Tutorial](https://mesh-adaptive-refinement-for-supercomputing-mars.readthedocs.io/en/latest/poiseuille_tutorial/)** - Incompressible Navier–Stokes, from the mesh to reading the output
- **[Taylor–Green Vortex (periodic)](https://mesh-adaptive-refinement-for-supercomputing-mars.readthedocs.io/en/latest/periodic_tgv_tutorial/)** - Canonical periodic validation case

**FEM**

- **[FEM Assembly](https://mesh-adaptive-refinement-for-supercomputing-mars.readthedocs.io/en/latest/FEM-Assembly/)** - Mesh → DOF map → sparsity → assembled CSR (GPU-native)
- **[CVFEM Kernels (GPU)](https://mesh-adaptive-refinement-for-supercomputing-mars.readthedocs.io/en/latest/CVFEM-Kernels/)** - Assembly kernel optimization variants

**Mesh & infrastructure**

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
    url = {https://github.com/dganellari/mars},
    year = {2018}
}
```

If you use the MARS Distributed backends (Kokkos, AMR and Unstructured) please use the following bibliographic entry


```
#!bibtex

@misc{mars_distributed,
    author = {Ganellari, Daniel and Zulian, Patrick and Rovi, Gabriele and Ramelli, Dylan},
    title = {{MARS} - {M}esh {A}daptive {R}efinement for {S}upercomputing. {G}it repository},
    url = {https://github.com/dganellari/mars},
    year = {2018}
}
```




