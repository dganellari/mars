[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause) [![Documentation](https://readthedocs.org/projects/mesh-adaptive-refinement-for-supercomputing-mars/badge/?version=latest)](https://mesh-adaptive-refinement-for-supercomputing-mars.readthedocs.io/en/latest/?badge=latest)


# M.A.R.S #
## Mesh Adaptive Refinement for Supercomputing ##

**[Read the Full Documentation](https://mesh-adaptive-refinement-for-supercomputing-mars.readthedocs.io/en/latest/)**

MARS is an open-source, GPU-native mesh management library for N-dimensional elements
(N <= 4). It is written in C++20; the element type, the floating-point precision and the
SFC key type are template parameters.

The main features of MARS are:

1. GPU-native unstructured meshes — the mesh is built and stored on the device
   (CUDA / HIP), with no host round-trips after load. Built on the cornerstone-octree
   library.

2. Space-filling-curve (SFC) domain decomposition — elements are identified by their
   lowest SFC corner key and distributed across ranks by cornerstone.

3. GPU-native finite-element and CVFEM assembly — element → DOF map → CSR sparsity →
   assembled matrix, all on the device, with several assembly kernel variants.

4. An incompressible Navier–Stokes solver for hexahedral meshes (CVFEM, projection
   method, Hypre BoomerAMG), validated on Poiseuille flow and the Taylor–Green vortex.

5. Distributed multi-rank execution via MPI, including a per-node halo for solver
   communication (CUDA-aware MPI) on top of the cornerstone element halo.

6. Experimental: high-order matrix-free CVFEM operators, and GPU-native adaptive mesh
   refinement (mark → refine → rebuild → solution transfer).

7. Lazy composition — adjacency, halo, and coordinate caches are built on first access
   to reduce GPU memory use and startup time.

MARS targets multi-core CPUs and GPUs (NVIDIA via CUDA, AMD via HIP); a Kokkos backend
covers the older structured-mesh path. Because the mesh stays on the device, libraries
built on MARS can run further operations directly on the GPU without going through the
host.

## Releases & status ##

Current release: **v0.1.0** — see [CHANGELOG.md](CHANGELOG.md). For the stable /
experimental / unsupported breakdown (what to rely on and what is still under
development), read [KNOWN_LIMITATIONS.md](KNOWN_LIMITATIONS.md). The major version is
`0`, so the public API may change between minor releases.

## Downloading MARS and its dependencies ##

Clone the repository. MARS has no git submodules: CMake fetches cornerstone-octree at
configure time, and googletest when tests are enabled and no system googletest is found.
A plain clone is all you need:

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
ranks (`-DMARS_RELEASE_TEST_RANKS=N`, default 4), that the CVFEM Poisson solve and the P1 Poisson
example `mars_ex1_poisson` reach the expected maximum on 1 and N ranks, and that the high-order
matrix-free operator passes its gates on one GPU and on N ranks. With `-DMARS_ENABLE_HYPRE=ON`
they also run the Navier–Stokes examples on 1 and N ranks: 10-step lid-driven cavity and channel
runs that must finish without a failed solve or NaN, and 100 steps of the Taylor–Green vortex,
whose kinetic energy must follow the viscous decay. They need python3 with numpy and an MPI
launcher, and take a few minutes on one GPU. ctest starts every GPU run through the MPI launcher
CMake found (`mpiexec`, or `srun` on Slurm), so on a Slurm cluster either run ctest inside an
allocation:

```bash
salloc -A <account> -N 1 -t 00:30:00      # add the partition/GPU flags your site needs
ctest -L release -V
```

or give the launcher your site's flags once at configure time and run ctest from the login node:

```bash
cmake -B build -DMPIEXEC_EXECUTABLE=$(which srun) \
  "-DMPIEXEC_PREFLAGS=--account=<account>;--time=00:10:00;--nodes=1"
```

The drivers pick GPU `rank % deviceCount` themselves. The Poiseuille validation
against the analytic profile (1500 steps) is opt-in: configure with
`-DMARS_ENABLE_VALIDATION_TESTS=ON`, then run `ctest -L validation` in a GPU allocation.

To use MARS from another CMake project, install it and point `CMAKE_PREFIX_PATH` at the
install prefix (see `examples/usage_from_external_cmake_project/`). The install holds the core
library and the mesh headers; the FEM and solver headers (`fem/`, `solvers/`) are not installed
yet, so code that uses the Navier–Stokes solver builds inside the source tree, like the examples:

```bash
cmake --install build --prefix <prefix>
```

```cmake
find_package(Mars REQUIRED)
target_link_libraries(my_app PRIVATE Mars::mars)
```

## MARS Kokkos requirements ##

The Kokkos backend (structured meshes, off by default) depends on the Kokkos and Kokkos
Kernels libraries. The unstructured GPU backend does not need Kokkos.

MARS finds Kokkos if it is installed on your system, standalone or as part of Trilinos. It
looks for KOKKOS_DIR or TRILINOS_DIR in the environment: with Trilinos it takes Kokkos from
$TRILINOS_DIR, otherwise it looks for Kokkos and Kokkos Kernels at $KOKKOS_DIR.

Use -DMARS_ENABLE_KOKKOS=ON to enable it. For more details check CMakeLists.txt.

The default when compiling MARS with Kokkos without specifying any other CMake flag is the
Kokkos/OpenMP execution space. Kokkos should also be compiled with OpenMP support. Otherwise the
default is the serial execution space.

To compile for CUDA the CMake flag MARS_ENABLE_CUDA=ON must be set. An example would be:
```
cmake -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_BUILD_TYPE=Release -DMARS_ENABLE_KOKKOS=ON -DMARS_ENABLE_CUDA=ON ..
```

If compiled for CUDA then Kokkos should also be compiled with CUDA (Kokkos_ENABLE_CUDA=ON) and CUDA_LAMBDA (Kokkos_ENABLE_CUDA_LAMBDA=ON) support.

## Unstructured Mesh Support ##

MARS supports GPU-native unstructured meshes through the Cornerstone library: the mesh is
partitioned along a space-filling curve (SFC) and managed on the GPU, on one rank or many.

### Key Features

- **GPU-Native Architecture**: All data structures live in device memory (`DeviceVector` via Cornerstone)
- **SFC-Based Partitioning**: Elements identified by space-filling curve keys and distributed along the curve
- **Lazy Composition**: Components (adjacency, halo, coordinates) allocated on demand to reduce GPU memory use
- **Thrust Algorithms**: CSR building, sorting, and reductions use GPU Thrust primitives
- **MPI Integration**: Multi-rank support via Cornerstone domain decomposition
- **Element Support**: Tetrahedra and hexahedra (triangle/quadrilateral tags exist but are not implemented yet)

### Quick Start

```cpp
#include "domain.hpp"

// Create GPU-native unstructured domain (read + partition + cstone sync).
// Template params: <ElementTag, RealType, KeyType, AcceleratorTag>.
ElementDomain<HexTag, double, uint64_t, cstone::execution::Gpu> domain("mesh_dir", rank, numRanks);

// Components built lazily on first access (all device-side):
const auto& offsets = domain.getNodeToElementOffsets();   // builds adjacency (CSR)
domain.cacheNodeCoordinates();                            // caches decoded coords
const auto& d_x = domain.getNodeX();                      // SoA node coordinates
const auto& d_conn = domain.getElementToNodeConnectivity(); // local node IDs per element
const auto& d_owner = domain.getNodeOwnershipMap();       // 1 = owned, 0 = ghost
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
matrix-free tests). Other optional add-ons:

- `-DMARS_ENABLE_HYPRE=ON` — Hypre BoomerAMG. The Navier–Stokes examples
  (`mars_poiseuille_flow`, `mars_tgv`, `mars_lid_driven_cavity`) and the segregated SIMPLE
  driver are built only with it.
- netCDF, when CMake finds it (pkg-config, or the `NETCDF_DIR` / `NETCDF_ROOT` environment
  variable), enables reading Exodus meshes such as the Poiseuille tutorial mesh. Without it
  the binary mesh directory format still works.
- `-DMARS_ENABLE_ADIOS2=ON`, `-DMARS_ENABLE_VTK=ON` — extra I/O backends.

See `cmake/MarsOptions.cmake` and `cmake/MarsDependencies.cmake` for the full
option list.

### Documentation

The [documentation](https://mesh-adaptive-refinement-for-supercomputing-mars.readthedocs.io/en/latest/)
has a quickstart, tutorials for v0.1 and reference pages:

**Getting started & tutorials**

- **[Quickstart](https://mesh-adaptive-refinement-for-supercomputing-mars.readthedocs.io/en/latest/Quickstart/)** - Clone, build, generate a mesh, run your first GPU assembly
- **[Poiseuille Channel Flow](https://mesh-adaptive-refinement-for-supercomputing-mars.readthedocs.io/en/latest/poiseuille_tutorial/)** - Incompressible Navier–Stokes, from the mesh to the validated result and weak scaling to 256 GPUs
- **[Taylor–Green Vortex (periodic)](https://mesh-adaptive-refinement-for-supercomputing-mars.readthedocs.io/en/latest/periodic_tgv_tutorial/)** - Periodic boundaries, one unknown per periodic point on any number of GPUs
- **[High-Order Matrix-Free Operator](https://mesh-adaptive-refinement-for-supercomputing-mars.readthedocs.io/en/latest/Matrix-Free-Tutorial/)** - The experimental high-order CVFEM operator apply on one and many GPUs

**Reference**

- **[FEM Assembly](https://mesh-adaptive-refinement-for-supercomputing-mars.readthedocs.io/en/latest/FEM-Assembly/)** - Mesh → DOF map → sparsity → assembled CSR (GPU-native)
- **[CVFEM Kernels (GPU)](https://mesh-adaptive-refinement-for-supercomputing-mars.readthedocs.io/en/latest/CVFEM-Kernels/)** - Assembly kernel variants
- **[Known limitations](KNOWN_LIMITATIONS.md)** and **[Changelog](CHANGELOG.md)**

### Implementation Details

The unstructured backend uses:
- **Lazy initialization** for memory efficiency (adjacency, halo, coordinates built on demand)
- **SFC keys as connectivity** for sparse global element identification
- **Thrust-based CSR** building via `sort_by_key`, `reduce_by_key`, `exclusive_scan`
- **Lowest SFC corner** representation (not centroids) for element identification

For more details, see the `backend/distributed/unstructured` directory and its testsuite.

# Contributors
Ganellari Daniel, Zulian Patrick, Ramelli Dylan and Rovi Gabriele.

Everyone who has committed to MARS since 2018:
[contributors graph](https://github.com/dganellari/mars/graphs/contributors?from=2018-06-01)
(GitHub's default view shows only the last three months).

# License
The software is released with NO WARRANTY and it is licensed under [BSD 3-Clause license](https://opensource.org/licenses/BSD-3-Clause)

# Copyright
Copyright (c) 2015 ETH-Z Eidgenössische Technische Hochschule Zürich, Institute of Computational Science - USI Università della Svizzera Italiana

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
