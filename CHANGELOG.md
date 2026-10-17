# Changelog

All notable changes to MARS are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/), and the project aims to follow
[Semantic Versioning](https://semver.org/). While the major version is `0`, the
public API may change between minor releases.

## [Unreleased]

## [0.1.1] — 2026-10-17

### Added
- HIP build for AMD GPUs, tested on MI300A: the tutorial examples pass their release checks on 1
  and 4 ranks, including the cavity, channel and Taylor–Green examples with Hypre built for HIP. The CUDA runtime,
  cuSPARSE/cuBLAS, Thrust and CUB calls map to HIP through one header layer on top of
  cornerstone's; warp-level kernels keep 32-lane groups inside AMD's 64-lane wavefronts. The
  README shows how to build Hypre 2.33 for HIP.

### Changed
- The matrix-free weak-scaling driver overlaps the halo exchange by default (`--no-overlap` to wait).
- GPU runtime calls in the domain setup (copies, allocations, synchronizations) are checked; a
  failure stops every rank with its location instead of continuing.

### Fixed
- AMR refinement levels (used by `mars_tgv` and the AMR drivers) read cornerstone's device copy of
  the tree leaves. They read the host copy, which only worked where the GPU can read CPU memory.
- The `team` assembly kernel indexes elements in 64 bits (its thread index wrapped above about
  134 million elements) and always launches the 256-thread block its shared memory is sized for
  (`--block-size` above 256 overflowed it).
- `HYPRE_DIR` is also searched in `lib64`.

## [0.1.0] — 2026-10-11

First tagged public release. MARS is a GPU-native mesh management and finite-element
assembly library for 3D tetrahedral and hexahedral meshes, built in C++20 on CUDA and the
cornerstone-octree library. The older serial and Kokkos backends, off by default and not
re-tested for this release, also handle simplices up to dimension 4.

### Added (stable)
- GPU-native unstructured meshes — mesh built and stored entirely on device, no host
  round-trips after load.
- Space-filling-curve (SFC) domain decomposition and load balancing via cornerstone.
- GPU-native finite-element / CVFEM assembly: element → DOF map → CSR sparsity →
  assembled matrix on device, with several optimized assembly kernels.
- Distributed multi-rank execution via MPI, including a per-node DOF halo for solver
  communication (CUDA-aware MPI) on top of the cornerstone element halo.
- The node exchanges keep one layer of ghost nodes current, 2.5 to 3 times less data than the whole
  cornerstone halo; `MARS_HALO_EXCHANGE=full` keeps the whole halo current.
- A mesh too small for the rank count stops every rank with a message instead of failing on a rank
  without elements. For small meshes the bucket of cornerstone's global tree shrinks, so every rank
  gets elements.
- Lazy composition of adjacency, halo, and coordinate caches (built on first access)
  to minimize VRAM and startup time.
- CMake install / `find_package(Mars)` packaging with the `Mars::mars` target
  (config installed to `<prefix>/lib/cmake/Mars`, found through `CMAKE_PREFIX_PATH`). The
  install holds the core library and the mesh headers, not yet the FEM and solver headers.
- A plain `cmake ..` on a CPU-only machine builds the core library; the unstructured
  backend is on by default in CUDA / HIP builds.
- Release checks: `ctest -L release` runs the documented drivers end to end on generated
  meshes and the shipped Poiseuille mesh (assembly rank-count invariance, Poisson solve, cavity/channel/Taylor–Green
  Navier–Stokes, high-order matrix-free gates on one GPU and on N ranks).
- The Poiseuille tutorial mesh ships in `tests/data/poiseuille/`; its validation against the
  analytic parabola is part of `ctest -L release` (1 and N ranks, needs netCDF).
- One incompressible Navier–Stokes solver for hex meshes, `fem/mars_navier_stokes.hpp`,
  runs `mars_poiseuille_flow`, `mars_tgv` and the new `mars_lid_driven_cavity`; the
  examples differ only in their boundary description (fixed velocity, p = 0,
  inlets/outlets, periodic pairs). CVFEM with equal-order nodes, Rhie–Chow face fluxes,
  a projection on the compact CVFEM Laplacian (the face fluxes are divergence-free to the
  solver tolerance on any rank count), skew-symmetric advection, BDF2. `DofSpace`
  (`fem/mars_dof_space.hpp`) maps node copies to unknowns: ghosts and periodic images of
  one point share one unknown, and every matrix is Pᵀ A P, formed by sending each copy's
  matrix row to the rank that owns its unknown. Both systems are solved with Hypre PCG +
  BoomerAMG (18 pressure iterations in the last step of the Poiseuille validation). Weak scaling on
  Alps GH200, 8M nodes per GPU: 72% of the one-node (4-GPU) step speed on 256 GPUs, Hypre-bound,
  with 20 to 23 pressure iterations at every size (Poiseuille tutorial, section 7).
- Weak scaling of MARS alone, the distributed p = 1 matrix-free operator apply at 16M elements per
  GPU, with the halo exchange overlapped (the default): 73% of the one-node rate on 1024 GPUs and
  69% on 2048 GPUs (32.8 billion DOFs); 65% and 45 to 62% without overlap (Matrix-Free tutorial).
- Validation: the Poiseuille 1500-step check (profile RMS error 4.551e-4 m/s, about 24 s on
  one GPU; `tests/reference/poiseuille/planar_validation.md`) passes on 1, 2, 4, 16 and 64 GPUs,
  and the Taylor–Green vortex kinetic energy follows the Stokes decay, identical on 1, 2, 4 and
  8 GPUs (8 over two nodes, `docs/periodic_tgv_tutorial.md`) and within the release band on 16
  and 64. `ctest -L release` runs the lid-driven
  cavity, the channel and the Taylor–Green vortex on 1 and N ranks.

### Changed
- Multi-rank node ownership comes from the SFC decomposition (single-block meshes; multi-block meshes keep the
  older vote ownership): a node belongs to the rank whose SFC range contains
  it, which every rank computes without communication. The domain sync sends each element to the owners of its
  corners, so every owner holds all elements around its nodes and owned rows are complete. Earlier versions relied
  on the distance-based halo reaching those elements, which failed at corner contacts between ranks. Periodic
  meshes use the same SFC ownership.

### Fixed
- Hex CVFEM kernels transform reference gradients with the inverse-transpose Jacobian. Earlier
  code used the inverse, which is wrong whenever an element's reference axes are not aligned with
  x, y, z (typical of meshes from mesh generators).
- `mars_ex1_poisson` numbers DOFs from the domain's node ownership and solves with CG and the node
  halo, like the Navier–Stokes solvers; before, its own DOF handler disagreed with the assembler's
  rows on more than one rank. The P1 tet assemblers now loop over every element a rank holds, halo
  elements included, and keep the ghost columns of owned rows; they used to drop both, which on more
  than one rank decoupled the ranks.
- `mars_cvfem_poisson` builds its DOF numbering, sparsity, boundary conditions and statistics on the
  GPU and runs on any number of ranks; before, it built them on the host, ran on one rank only, and
  handed the assembly kernels a host pointer that only unified-memory systems could read. The CG
  solver clips its Jacobi diagonal on the GPU instead of copying it to the host every solve.
- `mars_tgv` on several ranks: the old solver collapsed only the pressure at periodic
  points, and the multi-rank run lost the projection. Velocity and pressure now share one
  unknown per periodic point.

### Removed
- The cxxopts dependency and the legacy example programs that used it
  (`mars_examples` from `examples/mars_discretization.cpp`, `examples/mars_st_example.cpp`,
  `examples/serial/mars_serial_driver.cpp`); the examples parse their options directly.
- `examples/cvfem/`, a standalone STK-based assembly bench from another project.
- The 10.9k-line channel solver fork, `fem/mars_channel_flow.hpp`,
  `fem/mars_periodic_ns.hpp`, `fem/mars_periodic_space.hpp` and its host model
  `tests/periodic/`; the `--planar-ddt`, `--pressure-amg`, `--velocity-amg`, `--skew`,
  `--solver` and `--pressure-solve` options.

### Experimental
- High-order matrix-free CVFEM operators (hexahedra, p = 1 to 8), with DOF numbering on
  the device: the operator apply on one GPU and on several ranks, no solver yet
  (`docs/Matrix-Free-Tutorial.md`).
- Tetrahedral high-order operators: collapsed sum-factorization Galerkin and
  box-partition CVFEM.
- GPU-native adaptive mesh refinement (mark → refine → rebuild → transfer).
- Parallel AABB coarse search and a named per-interface ghost registry, each with a
  device path and a host reference.
- Segregated SIMPLE flow solver on the GPU (in development, built with `-DMARS_ENABLE_SEGREGATED=ON`).
- Multi-block Exodus side sets and a multi-state field history.
- MARSIR: an operator-spec → MLIR → tensor-core CUDA kernel generator
  (`marsir-compiler/`, `marsir-mlir/`). Research tooling, not part of the library
  build.

See [KNOWN_LIMITATIONS.md](KNOWN_LIMITATIONS.md) for the full stable / experimental /
unsupported breakdown.

[0.1.1]: https://github.com/dganellari/mars/releases/tag/v0.1.1
[0.1.0]: https://github.com/dganellari/mars/releases/tag/v0.1.0
