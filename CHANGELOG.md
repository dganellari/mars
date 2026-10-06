# Changelog

All notable changes to MARS are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/), and the project aims to follow
[Semantic Versioning](https://semver.org/). While the major version is `0`, the
public API may change between minor releases.

## [Unreleased]

## [0.1.0] — 2026-10-09

First tagged public release. MARS is a GPU-native mesh management and finite-element
assembly library for N-dimensional elements (N ≤ 4), built in C++20 on CUDA / HIP and
the cornerstone-octree library.

### Added (stable)
- GPU-native unstructured meshes — mesh built and stored entirely on device, no host
  round-trips after load.
- Space-filling-curve (SFC) domain decomposition and load balancing via cornerstone.
- GPU-native finite-element / CVFEM assembly: element → DOF map → CSR sparsity →
  assembled matrix on device, with several optimized assembly kernels.
- Distributed multi-rank execution via MPI, including a per-node DOF halo for solver
  communication (CUDA-aware MPI) on top of the cornerstone element halo.
- Lazy composition of adjacency, halo, and coordinate caches (built on first access)
  to minimize VRAM and startup time.
- CMake install / `find_package(Mars)` packaging with the `Mars::mars` target
  (config installed to `<prefix>/lib/cmake/Mars`, found through `CMAKE_PREFIX_PATH`). The
  install holds the core library and the mesh headers, not yet the FEM and solver headers.
- A plain `cmake ..` on a CPU-only machine builds the core library; the unstructured
  backend is on by default in CUDA / HIP builds.
- Release checks: `ctest -L release` runs the documented drivers end to end on generated
  meshes (assembly rank-count invariance, Poisson solve, cavity/channel/Taylor–Green
  Navier–Stokes, high-order matrix-free gates on one GPU and on N ranks).
- Hex CVFEM kernels transform reference gradients with the inverse-transpose Jacobian. Earlier
  code used the inverse, which is wrong whenever an element's reference axes are not aligned with
  x, y, z (typical of meshes from mesh generators).
- Multi-rank node ownership comes from the SFC decomposition: a node belongs to the rank whose SFC range contains
  it, which every rank computes without communication. The domain sync sends each element to the owners of its
  corners, so every owner holds all elements around its nodes and owned rows are complete. Earlier versions relied
  on the distance-based halo reaching those elements, which failed at corner contacts between ranks. Periodic
  meshes use the same SFC ownership.
- `mars_ex1_poisson` numbers DOFs from the domain's node ownership and solves with CG and the node
  halo, like the Navier–Stokes solvers; before, its own DOF handler disagreed with the assembler's
  rows on more than one rank. The P1 tet assemblers now loop over every element a rank holds, halo
  elements included, and keep the ghost columns of owned rows; they used to drop both, which on more
  than one rank decoupled the ranks.
- `mars_cvfem_poisson` builds its DOF numbering, sparsity, boundary conditions and statistics on the
  GPU and runs on any number of ranks; before, it built them on the host, ran on one rank only, and
  handed the assembly kernels a host pointer that only unified-memory systems could read. The CG
  solver clips its Jacobi diagonal on the GPU instead of copying it to the host every solve.
- The Poiseuille tutorial mesh ships in `tests/data/poiseuille/`; its validation run is
  opt-in with `-DMARS_ENABLE_VALIDATION_TESTS=ON`.
- One incompressible Navier–Stokes solver for hex meshes, `fem/mars_navier_stokes.hpp`,
  runs `mars_poiseuille_flow`, `mars_tgv` and the new `mars_lid_driven_cavity`; the
  examples differ only in their boundary description (fixed velocity, p = 0,
  inlets/outlets, periodic pairs). CVFEM with equal-order nodes, Rhie–Chow face fluxes,
  a projection on the compact CVFEM Laplacian (the face fluxes are divergence-free to the
  solver tolerance on any rank count), skew-symmetric advection, BDF2. `DofSpace`
  (`fem/mars_dof_space.hpp`) maps node copies to unknowns: ghosts and periodic images of
  one point share one unknown, and every matrix is Pᵀ A P, formed by sending each copy's
  matrix row to the rank that owns its unknown. Both systems are solved with Hypre PCG +
  BoomerAMG (18 pressure iterations per step in the Poiseuille validation). Weak scaling on
  Alps GH200, 8M nodes per GPU: 52% of the 1-GPU step speed on 256 GPUs, Hypre-bound, with
  21 to 23 pressure iterations at every size (Poiseuille tutorial, section 7).
- Validation on 1, 2 and 4 GPUs: the Poiseuille 1500-step check (profile RMS error
  4.551e-4 m/s, about 24 s on one GPU; `tests/reference/poiseuille/planar_validation.md`)
  and the Taylor–Green vortex kinetic energy against the Stokes decay, identical on 1, 2
  and 4 GPUs (`docs/periodic_tgv_tutorial.md`). `ctest -L release` runs the lid-driven
  cavity, the channel and the Taylor–Green vortex on 1 and N ranks.
- `mars_tgv` on several ranks: the old solver collapsed only the pressure at periodic
  points, and the multi-rank run lost the projection. Velocity and pressure now share one
  unknown per periodic point.
- Removed: the 10.9k-line channel solver fork, `fem/mars_channel_flow.hpp`,
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
- Segregated SIMPLE flow solver on the GPU, with a standalone public-channel driver.
- Multi-block Exodus side sets and a multi-state field history.
- MARSIR: an operator-spec → MLIR → tensor-core CUDA kernel generator
  (`marsir-compiler/`, `marsir-mlir/`). Research tooling, not part of the library
  build.

See [KNOWN_LIMITATIONS.md](KNOWN_LIMITATIONS.md) for the full stable / experimental /
unsupported breakdown.

[0.1.0]: https://github.com/dganellari/mars/releases/tag/v0.1.0
