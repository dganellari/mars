# Changelog

All notable changes to MARS are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/), and the project aims to follow
[Semantic Versioning](https://semver.org/). While the major version is `0`, the
public API may change between minor releases.

## [Unreleased]

### Changed
- `mars_tgv` runs on `PeriodicNavierStokes` (`fem/mars_periodic_ns.hpp`): velocity and
  pressure both keep one unknown per periodic point, and every operator is Pᵀ A P on
  that space (`fem/mars_periodic_space.hpp`), so `D u = 0` holds exactly on any rank
  count and the multi-rank guard is gone. The BDF2 pressure right-hand side uses
  `3ρ / (2 dt)`, matching the corrector. Skew-symmetric advection is the default;
  `--solver=hypre` and `--pressure-solve=K` are no longer options of `mars_tgv`.

### Added
- `tests/periodic/check_periodic_space.py`: host model of the multi-rank periodic
  projection (no GPU) and the reference numbers of the TGV GPU check.

## [0.1.0] — 2026-09-24

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
  (config installed to `<prefix>/lib/cmake/Mars`, found through `CMAKE_PREFIX_PATH`).
- A plain `cmake ..` on a CPU-only machine builds the core library; the unstructured
  backend is on by default in CUDA / HIP builds.
- Release checks: `ctest -L release` runs the documented drivers end to end on generated
  meshes (assembly rank-count invariance, Poisson solve, cavity/channel Navier–Stokes).
- Hex CVFEM kernels transform reference gradients with the inverse-transpose Jacobian. Earlier
  code used the inverse, which is wrong whenever an element's reference axes are not aligned with
  x, y, z (typical of meshes from mesh generators).
- Multi-rank node ownership comes from the SFC decomposition: a node belongs to the rank whose SFC range contains
  it, which every rank computes without communication. The domain sync sends each element to the owners of its
  corners, so every owner holds all elements around its nodes and owned rows are complete. Earlier versions relied
  on the distance-based halo reaching those elements, which failed at corner contacts between ranks. Periodic and
  multi-block meshes keep the previous ownership scheme (see KNOWN_LIMITATIONS.md).
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
- `mars_poiseuille_flow` is a short teaching example on a small solver module,
  `fem/mars_channel_flow.hpp`: the planar CVFEM projection with constrained pressure
  gradients, inlet lift and opening fluxes, BDF2. Both linear systems are assembled once
  and solved with Hypre PCG + BoomerAMG (15-17 pressure iterations per step instead of
  thousands of Jacobi-CG iterations). The 1500-step check passes on 1, 2 and 4 GPUs
  with profile RMS error 4.553e-4 m/s; the full run takes about 25 s on one GPU. The
  10.9k-line channel solver fork and the `--planar-ddt`, `--pressure-amg` and
  `--velocity-amg` options are gone.

### Experimental
- High-order matrix-free CVFEM operators (p ≥ 2), with DOF numbering on the device.
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
