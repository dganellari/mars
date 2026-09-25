# Changelog

All notable changes to MARS are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/), and the project aims to follow
[Semantic Versioning](https://semver.org/). While the major version is `0`, the
public API may change between minor releases.

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
- The Poiseuille tutorial mesh ships in `tests/data/poiseuille/`; its validation run is
  opt-in with `-DMARS_ENABLE_VALIDATION_TESTS=ON`.

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
