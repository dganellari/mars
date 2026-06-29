# Changelog

All notable changes to MARS are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/), and the project aims to follow
[Semantic Versioning](https://semver.org/). While the major version is `0`, the
public API may change between minor releases.

## [0.1.0] — 2026-06-29

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
- CMake install / `find_package(Mars)` packaging with the `Mars::mars` target.

### Experimental
- High-order matrix-free CVFEM operators (p ≥ 2).
- GPU-native adaptive mesh refinement (mark → refine → rebuild → transfer).

See [KNOWN_LIMITATIONS.md](KNOWN_LIMITATIONS.md) for the full stable / experimental /
unsupported breakdown.

[0.1.0]: https://github.com/dganellari/mars/releases/tag/v0.1.0
