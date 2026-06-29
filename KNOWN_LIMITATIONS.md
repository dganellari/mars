# Known Limitations (v0.1.0)

MARS v0.1.0 is the first public release. This document states plainly what is
validated, what is experimental, and what is not supported yet, so you can judge
whether MARS fits your use case. The major version is `0`: APIs may change.

## Stable — validated and supported
- Single-rank GPU-native unstructured mesh assembly pipeline
  (load → adjacency → DOF map → CSR sparsity → assembled matrix).
- Multi-rank distributed assembly and solve for non-periodic cases
  (e.g. lid-driven cavity, channel Navier–Stokes).
- Single-rank periodic Taylor–Green vortex.

## Experimental — usable, not yet hardened
- **High-order matrix-free CVFEM (p ≥ 2).** Validated single-rank and at scale for the
  operator action; not yet a turnkey solver path. Interfaces may change.
- **Adaptive mesh refinement (AMR).** Single-rank mark/refine/rebuild/transfer works;
  multi-rank AMR is under development.

## Not supported yet
- **Multi-rank periodic boundary conditions** (e.g. multi-rank periodic TGV). Periodic
  DOF collapse across rank boundaries is still under development; use single-rank for
  periodic cases.
- **Multi-rank Poiseuille channel.** Under investigation; single-rank works.

## Build / platform notes
- Primary supported build: CUDA (`-DMARS_ENABLE_CUDA=ON -DMARS_ENABLE_UNSTRUCTURED=ON`)
  on NVIDIA GPUs, architectures `70;80;90` by default (override with
  `-DCMAKE_CUDA_ARCHITECTURES=...`). HIP (AMD) is supported via `-DMARS_ENABLE_HIP=ON`.
- MPI is required by default (`-DMARS_ENABLE_MPI=ON`).
- Dependencies (cornerstone-octree, googletest, google/benchmark) are fetched by CMake
  at configure time, so a network connection is needed for a fresh configure.
- The test suite and FEM examples are GPU-oriented and most require a mesh input and/or
  MPI; there is not yet a CPU-only smoke test.
