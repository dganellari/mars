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

## Module status
The unstructured GPU backend (`backend/distributed/unstructured/`) is the active,
supported path. Other modules are present but build-OFF by default and not actively
developed for v0.1:
- `backend/serial/` — reference CPU backend, legacy (kept for debugging).
- Kokkos structured-mesh / SFC path — stable but not actively developed; prefer the unstructured backend.
- `backend/adios2/` — optional ADIOS2 I/O, off by default.
- `backend/vtk/`, `moonolith_adapter/` — optional/experimental extensions, off by default.

## Assembly kernel variants
`fem/` ships several CVFEM assembly kernels. The canonical paths are the **hex tensor**
kernel (`mars_cvfem_hex_kernel_tensor.hpp`) and the **graph** kernels (hex/tet). The other
hex variants (`_wmma`, `_perip`, `_aos`, `_colored`, `_shmem`, `_optimized`) are
hardware-targeted optimizations of the same math; `_core_example` is a reference/teaching
kernel, not for production. High-order matrix-free kernels are experimental (see above).
Unless you are benchmarking a specific GPU path, use the tensor or graph kernel.

## Build / platform notes
- Primary supported build: CUDA (`-DMARS_ENABLE_CUDA=ON -DMARS_ENABLE_UNSTRUCTURED=ON`)
  on NVIDIA GPUs, architectures `70;80;90` by default (override with
  `-DCMAKE_CUDA_ARCHITECTURES=...`). HIP (AMD) is supported via `-DMARS_ENABLE_HIP=ON`.
- MPI is required by default (`-DMARS_ENABLE_MPI=ON`).
- Dependencies (cornerstone-octree, googletest, google/benchmark) are fetched by CMake
  at configure time, so a network connection is needed for a fresh configure.
- The test suite and FEM examples are GPU-oriented and most require a mesh input and/or
  MPI; there is not yet a CPU-only smoke test.
