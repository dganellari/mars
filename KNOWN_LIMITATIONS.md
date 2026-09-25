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

The stable paths are validated on generated structured meshes (release checks: `ctest -L release`),
including element numberings that are not aligned with the coordinate axes.

## Experimental — usable, not yet hardened
- **High-order matrix-free CVFEM (p ≥ 2).** Validated single-rank and at scale for the
  operator action; not yet a turnkey solver path. Interfaces may change.
- **Adaptive mesh refinement (AMR).** Single-rank mark/refine/rebuild/transfer works;
  multi-rank AMR is under development.
- **Tetrahedral high-order operators** (`mars_ho_laplacian_tet.hpp`, collapsed
  sum-factorization). Interfaces may change.
- **Coarse search and ghost registry** (`mars_coarse_search.hpp`,
  `mars_ghost_registry.hpp`). The device paths are gated against the host references.
- **Segregated SIMPLE solver** (`fem/segregated/`). Converges on the single-GPU public
  channel case; multi-rank runs, general meshes and field-level parity with a reference
  code are not validated yet.
- **MARSIR** (`marsir-compiler/`, `marsir-mlir/`). Research code generator, off by
  default (`MARS_ENABLE_MARSIR`), not needed to build or use the library.

## Not supported yet
- **Multi-rank periodic boundary conditions** (e.g. multi-rank periodic TGV). Periodic
  DOF collapse across rank boundaries is still under development; use single-rank for
  periodic cases.
- **Poiseuille channel (`mars_poiseuille_flow`).** Does not reproduce the tutorial's validated
  profile at v0.1.0, on one or several ranks: the velocity stops developing after the first step.
  A projection repair for this channel solver is in progress; treat the tutorial as a description
  of the method until a release notes it as fixed.
- **Triangle and quadrilateral meshes.** `ElementDomain` supports `TetTag` and `HexTag` only;
  `TriTag`/`QuadTag` are rejected at compile time.
- **Example-level restrictions.** `mars_cvfem_poisson` is single-rank and refuses more ranks.
  `mars_ex1_poisson` applies u = 0 on the faces of the mesh's bounding box, so it is correct for
  box-shaped domains only.

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
hardware-targeted optimizations of the same math. High-order matrix-free kernels are experimental (see above).
Unless you are benchmarking a specific GPU path, use the tensor or graph kernel.

## Build / platform notes
- Primary supported build: CUDA (`-DMARS_ENABLE_CUDA=ON -DMARS_ENABLE_UNSTRUCTURED=ON`)
  on NVIDIA GPUs, architectures `70;80;90` by default (override with
  `-DCMAKE_CUDA_ARCHITECTURES=...`). HIP (AMD) is enabled with `-DMARS_ENABLE_HIP=ON`, but
  the HIP build was not re-verified for v0.1.0, and the FEM examples are CUDA-only.
- MPI is required by default (`-DMARS_ENABLE_MPI=ON`).
- Exodus mesh input and side sets need netCDF. Without it MARS still builds; reading an Exodus
  mesh then fails at runtime with a clear error, and the binary directory format still works.
- Multi-GPU runs: `mars_cvfem_graph`, `mars_cvfem_graph_tet`, `mars_ex1_poisson`,
  `mars_cvfem_poisson` and `mars_amr_ns_projection` select GPU `rank % deviceCount`. Other
  drivers (e.g. `mars_tgv`) expect the launcher to expose one GPU per rank (a binding wrapper
  or `CUDA_VISIBLE_DEVICES`); otherwise every rank uses GPU 0.
- Dependencies (cornerstone-octree, googletest, google/benchmark) are fetched by CMake
  at configure time, so a network connection is needed for a fresh configure.
- Without CUDA or HIP, `MARS_ENABLE_UNSTRUCTURED` defaults to OFF and a plain `cmake ..`
  builds only the core library. Its CPU tests are the MPI communication tests plus the
  install smoke test in `examples/usage_from_external_cmake_project/`.
- GPU builds: `ctest -L release` runs the documented drivers on generated meshes (see the
  README). The lower-level GPU domain tests still need a mesh directory in `MESH_PATH` and are
  skipped without one.

## HO DOF numbering: single-rank GPU path exists, but is not the default

`HODofHandler` used to split its build paths by rank count, not by device:
`build()` (host) for one rank, `buildDistributedGpu*()` (device) for many. So
single-rank drivers numbered their DOFs on the host and uploaded `elemDof`,
which bounded single-GPU problem size and setup time by host numbering (the
distributed path measured 80 s -> 8 s per rank at 625M DOF/GPU when it moved to
the device).

The device twin now exists — `buildGpu()` / `buildGpuDevice()` in
`mars_ho_dof_handler_gpu.hpp`. They feed `buildDistributedGpuCore` the
degenerate single-rank configuration (`myRank = 0`, every corner and element
owned by 0, no shared corners, global id == local id), so there is still only
one numbering implementation. Equivalence to host `build()` is gated by
`mars_cvfem_ho_matfree_test --dof-self-check` on the permutation-invariant
quantities (`numDof`/`nEdge`/`nFace`, the `DofKey` multiset, and the `elemDof`
identification classes) — the DOF ids themselves are a permutation, as on the
distributed path.

`mars_cvfem_ho_matfree_test` now numbers on the device only: `buildGpu()` with
`keepOwn`, and the apply reads `HoOwnershipDeviceData::elemDof` in place, so
there is no host build and no `elemDof` H2D. Measured on GH200 at E=32 (32768
hexes, up to 11.4M DOF), device vs host numbering: 7.8x at p=1, 11.4x at p=2,
3.3x at p=7. The speedup falls with p because the host cost is dominated by the
p-independent edge/face `std::map` work, which amortizes as p grows; 3.3x is the
steady-state per-DOF figure (47 ns host vs 14 ns device).

The host `build()` remains, as the oracle that `--dof-self-check` scores the
device numbering against. That gate, and only that gate, still builds on the
host.

## Tet HO DOF numbering

The tet apply was always device-side; only the numbering
(`HoCvfemTetDofHandler::build` / `HoTetDofHandler::build`, both `std::map` over a
`vector<pair<gid,weight>>` key) ran on the host. `buildGpu()` in
`mars_ho_dof_handler_tet_gpu.hpp` is the device twin, and every tet driver now
calls it; the host build survives only as the `--dof-self-check` oracle.

Measured on GH200 at p=3, Kuhn mesh, device vs host numbering: 0.45x at 3.4k DOF
(launch-bound), 8.65x at 466k, 8.72x at 1.56M. The ratio plateaus because both
sides are sort-dominated. Per DOF the host costs ~480 ns against hex's ~47 ns --
the tet key is a heap-allocated vector per node -- so the port is worth more here
than it was for hex.

Build memory is ~48 B/node while numbering (four uint64 key lanes plus the
permutation and scan buffers), freed before return. That, not correctness, is the
scale limit of the current key packing.
