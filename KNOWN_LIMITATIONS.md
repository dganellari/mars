# Known Limitations (v0.1.0)

MARS v0.1.0 is the first public release. This document states plainly what is
validated, what is experimental, and what is not supported yet, so you can judge
whether MARS fits your use case. The major version is `0`: APIs may change.

## Stable — validated and supported
- Single-rank GPU-native unstructured mesh assembly pipeline
  (load → adjacency → DOF map → CSR sparsity → assembled matrix).
- Multi-rank distributed assembly and solve for non-periodic cases
  (e.g. lid-driven cavity, channel Navier–Stokes).
- **Incompressible Navier–Stokes on hex meshes** (`fem/mars_navier_stokes.hpp`), the
  solver of `mars_poiseuille_flow`, `mars_tgv` and `mars_lid_driven_cavity`, including
  periodic boxes. Validated on 1, 2 and 4 GPUs:
  [Poiseuille](tests/reference/poiseuille/planar_validation.md) against the exact profile, and
  the [periodic Taylor–Green vortex](docs/periodic_tgv_tutorial.md) against the viscous decay,
  with the same kinetic energy on every rank count. The lid-driven cavity runs in the release
  tests on 1 and N ranks; it is not compared with a reference solution.

The stable paths are validated on generated structured meshes (release checks: `ctest -L release`),
including element numberings that are not aligned with the coordinate axes.

## Experimental — usable, not yet hardened
- **High-order matrix-free CVFEM (hexahedra, p = 1 to 8).** The release tests
  `marsReleaseHoMatfree` and `marsReleaseHoDistApply_npN` check the operator apply on one GPU
  and on several ranks; no linear solver uses it yet. Interfaces may change. See
  [the tutorial](docs/Matrix-Free-Tutorial.md).
- **Adaptive mesh refinement (AMR).** Single-rank mark/refine/rebuild/transfer works;
  multi-rank AMR is under development.
- **Tetrahedral high-order operators** (`mars_ho_laplacian_tet.hpp`, collapsed
  sum-factorization). Interfaces may change.
- **Coarse search and ghost registry** (`mars_coarse_search.hpp`,
  `mars_ghost_registry.hpp`). The device paths are gated against the host references.
- **Segregated SIMPLE solver** (`fem/segregated/`). Public-channel upwind and
  high-resolution fields agree with OpenAccel, with native 1/2/4-GPU rank parity.
  The [upwind duct study](tests/reference/openaccel/simple_duct/DAINT_RESULTS.md)
  passes refinement and rank comparisons. Accuracy on general meshes and multi-node
  scaling remain unvalidated. Native GPU SpMV is the default workaround for an
  unresolved residual mismatch in the Hypre/cuSPARSE path.
- **MARSIR** (`marsir-compiler/`, `marsir-mlir/`). Research code generator, off by
  default (`MARS_ENABLE_MARSIR`), not needed to build or use the library.

## Not supported yet
- **Navier–Stokes solver restrictions** (`fem/mars_navier_stokes.hpp`). Hex8 meshes
  only. On several ranks it needs SFC node ownership (the default); it stops under
  `MARS_OWNERSHIP=vote`. Planar mode (`mars_poiseuille_flow`) needs one layer of elements between two z
  planes. Meshes with hanging nodes are not supported, so `mars_tgv --adapt-every` gives
  wrong results: the solver does not constrain the hanging nodes that refinement leaves.
  Both systems are solved with PCG, which assumes a symmetric matrix; the CVFEM Laplacian
  is symmetric on the rectilinear meshes validated here but not in general on distorted
  hexes, which are not validated. On the 30k-node Poiseuille tutorial mesh more GPUs are
  slower, not faster; use `--cells` for scaling.
- **Triangle and quadrilateral meshes.** `ElementDomain` supports `TetTag` and `HexTag` only;
  `TriTag`/`QuadTag` are rejected at compile time.
- **Node ownership on multi-block meshes.** Multi-rank, single-block meshes, periodic ones included, give each
  node to the rank whose SFC range contains it and complete every owned node's element star during the domain
  sync, so owned rows are complete by construction. Multi-block (`MARS_BLOCK_NODE_IDENTITY`) meshes
  still use the previous scheme: the lowest claiming rank among halo peers owns a node, and the cornerstone halo
  search is widened by 1.5. That width is an empirical choice, not a guarantee; `MARS_ROW_DUMP` plus
  `tests/release/compare_rows.py` checks a mesh directly. `MARS_OWNERSHIP=vote` selects the previous scheme for
  every mesh.
- **Example-level restrictions.** `mars_cvfem_poisson` and `mars_ex1_poisson` apply u = 0 on the
  faces of the mesh's bounding box, so they are correct for box-shaped domains only. `mars_ex_beam_tet` and `mars_ex_beam_tet_distributed` are single-rank:
  their DOF handler (`UnstructuredDofHandler`) chooses node owners with its own rule, not the
  domain's. Multi-rank drivers number DOFs with `buildDofMappingGpu` from the domain's ownership,
  as `mars_ex1_poisson` and the Navier–Stokes solvers do.

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
- Multi-GPU runs: the multi-rank drivers of the release tests and the tutorials select GPU
  `rank % deviceCount`, except `mars_tgv`. `mars_tgv` and some research drivers expect the
  launcher to give each rank one GPU (a binding wrapper or `CUDA_VISIBLE_DEVICES`); otherwise
  every rank uses GPU 0.
- CMake fetches cornerstone-octree, and googletest when tests are on and no system googletest
  is found, at configure time, so a network connection is needed for a fresh configure.
- Without CUDA or HIP, `MARS_ENABLE_UNSTRUCTURED` defaults to OFF and a plain `cmake ..`
  builds only the core library. Its CPU tests are the MPI communication tests plus the
  install smoke test in `examples/usage_from_external_cmake_project/`.
- GPU builds: `ctest -L release` runs the documented drivers on generated meshes (see the
  README). The lower-level GPU domain tests still need a mesh directory in `MESH_PATH` and are
  skipped without one.

## High-order DOF numbering

The high-order DOF numbering runs on the device: `buildGpu()` and `buildDistributedGpu()` in
`mars_ho_dof_handler_gpu.hpp` for hexahedra, `buildGpu()` in `mars_ho_dof_handler_tet_gpu.hpp`
for tetrahedra. The drivers use these by default. The host builders (`HODofHandler::build()` and
`buildDistributed()`, and the tet `build()` of `mars_ho_dof_handler_tet.hpp`) remain as the
reference that the self-checks compare against (`mars_cvfem_ho_matfree_test --dof-self-check`,
`mars_ho_dist_apply_test --self-check`), and as the opt-in `--host-numbering` path of
`mars_ho_dist_apply_test`. The two numberings differ by a permutation, so the checks
compare permutation-invariant quantities: the DOF counts, the multiset of `DofKey`s and, for the
single-rank hex case, which element slots share a DOF.

The tet device numbering keeps about 48 bytes per element node while it runs (four 64-bit key
lanes plus the permutation and scan buffers, see `mars_ho_dof_handler_tet_gpu.hpp`) and frees
them before it returns. That memory, not correctness, limits the problem size of the current key
packing.
