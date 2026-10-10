# Known Limitations (v0.1.0)

MARS v0.1.0 is the first public release. This document states plainly what is
validated, what is experimental, and what is not supported yet, so you can judge
whether MARS fits your use case. The major version is `0`: APIs may change.

## Stable — validated and supported
- Single-rank GPU-native unstructured mesh assembly pipeline
  (load → adjacency → DOF map → CSR sparsity → assembled matrix).
- Multi-rank distributed assembly and solve, periodic boxes included
  (lid-driven cavity, channel and Taylor–Green Navier–Stokes).
- **Incompressible Navier–Stokes on hex meshes** (`fem/mars_navier_stokes.hpp`), the
  solver of `mars_poiseuille_flow`, `mars_tgv` and `mars_lid_driven_cavity`, including
  periodic boxes. Validated on 1 to 8 GPUs:
  [Poiseuille](tests/reference/poiseuille/planar_validation.md) against the exact profile, and
  the [periodic Taylor–Green vortex](docs/periodic_tgv_tutorial.md) against the viscous decay,
  with the same kinetic energy on every rank count. The release checks require the lid-driven cavity
  to give the same result on 1 and N ranks; it is not compared with a reference solution.

The stable paths are validated on generated structured meshes (release checks: `ctest -L release`),
including element numberings that are not aligned with the coordinate axes.

## Experimental — usable, not yet hardened
- **High-order matrix-free CVFEM (hexahedra, p = 1 to 8).** The release tests
  `marsReleaseHoMatfree` and `marsReleaseHoDistApply_npN` check the operator apply on one GPU
  and on several ranks; no linear solver uses it yet. On several ranks it needs the older vote
  node ownership, which its drivers select themselves (`MARS_OWNERSHIP=vote`). Interfaces may
  change. See [the tutorial](docs/Matrix-Free-Tutorial.md).
- **Adaptive mesh refinement (AMR).** Single-rank mark/refine/rebuild/transfer works;
  multi-rank AMR is under development.
- **Tetrahedral high-order operators** (`mars_ho_laplacian_tet.hpp`, collapsed
  sum-factorization). Interfaces may change.
- **Coarse search and ghost registry** (`mars_coarse_search.hpp`,
  `mars_ghost_registry.hpp`). The device paths are gated against the host references.
- **Segregated SIMPLE solver** (`fem/segregated/`). Under development and not built by
  default (`-DMARS_ENABLE_SEGREGATED=ON`).
- **MARSIR** (`marsir-compiler/`, `marsir-mlir/`). Research code generator, off by
  default (`MARS_ENABLE_MARSIR`), not needed to build or use the library.

## Not supported yet
- **Navier–Stokes solver** (`fem/mars_navier_stokes.hpp`): hexahedral meshes only, and conforming
  meshes only. Meshes refined by the AMR module have hanging nodes, which the solver does not
  constrain yet. The solver is validated on box-shaped hexahedra. On distorted hexahedra the CVFEM
  matrices are not symmetric, but both systems are solved with PCG, which assumes symmetry; the
  solver prints a warning when it finds such a matrix.
- **Triangle and quadrilateral meshes.** Only tetrahedra and hexahedra are supported. A 2D problem
  can run as one layer of hexahedra, as the Poiseuille example does.
- **Example-level restrictions.** `mars_cvfem_poisson` and `mars_ex1_poisson` set u = 0 on the
  faces of the mesh's bounding box, so they solve the intended problem only on box-shaped domains.
  `mars_ex_beam_tet` and `mars_ex_beam_tet_distributed` run on one rank only.

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
- Multi-rank node exchanges keep one layer of ghost nodes current: the nodes of the rank's own
  elements and of the elements that touch a node it owns. Code that reads ghost values deeper in
  cornerstone's halo needs `MARS_HALO_EXCHANGE=full`, which keeps the whole halo current (2.5 to 3
  times more data). The AMR drivers, the development drivers (pump, coupled solver, outlet gates) and
  the segregated SIMPLE solver select it themselves.
- A Hypre built with 32-bit global indices (its default) limits a Navier–Stokes run to
  2^31 − 1 ≈ 2.1 billion unknowns per system; the solver stops with a message above that. Larger
  runs need Hypre built with 64-bit global indices (`--enable-mixedint`, in Spack `hypre+mixedint`).
- Exodus mesh input and side sets need netCDF. Without it MARS still builds; reading an Exodus
  mesh then fails at runtime with a clear error, and the binary directory format still works.
- Multi-GPU runs: the drivers of the release tests and the tutorials select GPU
  `rank % deviceCount`. Some research drivers do not; give each of their ranks one GPU with the
  launcher (a binding wrapper or `CUDA_VISIBLE_DEVICES`), otherwise every rank uses GPU 0.
- CMake fetches cornerstone-octree, and googletest when tests are on and no system googletest
  is found, at configure time, so a network connection is needed for a fresh configure.
- Without CUDA or HIP, `MARS_ENABLE_UNSTRUCTURED` defaults to OFF and a plain `cmake ..`
  builds only the core library. Its CPU tests are the MPI communication tests; the
  `test_install` target builds `examples/usage_from_external_cmake_project/` against an install.
- GPU builds: `ctest -L release` runs the documented drivers on generated meshes (see the
  README). The lower-level GPU domain tests still need a mesh directory in `MESH_PATH` and are
  skipped without one.
