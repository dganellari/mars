# Known Limitations (v0.1.1)

MARS v0.1 is the first public release; v0.1.1 adds AMD GPUs (HIP). This document states plainly what is
validated, what is experimental, and what is not supported yet, so you can judge
whether MARS fits your use case. The major version is `0`: APIs may change.

## Stable — validated and supported
- Single-rank GPU-native unstructured mesh assembly pipeline
  (load → adjacency → DOF map → CSR sparsity → assembled matrix).
- Multi-rank distributed assembly and solve, periodic boxes included
  (lid-driven cavity, channel and Taylor–Green Navier–Stokes).
- **Incompressible Navier–Stokes on hex meshes** (`fem/mars_navier_stokes.hpp`), the
  solver of `mars_poiseuille_flow`, `mars_tgv` and `mars_lid_driven_cavity`, including
  periodic boxes. Validated:
  [Poiseuille](tests/reference/poiseuille/planar_validation.md) on 1, 2, 4, 16 and 64 GPUs against
  the exact profile, and the [periodic Taylor–Green vortex](docs/periodic_tgv_tutorial.md) against
  the viscous decay, with the same kinetic energy on 1 to 8 GPUs and within the release band on 16
  and 64. The release checks require the lid-driven cavity
  to give the same result on 1 and N ranks; it is not compared with a reference solution.

The stable paths are validated on generated structured meshes and the shipped Poiseuille mesh (release checks: `ctest -L release`),
including element numberings that are not aligned with the coordinate axes.

## Experimental — usable, not yet hardened
- **High-order matrix-free CVFEM (hexahedra, p = 1 to 8).** The release tests
  `marsReleaseHoMatfree` and `marsReleaseHoDistApply_npN` check the operator apply on one GPU
  and on several ranks; no linear solver uses it yet. On several ranks its high-order numbering
  (p > 1) needs the older vote node ownership, which `mars_ho_dist_apply_test` selects itself
  (`MARS_OWNERSHIP=vote`). Interfaces may change. See [the tutorial](docs/Matrix-Free-Tutorial.md).
- **Adaptive mesh refinement (AMR).** Single-rank mark/refine/rebuild/transfer works;
  multi-rank AMR is under development.
- **Tetrahedral high-order operators** (`mars_ho_laplacian_tet.hpp`, collapsed
  sum-factorization). Interfaces may change.
- **Coarse search and ghost registry** (`mars_coarse_search.hpp`,
  `mars_ghost_registry.hpp`). The device paths are gated against the host references.
- **Segregated SIMPLE solver** (`fem/segregated/`). Under development and not built by
  default (`-DMARS_ENABLE_SEGREGATED=ON`).
- **MARSIR** (`marsir-compiler/`, `marsir-mlir/`). Research code generator, off by
  default (`MARS_ENABLE_MARSIR`), not needed to build or use the library. CUDA only: its
  kernels are PTX loaded through the CUDA driver API.
  The BP1/BP3 example (`mars_marsir_bp`) runs only on its built-in deformed-cube mesh; its
  block-structured and generated DSS options run on one rank only; at scale the BoomerAMG coarse
  solve is the largest cost of a multigrid iteration.

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
- Kokkos structured-mesh / SFC path — optional, enabled with `MARS_ENABLE_KOKKOS`.
  Quad4/Hex8 generation, Morton/Hilbert ordering, Q1/Q2 scalar/two-component DOFs,
  scalar/two-component field gather/additive scatter and periodic cell halos have bounded CPU/MPI checks.
  Public-header and test-configuration failures from the baseline review are fixed.
  See [Structured Kokkos](docs/Structured-Kokkos.md) for the exact capability gates.
- `backend/adios2/` — optional ADIOS2 I/O, off by default.
- `backend/vtk/`, `moonolith_adapter/` — optional/experimental extensions, off by default.

Structured validation on 2026-11-04 uses macOS arm64, Open MPI 5.0.5 and Kokkos
4.4.01 at 1/2/4 ranks: final AppleClang 21 Debug builds use its Serial execution
space with direct and host-staged MPI; LLVM 19/libomp Release uses two OpenMP
threads with direct MPI. Installed consumers pass in both communication modes.
Each full configuration passes 11 CTests and 159 case/rank-count combinations.
Linux, CUDA and HIP remain
unvalidated. The [structured cluster kit](tests/structured/cluster/README.md)
prepares native CUDA and HIP builds with pinned Kokkos 4.4.01, direct and
host-staged MPI, 1/2/4-rank capability gates and an installed consumer. GPU source
cleanup and local CUDA syntax checks do not establish native GPU support.
Old PDE applications, staggered layouts, periodic DOF identification
and scaling are outside these checks. The shared `MARS_ENABLE_DISTRIBUTED_BACKEND`
option defaults to ON; Kokkos itself defaults to OFF.
Structured setup generates only local cell keys, selects exact balanced SFC
partitions with rank-sized metadata, and redistributes through sparse MPI peers.
DOF boundary scratch uses peer/DOF candidates. Device views contain no host context
or polymorphic owner, and data managers reuse communication buffers and requests.
These changes are locally validated, not measured scaling or native GPU results.
Plan creation is collective; start/wait is peer-only. Keep borrowed host contexts
alive and unmoved, use a new manager after changing a layout, and finish or destroy
plans before MPI/Kokkos finalize. The blocking cell-field and one-shot helper APIs
use the shared communicator and require consistently ordered calls.

## Assembly kernel variants
`fem/` ships several CVFEM assembly kernels. The canonical paths are the **hex tensor**
kernel (`mars_cvfem_hex_kernel_tensor.hpp`) and the **graph** kernels (hex/tet). The other
hex variants (`_wmma`, `_perip`, `_aos`, `_colored`, `_shmem`, `_optimized`) are
hardware-targeted optimizations of the same math. High-order matrix-free kernels are experimental (see above).
Unless you are benchmarking a specific GPU path, use the tensor or graph kernel.

## Build / platform notes
- Primary supported build: CUDA (`-DMARS_ENABLE_CUDA=ON -DMARS_ENABLE_UNSTRUCTURED=ON`)
  on NVIDIA GPUs, architectures `70;80;90` by default (override with
  `-DCMAKE_CUDA_ARCHITECTURES=...`).
- HIP (AMD, `-DMARS_ENABLE_HIP=ON`) is tested on MI300A (Beverin, ROCm 7.2) with Hypre 2.33 and
  netCDF 4.9.3 built for it: the production-integration run builds all examples and passes 48
  release tests plus 40 other tests; 12 other tests are skipped. The suite includes 1/4-rank checks.
  The FP64 tensor-core kernels (`wmma_tensor`, `wgmma_tensor` and the
  high-order DMMA operator) run on the AMD matrix cores (MFMA, two elements per wavefront). Not
  tested on AMD: more than one node, and performance beyond single kernels. Kernels tuned for
  NVIDIA can be slower. With the same P7 plane schedule on both, the high-order operator is
  1.3-1.6 times slower on MI300A than on H100, and 1.15-1.18 times without its atomic scatter
  (timing-only ablation): on MI300A every FP64 atomic goes past the L2. Its measured timings are
  operator results, not full-solver speedups. With its row map the `tensor` assembly, graph and full,
  on structured and unstructured hex meshes of 7.1e6 elements, runs within 3.4% between MI300A and H100
  (graph 4.82 and 4.88 ms on public hex192). The full assembly without a map, as callers that assemble
  once run it, takes 19.4 ms on MI300A and 9.6 ms on H100 on hex192: its global atomics all go past
  the L2 on MI300A. See `tests/performance/AMD-mfma.md` and `tests/performance/AMD-assembly.md` for
  the measured cases and validation scope.
  With Cray MPICH, multi-rank runs need `MPICH_GPU_SUPPORT_ENABLED=1`, and on MI300A
  `MPICH_GPU_IPC_ENABLED=0`: its on-node GPU IPC path crashes there.
- MPI is required by default (`-DMARS_ENABLE_MPI=ON`).
- Every rank must own part of the mesh. A mesh too small for the rank count stops on every rank with a
  message, instead of leaving a rank without elements.
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
  launcher (a binding wrapper, `CUDA_VISIBLE_DEVICES`, or `ROCR_VISIBLE_DEVICES` on AMD), otherwise
  every rank uses GPU 0.
- CMake fetches cornerstone-octree, and googletest when tests are on and no system googletest
  is found, at configure time, so a network connection is needed for a fresh configure.
- Without CUDA or HIP, `MARS_ENABLE_UNSTRUCTURED` defaults to OFF and a plain `cmake ..`
  builds only the core library. Its CPU tests are the MPI communication tests; the
  `test_install` target builds `examples/usage_from_external_cmake_project/` against an install.
- GPU builds: `ctest -L release` runs the documented drivers on generated meshes and the shipped Poiseuille mesh (see the
  README). The lower-level GPU domain tests still need a mesh directory in `MESH_PATH` and are
  skipped without one.
