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
- **Adaptive mesh refinement (AMR).** Mark/refine/rebuild/transfer works on one rank and runs
  distributed on cornerstone (hexahedra, red-refined tetrahedra). The refined solution has not
  been compared across rank counts, multi-rank tetrahedral AMR is only partly tested, and
  conforming tetrahedral bisection is not released.
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
- **General non-conformal flow coupling.** A scalar tetrahedral interface component
  pairs supplied nonmatching faces and forms diffusion stencils, checked on CUDA.
  Donor exchange, automatic CSR graph extension and owned-row insertion pass
  1/2/4-rank CUDA checks for the original dual-pass fixture. Independent passes do
  not conserve arbitrary nodal fields on nonmatching faces. A new coplanar shared-overlap
  path inserts equal and opposite scalar fluxes; host scalar, CUDA scalar and
  1/2/4-rank host/CUDA MPI checks pass. The new face-set search/coverage plan and complete
  scalar diffusion channel pass host and CUDA checks on 1/2/4 ranks, with equal and
  unequal diffusivities on procedural meshes.
  This path accepts stationary coplanar P1 triangular interfaces and supplied reciprocal
  candidate peers. Coverage checks compare the supplied faces; they cannot detect a region
  omitted from both side sets. The scalar mesh example accepts two disjoint Tet4 Exodus
  blocks through `prepare_nonconformal_mesh.py`, preserving file node and element indices.
  Input validation requires named inlet/outlet/wall and two interface sets covering all
  exterior faces. The file-based nonmatching Tet4 channel passes host and CUDA solves
  on 1/2/4 ranks, including analytic scalar, boundary-reaction and interface-flux checks.
  HIP validation remains pending. It replicates input staging
  for small meshes, then assembles owned rows and exchanges solution ghosts on device.
  It is not the scalable ElementDomain loader. Hex/tet interfaces, pyramids and
  general flow-driver integration remain unfinished. The separate scalar HIP checks remain unvalidated.
  An experimental `mars_nonconformal_flow` example adds steady constant-density,
  constant-viscosity channel flow on the PI's stationary 2 by 1 by 1 Tet4 mesh.
  It shares overlap mass fluxes between pressure correction and momentum advection,
  includes the coupled symmetric viscous stress, and uses physical pressure with a
  zero-pressure outlet. The inlet profile is 16 U y(1-y)z(1-z), with stationary no-slip
  walls; the plug-flow check instead prescribes uniform velocity on those walls.
  It requires the prepared input's `--flow-boundaries` face labels. Host flux,
  analytic plug-flow and MPI channel tests pass. The PI channel also passes CUDA
  on 1/2/4 ranks and the plug check on one GPU, with owned fields matching the host
  within 1.12e-16 absolute. The partitioned channel also passes user-run HIP
  comparisons on 1/2/4 ranks and the analytic plug check on one GPU.
  Its new partitioned-input mode reads only incident cells and required coordinates
  per rank, with persistent solution halos and explicit candidate peers. Host runs
  on 1/2/4 ranks reproduce the validated fields within 1.12e-16. User-supplied CUDA
  and HIP comparisons after the storage changes report maximum velocity difference
  2.78e-17 and pressure difference 1.67e-16 against the earlier CUDA reference.
  Raw files from these latest runs have not been independently retrieved.
  Offline preparation still reads the whole
  mesh and partitions by file-index ranges. This is not general ElementDomain
  integration or evidence of performance scaling. Its trace jump
  is weakly controlled, not an exact nodewise equality. On the current no-slip
  channel the maximum trace jump is 2.394e-3 and the normalized raw velocity
  divergence is 0.03835290 despite roundoff-level corrected-flux balance. Their
  accuracy on the original mesh is not established by conservation alone. A host
  procedural refinement study at n=2/4/8 now shows decreasing velocity and pressure
  differences against shared-node controls, including a fixed 2:1 nonmatching
  transverse resolution. Raw divergence decreases from n=4 to n=8, but is not
  monotone on the coarsest pair. See docs/Nonconformal-Channel-Flow.md.
  High-Reynolds flow, spatial
  convergence on general meshes, moving interfaces, mixed elements and pyramids
  are not validated by these checks.
- **Navier–Stokes solver** (`fem/mars_navier_stokes.hpp`): hexahedral meshes only, and conforming
  meshes only. Meshes refined by the AMR module have hanging nodes, which the solver does not
  constrain yet. The solver is validated on box-shaped hexahedra. On distorted hexahedra the CVFEM
  matrices are not symmetric, but both systems are solved with PCG, which assumes symmetry; the
  solver prints a warning when it finds such a matrix.
- **Triangle and quadrilateral meshes.** The unstructured GPU backend supports only tetrahedra and
  hexahedra; the structured Kokkos backend has Quad4. A 2D problem
  can run as one layer of hexahedra, as the Poiseuille example does.
- **Example-level restrictions.** `mars_cvfem_poisson` and `mars_ex1_poisson` set u = 0 on the
  faces of the mesh's bounding box, so they solve the intended problem only on box-shaped domains.
  `mars_ex_beam_tet` and `mars_ex_beam_tet_distributed` run on one rank only.

## Module status
The unstructured GPU backend (`backend/distributed/unstructured/`) is the active,
supported path. The Kokkos structured backend is also supported:
- Kokkos structured-mesh / SFC backend — optional, enabled with `MARS_ENABLE_KOKKOS`
  (OFF by default). It covers Quad4/Hex8 generation, Morton/Hilbert ordering,
  Q1/Q2 scalar/two-component DOFs, field gather/additive scatter and periodic cell
  halos, with direct and host-staged MPI. Its gates pass on Linux CPU, AMD MI300A
  and NVIDIA GH200 at 1/2/4 ranks over two nodes (2026-11-08); see
  [Structured Kokkos](docs/Structured-Kokkos.md#release-status). Not established:
  old PDE applications, staggered discretizations, periodic DOF identification,
  arbitrary polynomial orders and large-scale performance.

Other modules are present but build-OFF by default and not actively developed for v0.1:
- `backend/serial/` — reference CPU backend, legacy (kept for debugging).
- `backend/adios2/` — optional ADIOS2 I/O, off by default.
- `backend/vtk/`, `moonolith_adapter/` — optional/experimental extensions, off by default.
See the [release status](docs/Structured-Kokkos.md#release-status). The
[structured cluster kit](tests/structured/cluster/README.md) reproduces the native
builds and capability gates with pinned Kokkos 4.4.01.
The validated source uses the standard distributed partition for full Hilbert
meshes. All required structured correctness and installed-consumer checks pass
on CPU, AMD and NVIDIA. Four CUDA direct-MPI setup comparisons exceed the original
timing threshold and remain unresolved. They are accepted as non-blocking for
optional structured support, with the recorded timing rejections retained; no
performance-parity or scaling claim is made. See the release status for the
measured medians and the scope of the timing decision.
Old PDE applications, staggered layouts, periodic DOF identification
and scaling are outside these checks. The shared `MARS_ENABLE_DISTRIBUTED_BACKEND`
option defaults to ON; Kokkos itself defaults to OFF.
Structured setup generates only local cell keys, selects exact balanced SFC
partitions with rank-sized metadata, and redistributes through sparse MPI peers.
DOF boundary scratch uses peer/DOF candidates. Device views contain no host context
or polymorphic owner, and data managers reuse communication buffers and requests.
These changes pass the bounded CPU/CUDA/HIP gates above; scaling is not measured.
Plan creation is collective; start/wait is peer-only. Keep borrowed host contexts
alive and unmoved, use a new manager after changing a layout, and finish or destroy
plans before MPI/Kokkos finalize. The blocking cell-field and one-shot helper APIs
use the shared communicator and require consistently ordered calls.

The optional [conforming tetrahedral LEPP policy](docs/Conforming-Tet-AMR.md)
shares the legacy bisection core and adds sparse distributed closure, face-label
ancestry and dependent P1 transfer through cornerstone rebuilds. The policy
passes its CPU, NVIDIA GH200 and AMD MI300A gates at 1/2/4 ranks, including
two-node integration, with the published and pinned Cornerstone halo-depth
correction. Earlier procedural release checks passed: 10 CPU tests and all 50 GPU
release tests in both single-node and two-node layouts on each platform. It
supports single-block nonperiodic meshes with explicit triangle side sets;
node-only tags, multiblock/periodic meshes and coarsening are unsupported. GPU
validation uses 64-bit keys and double precision. Procedural labels and imported
Exodus triangle side sets pass on both GPU platforms. The Exodus fixture is a
single-block Tet4 cube with split coordinates in classic netCDF. The distributed
volume and boundary reader passes the combined AMR gate: thirteen checks on each
platform at 1/2/4 ranks, with CPU import/LEPP references and three GPU refinement
rounds per case. General Exodus layout coverage, scaling and cross-vendor
identical topology are not established.

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
  high-order DMMA operator) run on the AMD matrix cores (MFMA, two elements per wavefront). On AMD,
  two-node runs are correctness checks only (the 48 release tests in a two-node layout, and the
  structured suite); performance beyond single kernels is not tested. Kernels tuned for
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
