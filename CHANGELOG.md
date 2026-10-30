# Changelog

All notable changes to MARS are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/), and the project aims to follow
[Semantic Versioning](https://semver.org/). While the major version is `0`, the
public API may change between minor releases.

## [Unreleased]

## [0.1.1] — 2026-10-20

### Added
- HIP build for AMD GPUs, tested on MI300A: every example program builds, and all release checks
  (Navier–Stokes examples included, with Hypre and netCDF built for HIP) and tests pass on 1 and 4
  ranks. The CUDA runtime, cuSPARSE/cuBLAS, Thrust and CUB calls map to HIP through one header layer
  on top of cornerstone's; warp-level kernels keep 32-lane groups inside AMD's 64-lane wavefronts.
  The FP64 tensor-core kernels (`wmma_tensor`, `wgmma_tensor`, the high-order DMMA operator) run on
  AMD's matrix cores (MFMA). The README shows how to build Hypre 2.33 and netCDF for HIP.
- A step-by-step tutorial, ten steps and an appendix: the mesh, a first GPU loop,
  Poisson's equation with linear tetrahedra, the same with CVFEM on hexahedra, your own
  finite-element operator, the heat equation, the lid-driven cavity, Poiseuille flow, the
  Taylor–Green vortex, the high-order operator and adaptive refinement, and your own CMake project.
  The new programs (steps 1 to 6) have release checks on 1 and N ranks.
- `assembleElements` and `ElementValues` (`fem/mars_element_assembly.hpp`): assemble the element
  matrix and vector of any operator written as a small GPU functor, on linear tetrahedra and
  trilinear hexahedra. `assembleCvfemDiffusion` and `cvfemControlVolumes`
  (`fem/mars_cvfem_diffusion.hpp`): the CVFEM diffusion operator and control volumes in one call.

### Changed
- Device operations no longer carry redundant `Gpu` suffixes: for example,
  `buildDofMappingGpu` is now `buildDofMapping`. Geometry, DOF ownership and ACM helpers,
  examples and documentation use the same convention. High-order host reference builders
  use `buildHost`/`buildDistributedHost`; device builders use `buildDofs`/`buildDistributedDofs`,
  with `FromHost` wrappers for host inputs. Existing callers must update these names;
  the algorithms and execution paths are unchanged.
- The P7 FP64 Laplacian launchers use a plane schedule on AMD (MFMA) and on NVIDIA SM80+
  (`mma.sync`): each lane keeps u and the output along one direction in registers, and only the
  current plane of the element goes through shared memory. On the tested
  P7/E32 case a per-point and an affine apply take 0.43 and 0.21 ms on H100 (staged path 0.83 and
  0.69 ms) and 0.58 and 0.34 ms on MI300A (staged path 3.07 and 2.43-2.49 ms). Older NVIDIA GPUs keep the
  staged path. These are operator timings, not full-solver speedups.
- HIP `tensor` assembly scatters each row after its last face contribution. On public hex192,
  one MI300A rank, assembly time falls from 30.05 to 25.77 ms (14.2% less time), excluding
  output zeroing. Scratch storage falls from 1912 to 1084 bytes/lane; spills remain.
  Both optimizations retain the CUDA kernel paths. See `tests/performance/AMD-*.md` for scope.
- On AMD, reductions (sums, minima, maxima, counts) run on MARS kernels instead of rocPRIM, where
  each one compiled into about 130 kernels: `mars_amr_ns_projection` builds in under two minutes
  instead of over an hour. `-DMARS_HIP_ROCPRIM_REDUCTIONS=ON` goes back to rocPRIM; NVIDIA builds
  use Thrust as before.
- The `tensor` CVFEM assembly kernel takes the sub-control-surface nodes as compile-time values, so
  it indexes its per-node arrays with constants: its local memory per thread on H100 halves (1088 to
  512 bytes) and it runs 2.3 times faster; on MI300A it runs 6% faster but still spills registers.
  Same matrix on both.
- The Krylov solvers (CG, preconditioned CG, BiCGStab, GMRES) set up the sparse matrix-vector
  product once per solve instead of once per product. CG on 7.2 million unknowns takes 1.0 s
  instead of 19.9 s on MI300A, and 0.85 to 0.87 s instead of 0.97 s on H100.
- The matrix-free weak-scaling driver overlaps the halo exchange by default (`--no-overlap` to wait).
- GPU runtime calls in the domain setup (copies, allocations, synchronizations) are checked; a
  failure stops every rank with its location instead of continuing.

### Fixed
- AMR refinement levels (used by `mars_tgv` and the AMR drivers) read cornerstone's device copy of
  the tree leaves. They read the host copy, which only worked where the GPU can read CPU memory.
- The `team` assembly kernel indexes elements in 64 bits (its thread index wrapped above about
  134 million elements) and always launches the 256-thread block its shared memory is sized for
  (`--block-size` above 256 overflowed it).
- `HYPRE_DIR` is also searched in `lib64`.
- `smem_cache` assembly kernel: blocks with more distinct nodes than its cache holds read the rest
  from global memory, and its hash table holds every corner of a block (before, invalid reads or a
  hang). Its shared-memory alignment masks are 64-bit (they cleared the high address bits, which
  crashed on AMD), and it always launches the 256 threads it indexes with; so does
  `tensor_perip_lb2`, whose launch bounds allow no more.
- `wmma_tensor` zeroed twice the size of its scratch tiles, overwriting other per-element data.
- `ElementDomain::getNodeCount()` builds the local node numbering first: on several ranks it
  returned the count read from the file.
- `find_package(Mars)` on a HIP install recreates the ROCm library target and cornerstone's GPU
  library, which the installed targets link.
- The AMR bounding box reduced the device coordinates on the host.
- An installed MARS can build GPU programs against the unstructured backend: the headers are
  installed in the layout of the source tree, and `Mars::mars` carries the GPU definitions and nvcc
  flags; `find_package(Mars)` sets `MARS_GPU_LANGUAGE` and `MARS_GPU_ARCHITECTURES`. With CUDA,
  `libmars_unstructured` is device-linked when it is built. Tested on NVIDIA and AMD with the GPU
  program of the tutorial's appendix, built outside MARS and run on 4 ranks.

### Experimental
- MARSIR (CUDA only, off by default with `MARS_ENABLE_MARSIR`): in the element-local (cell-wise)
  solve every kernel now comes from the generator, including the DSS gather, the Jacobi
  preconditioner, the BiCGStab vector updates and the p-multigrid V-cycle. The V-cycle's p = 3 and
  p = 1 operators come from the same operator spec. At 64^3 on one GH200 the generated V-cycle
  is faster than the hand-written one (53 against 77 ms). Also new: a BoomerAMG coarse level,
  right-preconditioned BiCGStab with an FP32 V-cycle (one rank), and shared-memory layouts chosen
  against bank conflicts. Design notes and results: `docs/design/cellwise_unstructured_dss.md`.
- MARSIR (same option): the CEED benchmark problems BP1 and BP3 on element-local vectors
  (`mars_marsir_bp`), with the Galerkin operators generated from their specs (stored geometry, or
  trilinear geometry computed from the element corners) and communication-free flexible CG. BP3
  can be preconditioned by a p-multigrid V-cycle (degrees 7 -> 3 -> 1, Chebyshev smoothing,
  BoomerAMG on p = 1); at equal accuracy on one GH200 it solves 10x (p = 3) and 14x (p = 7) faster
  than DSS-Jacobi. The solve runs on several ranks (validated with BP3): 4 and 8 ranks reproduce
  the one-rank history to rounding (within one iteration with BoomerAMG), and with the same work
  per GPU, 64 GPUs keep 76% (p = 3) and 88% (p = 7) of the 8-GPU speed at equal accuracy. Also a block-structured DSS (one rank), hand-written and
  generated. Results: `docs/design/cellwise_bp_fcg.md`.

## [0.1.0] — 2026-10-11

First tagged public release. MARS is a GPU-native mesh management and finite-element
assembly library for 3D tetrahedral and hexahedral meshes, built in C++20 on CUDA and the
cornerstone-octree library. The older serial and Kokkos backends, off by default and not
re-tested for this release, also handle simplices up to dimension 4.

### Added (stable)
- GPU-native unstructured meshes — mesh built and stored entirely on device, no host
  round-trips after load.
- Space-filling-curve (SFC) domain decomposition and load balancing via cornerstone.
- GPU-native finite-element / CVFEM assembly: element → DOF map → CSR sparsity →
  assembled matrix on device, with several optimized assembly kernels.
- Distributed multi-rank execution via MPI, including a per-node DOF halo for solver
  communication (CUDA-aware MPI) on top of the cornerstone element halo.
- The node exchanges keep one layer of ghost nodes current, 2.5 to 3 times less data than the whole
  cornerstone halo; `MARS_HALO_EXCHANGE=full` keeps the whole halo current.
- A mesh too small for the rank count stops every rank with a message instead of failing on a rank
  without elements. For small meshes the bucket of cornerstone's global tree shrinks, so every rank
  gets elements.
- Lazy composition of adjacency, halo, and coordinate caches (built on first access)
  to minimize VRAM and startup time.
- CMake install / `find_package(Mars)` packaging with the `Mars::mars` target
  (config installed to `<prefix>/lib/cmake/Mars`, found through `CMAKE_PREFIX_PATH`). The
  install holds the core library and the mesh headers, not yet the FEM and solver headers.
- A plain `cmake ..` on a CPU-only machine builds the core library; the unstructured
  backend is on by default in CUDA / HIP builds.
- Release checks: `ctest -L release` runs the documented drivers end to end on generated
  meshes and the shipped Poiseuille mesh (assembly rank-count invariance, Poisson solve, cavity/channel/Taylor–Green
  Navier–Stokes, high-order matrix-free gates on one GPU and on N ranks).
- The Poiseuille tutorial mesh ships in `tests/data/poiseuille/`; its validation against the
  analytic parabola is part of `ctest -L release` (1 and N ranks, needs netCDF).
- One incompressible Navier–Stokes solver for hex meshes, `fem/mars_navier_stokes.hpp`,
  runs `mars_poiseuille_flow`, `mars_tgv` and the new `mars_lid_driven_cavity`; the
  examples differ only in their boundary description (fixed velocity, p = 0,
  inlets/outlets, periodic pairs). CVFEM with equal-order nodes, Rhie–Chow face fluxes,
  a projection on the compact CVFEM Laplacian (the face fluxes are divergence-free to the
  solver tolerance on any rank count), skew-symmetric advection, BDF2. `DofSpace`
  (`fem/mars_dof_space.hpp`) maps node copies to unknowns: ghosts and periodic images of
  one point share one unknown, and every matrix is Pᵀ A P, formed by sending each copy's
  matrix row to the rank that owns its unknown. Both systems are solved with Hypre PCG +
  BoomerAMG (18 pressure iterations in the last step of the Poiseuille validation). Weak scaling on
  Alps GH200, 8M nodes per GPU: 72% of the one-node (4-GPU) step speed on 256 GPUs, Hypre-bound,
  with 20 to 23 pressure iterations at every size (Poiseuille tutorial, section 7).
- Weak scaling of MARS alone, the distributed p = 1 matrix-free operator apply at 16M elements per
  GPU, with the halo exchange overlapped (the default): 73% of the one-node rate on 1024 GPUs and
  69% on 2048 GPUs (32.8 billion DOFs); 65% and 45 to 62% without overlap (Matrix-Free tutorial).
- Validation: the Poiseuille 1500-step check (profile RMS error 4.551e-4 m/s, about 24 s on
  one GPU; `tests/reference/poiseuille/planar_validation.md`) passes on 1, 2, 4, 16 and 64 GPUs,
  and the Taylor–Green vortex kinetic energy follows the Stokes decay, identical on 1, 2, 4 and
  8 GPUs (8 over two nodes, `docs/periodic_tgv_tutorial.md`) and within the release band on 16
  and 64. `ctest -L release` runs the lid-driven
  cavity, the channel and the Taylor–Green vortex on 1 and N ranks.

### Changed
- Multi-rank node ownership comes from the SFC decomposition (single-block meshes; multi-block meshes keep the
  older vote ownership): a node belongs to the rank whose SFC range contains
  it, which every rank computes without communication. The domain sync sends each element to the owners of its
  corners, so every owner holds all elements around its nodes and owned rows are complete. Earlier versions relied
  on the distance-based halo reaching those elements, which failed at corner contacts between ranks. Periodic
  meshes use the same SFC ownership.

### Fixed
- Hex CVFEM kernels transform reference gradients with the inverse-transpose Jacobian. Earlier
  code used the inverse, which is wrong whenever an element's reference axes are not aligned with
  x, y, z (typical of meshes from mesh generators).
- `mars_ex1_poisson` numbers DOFs from the domain's node ownership and solves with CG and the node
  halo, like the Navier–Stokes solvers; before, its own DOF handler disagreed with the assembler's
  rows on more than one rank. The P1 tet assemblers now loop over every element a rank holds, halo
  elements included, and keep the ghost columns of owned rows; they used to drop both, which on more
  than one rank decoupled the ranks.
- `mars_cvfem_poisson` builds its DOF numbering, sparsity, boundary conditions and statistics on the
  GPU and runs on any number of ranks; before, it built them on the host, ran on one rank only, and
  handed the assembly kernels a host pointer that only unified-memory systems could read. The CG
  solver clips its Jacobi diagonal on the GPU instead of copying it to the host every solve.
- `mars_tgv` on several ranks: the old solver collapsed only the pressure at periodic
  points, and the multi-rank run lost the projection. Velocity and pressure now share one
  unknown per periodic point.

### Removed
- The cxxopts dependency and the legacy example programs that used it
  (`mars_examples` from `examples/mars_discretization.cpp`, `examples/mars_st_example.cpp`,
  `examples/serial/mars_serial_driver.cpp`); the examples parse their options directly.
- `examples/cvfem/`, a standalone STK-based assembly bench from another project.
- The 10.9k-line channel solver fork, `fem/mars_channel_flow.hpp`,
  `fem/mars_periodic_ns.hpp`, `fem/mars_periodic_space.hpp` and its host model
  `tests/periodic/`; the `--planar-ddt`, `--pressure-amg`, `--velocity-amg`, `--skew`,
  `--solver` and `--pressure-solve` options.

### Experimental
- High-order matrix-free CVFEM operators (hexahedra, p = 1 to 8), with DOF numbering on
  the device: the operator apply on one GPU and on several ranks, no solver yet
  (`docs/Matrix-Free-Tutorial.md`).
- Tetrahedral high-order operators: collapsed sum-factorization Galerkin and
  box-partition CVFEM.
- GPU-native adaptive mesh refinement (mark → refine → rebuild → transfer).
- Parallel AABB coarse search and a named per-interface ghost registry, each with a
  device path and a host reference.
- Segregated SIMPLE flow solver on the GPU (in development, built with `-DMARS_ENABLE_SEGREGATED=ON`).
- Multi-block Exodus side sets and a multi-state field history.
- MARSIR: an operator-spec → MLIR → tensor-core CUDA kernel generator
  (`marsir-compiler/`, `marsir-mlir/`). Research tooling, not part of the library
  build.

See [KNOWN_LIMITATIONS.md](KNOWN_LIMITATIONS.md) for the full stable / experimental /
unsupported breakdown.

[0.1.1]: https://github.com/dganellari/mars/releases/tag/v0.1.1
[0.1.0]: https://github.com/dganellari/mars/releases/tag/v0.1.0
