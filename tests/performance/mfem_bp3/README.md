# BP3 hierarchy comparison

This driver solves the dimensionless steady scalar diffusion problem on a deformed
unit cube with homogeneous Dirichlet data. Its continuous tensor-product Qp space
uses GLL nodes and p+2 Gauss points per axis. The benchmark right-hand side is the
fine discrete operator applied to a prescribed nodal sine field, b=A u. The true
residual and solution error are checked independently of the CG stopping norm.
These discrete-system checks do not establish PDE convergence order.

`--h-levels H` (default 0, maximum 3) retains H uniform-refinement levels below
p=1. The polynomial levels remain 1 -> 3 -> 7, truncated at the requested degree.
The number of final cells per axis must be divisible by 2^H, with at least two
cells per axis on the geometric coarse mesh and enough cells for the MPI partition.

The hierarchy is refined in logical coordinates. Each retained mesh is then
independently deformed, so the final fine geometry, boundary conditions, RHS and
operator agree with h0. Deforming only the coarse mesh before refinement changes
the benchmark and is not equivalent. After independent deformation the physical
coarse and fine meshes need not be nested; interpolation retains MFEM's logical
refinement map. Coarse operators use MFEM's native rediscretization, not the
Galerkin products of the MARS geometric path. Equal h depth therefore does not
mean equal preconditioner.

The h0 default preserves the previous mesh construction and numerical settings.
An unmodified-driver control is still needed when comparing timings, because the
optional instrumentation adds host wrapper calls even when disabled. This driver
supports a native MFEM configuration comparison, not a claim to have found the
best possible MFEM settings.

The driver destroys its solver objects before explicitly finalizing Hypre and
MPI in `main`, so MPI cleanup does not wait for process-exit callbacks. Apply
the same lifecycle change to the original-driver timing control.

## Small CPU/MPI gate

Configure against an MPI-enabled MFEM 4.8 build with `BP3_CPU_REFERENCE=ON`.
The ordinary build still requires CUDA. Set `BP3_PHASE_TESTS=ON` to build the four
host timer-boundary checks plus a memory-placement/transpose wrapper check, and run them with CTest.

Run `check_hierarchy.py EXECUTABLE NEW_OUTPUT_DIRECTORY` locally. It uses
`mpiexec --oversubscribe` on 1, 4 and 8 ranks; `--launcher` replaces that prefix.
Its 15 cases compare h0/h1/h2 at p3 and h0/h1 at p7 on the same small deformed
problem, including detailed phase accounting. This is a CPU reference gate;
GPU and GPU-aware MPI validation are separate.

`--audit-prefix PATH` exports owned true-node coordinates, essential-boundary
flags, exact nodal values, RHS, an independent probe, its fine operator action,
and the computed solution. It runs after all solver timings and numerical checks,
and refuses meshes larger than 16 cells per axis. The comparison canonicalizes
physical-coordinate keys, rejects duplicate ownership, and checks every field.
It must reject a deform-before-refine hierarchy even if that hierarchy solves its
own changed system successfully.

## Optional completed-work measurements

`BP3_PHASE_PREFIX=/path/prefix` records the last solve in per-rank JSON files.
`BP3_PHASE_DETAIL=1` additionally records operators, smoothers and the V-cycle.
The shared implementation and parser tests live in `../bp3_phases/`.

Both drivers explicitly drain pending device work before entering the coarse
operation. That drain has its own interval. Coarse and Hypre-adapter regions end
only when device work has completed. The MARS coarse region includes its data
conversion/DSS around the Hypre adapter; MFEM's adapter includes its own vector
handling. Neither interval is a pure Hypre kernel timer.

The detailed MFEM operator includes true-DOF assembly/exchange. MARS records its
local operator and the outer DSS separately; DSS inside smoothing remains in the
smoother/cycle interval. Consequently the raw operator intervals alone are not an
apples-to-apples compiler comparison. Timers are nested: inclusive times must not
be added. Exclusive times partition each rank's measured solve interval.

Use one rank, selected by maximum completed solve time, for all subtractions.
Do not subtract independently maximized regions. Solve minus coarse still
includes entry waits, reductions, transfers, halo work and other solver work;
it is not fine-kernel time.

Compare original-driver, new-driver probe-off, coarse-probe and detailed-probe
runs, with equal accuracy and matching iterations. A predeclared 3% change in
median warm whole-solve time rejects attribution to the uninstrumented execution.
Retain intrusive measurements as diagnostics and use a timeline/event experiment
instead. Report uninstrumented hierarchy timings separately.

## Preserve the native execution paths

Detailed wrappers forward additive and batched actions as well as ordinary and
transpose actions. Multigrid calls ArrayAddMult; falling back to Operator's
implementation allocates a fresh temporary instead of reusing the constrained
operator's scratch. Fine operators and smoothers stay unwrapped when detailed
profiling is off.

With geometric levels, `--h-transfer assembled` (default) builds MFEM's own true
transfer in Hypre ParCSR form. Its applications and transpose use the device
backend. `--h-transfer native` retains MFEM 4.8's TrueTransferOperator reference,
whose geometric path uses a host element loop. Both use the same interpolation
map; the direct test checks actions, transpose, constants and the global adjoint
relation across MPI ranks. Mesh warping and rediscretized coarse forms are
unchanged. The transfer changes setup cost and storage, which must be reported
alongside warm solves. GPU performance acceptance is separate from CPU parity.
