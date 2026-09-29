# Prepare a user-owned OpenAccel case for native SIMPLE

`scripts/prepare_simple_deck.py` translates a restricted steady, laminar OpenAccel
YAML deck into native MARS command-line arguments. It is preparation only: Python
reads the deck, never the mesh. MARS's C++ Exodus reader, GPU topology checks and
GPU solver remain the execution path. Run the tool locally where the case lives;
keep private decks, generated arguments, logs and fields out of Git and review uploads.

The tool checks one fluid domain/material, constant positive density and dynamic
viscosity, zero initial velocity/pressure, fixed-frame physics, one normal-speed
inlet, one constant-static-pressure or average-static-pressure outlet and
stationary no-slip walls. It maps pseudo-time, relaxation, outlet blend and
upwind/high-resolution selection. Inlet/outlet `boundary_details.option` may be
omitted or explicitly `subsonic`; supersonic conditions remain unsupported.
Interpolation, pressure subiterations and expert settings must match the current
MARS implementation. In particular, OpenAccel defaults `relax_gradients` to true,
whereas this MARS path requires false; omission is rejected. Unknown physical
keys, expressions, moving walls, turbulence and ambiguous YAML are rejected.
Do not remove a rejected control just to get a passing preparation.

`velocity_interpolation_type: linear_linear` maps to
`--velocity-interpolation linear-linear`; `trilinear` remains the default.
This changes velocity sampling in fluxes, wall terms and gradient reconstruction,
not mesh coordinates. Pressure and both derivative interpolation settings still
require `linear_linear`. The numerical implementation requires a MARS rebuild;
see the [public shifted-interpolation check](velocity_interpolation.md) before
using the new mode on another case. Unknown expert keys remain rejected.

One explicit exception is `expert_parameters.coupled_pressure_velocity`: public
OpenAccel revision `0d69041` does not read this key anywhere in `src/`, so the
preparer accepts either boolean value without changing its arguments. This is
not a request for a monolithic pressure/velocity solver; native MARS still uses
SIMPLE. Non-boolean values and other unknown keys remain rejected. This exception
describes that reference revision, not versions that may implement the option.

`static_pressure` maps its constant `relative_pressure` in Pa to `--outlet-pressure`
and sets `--outlet-beta 1`. The existing trace law
`p_face = p_out + (1-beta)*(p_nearest-p_mean)` then fixes pressure at each open
outlet sample. The trace is set before the first pressure assembly; closed faces
retain that constant and use it again when reopening. At OpenAccel revision
`0d69041`, both pressure outlet modes share momentum/pressure assembly and
backflow selection. `pressure_profile_blend` is read only for
`average_static_pressure`; a leftover entry is ignored for `static_pressure`.
Time-dependent or spatially varying pressure input remains unsupported. This
mapping does not establish field parity on a new case.

Mesh coordinates must be in metres. The supplied mesh path must resolve to the
same file as the deck's `mesh.file_path` (relative to the saved deck directory).
Preparation follows path links but does not read or validate mesh bytes. Native
setup checks the single Tet4 block and complete side-set assignment before any
iteration; a preparation PASS alone is not mesh compatibility or field parity.
The optional string `mesh.automatic_decomposition_type` selects OpenAccel's
partitioner and is accepted without translation: MARS distributes the mesh with
Cornerstone. Other mesh controls, including transformations and decomposition
properties, remain rejected rather than silently omitted.

Boundary locations are resolved by the native reader using lowercase ASCII names
with spaces replaced by underscores, matching Ioss name normalization. It also
accepts `surface_<id>` and `sideset_<id>` aliases from the Exodus `ss_prop1` IDs;
IDs are not side-set sequence numbers. Blank or absent names require ID aliases.
As in Ioss, a stored `surface_<number>` name with a stale number is replaced by
the actual ID alias. Names or aliases that collide between side sets, repeated
selections of one side set, and unmatched selections are rejected. Every side
set must be selected, and GPU checks still require exact exterior-face coverage.
No names or counts are printed by these resolution errors. This native-reader
change requires rebuilding `mars_segregated_simple`; prepared arguments remain
usable. The rules follow Trilinos 16.2
[Ioex names and aliases](https://github.com/trilinos/Trilinos/blob/trilinos-release-16-2-0/packages/seacas/libraries/ioss/src/exodus/Ioex_DatabaseIO.C)
and [Ioss normalization](https://github.com/trilinos/Trilinos/blob/trilinos-release-16-2-0/packages/seacas/libraries/ioss/src/Ioss_Utils.C).
Split side-block aliases and arbitrary user aliases are not inferred.

Output scheduling and reference linear-solver settings are not translated. MARS
uses its own Hypre momentum/pressure solvers, true-residual checks and nonlinear
norms. `--reference-length` selects MARS's residual scale, not a physical model
parameter. Iteration budgets are specified at launch. Equal outer iteration counts
are useful diagnostic checkpoints, not equal physical time or a promise of equal
unconverged fields. The existing pinned public-channel comparator remains limited
to that public fixture; do not use it for other meshes.

OpenAccel allows named linear-solver definitions beside `solver_control` and
`output_control`, not just inline definitions inside
`solver_control.advanced_options.linear_solver_settings`. The preparer accepts
both forms. A named definition must be a mapping with a recognized `family`
(`petsc`, `hypre`, `trilinos`, `amgsolver` or `gmres`, case-insensitive).
Each `lookup` must resolve to such a definition; unused definitions are allowed.
Backend-specific options are not validated or copied into MARS. This follows
`linearSystem<N>::setupSolver` in public OpenAccel `0d69041` and does not claim
that MARS runs the reference's linear solver configuration.

`solver.restart_control` is different: its presence tells OpenAccel to load saved
fields. It is explicitly rejected, even if empty, because native SIMPLE currently
starts from zero fields. Unknown blocks are still rejected; they are not all
treated as output settings or linear-solver definitions.

## Preparation

Requires PyYAML in the preparation Python environment. The solver has no new
Python dependency. From the configured MARS CUDA/Hypre build directory, set
`SIMPLE_REFERENCE_DECK` to the actual saved deck and `MESH_BIG` to its existing
mesh path. These values remain on the user's machine. For example:

```bash
umask 077
simple_run=$(mktemp -d /capstor/scratch/cscs/gandanie/simple-private-XXXXXX)
MESH_BIG="$MESH_BIG" python3 ../scripts/prepare_simple_deck.py \
  --deck "$SIMPLE_REFERENCE_DECK" --mesh "$MESH_BIG" \
  --output "$simple_run/case" > "$simple_run/preparation.log" 2>&1
```

Stop on preparation failure. Validation collects independent compatibility issues
in one pass; checks that need a valid parent section cannot identify every error
inside a malformed section. Error messages use fixed public schema labels and
withhold custom key names, lookup names, values and file paths. File-access errors
give a generic message; check the paths and permissions locally. No argument files
are written when validation fails. `case.json`
records the deck hash and all arguments; `args.nul` provides arguments without
shell evaluation. Use a new output directory each time.

Preparation-only fixes do not require rebuilding the existing native executable.
The velocity-interpolation option changes numerical kernels and does require it:
`cmake --build . --parallel 4 --target mars_segregated_simple`.

## Direct short run

After successful preparation, in the same Bash shell:

```bash
mapfile -d '' -t simple_args < "$simple_run/case/args.nul"
unset CUDA_LAUNCH_BLOCKING MARS_OWNERSHIP MARS_HALO_FACTOR MARS_NODEHALO_ALLOW_INCONSISTENT
export MARS_HYPRE_SPMV_VENDOR=0 MARS_HYPRE_FLEXGMRES=0
export MARS_HYPRE_ABSTOL=0 MARS_HYPRE_VERBOSE=0 MARS_HYPRE_RESIDUAL_AUDIT=0
srun --account=csstaff --time=00:10:00 --nodes=1 --ntasks-per-node=4 \
  --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_segregated_simple \
  "${simple_args[@]}" --iterations 50 --report-every 10 \
  --residual-tol 1e-6 --mass-tol 1e-6 --change-tol 1e-6 \
  --linear-cache 1 --halo-overlap 1 --field-output distributed \
  --output-prefix "$simple_run/short" > "$simple_run/short.log" 2>&1
simple_status=$?
printf '%s\n' "$simple_status" > "$simple_run/short.exit"
```

Run the `srun` outside `set -e`, or temporarily disable it so exit 2 is recorded.
Zero means nonlinear convergence; 2 means the iteration budget was exhausted.
For a completed 50-step check with exit 2, require the exact final line beginning
`NOT CONVERGED: iteration limit iterations=50 ranks=4 `. Any other nonzero status
or missing final line is a failed short run. Do not label exit 2 convergence.
The solver checks finite fields, both linear residuals and global flux cancellation
throughout the run. Small continuity or mass imbalance is not required at startup.

Existing domain diagnostics print mesh metadata; the command stores all raw output
privately rather than echoing it to the terminal. Do not share those logs or field
files. No new mesh diagnostic or matrix dump is enabled. Only after compatibility
and short-run completion should a fresh run use the matched 2000-iteration budget;
this executable does not resume from the short run's CSV output.

## Shareable failure summary

Run the exporter on the machine holding the private log, after the job exits:

```bash
python3 ../scripts/simple_public_diagnostics.py \
  --log "$simple_run/short.log" --exit-file "$simple_run/short.exit" \
  --output "$simple_run/public-diagnostics.json"
```

Only `public-diagnostics.json` is intended for sharing or an exact-file transfer.
It contains fixed software labels and booleans: the failed stage, Krylov/SpMV
backend, residual-check verdicts, finite-value checks and completion status.
It contains no raw log lines, paths, mesh or boundary details, field values,
residual magnitudes, iteration counts, timings or input hashes. Full logs,
prepared arguments and fields remain private. No rebuild or simulation rerun is
needed; the exporter requires only Python 3.6 or later and its standard library.

Missing or conflicting diagnostics become `null` or `unknown`, never a pass.
The exporter combines observations without exposing their counts, refuses an
existing output file, and suppresses private values in its own error messages.
A successful export does not mean a successful solve. The absence of Hypre's
`false convergence` message does not rule out stagnation: that message depends
on verbosity. This summary identifies failure categories, not a root cause or
numerical parity. Add any new diagnostic through the fixed whitelist and its
redaction tests rather than uploading more of the raw log.

After rebuilding, `MARS_SIMPLE_PRESSURE_AUDIT=1` enables a failure-only device
audit of the original pressure matrix `A`, RHS `b` and halo-complete candidate
`phi`. The next export includes its booleans. The audit checks zero rows,
nonpositive diagonal entries, positive off-diagonals and whether every row
numerically annihilates the constant vector. This last check can reveal a missing
global pressure anchor; a negative result does **not** exclude an unanchored
disconnected component. Positive off-diagonals alone do not establish a bug in
a general CVFEM operator or prove AMG will fail.

It also compares the residual with a floating-point evaluation bound. For a row
with `m` stored entries, it computes
`gamma = (2*m+2)*u / (1-(2*m+2)*u)`, with double unit roundoff `u=2^-53`, and
`bound_i = gamma*(sum_j |A_ij|*|phi_j| + |b_i|)`.
The global L2 bound is compared with the residual and the unchanged application
limit. This is a cancellation warning, **not** a lower bound on attainable
accuracy, a condition-number estimate or permission to accept a rejected solve.
The `finite` flag also rejects overflow/underflow of the squared audit norms.

All row work and reductions use device buffers in CUDA builds, including the
MPI sums; only the fixed boolean report returns to the host for logging. Scratch
is allocated only after a rejected pressure solve. No field, row, norm, count or
geometry is printed, no successful iteration gains this work, and solve
acceptance is unchanged. A request on any rank enables the audit collectively.
Host references exercise the same row algebra on synthetic matrices. The
existing `mars_distributed_matrix_cuda_gate` includes these fixtures for
user-run GPU validation; host checks alone do not validate CUDA execution.

## Local checks

```bash
python3 -m unittest discover -s scripts -p test_prepare_simple_deck.py -v
python3 -m unittest discover -s scripts -p test_simple_public_diagnostics.py -v
```

The checked-in synthetic `preparation_fixture.i` tests cover extraction, static
and average pressure outlets, high-resolution selection, global
initialization, unsupported physics/interpolation/gradients/subiterations,
nonfinite values, duplicate YAML keys, unsafe tags, boundary ambiguity, literal
argument handling, path mismatch and existing-output rejection. An invalid dummy
mesh verifies that preparation never parses mesh contents. These are host checks,
not execution of a private case or a new CUDA validation.
Named/inline solver equivalence, missing/invalid lookups, unused definitions,
restart rejection, explicit subsonic boundaries, multiple simultaneous failures
and redacted CLI failures are covered separately.
