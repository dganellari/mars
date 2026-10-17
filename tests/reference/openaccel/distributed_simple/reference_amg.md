# First-step comparison with the reference AMG choices

The public `simple-reference-pressure-check-F87prv` report confirms that the
PETSc momentum predictor and both assembled equations match at the first step.
The pressure increments differ, although both satisfy the same saved reference
residual target on the common unscaled system. All twenty printed OpenAccel
pressure targets were met in the original-control history, without a cap or an
ambiguous rounded comparison. Those logged values are not independently
recomputed residuals for every step.

The normal GPU pressure profile deliberately substitutes PMIS for HMIS and
l1-Jacobi for the direct coarse solve. `prepare --reference-amg` removes those
two substitutions while keeping the original pressure and momentum targets,
restart, iteration cap, other supported AMG controls and PETSc momentum profile.
It reuses the saved OpenAccel first step and runs one new MARS step. Recovery
remains disabled. The default device profile is unchanged.

## Algebra and scope

The equation remains `A phi = b`, with the nodal Tet4 CVFEM physical pressure
increment `phi` in Pa. SIMPLE still applies `p += alpha_p phi` and
`u -= d grad(phi)`. Boundary contributions, shifted quadrature, volume scaling,
density, viscosity, pseudo-time and relaxation are unchanged. Only construction
and application of the pressure preconditioner change. Both independent residual
checks still use the original target `max(atol, rtol ||b||)`.

The executable option is `--pressure-amg-mode reference`; it requires an explicit
pressure profile and a one-step audit. The preparation wrapper also requires
`--both-linear-solvers`. It rejects history mode, unsupported settings and changed
baseline evidence before launch. There is no tolerance change or new acceptance
rule. A rebuilt MARS executable is permitted; the captured dependency identities,
rank count, launcher and solver environment remain checked.

**This mode permits CPU work inside Hypre.** In Hypre 2.33,
[HMIS](https://github.com/hypre-space/hypre/blob/v2.33.0/src/parcsr_ls/par_coarsen.c)
clones device matrix structure to host for its Ruge phase and returns the markers
to device; its subsequent PMIS phase follows the execution policy.
[Coarse solver 9](https://github.com/hypre-space/hypre/blob/v2.33.0/src/parcsr_ls/par_gauss_elim.c)
stages the coarse matrix and vectors on host for Gaussian elimination. These are
library-imposed transfers, not a GPU performance optimization. No library rebuild
is needed; this mode requires Hypre 2.33 or newer.

Keeping these choices does not make CPU and GPU solvers identical. Hypre versions,
ordering, ownership and implementation paths can still differ, including the
device PMIS phase inside HMIS. The experiment tests the two explicit substitutions;
it cannot guarantee field agreement or establish nonlinear convergence.

Each rank's post-setup settings must retain HMIS 10 and effective coarse solver 9.
Hypre can replace the coarse solver when coarsening stops early. Such a fallback
fails this comparison, because it would not test the requested direct solve.
`recorded_pressure_settings_verified: true` certifies the selected settings;
`recorded_gpu_pressure_settings_verified: false` in this mode records the explicit
permission for host work, rather than a failed settings check.

## Validation

The host-adapted production wrapper passes its numeric-refresh regression with
real sequential Hypre 2.33.0 and 3.1.0 under AddressSanitizer and
UndefinedBehaviorSanitizer. The added strong-diffusion cases exercise GMRES and
FlexGMRES with a multilevel hierarchy, actual HMIS/direct settings, a prepared
host coarse solve, two numeric updates and a known independent solution.
The existing device-profile and rejection cases remain included.

The Python comparison tests cover the unchanged baseline, original targets,
mixed PETSc/Hypre profile, host-work disclosure before launch, changed saved
arguments and rejected runtime substitutions. C++ option tests reject use outside
the first-step audit. These local checks do not validate CUDA execution or MPI.

The focused CUDA gate below checks the actual controls, more than one level,
host coarse preparation, independently known solutions, exchanged ghosts and
numeric refresh for both Krylov methods. It runs on 1, 2 and 4 GPUs before the
private comparison. A failed gate stops the recipe.

### Coarse-grid cutoff in the public gate

The first GPU run passed all 16 checks on one rank. On two ranks, GMRES and all
independent residual/solution checks passed, but the first FlexGMRES setup used
effective coarse smoother 6 instead of 9; its numeric refresh used 9 again.
The gate correctly rejected that setup. Under the supplied recipe, this stops
before the four-rank gate and private flow launch.

[Hypre 2.33 setup](https://github.com/hypre-space/hypre/blob/v2.33.0/src/parcsr_ls/par_amg_setup.c)
replaces a direct coarse solver with the generic smoother if another coarsening
attempt produces no coarse points or retains every point. Each new setup restores
the user-requested coarse solver first. A multilevel hierarchy therefore does not
prove that the direct solve remained active. The output identifies this fallback;
it does not identify which of those two coarsening conditions occurred.

The synthetic reference-AMG fixture now uses a cutoff of `N-1` rows (`119` for
its public 120-row system), rather than copying the reference cutoff of 9.
The connected fine graph must genuinely shrink once. The gate requires exactly
two levels, `0 < coarse_rows < N`, actual coarse solver 9 and prepared host direct
storage. This tests the direct-solve path before a small coarse graph can lose its
strong couplings. Residual, known-solution, ghost and numeric-refresh checks remain
unchanged. The private profile still uses its original cutoff, and its strict
post-setup check still rejects an effective solver other than 9.

Updated host-adapted wrapper regressions pass with Hypre 2.33 and 3.1. An all-weak
matrix deliberately reaches the fallback to 6; GMRES and FlexGMRES both solve it
accurately, but the unchanged profile check rejects it. A separate real Hypre 2.33
CPU probe passes the direct and fallback cases on 1/2/4 MPI ranks, including numeric
updates and independent residual/solution checks. The corrected CUDA gate then
passed all 16 checks on each of 1, 2 and 4 GPUs in the user-supplied results.
Those 48 passing checks cover the public fixture; private field parity is pending.

### Locating the private first-step timeout

The subsequent private run reached its scheduler time limit. The public file
checks show an initial snapshot and momentum settings on all four ranks, but no
pressure settings or first-step snapshot. Momentum settings are written after
the accepted momentum solve, so the unresolved interval starts after that solve
and ends before the saved pressure result. It includes the outlet update,
pressure assembly, packing, Hypre setup/solve and residual checks. The absence of
a pressure settings file does not distinguish those operations.

With `--first-step-audit`, the executable now flushes fixed progress labels from
each rank. The saved-launch wrapper displays only complete, allowlisted marker
lines, for example `MARS rank 0: hypre_setup`. The private raw log, coordinates,
boundary names and numerical values are not copied to the terminal. Normal runs
without that audit option do not emit these markers.

| Last entered stage | Work entered before the next marker |
| --- | --- |
| `outlet_update` | Outlet update, boundary moments and exchange |
| `pressure_assembly` | Pressure matrix and boundary assembly |
| `pressure_pack` | Pack owned matrix rows and RHS |
| `pressure_backend` / `hypre_initialize` | Enter the pressure backend and initialize Hypre |
| `hypre_matrix` / `hypre_vectors` | Prepare Hypre matrix or vector storage |
| `hypre_amg_configure` | Create and configure the preconditioner and Krylov solver |
| `hypre_setup` | Hypre setup, including AMG hierarchy construction |
| `hypre_solve` | Hypre Krylov solve |
| `hypre_solution` / `hypre_residual` | Retrieve the candidate and check its residual |
| `pressure_exchange` / `pressure_residual` | Exchange the candidate and check the MARS residual |
| `pressure_audit` / `pressure_correction` | Save the pressure audit or apply the correction |
| `first_step_complete` / `final_output` | Finish step bookkeeping or write final fields |

These are host API entry markers, not GPU completion or convergence reports.
No device synchronization, MPI collective or change to the numerical method is
added. A missing later marker localizes the interval; it does not by itself prove
a deadlock inside that call. The public failure summary also records
`simple_progress_last_stage_by_rank`. `hypre_complete` and
`flow_output_complete` mark return from the wrapper and the end of flow output;
neither overrides a failed residual or certifies successful process exit.

Local tests check partial and overlong log lines, private text rejection, failed
child exit preservation, and the real Hypre setup/solve callback sequence.

### Setup timeout: bounded stack capture

The public `simple-reference-amg-progress-lFwAtd/launch-public.json` reports a
scheduler time limit and last stage `hypre_setup` on all four ranks. Library
checks pass, with no reported CUDA error or linear rejection. This narrows the
interval to setup, before the pressure solve. It does not prove a deadlock or
identify a particular AMG operation. The marker also precedes a possible wrapper
barrier and includes the post-setup collective check before the next marker.
GMRES setup allocates Krylov workspace before calling AMG; HMIS and coarse solver
9 then include CPU work. Increasing the Krylov iteration cap cannot help here.

Use the following diagnostic instead of repeating the quiet one-step run.
It reuses the existing executable, case, controls, binding and saved launcher.
Only the output prefix changes, and GDB is inserted between the binding script
and the executable. It does not modify saved evidence or run OpenAccel. No build
or repeat of the public GPU gates is needed.

The script first verifies GDB's Python support, child debugging, timed interruption
and readable stopped frames on `/bin/sleep` without requesting GPU resources.
Unreadable frames or API failures refuse the GPU launch. A reported unwind stop
after reading named frames permits sampling, with the incomplete self-test stack
retained in the public report and a visible warning. The preflight checks frame
reading, not a complete backtrace of `/bin/sleep`.
Under the saved launcher, each rank stops at the first Hypre GMRES/FlexGMRES setup
call. It samples after 30 seconds in setup, or at setup return if that occurs first. If setup is
not reached, it samples after 120 seconds of execution. Each rank has a five-minute
watchdog, including debugger initialization and waiting for the other snapshots;
the saved scheduler limit is unchanged. Queue time is outside this watchdog.

**This diagnostic deliberately stops the solver.** A complete diagnostic is not
a completed CFD run. Only exact, allowlisted library function names and fixed
status flags are exported; no frame arguments, locals, addresses, source paths,
geometry or field values are exported. Raw debugger and solver output stays in
the private directory. The sampler uses GDB's
[frame-name API](https://sourceware.org/gdb/current/onlinedocs/gdb.html/Frames-In-Python.html)
and a [finish breakpoint](https://sourceware.org/gdb/current/onlinedocs/gdb.html/Finish-Breakpoints-in-Python.html).
Missing symbols or an incomplete unwind can limit what the snapshot establishes.
It samples host stacks, not GPU kernels, and a single sample cannot prove deadlock.

The first user-run sample (`simple-setup-stack-4Z9YlY`) reached setup on all four
ranks and did not observe its return before interruption. Its empty function lists
do not locate an internal call. The original sampler used
`FRAME_UNWIND_FIRST_ERROR`, which the
[GDB 13.2 Python exports](https://gnu.googlesource.com/binutils-gdb/+/refs/tags/gdb-13.2-release/gdb/python/py-frame.c)
omit, and discarded collected names when a later API call failed. The release API
reproduces the empty-list symptom locally; the cluster debugger version was not
recorded, so this does not establish its exact failure cause.

The sampler now uses the exported normal stop reasons and preserves partial
frames. `stack_errors` contains only fixed API-stage labels, never exception text.
`unwind_stop_reasons` distinguishes GDB's reported unavailable registers/memory,
frame-ID failures, missing saved PC and memory errors using fixed labels. Unknown
reasons remain `unknown`; the obsolete normal end marker `UNWIND_NULL_ID` is not
treated as an error.
Unreadable main-thread stacks produce `stack_unreadable`; usable but incomplete
unwinds produce `diagnostic_partial`. Neither is `diagnostic_complete`. An observed
setup return can also locate the stop when no allowlisted caller remains on the
stack. This corrects the earlier success label, not the unresolved setup slowdown.

The next user-run preflight (`simple-setup-stack-pqs4ki`) successfully interrupted
the sleep process and read stopped frame names, but reported `unwind_stopped`.
No GPU job was launched. Requiring a complete unwind rejected this useful frame
reading capability. The corrected preflight permits that specific limitation;
it still rejects frame/API failures and keeps partial solver captures distinct
from complete ones. The earlier report does not identify the exact unwind reason.

### Interpolation communication: two samples from one launch

The user-run `simple-setup-stack-VMpVen/trace/public.json` has usable, complete
host stacks on all four ranks. Each is inside
`hypre_BoomerAMGBuildExtPIInterpDevice`; ranks 0 and 1 are in
`hypre_MPI_Testall`, rank 3 in `hypre_MPI_Test`. Setup has not returned.
This locates the sampled delay before the pressure solve. It does not establish
an indefinite deadlock or identify rank 2 as its cause.

In the inspected Hypre 2.33 source, calls to those two polling wrappers occur in
[`hypre_DataExchangeList`](https://github.com/hypre-space/hypre/blob/v2.33.0/src/utilities/exchange_data.c).
Its contact, response and termination buffers are on the host. Interpolation
reaches this protocol while constructing communication packages for derived
matrices; see
[`new_commpkg.c`](https://github.com/hypre-space/hypre/blob/v2.33.0/src/parcsr_mv/new_commpkg.c)
and [`par_lr_interp_device.c`](https://github.com/hypre-space/hypre/blob/v2.33.0/src/parcsr_ls/par_lr_interp_device.c).
The source and stacks support a host communication discovery interval, not a
diagnosis of GPU-buffer transport failure. MARS's preceding halo exchange is
completed before pressure setup, and uses distinct tags. These observations do
not rule out corrupted communication metadata or another earlier defect.

The first filter omitted the internal discovery helpers and `MPI_Iprobe`.
The extended literal allowlist exposes them without exporting arguments or data.
`--communication-progress` also watches entry and return of
`hypre_DataExchangeList`. After the first sample, all ranks save their stacks
before any resumes. A second sample follows 15 seconds after resuming each rank.
If any rank has returned from setup, exited or failed to stop cleanly, the group
does not resume. The five-minute watchdog and original scheduler limit remain.

The public `communication_progress` report contains first-sample function names
and fixed flags; `rank_stacks` contains the final sample. Exchange counters stay
private. With usable samples and complete tracking:

- `exchange_advanced_between_samples=true` means an exchange entered or returned
  between samples, even when both stacks show the same function names.
- `false` with both `exchange_active_*` flags true means the watched exchange
  did not return during the interval. It does not prove that MPI made no internal
  progress, that other work stopped, or that the wait is permanent.
- `null` means the diagnostic cannot establish progress. Missing return watches,
  unreadable frames or a peer timeout remain incomplete evidence.

GDB changes timing, so this is not a performance measurement. Entry breakpoints
stop briefly; return watches are created after their callbacks have returned,
as required by the
[GDB breakpoint API](https://sourceware.org/gdb/current/onlinedocs/gdb.html/Breakpoints-In-Python.html).
The numerical method, original targets, executable, libraries and launcher are
unchanged. No MARS rebuild, OpenAccel run or full nonlinear run is needed.
The following command supersedes the single-sample recipe for this delay.

Run in the MARS terminal only, without an outer `srun`:

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
repo="$scratch/git/mars-v010-check"
test "$(git -C "$repo" branch --show-current)" = cstone
git -C "$repo" fetch origin cstone
git -C "$repo" merge --ff-only refs/remotes/origin/cstone
python3 -c 'import numpy, netCDF4, yaml; print("Python dependencies OK")'
umask 077
mkdir -p "$scratch/tmp"
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
run=$(mktemp -d "$scratch/simple-setup-communication-XXXXXX")
printf 'Private diagnostic directory: %s\n' "$run"
status=0
python3 "$repo/scripts/simple_setup_stack.py" \
  --pair "$scratch/simple-reference-amg-progress-lFwAtd/pair" \
  --output-dir "$run/trace" --communication-progress || status=$?
cat "$run/trace/public.json"
printf 'Share only: %s\n' "$run/trace/public.json"
exit "$status"
)
```

Local synthetic tests exercise setup entry, return, interruption, MPI-rank capture
coordination, watchdog cleanup, output filtering and preserved saved arguments.
They also cover the GDB 13.2 API without the missing sentinel, retained partial
frames, unreadable self-test refusal and incomplete-capture status. The user-run
samples verified debugger launch, setup entry, timed interruption and usable
Hypre stack names. Synthetic communication checks cover unchanged versus
advancing exchanges with identical stacks, peer coordination, early setup return,
missing return watches, watchdog deadlines and privacy. The subsequent user-run
coordinated capture verified entry/return tracking and two usable samples on all
four ranks. This is diagnosis, not a numerical fix.

## Inspect the saved exchange sequence without a new run

The user-run coordinated capture found all four ranks still inside
`hypre_DataExchangeList`, reached through `hypre_ParCSRMatrixGenerateFFFCDevice`
and communication-package creation. No exchange entry or return occurred between
samples on any rank. This locates an unfinished host exchange while building
interpolation; it does not establish its cause or indefinite deadlock.

Compare the already saved private counters before another GPU allocation:

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
repo="$scratch/git/mars-v010-check"
test "$(git -C "$repo" branch --show-current)" = cstone
git -C "$repo" fetch origin cstone
git -C "$repo" merge --ff-only refs/remotes/origin/cstone
umask 077
export PYTHONDONTWRITEBYTECODE=1
run=$(mktemp -d "$scratch/simple-setup-sequence-XXXXXX")
status=0
python3 "$repo/scripts/simple_setup_stack.py" \
  --inspect-dir "$scratch/simple-setup-communication-GEDu2q/trace" \
  --output-dir "$run/check" || status=$?
cat "$run/check/public.json"
printf 'Share only: %s\n' "$run/check/public.json"
exit "$status"
)
```

This reads only the saved stack/counter records and their public summary. It
does not inspect fields, load the solver libraries, start GDB or launch a solver.
It needs no build, GPU allocation or scientific Python packages. The original
records are unchanged; counters and arbitrary source text stay private.

`same_observed_exchange_ordinal: false` means the tracked call sequences differ
across ranks. It is a reason to inspect call ordering, not proof that particular
MPI tags or communicators differ. A true value rules out this ordinal difference
only: neither matching ordinals nor matching stacks prove identical protocol
phases or communicators. Missing records, inconsistent summaries, incomplete
tracking, nested active exchanges or progress during the interval remain
invalid or inconclusive. No solver settings should be changed from this flag alone.

## Inspect the protocol run that stopped on signals

The `simple-setup-protocol-ox40BU/trace/public.json` report records unexpected
signals on ranks 2 and 3, stopped in `hypre_ParCSRMatrixRAPKTDevice`. Ranks 0 and
1 were in communication; their waits may follow the other ranks stopping.
This instrumented run does not prove the original stall's cause or establish
that the interposer caused the new stops. Do not rerun it or ignore the signals
before classifying the existing evidence.

`--inspect-dir` now recognizes a saved protocol capture. It checks the saved
stack reports against the public capture and exports each rank's trace state
and a fixed parser rejection label where needed. A valid trace ending between
exchanges is no longer described as a broken sequence. It also extracts only
allowlisted signal names from the private GDB logs. Raw log lines, addresses,
peer lists, tags, counters and message contents are not exported.

Run this block in one terminal. It uses only the Python standard library and
the already saved files: no rebuild, GDB process, library loading, solver launch
or GPU allocation. The original capture is left unchanged.

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
repo="$scratch/git/mars-v010-check"
test "$(git -C "$repo" branch --show-current)" = cstone
git -C "$repo" fetch origin cstone
git -C "$repo" merge --ff-only refs/remotes/origin/cstone
umask 077
export PYTHONDONTWRITEBYTECODE=1
run=$(mktemp -d "$scratch/simple-protocol-inspection-XXXXXX")
status=0
python3 "$repo/scripts/simple_setup_stack.py" \
  --inspect-dir "$scratch/simple-setup-protocol-ox40BU/trace" \
  --output-dir "$run/check" || status=$?
cat "$run/check/public.json"
printf 'Share only: %s\n' "$run/check/public.json"
exit "$status"
)
```

`inspection_complete` means the saved records were processed, not that setup
succeeded. `rank_signal_checks` retains missing logs and unknown signal names
explicitly. `protocol_checks.rank_trace_checks` distinguishes valid active,
between-exchange, pre-exchange and returned setup states from malformed data.
Signal names come from debugger text; they do not identify the faulting
instruction or establish the reason for the original uninstrumented wait.

## Identify the unfinished protocol stage

The offline inspection above supersedes another protocol launch for the
reported `ox40BU` result. Keep the following recipe for a future justified
capture after its unexpected stops have been explained.

The saved sequence inspection reports the same observed exchange ordinal on all
four ranks, with no entry or return between samples. It does not identify the
exchange's `rnum`, communicator or pending requests. Do not change the pressure
algorithm on that evidence alone.

The opt-in `--protocol-trace` records host message-call metadata during Hypre
setup. It distinguishes contact replies, child completion notices and parent
release, and checks whether matching sends were observed and returned. The
tracer forwards the original arguments to the already loaded, captured Hypre
library. Message buffers are never read or modified. Peer identifiers, tags,
request counters and raw traces stay private; only fixed checks are exported.

This diagnostic adds noncollective communicator queries and private file writes,
so it changes timing. A completed send does not prove receipt, and an absent
observed call is not proof of a lost message. Matching call contexts require the
same ordinal and `rnum`, with all traced communication on `MPI_COMM_WORLD`.
Untracked requests, unsupported contexts, missing hooks and incomplete records
remain inconclusive. This is not a performance measurement or deadlock proof.

Run this once in the **MARS terminal**. It compiles only a small host interposer
from the captured Hypre 2.33 headers and the existing build's MPI configuration;
MARS does not need rebuilding. The generated exchange declaration comes from
the installed header, without including CUDA internals in the host compiler.
The interposer's MPI libraries must match the saved case. An existing
`LD_PRELOAD` is refused. GDB itself and its sleep self-test are not preloaded.

The saved launcher, binding, executable and solver controls are preserved. The
existing timed setup sample, coordinated stop and five-minute debugger watchdog
bound the launch. No outer `srun`, OpenAccel rerun or repeat of public GPU gates
is needed.

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
repo="$scratch/git/mars-v010-check"
test "$(git -C "$repo" branch --show-current)" = cstone
git -C "$repo" fetch origin cstone
git -C "$repo" merge --ff-only refs/remotes/origin/cstone
python3 -c 'import numpy, netCDF4, yaml; print("Python dependencies OK")'
umask 077
mkdir -p "$scratch/tmp"
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
run=$(mktemp -d "$scratch/simple-setup-protocol-XXXXXX")
printf 'Private diagnostic directory: %s\n' "$run"
status=0
python3 "$repo/scripts/simple_setup_stack.py" \
  --pair "$scratch/simple-reference-amg-progress-lFwAtd/pair" \
  --output-dir "$run/trace" \
  --protocol-trace --build-cache "$repo/build-hypre/CMakeCache.txt" || status=$?
cat "$run/trace/public.json"
printf 'Share only: %s\n' "$run/trace/public.json"
exit "$status"
)
```

`protocol_checks.rank_waits` locates each rank's observed wait. For contact
replies, `matching_contacts_received` separates contacts that reached their
destination from replies whose send call was not observed. The send-call flags
describe all receives in the current group; they do not identify which
individual request remains incomplete. Numerical controls and message contents
are outside this diagnostic. If setup returns during sampling, the stack report
retains that fact; no unfinished protocol is claimed.

Local tests cover missing replies, unreceived contacts, unreturned sends,
different rounds/communicators, completed waits, truncation, privacy, build
refusal and inferior-only preload. The standalone public
`simple_performance/setup_protocol_check.cpp` fixture builds a tridiagonal host
matrix and runs HMIS/interpolation setup. The interposer was compiled with
warnings as errors and exercised with real CPU Hypre 2.33/OpenMPI on 1, 2 and 4
local ranks. Every setup completed and its exchange records passed parsing.
This does not validate the Linux preload or CUDA/Cray execution; the bounded
user-run diagnostic supplies that evidence.

## Earlier one-step progress recipe

The setup-stack diagnostic above supersedes this recipe for the reported timeout.
Keep this recipe for a normal first-step comparison once its cause is addressed.

Use the working MARS uenv/venv and existing PETSc-enabled build. Run the Git update
in this terminal only. Rebuild only SIMPLE for the new markers; the 48 passing
public GPU checks do not need repeating. The `run` command executes the saved
private launcher once with its original resource and time limits; do not wrap it
in another `srun`. OpenAccel is not relaunched. Progress labels appear while the
run is active; only the public JSON should be shared after a failure.

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
repo="$scratch/git/mars-v010-check"
test "$(git -C "$repo" branch --show-current)" = cstone
git -C "$repo" fetch origin cstone
git -C "$repo" merge --ff-only refs/remotes/origin/cstone
python3 -c 'import numpy, netCDF4, yaml; print("Python dependencies OK")'
umask 077
mkdir -p "$scratch/tmp"
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
cd "$repo/build-hypre"
cmake --build . --parallel 8 --target mars_segregated_simple
run=$(mktemp -d "$scratch/simple-reference-amg-progress-XXXXXX")
printf 'Private results: %s\n' "$run"

if ! python3 "$repo/scripts/simple_pressure_profile.py" prepare \
  --baseline-pair "$scratch/simple-first-step-MPQGzD/pair" \
  --defaults-probe "$scratch/simple-hypre-defaults-PrkNXY/probe" \
  --executable "$repo/build-hypre/examples/distributed/unstructured/mars_segregated_simple" \
  --output-dir "$run/pair" --both-linear-solvers --reference-amg \
  --output "$run/preparation-public.json"; then
  cat "$run/preparation-public.json"
  printf 'Share only: %s\n' "$run/preparation-public.json"
  exit 1
fi

if ! python3 "$repo/scripts/simple_pressure_profile.py" run \
  --pair "$run/pair" --output "$run/launch-public.json"; then
  cat "$run/launch-public.json"
  printf 'Share only: %s\n' "$run/launch-public.json"
  exit 1
fi
if ! python3 "$repo/scripts/simple_pressure_profile.py" compare \
  --pair "$run/pair" --detail-dir "$run/private-comparison" \
  --output "$run/comparison-public.json"; then
  cat "$run/comparison-public.json"
  printf 'Share only: %s\n' "$run/comparison-public.json"
  exit 1
fi

python3 - "$run/comparison-public.json" <<'PY'
import json, sys
with open(sys.argv[1]) as stream:
    report = json.load(stream)
keys = ("comparison_status", "failed_check", "pressure_amg_mode",
        "hypre_host_staging_permitted", "recorded_pressure_settings_verified",
        "original_pressure_targets_preserved", "original_momentum_targets_preserved",
        "gpu_adaptations", "identical_linear_solvers_verified")
summary = {key: report.get(key) for key in keys}
comparison = report.get("comparison", {})
summary["first_step"] = {key: comparison.get(key) for key in
    ("all_snapshots_match", "first_differing_iteration", "first_differing_fields",
     "first_differing_stage", "first_step_stage_matches", "pressure_common_system_checks")}
print(json.dumps(summary, indent=2))
PY
printf 'Share only: %s\n' "$run/comparison-public.json"
)
```

Agreement would show that retaining these two choices is sufficient for this
first-step comparison. If the pressure increment still differs with verified
settings and passing original residual targets, the substitutions are not the
complete explanation. In either case, keep the comparison bounded; completion
alone does not justify another long nonlinear run.
