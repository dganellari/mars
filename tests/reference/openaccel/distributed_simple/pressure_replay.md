# Replay the rejected pressure equation

## Inspect the completed recovery without another solve

The public `simple-pressure-recovery-fuUIJf/comparison-public.json` result verifies
both replay manifests and the common frozen system. The initial MARS candidate,
the recovered candidate and the reference candidate all fail the independent
residual target. Recovery exhausted its three-correction budget. In this code
path every correction was accepted only after the original residual interval
decreased; the summary alone does not measure how much it decreased. This result
does not justify resuming the full-flow history or claim an unattainable target.

`compare --recovery-progress` reads the saved, hash-bound `.recovery` traces on
the user's machine. It exports fixed bands for accepted residual reductions and
the remaining upper-bound distance to the original target. It also reports
whether a correction hit its cap, whether all reported correction residuals
met 0.1, and whether the final evaluation interval width is below the target.
These are correction-boundary observations, not a Krylov iteration history.
The independent final residual check remains the convergence decision. No raw
norms, iteration counts, matrix data or private paths are added to public JSON.

The exporter requires complete accepted-step traces (target reached or budget
exhausted), matching rank intervals, report/trace iteration counts and the
recorded stop reason. Malformed, incomplete or nonfinite traces fail. Exact
binary64 rational comparisons choose bands without overflow at tiny targets.

Run this **in the MARS terminal**. It does not rebuild, allocate GPUs or launch a
solver. It rechecks the existing result with the existing CPU residual checker.

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
repo="$scratch/git/mars-v010-check"
git -C "$repo" fetch origin cstone
git -C "$repo" merge --ff-only refs/remotes/origin/cstone
umask 077
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
mkdir -p "$TMPDIR"
reference="$scratch/simple-reference-replay-URu8yA/reference"
IFS= read -r archive < "$reference/input-archive-current.txt"
summary=$(mktemp -d "$scratch/simple-recovery-progress-XXXXXX")/public.json
status=0
python3 "$repo/scripts/simple_pressure_replay.py" compare \
  --capture-run "$scratch/simple-pressure-frozen-RlGQaI/capture" \
  --mars-run "$scratch/simple-pressure-recovery-fuUIJf/mars" \
  --reference-run "$reference" --reference-input-archive "$archive" \
  --checker "$repo/build-hypre/examples/distributed/unstructured/mars_simple_pressure_residual_check" \
  --recovery-progress --output "$summary" || status=$?
cat "$summary"
printf 'Share only: %s\n' "$summary"
exit "$status"
)
```

## Four-correction recovery on the frozen GPU system

The public `simple-recovery-progress-N9uuXM/public.json` report verifies that each
of the three corrections reduced the residual by at least 10 times. None reached
its iteration cap, and the final residual upper bound is between one and ten
times the original target. Its lower bound exceeds the target, so this remains
a definite failure, not an inconclusive check. The evaluation interval width is
smaller than the target; this alone does not establish an attainable error floor.

This measured progress supports allowing one more correction with
`--recovery-rounds 4`. Another tenfold reduction would suffice, but repeating
that reduction is not guaranteed. Three-correction runs remain supported and
verifiable. Production SIMPLE and the saved reference are unchanged.

For the same pressure-correction equation `A p' = b`, form the compensated defect
`r = b - A p'`, solve `A delta = r` from zero, and try `p' + delta`. Both unknowns
have pressure units. Boundary elimination, pressure reference, matrix, RHS,
partition and the original acceptance target are unchanged. This changes how
the linear system is solved; it is not an identical OpenAccel iteration path.

There are at most four additional solves, each with the original iteration cap,
relative tolerance 0.1, absolute tolerance zero and minimum iterations zero.
The initial solve has its own original budget. The same Krylov object and AMG
hierarchy are reused without another setup. This costs at most five times the
original Krylov iteration budget; it does not promise a speedup. A capped
correction can be used if the original residual decreases. Nonfinite values,
fatal backend errors and failure to establish a decrease stop recovery. Original
Krylov controls are restored afterward.

The CUDA path computes defects, updates and norms on device and uses Hypre's
device communication map. Recovery explicitly enables GPU-aware MPI; the runtime
must support device buffers. Only scalar reductions/control and diagnostic file
I/O use host data. The initial candidate is saved under `result/initial`; private
per-rank `.recovery` files record bounds and correction exit metadata. The initial
`result_*` backend fields retain their original meaning and are labelled as such
in public comparisons. The new files are covered by the output manifest.

The internal decrease test includes compensated row-error bounds and a
conservative FP64 norm-reduction margin. It rejects subnormal nonzero RHS squares
and stops if residual squares underflow; it does not certify those scales. With
zero RHS and zero absolute tolerance, the evaluation bound can leave an exact
zero candidate inconclusive. The independent checker of the original captured
rows is the final authority, regardless of the internal stopping reason.

Run this **once in the MARS terminal**. It uses the saved rank count and launcher;
no new OpenAccel run or full-flow run is launched. Detailed output stays private.

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
repo="$scratch/git/mars-v010-check"
test "$(git -C "$repo" branch --show-current)" = cstone
git -C "$repo" fetch origin cstone
git -C "$repo" merge --ff-only refs/remotes/origin/cstone
python3 -c 'import numpy, netCDF4, yaml'
umask 077
export TMPDIR="$scratch/tmp" XDG_CACHE_HOME="$scratch/.cache"
export PYTHONDONTWRITEBYTECODE=1 MPICH_GPU_SUPPORT_ENABLED=1
mkdir -p "$TMPDIR" "$XDG_CACHE_HOME"
cmake --build "$repo/build-hypre" --parallel 4 --target \
  mars_simple_pressure_replay mars_simple_pressure_residual_check

capture="$scratch/simple-pressure-frozen-RlGQaI/capture"
profile="$scratch/simple-original-pressure-1xGpr2/pair"
reference="$scratch/simple-reference-replay-URu8yA/reference"
IFS= read -r archive < "$reference/input-archive-current.txt"
test -f "$archive/archive.json"
run=$(mktemp -d "$scratch/simple-pressure-recovery4-XXXXXX")
printf 'Private results: %s\n' "$run"
status=0
python3 "$repo/scripts/simple_pressure_replay.py" replay \
  --capture-run "$capture" --backend mars --profile gpu-reference \
  --gpu-profile-pair "$profile" --recovery-rounds 4 \
  --executable "$repo/build-hypre/examples/distributed/unstructured/mars_simple_pressure_replay" \
  --output-dir "$run/mars" || status=$?
cat "$run/mars/public.json"
if (( status != 0 )); then exit "$status"; fi
python3 "$repo/scripts/simple_pressure_replay.py" compare \
  --capture-run "$capture" --mars-run "$run/mars" \
  --reference-run "$reference" --reference-input-archive "$archive" \
  --checker "$repo/build-hypre/examples/distributed/unstructured/mars_simple_pressure_residual_check" \
  --output "$run/comparison-public.json" || status=$?
cat "$run/comparison-public.json"
printf 'Share only: %s\n' "$run/comparison-public.json"
exit "$status"
)
```

Success requires `residual_checks.mars.residual_passed=true`. The comparison also
reports `initial_residual_checks.mars` and `recovery_checks.mars`, so improvement
cannot erase the original failure. `stopping_checks.mars` describes the initial
solve, not the corrected candidate. `replay_complete` alone is not convergence.
If the final independent check fails or is inconclusive, do not resume the long
flow history. A pass permits evaluating this recovery in SIMPLE; it does not
establish pump convergence or OpenAccel field parity.

Local checks cover recovery after a capped solve, budget exhaustion at a stricter
target, already-converged, zero-RHS and inconsistent systems, original/final
metadata separation and tampered evidence. Real CPU Hypre 2.33/3.1 builds run
under ASan/UBSan. CPU MPI Hypre 2.32 exercises 1/2/4 ranks with separate owned rows
and exchanged off-rank values. The four-correction regression includes a target
that fails after three corrections and passes after four in the sequential
builds, plus distributed success and budget-exhaustion cases. Correction counts
can depend on the partition. All 70 replay tests pass with these local builds;
this does not replace the CUDA build/run above.
The optional local test executables are selected with
`MARS_TEST_PRESSURE_REPLAYS` (colon-separated), `MARS_TEST_PRESSURE_MPI_REPLAY`
and `MARS_TEST_PRESSURE_CHECKER`; run `scripts/test_simple_pressure_replay.py`
with `PYTHONPATH=scripts`. `MPIEXEC` selects the local MPI launcher.

## Replay the failing system with the verified GPU profile

The original-pressure first-step comparison found matching unscaled pressure
matrix/RHS and distinct increments that both meet the reference target. Earlier
tighter first-step results improved field agreement, but the tighter history
failed a pressure solve. This experiment applies the verified GPU-compatible
profile to that saved failure; it does not launch another flow run.

`--profile gpu-reference --gpu-profile-pair PAIR` verifies the saved profile,
its first-step runtime settings, baseline/reference identity and library hashes.
It copies the profile's algorithm controls into a new settings file and replaces
only its rtol/atol with the frozen system's captured target. The original profile,
failure capture and reference replay remain unchanged. Rank count, MARS libraries
and the reference's explicit algorithm controls must match the capture.

For the frozen pressure equation A p' = b, p' remains in pressure units. The CSR,
RHS, ownership and zero initial guess are unchanged. Acceptance remains
max(atol, rtol ||b||), evaluated by the independent compensated residual checker
with its evaluation bound. The replay uses the production profile setter order,
including cycle-specific smoothers, and checks actual settings after AMG setup.
AMG builds a fresh hierarchy; this is not a replay of the original hierarchy or
an identical CPU/GPU solver. Native GPU SpMV stays selected. No correction solves
or relaxed acceptance thresholds are added.

Run this **once in the MARS terminal**, with its working uenv/Python environment.
It rebuilds only the replay and checker, reuses the captured launcher, and compares
against the archived reference replay. No new OpenAccel job is needed. Keep the
source profile pair available for later identity checks, even if runtime
libraries have been archived. All detailed files remain private on scratch.

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
repo="$scratch/git/mars-v010-check"
test "$(git -C "$repo" branch --show-current)" = cstone
git -C "$repo" fetch origin cstone
git -C "$repo" merge --ff-only refs/remotes/origin/cstone
python3 -c 'import numpy, netCDF4, yaml'
umask 077
mkdir -p "$scratch/tmp"
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
cmake -S "$repo" -B "$repo/build-hypre" \
  -DMARS_ENABLE_CUDA=ON -DMARS_ENABLE_MPI=ON \
  -DMARS_ENABLE_UNSTRUCTURED=ON -DMARS_ENABLE_HYPRE=ON \
  -DMARS_ENABLE_FEM_EXAMPLES=ON -DMARS_ENABLE_SEGREGATED=ON
cmake --build "$repo/build-hypre" --parallel 4 --target \
  mars_simple_pressure_replay mars_simple_pressure_residual_check

capture="$scratch/simple-pressure-frozen-RlGQaI/capture"
profile="$scratch/simple-original-pressure-1xGpr2/pair"
reference="$scratch/simple-reference-replay-URu8yA/reference"
IFS= read -r archive < "$reference/input-archive-current.txt"
test -f "$archive/archive.json"
run=$(mktemp -d "$scratch/simple-pressure-gpu-profile-XXXXXX")
printf 'Private results: %s\n' "$run"
status=0
python3 "$repo/scripts/simple_pressure_replay.py" replay \
  --capture-run "$capture" --backend mars --profile gpu-reference \
  --gpu-profile-pair "$profile" \
  --executable "$repo/build-hypre/examples/distributed/unstructured/mars_simple_pressure_replay" \
  --output-dir "$run/mars" || status=$?
cat "$run/mars/public.json"
if (( status != 0 )); then exit "$status"; fi

python3 "$repo/scripts/simple_pressure_replay.py" compare \
  --capture-run "$capture" --mars-run "$run/mars" \
  --reference-run "$reference" --reference-input-archive "$archive" \
  --checker "$repo/build-hypre/examples/distributed/unstructured/mars_simple_pressure_residual_check" \
  --output "$run/comparison-public.json" || status=$?
cat "$run/comparison-public.json"
printf 'Share only: %s\n' "$run/comparison-public.json"
exit "$status"
)
```

`replay_complete` verifies execution and applied settings, not convergence.
In the comparison, inspect `residual_checks.mars.residual_passed` and
`stopping_checks.mars`: a failed or inconclusive independent check does not pass
because Hypre reports convergence. The old reference can remain a recorded
failure. A new MARS pass justifies a tighter-history test; it does not establish
first-step field parity or nonlinear pump convergence. If it still fails, retain
this frozen comparison for further solver analysis rather than repeat a long run.

Local validation covers synthetic profile provenance, target preservation,
wrong cycle settings, independent rejection of bad solutions, and real CPU
Hypre 2.33/3.1 GMRES/FlexGMRES with one and multiple AMG levels. CPU profile tests
exercise configuration and residual checks; CUDA compilation and execution remain
a user-run gate.


This is an opt-in diagnostic, not a new pressure algorithm. It freezes the first
rejected pressure equation `A phi = b` before the existing error stops SIMPLE.
It saves owned CSR rows, the local-to-global solver map, RHS, halo-complete
candidate, tolerance and live Krylov/AMG controls. Nothing is captured on success,
on a momentum rejection, or when the linear backend cannot supply a candidate.
Experimental pressure refinement must remain off.

The capture contains **private mesh-derived data**, even though it has no mesh
coordinates. Matrices, vectors, settings, logs and detailed reports must remain
on the cluster. The option explicitly enables device-to-host copies for this
diagnostic file output only. Normal solver iterations have no added field copies.
Directories must be new; parts have mode 0600 and directories mode 0700.

The replay starts from zero, with the identical captured rows, RHS, rank count,
ownership and stopping target. There are two configurations:

* MARS: its recorded Krylov and AMG controls, native GPU SpMV and GPU relaxation18.
* Reference: the saved OpenAccel pressure configuration, interpreted using the
  inspected `ContextHYPRE` calls and the reference's captured CPU Hypre library.
  Omitted controls remain omitted, so that library supplies its defaults.
  OpenAccel forces AMG to one cycle with tolerance zero and does not forward its
  generic minimum-iteration option; the replay does the same.

Both create a **fresh** hierarchy. They do not reproduce the original cached
hierarchy or OpenAccel's original mesh partition. This isolates configuration
and library behavior on one MARS system; it is not a complete application replay.
The optional reference `--profile captured` instead applies the recorded MARS
controls to the CPU library. That is a separate controlled comparison, not the
OpenAccel profile.

The independent checker evaluates the original owned rows using compensated
products/sums and their existing Dot2Err evaluation bound. It encloses the global
residual and RHS norms with outward rounding. A pass requires the residual upper
bound to meet the lower bound of `max(atol, rtol*||b||)`; a failure requires the
residual lower bound to exceed the upper limit. Overlap is inconclusive, and
nonfinite values never pass. This certifies the saved binary64 algebra under the
checker’s floating-point assumptions, not the PDE error or nonlinear convergence.

Backend return/error/converged values and effective settings stay in the private
report; only fixed booleans and profile labels enter the public comparison.
The code does not instrument Hypre's actual early-exit branch. A convergence flag
with a failed common residual shows false acceptance, not which branch caused it.

## MARS terminal: build and capture

Keep the working MARS uenv and Python environment. The history pointer must name
the previously prepared twenty-step pair, whose reference capture completed and
whose MARS run reached the pressure rejection. No new deck or mesh is prepared.
This block enables the segregated targets in the existing build and rebuilds
them, then reuses that pair's launcher,
controls and restored solver environment. It disables field/snapshot output and
adds the private failure capture. It preserves the original iteration count.

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
repo="$scratch/git/mars-v010-check"
test "$(git -C "$repo" branch --show-current)" = cstone
git -C "$repo" fetch origin cstone
git -C "$repo" merge --ff-only refs/remotes/origin/cstone
cmake -S "$repo" -B "$repo/build-hypre" \
  -DMARS_ENABLE_CUDA=ON -DMARS_ENABLE_MPI=ON \
  -DMARS_ENABLE_UNSTRUCTURED=ON -DMARS_ENABLE_HYPRE=ON \
  -DMARS_ENABLE_FEM_EXAMPLES=ON -DMARS_ENABLE_SEGREGATED=ON
cmake --build "$repo/build-hypre" --parallel 4 --target \
  mars_segregated_simple mars_simple_pressure_replay mars_simple_pressure_residual_check
python3 -c 'import numpy, netCDF4, yaml'
umask 077
mkdir -p "$scratch/tmp"
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
run=$(mktemp -d "$scratch/simple-pressure-frozen-XXXXXX")
printf '%s\n' "$run" > "$scratch/simple-pressure-frozen-current.txt"
printf 'Private results: %s\n' "$run"
pair=$(cat "$scratch/simple-pressure-history-current.txt")
status=0
python3 "$repo/scripts/simple_pressure_replay.py" capture \
  --pair "$pair" \
  --executable "$repo/build-hypre/examples/distributed/unstructured/mars_segregated_simple" \
  --output-dir "$run/capture" || status=$?
cat "$run/capture/public.json"
exit "$status"
)
```

`capture_complete` means the rejected equation was saved and the original failure
was retained. It does **not** mean that SIMPLE passed. A missing capture, a signal
before completion, an unsupported profile, or changed inputs fails the diagnostic.
Vendor SpMV captures are rejected: this experiment requires the native path.

## MARS terminal: replay that equation on the GPU

This uses the captured MARS launcher, including its rank count and GPU binding.
It performs one linear solve, not another flow history.

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
repo="$scratch/git/mars-v010-check"
umask 077
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
run=$(cat "$scratch/simple-pressure-frozen-current.txt")
status=0
python3 "$repo/scripts/simple_pressure_replay.py" replay \
  --capture-run "$run/capture" --backend mars \
  --executable "$repo/build-hypre/examples/distributed/unstructured/mars_simple_pressure_replay" \
  --output-dir "$run/mars" || status=$?
cat "$run/mars/public.json"
exit "$status"
)
```

## OpenAccel terminal: replay with the reference library and controls

Use the working OpenAccel environment. Compilation uses its recorded plain
compiler and explicit MPI flags from `prgenv/CMakeCache.txt`. Installed internal
Hypre headers must match the library. The launcher retains the captured MARS
rank count/partition. Each rank records its actual loaded Hypre and MPI paths;
their hashes must match the corresponding application capture. Dependencies and
input/output hashes are checked before comparison. The full launcher child
environment and original application effective AMG hierarchy are not certified.
This automated identity check requires shared Hypre and MPI libraries; static
builds are rejected rather than silently losing the library binding.

The reference probe is compiled with `-fPIC`. Without it, a function pointer can
refer to an executable PLT stub, making `dladdr` name the executable instead of
the shared library ([documented limitation](https://man7.org/linux/man-pages/man3/dladdr.3.html)).
This changes library identification, not solver settings or residual acceptance.

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
repo="$scratch/git/mars-v010-check"
umask 077
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
run=$(cat "$scratch/simple-pressure-frozen-current.txt")
status=0
python3 "$repo/scripts/simple_pressure_replay.py" replay \
  --capture-run "$run/capture" --backend reference \
  --build-cache "$scratch/git/OpenAccel/prgenv/CMakeCache.txt" \
  --output-dir "$run/reference" || status=$?
cat "$run/reference/public.json"
exit "$status"
)
```

## Inspect an existing failed replay

Older summaries use `replay_launch` for both execution failures and checks after
execution. New runs distinguish `launcher_exit`, `completion_marker` and
`loaded_library_identity`. Do not rerun a flow capture to diagnose that label.

In the same terminal environment that launched the replay, `inspect` reads its
saved log, exit status and metadata. It checks input hashes and per-rank library
identities and reports which output parts exist. It neither launches a program
nor evaluates a residual. `inspection_complete` only means the inspection ran;
`failed_check` identifies the failed check. Unknown messages and private paths,
rank counts and solver values are never copied into its public output.

`library_identity_checks` gives separate Hypre and MPI verdicts. Fixed failure
labels distinguish a wrong hash, an executable recorded instead of a library,
an unreadable file, a relative path, and a missing or malformed identity record.
`executable_instead_of_library` is consistent with the PLT issue; it does not
establish which shared library that old process actually used. The old result
remains unverified. Neither matching `ldd` output nor an exit code of zero
overrides the per-rank hash check.

For the saved reference attempt:

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
repo="$scratch/git/mars-v010-check"
test "$(git -C "$repo" branch --show-current)" = cstone
git -C "$repo" fetch origin cstone
git -C "$repo" merge --ff-only refs/remotes/origin/cstone
umask 077
mkdir -p "$scratch/tmp"
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
run=$(cat "$scratch/simple-pressure-frozen-current.txt")
summary=$(mktemp -d "$scratch/simple-replay-inspection-XXXXXX")/public.json
status=0
python3 "$repo/scripts/simple_pressure_replay.py" inspect \
  --replay-run "$run/reference" --output "$summary" || status=$?
cat "$summary"
printf 'Share only: %s\n' "$summary"
exit "$status"
)
```

## Retry an executable-address identity failure

Run this block only in the OpenAccel terminal. It first inspects the old result
without launching anything. It retries only if the sole identity issue is an
executable recorded instead of a library. Other failures stop with the public
inspection report. The retry compiles the small reference replay with `-fPIC`
and solves the already captured equation once; no MARS build, new flow capture
or full pump run is needed. Original results remain untouched.

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
repo="$scratch/git/mars-v010-check"
test "$(git -C "$repo" branch --show-current)" = cstone
git -C "$repo" fetch origin cstone
git -C "$repo" merge --ff-only refs/remotes/origin/cstone
umask 077
mkdir -p "$scratch/tmp"
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
run=$(cat "$scratch/simple-pressure-frozen-current.txt")
retry=$(mktemp -d "$scratch/simple-reference-replay-XXXXXX")
python3 "$repo/scripts/simple_pressure_replay.py" inspect \
  --replay-run "$run/reference" --output "$retry/previous-public.json"
cat "$retry/previous-public.json"
python3 - "$retry/previous-public.json" <<'PY'
import json, sys
with open(sys.argv[1]) as stream:
    report = json.load(stream)
issues = {reason for check in report.get('library_identity_checks', {}).values()
          for reason in check['failures']}
if (report['failed_check'] != 'loaded_library_identity'
        or not report['launch_inputs_unchanged']
        or issues != {'executable_instead_of_library'}):
    sys.exit('STOP: saved evidence does not identify the executable-address issue. Share this public JSON.')
PY
status=0
python3 "$repo/scripts/simple_pressure_replay.py" replay \
  --capture-run "$run/capture" --backend reference \
  --build-cache "$scratch/git/OpenAccel/prgenv/CMakeCache.txt" \
  --output-dir "$retry/reference" || status=$?
cat "$retry/reference/public.json"
if test "$status" -eq 0; then
  printf '%s\n' "$retry/reference" > "$run/reference-retry-current.txt"
fi
printf 'Share only: %s\n' "$retry/reference/public.json"
exit "$status"
)
```

The comparison below uses the verified retry when its pointer exists. Its library
identities apply only to the new process; the old result remains unverified.
The common residual check still decides numerical acceptance.

## Separate uenvs: archive the reference dependencies once

If comparison reports missing reference MPI libraries, other runtime libraries
or build inputs, those paths may belong to the OpenAccel uenv. Run the following
block **in the OpenAccel terminal**, with that environment active. It reads the
existing successful replay, checks all original input and output hashes and the
recorded per-rank loaded-library identities, then copies the executable,
libraries and source/build inputs into a new private archive under scratch.
Each copied file must match its recorded launch hash. No build, allocation or
solver is launched, and saved results are not rewritten.

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
run=$(cat "$scratch/simple-pressure-frozen-current.txt")
reference="$run/reference"
if test -f "$run/reference-retry-current.txt"; then
  IFS= read -r reference < "$run/reference-retry-current.txt"
fi
archive=$(mktemp -d "$scratch/simple-replay-inputs-XXXXXX")/inputs
status=0
python3 "$repo/scripts/simple_pressure_replay.py" archive-inputs \
  --capture-run "$run/capture" --replay-run "$reference" \
  --output-dir "$archive" || status=$?
cat "$archive/public.json"
if test "$status" -eq 0; then
  printf '%s\n' "$archive" > "$reference/input-archive-current.txt"
fi
printf 'Share only: %s\n' "$archive/public.json"
exit "$status"
)
```

The archive is bound to the exact replay record and capture. Comparison hashes
the archived bytes against the original input manifest before and after residual
evaluation; it still checks the live capture and output files. An archive from
another replay, missing or changed archived bytes, changed results, and a public
summary supplied in place of an archive all fail. Source dependencies need only
be visible when archiving. The comparison labels this scope explicitly as
`archived_dependencies_and_live_capture_inputs`; it does not claim those original
uenv paths remain unchanged or mounted at comparison time. Without archive
arguments, live input checks remain the default.

This uses the same local manifest trust model as replay: hashes detect changes
relative to the saved records, not coordinated forgery of records and files.
It is not a signed execution attestation or a convergence certificate. Keep the
whole archive private; filenames and manifests can identify the private run.
Only its fixed `public.json` is shareable.

## MARS terminal: compare all three candidates

The host checker reads only the frozen diagnostic files. Production assembly,
solving and halo exchange remain on the GPU. Share only the final public JSON.
The comparison also rechecks both replays' recorded inputs and outputs. Use the
archive step above if the reference runtime and build inputs are not visible in
this terminal. The block uses that archive when its pointer exists; otherwise
both replays' original inputs must be visible. `--mars-input-archive` supports
the same mechanism for a MARS replay archived in its own environment.

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
repo="$scratch/git/mars-v010-check"
umask 077
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
run=$(cat "$scratch/simple-pressure-frozen-current.txt")
reference="$run/reference"
if test -f "$run/reference-retry-current.txt"; then
  IFS= read -r reference < "$run/reference-retry-current.txt"
fi
archive_args=()
if test -f "$reference/input-archive-current.txt"; then
  IFS= read -r archive < "$reference/input-archive-current.txt"
  archive_args=(--reference-input-archive "$archive")
fi
summary=$(mktemp -d "$scratch/simple-pressure-comparison-XXXXXX")/public.json
status=0
python3 "$repo/scripts/simple_pressure_replay.py" compare \
  --capture-run "$run/capture" --mars-run "$run/mars" --reference-run "$reference" \
  --checker "$repo/build-hypre/examples/distributed/unstructured/mars_simple_pressure_residual_check" \
  "${archive_args[@]}" --output "$summary" || status=$?
cat "$summary"
printf 'Share only: %s\n' "$summary"
exit "$status"
)
```

`completed` says the evidence was checked, not that either solver passed. Read
`residual_checks` for the original rejected candidate and the two new candidates.
If the MARS replay succeeds where the original failed, fresh setup/replay did not
reproduce the failure; do not attribute that difference to OpenAccel controls.
If both fail the common test, changing to the reference profile alone has not
resolved the accuracy problem. An inconclusive interval requires better residual
resolution, not tolerance relaxation.

### Read the stopping behavior from the same saved results

The completed comparison also exports `stopping_checks` for MARS and reference.
These use the iteration counts, effective limit, method, error flags and reported
relative residual already saved by the replay executable. No rebuild or solver
run is needed. All numeric values remain private; only fixed labels and booleans
are exported. Missing or malformed stopping fields fail the comparison at
`<candidate>_replay_stopping_report`.

* `iteration_limit_relation`: below, at, above, or mixed across ranks.
* `assessment`: distinguishes a definite residual failure before the limit,
  at/above the limit, or with zero iterations. Fatal errors, nonfinite results,
  disagreement between ranks, passed residuals and inconclusive residuals have
  separate labels.
* `convergence_claim_contradicted`: at least one backend convergence flag is set,
  but the independent bounded residual definitely fails.
* `solve_return_nonzero` and `global_error_nonzero`: retain nonconvergence errors
  even when `fatal_backend_error_seen` is false.
* `reported_relative_residual_below_rtol`: compares Hypre's reported value with
  its recorded relative tolerance only. This is not the full stopping test;
  `absolute_tolerance_enabled` flags a potentially larger absolute allowance.
  Independent residual acceptance remains authoritative.

A failure below the limit shows that the iteration cap did not end that solve;
it does not prove stagnation, residual drift or an unattainable tolerance. A
failure at the limit does not guarantee that more iterations will help.
`actual_exit_branch_verified` remains false. The original failed flow candidate
has no replay stopping report; these labels describe the two fresh solves only.

To update an existing successful comparison, run the block below in the **MARS
terminal only**. It reuses the reference dependency archive if present. Original
captures, replays and reports remain untouched. Share its new public JSON.

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
repo="$scratch/git/mars-v010-check"
test "$(git -C "$repo" branch --show-current)" = cstone
git -C "$repo" fetch origin cstone
git -C "$repo" merge --ff-only refs/remotes/origin/cstone
umask 077
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
mkdir -p "$TMPDIR"
run=$(cat "$scratch/simple-pressure-frozen-current.txt")
reference="$run/reference"
if test -f "$run/reference-retry-current.txt"; then
  IFS= read -r reference < "$run/reference-retry-current.txt"
fi
archive_args=()
if test -f "$reference/input-archive-current.txt"; then
  IFS= read -r archive < "$reference/input-archive-current.txt"
  archive_args=(--reference-input-archive "$archive")
fi
summary=$(mktemp -d "$scratch/simple-pressure-stopping-XXXXXX")/public.json
status=0
python3 "$repo/scripts/simple_pressure_replay.py" compare \
  --capture-run "$run/capture" --mars-run "$run/mars" --reference-run "$reference" \
  --checker "$repo/build-hypre/examples/distributed/unstructured/mars_simple_pressure_residual_check" \
  "${archive_args[@]}" --output "$summary" || status=$?
cat "$summary"
printf 'Share only: %s\n' "$summary"
exit "$status"
)
```

### If comparison stops before a verdict

Older versions used `replay_identity` for every failure after locating the
checker, including checker execution and report parsing. That label alone does
not identify a stale replay or a numerical failure.

The current comparison checks both replay records before running the residual
checker. `replay_evidence_checks` distinguishes missing, changed and unreadable
files using fixed categories: executable, capture input, Hypre library, MPI
library, other runtime library, source/build input and replay output. No paths,
hashes, counts or private values are exported. A missing library may be outside
the current uenv mount; it does not establish that the library was deleted.

`failed_candidate` identifies the original, MARS or reference candidate.
`failed_check` distinguishes record/binding/input/output checks, checker launch,
checker exit, checker output parsing and replay report parsing. Checker failures
also report a fixed launch reason or exit code and fixed loader/error flags.
All existing input/output hashes remain required, including checks after
evaluation. File problems are never converted into numerical passes.

For an old `replay_identity` result, fetch the updated script in one terminal
only, then repeat the comparison block with the same capture/replays and a fresh
summary. This is saved-data checking only: no solver launch, capture, build or
allocation. If inputs belong to the other uenv, archive them in that environment
using the step above. If they are unavailable there too, restore the exact
recorded files; do not alter their hashes to bypass the check.

## Local verification

The tests use synthetic matrices only. Run
`python3 -m unittest discover -s scripts -p test_simple_pressure_replay.py` for
format, provenance and independent-residual checks. It compiles the checker into
`test-scratch/`. Set `MARS_TEST_PRESSURE_REPLAYS` to a colon-separated list of
locally built CPU replay executables to exercise real GMRES/FlexGMRES solves and
settings round trips. The existing `marsSimpleLinearRejection{1,2,4}` host gates
also check capture timing, unchanged rows/RHS and preservation of rejection.

Alps CUDA compilation, GPU replay and private pressure results remain separate
validation gates; host tests do not establish them.
