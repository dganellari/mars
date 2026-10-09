# Replay the rejected pressure equation

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
