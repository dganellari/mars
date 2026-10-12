# Twenty steps at the original pressure targets

The retained-pressure and early-recovery experiments followed the tightened
pressure-accuracy history. The early-budget GPU comparison rejected the faster
candidate because it missed that experiment's original captured target. It is
not an accepted optimization or an explanation of nonlinear convergence.

This comparison returns to the original `MPQGzD` deck and its pressure targets.
It runs both solvers for up to twenty iterations from declared zero fields and
compares states zero through twenty. The existing original first-step capture
and all recovery results remain intact. A new reference run is required: the
saved one-step outputs cannot be used as a twenty-step history.

## Controls and evidence

`simple_pressure_profile.py prepare --history` derives a new private deck from
the verified original one-step pair and resolves the same GPU profile as the
original-pressure comparison. Tightened or already profiled baselines fail.
Only the iteration/output controls and matrix-dump switches change in the
reference deck. Explicit pressure targets, momentum settings, boundary
conditions, density, viscosity, pseudo-time and relaxation are preserved.

The nodal CVFEM pressure increment still solves `A phi = b` in pressure units,
with `max(atol, rtol*||b||)` as its original acceptance target. No assembly,
gradient, flux, velocity update or stopping tolerance is changed. MARS uses
ordinary FP64 pressure with `--pressure-expansion 0 --pressure-refinement 0`.
Its existing independent runtime residual checks still reject failed solves.

As in the one-step GPU profile, CPU HMIS is replaced by GPU PMIS and coarse
direct relaxation by GPU l1-Jacobi. The library versions also differ. Identical
linear algorithms or bitwise iteration histories are not claimed.

The driver permits the pressure profile with full distributed snapshots. After
each accepted pressure solve it writes that rank's actual settings to a private
`flow-pressure-step-N-rank-R.settings` file. It writes no matrix audit in this
mode. The comparator requires the complete twenty-step/rank coverage and checks
every record against the requested settings. The ordinary production default
and expansion paths remain unchanged.

Both fresh launches must match their respective saved executable/library,
launcher, rank and numerical-environment records; only rebuilding MARS is
allowed. Known MARS scalar overrides are restored in its child process. The
OpenAccel run requires its own saved runtime and executable, including PETSc
support where the deck needs it. No executable is guessed from the shell's PATH.

The comparison checks matching initial fields, node mapping, complete outputs,
all recorded controls and disabled recovery. It reports velocity and pressure
max/RMS error bands at each state, using the existing `U` and `rho U^2` scales
and `1e-5` field threshold. Exact errors, fields, settings and paths stay private.
`twenty_step_fields_match` proves only this bounded field comparison. A mismatch
does not by itself prove a discretization error; different finite-accuracy
linear solutions can produce different trajectories. The history does not save
the matrices needed to recompute every pressure residual independently, and
does not certify OpenAccel's backend convergence flags.

## Run both blocks, in order

Use the working MARS uenv/venv for the first block. It updates/builds only MARS,
checks the CLI, prepares the new pair, then runs MARS once. If it reports a
failure, share its public JSON and stop; do not run the second block yet.
The wrappers reuse each recorded complete `srun` command, including rank count
and GPU/NUMA binding. Do not wrap them in another `srun`.

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
repo="$scratch/git/mars-v010-check"
test "$(git -C "$repo" branch --show-current)" = cstone
git -C "$repo" fetch origin cstone
git -C "$repo" merge --ff-only origin/cstone
python3 -c 'import numpy, netCDF4, yaml'
cmake --build "$repo/build-hypre" --parallel 4 \
  --target mars_segregated_simple mars_simple_options_test
"$repo/build-hypre/mars_simple_options_test"
umask 077
mkdir -p "$scratch/tmp"
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
run=$(mktemp -d "$scratch/simple-original-history-XXXXXX")
probe="$repo/scripts/simple_pressure_profile.py"
python3 "$probe" prepare --history \
  --baseline-pair "$scratch/simple-first-step-MPQGzD/pair" \
  --defaults-probe "$scratch/simple-hypre-defaults-PrkNXY/probe" \
  --executable "$repo/build-hypre/examples/distributed/unstructured/mars_segregated_simple" \
  --output-dir "$run/pair"
printf '%s\n' "$run" > "$scratch/simple-original-history-current.txt"
printf 'Private results: %s\n' "$run"
status=0
python3 "$probe" run --pair "$run/pair" --solver mars \
  --output "$run/mars-public.json" || status=$?
cat "$run/mars-public.json"
printf 'Share only: %s\n' "$run/mars-public.json"
exit "$status"
)
```

After `capture_complete`, run this in the **OpenAccel terminal**, with its
working environment and Python dependencies. It runs the saved reference
executable once and compares the two histories. Do not rebuild OpenAccel or pull
the same checkout concurrently from both terminals.

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
repo="$scratch/git/mars-v010-check"
python3 -c 'import numpy, netCDF4, yaml'
umask 077
mkdir -p "$scratch/tmp"
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
IFS= read -r run < "$scratch/simple-original-history-current.txt"
probe="$repo/scripts/simple_pressure_profile.py"
status=0
python3 "$probe" run --pair "$run/pair" --solver openaccel \
  --output "$run/reference-public.json" || status=$?
cat "$run/reference-public.json"
if test "$status" -ne 0; then exit "$status"; fi
python3 "$probe" compare --pair "$run/pair" \
  --detail-dir "$run/private-comparison" \
  --output "$run/comparison-public.json" || status=$?
cat "$run/comparison-public.json"
printf 'Share only: %s\n' "$run/comparison-public.json"
exit "$status"
)
```

The complete direct `srun` line is saved locally in each of
`pair/mars-launch-command.txt` and `pair/openaccel-launch-command.txt`, immediately
before its launch. These contain private case arguments; inspect them locally
only. They document the exact command, while the wrapper additionally supplies
the checked environment and records. Do not rerun those lines in a used output
directory or share them.

`capture_complete` alone is not a field pass. Exit 2 with the required iteration
limit is a valid MARS capture. A rejected solve or missing step is a failed
comparison. The driver retains its existing early-convergence behavior; an
early stop cannot be relabelled a complete twenty-step comparison. Both wrappers
are quiet while collecting private logs. No long convergence run or performance
measurement is requested by these blocks.

## Local checks

```bash
PYTHONPATH=scripts python3 -m unittest test_simple_pressure_history \
  test_simple_pressure_profile test_simple_startup_probe test_simple_pressure_probe
```

Synthetic tests cover full histories, later field differences, missing or changed
settings, missing snapshots, reused reference records, altered controls,
initial-field mismatch, enabled recovery, runtime selection and failed launches.
The C++ options gate checks the supported snapshot/profile combination and
rejects partial coverage, gathered output and failure-capture combinations.
CUDA compilation and the real twenty-step runs remain user-executed gates.
