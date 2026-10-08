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
This block rebuilds the three affected targets, then reuses that pair's launcher,
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

## MARS terminal: compare all three candidates

The host checker reads only the frozen diagnostic files. Production assembly,
solving and halo exchange remain on the GPU. Share only the final public JSON.

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
repo="$scratch/git/mars-v010-check"
umask 077
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
run=$(cat "$scratch/simple-pressure-frozen-current.txt")
status=0
python3 "$repo/scripts/simple_pressure_replay.py" compare \
  --capture-run "$run/capture" --mars-run "$run/mars" --reference-run "$run/reference" \
  --checker "$repo/build-hypre/examples/distributed/unstructured/mars_simple_pressure_residual_check" \
  --output "$run/comparison-public.json" || status=$?
cat "$run/comparison-public.json"
printf 'Share only: %s\n' "$run/comparison-public.json"
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
