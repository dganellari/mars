# Earlier pressure recovery on the saved system

The longer retained-pressure flow was interrupted on 2026-10-12. The user-side
sample showed outer iterations advancing, recovery in most recent iterations,
mixed momentum history and flat continuity history. This establishes expensive
extra work, but does not measure its share of wall time or explain nonlinear
convergence. Do not repeat that long run unchanged.

This experiment changes only when recovery starts on the frozen pressure system.
It makes no production SIMPLE change. It is not matched OpenAccel stopping
behavior, a nonlinear convergence fix, or a performance result.

## Algebra and limits

Both arms solve the same saved nodal pressure-increment equation `A x = b` from
zero, using the same executable, GPU AMG profile, rank count and launch controls.
Assembly, boundary conditions, physical controls and the pressure target stay
unchanged. The ordinary arm has its saved initial cap `M`. The early arm uses
`K = kdim`, one GMRES restart cycle, with `miniter <= K < M`. Both then allow the
existing four correction solves, each with its original cap `M`. Thus the early
arm's Krylov work is bounded by `K + 4 M`, not `K` alone.

`--initial-restart-cycles 1` is restricted to the GPU reference profile with
retained recovery enabled. The replay reads back the actual initial cap and
restores `M` before corrections. It records both; stopping diagnostics compare
the first solve against its actual cap. No convergence flag overrides the final
independent residual check on the captured rows:

`||b - A x_high - A x_low|| <= max(atol, rtol ||b||)`.

The existing compensated evaluator also accounts for its evaluation bound.
If no retained update is produced, the ordinary final candidate is checked.
Neither arm is accepted when that check fails or is inconclusive, or when a
fatal backend error occurred. Earlier recovery may need more work or may fail.

## One GPU comparison, using existing data

Run the following **once, in the MARS terminal**. It builds only the replay and
independent checker, then performs two pressure-system replays. It does not
launch either full flow solver. The wrapper uses the saved complete `srun`
command, including binding and rank count; do not wrap it in another `srun`.

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
  --target mars_simple_pressure_replay mars_simple_pressure_residual_check
umask 077
mkdir -p "$scratch/tmp"
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
run=$(mktemp -d "$scratch/simple-pressure-budget-XXXXXX")
printf 'Private results: %s\n' "$run"
for arm in baseline early; do
  extra=()
  if test "$arm" = early; then extra=(--initial-restart-cycles 1); fi
  status=0
  python3 "$repo/scripts/simple_pressure_replay.py" replay \
    --capture-run "$scratch/simple-pressure-frozen-RlGQaI/capture" \
    --backend mars --profile gpu-reference \
    --gpu-profile-pair "$scratch/simple-original-pressure-1xGpr2/pair" \
    --recovery-rounds 4 --retain-update-expansion "${extra[@]}" \
    --executable "$repo/build-hypre/examples/distributed/unstructured/mars_simple_pressure_replay" \
    --output-dir "$run/$arm" || status=$?
  if test "$status" -ne 0; then cat "$run/$arm/public.json"; exit "$status"; fi
done
status=0
python3 "$repo/scripts/simple_pressure_replay.py" compare-budget \
  --capture-run "$scratch/simple-pressure-frozen-RlGQaI/capture" \
  --baseline-run "$run/baseline" --early-run "$run/early" \
  --checker "$repo/build-hypre/mars_simple_pressure_residual_check" \
  --output "$run/comparison-public.json" || status=$?
cat "$run/comparison-public.json"
printf 'Share only: %s\n' "$run/comparison-public.json"
exit "$status"
)
```

For direct-launch inspection, each arm's private `launch-start.json` contains
the entire command as an argument list. Its tail is the replay executable,
`capture/system`, that arm's `gpu-reference.settings` and its new `result`
directory. Do not rerun it in the same directory: replay output is write-once.

## Reading the result

`compare-budget` verifies both manifests, loaded libraries, source profile,
same executable and saved launcher, captured system, recorded environment,
correction budget and final target. Only the initial cap may differ.
It counts initial plus all correction Krylov iterations, including rejected
corrections. Exact counts and times remain private.

The measured time is the maximum over ranks of initial solve plus recovery
wall time. Device work is synchronized at those boundaries. It includes
recovery residuals and audit work, and excludes AMG setup, disk output, queue
wait and the independent checker. Each replay rebuilds its AMG hierarchy.
The two launches may use different nodes; one pair is an exploratory screen,
not a repeatability or production speedup measurement.

`less_work_and_shorter_measured_solve` requires both candidates to pass,
fewer total Krylov iterations and at least a 10% shorter measured solve.
This threshold screens the observation; it is not a statistical confidence
bound. `candidate_rejected` forbids claiming a benefit even when the rejected
arm ran faster. `no_clear_work_and_time_benefit` leaves production unchanged.
All outcomes leave nonlinear convergence and OpenAccel field parity unverified.
