# Original-pressure first-step comparison

This experiment returns to the original OpenAccel pressure targets and saved
first step. It runs only MARS. It does not use the tightened-pressure experiment
as the reference, change momentum controls, enable correction-equation recovery,
or accept a failed pressure residual.

The first-step assembly audit already found matching momentum and pressure
equations but different pressure increments. The purpose here is to test the
remaining pressure-solver configuration difference at the original tolerances.
For the assembled equation `A phi = b`, `phi` is a physical pressure increment in
Pa. SIMPLE still applies `p += alpha_p phi` and `u -= d grad(phi)`, with the same
Tet4 CVFEM quadrature, boundary contributions, pseudo-time term, density and
relaxation. A changed preconditioner can change a finite-accuracy iterate even
when the equation and stopping target agree. No discretization term is changed.

`simple_pressure_profile.py` resolves explicit OpenAccel options against the
previously measured CPU Hypre 3.1.0 defaults. It checks that the saved defaults
probe used the same reference executable and library hashes. The copied probe
log must agree with its saved public report and loaded-library record. This
reuses recorded evidence; it does not inspect the old process or certify its
matrix-dependent hierarchy.

Two explicit substitutions keep matrix/vector work on the GPU:

| Reference choice | GPU profile |
|---|---|
| HMIS coarsening (10) | PMIS (8) |
| Coarse Gaussian elimination (9) | l1-Jacobi (18) |

The source basis is OpenAccel `src/solver/HYPRE/ContextHYPRE.h` at public revision
`0d69041ba1afda63e9e4328d9e0d9834bba37756`, and Hypre 2.33.0/3.1.0:
`parcsr_ls/par_coarsen.c` stages HMIS data on host;
`parcsr_ls/par_gauss_elim.c` stages coarse matrix/vector work on host, including
the device-LU variants. Those variants are therefore not used as a purported
GPU-only replacement. GPU down/up relaxation 13/14 is retained when the reference
omits a generic relaxation option. The one-level fallback is explicitly 6 and
one sweep, following `parcsr_ls/par_cycle.c`. Explicit generic relaxation/sweep
options are applied before the cycle overrides. Unsupported combinations fail.

The reference restart, method, iteration cap, default minimum iterations and
original rtol/atol are applied to the pressure instance. OpenAccel does not
forward its generic `min_iterations` to Hypre; the captured Hypre default is
used. AMG remains one V-cycle with tolerance zero. The production default path
and momentum instance retain their existing controls. Native GPU SpMV and both
independent acceptance checks remain in effect.

`--pressure-solver-profile` requires explicit pressure tolerances and the existing
one-step audit mode. Rank zero reads the small settings file and broadcasts it;
there are no added host matrix/vector computations. The audit already enables
private field/matrix output and its device-to-host copies. Each rank also writes
its actual Krylov and AMG settings after setup. The comparison requires these
to match the requested profile, including the down/up/coarse relaxation values.

The new pair retains an immutable reference to the original pair. It does not
copy an old launch record and relabel it as fresh. The MARS executable may be
rebuilt, but its runtime library hashes, rank count, launcher and recorded
numerical environment must match the baseline. Known scalar environment controls
are restored only in the child process. Unknown or device-binding changes stop
before launch. Both baseline and new output hashes are checked again during
comparison. OpenAccel runtime mounts need not be present in the MARS terminal.

## One block, in the MARS terminal

Run this once, in only one tmux pane. It builds the two affected executables,
runs public synthetic gates on 1/2/4 GPUs, and only then prepares and launches
the private MARS first step. It reuses the original `MPQGzD` reference capture
and the successful `PrkNXY` defaults probe. It does not launch OpenAccel.
The public gates exercise one-level/multilevel hierarchies, GMRES/FlexGMRES,
two numerical updates and the exchanged independent residual.

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
repo="$scratch/git/mars-v010-check"
umask 077
mkdir -p "$scratch/tmp"
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
python3 -c 'import numpy, netCDF4, yaml'
test "$(git -C "$repo" branch --show-current)" = cstone
git -C "$repo" fetch origin cstone
git -C "$repo" merge --ff-only refs/remotes/origin/cstone
cmake -S "$repo" -B "$repo/build-hypre" \
  -DMARS_ENABLE_CUDA=ON -DMARS_ENABLE_MPI=ON -DMARS_ENABLE_TESTS=ON \
  -DMARS_ENABLE_UNSTRUCTURED=ON -DMARS_ENABLE_HYPRE=ON \
  -DMARS_ENABLE_FEM_EXAMPLES=ON -DMARS_ENABLE_SEGREGATED=ON \
  -DCMAKE_PROJECT_mars_INCLUDE="$repo/tests/reference/openaccel/distributed_matrix/inject.cmake"
cmake --build "$repo/build-hypre" --parallel 4 --target \
  mars_segregated_simple mars_distributed_matrix_cuda_gate
run=$(mktemp -d "$scratch/simple-original-pressure-XXXXXX")
mkdir "$run/public-gates"
printf 'Results: %s\n' "$run"
git -C "$repo" rev-parse HEAD > "$run/public-gates/mars-revision.txt"
exe="$repo/build-hypre/examples/distributed/unstructured/mars_segregated_simple"
gate="$repo/build-hypre/mars_distributed_matrix_cuda_gate"
sha256sum "$exe" "$gate" > "$run/public-gates/executables.sha256"
for np in 1 2 4; do
  srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node="$np" \
    --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
    "$HOME/affinity/bind_numa.sh" "$gate" \
    2>&1 | tee "$run/public-gates/matrix-$np.log"
done
python3 "$repo/scripts/simple_pressure_profile.py" prepare \
  --baseline-pair "$scratch/simple-first-step-MPQGzD/pair" \
  --defaults-probe "$scratch/simple-hypre-defaults-PrkNXY/probe" \
  --executable "$exe" --output-dir "$run/pair"
status=0
python3 "$repo/scripts/simple_pressure_profile.py" run \
  --pair "$run/pair" --output "$run/launch-public.json" || status=$?
cat "$run/launch-public.json"
if (( status != 0 )); then exit "$status"; fi
python3 "$repo/scripts/simple_pressure_profile.py" compare \
  --pair "$run/pair" --detail-dir "$run/private-comparison" \
  --output "$run/comparison-public.json"
cat "$run/comparison-public.json"
printf 'Share only: %s\n' "$run/comparison-public.json"
)
```

All new files, logs and temporary data are on capstor scratch. The existing
affinity script is only executed from home; it is not used as an output location.
The pair, settings, raw private logs, matrices and fields stay on the cluster.
Share only `launch-public.json`, `comparison-public.json` and public gate results.
The script refuses to overwrite existing pairs or public reports.

`capture_complete` means the one-step run finished, not that it matched.
`comparison_status=completed` means the evidence was checked; inspect the nested
`first_step_stage_matches` and `all_snapshots_match` for the result. A mismatch
remains a mismatch. Equal first-step results would justify the subsequent
original-control twenty-step comparison, not a pump-convergence claim.
The CPU/GPU library versions, ordering and two AMG substitutions remain different;
`identical_linear_solvers_verified` stays false.

## Local validation

Synthetic Python tests exercise original-target preservation, generic/cycle
setter ordering, explicit settings, saved defaults identity, reference reuse,
tampered inputs, missing/ignored actual controls and environment isolation.
The real CPU Hypre wrapper regression exercises GMRES/FlexGMRES with one and
multiple levels, a fresh momentum instance and two cached pressure updates under
ASan/UBSan. These checks validate configuration and lifecycle, not CUDA execution
or performance. Run the public GPU gates above before the private comparison.
