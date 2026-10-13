# Original-pressure first-step comparison

For the follow-on bounded run at these original targets, use the
[twenty-step original-control comparison](pressure_original_history.md).
It runs both solvers afresh with recovery disabled and checks every snapshot.

This experiment returns to the original OpenAccel pressure targets and saved
first step. It runs only MARS. It does not use the tightened-pressure experiment
as the reference, change momentum controls by default, enable correction-equation recovery,
or accept a failed pressure residual.

The first-step assembly audit found momentum and pressure equations matching
within the positive-row-scaling tolerance, but different pressure increments. The purpose here is to test the
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
and, unless explicitly selected below, momentum instance retain their existing controls. Native GPU SpMV and both
independent acceptance checks remain in effect.

`--pressure-solver-profile` requires explicit pressure tolerances and either the
one-step audit mode, full distributed snapshots, or opt-in
[retained pressure expansion](pressure_expansion.md).
The recipe below uses the audit. Rank zero reads the small settings file and broadcasts it;
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

## Momentum and pressure controls

Use `prepare --both-linear-solvers` to resolve both `coupled_navier_stokes`
and `pressure_correction` from the original deck. Each uses the reference's
lookup precedence: equation, `segregated_flow`, then `default`. The saved
reference deck and fields stay unchanged. Separate files configure the two
MARS instances; separate settings records verify every selected rank and step.
The public result reports `original_momentum_targets_preserved` and
`recorded_gpu_momentum_settings_verified` only after those checks pass.

The supported profile is unscaled Hypre GMRES or FlexGMRES with the restricted
BoomerAMG options above. A PETSc, Trilinos, MGR, matrix-scaling or other unsupported
configuration stops preparation with `unsupported_momentum_profile`; it is not
silently translated into Hypre. The previously measured reference defaults are
reused only after their executable/library identities are checked. Identical
CPU/GPU solvers are still not claimed: Hypre versions, row order, coarsening and
coarse relaxation can differ.

If preparation reports `unsupported_momentum_profile`, inspect the saved
configuration before rebuilding or allocating GPUs:

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
repo="$scratch/git/mars-v010-check"
umask 077
export PYTHONDONTWRITEBYTECODE=1
run=$(mktemp -d "$scratch/simple-linear-support-XXXXXX")
run_status=0
python3 "$repo/scripts/simple_pressure_profile.py" inspect-profiles \
  --baseline-pair "$scratch/simple-first-step-MPQGzD/pair" \
  --defaults-probe "$scratch/simple-hypre-defaults-PrkNXY/probe" \
  --output "$run/public.json" || run_status=$?
cat "$run/public.json"
printf 'Share only: %s\n' "$run/public.json"
exit "$run_status"
)
```

The report covers both equations. It exports known backend/method names, scaling
flags and the first rejected component; GPU restrictions list every failing
control in the checked profile. Unknown names, numeric values, lookup names and
exception messages are withheld. The same GPU validator is used by preparation
and inspection. Saved case, launch and defaults evidence is checked before and
after inspection; neither solver nor the runtime library loader is launched.
`inspection_complete` means the check finished. `both_profiles_supported` states
whether these settings can be mapped; neither field establishes solver convergence
or field agreement. This command works in the MARS terminal without mounting the
OpenAccel environment.

For `family: PETSc`, the inspector reads the flat `ksp_*`, `pc_*`, `sub_*`,
`mat_type` and `vec_type` options. It does not read Hypre's nested
`options.solver`/`options.precond` mappings as PETSc settings. Omitted PETSc method
and preconditioner types are reported as `unspecified`, not as GMRES/none.
Known declared option names and categorical values are exported; numeric values,
viewer paths and unrecognized names remain private. `runtime_settings_verified`
stays false: this is a source-based inspection, not a query of a live PETSc object.
The earlier generic inspector's GMRES/none labels for non-Hypre families were
not valid evidence and must not be used to choose a momentum algorithm.

In the inspected OpenAccel revision, `ContextPETSc.h::setOptions_` lowercases
the flat options before passing them to PETSc. `createKSP_` calls
`KSPSetFromOptions`, then overrides rtol, atol and the iteration cap with the
equation context's values; its generic minimum-iteration setting is not forwarded.
PETSc's stopping norm may be preconditioned, depending on KSP configuration;
see [KSPSetTolerances](https://petsc.org/release/manualpages/KSP/KSPSetTolerances/).
Matching tolerance numbers alone would therefore not establish the same stopping
test as MARS's unpreconditioned check. MARS SIMPLE currently has no PETSc momentum
backend. Inspection continues to reject that family instead of substituting Hypre.
The confirmed mixed PETSc-momentum/Hypre-pressure reference explains preparation's
family rejection; it does not establish the cause of the observed field difference.

For momentum `A du = b`, `du` is a velocity increment; matrix coefficients have
units kg/s and the RHS has units N. Both its Krylov target and the two independent
original-system checks use `max(atol, rtol ||b||)`. Pressure retains its own target
with RHS units kg/s. Neither target is tightened or relaxed by this mapping.
SIMPLE equations, relaxation, update order and the zero initial fields are unchanged.
Recovery remains disabled. MARS continues to reject failed linear solves even
though the inspected reference wrapper discards Hypre's return code.

This is an iteration-state comparison. Nonlinear MARS norms and stopping criteria
are unchanged; the reference's diagonally normalized component RMS is a different
quantity. Neither completion nor a profile match establishes nonlinear convergence.
The wall coefficient now uses the reference minimum normal distance:
`4 mu A_sample / max(h, epsilon)`, with `h = V_tet / A_sample`. The original
arithmetic is retained for `h >= epsilon`, and the first outer step still has zero
wall coefficient. This extreme-scale regularization difference has not been linked
to the observed pump discrepancy.

After building and passing the public gates in the block below, select
`--both-linear-solvers` on its preparation command. Start with the saved original
one-step reference, rather than repeating the interrupted strict-target long run:

```bash
python3 "$repo/scripts/simple_pressure_profile.py" prepare \
  --baseline-pair "$scratch/simple-first-step-MPQGzD/pair" \
  --defaults-probe "$scratch/simple-hypre-defaults-PrkNXY/probe" \
  --executable "$exe" --output-dir "$run/pair" --both-linear-solvers
```

`--history --both-linear-solvers` additionally requires fresh OpenAccel and MARS
captures as in [the original history recipe](pressure_original_history.md).
Do not reuse an old settings file or launch record to represent the new mapping.
Local validation includes the real CPU Hypre 2.33.0 and 3.1.0 regressions,
independent algebra and wall-distance checks, and collective acceptance tests
on 1/2/4 host ranks. These do not validate CUDA execution. The public GPU gate
now includes coupled three-component momentum profiles and the wall-distance
threshold check; execution of those new cases is pending.

## One block, in the MARS terminal

Run this once, in only one tmux pane. It builds the two affected executables,
runs public synthetic gates on 1/2/4 GPUs, and only then prepares and launches
the private MARS first step. It reuses the original `MPQGzD` reference capture
and the successful `PrkNXY` defaults probe. It does not launch OpenAccel.
The public gates exercise one-level/multilevel hierarchies, GMRES/FlexGMRES,
two numerical updates and the exchanged independent residual.

The pressure-profile cases use a seven-point diffusion operator on a synthetic
6 by 5 by 4 grid, with zero exterior Dirichlet values and a positive reaction.
For spacing one, the diagonal is `2*(wx+wy+wz)+0.125` and axial neighbors have
coefficients `-wx`, `-wy`, `-wz`. Positive axis weights change between updates.
The known solution, shuffled source IDs and uneven ownership are shared with
the CPU regression. Requested controls and actual hierarchy depth are checked
separately; the multilevel cases still require more than one level.

The initial GPU gate reported four combined settings/depth failures while its
linear solves passed. Its reused random matrix was unsuitable for multilevel
coverage: all 120 rows in both updates satisfy Hypre's all-weak row test at
`maxrowsum=0.9`. Local CPU Hypre 2.33.0 reproduces a single level with the requested
relaxation values intact. Hypre's GPU strength kernel applies the same condition.
The old combined output did not identify individual mismatched settings; the
corrected gate prints these and the measured depth. Solver controls, tolerances
and acceptance tests are unchanged.

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

## Recheck the saved original-pressure result

The reported `1xGpr2` comparison passed its profile and evidence checks but still
differed at the first pressure increment. This profile has not established field
parity. Both candidates passed their own original pressure residual targets.
That alone does not say whether they pass the same target on the same matrix:
independent positive row scaling preserves the exact solution but can change a
residual norm and an iterative stopping decision.

`pressure_common_system_checks` applies both saved pressure increments to the
reference matrix and RHS, using the single reference limit
`max(atol, rtol*||b_reference||₂)`. It also applies the reference increment to the
MARS matrix/RHS against the recorded MARS runtime limit. Unscaled row comparisons
use a common divisor for each pair of rows. The report distinguishes unit,
uniform nonunit and nonuniform row scaling within `1e-10`; these are approximate
classifications, not exact identities. The numerical details remain private.

- `distinct_increments_pass_common_fp64_check`: both differing candidates pass
  the same recomputed reference check. This is consistent with the original
  tolerance allowing different finite-accuracy iterates.
- `mars_candidate_fails_reference_target`: passing each solver's own check did
  not give acceptance on the common reference system. Inspect the scaling flags
  before attributing this to a backend.
- Assembly or referenced-copy discrepancies are reported ahead of those
  interpretations. A reference candidate failing its own check is reported
  separately when increments differ.

These checks multiply in FP64 and sum with `math.fsum`; they have no certified
evaluation-error bound. Cancellation near a stopping limit can therefore make
the verdict uncertain. They neither prove forward accuracy nor identify an AMG
or Krylov defect, and they do not change the existing field verdict.

Run this in the MARS terminal only. It reads the existing capture locally; no
rebuild, OpenAccel environment, allocation or flow run is required. Only the
new public checks are printed, with the complete public report saved alongside.

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
run=$(mktemp -d "$scratch/simple-pressure-common-XXXXXX")
status=0
python3 "$repo/scripts/simple_pressure_profile.py" compare \
  --pair "$scratch/simple-original-pressure-1xGpr2/pair" \
  --detail-dir "$run/private" --output "$run/public.json" || status=$?
python3 - "$run/public.json" <<'PY'
import json, sys
with open(sys.argv[1]) as stream:
    report = json.load(stream)
comparison = report.get('comparison', {})
print(json.dumps({
    'comparison_status': report.get('comparison_status'),
    'failed_check': report.get('failed_check'),
    'first_differing_stage': comparison.get('first_differing_stage'),
    'pressure_common_system_checks': comparison.get('pressure_common_system_checks')
}, indent=2))
PY
printf 'Share only: %s\n' "$run/public.json"
exit "$status"
)
```

## Local validation

Synthetic Python tests exercise original-target preservation, generic/cycle
setter ordering, explicit settings, saved defaults identity, reference reuse,
tampered inputs, missing/ignored actual controls and environment isolation.
The real CPU Hypre wrapper regression uses the same public matrix generator as
the GPU gate. It retains the old random case to check legitimate one-level
collapse and exercises diffusion with one and multiple levels, GMRES/FlexGMRES,
a fresh momentum instance and two cached pressure updates under ASan/UBSan.
These checks validate configuration and lifecycle, not CUDA execution or
performance. Run the public GPU gates above before the private comparison.
