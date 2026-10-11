# Controlled first-step pressure accuracy

The saved MPQGzD first-step audit matches momentum assembly, the predictor and
influence, and pressure assembly within its tolerances. The first differing
field is the raw pressure increment. The user's subsequent saved-file check
(`simple-pressure-audit-3sgMma/public.json`) reports that both pressure residuals
satisfy the same declared target, and referenced MARS ghost copies equal their
owners. Both targets are looser than the common relative `1e-8` diagnostic.
This does not establish a failed solve or explain the field difference yet.

For the same discrete pressure equation `A_p phi = b_p`, this experiment requests
`||b_p - A_p phi||_2 <= rtol * ||b_p||_2`, with
`rtol = min(original_rtol, 1e-10)` and `atol = 0` in both codes. It tests whether
reducing the permitted residual brings the first-step fields into agreement.
No equation, boundary condition, relaxation or nonlinear criterion changes.
Small residuals alone do not bound solution error without conditioning evidence.

`scripts/simple_pressure_probe.py` creates a fresh one-step pair from a completed
first-step capture with `pressure_linear_policy=reference`. It deep-copies the
resolved OpenAccel pressure configuration into a dedicated `pressure_correction`
block, preserving its backend, preconditioner, restart dimension, iteration cap
and other options. Only the two pressure targets change. Shared fallback or
named configurations used by momentum are not edited; the resolved momentum
configuration is checked explicitly. The reference policy emits the same
explicit pressure targets for MARS's Krylov and original-CSR residual checks.

This policy does not copy the reference preconditioner, restart dimension,
iteration limits or AMG options into MARS. Use the read-only
[saved pressure settings check](pressure_solver_settings.md) to identify those
differences before attributing a pressure mismatch to the equations.

The helper reuses each solver's captured executable, rank count and launcher
arguments. It rejects changed executable/library hashes or recorded solver
environment overrides before launching. Use the same OpenAccel and MARS runtime
environments as the baseline. **Do not rebuild either executable for this test.**
Pulling these Python scripts does not change the captured executable.

Runtime identity covers the recorded binaries, libraries, launchers and
`MARS_`, `HYPRE_`, `CUDA_`, `MPICH_`, `OMP_` environment values. It does not
prove identical hardware or every possible external setting. Older captures did
not record `PETSC_OPTIONS` or `PETSC_OPTIONS_YAML`; this helper rejects those
overrides if present rather than guessing their historical values. A changed
runtime needs an explicitly reviewed new baseline, not bypassing the checks.

An environment mismatch now reports only allowlisted variable names and whether
other names changed; no values or unknown names are exported. The opt-in
`run --restore-solver-environment` option restores known scalar MARS solver,
halo and execution controls from the verified baseline, including removing a
current override that was absent there. It applies only to the MARS child
process; the interactive shell stays unchanged. It does not restore GPU binding,
library paths, geometry-debug switches, file-output paths or unknown variables.
The full environment equality check still runs after restoration, and the launch
record stores the actual child environment. Binary and library checks stay strict.

All decks, matrices, fields, identifiers, logs and detailed errors stay private
on capstor. Only the fixed-label public JSON reports may be shared. The scripts
run where the user owns those files; no private artifact is required locally.

## OpenAccel terminal

Use the existing working OpenAccel environment. No build is needed. The stored
MPI launcher is reused exactly, including its account, rank count and time limit.
The run command prints no private log or solver command. Save the printed pair
directory for the MARS terminal.

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
repo="$scratch/git/mars-v010-check"
git -C "$repo" pull --ff-only
python3 -c 'import numpy, netCDF4, yaml'
umask 077
mkdir -p "$scratch/tmp"
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
run=$(mktemp -d "$scratch/simple-pressure-accuracy-XXXXXX")
pair="$run/pair"
probe="$repo/scripts/simple_pressure_probe.py"
python3 "$probe" prepare \
  --baseline-pair "$scratch/simple-first-step-MPQGzD/pair" \
  --output-dir "$pair"
printf 'Private pair: %s\n' "$pair"
status=0
python3 "$probe" run --pair "$pair" --solver openaccel \
  --output "$run/reference-public.json" || status=$?
cat "$run/reference-public.json"
printf 'Private pair: %s\nShare only: %s\n' "$pair" "$run/reference-public.json"
exit "$status"
)
```

## MARS terminal

Use the existing MARS uenv and the unchanged executable used for MPQGzD. Enter
the exact pair path printed above. The recorded one-step launcher permits exit
2 (iteration limit); a pressure-solve rejection is a failed capture, not a usable
field comparison. Both original launch histories remain intact.

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
repo="$scratch/git/mars-v010-check"
python3 -c 'import numpy, netCDF4, yaml'
umask 077
mkdir -p "$scratch/tmp"
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
printf 'Fresh pressure pair directory: '
IFS= read -r pair
test -f "$pair/pair.json"
probe="$repo/scripts/simple_pressure_probe.py"
status=0
python3 "$probe" run --pair "$pair" --solver mars \
  --output "$pair/mars-public.json" || status=$?
cat "$pair/mars-public.json"
if test "$status" -ne 0; then exit "$status"; fi
python3 "$probe" compare --pair "$pair" \
  --detail-dir "$pair/accuracy-private" \
  --output "$pair/accuracy-public.json" || status=$?
cat "$pair/accuracy-public.json"
printf 'Share only: %s\n' "$pair/accuracy-public.json"
exit "$status"
)
```

## Interpreting the result

- `first_step_matches`: both independently recomputed pressure residuals meet
  the tighter targets and all audited stages match. If the baseline pressure
  differed, this supports permitted linear-solve error as the explanation for
  that first step. It does not establish full pump convergence or a general fix.
- `first_step_still_differs`: the tighter residuals pass and assembly/predictor
  comparisons still match, but at least one corrected field or pressure
  intermediate differs. Conditioning and the remaining differences need analysis;
  this alone does not prove an assembly defect.
- `upstream_stage_mismatch`: momentum or pressure assembly/predictor agreement
  failed in the tightened pair; do not attribute that result to pressure accuracy.
- `pressure_target_not_met`, a failed capture, or invalid evidence: inconclusive.
  Preserve the attempt. Do not raise iteration caps, loosen acceptance, or accept
  an incomplete binary audit automatically. MARS can reject before writing the
  pressure observer; the public launch JSON then contains safe failure diagnostics.

The comparison recomputes both old and new audits into fresh private directories.
It neither overwrites prior reports nor claims an automatic root-cause verdict.
This is a one-step experiment; no new GPU result has been established by preparing
the helper. Synthetic tests check shared solver lookup, unchanged momentum/caps,
strict target mapping, missing/tampered evidence, runtime drift, residual failure,
upstream mismatch, failed launches and private-output suppression.

## Retry W3SoRZ after the environment preflight stopped

The OpenAccel capture completed; the MARS `solver_environment_changed` preflight
stopped before launching a solver. The specific difference was not identified by
that older report. Keep the completed reference and the old failure report.
Run the following in the MARS terminal with its working runtime and Python.
No C++ rebuild or OpenAccel rerun is needed. This explicitly restores only the
allowlisted controls for the child process; any remaining mismatch is reported
safely and still prevents launching. Reports go into a new scratch directory.

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
repo="$scratch/git/mars-v010-check"
git -C "$repo" pull --ff-only
python3 -c 'import numpy, netCDF4, yaml'
umask 077
mkdir -p "$scratch/tmp"
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
pair="$scratch/simple-pressure-accuracy-W3SoRZ/pair"
probe="$repo/scripts/simple_pressure_probe.py"
retry=$(mktemp -d "$scratch/simple-pressure-retry-XXXXXX")
status=0
python3 "$probe" run --pair "$pair" --solver mars \
  --restore-solver-environment --output "$retry/launch-public.json" || status=$?
cat "$retry/launch-public.json"
if test "$status" -ne 0; then
  printf 'Share only: %s\n' "$retry/launch-public.json"
  exit "$status"
fi
python3 "$probe" compare --pair "$pair" --detail-dir "$retry/private" \
  --output "$retry/comparison-public.json" || status=$?
cat "$retry/comparison-public.json"
printf 'Share only: %s\n' "$retry/comparison-public.json"
exit "$status"
)
```

```bash
PYTHONPATH=scripts python3 -m unittest test_simple_pressure_probe test_simple_first_step_audit test_simple_startup_probe test_simple_snapshot_compare test_prepare_simple_deck
```

## Saved gradient replay after hkkN7u

The public hkkN7u report verifies the pressure-only experiment and recorded runtime
identity. Both pressure residual checks now pass; the raw pressure increment and
corrected velocity and pressure match the comparison tolerances. Only the
pressure-increment gradient still differs. This is **one step**, not a converged
pump comparison. A scalar passing its tolerance is not an identical scalar:
differentiation can amplify the remaining difference.

For the supported shifted, incremental Tet4 reconstruction, both source paths use

```text
V_i = sum(t containing i) V_t/4
A_ij,t = (V_t/4) (grad N_j - grad N_i)
(G phi)_i = sum(t containing i, j != i) 0.5*(phi_j-phi_i)*A_ij,t / V_i
```

Here `phi` is the raw pressure correction in Pa, `G phi` has units Pa/m, and
`N_i` is the affine tetrahedral basis. The shifted boundary sample equals its
node value, so the incremental boundary contribution is zero for the supported
wall/inlet/outlet path. This operator is not affine-exact at every boundary
node. Relevant implementations are MARS `mars_segregated_simple.hpp`
(`SimpleGradientInterior/Boundary/Finish`) and pinned OpenAccel `nodeField.hpp`
(`updateGradientField`), with shifted Tet4/Tri3 master elements. MARS sums before
dividing by dual volume; OpenAccel divides each term before accumulation.
Their geometric evaluation paths also differ, so exact floating-point identity
is not assumed.

`compare --gradient-audit` reads the already hashed input mesh and captures on
the user's machine, reconstructs this common operator, and checks

```text
g_M - G(phi_M)
g_R - G(phi_R)
(g_M - g_R) - G(phi_M - phi_R)
```

It applies `G` directly to the pressure difference. The reference intermediates
must be stored as float64, and duplicate reference copies must agree exactly.
Connectivity is streamed by block; only Tet4 is accepted. The production solver,
captured binaries, original reports, and launch records are unchanged. This is
optional CPU postprocessing of existing files, not a CPU solver path.

Each reconstruction uses consistency tolerance `1e-10*(P0/L + B)`, where `B`
is the absolute edge-area/pressure accumulation divided by dual volume. The
closure tolerance sums the two individual allowances and the delta-field
allowance. These are diagnostic tolerances, **not certified roundoff bounds**.
An explanation additionally requires every allowance to be below 1% of the
observed maximum gradient difference. Numerical values remain in the fresh
private report; the public report contains only fixed labels and booleans.

- `input_difference_explains_gradient_within_replay_tolerance`: both individual
  gradients reconstruct, the difference closes, and the check resolves that
  difference. This supports amplification of the remaining scalar difference;
  it does not prove backend convergence or nonlinear pump agreement.
- `reconstruction_mismatch`: at least one reconstruction/closure fails. Further
  geometry, mapping, capture or implementation checks are needed; it does not
  uniquely identify a kernel bug.
- `insufficient_replay_resolution`: the allowed replay error is too large for
  that attribution.
- `gradients_within_field_tolerance`: the saved gradients already meet the
  original `1e-5` scaled check and reconstruct successfully.

The original stage verdict is preserved even if the difference is explained.
No tolerance is relaxed to turn `first_step_still_differs` into a pass.

### Exact checks after the unresolved ffIlzR replay

The user's `simple-gradient-replay-ffIlzR/public.json` reports
`insufficient_replay_resolution`: the allowance above cannot discriminate the
observed difference. That does not measure the replay's actual arithmetic error
or establish a gradient implementation defect.

For gradients that fail the original field tolerance, `--gradient-audit` now
also computes an exact rational reconstruction at selected nodes. Every stored
binary64 coordinate, pressure and captured gradient is converted to a rational
number **before subtraction**. With `a=x1-x0`, `b=x2-x0`, `c=x3-x0`, let
`C1=b cross c`, `C2=c cross a`, `C3=a cross b`, `C0=-(C1+C2+C3)`, and
`d=a dot C1`. Each incident tetrahedron contributes

```text
numerator_i += sum(j != i) sign(d)*(phi_j-phi_i)*(Cj-Ci)/48
volume_i += abs(d)/24
```

Dividing the accumulated numerator by the accumulated volume gives the same
operator above, without replay rounding. This is exact for the saved inputs;
it does not reproduce either solver's floating-point geometry or accumulation.

Selection takes up to eight nodes for each of four scores: observed gradient
difference, MARS replay error, reference replay error, and difference-closure
error. Only nodes failing the original field tolerance are candidates. Scores
come from the floating-point replay, with source-row order breaking ties. The
union has at most 32 nodes; it is not a certified set of global worst exact
errors. Every tetrahedron touching a selected node is included. More than
20,000 such elements yields `inconclusive_work_limit`, never a partial-star pass.

The public `exact_selected_node_checks` compares each individual reconstruction
error and the closure error against **1% of that node's own maximum-component
captured gradient difference**, using exact comparisons. It reports:

- `input_difference_explains_checked_nodes`: all three errors meet that check.
- `captured_gradient_difference_not_explained_at_checked_nodes`: closure fails.
- `individual_reconstruction_discrepancy_at_checked_nodes`: closure passes but
  at least one individual reconstruction fails; shared errors must not hide.

These are selected-node conclusions. `all_field_tolerance_failures_checked`
states whether selection included every node failing the original tolerance;
even then, nodes within that tolerance are not audited exactly. Global operator
equivalence, runtime roundoff bounds, and the cause of the pressure difference
remain unproved. The original global assessment, stage verdicts and pressure
residual checks are unchanged. Exact ratios and node identifiers are written
only to the fresh private report; the public report contains labels and booleans.

Run **once in one terminal**, with the working NumPy/netCDF4/PyYAML environment.
Both terminals share the same repository; do not fetch or merge concurrently.
**No build, allocation, OpenAccel rerun or MARS rerun is required.**

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
repo="$scratch/git/mars-v010-check"
if test "$(git -C "$repo" branch --show-current)" != cstone; then
  echo 'Stop: this checkout is not on cstone.'
  exit 1
fi
git -C "$repo" fetch origin cstone
git -C "$repo" merge --ff-only refs/remotes/origin/cstone
python3 -c 'import numpy, netCDF4, yaml'
umask 077
mkdir -p "$scratch/tmp"
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
run=$(mktemp -d "$scratch/simple-gradient-exact-XXXXXX")
status=0
python3 "$repo/scripts/simple_pressure_probe.py" compare \
  --pair "$scratch/simple-pressure-accuracy-W3SoRZ/pair" \
  --gradient-audit --detail-dir "$run/private" \
  --output "$run/public.json" || status=$?
cat "$run/public.json"
printf 'Share only: %s\n' "$run/public.json"
exit "$status"
)
```

Local validation adds public synthetic tests for the analytic single-tet action,
constant fields, orientation/translation/scaling, block connectivity, amplified
scalar errors, false explanations, storage/copy rejection, and preserved capture
provenance. Exact checks cover cancellation before subtraction, complete stars
across blocks, work-limit rejection, deterministic selection, and preservation
of a failed pressure target even when the gradient difference is explained.
A compiled host gate checks the replay against the production C++
geometry kernels on shared, skew tetrahedra; this does not validate GPU execution.

## Twenty steps with the verified pressure targets

The user's `simple-gradient-exact-wnbVQ3/public.json` reports that exact replay
explains the gradient difference at every selected node to within 1% of its
local difference. Not all failing nodes were selected. The corrected first-step
velocity and pressure already match; this supports moving to a longer history
without changing a gradient kernel. It does not establish global gradient
agreement, a complete root cause, or nonlinear convergence.

`simple_pressure_probe.py extend` prepares a fresh 20-step run **from the same
zero initialization**, not a restart from the saved solution. Its baseline is the
completed tighter one-step pair W3SoRZ. Preparation rechecks that pair's pressure
residuals, assembly, predictor, influence, raw increment and corrected fields.
The remaining gradient mismatch does not block this history experiment.

The helper preserves the tight pressure targets, all other numerical controls,
and the recorded binaries, libraries, rank counts, launchers and solver
environment. Only iteration limits and diagnostic output change. The reference
momentum and pressure configurations are copied into dedicated blocks with
`write_system=false`; the original shared/named configurations remain intact.
This disables the inherited one-step matrix dumps without altering the solves.
The new MARS launch also omits `--first-step-audit`.

The existing `run` preflight rejects runtime drift against W3SoRZ, which is still
checked against its original baseline. The MARS-only environment restoration
option has the same restricted scope as above. **Do not rebuild either binary.**
Preparation and comparison read the saved private matrices locally; only the
fixed-label public reports may be shared.

Run the following first in the **OpenAccel terminal**, with its working runtime
and Python. Fetch only in this terminal; wait for this block before using the
MARS block. Each `run` reuses that solver's recorded `srun` command, including
the original allocation options and MARS GPU binding.

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
run=$(mktemp -d "$scratch/simple-pressure-history-XXXXXX")
pair="$run/pair"
probe="$repo/scripts/simple_pressure_probe.py"
python3 "$probe" extend \
  --baseline-pair "$scratch/simple-pressure-accuracy-W3SoRZ/pair" \
  --output-dir "$pair"
printf '%s\n' "$pair" > "$scratch/simple-pressure-history-current.txt"
printf 'Private pair: %s\n' "$pair"
status=0
python3 "$probe" run --pair "$pair" --solver openaccel \
  --output "$run/reference-public.json" || status=$?
cat "$run/reference-public.json"
exit "$status"
)
```

After `capture_complete`, use the **MARS terminal**, with its working runtime and
Python. This reads the new pair's scratch pointer; it neither fetches nor builds.

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
probe="$scratch/git/mars-v010-check/scripts/simple_pressure_probe.py"
python3 -c 'import numpy, netCDF4, yaml'
umask 077
mkdir -p "$scratch/tmp"
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
pair=$(cat "$scratch/simple-pressure-history-current.txt")
test -f "$pair/reference/launch.json"
status=0
python3 "$probe" run --pair "$pair" --solver mars \
  --restore-solver-environment --output "$pair/mars-public.json" || status=$?
cat "$pair/mars-public.json"
if test "$status" -ne 0; then exit "$status"; fi
python3 "$probe" compare --pair "$pair" \
  --detail-dir "$pair/history-private" \
  --output "$pair/history-public.json" || status=$?
cat "$pair/history-public.json"
printf 'Share only: %s\n' "$pair/history-public.json"
exit "$status"
)
```

`twenty_step_fields_match` means states 0–20 meet the unchanged velocity/U and
pressure/(rho U²) tolerance `1e-5`, without a pressure offset. Otherwise
`twenty_step_fields_differ` identifies the first differing state and field.
A completed comparison returns exit 0 in either case; read its `outcome`.
Failed launches, missing states, or runtime changes cannot pass. Early MARS
convergence before step 20 is incomplete evidence for this experiment; do not
loosen convergence criteria to force the history.

No matrices or intermediate gradients are captured in this history, so the
comparator does not independently recompute pressure residuals at every step.
`pressure_residuals_recomputed_at_every_step` and
`nonlinear_convergence_verified` remain false even if all saved fields match.
The original one-step captures and reports remain unchanged.

### Inspecting a rejected pressure solve

If the history stops at a pressure solve, re-export its saved log before changing
solver controls. Run this in **one terminal only**; it needs neither a rebuild
nor another GPU run. The existing capture and comparison inputs stay unchanged.

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
pair=$(cat "$scratch/simple-pressure-history-current.txt")
test -f "$pair/mars/run.log"
test -f "$pair/mars/run.exit"
summary=$(mktemp -d "$scratch/simple-pressure-stop-XXXXXX")/public.json
python3 "$repo/scripts/simple_public_diagnostics.py" \
  --log "$pair/mars/run.log" --exit-file "$pair/mars/run.exit" \
  --output "$summary"
cat "$summary"
printf 'Share only: %s\n' "$summary"
)
```

The extra flags distinguish a first-step rejection from a later one, zero from
positive Krylov iterations, and zero from positive RHS norms and acceptance
limits. Missing, invalid, or conflicting records remain unknown (`null`). Raw
counts and norms stay private. These flags do not identify a stopping cause.

`krylov_work_residual_within_limit` checks the logged work vector, which need not
contain the final `b-Ax` after an early stop. It cannot override either explicit
residual check or turn a failed run into a pass. Likewise, being within the
pressure audit's roundoff upper bound does not certify an accurate solution or
prove that a tighter solve is impossible.

## Experimental correction is excluded from this comparison

`--pressure-refinement` defaults to `0` and must remain disabled here. OpenAccel's
reference Hypre wrapper does not have MARS's new bounded defect-correction loop.
Preserving the equation and target alone is insufficient to keep this controlled
comparison of iteration histories unchanged. The comparator rejects enabled,
ambiguous or attempted recovery with `pressure_refinement_not_comparable` (or
the enclosing capture/control failure label).

The optional recovery mechanism, its tests and limitations are documented in
[pressure refinement](../simple_performance/pressure_refinement.md). They are a
separate robustness experiment and do not resolve the reference comparison.
The default still rejects pressure solutions that miss the original checks.
After rebuilding, prepare a fresh capture with recovery disabled; do not replace
a recorded executable or rewrite an old pair's provenance.
