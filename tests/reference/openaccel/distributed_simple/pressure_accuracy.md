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

```bash
PYTHONPATH=scripts python3 -m unittest test_simple_pressure_probe test_simple_first_step_audit test_simple_startup_probe test_simple_snapshot_compare test_prepare_simple_deck
```
