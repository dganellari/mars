# Pressure solver settings in saved comparisons

The `reference` pressure policy copies **only `rtol` and `atol`** into MARS.
`mapped_controls_verified` does not certify equal linear solvers. The existing
comparison explicitly reports `identical_linear_solvers_verified=false`.

`scripts/simple_pressure_settings.py` compares the resolved pressure block in a
saved OpenAccel deck with MARS's source settings and **recorded** launch environment.
It reads the deck, case, pair, launch records, logs and exit files on the user's
machine. It never opens the mesh, matrices or field files, probes the current
runtime, or launches either solver. An unsuccessful history can still be inspected
if both complete launch records exist; inspection success is not solver success.

The output lists only fixed setting labels under `known_matches`,
`known_differences` and `unresolved_settings`. No private option values, unknown
keys, paths, boundary names or exception messages are exported. Enabled or
attempted experimental pressure refinement is rejected.

## Source contract

Inspected on 2026-10-08: OpenAccel
`0d69041ba1afda63e9e4328d9e0d9834bba37756`, MARS `2c0c669c`.
The model covers calls made by `ContextHYPRE`, `HypreSimpleSolve<1>` and
`HypreGMRESSolver`. It does not infer the linked Hypre library's omitted defaults
or prove that a captured binary was compiled from these revisions.

| Control | OpenAccel source | MARS SIMPLE source |
| --- | --- | --- |
| Krylov method | `options.solver.type`; GMRES if `options` is absent | GMRES; `MARS_HYPRE_FLEXGMRES` selects FlexGMRES |
| Preconditioner | None unless `options.precond` specifies one | BoomerAMG |
| Maximum iterations | Top-level `max_iterations`, default 20 | 2000 |
| Restart dimension | `options.solver.kdim`; otherwise Hypre default | 100 |
| Minimum iterations | Top-level `min_iterations` is **not forwarded** by `ContextHYPRE` | `MARS_HYPRE_MINITER`, default 3; nonpositive means no setter |
| Pressure tolerances | Top-level `rtol`, `atol` | Explicit `--pressure-linear-rtol`, `--pressure-linear-atol` override the wrapper defaults |
| BoomerAMG controls | Recognized `options.precond` entries; otherwise library defaults | Explicit wrapper choices, with eight `MARS_AMG_*` overrides |
| AMG application | Preconditioner `maxiter=1`, `tol=0` enforced after parsing | `maxiter=1`, `tol=0` |
| Initial linear solution | Zeroed for every solve | Zeroed for every solve |
| Acceptance | Reference backend behavior | Additional Hypre and MARS explicit residual checks |

Relevant code:

- OpenAccel `src/solver/HYPRE/ContextHYPRE.h`: `setupSolver_`, `setSolver_`,
  `setPreconditioner_`, `setGMRES_`, `setFlexGMRES_`, `setBoomerAMG_`.
- OpenAccel `src/solver/linearSolverContext.hpp`: `setupDefaults_`, `setOptions_`.
- MARS `backend/distributed/unstructured/fem/segregated/mars_segregated_simple_distributed.hpp`:
  `HypreSimpleSolve` constructor and `set_tolerances`.
- MARS `backend/distributed/unstructured/solvers/mars_hypre_gmres_solver.hpp`:
  AMG and Krylov setup, environment parsing and independent residual checks.

Reference option names inside solver/preconditioner blocks are case insensitive,
but the required `type` key must use that spelling. Unsupported entries can be
silently ignored by the reference; the report flags their presence without
disclosing their names. For example, writing `pmaxelmts` in that reference block
does not call `HYPRE_BoomerAMGSetPMaxElmts` in the inspected source.

An explicit MARS setting and an omitted reference setting remain **unresolved**,
even when a familiar Hypre default seems likely. Two omitted settings also remain
unresolved because library versions and builds can differ. Malformed numeric
environment overrides are not certified. MGR is reported as a different
preconditioner; its internal options are outside this audit.

The checks bind the saved deck and case to `pair.json`, compare start/end launch
metadata and the exact solver argument suffix, and verify log/exit hashes.
Executable/library hashes are checked for consistency in the saved records, not
against present files. Mesh/field provenance is intentionally not checked here;
the original startup/pressure comparison retains that responsibility. Launcher
scripts can alter the child environment, and runtime-effective settings are not
queried. The output never claims identical linear solvers or attributes the
observed field difference to a setting merely because that setting differs.

## Inspect the existing pressure pair

Run once in either cluster terminal with the working Python environment. Only
PyYAML is required. No C++ rebuild, environment switch or allocation is needed.
Do not run simultaneous Git updates in both terminals.

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
repo="$scratch/git/mars-v010-check"
test "$(git -C "$repo" branch --show-current)" = cstone
git -C "$repo" fetch origin cstone
git -C "$repo" merge --ff-only refs/remotes/origin/cstone
python3 -c 'import yaml'
umask 077
mkdir -p "$scratch/tmp"
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
run=$(mktemp -d "$scratch/simple-pressure-settings-XXXXXX")
status=0
python3 "$repo/scripts/simple_pressure_settings.py" \
  --pair "$scratch/simple-pressure-accuracy-W3SoRZ/pair" \
  --output "$run/public.json" || status=$?
cat "$run/public.json"
printf 'Share only: %s\n' "$run/public.json"
exit "$status"
)
```

Use the fixed setting differences to choose the next controlled comparison.
Keep experimental recovery disabled and retain the original residual checks.
Matching these source settings alone will not establish history parity or pump
convergence; those remain numerical validation steps.
