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

## Read defaults from the captured reference library

`scripts/simple_hypre_defaults.py` builds a small parameter probe against the
Hypre shared library recorded in the saved OpenAccel capture. It locates the
installed headers beside the resolved library, checks the current executable
and library hashes against that capture, and checks the probe's Hypre/MPI links.
The running probe reports its actual Hypre and MPI library identities through
`dladdr`; these must also match. A failed build or launch produces a failed
public report, even if some output exists.

The probe creates fresh GMRES, FlexGMRES and BoomerAMG objects. It reads restart
dimensions, minimum iterations, coarsening, interpolation, relaxation, sweeps,
coarse-grid limits and the other defaults left unresolved above. It creates no
matrix and performs no setup or solve. Header/runtime versions and available
getter-versus-header values must agree before results are exported. A few AMG
parameters lack getters; those use the matching installation's internal headers.
This checks common layout fields, not a complete ABI proof.

Only `public.json` is shareable. It contains fixed keys, version/build metadata
and **fresh-object library defaults**, never deck values, paths or log fragments.
The build log, launch log and provenance remain in the private output directory.
It does not open mesh, matrix or field files; saved input/log files are hashed
only to verify capture identity. Existing captures and solver settings are not
modified.

These are library initialization defaults, not proof of the original executable's
effective solver settings. Explicit reference setters still take precedence;
matrix-dependent AMG setup choices are outside this probe. In particular, report
down/up/coarse relaxation and sweeps separately: their defaults need not agree.
CPU and GPU execution policies can select different defaults. Matching library
defaults will not certify equal preconditioners or convergence histories.

Run this block **only in the OpenAccel terminal**, with its working MPI/compiler
environment. No MARS or OpenAccel rebuild is needed; the script compiles only the
small probe. Use one task on a compute node because Hypre initialization may use
the GPU in a CUDA build. All files stay in capstor scratch. The default compiler
is `mpicxx`; `--cxx` selects another MPI C++ wrapper, and `--include-dir` can supply
an additional dependency header directory if the installation requires one.

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
run=$(mktemp -d "$scratch/simple-hypre-defaults-XXXXXX")
status=0
python3 "$repo/scripts/simple_hypre_defaults.py" \
  --pair "$scratch/simple-pressure-accuracy-W3SoRZ/pair" \
  --output-dir "$run/probe" \
  --launcher srun --account=csstaff --time=00:03:00 --nodes=1 \
    --ntasks-per-node=1 --export=ALL --kill-on-bad-exit=1 \
    "$HOME/affinity/bind_numa.sh" || status=$?
cat "$run/probe/public.json"
printf 'Share only: %s\n' "$run/probe/public.json"
exit "$status"
)
```

Local validation: the probe compiles and reads defaults with CPU Hypre 2.32.0,
2.33.0 and 3.1.0; synthetic script tests cover stale captures, changed libraries,
MPI mismatches, failed builds/launches and public-output restrictions. These
local results do not establish which defaults the captured Alps library uses;
that requires the command above. The script requires a shared Hypre installation;
it rejects missing or ambiguous library records instead of guessing a prefix.
