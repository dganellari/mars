# Compare saved steady SIMPLE snapshots privately

`scripts/simple_snapshot_compare.py` runs on the user's machine holding both
results. It reads existing files; it launches no simulation and changes no solver
or tolerance. It accepts a completed MARS run that converged **or reached its
iteration limit**. A stalled reference can therefore be compared without first
changing the calculation to force convergence.

The detailed JSON contains private numerical results, filenames and hashes. Only
the separate public JSON may be shared. It contains fixed labels and booleans:
no geometry, node counts, boundary names, input hashes, iterations or physical
values. Exceptions are redacted, including errors while reading a private deck.
New report files use owner-only permissions and never overwrite an existing file.

## Inputs and evidence

- A native MARS prefix with `-metrics.csv` and either `-fields.csv` or a
  `mars-simple-fields-v1` distributed `-fields.json` manifest and all its parts.
- Its completed log and exit file. The log reports must match the metrics, the
  rank count must be consistent, and the final iteration must equal `--iteration`.
  Exit 0 means convergence; exit 2 must have an iteration-limit completion.
  Supply the original nonlinear targets if they differ from the defaults 1e-6.
- The `case.json` from `prepare_simple_deck.py` used for that MARS run, and the
  native input Exodus mesh. Its path defaults to the saved preparation's path;
  an explicit `--mesh` must resolve to that same file. The printed physical/relaxation/interpolation
  controls and explicit pressure targets must agree with the prepared arguments.
- An OpenAccel result directory containing one `results.e*` family. Every MPI
  piece is required; mixed families are rejected. Each piece must contain the
  requested `time_whole` coordinate and the same strictly increasing saved
  coordinates. For the supported **steady** deck this selects an iteration,
  not elapsed physical time. This tool does not interpolate snapshots.
- Optional saved reference deck in that directory (`.yaml`, `.yml`, `.i`). An
  exact match to the preparation's deck hash is preferred. Otherwise a sole
  candidate is translated and compared by supported controls. Missing,
  ambiguous, unsupported or different decks are reported explicitly; field
  differences are still computed, without claiming matching reference settings.

Both Exodus layouts (combined and separate nodal variables) are supported, with
`velocity_x`, `velocity_y`, `velocity_z`, `pressure`. MARS source-row IDs map
through the input mesh's `node_num_map` (or Exodus's implicit 1-based numbering).
Partitioned reference outputs must carry their global node maps. Every source
node must occur, MARS rows must occur exactly once, and reference ghost values
must agree to roundoff. Coordinates are checked at
`64 * epsilon * max(1 metre, max |coordinate|)`. No coordinate fitting or nearest
neighbour matching is performed.

These checks bind the selected fields to the input's **nodes**, not its element
connectivity. Existing outputs have no runtime signature binding fields, boundary
selection and executable to a launch; `full_run_provenance_verified` stays false.
The final MARS peak is cross-checked against its metrics, but this alone cannot
detect all swapped field files. The user must select the original run artifacts.
Input hashes in the private report record what was compared, not when it ran.

## Interpretation

`comparison_status=completed` means the comparison executed, not that the fields
agree or either nonlinear solve converged. `invalid_evidence` gives a fixed
`failed_check` label. It never treats missing evidence as a numerical pass.

Velocity errors use the nodal **vector** difference divided by the prepared inlet
speed U. Pressure errors use absolute pressure difference divided by rho U².
Reported RMS is an unweighted nodal RMS, not a volume integral. No mean pressure
shift is subtracted. Private output also records the signed mean difference.

Each error gets the tightest band it meets: `within_1e_minus_5`,
`within_1_percent`, `within_5_percent`, or `over_5_percent`. The last three are
diagnostic bands, not CFD accuracy criteria. `snapshot_fields_within_tolerance`
requires **both maximum errors <= 1e-5** in those scales and a hash-matched saved
deck requesting raw solver values, as described below. The peak-speed band
instead uses the reference peak as denominator and reports `zero_reference_peak`
when it is zero. Equal peaks do not imply equal velocity fields.

`reference_settings_status=mapped_controls_match` checks only the subset supported
by the deck bridge, not the reference launch or all solver settings. Linear
solvers/preconditioners and nonlinear norms still differ; therefore
`identical_linear_solvers_verified` stays false. Same-iteration mismatch on an
unconverged calculation does not by itself identify an implementation bug: the
nonlinear trajectories can differ. Inspect the private errors before changing
relaxation, physical timescale or boundary conditions.

## Boundary output and localization

OpenAccel `0d69041` has a separate **output** setting:
`simulation.solver.output_control.corrected_boundary_values` (default `false`).
In `src/simulation/simulationIO.cpp`, `writeResults()` temporarily replaces
boundary nodal values with their boundary-field values, writes Exodus, and
restores the solver values. Comparing that corrected output to MARS's raw nodal
fields is not a solver-state comparison. This setting is independent of the
physical controls already checked by the preparation tool.

The public summary now reports `reference_deck_output_values` as `solver_values`,
`boundary_corrected` or `unknown`. `solver_field_comparison_supported` requires
`solver_values` and an exact match to the preparation's deck hash. A corrected,
missing, malformed or unverified deck cannot produce a snapshot parity pass,
even if the selected numerical fields happen to agree. Numerical error bands
are still reported. The saved deck is evidence of requested settings, not proof
of the actual reference launch; full provenance remains unverified.

To localize a mismatch, the comparator also forms the union of the input mesh's
stored side-set and node-set nodes. Tet4 side sets use Exodus element-row and
side numbering across blocks, not global node or element IDs. The public report
contains the same four error bands separately for `tagged_boundary` and
`other_nodes`. It exports no names, counts, coordinates or numerical errors.
Unsupported side topology or absent tags gives `boundary_localization_status`
`unavailable`; malformed indexing fails the evidence check. Empty groups are
reported explicitly rather than treated as an agreement.

These labels do not certify complete physical boundary coverage: extra node sets
can include interior nodes, and missing tags can omit boundary nodes. Consequently
`other_nodes` is not automatically a certified interior. Regional agreement never
replaces the full-field check or excuses an error. With complete boundary tags,
an error away from the boundary cannot be explained solely by output correction.
No pressure offset or boundary values are fitted or removed.

Rerun the comparison on existing results to obtain these fields. No rebuild,
solver launch, environment change or tolerance change is needed. This closes an
unchecked comparison assumption; it does not establish the cause of a particular
private mismatch.

## User-side command

Requires Python 3.6+, numpy, netCDF4 and PyYAML. All arrays are post-processing
data on the host; the solver is not involved. Memory grows linearly with node
count. CSV parts and reference files are read sequentially; Exodus fields use
bounded chunks. Do not upload any of these inputs or the private report.

From the MARS checkout, with variables pointing at the saved artifacts:

```bash
(
set -euo pipefail
umask 077
: "${MARS_RUN:?}" "${REFERENCE_RUN:?}" "${CASE_JSON:?}" "${MESH_BIG:?}"
report_dir=$(mktemp -d /capstor/scratch/cscs/gandanie/simple-comparison-XXXXXX)
printf 'Comparison directory: %s\n' "$report_dir"
python3 scripts/simple_snapshot_compare.py \
  --reference-dir "$REFERENCE_RUN" --mesh "$MESH_BIG" --case "$CASE_JSON" \
  --mars-prefix "$MARS_RUN/flow" --mars-log "$MARS_RUN/run.log" \
  --mars-exit "$MARS_RUN/run.exit" --iteration 2000 \
  --private-report "$report_dir/private.json" --output "$report_dir/public.json"
cat "$report_dir/public.json"
)
```

Change the prefix and iteration only to those of the saved runs being compared.
The public JSON is also written on an evidence failure (exit 1); inspect that
file after failure. Exit 0 means export completed even if fields differ.

If dependencies are absent, install them into a **scratch** directory, leaving
home unchanged, then rerun the command with that directory on `PYTHONPATH`:

```bash
deps=$(mktemp -d /capstor/scratch/cscs/gandanie/simple-compare-python-XXXXXX)
python3 -m pip install --no-cache-dir --only-binary=:all: --target "$deps" numpy netCDF4 PyYAML
export PYTHONPATH="$deps${PYTHONPATH:+:$PYTHONPATH}"
```

## Localize disagreement at an early saved state

`scripts/prepare_simple_snapshot_probe.py` prepares a short run; it does not
launch one. Use it after the late-state comparison shows disagreement. It requires
the completed baseline directory with `run.log`, `run.exit`, `flow-metrics.csv`
and `executable.sha256`, plus the original preparation's `case.json` and `args.nul`.

By default the executable must match the baseline SHA-256. The original arguments must match
the preparation and printed controls, and the reference deck hash must still
match. Every reference shard must have the same saved coordinates. The probe
selects the **earliest positive** coordinate, requires an integer steady iteration
strictly before the baseline's final iteration, and refuses anything beyond 100
iterations (`--max-iteration` can lower this cap). Initialization is not selected.

The generated arguments preserve physical controls, explicit pressure targets,
cache/halo settings and field-output mode. The recipe preserves rank count, starts
from the same zero initialization, and changes only the iteration cap and report
frequency. Supply the baseline's nonlinear targets if they differ from 1e-6.
Only one-node baselines with 1–4 ranks and profiling disabled are accepted. The
existing executable is reused; a source pull does not require a solver rebuild.

Run preparation on the machine holding the data. With `baseline`, `case_file`,
`reference`, `exe` and a fresh private `run` directory set to the relevant paths:

```bash
python3 ../scripts/prepare_simple_snapshot_probe.py \
  --baseline "$baseline" --case "$case_file" --reference-dir "$reference" \
  --executable "$exe" --output-dir "$run/probe" --output "$run/preparation-public.json"
```

If preparation succeeds, `probe/args.nul`, `iteration.txt`, `ranks.txt` and
`probe.json` are private launch metadata. If it fails, share only
`preparation-public.json`; do not bypass stale arguments.
The one-node Alps launch, from the configured MARS build directory, is:

```bash
mapfile -d '' -t args < "$run/probe/args.nul"
iteration=$(cat "$run/probe/iteration.txt")
np=$(cat "$run/probe/ranks.txt")
unset CUDA_LAUNCH_BLOCKING MARS_OWNERSHIP MARS_HALO_FACTOR MARS_NODEHALO_ALLOW_INCONSISTENT
export MARS_HYPRE_SPMV_VENDOR=0 MARS_HYPRE_FLEXGMRES=0 MARS_HYPRE_ABSTOL=0
export MARS_HYPRE_VERBOSE=0 MARS_HYPRE_RESIDUAL_AUDIT=0 MARS_SIMPLE_PRESSURE_AUDIT=1
python3 ../scripts/run_simple_diagnostics.py \
  --log "$run/run.log" --exit-file "$run/run.exit" --output "$run/run-public.json" -- \
srun --account=csstaff --time=00:15:00 --nodes=1 --ntasks-per-node="$np" \
  --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
  "$HOME/affinity/bind_numa.sh" "$exe" "${args[@]}" --output-prefix "$run/flow"
```

That environment matches the recorded baseline recipe for the current comparison;
it is not a general proof of shared-library or unrecorded environment identity.
Tracing is omitted; do not use this check as a performance benchmark. Raw stdout
and stderr remain in the private log. Handle the expected iteration-limit exit 2
before proceeding in a `set -e` shell; other failures require the public run
diagnostic, not a field comparison.

Then run `simple_snapshot_compare.py` with the same `reference` and `case_file`,
`--iteration "$iteration"` and the new `flow` prefix/log/exit. An early nonlinear
convergence before that selected coordinate is rejected as an iteration mismatch;
the script does not silently change stopping criteria to reach it.

Early agreement and late disagreement locate growth between those checkpoints.
An early mismatch means disagreement is already present **by the first saved
reference state**, which need not be the first iteration. Neither outcome alone
identifies whether assembly, boundary updates or differing linear solves caused
it. Keep solver controls unchanged until that distinction is investigated.

### Comparing an implementation repair

After intentionally rebuilding a changed solver, add `--allow-executable-change`
to the preparation command. The default still rejects a changed executable. This
opt-in preserves all baseline, argument, deck, runtime-option and iteration checks;
it only permits the binary hash to differ. Both hashes are recorded privately in
`probe.json`. The public `binary_matches_baseline` remains false for a changed
binary and `executable_change_allowed` records the opt-in. Do not replace the
baseline hash file or describe this comparison as a same-binary experiment.

The [nodal inlet correction](inlet_normal_speed.md) needs this opt-in for its short
comparison after the public GPU checks. Use the rebuilt executable, the original
case and reference, and a fresh output directory. No OpenAccel rebuild is needed.

## Local synthetic checks

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=scripts \
  python3 -m unittest test_simple_snapshot_compare test_simple_convergence_summary \
    test_simple_public_diagnostics test_prepare_simple_deck test_prepare_simple_snapshot_probe
```

These test synthetic Exodus shards and field files, large global IDs, reordering,
identity, missing/duplicate data, private-error redaction, evidence mismatches and
unchanged absolute-pressure semantics. They are not results from a private mesh.
