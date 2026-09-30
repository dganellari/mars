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
requires **both maximum errors <= 1e-5** in those scales. The peak-speed band
instead uses the reference peak as denominator and reports `zero_reference_peak`
when it is zero. Equal peaks do not imply equal velocity fields.

`reference_settings_status=mapped_controls_match` checks only the subset supported
by the deck bridge, not the reference launch or all solver settings. Linear
solvers/preconditioners and nonlinear norms still differ; therefore
`identical_linear_solvers_verified` stays false. Same-iteration mismatch on an
unconverged calculation does not by itself identify an implementation bug: the
nonlinear trajectories can differ. Inspect the private errors before changing
relaxation, physical timescale or boundary conditions.

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

## Local synthetic checks

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=scripts \
  python3 -m unittest test_simple_snapshot_compare test_simple_convergence_summary \
    test_simple_public_diagnostics test_prepare_simple_deck
```

These test synthetic Exodus shards and field files, large global IDs, reordering,
identity, missing/duplicate data, private-error redaction, evidence mismatches and
unchanged absolute-pressure semantics. They are not results from a private mesh.
