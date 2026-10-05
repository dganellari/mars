# Independent geometry and outlet probes

The straight public channel establishes agreement for its geometry and controls.
The bent-inlet host/MPI gates compare MARS implementations with each other; they
are not an independent OpenAccel check. These four short references separate the
geometry and outlet choices:

| Case | Mesh | Outlet |
| --- | --- | --- |
| straight-average | Original pinned public channel | Average pressure, beta 0.05 |
| straight-static | Same channel | Constant static pressure, beta 1 |
| warped-average | Warped, rotated public channel | Average pressure, beta 0.05 |
| warped-static | Same warped channel | Constant static pressure, beta 1 |

All use high-resolution advection, `linear_linear` velocity interpolation,
rho=1 kg/m³, mu=0.1 Pa s, inward speed 0.1 m/s, outlet pressure 0 Pa, and zero
initial velocity/pressure. Relaxation and pseudo-time remain 0.3/0.3/0.75 and
0.01 s. Each runs exactly 20 outer iterations and saves every state. This is a
matched-iteration test, not a convergence test.

For the public coordinates (x,y,z), the warp is
`(x + 0.2*y*y + 0.1*z*z, y + 0.05*x*x, z)`, followed by the proper rotation
`[[.8,-.6,0],[.48,.64,-.6],[.36,.48,.8]]`. Coefficients multiplying squared
coordinates have units of inverse metres. Connectivity, node IDs and side sets
are unchanged. The runner rejects any inverted or collapsed tetrahedron.

`scripts/openaccel_simple_geometry_probe.py` accepts only the checksum-pinned
public channel and deck from the completed public capture, and requires that
capture's exact executable. It derives the four cases itself; it accepts no
user-supplied replacement geometry or physical parameters. Legacy instrumentation
is disabled, so the changed cases cannot be mistaken for frozen captures of the
original fixture. Each `probe.json` records input, executable and output hashes,
the command, revision, iteration coverage and exit status.

## Daint reference command

Use the existing **OpenAccel environment** that ran the previous public
reference. Keep the MARS build environment in its separate terminal. No build is
needed. The library check must succeed before launching the single-rank reference
job. The Python environment needs numpy and netCDF4.

```bash
(
set -euo pipefail
root=/capstor/scratch/cscs/gandanie/git
git -C "$root/mars-v010-check" pull --ff-only
python3 -c 'import numpy, netCDF4'
exe="$root/OpenAccel-reference-updates-IxgJIp/source/build/openaccel-3D.exe"
links=$(ldd "$exe")
if printf '%s\n' "$links" | grep -q 'not found'; then
  printf '%s\n' "$links"
  echo 'Restore the OpenAccel runtime environment before launching.'
  exit 1
fi
run=$(mktemp -d /capstor/scratch/cscs/gandanie/simple-geometry-reference-XXXXXX)
mkdir "$run/tmp"
export TMPDIR="$run/tmp"
printf 'Public results: %s\n' "$run"
git -C "$root/mars-v010-check" rev-parse HEAD > "$run/mars-revision.txt"
printf '%s\n' "$links" > "$run/reference-libraries.txt"
srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=1 \
  --cpus-per-task=1 --cpu-bind=cores --export=ALL --kill-on-bad-exit=1 \
  env OMP_NUM_THREADS=1 OMP_PROC_BIND=close OMP_PLACES=cores \
  python3 "$root/mars-v010-check/scripts/openaccel_simple_geometry_probe.py" \
  --capture "$root/OpenAccel-reference-updates-IxgJIp/run" \
  --executable "$exe" --output "$run/reference" \
  2>&1 | tee "$run/launch.log"
printf 'Public results: %s\n' "$run"
)
```

The four result directories contain public synthetic data only. Comparing every
saved state against the same MARS kernels locates the earliest discrepancy. A
failure restricted to the static cases points toward outlet semantics; one
restricted to warped cases points toward geometry-sensitive operations. These
are localization clues, not proofs of cause. Passing all four leaves other
controls, longer histories and reference-run provenance to investigate.

## Current evidence

Ten local harness checks pass, including unknown input rejection, pinned binary
identity, independent outlet/geometry controls, launch failure and missing states.
Generation against the actual pinned public mesh preserves every non-coordinate
variable. The minimum positive element determinant is 0.015625 m³ before the warp
and 0.013505859375 m³ afterward. Neither determinant is the element volume, which
is one sixth of it.

All four cases also complete 20 local steps through the production arithmetic in
the host `SimpleRunner` with its test-only dense linear solver. The control case
matches all first 20 states of the existing independent public OpenAccel run:
maximum velocity difference/U is 6.46e-13, and maximum absolute pressure
difference/(rho U²) is 1.04e-10. This checks the local comparison path; it does
not establish agreement for the other three cases.

## Independent reference results, 2026-10-01

The user completed all four references on Daint in
`/capstor/scratch/cscs/gandanie/simple-geometry-reference-MUYE1J/reference`,
using harness revision `ac63b2c27e5c4764c8b5dc0e3158a5f9cc8fcfaf`.
Retrieved manifests, decks, meshes, logs, exits and field files were checked:
all four exited zero, saved iterations 1 through 20, and match their recorded
SHA-256 hashes. Local and reference decks agree exactly; connectivity, side sets
and node IDs agree, and coordinate differences are below 2e-15 m.

Comparison by source node ID covers every node at every saved iteration. Velocity
uses the vector difference divided by U=0.1 m/s; pressure uses the absolute
difference divided by rho*U²=0.01 Pa. No pressure offset is removed.

| Case | Maximum velocity difference / U | Maximum pressure difference / (rho U²) |
| --- | --- | --- |
| straight-average | 6.46e-13 | 1.04e-10 |
| straight-static | 5.91e-13 | 8.81e-11 |
| warped-average | 7.84e-13 | 1.96e-10 |
| warped-static | 7.32e-13 | 1.96e-10 |

These are maxima across all 20 states. They validate the host production
arithmetic against the independent reference for these cases, including the
nonplanar inlet and both outlet types. They do not identify the remaining private
snapshot discrepancy. No production kernel, solver tolerance or private input
changed. CUDA ingestion and MPI execution of these four references were
unverified at this host milestone; the GPU results below extend the evidence to
the production driver at iteration 20.

## Production GPU results, 2026-10-05

All four cases pass at iteration 20 on one and four GPUs. The one-GPU
straight-average result is in `simple-geometry-mars-fUCfA0`; the other seven
are in `simple-geometry-retry-FJGrnL`, both under
`/capstor/scratch/cscs/gandanie/`. They use the same executable SHA-256
`385c99163d1930c6462ed38090b358f4d66b232efefe76db3848ecea22052fea`,
with recorded source revision `ac63b2c27e5c4764c8b5dc0e3158a5f9cc8fcfaf`.

The saved OpenAccel comparison reports give the following maximum nodal
differences, using the velocity and pressure scales defined above. All 64
input hashes in those reports match the retrieved files, including the reference
meshes, decks and Exodus results. Each completed MARS run exits 2 and records
iteration 20; this expected exit is not a convergence failure for this test.

| Case | GPUs | Velocity difference / U | Pressure difference / (rho U²) |
| --- | --- | --- | --- |
| straight-average | 1 | 5.39e-13 | 3.67e-11 |
| straight-average | 4 | 5.64e-13 | 4.25e-11 |
| straight-static | 1 | 2.90e-13 | 2.70e-11 |
| straight-static | 4 | 3.19e-13 | 3.16e-11 |
| warped-average | 1 | 1.07e-12 | 2.05e-10 |
| warped-average | 4 | 1.08e-12 | 2.12e-10 |
| warped-static | 1 | 5.72e-13 | 5.78e-11 |
| warped-static | 4 | 5.79e-13 | 5.68e-11 |

An independent comparison of the retrieved one- and four-GPU CSV fields covers
all 425 nodes per case, with identical coordinates. Across the four cases, the
largest rank differences are 6.92e-14 for velocity/U and 3.26e-11 for pressure/(rho U²).
These results validate the native ingestion and distributed SIMPLE path for
these four public cases. They do not establish nonlinear convergence, performance
scaling or agreement on a private pump mesh; the remaining private mismatch is
not explained by this test.

## Production CUDA/MPI comparison

Run this from `mars-v010-check/build-hypre` in the **MARS uenv**. The OpenAccel
environment is not needed. Reuse the current executable containing the inlet and
outlet repairs; this documentation change needs no rebuild. The Python environment
needs numpy, netCDF4 and yaml. Outputs and temporary files stay on capstor.

The one- and four-rank runs use the exact saved reference decks and meshes.
Exit 2 is accepted only as a completed iteration-limited run; the comparator must
independently verify completion at iteration 20 and field agreement. This is not
a nonlinear convergence test or a scaling measurement. The existing comparison
script calls its detailed output a private report, but all data in this recipe
are synthetic public data and the full report may be shared.

These deliberately short MARS runs use `--kill-on-bad-exit=0` so every rank can
finish and return the expected exit code 2. With `--kill-on-bad-exit=1`, Slurm
can terminate a remaining rank during shutdown and return 143 instead. The
exit-status and field checks below still reject unsuccessful comparisons.

```bash
(
set -euo pipefail
python3 -c 'import numpy, netCDF4, yaml'
exe=./examples/distributed/unstructured/mars_segregated_simple
test -x "$exe"
links=$(ldd "$exe")
if printf '%s\n' "$links" | grep -q 'not found'; then
  printf '%s\n' "$links"
  echo 'Restore the MARS runtime environment first.'
  exit 1
fi
reference=/capstor/scratch/cscs/gandanie/simple-geometry-reference-MUYE1J/reference
run=$(mktemp -d /capstor/scratch/cscs/gandanie/simple-geometry-mars-XXXXXX)
mkdir "$run/tmp"
export TMPDIR="$run/tmp"
printf 'Public results: %s\n' "$run"
git rev-parse HEAD > "$run/mars-revision.txt"
sha256sum "$exe" > "$run/executable.sha256"
printf '%s\n' "$links" > "$run/libraries.txt"
unset CUDA_LAUNCH_BLOCKING MARS_OWNERSHIP MARS_HALO_FACTOR MARS_NODEHALO_ALLOW_INCONSISTENT
export MARS_HYPRE_SPMV_VENDOR=0 MARS_HYPRE_FLEXGMRES=0 MARS_HYPRE_ABSTOL=0
export MARS_HYPRE_VERBOSE=0 MARS_HYPRE_RESIDUAL_AUDIT=0 MARS_SIMPLE_PRESSURE_AUDIT=0
for case in straight-average straight-static warped-average warped-static; do
  ref="$reference/$case"
  python3 ../scripts/prepare_simple_deck.py --deck "$ref/input.i" \
    --mesh "$ref/channel.exo" --output "$run/$case/case"
  mapfile -d '' -t args < "$run/$case/case/args.nul"
  for np in 1 4; do
    out="$run/$case/np$np"
    mkdir "$out"
    set +e
    srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node="$np" \
      --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=0 \
      ~/affinity/bind_numa.sh "$exe" "${args[@]}" \
      --output-prefix "$out/channel" --iterations 20 --report-every 10 \
      --field-output gathered --residual-tol 1e-6 --mass-tol 1e-6 --change-tol 1e-6 \
      2>&1 | tee "$out/run.log"
    statuses=("${PIPESTATUS[@]}")
    set -e
    printf '%s\n' "${statuses[0]}" > "$out/run.exit"
    if (( statuses[1] != 0 || (statuses[0] != 0 && statuses[0] != 2) )); then exit 1; fi
    python3 ../scripts/simple_snapshot_compare.py \
      --reference-dir "$ref" --case "$run/$case/case/case.json" --iteration 20 \
      --mars-prefix "$out/channel" --mars-log "$out/run.log" --mars-exit "$out/run.exit" \
      --private-report "$out/field-errors.json" --output "$out/summary.json"
    python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); print(d["errors"]); s=json.load(open(sys.argv[2])); assert s["snapshot_fields_within_tolerance"] and s["reference_settings_status"] == "mapped_controls_match"' \
      "$out/field-errors.json" "$out/summary.json"
  done
done
printf 'Public results: %s\n' "$run"
)
```
