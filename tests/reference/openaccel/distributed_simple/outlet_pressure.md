# Outlet pressure at shared nodes

The pressure correction and boundary mass-flux update use an area-weighted nodal
boundary pressure. They must not use the face trace directly when neighbouring
faces carry different values. This matters after outlet reversal closes a face.

The pinned public OpenAccel source is `0d69041ba1afda63e9e4328d9e0d9834bba37756`:

- `flowModel::updatePressureBoundarySideFieldAverageStaticPressure_` refreshes
  open-face traces and retains closed-face traces, then calls
  `nodeSideField<scalar,1>::interpolate`.
- That interpolation averages **all** adjacent outlet samples onto their nearest
  nodes, including closed faces. It uses sample-area magnitudes as weights.
- `pressureCorrectionAssembler::assembleElemTermsBoundaryOutletSpecifiedPressure_`
  and `flowModel::updateMassFlowRateBoundaryFieldOutletSpecifiedPressure_` read
  this nodal field when forming the boundary pressure gradient.
- The reversal decision still reads the original face trace.

For outlet samples attached to node `i`, with face pressure `t_fi` in Pa and sample
area `a_fi` in square metres,

```
p_boundary[i] = sum_f(a_fi * t_fi) / sum_f(a_fi).
```

Closed faces are excluded from the open-outlet mean used to refresh `t_fi`, but
are included in this separate interpolation. Before closure, all incident traces
at a node agree, so using the face value directly hides the missing operation.

MARS retains both representations. `SimpleOutletTraceSum` accumulates the area
once at setup and the weighted pressure each iteration. `SimpleOutletTraceFinish`
normalizes it. Both run on the GPU in production. Complete element stars provide
all incident boundary faces at owned nodes; owners publish the resulting scalar
through the existing device halo. Ghost partial sums are never used. This adds one
scalar per ghost and one exchange round per iteration, overlapping interior
pressure assembly when `--halo-overlap 1`. No field host staging is added.

## Regression evidence, 2026-09-30

The independent reference uses the pinned public 425-node channel, high-resolution
advection and `linear_linear` velocity interpolation. Its manifest reports 1431
iterations, with mesh SHA-256
`3b3a2d8101291955d496f9c696ce250cb5acc5c2e1e34eae05c9db5ae5f861fd`
and result SHA-256
`e94fa17fe753f679bb6cddcde3a4852febd94e7f63e7c5c0b5c7880578af685f`.
The files were retrieved and both hashes checked before comparison.

Running the production kernels through the local host `SimpleRunner` on that
mesh isolates the first discrepancy: old MARS agrees through iteration 53, when
two outlet faces close. At iteration 54 its maximum absolute velocity difference
is `1.009e-4 m/s` and pressure difference is `6.705e-3 Pa`. With nodal averaging,
those differences fall to `1.595e-13 m/s` and `1.189e-12 Pa`.

At iteration 1431, maximum vector velocity difference divided by `U=0.1 m/s` is
`1.485e-13`; maximum absolute pressure difference divided by `rho*U^2=0.01 Pa` is
`1.066e-11`. No pressure offset is removed. This is a host-kernel comparison with
an independently executed reference. The separate GPU results follow below.

The algebra regression uses unequal sample areas, a retained closed-face trace,
and an updated open-face trace. It checks the shared nodal values and their use in
both pressure assembly and flux updates. Distributed gates compare the nodal
boundary pressure explicitly and poison unexchanged ghost entries.
The final local ASan/UBSan build passes 2139 algebra checks and 15 shifted/sheared
host MPI tests, covering 1/2/4 ranks and overlapping/synchronous assembly. Bypassing
the nodal pressure in the boundary kernel makes the algebra regression fail.

## Daint CUDA/Hypre validation, 2026-09-30

The user ran revision `02459212cd81c0497e3b76d96a69b894c273d7b3` on 1/2/4 GPUs
with Hypre 2.33.0 and native GPU SpMV. Public results are saved at
`/capstor/scratch/cscs/gandanie/simple-outlet-pressure-OVYNJC`.
The recorded executable SHA-256 is
`385c99163d1930c6462ed38090b358f4d66b232efefe76db3848ecea22052fea`.
Retrieved fields, metrics, logs and exit files match their comparison hashes;
recomputing differences against the pinned reference reproduces all three reports.

| GPUs | Max velocity difference / U | Max pressure difference / (rho U^2) | Exit |
| --- | --- | --- | --- |
| 1 | 1.746e-13 | 1.101e-11 | 2 |
| 2 | 1.608e-13 | 9.770e-12 | 2 |
| 4 | 1.858e-13 | 1.208e-11 | 2 |

Every run reaches iteration 1431 and passes the matched-snapshot field tolerance
of `1e-5`, without removing a pressure offset. The four-rank log records SFC
ownership and element-star completion. This validates the repaired CUDA path and
rank agreement on this public case; full runtime provenance and identical linear
solvers are not attested by the comparator.

Iteration wall time was 30.84 seconds on one GPU (`nid005365`), versus
69.18/80.22 seconds on 2/4 GPUs (`nid006489`). These individual runs of a 425-node
test show no speedup and do not establish scalability.

## Conservation is a separate result

This shifted reference reports native convergence despite a roughly 62.96% mass
imbalance. The corrected local run reproduces it. At iteration 1431 the raw
Rhie--Chow outlet flux sums to `0.0999999999741 kg/s`. Relaxation against the stored
history gives `0.115739259302 kg/s`; clipping negative outlet samples adds
`0.0472177779572 kg/s`, leaving `0.162957037260 kg/s`, versus inlet `-0.1 kg/s`.

The mass discrepancy is therefore reproduced by the reference's flux processing;
it is not cured by this pressure-interpolation repair. MARS retains its continuity
and mass-balance convergence checks and must report an iteration limit for this
case. Compare matched snapshots; do not use the converged-field checker or claim
physical convergence, private-pump parity or general scaling from this result.

All three GPU runs retain the iteration-limit status and reproduce the roughly
62.96% imbalance. Existing reference binaries and results can be reused;
OpenAccel need not rebuild.

## Interactive Daint verification

Use the existing **MARS uenv**, not the OpenAccel session. This reuses the saved
public reference and runs MARS on 1/2/4 GPUs. Exit 2 is expected from the mass
balance gate; the separate snapshot comparison must pass at iteration 1431.
All generated files stay on capstor scratch.

```bash
(
set -euo pipefail
cd /capstor/scratch/cscs/gandanie/git/mars-v010-check/build-hypre
git pull --ff-only
cmake --build . --parallel 4 --target mars_segregated_simple mars_segregated_simple_algebra_check
./examples/distributed/unstructured/mars_segregated_simple_algebra_check
python3 -c 'import numpy, netCDF4, yaml'
exe=./examples/distributed/unstructured/mars_segregated_simple
links=$(ldd "$exe")
if printf '%s\n' "$links" | grep -q 'not found'; then
  printf '%s\n' "$links"
  echo 'Stop: use the MARS uenv shell.'
  exit 1
fi
reference=/capstor/scratch/cscs/gandanie/simple-shifted-bWj8Qu/openaccel
run=$(mktemp -d /capstor/scratch/cscs/gandanie/simple-outlet-pressure-XXXXXX)
mkdir "$run/tmp"
export TMPDIR="$run/tmp"
printf 'Public results: %s\n' "$run"
git rev-parse HEAD > "$run/mars-revision.txt"
sha256sum "$exe" > "$run/executable.sha256"
python3 ../scripts/prepare_simple_deck.py --deck "$reference/input.i" \
  --mesh "$reference/channel.exo" --output "$run/case"
mapfile -d '' -t args < "$run/case/args.nul"
unset CUDA_LAUNCH_BLOCKING MARS_OWNERSHIP MARS_HALO_FACTOR MARS_NODEHALO_ALLOW_INCONSISTENT
export MARS_HYPRE_SPMV_VENDOR=0 MARS_HYPRE_FLEXGMRES=0 MARS_HYPRE_ABSTOL=0
export MARS_HYPRE_VERBOSE=0 MARS_HYPRE_RESIDUAL_AUDIT=0 MARS_SIMPLE_PRESSURE_AUDIT=0
for np in 1 2 4; do
  out="$run/np$np"
  mkdir "$out"
  set +e
  srun --account=csstaff --time=00:20:00 --nodes=1 --ntasks-per-node="$np" \
    --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
    ~/affinity/bind_numa.sh "$exe" "${args[@]}" \
    --output-prefix "$out/channel" --iterations 1431 --report-every 100 \
    --field-output gathered --residual-tol 1e-6 --mass-tol 1e-6 --change-tol 1e-6 \
    2>&1 | tee "$out/run.log"
  statuses=("${PIPESTATUS[@]}")
  set -e
  printf '%s\n' "${statuses[0]}" > "$out/run.exit"
  if (( statuses[1] != 0 || (statuses[0] != 0 && statuses[0] != 2) )); then exit 1; fi
  python3 ../scripts/simple_snapshot_compare.py \
    --reference-dir "$reference" --case "$run/case/case.json" --iteration 1431 \
    --mars-prefix "$out/channel" --mars-log "$out/run.log" --mars-exit "$out/run.exit" \
    --private-report "$out/field-errors.json" --output "$out/summary.json"
  python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); print(d["errors"]); s=json.load(open(sys.argv[2])); sys.exit(0 if s["snapshot_fields_within_tolerance"] else 1)' \
    "$out/field-errors.json" "$out/summary.json"
done
printf 'Public results: %s\n' "$run"
)
```
