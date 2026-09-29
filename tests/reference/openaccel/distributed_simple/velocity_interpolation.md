# Velocity interpolation in Tet4 SIMPLE

`--velocity-interpolation trilinear` is the existing default.
`--velocity-interpolation linear-linear` implements OpenAccel's
`velocity_interpolation_type: linear_linear`. Advection selection is independent;
both modes support upwind and high-resolution advection.

The equations remain steady incompressible laminar momentum and continuity,
with velocity in m/s, physical pressure in Pa, density in kg/m³ and dynamic
viscosity in Pa s. The velocity inlet, pressure outlet, stationary no-slip wall,
pseudo-time and relaxation conventions are unchanged.

## Discrete correspondence

For the median-dual SCS associated with Tet4 edge (L,R), standard field weights
are 13/36 on L,R and 5/36 on the other two nodes. Shifted weights are 1/2 on L,R
and zero elsewhere. On a Tri3 boundary, standard weights are 11/18 at the nearest
node and 7/36 at each other node; shifted sampling uses only the nearest node.

`native_interior` applies this selection to velocity and pressure-influence
interpolation. `native_boundary` applies it to face sampling, including the
no-slip momentum coefficient. Constant density/viscosity and prescribed uniform
boundary data retain their values under either set of weights.

Velocity gradients use the same field weights in their median-dual surface sum:
`grad(u)_n = sum[(u_sample - u_n) A_sample] / V_n`.
The incremental subtraction preserves constants; shifted boundary samples
contribute zero because their sample is exactly the nearest nodal value.
Both the single-rank reference runner and distributed production runner select
these weights before high-resolution limiting. Pressure reconstruction remains
shifted. Tet4 derivative shapes are constant and do not move with sample location.

Coordinate interpolation always uses the standard shapes. Neither area vectors
nor limiter reconstruction points change. The same oriented mass flux still
enters the two incident continuity rows with opposite signs. Outlet reversal
uses face means; equal Tri3 sample areas make both interpolation choices give
the same mean from the same nodal values. Changes in the solved fields can still
change which outlets close.

Public source correspondence at OpenAccel `0d69041`: `nodeField::isShifted`,
`nodeField::updateGradientField`, the segregated momentum/pressure element and
boundary assemblers, and Nalu `TetSCS::shifted_shape_fcn`/Tri3 face shapes.
MARS changes only sampling selectors in existing device kernels; no new halo
messages, field transfers or host numerical work are introduced.

## Verification status

Local checks on 2026-09-29: 1,985 independent SIMPLE algebra checks include shifted
endpoint fluxes, pairwise continuity cancellation, gradient surface sums,
unchanged geometric weights and nodal wall blocks on an oblique tetrahedron.
The focused 24-test host suite passes with ASan/UBSan: CLI validation, existing
standard high-resolution regressions and shifted upwind/high-resolution parity
on 1/2/4 ranks through 12 iterations with oblique inlet and outlet reversal.
The preparer passes 20 tests; the public reference comparator passes 34 tests,
including rejecting missing, duplicated or mismatched interpolation provenance.

These are host tests. CUDA execution, nonlinear convergence and independent
OpenAccel field parity for shifted velocity interpolation remain pending. Earlier
channel/duct results use standard sampling and do not validate this new option.

## Public GPU check

From the existing MARS CUDA/Hypre build on Daint, in Bash:

```bash
set -euo pipefail
cd /capstor/scratch/cscs/gandanie/git/mars-v010-check/build-hypre
git pull --ff-only
cmake --build . --parallel 4 --target mars_segregated_simple
shifted_run=$(mktemp -d /capstor/scratch/cscs/gandanie/simple-shifted-XXXXXX)
export shifted_run
channel_mesh=/capstor/scratch/cscs/gandanie/git/OpenAccel-simple-converged-20260925-151406/channel.exo
printf 'Results: %s\n' "$shifted_run"
git rev-parse HEAD > "$shifted_run/mars-revision.txt"
sha256sum ./examples/distributed/unstructured/mars_segregated_simple > "$shifted_run/executable.sha256"
unset CUDA_LAUNCH_BLOCKING MARS_OWNERSHIP MARS_HALO_FACTOR MARS_NODEHALO_ALLOW_INCONSISTENT
export MARS_HYPRE_SPMV_VENDOR=0 MARS_HYPRE_FLEXGMRES=0 MARS_HYPRE_ABSTOL=0
export MARS_HYPRE_VERBOSE=0 MARS_HYPRE_RESIDUAL_AUDIT=0
for np in 1 2 4; do
  mkdir "$shifted_run/np$np"
  set +e
  srun --account=csstaff --time=00:20:00 --nodes=1 --ntasks-per-node="$np" \
    --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
    ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_segregated_simple \
    --mesh "$channel_mesh" --output-prefix "$shifted_run/np$np/channel" \
    --advection high-resolution --velocity-interpolation linear-linear \
    --iterations 5000 --report-every 100 \
    --residual-tol 1e-6 --mass-tol 1e-6 --change-tol 1e-6 \
    2>&1 | tee "$shifted_run/np$np/run.log"
  statuses=("${PIPESTATUS[@]}")
  set -e
  printf '%s\n' "${statuses[0]}" > "$shifted_run/np$np/run.exit"
  if (( statuses[0] != 0 || statuses[1] != 0 )); then exit 1; fi
done
for np in 2 4; do
  python3 ../tests/reference/openaccel/distributed_simple/compare_fields.py \
    "$shifted_run/np1/channel-fields.csv" "$shifted_run/np$np/channel-fields.csv"
done
```

In the working OpenAccel environment, reuse the recorded public executable;
no OpenAccel rebuild is needed. Keep `shifted_run` set to the printed directory
if changing shells. This uses only the pinned public channel, never a private deck.
Python needs NumPy/netCDF4 for the final comparison.

```bash
set -euo pipefail
: "${shifted_run:?Set shifted_run to the public MARS results directory}"
accel_root=/capstor/scratch/cscs/gandanie/git
reference="$shifted_run/openaccel"
srun --account=csstaff --time=00:10:00 --nodes=1 --ntasks-per-node=1 \
  --cpus-per-task=1 --cpu-bind=cores --export=ALL --kill-on-bad-exit=1 \
  env OMP_NUM_THREADS=1 OMP_PROC_BIND=close OMP_PLACES=cores \
  python3 "$accel_root/mars-v010-check/scripts/openaccel_simple_convergence.py" run \
  --capture "$accel_root/OpenAccel-reference-updates-IxgJIp/run" \
  --executable "$accel_root/OpenAccel-reference-updates-IxgJIp/source/build/openaccel-3D.exe" \
  --advection high-resolution --velocity-interpolation linear-linear --output "$reference"
python3 "$accel_root/mars-v010-check/scripts/openaccel_simple_convergence.py" compare \
  --reference "$reference" --mars "$shifted_run/np1" \
  --native-mesh "$accel_root/OpenAccel-simple-converged-20260925-151406/channel.exo" \
  --output "$shifted_run/openaccel-comparison.json"
```

Require the existing tolerances: 1e-6 for rank parity, 1e-5 for scaled reference
field differences and 1e-6 for final reference changes. The comparator checks the
new mode in the reference manifest, reference input and MARS log. Historical
unlabelled logs are accepted only for the old standard mode. A successful public
check does not establish compatibility with additional private expert settings.
