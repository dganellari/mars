# High-resolution SIMPLE advection

`--advection high-resolution` selects the ordinary velocity limiter from
OpenAccel `0d69041ba1afda63e9e4328d9e0d9834bba37756`. Upwind remains the default.
This is steady incompressible laminar flow: density in kg/m³, dynamic viscosity
in Pa s, velocity in m/s and physical pressure in Pa. It retains the existing
Tet4 median-dual quadrature, boundary equations and SIMPLE pseudo-time terms.
The inlet has prescribed inward-normal velocity, walls are no-slip, and the
outlet uses the existing mean-pressure trace and close/reopen treatment.

## Discrete contract

For a component u at a node, take the minimum and maximum over that node and
its Tet4 edge neighbours. For every adjacent interior SCS point and every
physical boundary sample assigned to the node, compute

```
r = grad(u) dot (x_sample - x_node)
r_hat = r + epsilon if r > 0, otherwise r - epsilon
y = min(u_max - u, u - u_min) / abs(r_hat)
candidate = min(1, (y*y + 2*y)/(y*y + y + 2))
```

Use double machine epsilon; samples with `abs(r_hat) <= 100*epsilon` do not
limit. Standard coordinate shapes apply: interior weights 13/36 at the two
edge nodes and 5/36 elsewhere; boundary weights 11/18 at the nearest node and
7/36 elsewhere. All boundaries constrain the limiter, including closed outlets.
Boundary prescribed values do not enlarge the nodal extrema stencil.

Take the minimum candidate over all samples, then update persistent history:

```
beta = min(candidate, 0.75*beta_previous + 0.25*candidate)
```

Beta starts at zero, updates once at initialization and once per corrected
iterate, and remains fixed during repeated assembly of that iterate. The
limiter is componentwise and two-sided; it is not a sign-selected slope ratio.
Affine fields need not retain beta=1 near stencil extrema.

For a stored mass flux m [kg/s], the interior momentum flux is
`m * (u_upwind + beta_upwind * grad(u_upwind) dot dx)` [N]. The two nodal RHS
contributions are equal and opposite. The implicit upwind matrix is unchanged;
the reconstruction is a deferred RHS correction. Diffusion and pressure
correction are unchanged. The advection coefficient beta here is distinct from
the outlet `--outlet-beta` control.

Source correspondence: OpenAccel `nodeField.hpp` lines 1957–2128 (extrema),
2277–2607 (interior limiter), 3182–3466 (boundary limiter and relaxation),
`navierStokesEquation.cpp` lines 97–100 (initialization), and
`segregatedFlowEquations.cpp` lines 213–225 (corrected-velocity update).
MARS uses `mars_segregated_high_resolution.hpp` for the limiter and the existing
`tet_interior` reconstruction through `SimpleInterior<3>`.

Scope matches the public reference deck: cap 1, unrelaxed/unlimited gradients,
constant density, no shock sensor, no interfaces. Other OpenAccel high-resolution
variants and turbulent models are not implied by this option.

## GPU and MPI

The existing device CSR supplies each owned node's complete neighbour stencil.
Bounds, sample minima and history updates run on the GPU with persistent scratch.
Only owned rows are evaluated; neighbour velocities are already halo-complete.
Complete element stars provide all required samples. The existing first halo round publishes pressure gradient,
velocity gradient and beta together (15 doubles). The full iteration still has
four field rounds, now 26 doubles per ghost. No mesh/field host staging or reverse
limiter reduction is added. File I/O, MPI/Hypre API calls and scalar reports remain
host-controlled.

## Validation and interactive Daint run

Local validation on 2026-09-28: 3772 host algebra checks cover two-sided bounds, extrema, small slopes, upward
relaxation, boundary constraints and conservative reconstruction. Distributed
tests compare matrix rows, fields and limiter histories on 1/2/4 ranks, including
an oblique sheared initial field, outlet reversal and poisoned ghost gradients.
All 62 existing CPU/MPI regression cases passed. The first high-resolution MPI
checks exposed an undersized halo stride capacity; after increasing it to 15,
the focused 10-case CPU/MPI suite passed under a strict C++20 build. The 81-node
host channel converged at iteration 2907, and 25 reference-comparator tests passed.
Subsequent user-reported Daint runs of the 425-node public channel converged at
1318 iterations on 1/2/4 GPU ranks. OpenAccel converged at 1315; the saved-field
comparison passed with maximum scaled errors of 3.1101e-9 in velocity and
1.2877e-8 in absolute pressure, below 1e-5. Results are in
`/capstor/scratch/cscs/gandanie/simple-highres-yQnE5j`; the reference is
`/capstor/scratch/cscs/gandanie/openaccel-highres-yQnE5j/run`.
The later optimized channel regression also passed on 1/2/4 ranks. These are
user-reported GPU results, separate from the host checks above and from the
[independently rechecked upwind duct study](../simple_duct/DAINT_RESULTS.md).
High-resolution duct refinement remains unvalidated. The commands below
reproduce the public-channel checks.

From the configured MARS CUDA/Hypre build, after pulling:

```bash
cmake -S .. -B .
cmake --build . --parallel 4 --target mars_segregated_simple mars_simple_high_resolution_check
./mars_simple_high_resolution_check
highres_run=$(mktemp -d /capstor/scratch/cscs/gandanie/simple-highres-XXXXXX)
channel_mesh=/capstor/scratch/cscs/gandanie/git/OpenAccel-simple-converged-20260925-151406/channel.exo
git rev-parse HEAD > "$highres_run/mars-revision.txt"
set -o pipefail
for np in 1 2 4; do
  mkdir "$highres_run/np$np"
  srun --account=csstaff --time=00:20:00 --nodes=1 --ntasks-per-node="$np" \
    --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
    ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_segregated_simple \
    --mesh "$channel_mesh" --output-prefix "$highres_run/np$np/channel" \
    --advection high-resolution --iterations 5000 --report-every 100 \
    --residual-tol 1e-6 --mass-tol 1e-6 --change-tol 1e-6 \
    2>&1 | tee "$highres_run/np$np/run.log" || break
done
python3 ../tests/reference/openaccel/distributed_simple/compare_fields.py \
  "$highres_run/np1/channel-fields.csv" "$highres_run/np2/channel-fields.csv"
python3 ../tests/reference/openaccel/distributed_simple/compare_fields.py \
  "$highres_run/np1/channel-fields.csv" "$highres_run/np4/channel-fields.csv"
printf 'Saved: %s\n' "$highres_run"
```

Run the independent reference from the working OpenAccel environment. Reuse the
existing executable; no STK rebuild or fresh capture is needed. The script verifies
its checksum against the earlier public capture, changes only the declared public
advection/convergence controls, and disables export instrumentation.

```bash
highres_reference=$(mktemp -d /capstor/scratch/cscs/gandanie/openaccel-highres-XXXXXX)/run
srun --account=csstaff --time=00:10:00 --nodes=1 --ntasks-per-node=1 \
  --cpus-per-task=1 --cpu-bind=cores --export=ALL --kill-on-bad-exit=1 \
  env OMP_NUM_THREADS=1 OMP_PROC_BIND=close OMP_PLACES=cores \
  python3 /capstor/scratch/cscs/gandanie/git/mars-v010-check/scripts/openaccel_simple_convergence.py run \
  --capture /capstor/scratch/cscs/gandanie/git/OpenAccel-reference-updates-IxgJIp/run \
  --executable /capstor/scratch/cscs/gandanie/git/OpenAccel-reference-updates-IxgJIp/source/build/openaccel-3D.exe \
  --advection high-resolution --output "$highres_reference"
```

With NumPy/netCDF4 available, compare the saved one-rank fields (use the paths
printed above if changing shells):

The current distributed driver does not print the mesh path. The comparator
accepts its exact single-rank scheme banner and checks the pinned mesh hash,
source-row/global-ID mapping and output coordinates. The report records
`native_input_path_recorded: false`; this does not attest which file path was
opened. Legacy logs that include a path must still match it exactly. Existing
saved fields can be compared without rebuilding or rerunning either solver.

```bash
python3 /capstor/scratch/cscs/gandanie/git/mars-v010-check/scripts/openaccel_simple_convergence.py compare \
  --reference "$highres_reference" --mars "$highres_run/np1" \
  --native-mesh "$channel_mesh" --output "$highres_run/openaccel-comparison.json"
```

Require convergence in both solvers and the existing 1e-5 scaled field tolerance,
including absolute pressure, plus the reference final-change tolerance of 1e-6.
Do not compare high-resolution results against the old upwind fields. The Python
reader is solely for reference validation; MARS continues to read Exodus natively.
