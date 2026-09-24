# Converged public SIMPLE field comparison

The standalone CUDA/Hypre driver passed the public channel at iteration 1277
(user-provided Daint result, `simple-channel-4VjvlT`). Momentum was 9.90686e-7,
continuity 3.25978e-11, relative mass imbalance 1.98286e-12 and peak speed
0.188288 m/s. This is nonlinear convergence, not yet final OpenAccel field parity.

Reuse the existing instrumented OpenAccel executable and its completed public
capture. No STK compilation, MARS rebuild or repeated MARS run is needed.
`scripts/openaccel_simple_convergence.py run` checks the pinned deck, public mesh
and captured executable hashes before preparing a fresh output directory. It
turns off all capture hooks, so no per-element JSON is written during the long run.
The input retains uncorrected nodal output (`corrected_boundary_values: false`),
which is the state that MARS exports.

Only four controls change from the two-iteration reference deck:

| Control | Capture | Convergence run |
|---|---:|---:|
| Maximum outer iterations | 2 | 5000 |
| Native RMS residual target | 1e-6 | 1e-10 |
| Linear relative tolerance | 1e-8 | 1e-12 |
| Linear absolute tolerance | 1e-12 | 1e-14 |

Geometry, SIMPLE ordering, physical properties, upwind advection, pseudo-time,
relaxations, outlet beta and boundary conditions stay fixed. Each iteration is
saved to Exodus. The native residual target is deliberately tighter to reduce
reference stopping error; its normalization is not MARS's normalization.
An iteration-limit exit is rejected even if OpenAccel returns zero.

## Interactive Daint commands

Use the working **OpenAccel environment**. From `mars/mlir`, pull the MARS scripts;
no CMake/build command is needed:

```bash
git pull --ff-only
reference_run=/capstor/scratch/cscs/gandanie/git/OpenAccel-simple-converged-$(date +%Y%m%d-%H%M%S)
printf 'Reference results: %s\n' "$reference_run"
srun --account=csstaff --time=00:30:00 --nodes=1 --ntasks-per-node=1 \
  --cpus-per-task=1 --cpu-bind=cores --export=ALL --kill-on-bad-exit=1 \
  env OMP_NUM_THREADS=1 OMP_PROC_BIND=close OMP_PLACES=cores \
  python3 ../scripts/openaccel_simple_convergence.py run \
  --capture /capstor/scratch/cscs/gandanie/git/OpenAccel-reference-updates-IxgJIp/run \
  --executable /capstor/scratch/cscs/gandanie/git/OpenAccel-reference-updates-IxgJIp/source/build/openaccel-3D.exe \
  --output "$reference_run"
```

Then compare on the login node; no GPU allocation or MARS environment switch:

```bash
comparison_python=/capstor/scratch/cscs/gandanie/git/OpenAccel/spack/var/spack/environments/accel-clean-20260914/.spack-env/view/bin/python3
"$comparison_python" ../scripts/openaccel_simple_convergence.py compare \
  --reference "$reference_run" \
  --mars /capstor/scratch/cscs/gandanie/git/mars/mlir/simple-channel-4VjvlT \
  --output "$reference_run/field-comparison.json"
```

That is the Spack-view Python previously used for Exodus statistics. `compare`
requires numpy and netCDF4; `run` uses only the standard library. Keep the printed
reference directory if changing shells. Both commands preserve existing outputs.
The runner checks the binary hash against the old capture; if it changed, stop
and inspect provenance rather than removing the check or recompiling blindly.

## Comparison contract

The equations and physical units are those in [README.md](README.md): steady
laminar incompressible Navier–Stokes, physical pressure in Pa, density 1 kg/m³,
U=0.1 m/s. Both fields use the same absolute pressure boundary condition.

- Require completed run metadata, unmodified result/deck/mesh/log hashes and the
  pinned OpenAccel revision. Require MARS's convergence banner to match the last
  metrics row and all original convergence gates to pass.
- Match all 425 nodes using the Exodus `node_num_map` and the MARS `channel.json`
  global-ID mapping, then check coordinates to 1e-12 m. Never compare by file order.
- Read the last two reference states, requiring their saved iteration values to
  be the reported final iteration and its predecessor. Reject incomplete or
  nonfinite fields and unsupported output layouts.
- Report maximum and nodal RMS errors of velocity scaled by U and pressure scaled
  by rho*U²=0.01 Pa. These comparison RMS norms are unweighted nodal norms, not
  either solver's residual norm. Do not subtract a pressure mean; report the mean
  difference separately so a pressure-level defect remains visible.
- Require maximum vector velocity error and absolute pressure error <=1e-5 in
  these scales (1e-6 m/s and 1e-7 Pa). Require the reference's last-step maximum
  field changes <=1e-6 in the same scales. These are fixed engineering comparison
  tolerances, not bitwise tests or a proof that stopping error is bounded by the
  last change. Tight linear/reference solves reduce contamination by stopping
  error; a failed gate needs investigation, not automatic tolerance relaxation.

The JSON records metrics, tolerances and evidence hashes, including on a numerical
comparison failure. Raw velocity/pressure parity does not replace stabilized-flux
conservation checks, analytic spatial-accuracy tests, MPI validation or pump work.
The already executed two-iteration gate separately checks stored flux and continuity.

Local tests exercise real NetCDF readers on explicitly synthetic fields (both
combined and split coordinate/nodal storage), permuted IDs, pressure offsets,
state drift, stale output, duplicate IDs, changed files and failed launches. They
are transport/comparison tests, not an executed OpenAccel convergence result:

```bash
python3 -m unittest discover -s scripts -p test_openaccel_simple_convergence.py -v
```
