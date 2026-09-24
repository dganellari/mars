# Two native SIMPLE iterations on the public channel

The standalone single-GPU run now has Daint convergence evidence at iteration
1277. Continue with the [converged OpenAccel field comparison](convergence.md),
which reuses the existing reference executable and the completed MARS run.

This gate starts from zero velocity, pressure, stored flux, mass divergence and
wall coefficients. MARS computes two complete steady SIMPLE iterations using
native geometry, reconstruction, assembly, linear solves and state updates.
The captured OpenAccel fields are read only as expected answers. This extends
the preceding frozen assembly gate; it is still a restricted single-rank
integration executable, not a general mesh-based solver or convergence test.

The public fixture is the same 425-node, 1,536-Tet4 channel. Its pinned deck uses
rho=1 kg/m^3, dynamic viscosity=0.1 Pa s, inlet normal speed=0.1 m/s, no-slip
walls and one average-pressure outlet with beta=0.05 and reference pressure=0.
No pump inputs are used. The packer verifies the existing capture, deck and mesh
hashes. OpenAccel stays at 0d69041 with the pinned e351ba5 solver submodule.

## Equations and ordering

The target equations are steady incompressible Navier–Stokes,
`div(rho*u)=0` and `rho*(u.grad)u - div(mu*(grad(u)+grad(u)^T)) + grad(p)=0`.
Pressure is physical pressure in Pa. The discretization is the previously
verified median-dual CVFEM with upwind advection and Rhie–Chow mass flux.
The pseudo-time of 0.01 s contributes rho*V/pseudo_dt to momentum's increment
matrix; it does not advance physical time or BDF history.

1. Reconstruct current velocity and pressure gradients and assemble the full
   three-component momentum system `A_u*delta_u=b_u`. Apply alpha_u=0.3 only
   through its existing diagonal relaxation; multiply boundary-node RHS by
   0.75 once per node. Compute `d=V/(relaxed diagonal+SMALL)` and add delta_u.
2. Refresh the outlet trace using the area-weighted mean of old pressure.
   Assemble `A_p*phi=-net_mass_flux` using the computed d and predicted velocity.
   The outlet partial supplies the pressure-level anchor; no pin or mean removal.
3. Update `p += 0.3*phi` and reconstruct the gradient of **raw phi**.
4. Update interior and boundary mass flux using new p, predicted u and the old
   reconstructed grad(p). Apply relax_mass=0.75 to stored history. Apply the
   outlet flux clipping/reversal decision, then accumulate nodal mass divergence.
5. Correct `u -= d*grad(phi)`. Neither pressure relaxation nor another velocity
   relaxation multiplies this correction. The next iteration reconstructs the
   updated fields and uses the evolved flux/mass-divergence history.

The reference fills laminar wall coefficients in postSolve, so iteration 1 uses
zero coefficients. Later iterations use `4*mu*|A_sample|^2/V_tet`, equivalent to
its quarter-normal-edge formula. Changing this startup behavior would break
reference parity even if it were a desirable separate physics change.

The pressure matrix and momentum matrix are not assumed symmetric. CUDA uses
the existing MARS Hypre GMRES/BoomerAMG adapter with node-major component DOFs,
zero initial increments and rtol=1e-12. Momentum requests three-function AMG.
Block CSR expands on-device into sorted scalar CSR without transposing blocks.
Each solve is checked using a separately evaluated true `A*x-b` residual:
`||r||_2 <= 1e-13 + 1e-10*||b||_2`. No Hypre setup is reused across changed matrices.

## Bounds and evidence

Only the public fixed-frame, constant-property, upwind SIMPLE profile is enabled.
The integrated reversal path is described below. SIMPLEC, moving meshes, transient terms, forces,
turbulence, other boundary types, MPI ownership and general input decks remain
outside this executable. The two iterations do not establish nonlinear convergence.

Local strict C++ host builds pass 457 independent checks and 35,544 comparisons
against actual captured update states. The independent host oracle uses dense
partial-pivot elimination, not Hypre. Maximum scaled state difference is 1.43e-8
for iteration 2's raw pressure increment; corrected velocity differs by at most
5.89e-10 including the predictor. Captured linear solves used rtol=1e-8, so the
state gate allows `2e-11 + 1e-7*max(abs(expected))` for each field. True host solve
relative residuals are below 1e-13. Flux cancellation is checked independently;
the net balances (-0.0525 and approximately -0.00231336 kg/s) remain nonzero,
just as in the unconverged two-iteration reference.

The original two-iteration CUDA/Hypre path passed on Daint (35,544 comparisons).
The shared runtime, new reversal gate and convergence driver require new GPU validation.
Full array downloads exist only for these validation comparisons. Geometry,
assembly, gradients, scalar CSR, solution updates and history stay on-device
between comparisons; solver acceptance downloads scalar reductions only.
Persistent arrays and graph storage are reused. No performance claim is made.

## Interactive Daint gate

After these changes have been reviewed and published, from the configured MARS
CUDA/Hypre build directory (`mars/mlir`), run:

```bash
(
set -euo pipefail
git pull --ff-only
simple_work=$(mktemp -d "$PWD/segregated-simple-XXXXXX")
python3 ../scripts/prepare_openaccel_simple.py \
  /capstor/scratch/cscs/gandanie/git/OpenAccel-reference-updates-IxgJIp/run/exports \
  --boundary /capstor/scratch/cscs/gandanie/git/OpenAccel-reference-boundary-TljPP1/run/exports/boundary \
  --updates /capstor/scratch/cscs/gandanie/git/OpenAccel-reference-updates-IxgJIp/run/exports/updates \
  --output "$simple_work/inputs.txt"
cmake -S .. -B .
cmake --build . --target mars_segregated_simple_check mars_segregated_simple_algebra_check -j4
./examples/distributed/unstructured/mars_segregated_simple_algebra_check
srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=1 \
  --export=ALL --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_segregated_simple_check \
  "$simple_work/inputs.txt" 2>&1 | tee "$simple_work/simple.log"
git rev-parse HEAD > "$simple_work/mars-revision.txt"
sha256sum ./examples/distributed/unstructured/mars_segregated_simple_check \
  "$simple_work/inputs.txt" > "$simple_work/sha256.txt"
printf 'Saved: %s\n' "$simple_work"
)
```

No OpenAccel rebuild or rerun is needed. A pass closes the single-GPU two-iteration
integration gate. The standalone public driver below extends that gate; distributed ownership/halos,
general meshes and converged reference-field comparisons remain pending.

## Standalone public-channel iteration driver

`mars_segregated_simple` reads coordinates, connectivity and boundary tags only.
It uses the same `SimpleRunner` and CUDA/Hypre solves as the two-iteration gate;
no captured stage fields or expected answers are loaded. Material, relaxation,
pseudo-time, upwind interpolation and boundary settings remain the pinned profile
above. This is a one-rank public-channel driver, not an arbitrary-mesh pump driver.

Outlet reversal flags now persist between iterations. The pinned reference's
artificial slip wall omits the opening's momentum/pressure boundary blocks and
sets its mass flux to zero. Closed faces are excluded from the pressure mean and
retain their old trace. Reopening requires outward face-average velocity and
interior mean pressure at least as large as that trace. The current iteration's
flux stays zero when reopening; the next iteration restores the opening blocks.
All-closed outlets fail explicitly before a pressure solve: the pressure anchor
would be absent with prescribed inflow. No pressure pin or fallback wall model
is silently introduced. Source: pinned `flowModel.cpp` outlet flag/trace updates
and segregated momentum/pressure boundary assemblers.

At the start of each iteration, the driver assembles the momentum defect for the
current velocity, pressure and stored flux history. Let `R_u` be that nodal RHS
before the boundary factor 0.75, `R_m` the stored nodal mass-flux sum and `V` the
nodal dual volume. It reports

- momentum: `sqrt(sum(|R_u|^2/V)/sum(V))*L/(rho*U^2)`;
- continuity: `sqrt(sum(R_m^2/V)/sum(V))*L/(rho*U)`;
- mass balance: `abs(Q_in+Q_out)/(rho*U*A_in)`;
- velocity/pressure changes: volume-weighted RMS changes divided by `U` and
  `rho*U^2`, respectively;
- flux-history change: the largest sample-flux change divided by `rho*U*A_in`.

Here `L=1 m` is the channel height, `U=0.1 m/s`, and `Q` denotes mass flow in
kg/s. These are dimensionless MARS norms, not OpenAccel's displayed RMS columns.
Pseudo-time and diagonal relaxation affect the increment matrix, not this steady
momentum RHS. The boundary RHS factor is divided out in the diagnostic. Continuity
includes every node and uses the same stored stabilized flux as the iteration.
`sum(R_m)=Q_in+Q_out` is checked to `1e-10` relative to prescribed inflow.

Success requires momentum and continuity below `--residual-tol`, mass balance
below `--mass-tol`, all three state-change measures below `--change-tol`, no flag
changes and at least two completed iterations. Reaching the iteration limit prints
`NOT CONVERGED` and returns 2; a runtime/linear/consistency failure returns 1.
Final nodal fields and per-iteration metrics are CSV files. Field `node` is the
zero-based packed row; the mesh JSON records `node_global_ids` for reference matching.
Full field downloads
occur only for the explicit final export; scalar reductions are used during
iteration. The host/direct build is a validation oracle, not a production CPU path.

The host/direct standalone run converged at iteration 1,277 with the default
`1e-6` thresholds: momentum `9.90686e-7`, continuity `3.25978e-11`, relative mass
imbalance `1.98244e-12`, velocity change `2.56626e-10`, pressure change `1.14102e-9`
and flux-history change `4.44996e-12`. Peak speed was `0.188288 m/s`; no face reversed
in that run. This is internal convergence on the public mesh, not a converged
OpenAccel field comparison or GPU result. The new synthetic reversal gate supplies
the missing close/reopen branch coverage.

The extended `mars_segregated_simple_check` additionally exercises native kernels
on a unit tetrahedron through closure, retained trace, pressure-qualified reopening,
restored flux/pressure anchor and device diagnostics reductions. Its two-iteration
reference comparisons remain unchanged. CPU algebra and host checks do not certify
the new CUDA branches; run the following commands on Daint.

From the existing MARS CUDA/Hypre build directory, after pulling this change:

```bash
(
set -euo pipefail
git pull --ff-only
cmake -S .. -B .
cmake --build . --parallel 4
simple_run=$(mktemp -d "$PWD/simple-channel-XXXXXX")
python3 ../scripts/prepare_openaccel_simple.py \
  /capstor/scratch/cscs/gandanie/git/OpenAccel-reference-updates-IxgJIp/run/exports \
  --boundary /capstor/scratch/cscs/gandanie/git/OpenAccel-reference-boundary-TljPP1/run/exports/boundary \
  --mesh-only --output "$simple_run/channel.txt"
./examples/distributed/unstructured/mars_segregated_simple_algebra_check
srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=1 \
  --export=ALL --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_segregated_simple_check \
  "$PWD/segregated-simple-ZH466h/inputs.txt" 2>&1 | tee "$simple_run/gate.log"
git rev-parse HEAD > "$simple_run/mars-revision.txt"
sha256sum ./examples/distributed/unstructured/mars_segregated_simple \
  ./examples/distributed/unstructured/mars_segregated_simple_check \
  "$simple_run/channel.txt" > "$simple_run/sha256.txt"
# Save the path before launching so a failed or unconverged run is easy to locate.
printf 'Results: %s\n' "$simple_run"
srun --account=csstaff --time=00:20:00 --nodes=1 --ntasks-per-node=1 \
  --export=ALL --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_segregated_simple \
  --mesh "$simple_run/channel.txt" --output-prefix "$simple_run/channel" \
  --iterations 2000 --report-every 20 \
  --residual-tol 1e-6 --mass-tol 1e-6 --change-tol 1e-6 \
  2>&1 | tee "$simple_run/run.log"
)
```

No OpenAccel build or new capture is needed for this gate/run. A later converged
field comparison needs a longer OpenAccel run with this same public deck; two
captured iterations cannot supply that final reference. MPI and general input
meshes remain outside this milestone.
