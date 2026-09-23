# Two native SIMPLE iterations on the public channel

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
A transition to a reversed outlet fails explicitly: selecting new artificial-wall
blocks is not integrated here. SIMPLEC, moving meshes, transient terms, forces,
turbulence, other boundary types, MPI ownership and general input decks remain
outside this executable. The two iterations do not establish nonlinear convergence.

Local strict C++ host builds pass 230 independent checks and 35,544 comparisons
against actual captured update states. The independent host oracle uses dense
partial-pivot elimination, not Hypre. Maximum scaled state difference is 1.43e-8
for iteration 2's raw pressure increment; corrected velocity differs by at most
5.89e-10 including the predictor. Captured linear solves used rtol=1e-8, so the
state gate allows `2e-11 + 1e-7*max(abs(expected))` for each field. True host solve
relative residuals are below 1e-13. Flux cancellation is checked independently;
the net balances (-0.0525 and approximately -0.00231336 kg/s) remain nonzero,
just as in the unconverged two-iteration reference.

The CUDA path is implemented but has not been compiled or executed locally.
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
integration gate. The following work is general driver integration, reversal
selection, distributed ownership/halos and converged public-case comparisons.
