# Outlet correction: true-J FGMRES

**Current validation, GPT/Codex, 2026-09-11:** public 1/2/4-rank runs report PASS.
The signed-cut diagnostic fix and separate empty-opening fixture now need their
focused Daint gates. See [the current record and commands](gpt_outlet_cut_validation_2026-09-11.md).
The dated implementation/run notes below are historical; their pending statuses
refer to the evidence available at that time.

Author: GPT/Codex. Date: 2026-09-10.
Status: implementation and host validation; CUDA/MPI execution pending.

## Trigger and scope

The user's public-channel Daint log, job 4642390, exercises the corrected Armijo
rule. Corrections 5–7 pass while the old rule would reject them. Predicted and
measured residual ratios agree to roundoff along the tested directions. However,
RMS stalls at 0.021267878463389435 in the first physical step: damping falls from
0.0088508848 to 3.5493162e-10. Correction 8 cannot produce measurable decrease.
The line-search repair was necessary but insufficient; no full timestep passed.

The change replaces the approximate single-direction solve with restarted flexible
GMRES on the actual frozen correction operator. `Apre` becomes its right
preconditioner. The physical Chorin/BDF2 step, trace update, momentum relaxation,
and final continuity/boundary acceptance tolerances remain the same.

This is not SIMPLE. SIMPLE needs an outer momentum/pressure/flux iteration and a
separately derived momentum response. The flux evaluator, outlet quadrature, and
conservation gates can support it; the Chorin Jacobian must not be relabeled as its
Schur complement. The user approved this bounded correction-solver step, not an
indefinite expansion of the Chorin method.

## Discrete action and scaling

The PDE remains constant-density incompressible Navier–Stokes with physical
pressure, kinematic viscosity `nu`, and dynamic viscosity `rho*nu`. The public
duct's velocity inlet, no-slip walls, free outlet velocity, and frozen mean-pressure
trace are unchanged. Unknowns are nodal Tet4 velocity and pressure. Interior SCS
fluxes and the three area/3 samples per opening triangle define integrated volume
continuity `R`, in volume/time. Lumped nodal volume is `V_i`.

Using the boundary-aware production gradient `G_v` and prescribed-velocity mask Q:

```
h = dtEff/rho                       (dtEff=dt at BDF1, 2*dt/3 at BDF2)
delta_u = -h Q G_v phi
h J phi = R(delta_u, phi; delta_trace=0, delta_Gbar=0)
J phi = -R_current/h
```

The action evaluates the homogeneous flux directly, with zero trace and no smooth
reconstructed-gradient perturbation. It does not subtract two large residuals.
Interior and opening scatters are the same kernels as the physical residual. The
pressure boundary derivative is therefore included by evaluation, with no CSR
support assumption and no gauge removal. The action changes only scratch pressure
increments, their gradients, and scratch velocity/residual arrays.

Let `S_ii=1/sqrt(V_i)`. FGMRES solves `S J phi = -S R_current/h`, applying the right
preconditioner `Apre^-1 S^-1`. Its Euclidean residual norm equals the continuity
RMS times `sqrt(sum V)/h`, so its relative residual uses the physical acceptance
weights. This is left row scaling, not a change to the pressure unknown or PDE.

`mars_outlet_fgmres.hpp` holds the production controller. It stores each actual
preconditioned basis vector Z, uses two global Gram-Schmidt passes, applies Givens
rotations before handling breakdown, and recomputes the true residual before
accepting convergence. It respects an iteration cap shorter than a restart and
handles partial final restarts. Zero RHS is exact-zero tested; no absolute `1e-14`
cutoff silently accepts a small nonzero equation.

The existing CSR-only flexible GMRES implementation is not directly reused: it
lacks a matrix-free action callback. The new controller has device operations in
`mars_outlet_krylov.hpp` and an independent dense host adapter for regression tests.

## GPU and MPI behavior

Only owned DOFs enter the Krylov vectors and reductions. Each action publishes phi,
reverse-adds three gradient accumulators, publishes three masked velocity responses,
and reverse-adds continuity once. Frozen nodal coefficients and gradient samples
already have current halos. Unique opening-facet ownership is unchanged.

Restart depth is 40. Persistent workspace holds 85 owned-DOF vectors, three nodal
velocity responses, one nodal residual, and 41 projection scalars, in addition to
the existing solver scratch. Each orthogonalization pass batches its column dot
products into one device kernel and one small host MPI reduction. Full Krylov and
flow fields stay on device. No allocation of a full Jacobian or global gather is
introduced. Empty local projection work still contributes zeros to the reductions.

The existing Hypre helper solves `Apre` inside each preconditioning application and
may take different iteration counts, which is why the outer method is flexible.
That helper rebuilds Hypre setup each call; this implementation makes no performance
claim. Zero-owned-DOF support of the existing Hypre wrapper is not established by
the empty-projection gate. Empty-opening-facet ranks remain valid.

Outlet gradient launches now check launch errors without adding device-wide fences
to each action; halo operations and reductions provide the required dependencies.
The legacy gradient path retains its previous synchronization. The physical step
still ends with checked device synchronization.

## Gates and execution status

- Executed locally: 46 checks of the production FGMRES controller, including varying
  preconditioning, nonnormal stagnation, happy/singular breakdown, restart/cap
  handling, exact-zero and scaled RHS, nonfinite actions, and rejection of a false
  Hessenberg convergence claim. Clang C++20, `-Wall -Wextra -Werror`, ASan and UBSan.
- Executed locally: all 1,435 existing production outlet evaluator/algebra checks,
  with the same compiler warnings and sanitizers.
- Added, not GPU-executed: eight batched projection checks in
  `mars_outlet_kernel_gate`, including padded basis storage and ranks with zero
  local projection work. The kernel gate now has 393 checks per rank.
- Added, not GPU-executed: `--outlet-channel-check` tests `J*0`, homogeneity, unchanged
  physical/frozen/history fields, and finite differences along BOTH pressure and
  velocity responses for epsilon `1e-3,1e-4,1e-5`. Relative action error must be below
  `1e-8`. It runs each physical step, including BDF1 and BDF2, with the actual domain
  halos. Perturbation snapshots stay on device and are restored before the solve.

Next: build `mars_outlet_kernel_gate` and `mars_pump` on Daint. Run the kernel gate
on one rank, then the existing public one-rank channel command with
`MARS_SOLVE_TRACE=1`. Look for `[outlet-jacobian]`, `[outlet-krylov]`, and the final
eight-step channel PASS. Only after one-rank success proceed to two/four-rank gates
and comparisons. Failure is evidence to diagnose, not permission to relax tolerances.

The complete channel commands are in
[the integration instructions](gpt_outlet_channel_integration_2026-09-10.md).
