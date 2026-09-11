# Fully implicit coupling does not repair the blended outlet

Author: GPT/Codex. Date: 2026-09-11.
Status: derivation and executed public host checks; no production implementation.
Source baseline: `f82222ba`, including the independently assembled operators in
`scripts/outlet_temporal_audit.py`.

Updating the trace, reconstructed gradient, and momentum response consistently
at the new pressure does not stabilize the beta=0.05 boundary formulation in
this model. The resulting BDF2 amplification factor is 1.63973. The beta=1
counterpart has no growing eigenmodes in the two parameter sets tested.

This narrows the earlier [temporal audit](gpt_outlet_temporal_audit_2026-09-11.md):
lagged feedback is implicated in the observed GPU failure, but removing its lag
does not suffice. The spatial boundary closure also needs examination. Do not
implement the fully implicit beta=0.05 candidate as a fix.

## Equations and conventions

The public procedural channel, Tet4 nodal spaces, volume and opening quadrature,
pressure-gradient signs, prescribed velocity mask, and nodal stabilization
coefficient are inherited from the previous audit. The continuum equations are
incompressible momentum with physical pressure p, density rho, kinematic
viscosity nu, and zero body force. Advection is omitted to isolate the linear
coupling; this is not a nonlinear CFD validation.

Let a=dt_eff and z=a*p/rho. B maps velocity to integrated continuity flux.
The pressure trace has p_ref=0 and map F=(1-beta)*(E-1*w^T*E), with scalar-area
weights w and outlet restriction E. Write G=G_v+G_t F for the complete pressure
gradient. C and T are the compact volume-pressure and trace-flux maps divided
by a/rho; W maps a reconstructed scaled gradient to flux. With every pressure
quantity evaluated at z_new, the pressure part of continuity is

```
S = C + T F + W G.
```

All these maps include the original masks and boundary quadrature. In particular,
the interior reconstructed-gradient mask and the boundary face/opposite blend
are preserved. The host coefficient is D/(a/rho), with relax_u=0.3 folded in
once. It is frozen within each step, using that step's BDF momentum diagonal.

On free velocity DOFs define

```
A_hat = M^-1 * a * (M/a + nu*K),   H = A_hat^-1.
```

G below is restricted to these DOFs, and B_f is B restricted to their columns.
After eliminating prescribed velocities, the proposed fully implicit step is

```
[ A_hat   G ] [ u_new ] = [ r ]
[ B_f     S ] [ z_new ]   [ b ]
```

Here r contains w=u_old for BDF1 or w=(4*u_old-u_older)/3 for BDF2, plus the
diffusion Dirichlet lift. The continuity source b=-B_fixed*u_prescribed includes
the inlet flux. The pressure level is set by the outlet closure, with no added
pin or mean subtraction. a=dt for startup and 2*dt/3 thereafter.

The exact Schur elimination is

```
J = S - B_f H G
z_new = J^-1 (b - B_f H r)
u_new = H r - H G z_new.
```

This replaces the mass-only pressure correction response with H G and refreshes
both trace and reconstructed gradient. Pressure has no independent lagged state
in this counterfactual. It is a different method from production, not another
execution of its existing correction loop.

For homogeneous boundary data, the velocity history map is

```
V = H + H G J^-1 B_f H
u_new = V (4*u_old-u_older)/3.
```

An eigenvalue mu of V gives the two BDF2 roots satisfying
`lambda^2 - (4*mu/3)*lambda + mu/3 = 0`. The script checks the leading root
against the full two-history-state action, not just the polynomial formula.

## Executed checks and results

The reproducible host check is
[`scripts/outlet_coupled_time_check.py`](../../scripts/outlet_coupled_time_check.py).
It reuses the public host spatial assembly; it does not execute CUDA/HIP or MPI.

| Fully implicit host experiment | BDF2 spectral radius | Growing modes |
|---|---:|---:|
| beta=0.05, dt=2e-6, nu=1e-4 | 1.6397305245 | 8 |
| beta=1, dt=2e-6, nu=1e-4 | 0.9999999934 | 0 |
| beta=0.05, dt=0.01, nu=0.1 | 1.5279805481 | 8 |
| beta=1, dt=0.01, nu=0.1 | 0.9711307592 | 0 |
| beta=0.05, small dt/nu, reconstructed-gradient flux removed | 1.5038528634 | 8 |

The last row is an isolation experiment, not a proposed stabilization. Growth
persists after removing W G, so refreshing or removing that term alone does not
resolve the problem. No parameter sweep was used to choose a damping value.

Each of the five cases passed seven algebra checks:

- Schur elimination versus a separately solved full block system: relative
  difference <=2.46e-15.
- Block backward error <=5.59e-18; momentum and continuity residual checks
  <=1.09e-15.
- Centered finite difference of the flux after applying the momentum response
  versus J: relative difference <=3.45e-12.
- Homogeneous time-map action versus direct stepping: relative difference
  <=1.51e-15.
- Leading BDF2 eigenpair residual <=4.68e-15.

The checks use arbitrary nonzero prescribed velocities for the block comparison,
then homogeneous data for the spectrum. Thresholds are 1e-10 for algebra/eigenpair
checks and 1e-7 for the finite difference, comfortably above observed roundoff.
The classification uses |lambda|>1+1e-8; it is not a bound on transient growth.

The small-dt leading unstable velocity mode has 86.4% of its lumped-mass-weighted
energy on outlet nodes. A separate forced 200-step host trajectory, with the
public inlet ramp, reaches speed 1.326e40 for beta=0.05. The beta=1 counterpart
ends at 1.23276. These are host results for the counterfactual method. The previous
production beta=1 GPU run ended at 1.233, but that numerical proximity is not a
device parity check for the new implicit model.

Reproduce with NumPy from the repository root:

```bash
python3 scripts/outlet_coupled_time_check.py --beta=.05 --expect=unstable
python3 scripts/outlet_coupled_time_check.py --beta=1 --expect=stable
python3 scripts/outlet_coupled_time_check.py --beta=.05 --dt=.01 --nu=.1 --steps=0 --expect=unstable
python3 scripts/outlet_coupled_time_check.py --beta=1 --dt=.01 --nu=.1 --steps=0 --expect=stable
python3 scripts/outlet_coupled_time_check.py --beta=.05 --compact-only --steps=0 --expect=unstable
```

## Architectural consequence

On the free velocity rows, the current maps obey
`G_v=-M^-1*B_f^T`. With feedback, `G=G_v+G_t F`: the trace adds boundary pressure
work. In the exact Schur complement it contributes `-B_f H G_t F`, alongside
changes in the pressure-flux map S. Their combined stability is not guaranteed
by satisfying continuity or by selecting an implicit solver. The measured growing
mode is a concrete failure of this combination on the public fixture; this is
not a general theorem about all blended pressure outlets.

For a frozen linear system, converged SIMPLE iterations would solve the block
equations they target. If they target the block above at beta=0.05, the growing
physical-step mode remains even with exact solves. Thus introducing SIMPLE alone
cannot repair this particular target system. SIMPLE is still not implemented
in production MARS.

Keep the executed beta=1 GPU case as the public development baseline. The exact
block above at beta=1 provides a host reference for a future segregated iterative
solver. Before restoring an average-pressure feedback boundary, derive its
momentum boundary work and continuity flux together and require a stable time-map
gate. A guessed trace damping or another update-order change is not supported
by these results. No production default, boundary kernel, or solver was changed
in this review.
