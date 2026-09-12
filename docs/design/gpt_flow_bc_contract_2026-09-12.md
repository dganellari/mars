# Velocity inlet, pressure outlet: source and units check

Author: GPT/Codex. Date: 2026-09-12.
Status: source inspection and user-provided physical-property clarification.
No confidential geometry, images or field outputs inspected; no production edit.
The user clarified that the PI's reference calculation uses OpenAccel. Treat
those observations as attributed reference results, not as measurements of MARS.

The specified physical target is an inward-normal velocity inlet, a pressure
outlet and no-slip walls elsewhere. No outlet resistance was specified. The
separately explored passive-resistance host model is not selected for production.
The clarification does not establish the precise pointwise/mean pressure law at
the outlet or the normalization of the reference solver's convergence residual.

## Fluid-property correction

The user supplied rho=1000 kg/m^3 and dynamic viscosity mu=1e-3 Pa s.
MARS's `--nu` is kinematic viscosity, as declared at `mars_pump.cu:133`:

```
nu = mu/rho = 1e-6 m^2/s
--rho=1000 --nu=1e-6
```

The earlier public-channel commands explicitly selected `--nu=1e-4` and thus
mu=0.1 Pa s at this density, 100 times the specified dynamic viscosity. They
remain valid performance and solver comparisons at their stated parameters;
they are not physical-property matches to the newly clarified target. Do not
rewrite their reported parameters or infer new GPU validation at nu=1e-6.
The user clarified that this higher viscosity was a workaround because the
water-viscosity case could not previously be run successfully; it was not an
unintentional units conversion. The next public GPU control changes only nu
to 1e-6, retaining beta=1, AMG-cycle and the previous acceptance tolerances.
The beta=.05 instability established on the public MARS replica does not show
that OpenAccel's own outlet implementation has the same instability.

At nu=1e-6, the executed host lagged-feedback spectrum still has eight growing
modes, radius 1.2340702145, with leading eigenpair residual 5.43e-15. The host
fixed-trace spectrum did not pass its existing eigenpair-residual gate; no
spectral stability claim is made for that new parameter set. This is an
eigensolve verification failure, not evidence by itself of a growing mode.
The separate 200-step fixed-trace host trajectory at nu=1e-6 stays bounded,
ending at speed 1.2327604 and continuity RMS 1.34e-16 /s. It omits advection;
the user-run GPU water test remains necessary.

## Inlet implementation

`inletPernodeNormal` defaults to true (`mars_pump.cu:122`). For a nondegenerate
inlet node, the driver constructs its area-weighted normal and assigns

```
a_i = sum_incident_inlet_faces A_f/3
u_target,i = -U(t)*a_i/|a_i|
```

The default minus sign assumes outward face winding. `--inlet-flip-normal`
reverses it. `--no-inlet-pernode-normal` uses one common averaged inlet
direction instead. Construction is at `mars_pump.cu:1043,1080–1106` and the
velocity target is applied at `mars_ns_pump_solver.hpp:5421–5445` when neither
pressure-driven inlet nor flux-Neumann inlet is selected.

For `--inlet-velocity=.5`, the nodal target magnitude is .5 m/s after the source
ramp. During a ramp of N steps it is `.5*min(1,step/N)`; the driver updates both
Uinf and the target at `mars_pump.cu:1767–1771`.

This is normal to the averaged nodal normal, not necessarily to every incident
triangle on a curved opening. Degenerate nodal vector sums fall back to the
global direction. Inlet setup still gates on the norm of the global summed
area vector, so a patch with canceling normals can skip per-node setup. These
are source limitations; no claim is made that either occurs on private data.

Opening nodes override wall membership at shared edges (`mars_pump.cu:814–823`).
Thus "all other sides no-slip" does not mean every node shared with a wall is
velocity-zero. This precedence requires an explicit boundary convention when
comparing codes.

Source inspection establishes the intended target, not the field actually
stored in a particular run. The geometry owner can inspect inlet vectors
locally after ramp completion: speed, inward orientation relative to outward
normals, and agreement with the prescribed target. That inspection remains
outside GPT's access to confidential artifacts.

## Convergence and Reynolds-number conventions

`--tol` controls the linear solve. The outlet path separately accepts full
stabilized continuity and boundary balance. `--steady-tol` currently tests
`abs(u_rms_new-u_rms_old)/u_rms_old` at reporting times for three consecutive
reports (`mars_pump.cu:1828–1832`), with a health check. It does not test the norm
of the velocity-field change or the steady momentum residual. Equal RMS speeds
can hide changing flow patterns. Setting it to 1e-6 alone does not reproduce
another solver's unspecified 1e-6 convergence criterion.

Local Reynolds number needs a matching velocity and length scale:
`Re = rho*U_local*D_local/mu = U_local*D_local/nu`. A bounding-box diagonal is
not a local passage diameter. No private diameter or local Reynolds number was
measured or inferred here.
