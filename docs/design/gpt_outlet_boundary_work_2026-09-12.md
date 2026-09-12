# Outlet pressure work and a passive mean-pressure candidate

Author: GPT/Codex. Date: 2026-09-12.
Status: derivation, source inspection, and executed public host checks.
No production CUDA or boundary defaults changed. The candidate below changes
the physical boundary law; it is not a beta=0.05 fix ready to deploy.
The user's subsequent BC clarification does not prescribe outlet resistance.
This candidate is therefore not selected for production; the reproducible host
reference below records the explored alternative, not an implementation plan.

## Decision and scope

Keep the GPU-validated beta=1 path with `MARS_OUTLET_PRECOND=amg-cycle` as the
working baseline. Its [recorded timing comparison](gpt_outlet_preconditioner_reuse_2026-09-12.md)
measures 7.24x lower runtime on the public channel with tracing/profiling off.

The old mean-pressure feedback can inject energy even with exact coupled
pressure/velocity solves. A passive alternative admits nonuniform pressure while
preserving its prescribed area mean, but ties the fluctuations to normal velocity
through a dimensional resistance. There is no justified conversion from the old
dimensionless beta to that resistance. Preserve the specified velocity inlet,
pressure outlet and no-slip walls; see the [BC and units check](gpt_flow_bc_contract_2026-09-12.md).

## Existing pressure feedback: an explicit energy witness

Use the existing public 425-node/1,536-tet procedural channel only. The PDE is
incompressible momentum with physical pressure p [Pa], density rho [kg/m^3],
kinematic viscosity nu [m^2/s], and prescribed inlet/wall velocity. Advection is
omitted in the host model. Viscosity uses the existing stiffness matrix and
natural treatment on free outlet velocity rows. Velocity and volume pressure
are nodal; the opening trace is a separate field. Opening quadrature is A_f/3
per vertex and interior velocity flux uses the existing SCS average.

Let M be lumped control-volume mass (volume, without density), B the integrated
velocity-divergence map restricted to free velocity columns, N the raw outlet
normal-flux map, and t the pressure trace. In scaled variables z=a*p/rho and
t_z=a*t/rho, a=dt for BDF1 and 2*dt/3 for BDF2, the gradient obeys

```
M g = -B^T z + N^T t_z.
u^T M g = -z^T B u + t_z^T N u.
```

This is the actual volume-gradient/trace split used by
`addOutletGradientTermKernel` in `mars_ns_pump_solver.hpp:11303`, called from
the predictor and pressure-increment gradient. The positive sign of the trace
term follows the outward normal and the momentum force `-grad(p)/rho`.
Positive pressure work in this formula removes kinetic energy; negative work
adds it. Prescribed velocities introduce separate source work, excluded in the
homogeneous eigenmode experiment.

The current `outletTraceValue` at `mars_ns_pump_solver.hpp:11265` sets
`t_z=(1-beta)*P*E*z`, where E restricts to outlet nodes and P removes the
scalar-area-weighted mean. No condition makes `t_z^T N u` nonnegative.

An exact solve of the existing fully implicit host block at beta=.05,
dt=2e-6, nu=1e-4 gives the following on its leading velocity-map eigenvector.
The history eigenvector has unit Euclidean norm; these are scaled algebraic
work values, not measured physical watts:

| Contribution | Value |
|---|---:|
| Volume pressure work `-z^T B u` | +0.0004644920975 |
| Outlet trace work `t_z^T N u` | -0.002516443069 |
| Total pressure work | -0.002051950971 |

The velocity-map eigenvalue is 1.451027454, yielding the previously verified
BDF2 amplification 1.639730525. The volume contribution dissipates energy on
this mode; the trace injects more. This is an explicit witness for this discrete
system, not a claim that all pressure-feedback outlets are unstable.
It supports examining the boundary equations before introducing SIMPLE: exact
solutions of these same block equations still grow.

## Passive alternative and its units

For every outlet node, accumulate scalar area `w_i=sum_f |A_f|/3` and vector
area `a_i=sum_f A_f/3` from uniquely owned facets. Define

```
v_n,i = (a_i dot u_i)/w_i
mean_A(v_n) = sum_i w_i*v_n,i / sum_i w_i
t_i = p_ref + R*(v_n,i - mean_A(v_n)),       R >= 0 [Pa s/m].
```

On a bent patch this is the quadrature average of normal velocity, not velocity
dotted with a normalized summed normal. Nodes with canceling normals retain
their scalar area. The formula preserves `mean_A(t)=p_ref` exactly.

For the fluctuating part of the trace, the boundary work is

```
sum_i (t_i-p_ref)*(a_i dot u_i)
    = R*sum_i w_i*(v_n,i-mean_A(v_n))^2 >= 0.
```

The mean contribution `p_ref*sum_i a_i dot u_i` is external pressure work and
is not claimed to be nonnegative. The resistance damps deviations from mean
normal flow; it does not impose a resistance on total discharge. This is a
nonlocal mean-pressure impedance law, different from keeping 95% of interior
pressure fluctuations. R=0 is uniform prescribed trace (the beta=1 reference).

Let W=diag(w), P=I-1*w^T/sum(w), and zeta=a*R/rho [m]. In scaled variables,
the homogeneous trace map is `t_z=F_u*u`, `F_u=zeta*P*W^-1*N`.
`P*W^-1=W^-1-1*1^T/sum(w)` is symmetric positive semidefinite. Therefore the
added integrated momentum operator `N^T F_u` is also symmetric positive
semidefinite. This proves passivity of this boundary term only. The existing
Rhie–Chow reconstruction has additional cross terms; this is not a proof of
energy stability of the complete discretization.

## Coupled equations and the derivative that must accompany the law

Retain the existing stabilized flux operators from the public replica: C is
the compact volume-pressure map, T the compact trace map, and W_g the
reconstructed-gradient flux map (all scaled as in the earlier temporal audit).
Set

```
S0 = C + W_g G_v
L  = T + W_g G_t
A_R = A + G_t F_u
B_R = B + L F_u

[ A_R   G_v ] [ u_new ] = [ r ]
[ B_R    S0 ] [ z_new ]   [ b ]
```

Here A is the BDF-scaled viscous momentum operator. Prescribed velocity columns
are eliminated into r and b, including their F_u terms. The host gate evaluates
the original flux expression separately from these block matrices. Its exact
Schur complement is `J_R=S0-B_R*A_R^-1*G_v`.

Both `G_t F_u` and `L F_u` are essential. Increments satisfy
`delta t_z=F_u*delta u`. Freezing this trace while changing velocity would solve
a different system. The reference keeps physical R constant across BDF startup;
zeta changes by 2/3 at the BDF1-to-BDF2 transition.

This candidate uses the exact viscous momentum response and fully current
pressure reconstruction. It is not the production Chorin correction with one
trace formula substituted. An iterative segregated solver could target this
block, but no SIMPLE implementation is supplied here. The existing FGMRES/AMG
infrastructure is reusable only after the new true action and its acceptance
residual have been implemented and checked.

## Executed host gates

`scripts/outlet_passive_trace_check.py` reuses the independent public spatial
assembly, not production CUDA kernels. Six cases passed 80 algebra checks in
total. The cases were chosen as zero-resistance controls, a moderate dimensional
example, and stiff-resistance stress tests, not as a sweep to tune a stable R.

| dt | nu | rho | R [Pa s/m] | BDF2 spectral radius | Growing modes |
|---:|---:|---:|---:|---:|---:|
| 2e-6 | 1e-4 | 1000 | 0 | 0.9999999934 | 0 |
| 2e-6 | 1e-4 | 1000 | 500 | 0.9999999916 | 0 |
| 2e-6 | 1e-4 | 1000 | 187500000 | 0.9999999916 | 0 |
| .01 | .1 | 1 | 0 | 0.9711307592 | 0 |
| .01 | .1 | 1 | 37.5 | 0.9602814500 | 0 |
| .01 | .1 | 1 | 712.5 | 0.9602098250 | 0 |

Checks cover Schur versus direct block solve, backward error, independently
evaluated continuity, trace mean and pressure work, mass-weighted adjointness,
impedance symmetry, constant-pressure/reference-shift invariance, finite-
difference Schur action, direct time-map action, and the leading eigenpair.
R=0 also reproduces the existing fully implicit fixed-trace reference.
Maximum non-FD check error across cases is 1.16e-14; finite-difference errors
are at most 8.12e-12. Acceptance limits are 1e-10 and 1e-7 respectively.
The spectral classification uses modulus greater than 1+1e-8; it does not bound
nonnormal transient growth or nonlinear advection.

Two 200-step forced histories at the small dt/nu stay bounded: R=500 ends with
speed 1.2327585, R=187500000 with speed 1.2329043. Their final integrated-
continuity RMS values are below 6.2e-16 /s. These are exact linear host solves,
not new GPU measurements or physical validation of either resistance.

Reproduce the main gates:

```bash
python3 scripts/outlet_passive_trace_check.py --resistance=0 --steps=0 --audit-pressure-feedback
python3 scripts/outlet_passive_trace_check.py --resistance=500
python3 scripts/outlet_passive_trace_check.py --resistance=187500000
python3 scripts/outlet_passive_trace_check.py --dt=.01 --nu=.1 --rho=1 --resistance=0 --steps=0
python3 scripts/outlet_passive_trace_check.py --dt=.01 --nu=.1 --rho=1 --resistance=37.5 --steps=0
python3 scripts/outlet_passive_trace_check.py --dt=.01 --nu=.1 --rho=1 --resistance=712.5 --steps=0
```

## Corner control and remaining limitations

An additional affine probe with p=2x-3y+5z and the exact outlet trace reproduces
the gradient on free interior nodes to 1.34e-14 in the existing spatial replica,
but fails at free wall/outlet corner nodes: at (4,0,1), the gradient is
(4,-8.25,6.25), rather than (2,-3,5). The driver explicitly lets opening nodes
override wall nodes (`mars_pump.cu:814`); the source uses the matching divergence-
transpose gradient (`mars_ns_pump_solver.hpp:8790` region). This is a boundary
accuracy concern in the inherited model, not a newly executed production-kernel
gate or proof of the beta instability's cause. The existing compact Tet4 affine
gate does not test this complete assembled nodal gradient.

The new host-only `--wall-precedence` control keeps those intersection nodes
velocity-fixed and rebuilds every dependent operator. At beta=.05 it still has
eight growing modes: lagged-map radius 1.258542323, fully implicit radius
1.661563096. Their time-map and leading eigenpair residuals are below 5.1e-15.
Changing corner precedence alone is therefore not a supported stability repair.
Both counterfactuals preserve the current default host model unless opted in:

```bash
python3 scripts/outlet_temporal_audit.py --beta=.05 --wall-precedence --steps=100 --spectrum
python3 scripts/outlet_coupled_time_check.py --beta=.05 --wall-precedence --steps=0 --expect=unstable
```

Before production use of the passive candidate: choose the physical boundary
target and resistance; settle/validate the corner boundary treatment; add real
kernel Jacobian and conservation gates; implement the complete coupled response;
then run public 1/4-rank and longer-time tests. Do not change beta's meaning or
promote a stable host spectrum into a production-accuracy claim.
