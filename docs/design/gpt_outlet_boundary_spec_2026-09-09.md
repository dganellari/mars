# Average-pressure outlet: discrete contract and acceptance gates

Author: GPT/Codex. Date: 2026-09-09.
Status: mathematical specification and independent review; no CUDA implementation.
Source snapshot: MARS `018807b4f58bed851e33c72c9e1c808dc56c69aa`.

This develops the outlet proposal in
[GPT's architecture review](openaccel_math_review_2026-09-09.md), following
[Claude's implementation status](outlet_trace_status.md). Claude authored the
committed trace scaffold, collective fix, and existing host trace test. The equations,
additional gates, and numerical counterexample below are GPT's work. Sources were
read from an immutable copy of the named commit while Claude worked in the live tree.
Later changes are outside this snapshot review. No confidential geometry or results
were inspected.

## Decision and first supported path

Claude correctly stopped before freeing pressure rows without their boundary operator.
Complete the mask, continuity residual, pressure derivative, pressure gradient/corrector,
and gauge together. The gradient/corrector is an essential part of this change.

Use the existing **SCS → Hypre → Apre** entry for the first implementation, with a
declared approximate correction operator and a full-residual convergence test.
Keep Chorin/BDF2; this does not require a SIMPLE momentum outer loop.

The explicit target is `--solver=hypre --vms-stab --rc-implicit`, with
`--outlet=do-nothing --outlet-beta=0.05`. The driver otherwise defaults to a prescribed
mass-conserving velocity outlet. For the first synthetic flow gate use prescribed inlet velocity,
free outlet velocity, `pumpDp=0`, fixed timestep, and `relaxMass=1`. Retain the accepted
momentum relaxation. Keep `rotationalPressureCorrection=false`: its extra pressure
update (`P:10610`) is outside the `p_new=p+phi` map derived here.
These are proposed gate settings, not a reconstruction of an
undisclosed production command.

Exclude `--pressure-k`, `--rc-only`, `--rc-blend`, `--pspg`, `--flux-neumann`, prescribed
opening-flux sources, and `MARS_HYPRE_USE_DDT` from this first implementation. Keep the
full VMS gradient difference, with `MARS_VMS_COMPACT_ONLY` unset. Reject unsupported
combinations explicitly when the new outlet mode is enabled. Pressure-driven inlets,
additional velocity projections, and other pressure paths need their own extension.

| Source fact at this snapshot | Consequence |
|---|---|
| Driver `:435–438` selects enum `DDT` without `--pressure-k` | The enum alone does not identify the solved matrix. |
| `P:10028–10037` selects `Apre` when `MARS_HYPRE_USE_DDT` is absent | The selected Hypre path solves a modified stiffness approximation. |
| `P:4978–5054` assembles compact RC sensitivity into `d_valuesPre` | This is not the full derivative of the temporal velocity correction. |
| `P:11180`, driver `:1574` refresh the trace before the timestep | This supplies a useful once-per-step trace lifetime. |
| `P:11265` is the only trace consumer, via boundary reporting | The trace does not yet enter the predictor, pressure solve, or corrector. |

Here `P` is `backend/distributed/unstructured/fem/mars_ns_pump_solver.hpp`, `VMS` is
`backend/distributed/unstructured/fem/mars_vms_pressure_stab.hpp`, and `driver` is
`examples/distributed/unstructured/mars_pump.cu`. All line numbers refer to the snapshot.

## 1. Unknowns, signs, and boundary quadrature

For this constant-density, fixed-geometry formulation,

\[
\partial_t u+(u\cdot\nabla)u-\nu\Delta u+\rho^{-1}\nabla p=f/\rho,
\qquad \nabla\cdot u=0.
\]

`p` and its increment `phi` are physical pressure. `M=diag(V_i)`, repeated for the
three velocity components, contains geometric dual volumes. Let
\(h=\Delta t_{\rm eff}/\rho\): BDF1 uses `dt`, active constant-step BDF2 uses `2*dt/3`
(`P:10311`). The influence coefficient `D_f` in the volumetric RC flux has units
time/density. It is not necessarily `h`: momentum diagonal and relaxation matter.

Let `Q` be the velocity correction mask: zero on prescribed velocity components,
one on free components. The first implementation excludes extra opening-normal
projections. If such a projection is later enabled, its actual action belongs in `Q`.

An interior SCS flux is positive from its left to right node:
\(q_f=\tfrac12(u_L+u_R)\cdot A_f\), scattered `+q_f` to L and `-q_f` to R.
`B_i` denotes that integrated velocity-divergence matrix (`P:2396`). Exterior area
vectors point outward. Continuity accumulators have units volume/time; divide by
dual volume only for normalized diagnostics, not during conservative assembly.

**Choose nodal, lumped triangle quadrature for this SCS extension.** For each outlet
triangle `f=(a,b,c)`, use three samples `b=(f,r)`, one at each vertex `r`, each with
vector area `A_b=A_f/3` and scalar area `a_b=|A_f|/3`. A boundary sample scatters its
flux to its own vertex control volume. Then

\[
(B_o u)_i=u_i\cdot a_{o,i},\qquad
a_{o,i}=\sum_{f\ni i} A_f/3.
\]

This selects the local continuity scatter; a total triangle flux alone does not
specify it. It agrees with the existing raw SCS opening term (`P:2887,9192`).
Prescribed inlet/wall exterior fluxes are separate affine sources; their velocity
correction is zero. Do not add the old raw outlet term again after the new flux.

This is a defined MARS quadrature, not a claim of complete OpenAccel quadrature parity.
The FEM opening uses a different consistent surface mass matrix; see section 8.

## 2. Trace map and the two kinds of surface area

For one selected outlet patch, sample the **free volume pressure** at those same
vertices, `q=E p`. Thus `E` is restriction to the corresponding volume DOF, not
restriction to a clamped target or uncorrected opposing-node pressure. Require
`E*1=1` and affine-exact sampling at the declared sample positions.

\[
w_b=\frac{a_b}{\sum_c a_c},\qquad
F(p)=p_{\rm ref}\mathbf1+(1-\beta)(I-\mathbf1w^T)Ep,\qquad t=F(p).
\]

Consequently `w^T t=p_ref` and `t_b-t_c=(1-beta)(q_b-q_c)`. The entire mean is
prescribed; beta=0.05 retains 95% of the spatial fluctuation. It is not temporal
relaxation. The enabled mode accepts `0 <= beta <= 1`; the existing negative sentinel
continues to disable it. Do not merge independently prescribed outlet patches' means.

For nodal storage, aggregate the scalar weights as
\(a_i=\sum_{f\ni i}|A_f|/3\), separately from the vector `a_o,i` used for flux.
The current `OutletMomentFunctor` (`P:10929`) uses
\(|\sum_{f\ni i} A_f/3|\). That equals the scalar area only for aligned facet normals.
Its host trace tests supply areas directly and do not test this geometric construction.

Build scalar areas on device from uniquely owned facets; reverse-add to node owners,
then publish owner values if ghost readers need them. Reduce pressure moments over
owned nodes only. Empty local ranks contribute `(0,0)` and join every collective.
A globally absent/zero-area **requested** outlet is a collective setup error: a
fallback mean prevents division by zero but supplies no pressure boundary condition.
The scalar helper may retain its finite fallback for isolated algebra tests.

Different storage alone does not guarantee a nonuniform converged trace. If some
remaining enforcement imposes `Ep=t`, the fixed-point equation with beta>0 forces a
uniform trace. Remove all old outlet pressure masks/relocks in the enabled mode and
test a nontrivial response to tangential pressure forcing. Do not promise that the
profile survives merely because `d_pTraceOutlet` and `d_p` are separate arrays.

## 3. The boundary-aware nodal pressure gradient

The current transpose scatter (`P:1579`), volume division, and sign flip give
\(G_0p=-M^{-1}B_i^Tp\) in both predictor (`P:8503–8535`) and corrector (`P:10354–10383`).
Adding an outlet pressure trace requires the trace **difference**, not just `t*a/V`.

Let `B=B_i+B_o`. Let `Z` map sample trace values to nodal vector force:
\((Zt)_i=\sum_{b\to i}t_b A_b\). Use

\[
\boxed{g(p,t)=M^{-1}(-B^Tp+Zt)=G_vp+G_t t.}
\]

With nodal trace storage the concrete change is

\[
\boxed{g_i(p,t)=(G_0p)_i+(t_i-p_i)a_{o,i}/V_i.}
\]

Apply `Q` when this gradient changes velocity. Other prescribed-velocity boundary
fluxes remain affine sources. Outlet area alone is not the entire exterior area of
a dual cell touching an inlet or wall.

Two exact algebra gates, using this same quadrature, are

\[
Z\mathbf1=B_o^T\mathbf1,\quad g(c\mathbf1,c\mathbf1)=0,
\qquad u^TMg=-p^TBu+t^TZ^Tu.
\]

The last identity checks pressure work and the boundary adjoint. It does not prove
full time-discrete energy stability or affine reproduction by the existing midpoint
SCS nodal gradient. The affine compact-reconstruction test below is a separate test.

### What trace freezing does and does not remove

For `p_new=p+phi` and trace change `delta t`,

\[
\delta g=G_v\phi+G_t\delta t,\qquad \delta u=-hQ\delta g.
\]

If the same trace is used before the predictor and through the corrector,
`delta t=0`. The frozen-trace increment gradient still contains

\[
\boxed{(G_v\phi)_i=(G_0\phi)_i-\phi_i a_{o,i}/V_i.}
\]

Thus no known **trace-change** lift is needed within that step, but this unknown
boundary contribution remains. `phi=0` at outlet volume nodes would remove it
incorrectly. In the current scaffold neither gradient reads the trace: identical
trace use in predictor/corrector is a contract to implement, not an existing result.

## 4. Complete flux residual and its exact correction derivative

At an interior SCS retain the existing VMS flux, coefficient interpolation, and
velocity-BC masking of the reconstructed-gradient term (`VMS:85–209`). Freeze its
predictor gradient `bar g` and coefficients during the inner correction iteration.

For the selected outlet triangle, define a concrete compatible local scatter:

\[
q_{f,r}=u_r\cdot(A_f/3)+D_f\left(\bar g_f-\nabla p_f^{\rm mix}\right)\cdot(A_f/3),
\]
\[
\nabla p_f^{\rm mix}=p_o\nabla N_o+\sum_{r\in f}t_r\nabla N_r,\qquad
\bar g_f=\tfrac12\left(\tfrac13\sum_{r\in f}\bar g_r+\bar g_o\right).
\]

Here `o` is the opposite tetrahedron node. Choose `D_f` as the average of the three
face-node coefficients, without an opposing-node coefficient. This corrects the
face/opposite diffusivity blend in the current diagnostic (`P:11093–11099`).
The reconstructed-gradient blend above matches its existing boundary interpolation;
it is a declared boundary interpolation, distinct from the interior velocity-BC mask.
Sum `q_f,r` to obtain the triangle flux; scatter each to row `r`. The advective sum
is exactly `mean(u_a,u_b,u_c) dot A_f`. Retain this same local flux for diagnostics.

`D_f` is written as scalar for the current shared-diagonal case. For componentwise
coefficients each component is multiplied before its area contraction. Use the
same coefficient values and epoch in residual and pressure derivative, including
BDF startup. Do not use a boundary identity diagonal as a physical momentum diagonal.

Stack interior and outlet samples. Let `C` scatter sample fluxes to continuity rows,
`F_u` map velocities to oriented sample flux, `W` map reconstructed gradients to
oriented samples, `T=diag(D_f)`, and `L_v p+L_t t` be their compact pressure-gradient
area contractions. Then `C F_u=B` and the integrated volumetric residual is

\[
R(u,p;\bar g,t)=C\{F_u u+T[W\bar g-L_vp-L_tt]\}+b_{\rm prescribed}.
\]

Interior samples have no trace columns. Outlet compact `L_v` retains only the
opposite-node pressure; `L_t` contains the face trace values. The **absolute** term
`-C T L_t t` stays in the residual even when `delta t=0`.

With `bar g`, `t`, `T`, geometry, boundary masks, and flux history fixed,

\[
R^+=R+hJ\phi,\qquad
\boxed{J=H_\Gamma+K_D/h,\quad H_\Gamma=BQM^{-1}B^T,\quad K_D=-CTL_v,}
\]
\[
\boxed{J\phi=-R/h.}
\]

The sign agrees with `buildPressureRhsKernel` (`P:2787`). `J` has the dimensions
of a stiffness matrix; pressure increments have pressure units. This is the exact
Jacobian of the declared **Chorin correction map**, not an exact Schur complement
using the inverse implicit momentum matrix.

The compact outlet partial alone is
\(\partial q_{f,r}/\partial p_j=-D_f(A_f/3)\cdot\nabla N_o\,\delta_{jo}\).
It fits tetrahedral stiffness adjacency. The full derivative additionally includes
the velocity response in `H_Gamma`; cross-element nodal gradients can enlarge its
support. Existing compact CSR entries do not prove the full derivative fits.

## 5. Practical correction iteration, gauge, and time lag

Use `Apre` as an explicitly **approximate** operator `A_tilde`, extended with
`K_D,out/h = -C_out T_out L_v,out/h` and with its obsolete outlet constraints removed.
The flux partial `dq/dp` must therefore be scaled by `rho/dtEff` before this matrix addition.
It approximates `J`; do not claim equality or accept its linear residual as continuity.
The existing VMS repeated-correction loop (`P:11885–11913`) is an entry point to adapt,
not proof that a fixed number of passes converges.

For fixed timestep data:

1. Build `t=F(p_old)` before the predictor. Evaluate the new boundary-aware predictor
   gradient and keep its halo-complete `bar g` for the inner corrections.
2. After implicit diffusion, initialize `u_iter=u**`, `p_iter=p_old`. Build the full
   residual above from the current iterate, including all outlet continuity rows.
3. Solve `A_tilde * delta_phi = -R/h` using the existing Hypre GMRES interface.
4. Apply the same damping to both updates:
   `p_iter += omega*delta_phi`, `u_iter -= h*Q*G_v*(omega*delta_phi)`.
   Preserve prescribed velocities and update owner/ghost state before reading fluxes.
5. Rebuild the actual frozen-state residual. Continue only with measured residual
   contraction; terminate successfully only at the specified residual/balance tolerance.
   Preserve BDF history once per physical step, not once per correction.

For fixed data the residual iteration is
\(R_{k+1}=(I-\omega J\widetilde A^{-1})R_k\).
Require small synthetic spectral/contraction gates before the flow gate. Damping
is an algorithmic correction parameter, distinct from beta, `relaxU`, and `relaxMass`.
Do not choose it from the outlet profile blend or assume damping guarantees stability.
If this iteration cannot contract, use a verified true-`J` operator/Krylov path before
running the new outlet; the incomplete assembled Gram or experimental coupled FGMRES
is not a shortcut. This is a linear pressure-correction iteration at frozen momentum
data; it does not reassemble SIMPLE momentum equations. Record correction count,
damping, and stopping tolerances in comparisons. If those also change from the
historical run, do not attribute the entire difference to the outlet formula alone.

### Pressure level and solver assumptions

For the full boundary-paired gradient,

\[
\mathbf1^TH_\Gamma\mathbf1=
\|M^{-1/2}QB^T\mathbf1\|^2.
\]

A free outlet with nonzero normal correction therefore anchors the constant mode
through the temporal correction itself. The RC partial is not necessarily its only
anchor. With an interior-only gradient, the compact boundary term can instead remove
the constant vector. Neither argument proves absence of other null modes.

On a tetrahedron with outward face area, `A_f dot grad(N_o) < 0`, so the compact
outlet derivative has positive entries in the **opposite-node column** of face-node
rows. It is generally nonsymmetric, not simply a positive diagonal Robin penalty.
Test rank/small singular values of both the true `J` and `A_tilde` on the small gates.
Do not silently add a gauge pin or subtract the pressure/RHS mean: prescribed exterior
pressure already supplies a physical level. Any retained genuine pressure-inlet
condition requires its own elimination/lifting derivation outside the first gate.

Use Hypre GMRES for the selected path; disallow a forcing `MARS_HYPRE_KRYLOV=pcg`
unless symmetry, definiteness, and preconditioner assumptions are proved for that mode.
Refresh wrapped matrix storage and preconditioner setup when coefficients change.
If scaling is introduced, apply it consistently to operator, RHS, inverse solution map,
and preconditioner; it cannot repair a missing term or incorrect pressure mask.

### What a converged inner residual establishes

Once-per-step trace and reconstructed-gradient freezing establish only
`R(u_new,p_new;bar_g_old,t_old)=0`. A refreshed-gradient residual differs by
`C T W [g(p_new,t_old)-bar_g_old]`; the trace defect is `t_old-F(p_new)`.
Report those defects separately. BDF2 velocity coefficients do not automatically
make this explicit boundary lag second order; demonstrate timestep convergence.

If the intended end state requires a fully current trace and reconstructed gradient,
add an outer boundary/reconstruction refresh and converge both defects. Freeze them
again during each inner solve. With a known trace change and frozen `bar g`,

\[
J_t=-BQG_t-CTL_t/h,\qquad J\phi=-R/h-J_t\delta t.
\]

Equivalently apply `-h Q G_t delta_t` to velocity, switch the trace, rebuild `R`, and
solve; do not apply both lifts. If the reconstructed gradient changes during this
map, include its known change, or its derivative when treated implicitly:
`J_fresh = H_Gamma + C T (W G_v-L_v)/h`.
The fully implicit trace derivative `(1-beta)(I-1*w^T)E` introduces global coupling;
it is outside this first frozen-trace implementation.

## 6. Small gates with explicit expected results

The following are acceptance requirements, not claims about unrun CUDA kernels.

### Trace values and bent surface

Use samples `(2,4,8)` Pa, scalar areas `(1,2,1)`, `p_ref=10` Pa. The mean is `4.5`.

| beta | Expected trace |
|---|---|
| 0.05 | `(7.625, 9.525, 13.325)` |
| 0 | `(7.5, 9.5, 13.5)` |
| 1 | `(10, 10, 10)` |

Adding 100 Pa to every sample leaves the trace unchanged. Verify that building it
leaves volume pressure unchanged. Beta=1 recovers a uniform **side trace**, not
necessarily the legacy strongly constrained nodal solver.

For a bent-surface area test, give one sample two perpendicular unit-area contributions
and another one unit-area contribution. With sample pressures `(0,3)`, true weights
`(2,1)` give mean `1`; vector-norm weights `(sqrt(2),1)` give `1.242640687119285`.
At beta=0.05 the latter trace has physical mean `p_ref-0.230508652763321`.
This is a geometric-weight failure that the existing supplied-area script cannot catch.

### Compact boundary reconstruction and signs

Use the unit tetrahedron `x0=(0,0,0)`, `x1=(1,0,0)`, `x2=(0,1,0)`, `x3=(0,0,1)`.
Outlet face `(1,2,3)` has outward area `(1/2,1/2,1/2)`, opposite node `0`.
For `p(x)=7+2x-3y+5z`, nodal pressure is `(7,9,4,12)`.
Feed the **manufactured** trace `(9,4,12)` directly and supply reconstructed gradient
`bar_g=(2,-3,5)`. The mixed element gradient must be `(2,-3,5)`; with uniform velocity
`(1,2,3)`, the triangle volumetric flux is `3` and the RC difference is zero.

With velocity, reconstructed gradient, coefficients, and other pressures fixed:
an opposite-node pressure perturbation `delta` changes total flux by `+1.5*D*delta`;
a common trace perturbation changes it by `-1.5*D*delta`. For `D=2`, `delta=0.2`,
the answers are `+0.6` and `-0.6`. Repeat on translated, scaled, and skewed but
well-conditioned tetrahedra. Supply the analytic nodal gradient in this local gate;
test the actual reconstructed nodal-gradient operator separately.

Do not feed this affine tangential trace through beta=0.05 and demand the same affine
solution: that map deliberately modifies its variation.

### Adjoint, Jacobian, and a correction-convergence counterexample

Check `u^T M g + p^T B u - t^T Z^T u = 0` and `g(c,c)=0` on the assembled small
operator. With `bar_g,t,D,Q` frozen, finite-difference **both** pressure and velocity
along `(delta_p,delta_u)=(v,-h Q G_v v)`. Compare with `h J v`. A pressure-only
perturbation tests the compact partial, not the Chorin correction Jacobian.
Use `epsilon=1e-3,1e-4,1e-5` after nondimensionalizing a small linear fixture; seek
FP64 relative error below `1e-8`, using a nonzero operator-action scale.

Here is an independently calculated small counterexample to automatic contraction
of a compact approximate solve. It is a matrix fixture, **not a physical channel run**.
Use the unit tet, midpoint SCS areas from `mars_cvfem_tet_area.hpp`, outlet `(1,2,3)`,
the nodal surface scatter above, `M=I/24`, and `Q=I`. Then

\[
H_\Gamma=\frac1{48}
\begin{bmatrix}33&5&5&5\\5&77&-1&-1\\5&-1&77&-1\\5&-1&-1&77\end{bmatrix},\quad
K=\frac1{6}
\begin{bmatrix}3&-1&-1&-1\\-1&1&0&0\\-1&0&1&0\\-1&0&0&1\end{bmatrix}.
\]

Set nondimensional `h=1`, `D=0.3`, and `K_b[1:4,0]=0.15`, all other entries zero.
`K_D=0.3*K+K_b`, `J=H_Gamma+K_D`, while the compact approximation is
`A_tilde=1.3*K+K_b`. Both are nonsingular here, but eigenvalues of `A_tilde^-1 J` are
approximately `(1.0354710585,13.0478622748,7.7307692308,7.7307692308)`.
The undamped iteration has spectral radius `12.0478622748`; even omega=0.25 gives
`2.2619655687`. Thus neither removing the constant mode nor adding the correct compact
boundary partial guarantees contraction. These numbers do not prescribe an omega
for another mesh. Non-normal operators also require measured residual histories.
This counterexample concerns stationary correction with `A_tilde^-1`; it does not
show that `A_tilde` is unusable as a preconditioner for GMRES solving the actual `J`.

### Ownership and continuity

Partition the three-sample trace gate into moments `(10,2)`, `(8,2)`, and `(0,0)` on
three ranks. Every rank must obtain mean `4.5`. Test one/two/four-rank meshes with
shared outlet nodes, an empty-facet rank, and valid zero-local-work cases.

Require, before any pressure equation replacement,
`sum_owned(R_i) = sum_unique_exterior_samples(q_b)`.
Use the same accepted stabilized flux, coefficient epoch, and prescribed inlet source
on both sides. For small FP64 fixtures, target an error below `1e-11` times the sum
of absolute flux contributions, with a nonzero physical scale for the zero-flux case.
Compare local residual vectors across partitions after matching global DOFs, not only
their global sum. Check constant-mode response and smallest singular values too.

For scalar trace and local affine algebra on unit, well-conditioned fixtures, use
FP64 `1e-12` or FP32 `1e-5` times the pressure/term scale. Scale tolerances for geometry
conditioning and accumulation count; these are proposed acceptance thresholds.

### Manufactured channel and temporal gates

For `0<x<L`, `-H<y<H`, compatible periodic/slip span boundaries, and Laplacian viscosity,
use `u_x=G*(H^2-y^2)/(2*mu)`, `u_y=u_z=0`, `p=p_out+G*(L-x)`.
Per unit span, `Q=2*G*H^3/(3*mu)`. Choose `H=1`, `L=4`, `rho=mu=G=1`, `p_out=0`:
`Umax=0.5`, `Q=2/3`, and inlet pressure `4`. Prescribe its parabolic inlet velocity
and no-slip walls, with the new pressure outlet. Use physical viscosity explicitly.

If viscosity uses full symmetric stress, prescribe its compatible tangential opening
traction too; pressure-only traction is not that exact channel solution. For the
componentwise Laplacian, `partial_n u=0` at the opening is compatible.

P1 velocity cannot represent the quadratic profile exactly. Check mesh convergence,
not roundoff equality to the continuum solution; propose less than 1% flow and velocity
L2 error by a resolved shape-regular mesh with `h_mesh/H <= 1/32`, decreasing over the
last two refinements. Report achieved orders rather than assuming them. Converge
the discrete inner residual and boundary/reconstruction defects before interpreting
spatial errors. For a steady gate target normalized fixed-point defects below `1e-8`.
Add a manufactured tangential-pressure forcing gate for a nonuniform trace response.
For temporal order, use a **time-dependent** manufactured state with a changing
mean-free outlet profile on a fixed, sufficiently resolved mesh. For example choose
sampled volume pressure `p_i(t)=p0_i+a*sin(t)*y_i`, build its exact current-time trace
with `F`, and derive compatible forcing from the fresh-trace equations. If using a
semidiscrete manufactured problem with an artificial continuity source, label it and
include that source in its balance identity. Do not construct the forcing using the
same stale trace as the algorithm, which would hide the lag being tested. Refine
`dt,dt/2,dt/4` at the same physical final time. A steady-channel timestep sweep alone
cannot establish the order of the boundary lag.

## 7. CUDA/MPI implementation boundaries

- Cache unique facet ownership, scalar areas, vector areas, opposite-node mapping,
  trace values, coefficient fields, and accepted sample flux on device. Do not add
  full host copies or default mesh-identifying output.
- Interior owned-element and boundary owned-facet scatters must reverse-add their
  contributions exactly once. The existing owner-complete nodal opening addition
  runs after reverse-add; do not also execute it for facets already scattered.
- Both mean reduction and fatal setup status must follow the same collective sequence
  on every rank. Review the local storage-size early return in `updateOutletPressureTrace`
  (`P:11183`) when supporting empty ranks; guard zero-block launches without skipping MPI.
- Remove outlet identity rows, column elimination, RHS zeroing, and post-solve pressure
  relocking together, including halo copies of masks. Keep genuine other BCs unchanged.
- Assemble the approximate compact boundary derivative before copying values into
  the active wrapped matrix; do not update only the source allocation. Initialize/rebuild
  the Hypre operator and preconditioner from the resulting values.
- Residual evaluation must read an explicit iteration state and must not advance flux
  history. `divMaxVmsOwned` is not the acceptance residual (`P:10178` onward); its existing
  scalar-coefficient/RMS defects remain a separate implementation issue to address where
  this new mode obtains its convergence test.
- Review `runCorrectorStep` and the VMS repetition loop together: cumulative pressure,
  the corrected velocity fed into the next residual, damping, owner/ghost publication,
  and BDF history must all use the same iteration. A pass count is not convergence.
- Follow the ultracode gate before implementing these interacting changes. This
  document does not change the model/effort setting or authorize claims of GPU validation.

## 8. Why the FEM opening helper is not the SCS boundary operator

`mars_fem_projection.hpp:125` uses
\((E_f u)_i=A_f\cdot(2u_i+u_j+u_k)/12\).
Its corresponding frozen-trace gradient addition is
\(-A_f(2\phi_i+\phi_j+\phi_k)/(12V_i)\), not the SCS nodal `-phi_i*a_o,i/V_i`.
`femOpeningSurfaceActive` (`P:9131`) requires FEM projection, tet elements, no prescribed
opening source, and `nnzFemGram==0`. Turning that helper on in the selected SCS path
would mix quadratures. Extend each family with its own paired operators when requested.

## Validation performed for this specification

GPT reran Claude's committed `scripts/outlet_trace_check.py` from the frozen snapshot:
all printed host gates passed, exit code 0. It is an algebra replica, not execution
of `outletTraceValue` or a GPU/MPI test.

GPT evaluated the explicit trace, bent-area, affine boundary, and small-matrix numbers
above with NumPy arithmetic. A separate GPT review independently reconstructed the
six SCS area vectors and confirmed `H_Gamma` and the counterexample eigenvalues.
The source, units, adjoint derivation, and frozen-state Jacobian received independent
reviews. No new solver/test source was written; no GPU, MPI, confidential-case, or
performance run was made. The stated acceptance gates remain implementation work.
