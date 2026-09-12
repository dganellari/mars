# OpenAccel-to-MARS numerical contract, version 1

Author: GPT/Codex. Date: 2026-09-12.
Status: source-derived architecture; ready for the reference-harness task below.
No OpenAccel simulation or new GPU kernel was executed for this document.
Implementer: Claude. GPT reviews the numerical interfaces and bounded diffs.

This implements Stage 0 of the [port plan](gpt_openaccel_gpu_port_plan_2026-09-12.md).
It defines a separate SIMPLE solver, followed by the SIMPLEC variant. It does
not redefine the existing Chorin/BDF2 solver. All fixtures and source material
here are public; private application results are outside this contract.

## 1. Reference and first configuration

Reference: [OpenAccel at 0d69041](https://github.com/CCFNUM/OpenAccel/tree/0d69041ba1afda63e9e4328d9e0d9834bba37756),
full revision `0d69041ba1afda63e9e4328d9e0d9834bba37756`.
The local checkout was clean during inspection. Its `src/solver` gitlink is
[LibLinSolve at e351ba5](https://github.com/CCFNUM/LibLinSolve/tree/e351ba5eeaf3537dcc53d0aba09a8347f0a44cd0),
full revision `e351ba5eeaf3537dcc53d0aba09a8347f0a44cd0`.
That dependency was fetched separately for the residual audit; the OpenAccel
checkout was not modified. Pin both, including the vendored Nalu geometry.

The [machine-readable contract](../../tests/data/public_openaccel_reference/contract_v1.json)
records profiles, analytic fixtures, required exports and inspected source hashes.
It is harness input, **not OpenAccel's native input syntax**. Claude must translate
it through the pinned parser and export the effective settings from the running
reference. A native deck that silently ignores an option fails the harness.

The PI's exact revision, SIMPLE/SIMPLEC choice, steady/transient settings,
advection and residual criterion remain unknown. The following is our explicit
public reference profile, not an assertion about the PI's deck:

| Setting | Initial profile |
|---|---|
| Domain | One fixed, nonrotating 3D Tet4 domain; no interfaces, AMR or periodicity |
| Fluid | Incompressible, laminar, constant density and dynamic viscosity |
| Primary time mode | Steady, physical timescale 0.01 s used only as the reference's false-transient diagonal |
| Algorithm | `consistent_=false`, `fractionalStepMethod_=false`; one momentum and one pressure solve per outer iteration |
| Advection | `upwind`; high-resolution blending zero; NSO off |
| Interpolation | Velocity `trilinear` (unshifted); pressure `linearLinear` (shifted); both gradient-interpolation settings `linearLinear` |
| Relaxation | Velocity 0.3; pressure 0.3; stored mass flux 0.75; `relaxGradients_=false` |
| Gradient switches | Incremental gradient change on; limiter off; no symmetry surfaces |
| Other controls | No field-update clipping or acceleration, body forces, sources, turbulence, buoyancy or compressible extensions; retain reference boundary-flux clipping |
| BCs | Positive inward normal-speed inlet; stationary no-slip walls; average static-pressure outlet, reference pressure 0 Pa, blend 0.05 |
| Diagnostic variants | Outlet blend 1; SIMPLEC; gradient relaxation enabled; BDF1 and BDF2 are separately identified comparisons |

Use two property cases: `rho=1, mu=0.1` for elementary checks and water
`rho=1000 kg/m^3, mu=0.001 Pa s`, hence `nu=mu/rho=1e-6 m^2/s`.
The viscosity passed to the reference stress is **mu**, not nu. Use inlet speed
0.1 m/s for the first public fixture; water application parameters are a later
explicit profile. No mesh-independent Reynolds number is assigned.

The first frozen-state and single-iteration gates need no claim of nonlinear
convergence. The implementation must eventually support a converged steady
public case and physical-time iterations. Upwind parity does not certify the
reference's high-resolution limiter or the PI's exact configuration.

## 2. Unknowns, units and quadrature

The target continuum equations, with no sources, are

\[
\partial_t(\rho u)+\nabla\cdot(\rho u\otimes u)
 =-\nabla p+\nabla\cdot[\mu(\nabla u+\nabla u^T)],\qquad
\nabla\cdot(\rho u)=0.
\]

Pressure is physical pressure in Pa. Node values are `u_i` in m/s, `p_i` in Pa,
and reconstructed gradients `g_i=grad(p)_i` in Pa/m. The nodal control volume
`V_i` is in m^3. Each tetrahedron owns six interior subcontrol-surface samples;
each exterior triangular facet has three samples. Store **mass** flux `m_s`
in kg/s. A momentum matrix entry has units kg/s and its RHS units N. Pressure
matrix entries have units (kg/s)/Pa. The influence coefficient `d_ij` has
units m^3 s/kg; `d_ij g_ij` is a velocity.

Use linear tet shape functions and constant physical derivatives `b_a=grad N_a`.
The finite-volume unknowns are collocated at the vertices of a median-dual
control volume; this is not an instruction to substitute a Galerkin mass/stiffness
pair. The local dual volumes are V_tet/4 for an affine Tet4. Assemble nodal
volumes over all incident elements exactly once.

The interior directed pairs in Nalu `TetSCS` are `(0,1),(1,2),(0,2),(0,3),(1,3),(2,3)`.
For a sample from L to R, scatter `+m_s` to L and `-m_s` to R. Exterior area
vectors point outward and scatter to their nearest face node. Derive area vectors
from the actual reference geometry routines, then compare MARS by oriented
sample identity. Coordinate order and a comment saying "edge midpoint" are
insufficient evidence.

| Location | Unshifted shape weights | Shifted weights |
|---|---|---|
| Interior SCS joining L,R | 13/36 at L,R; 5/36 at each other node | 1/2 at L,R; zero elsewhere |
| Triangular boundary sample nearest r | 11/18 at r; 7/36 at the other two face nodes | 1 at r; zero elsewhere |

On a planar Tet4 face each boundary area is A_face/3. Its **sample interpolation**
is still the table above. In particular, velocity and `d` at an unshifted boundary
sample are not the mean of the three face-node values. Coordinate interpolation
and field interpolation have separate switches in the reference; export both.
Tet derivatives are constant, but this does not erase differences in sampled fields.

The reference geometry is in `external/nalu/src/master_element/Tet4CVFEM.C`
(`TetSCV`, `TetSCS`) and `MasterElement.C` (`Tri3DSCS`, lines 610–714).

## 3. Reconstructed gradient is a separate operator

For a scalar nodal field q, no symmetry or interfaces, and incremental assembly
enabled, write the reference reconstruction as

\[
(Gq)_i={1\over V_i}\left[
\sum_{s\ni i}\sigma_{is}(I_s q-q_i)A_s
+\sum_{b:r(b)=i}(I_b q-q_i)A_b\right].
\]

Here `I` uses **that field's interpolation scheme**, and boundary q is gathered
from the volume nodal field. It is not the independently prescribed pressure
trace. With shifted pressure interpolation the boundary increment is zero at
each nearest-node sample; the interior contributions remain. Constant q must
give zero. Do not require this reconstruction on a boundary-truncated control
volume to reproduce every affine gradient: export and match the actual operator.
Compact element differentiation of affine q is exact and is a different gate.

`nodeField::updateGradientField` (`nodeField.hpp:3470–4896`) first uses the full
reconstruction. On later calls it applies `(1-eta) g_old + eta Gq`; `eta=1` in
the primary profile. The optional reference default is `eta=0.5625`. Pressure
correction always has `eta=1` (`fieldBroker::setupPressureCorrection`). Its
interpolation matches pressure. Symmetry projection and limiting must be off or
separately implemented; neither may be inferred from this formula.

**Do not insert the Chorin outlet boundary-gradient lift into this operator.**
OpenAccel's mixed trace enters the compact boundary flux in section 5. The
pressure source and correction use this separately reconstructed nodal gradient.

## 4. Momentum increment, full stress and influence coefficients

Let `m^k` be the stored transport flux entering outer iteration k, and
`c_i^k=sum sigma*m^k + sum boundary*m^k`. This is an integrated mass imbalance
in kg/s, not a volume-normalized divergence. For fixed density/mesh the
reference `updateMassDivergenceField` computes precisely this sum.

The momentum equation solves an **increment**: `A_hat delta_u = b_u`, then
`u_star=u^k+delta_u`. The solution vector starts at zero. The momentum assembly
holds `m^k`, reconstructed gradients, properties and boundary state fixed.

For each component the unrelaxed residual contains:

\[
R_{u,i}=T_i+\sum_s\sigma_{is}
 \{m_s^k u_{up,s}-\mu_s(\nabla u+\nabla u^T)_s A_s\}
+R_{u,i}^{boundary}+V_i g_i^k-c_i^k u_i.
\]

This last `-c_i u_i` is present in the pinned reference. It converts the
transport operator when the stored flux is not locally conservative. Its RHS
contribution is `+c_i u_i`; its diagonal addition is `max(-c_i,0)`. Do not omit
it merely because the final solution should be conservative.

For upwind, `u_up=u_L` if m>=0, otherwise u_R. Frozen convection adds
`max(m,0)` in the L column and `min(m,0)` in the R column of the L row, with
the opposite contributions to the R row. In the high-resolution extension the
RHS uses `u_up+blend_up*(grad u_up dot (x_s-x_up))`, while the implicit part
remains the upwind coefficients. The limiter itself is outside version 1.

The viscous L-row block for node a is

\[
(A_{La}^{visc})_{ij}=-\mu_s[\delta_{ij}(b_a\cdot A_s)+b_{a,i}A_{s,j}],
\qquad A_{Ra}^{visc}=-A_{La}^{visc}.
\]

Both terms are implicit in `navierStokesAssemblerElemTerms.cpp`. This requires
3x3 component blocks, or their mathematically equivalent scalar-CSR expansion.
Three independent scalar diffusion solves do **not** reproduce this matrix.
The compressible `2/3 div(u)` stress branch is inactive here.

For steady operation, add `rho*V_i/t_pseudo` to the diagonal only; do not add a
physical-time term to the RHS. For transient operation,

\[
T_i={\rho V_i\over\Delta t}
(\gamma_0 u_i+\gamma_1 u_i^n+\gamma_2 u_i^{n-1}).
\]

BDF1 has `(1,-1,0)`. BDF2 uses the pinned `BDF2::coeff(dt,dt_previous)`;
for constant dt it is `(3/2,-2,1/2)`. Export the reference startup dispatch
and histories before enabling BDF2; do not copy MARS's startup state machine
without that comparison. Histories move once per physical step, not per outer
iteration. Momentum diagonal, RHS and d must use the same time coefficients.

Assembly order is mandatory:

1. Form interior, node and boundary blocks and `b_u=-R_u`.
2. Apply reference constraints (none in the first single-domain profile).
3. Divide each component's scalar diagonal by `alpha_u`. Other block entries
   are unchanged (`phiAssembler::postAssemble`, `assembleRelaxation_`).
4. Compute, from this diagonal `a_hat_ij`, `d_ij=V_i/(a_hat_ij+SMALL)`.
5. Multiply boundary-node RHS components by 0.75 **once per node**. In steady
   mode all the selected physical boundary nodes participate; in transient
   mode only no-slip wall nodes do. This operation does not change the matrix
   (`assembleBoundaryRelaxation_`). Then apply symmetry conditions if enabled.
6. Solve the full velocity increment with explicit field relaxation 1.

Do not multiply d by alpha_u again. Do not add the usual absolute-variable
under-relaxation RHS formula to an already assembled increment residual.

For SIMPLEC, the reference additionally forms

\[
\widetilde d_{ij}=V_i/[a\_hat_{ij}+\sum_{a\ne i}(A_{ia})_{jj}+SMALL].
\]

The sum excludes the entire diagonal node block; cross-component off-diagonal
entries do not enter this coefficient. SIMPLEC uses d_tilde in the pressure
**matrix** and velocity correction, while fresh/stored flux evaluation still
uses d. This is a distinct parity profile. Export denominators and diagnose
nonfinite values; no silent coefficient clipping.

## 5. Continuity, pressure matrix and boundaries

Use `S` for signed sample-to-node scatter, `I_u` for velocity interpolation,
and `B` for compact element differentiation. For an interior sample,

\[
\widehat m_s=\rho\sum_j A_{s,j}
\left[(I_u u^\star)_{s,j}-d_{s,j}
\left(\sum_a b_{a,j}p_a-\tfrac12(g_{L,j}^k+g_{R,j}^k)\right)\right],
\quad d_{s,j}=\sum_a N_a^u(s)d_{a,j}.
\]

The arithmetic adjacent-node gradient average is the constant-density path.
Do not replace it with interpolating all four nodal gradients. Build the
pressure RHS from this newly evaluated flux, not from relaxed `m^k`:

\[
R_c=S\widehat m+q_{boundary},\qquad H\phi=-R_c.
\]

For SIMPLE, an interior L-row entry is

\[
H_{La}=-\rho\sum_j A_{s,j}d_{s,j}b_{a,j},\qquad H_{Ra}=-H_{La}.
\]

Use interpolated d_tilde in this entry for SIMPLEC. There is no extra rho/dt
factor: momentum coefficients already determine d. Pressure post-assembly
applies constraints but **skips the base diagonal under-relaxation**.

This H is the reference's compact SIMPLE/SIMPLEC pressure approximation. It is
not the exact Jacobian of the full iteration, and is not MARS's Chorin true-J.
For an actual frozen block system `[A G; D C]`, exact elimination gives
`C-D A^-1 G`; the reference's component-diagonal approximation, compact gradient,
interpolations and lagged reconstructed gradient must be kept distinct. Do not
assert `G=-D^T`, pressure symmetry or positive definiteness from continuum
integration by parts. Export matrix actions after BCs; use GMRES for a
nonsymmetric system and FGMRES if preconditioner actions vary.

### Inlet and wall

The reference normal-speed inlet forms a nodal side vector by summing **inward**
sample area vectors and normalizing it to the specified speed. It then
interpolates those nodal side vectors to face samples using velocity's shape
functions (`velocity::updateBoundarySideFieldNormalSpeed`). On a planar inlet
this is exactly inward normal speed. On a bent patch a shared node has one
averaged direction, so do not assert every face-sample vector has the specified
normal direction and magnitude. Retain the separate side and volume fields.

At an active specified-velocity inlet the continuity source is
`rho*u_bc,b dot A_b`; it has no pressure derivative. The momentum boundary
kernel uses prescribed side velocity for advective RHS, substitutes nodal-side
velocity at the face nodes of its local viscous workspace, and zeroes those
face columns of the boundary derivative. It does not globally replace all
volume velocity rows with identity rows.

At a stationary wall the normal mass flux is zero. The laminar momentum wall
term is a tangential coefficient law, not MARS's hard nodal no-slip mask:

\[
R_b^{wall}=c_b(I-nn^T)(I_b^u u-u_{bc,b}),\quad
c_b={\mu_b|A_b|\over 0.25 y_b},\quad
y_b=\max[-n\cdot(x_o-x_r),SMALL].
\]

Here o is the tet opposite node and r the nearest face node. The matrix is
the derivative with respect to face u. The reference updates wall coefficients
in initialization/post-solve; they are constant for this fixed laminar profile.

At a corner, add the relevant boundary-facet contributions to the same volume
row. Do not impose a new "opening overrides wall" global mask. Normalized side
values are boundary data, while volume values remain unknowns. The reference
also has a shared nodal-side buffer: the harness must export it after each
boundary update, along with boundary order, to expose any write-order effect
at patch intersections. Version 1 orders `walls`, then `inlet`, then `outlet`.
Do not silently replace an observed order dependence with a new rule; report it
before the production boundary implementation. All facets still contribute.

### Average-pressure outlet

At pressure pre-solve, estimate the pressure mean over **unflagged** outlet
samples using pressure interpolation and scalar sample areas:

\[
\bar p={\sum_{b\ active}|A_b|(I_b^p p)\over\sum_{b\ active}|A_b|},\qquad
t_b=p_{ref}+(1-\beta)(p_{r(b)}-\bar p).
\]

The estimate is interpolated; the trace update uses the nearest node. The
nodal-side pressure is the scalar-area weighted average of incident trace
samples. Source: `updatePressureBoundarySideFieldAverageStaticPressure_`
in `flowModel.cpp:20761–21060` and `nodeSideField<scalar,1>::interpolate`.
The primary shifted-pressure profile makes the estimate a nearest-node
quadrature. Do not replace scalar area weights by the norm of summed normals.

The boundary compact gradient gathers volume p at the opposite node and
**nodal-side pressure** at face nodes:

\[
(Bp)_{b,j}=p_o b_{o,j}+\sum_{a\in face}t^{node}_a b_{a,j}.
\]

For boundary sample r, the reconstructed contribution is
`(g_r+g_o)/2`, not `(mean_face(g)+g_o)/2`. Both velocity and d are interpolated
with the three **sample-specific** velocity shape weights, using face nodes
only. Thus

\[
\widehat m_b=\rho\sum_j A_{b,j}
[(I_b^u u^\star)_j-d_{b,j}((Bp)_{b,j}-(g_{r,j}^k+g_{o,j}^k)/2)].
\]

With the nodal-side trace frozen for this pressure solve, the compact derivative
in its nearest-node row and opposite-node column is

\[
H_{r,o}\mathrel{+}=-\rho\sum_j A_{b,j}d^{LHS}_{b,j}b_{o,j}.
\]

The boundary face columns are removed by the workspace BC multiplier. Do not
zero the entire pressure row. The pressure-outlet momentum term uses stored
outflow mass flux times nearest-node velocity, plus the tangential projection
of the full viscous traction: `-mu*(I-nn^T)*(grad u+grad u^T)*A_b` in the residual.
This is `IP_ZERO_NORMAL_STRESS__`, not a blanket zero viscous flux.

### Reversal, active area and gauge

Reversal is part of the pinned reference method, not a logging convention.
After updating stored flux, `updateFlowReversalFlag` first saves the net face
mass flux and clips individual outlet samples to `max(m_b,0)`. It treats a face
with negative saved net mass flux (`sum_b m_b < -SMALL`) as a slip wall,
sets its sample flags and zeroes all its stored boundary fluxes. A flagged face
is reopened when the arithmetic face-average velocity has nonnegative outward
normal component and the average side pressure is no larger than the average
volume face pressure. Otherwise its flux stays zero. These tests use the
pre-velocity-correction state at that point in the iteration. Pressure-specified
inlets have a corresponding reversal branch with inflow-only sample clipping;
our specified-velocity inlet has zero-gradient pressure and does not enter that
branch. Do not add that clipping to the selected inlet type.

Flagged outlet samples contribute neither the open pressure-boundary matrix/RHS
nor the open momentum terms; their trace update is skipped. The trace average
also excludes them. Do not allow negative outlet flux without this policy and
call it reference parity. Version 1's assembled fixtures start unflagged;
transition tests explicitly exercise flagging and reopening.

In an open reference domain no extra pin/mean subtraction is introduced. The
compact boundary entries break the constant nullspace. In a closed domain the
reference pressure-level routine zeros a selected row except for its original
diagonal and sets its RHS to zero. Gate closed domains separately.

Two defensive MARS behaviors are required and must be identified as such:
collectively reject a missing required CSR entry, and stop clearly if no active
pressure outlet remains in this open-domain profile. The source mean divides
by active area without a zero-area fallback here; do not reproduce a NaN or
invent a pressure pin/another BC to continue. Extra model handling requires
a separate decision. A nonzero constant response is necessary, not proof
that every possible assembled pressure matrix is nonsingular.

## 6. One complete iteration and field lifetimes

`segregatedFlowEquations::solve` defines this order. A pressure sub-iteration
is not a physical timestep. The first profile has one of each per outer iteration.

| Step | Reads | Writes / ordering requirement |
|---|---|---|
| 0. Start outer iteration | Accepted u,p,m, their gradients and side state | Snapshot previous-iteration u,p; refresh velocity BC side fields and properties |
| 1. Momentum | Stored m and c; p gradient; wall/side state | A_raw,b_raw -> relaxed A,b; d and optionally d_tilde; solve delta_u; u_star=u+delta_u |
| 2. Pressure pre-solve | Current volume p, previous reversal flags | Refresh outlet side trace and nodal-side trace; freeze them for this pressure solve and ensuing flux update |
| 3. Pressure solve | u_star,current p, old reconstructed g, current d, trace | Fresh sample flux and R_c; assemble H; solve H*phi=-R_c; p_new=p+alpha_p*phi |
| 4. Correction gradient | Raw, unrelaxed phi | Publish phi; construct G(phi) with eta=1; keep separate from pressure-gradient history |
| 5. Stored flux | u_star,p_new, old reconstructed g,d, frozen trace,m_old | m_new=alpha_m*m_hat(u_star,p_new,g_old,d,trace)+(1-alpha_m)*m_old; apply the actual reference boundary-specific updates |
| 6. Boundary status / divergence | Newly stored m and pre-correction u,p | Update reversal flags and zero blocked flux; form stored nodal c from the resulting flux |
| 7. Velocity correction | u_star,G(phi),d (d_tilde for SIMPLEC) | u_new=u_star-d*G(phi), without multiplying by alpha_p or alpha_m; no imported Chorin boundary mask |
| 8. Derived state | p_new,u_new | Refresh g_p, field scales, grad u, then velocity blending; post-solve wall coefficients and reports |

The specified-velocity inlet's stored flux also uses the relaxation formula,
with `m_hat=rho*u_bc dot A`, while its pressure RHS uses the fresh prescribed
source. The element field passes its URF to the side field when allocating it;
export both effective URFs. With constant BCs and a consistently initialized
inlet flux, this relaxation leaves the value unchanged. A ramp/history gate
must expose the difference. Stationary wall flux is initialized to zero and
the boundary update dispatcher leaves it unchanged. No second flux
reconstruction is inserted after step 7.

Pressure, corrected nodal velocity and stored mass flux are therefore different
iteration states. Demanding one-pass exact continuity of raw corrected velocity
would change the algorithm. For frozen trace/gradient/coefficients, the SIMPLE
interior and active-outlet fresh flux change under `p+=alpha_p*phi` is its
compact pressure sensitivity applied to that increment. This identity is a
useful gate; it does not include reversal, trace refresh or gradient refresh.

Every cache has a revision: topology/map, properties, time coefficients,
advection mass flux, momentum matrix/relaxation, d, boundary flags and pressure
matrix. SIMPLE momentum and pressure coefficients can change on **every outer
iteration**. Hypre hierarchy reuse across changed values is not authorized by
the existing Chorin physical-step cache. Publish changed node values before
any peer uses them. Preserve mass-flux iteration history across restart;
reinitializing it from raw velocity changes the resumed iteration.

## 7. Residuals and acceptance

The pinned LibLinSolve calls `CRSResidual` (r=b-Ax) before a solve and computes,
for component j,

\[
e_j={1\over s_j+SMALL}
\sqrt{{1\over N_{owned,global}}\sum_{i\ owned}(r_{ij}/a_{ij,ij})^2}.
\]

The increment starts at zero, so its initial r is the assembled RHS. The
diagonal is the **final solved matrix diagonal**, including relaxation and
boundary treatment. The source asserts positive scalar diagonals. There is
no volume weighting in this RMS. `linearSystem::convergenceReport_` uses the
scaled **initial** residual, not the final Krylov residual.

Velocity scale is the global maximum speed, with the generic near-zero fallback
to 1 when its range is below 1e-9. Pressure scale is overridden by
`flowModel::updatePressureScale` to `max(p_max-p_min,0.5*rho_scale*U_scale^2)`.
The linear solve reports before the field update. Momentum then updates its
velocity scale; pressure correction temporarily calls generic `p.updateScale()`;
the coupled loop restores the dynamic-pressure-based scale after velocity
correction. Export the actual scale at each report rather than recomputing it
from the final state. For constant positive density
in our profiles, rho_scale is rho. `SMALL` is scalar machine epsilon in the
pinned reference, not a user-selected physical regularization.

The reference's `is_converged_` can also become true on its iteration limit.
MARS must report `residual_converged` and `stopped_at_limit` separately and
must not call the latter successful convergence. Record linear final residual,
initial normalized outer residual and update norms as separate quantities.
The PI's requested 1e-6 has not been identified with this source criterion yet.

Additionally report, over all owned continuity rows, both fresh-RHS and stored
flux residuals and their boundary sums. With c in kg/s, the full normalized
divergence RMS is `sqrt(sum(c_i/rho)^2/V_i / sum(V_i))` in 1/s. Verify
`sum_owned c_i = sum_unique_boundary m_b` for fixed density/mesh. Raw velocity
divergence remains a separate diagnostic. A steady solution needs converged
momentum and pressure outer residuals and mass balance; a stable scalar speed
RMS alone is insufficient. Energy/pressure-work balance and spatial refinement
remain physical validation beyond the algebraic parity gate.

## 8. MARS implementation boundaries and typed reuse

All names under `fem/segregated/` below are **proposed new interfaces**. No such
production algorithm or executable is claimed to exist by this document.
Keep them outside `mars_ns_pump_solver.hpp`.

| Proposed operation | Contract/source | Reuse classification | First gate |
|---|---|---|---|
| `build_segregated_geometry` | Section 2; TetSCS/Tri3DSCS | Adapt existing topology; new sample tables and canonical IDs | Per-sample reference geometry, orientation and dual-volume comparison |
| `reconstruct_gradient` | Section 3; `nodeField::updateGradientField` | New reference-compatible assembly; reuse halos | Constant field and full nodal action parity, including boundary rows |
| `assemble_momentum` | Section 4; NS element/node/BC assemblers | New 3x3-block values; reuse device CSR storage | Every local block/RHS and assembled action before/after relaxation |
| `build_influence_coefficients` | `computeDUCoefficients` | New componentwise coefficient kernel | d,d_tilde and exact diagonal stage |
| `update_boundary_state` | Sections 5–6; flowModel and side fields | Adapt ownership/scratch, new reference policy | Trace, nodal-side interpolation, corners and reversal transitions |
| `evaluate_mass_flux` | Sections 5–6; pressure assembly and flowModel | New evaluator with fresh/stored variants | Every sample, RHS and stored history |
| `assemble_pressure_correction` | Compact H and frozen trace | New values, reuse CSR/Hypre | Finite-difference compact action, rows/columns, gauge checks |
| `advance_segregated_iteration` | Section 6 | New controller | Ordered state dump and single-iteration parity |

Inspected existing signatures permit the following reuse (these are interface
sketches, not compiled new call sites):

```cpp
using Solver = mars::fem::HypreGMRESSolver<double, int, cstone::GpuTag>;
using Matrix = Solver::Matrix;
using Vector = Solver::Vector;
Matrix H;
// Allocate once; fill sorted CSR on device using the new reference assembly.
H.allocate(n_owned, n_local_columns, nnz);
Vector rhs, phi;
thrust::device_vector<HYPRE_BigInt> d_local_to_global;
Solver pressure_solver(comm, max_iter, tolerance);
// Maps and vectors must be sized/populated before this call.
bool solved = pressure_solver.solve(H, rhs, phi, global_row_begin,
    global_row_end, 0, global_node_count, d_local_to_global);
```

`SparseMatrix<int,double,cstone::GpuTag>` exposes `rowOffsets()`, `colIndices()`,
`values()` and their `*Ptr()` accessors. Do not call its current `sortColumns()`
in the production loop: it sorts through host arrays. The new block momentum
can be expanded into this scalar CSR using component row IDs `3*g_node+j`.
That expansion/map is new work and must preserve all cross-component blocks;
the scalar solver API is reusable, a block-aware AMG configuration is not yet
validated. Do not assume a one-cycle AMG action is an adequate complete solve.

The prepared wrapper exposes `enable_reuse(bool amg_cycle=false)` and
`invalidate_setup()`. Its pointer checks cannot detect in-place value changes;
the new owner must invalidate collectively on a matrix/map/configuration
revision. Reuse within a frozen linear solve is permitted. Hierarchy lagging
across outer iterations needs a later, measured validation package.

`ElementDomain` supplies these existing vector-templated member calls:

```cpp
domain.reverseExchangeNodeHaloAdd(d_scalar_sum); // sum ghost contributions to owners
domain.exchangeNodeHalo(d_scalar_sum);          // publish completed owner values
domain.exchangeNodeHaloBlock(d_velocity, 3);    // node-interleaved overwrite halo
```

The scalar vector element type must equal the domain RealType. Reverse-add
leaves ghost slots untouched. Use zeroed scratch and unique owned-element/facet
scatters, then reverse-add and publish for complete nodal sums. Do not reverse-add
again after publishing or mix this route with full-star owned-row gathering.
For matrix rows, gather a complete incident-element star for each owned row;
do not also add halo rows. Assert missing connectivity and missing CSR entries.
All ranks enter reductions even when they own no relevant boundary samples.
No full nodal/CSR arrays return to the host in the production iteration.

`mars::Tet4CVFEM::jacobian_and_dNdx<RealType>(coords,det,dNdx)` and
`get_scs_nodes(int,int&,int&)` are candidates for direct algebraic reuse after
reference parity. **Do not reuse `Tet4CVFEM::scs_coords`**: its current table
does not implement the reference unshifted sample locations; even its first
entry `(0.5,0.5,0)` does not correspond to declared pair `(0,1)` under these
shape functions. This contract does not change that unrelated API.

## 9. Claude's next coding task: actual reference exports

**Recommended effort: Medium; one implementation agent.** The task is bounded
to the public oracle and comparator. If linking/extracting the reference needs
difficult architectural work, state the obstacle and use High for that named
step. No Ultra or additional agents are requested. GPT's subsequent focused
numerical review is High; documentation/commit work is Low.

Deliver this task before writing the production SIMPLE kernels:

1. Read this contract and its JSON manifest. Build a standalone reference
   validation target against the pinned OpenAccel code and pinned dependency,
   or instrument a separate checkout of that exact reference. Keep a precise
   patch and build instructions. Preserve source notices. Do not substitute
   independently retyped equations for the reference calls.
2. Create the native reference deck from the explicit profile. Export actual
   effective controls, boundary enumeration, gradient relaxation and histories.
   Use only the analytic fixtures in the manifest and public channel geometry.
3. Export actual TetSCS/Tri3DSCS geometry; momentum element/node/boundary
   matrices and RHS; relaxed matrix/RHS and d; pressure matrix/RHS; every fresh
   and stored flux; side/nodal-side traces and reversal flags; gradients;
   solved increments and the ordered states in section 6. Include raw and
   scaled residuals, scales, revisions, precision and output hashes.
4. Compare by global node ID, component, parent element and oriented sample
   ID, not array offset. Require local per-sample and matrix/action agreement
   before a single-iteration or final global balance comparison. Mark algebra
   fixtures that cannot represent a closed global solve as local-only.
5. Add the sensitivity and transition gates below. Export from real production
   reference assembly/evaluation; a Python model is a useful independent check,
   not the source of a claimed OpenAccel result.
6. Commit the harness, public fixtures and compact provenance manifest as a
   separate reviewable change. Report executed versus pending tests and exact
   commands. Large oracle dumps belong in reproducible artifacts, not private
   result directories. Ask GPT to review the diff and sample outputs before
   proceeding to the device-state stage.

Required falsifiable gates:

- On the unit tet `(0,0,0),(1,0,0),(0,1,0),(0,0,1)`, outlet `(1,2,3)` has
  `A=(1/2,1/2,1/2)`, `b_0=(-1,-1,-1)`. With rho=1 and constant d=2, the
  opposite-pressure derivative is +1 per boundary sample, +3 over the face.
  Opposite perturbation +0.2 changes mass flux by +0.6; common nodal-side
  trace perturbation +0.2 changes it by -0.6. No Chorin rho/dt multiplier.
- For constant u=(1,2,3) and matching affine compact/reconstructed pressure
  gradient `(2,-3,5)`, stabilization is zero and samples are `(1,1,1)`, total 3.
  This supplies reconstructed gradients directly; it does not assert an
  affine-exact reconstructed gradient at an arbitrary truncated boundary node.
- Set compact p and u to zero, d=2, opposite reconstructed gradient zero,
  face reconstructed x-gradients `(0,3,9)`, other components zero. Samples
  must be `(0,1/2,3/2)`. A face-mean gradient gives `(2/3,2/3,2/3)` with the
  **same total**, so a face-total-only gate cannot distinguish the operators.
- Vary face d as `(2,4,8)` with the same gradient fixture. Check each sample's
  unshifted d `(32/9,79/18,109/18)` and flux `(0,79/72,109/24)`, then perturb
  the opposite d: no boundary coefficient may change.
  Distinguish the derivative holding d and reconstructed gradients fixed from
  the Jacobian of the whole coefficient/gradient/trace update.
- Check momentum cross-component blocks on a skew tet. Vary stored mass flux
  without changing u to expose both upwind and the nodal `-c*u` term. Test
  alpha_u only once and boundary RHS damping only once at a shared corner.
- Compare d_tilde to the diagonal-plus-same-component-off-diagonal sum. SIMPLEC
  pressure matrix changes while its fresh-flux evaluator still uses d.
- Exercise nonuniform pressure on a bent outlet, scalar-area trace aggregation,
  each reversal transition, all-outlet-blocked failure, and closed-domain gauge.
- Use two successive outer iterations to detect mass history reset, premature
  g refresh or a post-correction flux rewrite. Test first versus later gradient
  reconstruction when gradient relaxation is enabled.
- Before production exposure, compare 1/2/4-rank owned results and global
  balance, including empty-opening ranks, shared corners and a partition that
  places opposite/face nodes on different ranks. This is later GPU validation;
  no CUDA/MPI execution is claimed by Stage 0.

Start double-precision local algebra comparisons at
`|a-b| <= 1e-12*max(1,max|reference block|)` in their stated SI units for the
well-scaled fixtures. Re-evaluate this tolerance for mesh/property scaling and
conditioning; never loosen it to hide a different stencil. For finite differences
use multiple perturbation sizes and show the roundoff/truncation trend. Linear
solutions use achieved residual and conditioning, not the matrix-entry tolerance
as a universal solution-error bound. Distributed reductions allow ordering drift.

## 10. Source map and limits of this handoff

Paths below are relative to the pinned OpenAccel repository unless marked MARS.
The JSON manifest hashes the inspected files so line references are tied to a
specific version. Numerical claims above are source inspection and derivation;
the harness must still demonstrate execution.

| Term/state | Source symbol / file |
|---|---|
| Outer order and full velocity correction | `segregatedFlowEquations::{solve,postSolve,preTimeStep}`, `src/equation/flow/segregatedFlowEquations.cpp` |
| Increment solve, BC refresh, zero start | `navierStokesEquation::{preSolve,solve}`, `src/equation/flow/navierStokes/navierStokesEquation.cpp`; `equation::correctField_`, `src/equation/equation.h` |
| Momentum transport and symmetric stress | `assembleElemTermsInterior_`, `src/assemble/flow/segregatedFlow/navierStokes/navierStokesAssemblerElemTerms.cpp` |
| Time and divergence terms | Steady/BDF node assemblers, `navierStokesAssemblerNodeTerms.cpp` in that directory |
| Relaxation and d | `phiAssembler::postAssemble`, `src/assemble/phiAssembler/phiAssembler.h`; `navierStokesAssembler::{computeDUCoefficients,postAssemble,assembleBoundaryRelaxation_}` |
| Momentum inlet/wall/outlet | `navierStokesAssemblerElemBoundaryConditions.cpp`, including `IP_FULL_STRESS_FIXED_VEL__` and `IP_ZERO_NORMAL_STRESS__` |
| Pressure matrix and fresh RHS | `src/assemble/flow/segregatedFlow/pressureCorrection/pressureCorrectionAssemblerElemTerms.cpp` and `pressureCorrectionAssemblerElemBoundaryConditions.cpp` |
| Pressure post-assembly/gauge | `pressureCorrectionAssembler::{postAssemble,adjustMatrixForPressureReference}` |
| p relaxation and raw phi | `pressureCorrectionEquation::{solve,correctField_}`, `src/equation/flow/pressureCorrection/pressureCorrectionEquation.{cpp,h}` |
| Stored flux and boundary flux | `flowModel::{updateMassFlowRate,updateMassFlowRateBoundaryFieldOutletSpecifiedPressure_,updateMassDivergenceField_}`, `src/model/flow/flowModel.cpp` |
| Trace, reversal, scales and wall coefficient | `flowModel` methods in sections 5 and 7; same file |
| Normal inlet and side transfer | `src/field/nodeField/vector/transport/flow/velocity.cpp`; `src/field/nodeSideField/nodeSideField.cpp`; `src/field/sideField/sideField.cpp` |
| Gradients and scales | `src/field/nodeField/nodeField.hpp`; `src/realm/fieldBroker.cpp` |
| Geometry and normal wall distance | Vendored Nalu master elements; `src/mesh/meshGeometry.cpp` |
| Residual and stopping | `src/equation/linearSystem.h`; pinned LibLinSolve `residual.h`, `linearSolverContext.h`, `matrix/operators/distributed/CRSResidual.h` |
| Existing GPU matrix/solver API | MARS `backend/distributed/unstructured/fem/mars_sparse_matrix.hpp`, `solvers/mars_hypre_gmres_solver.hpp` |
| Existing nodal communication | MARS `backend/distributed/unstructured/domain.hpp` |

The first implementation checkpoint is a working **reference export and
comparison harness**, not an exposed MARS solver. Remaining runtime questions
are deliberately assigned to it: native-deck interpretation, shared nodal-side
write ordering, reference BDF startup, achieved linear-solve accuracy and
distributed assembly equivalence. Findings that change this contract require
a short attributed correction before dependent numerical kernels are written.
