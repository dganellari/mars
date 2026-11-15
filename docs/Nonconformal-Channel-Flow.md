# Stationary non-conformal channel flow in MARS

**Status: experimental implementation with host, CUDA and HIP channel validation, 21 October 2026.**

This note explains the stationary Tet4 interface implementation in the local
`feature/nonconformal-scalar` branch. It describes the experimental source;
it is not a description of a released MARS capability.
The scalar and momentum/continuity examples have passed host and CUDA tests,
including the PI's mesh on 1, 2 and 4 MPI ranks. The new flow example also passes
the analytic plug check on one GPU. The later partitioned-input and device block-graph
changes also pass user-run CUDA and HIP field comparisons on 1/2/4 ranks.
The latest GPU comparison reports were supplied by the user; their raw fields and
build manifests have not been independently retrieved for this update.

## 1. What we are trying to achieve

The PI wants to mesh two parts of a fluid domain separately and solve one flow
problem across their common interface. Our first case is a stationary channel
split into two tetrahedral blocks. The triangular faces on one side do not match
the triangular faces on the other side.

MARS must connect the **equations** across that interface without requiring the
two meshes to share nodes. Fluid leaving one block must enter the other, and
momentum must be transferred consistently with that mass flow.

The work has progressed through three levels:

1. **Interface components:** geometric pairing, communication, sparse-matrix
   connectivity and conservative scalar fluxes.
2. **A complete scalar problem:** diffusion through the PI's nonmatching Tet4
   channel, checked against an analytic solution.
3. **The current flow example:** three velocity components, pressure, momentum
   transport and a pressure-correction equation enforcing discrete mass balance.

The third level now has host and CUDA execution checks. The first two alone
would not give the PI a channel-flow solver; the new checks still do not establish
spatial convergence or general flow-driver support.

## 2. What makes the mesh non-conformal?

A conforming connection has compatible faces, edges and vertices. In a usual
continuous discretization, adjacent elements refer to the same unknown at a
shared node.

At our non-conformal interface, the two blocks occupy adjacent regions and their
boundary surfaces coincide geometrically, but their surface triangulations and
node identities are independent. One triangle can overlap several triangles on
the other side. Even coincident corner coordinates retain distinct node IDs.

```text
inlet                 internal interface                   outlet
x = 0                       x = 1                           x = 2
  |                           |                               |
  |       Tet4 block A        |        Tet4 block B            |
  |           ---- fluid flow ---->                           |
  |                           |                               |
                   triangulation A != triangulation B
                   separate nodes and unknowns on each side
```

This interface is an artificial mesh boundary inside one fluid. It is not a wall
or a material membrane. It must transmit the flow rather than impose no-slip or
zero flux.

The two values of a field evaluated at the same interface location are called
its **traces**: for example, velocity `u_A` interpolated from block A and `u_B`
interpolated from block B. The numerical coupling controls their difference
through fluxes and penalty terms. It does not merge the unknowns or enforce exact
node-by-node equality.

### Why this is useful

Independent meshes let different parts use different resolutions and meshing
strategies. They avoid forcing a change on one side to propagate through the
entire volume mesh. They also provide an interface-coupling foundation for future
relative mesh motion.

For a conforming transition between hexahedra and tetrahedra, pyramids can connect
quadrilateral and triangular faces. A non-conformal connection can avoid that
particular conforming transition requirement, but each side still needs supported
element operators. **Our current implementation supports Tet4/Tet4 interfaces;
it does not yet make Hex8/Tet4 interfaces or pyramids usable.**

Ghosting is a separate concern: it supplies remote data when the two sides belong
to different MPI ranks. It is needed even for a stationary interface. Moving
meshes would additionally require updating the pairing and communication when
the overlap changes, plus the appropriate moving-grid or rotating-frame equations.
None of that mesh-motion physics is implemented by this channel example.

## 3. The mesh and physical problem used here

The PI's `non_conformal_tet_tet.exo` describes a channel with
`0 <= x <= 2`, `0 <= y,z <= 1`, split at `x = 1`. Each block is internally
conforming; the two blocks have separate node identities at the interface.

The example solves steady incompressible flow with constant density and dynamic
viscosity, no body force and no turbulence model:

$$
\nabla\cdot(\rho\mathbf u)=0,
\qquad
\nabla\cdot(\rho\mathbf u\otimes\mathbf u)
+\nabla p-\nabla\cdot\left[\mu(\nabla\mathbf u+\nabla\mathbf u^T)\right]=0.
$$

| Quantity | Meaning | Units |
|---|---|---|
| `u` | Velocity, three components per node | m/s |
| `p` | Physical pressure, not pressure divided by density | Pa |
| `rho` | Density | kg/m³ |
| `mu` | Dynamic viscosity, not kinematic viscosity | Pa·s |
| `m_dot` | Integrated mass flux through a surface patch | kg/s |

Both velocity and pressure use nodal linear Tet4 interpolation. The conservation
equations use **CVFEM: control-volume finite elements**. Each tetrahedron is split
into subvolumes associated with its nodes; the node's control volume is the union
of its incident subvolumes. Fluxes are assembled through their surfaces.

The default channel settings are:

- Inlet: `u_x = 16 U y(1-y) z(1-z)`, `u_y = u_z = 0`, with `U = 0.1`.
- Walls: stationary, no-slip velocity.
- Outlet: prescribed pressure trace `p = 0`, with the existing open-boundary
  momentum treatment. This pressure condition fixes the pressure gauge.
- `rho = 1`, `mu = 0.1`. Taking the channel height `L = 1` and the inlet peak
  speed `U` as reference scales gives `Re = rho U L / mu = 1`.

The polynomial inlet is a prescribed profile, not a claim that we impose the exact
fully developed square-duct solution. Boundary fluxes are evaluated by the discrete
face quadrature, so their numerical integral need not equal the exact integral of
the continuous polynomial.

There is also an analytic **plug-flow check**. It prescribes uniform axial velocity
on the inlet **and on the side walls**, so `u = (U,0,0)` and `p = 0` are compatible
exact fields. It tests the coupling; it is not the stationary no-slip-wall case.

Source: [flow driver](https://github.com/dganellari/mars/blob/v0.2.0/examples/distributed/unstructured/mars_nonconformal_flow.cpp)
and [flow runtime](https://github.com/dganellari/mars/blob/v0.2.0/examples/distributed/unstructured/mars_nonconformal_flow_runtime.hpp).

## 4. Turning two surface meshes into one conservative interface

### 4.1 Read and validate the two blocks

[The preparation script](https://github.com/dganellari/mars/blob/v0.2.0/scripts/prepare_nonconformal_mesh.py) reads Exodus
coordinates, Tet4 connectivity and the five named side sets: inlet, outlet, wall,
interface left and interface right. Their names can be mapped through the script's
options. It preserves the file's node and element indexing and does not weld
coincident nodes or convert unsupported element types.

Preparation checks finite, nondegenerate geometry, separate block node identities,
connected blocks and exact assignment of exterior faces to the supplied side sets.
Duplicate faces, internal faces incorrectly labelled as boundaries and unsupported
Hex8 or pyramid elements are rejected.

The `--flow-boundaries` option retains explicit inlet/outlet/wall face labels in
the prepared `MARSNC02` input. The flow example needs those labels because a
velocity condition and a pressure outlet cannot be reconstructed from the earlier
scalar boundary values alone. The older scalar input format remains supported
by the scalar example.

### 4.2 Search for candidate faces

Each interface face gets an axis-aligned bounding box. MARS's coarse search finds
potentially overlapping faces locally and on supplied candidate MPI peers.
Bounding-box overlap only identifies candidates; it does not establish a valid
fluid connection.

The interface plan then obtains opposing tetrahedron data and checks the actual
surface geometry. Faces must be coplanar within tolerance, have compatible
opposite outward normals, and cover each other without positive-area duplicates.
The coverage check tests overlaps as well as total area: a duplicated region must
not be allowed to compensate for a hole elsewhere.

There is a limit to this check: it validates the supplied surfaces. If both side
sets omit the same region, the search alone cannot discover that missing region.
The input's complete exterior-face labelling is a separate safeguard.

Source: [`NonconformalInterfacePlan`](https://github.com/dganellari/mars/blob/v0.2.0/backend/distributed/unstructured/fem/mars_nonconformal_search.hpp).

### 4.3 Build shared overlap patches

A triangular boundary face is divided into six small median-dual triangles: two
belong to each face vertex. For each candidate pair of faces, the implementation
intersects these small triangles in a face-local two-dimensional coordinate system.
Every nonempty intersection is a convex polygon with a physical area and centroid.

One resulting patch belongs to a control-volume row on A and a row on B. It stores:

- the two receiving row IDs;
- the eight nodal columns from the two adjacent tetrahedra;
- its area, common normal and two normal length scales;
- interpolation weights at the common centroid and the two Tet4 shape gradients.

The tetrahedron's off-face vertex has zero face interpolation weight but generally
a nonzero gradient coefficient. We therefore need the **whole opposing
tetrahedron**, not just its three face vertices.

The patch uses one centroid evaluation. This exactly integrates an affine
integrand over a planar polygon; it is not exact quadrature for arbitrary nonlinear
products in the flow flux. Conservation follows from using the same patch rule
on both sides, independently of that integration accuracy.

Source: [`InterfacePatch` and `nonconformalConservativePair`](https://github.com/dganellari/mars/blob/v0.2.0/backend/distributed/unstructured/fem/mars_nonconformal_overlap.hpp).

### 4.4 Insert equal and opposite contributions

Choose a canonical orientation with normal `n` pointing out of A. A patch flux
`F_q` contributes

$$
R^A_{i(q)}\mathrel{+}=F_q,
\qquad
R^B_{j(q)}\mathrel{-}=F_q.
$$

Thus the patch's contribution to the sum over both blocks is zero. This applies
to arbitrary nodal values, not only a constant or affine test field.

The two row owners can evaluate the patch on different ranks. They use the same
canonical face ordering, geometry and refreshed field values, then opposite signs.
There is no requirement to evaluate the floating-point flux physically once on
one rank; the requirement is the same discrete flux definition on both owners.

This distinction mattered during implementation. Integrating A and B separately
with unrelated surface quadrature can produce two different flux totals on a
nonmatching interface. Our current shared-overlap path avoids that mismatch.

## 5. From scalar diffusion to momentum and continuity

### 5.1 The scalar problem established the interface machinery

For scalar diffusion, `q = -k grad(phi)`. On a patch with area `S`, the flux is
an average of the two physical normal fluxes plus a jump penalty:

$$
F_\phi=S\left[
-\tfrac12(k_A\nabla\phi_A+k_B\nabla\phi_B)\cdot\mathbf n
+\lambda_k(\phi_A-\phi_B)\right],
\qquad
\lambda_k=\tfrac12\left(\frac{k_A}{h_A}+\frac{k_B}{h_B}\right).
$$

Here `h_A` and `h_B` are the tetrahedral altitudes normal to the interface faces.
The penalty discourages a discontinuity while the gradient terms transmit flux.
For constant fields the flux vanishes; for a common affine field with equal
diffusivities the jump vanishes and the physical diffusion flux is recovered.

The PI scalar test uses the analytic field `phi = x/2`. It established that mesh
reading, geometric overlap, cross-interface matrix entries, ghost exchange and
the linear solve work together before introducing pressure/velocity coupling.

### 5.2 The flow interface needs a mass flux

For collocated nodal velocity and pressure, simply averaging velocity at a face
does not provide the pressure coupling used by this discretization. The flow
example includes a pressure-gradient correction and a pressure-jump term.

The following notation describes the actual kernel:

- `bar(u)` and `bar(p)` are the arithmetic means of the two interpolated traces.
- `bar(g)` is the mean of the two element pressure gradients.
- `bar(G)` is the mean interpolated, reconstructed nodal pressure gradient.
- `D = diag(d_x,d_y,d_z)` is the mean interpolated momentum influence.
- At an unconstrained node, `d_i = V / a_ii`, using the relaxed momentum diagonal.
  Its units are m³·s/kg, so `D grad(p)` has velocity units.
- `H = (1/h_A + 1/h_B)/2` and `d_n = n^T D n`.

Then the patch's outward mass flux is

$$
\dot m_q=\rho S\left[
\overline{\mathbf u}\cdot\mathbf n
+\bigl(D(\overline{\mathbf G}-\overline{\mathbf g})\bigr)\cdot\mathbf n
+d_n H(p_A-p_B)\right].
$$

The gradient difference is a Rhie–Chow-type pressure stabilization. For a common
affine pressure represented consistently on both sides, that difference and the
pressure jump vanish. A pressure discontinuity drives a flux from high pressure
to low pressure with the chosen outward-normal convention.

The kernel also returns each coefficient `c_qk = d(m_dot_q)/d(p_k)`, holding the
reconstructed gradients and influence fixed. The pressure-correction assembly
and later flux correction use these same coefficients.

Source: [`nonconformalFlowFlux`](https://github.com/dganellari/mars/blob/v0.2.0/backend/distributed/unstructured/fem/mars_nonconformal_flow.hpp).

### 5.3 Momentum uses the same transported mass

The interface momentum flux combines upwind advection, the full symmetric viscous
stress and a velocity-jump penalty. In compact notation, its advective/viscous part is

$$
\mathbf F_{m,q}=
\dot m_q^{\mathrm{stored}}\mathbf u_{\mathrm{upwind}}
-\frac{S}{2}\left[
\mu_A(\nabla\mathbf u_A+\nabla\mathbf u_A^T)
+\mu_B(\nabla\mathbf u_B+\nabla\mathbf u_B^T)\right]\mathbf n
+S\lambda_\mu(I+\mathbf n\mathbf n^T)(\mathbf u_A-\mathbf u_B),
$$

where `lambda_mu = (mu_A/h_A + mu_B/h_B)/2`. Upwind selects A for nonnegative
stored mass flux and B for negative flux. The stored flux is the corrected flux
from the preceding pressure-correction step.

The full symmetric stress couples the three velocity components; it is not
implemented as three unrelated scalar diffusion equations. The jump term acts
on tangential and normal velocity differences.

The pressure traction is `S bar(p) n`. In the assembled momentum equation it is
represented through the reconstructed pressure gradient: the interface correction
replaces the one-sided trace by the common trace. It must not also be added as a
second explicit pressure force. The kernel's separate `pressureForce` output is
used to check the combined interface momentum balance.

Both the mass flux and momentum flux receive opposite row signs across the
interface. Copying velocity values between meshes alone would not give this
conservative coupling.

## 6. The channel iteration

The new runner uses a **SIMPLE-style iteration**: Semi-Implicit Method for
Pressure-Linked Equations. It separates momentum and pressure correction while
coupling them through the mass flux and momentum influence.

For each outer iteration it:

1. Refreshes ghost states and reconstructs the pressure gradient, including the
   common interface trace and prescribed outlet pressure trace.
2. Assembles the three-component momentum residual and matrix, using the stored
   conservative mass flux for advection. It solves for a velocity increment and
   forms the predictor velocity.
3. Exchanges the predictor and momentum influence, then assembles the continuity
   residual and its pressure-correction matrix.
4. Solves `B p' = -R_c`. The interface part of `B` consists of the coefficients
   returned by the same mass-flux kernel.
5. Corrects and stores each flux using the **full** pressure increment:
   `m_dot_new = m_dot_predictor + sum_k(c_qk p'_k)`.
6. Updates pressure with relaxation, `p <- p + 0.3 p'`, and unconstrained velocity
   with `u <- u_predictor - D G(p')`.
7. Checks outer residuals, and on convergence checks mass balance, interface
   momentum balance, finite fields and any requested analytic comparison.

The momentum diagonal uses relaxation `alpha_u = 0.7` and an artificial
pseudo-time contribution with `delta_tau = 1`. This stabilizes a steady iteration;
these iterations are not physical time steps or a time-accuracy test.

Pressure relaxation does not scale down the stored mass-flux correction. Applying
the `0.3` factor there would retain part of the continuity defect. Recomputing an
unrelated velocity-only flux afterwards would likewise lose the flux consistency.

Both linear systems use the existing MARS
[`outlet_fgmres`](https://github.com/dganellari/mars/blob/v0.2.0/backend/distributed/unstructured/fem/mars_outlet_fgmres.hpp)
controller with Jacobi preconditioning, restart 100, a 3000-iteration cap and a
relative residual target of `1e-10`. It recomputes `b-Ax` before accepting a solve.
The operators are not assumed symmetric positive definite, so the example does
not select conjugate gradient. It does not use Hypre or PETSc for these solves.

The outer target is `1e-7` relative to the initial momentum and continuity residual
norms, with small physical-scale absolute floors. A successful linear solve alone
does not establish outer convergence or spatial accuracy.

## 7. How this uses MARS

This is a new example assembled from existing MARS numerical and communication
components, plus the new interface operators. It is not an external solver that
only uses MARS to read a mesh.

| Component | Role in this implementation |
|---|---|
| [Coarse search](https://github.com/dganellari/mars/blob/v0.2.0/backend/distributed/unstructured/fem/mars_coarse_search.hpp) and [fine search](https://github.com/dganellari/mars/blob/v0.2.0/backend/distributed/unstructured/fem/mars_fine_search.hpp) | Candidate face discovery and triangular interpolation geometry |
| [Face extraction](https://github.com/dganellari/mars/blob/v0.2.0/backend/distributed/unstructured/fem/mars_nonconformal_mesh_faces.hpp) | Builds face records from Tet4 connectivity and side-set tags |
| [Interface plan](https://github.com/dganellari/mars/blob/v0.2.0/backend/distributed/unstructured/fem/mars_nonconformal_search.hpp) and [overlap geometry](https://github.com/dganellari/mars/blob/v0.2.0/backend/distributed/unstructured/fem/mars_nonconformal_overlap.hpp) | Coverage checks, common patches and receiving rows |
| [Ghost registry](https://github.com/dganellari/mars/blob/v0.2.0/backend/distributed/unstructured/fem/mars_ghost_registry.hpp) and [donor exchange](https://github.com/dganellari/mars/blob/v0.2.0/backend/distributed/unstructured/fem/mars_nonconformal_donor_exchange.hpp) | Ownership, remote element import and repeated field refresh |
| [Interface graph builder](https://github.com/dganellari/mars/blob/v0.2.0/backend/distributed/unstructured/fem/mars_nonconformal_scalar_graph.hpp) | Adds opposing-element columns to the owned sparse rows |
| [Tet4 bulk operator](https://github.com/dganellari/mars/blob/v0.2.0/backend/distributed/unstructured/fem/segregated/mars_segregated_tet_interior.hpp), [geometry](https://github.com/dganellari/mars/blob/v0.2.0/backend/distributed/unstructured/fem/segregated/mars_segregated_geometry.hpp) and [boundary operator](https://github.com/dganellari/mars/blob/v0.2.0/backend/distributed/unstructured/fem/segregated/mars_segregated_boundary.hpp) | Existing CVFEM volume and external-boundary algebra |
| [Flow flux](https://github.com/dganellari/mars/blob/v0.2.0/backend/distributed/unstructured/fem/mars_nonconformal_flow.hpp) | Shared mass flux, coupled momentum flux and their linear coefficients |
| [Shared example runtime](https://github.com/dganellari/mars/blob/v0.2.0/examples/distributed/unstructured/mars_nonconformal_scalar_runtime.hpp) and [flow runtime](https://github.com/dganellari/mars/blob/v0.2.0/examples/distributed/unstructured/mars_nonconformal_flow_runtime.hpp) | Device arrays, matrix actions, reductions, iteration and diagnostics |

The older [`mars_nonconformal.hpp`](https://github.com/dganellari/mars/blob/v0.2.0/backend/distributed/unstructured/fem/mars_nonconformal.hpp)
contains an unfinished integration sketch. It is not the implementation exercised
by this example. Similarly, the existing
[search-and-ghosting guide](https://github.com/dganellari/mars/blob/v0.2.0/docs/Interface-Search-And-Ghosting.md) explains enabling
components; their existence alone does not establish a complete flow coupling.

### Sparse algebra and MPI ownership

The pressure graph has one scalar row per owned node. The momentum graph expands
each node connection into a 3-by-3 velocity block. Interface rows include columns
from both adjacent tetrahedra, so the linear systems couple the two blocks rather
than solving each block independently with frozen interpolated boundary values.

Each row has one owner. Remote column values are ghost copies refreshed from that
owner; they are not additional independent unknowns. Global norms and balances
count owned rows only. The example assembles the receiving rows on their owners,
so it does not rely on adding duplicate ghost-row contributions back to owners.

The mesh is stationary. Pairing, graph connectivity and communication lists are
built once. Field values are exchanged repeatedly as velocity, pressure and
pressure correction change.

### What stays on the GPU, and what currently does not

The device path keeps numerical fields, matrix values, matrix-vector operations,
interface evaluation and packed halo payloads on the GPU. Ghost exchange uses
GPU-aware MPI with persistent buffers and stream completion before communication.
Norms and the small Krylov Hessenberg problem involve host scalars; CSV output
copies owned fields to the host explicitly.

The original `--mesh` mode replicates input staging. The new `--partitioned-mesh`
mode reads one file per rank containing only cells incident to owned rows, their
local coordinates, and boundary/interface tags. `partition_nonconformal_mesh.py`
prepares these files offline, preserves the original global node identities, and
includes off-face corners needed by gradient reconstruction. Coincident interface
nodes remain independent. Candidate peer lists combine volume dependencies and
conservative interface envelopes; runtime searches prune those candidates.

The partitioned path retains global row IDs in the algebra but uses local node
indices for geometry and boundary data. Donor columns map into persistent state
slots, and the momentum block graph is expanded on device without copying it to
the host. GhostRegistry still needs host identity metadata during setup.

This removes full-mesh replication from the solver processes; it does not establish
large-scale performance. Offline partition preparation still reads the whole mesh
and uses file-index row ranges rather than a geometry-aware partitioner. This is
not an adapter for `ElementDomain`/Cornerstone or automatic interface support in
other Navier–Stokes drivers.

## 8. What the checks establish

These are different questions and must remain separate:

| Check | What it establishes |
|---|---|
| Constant and affine fields | Basic consistency of interpolation, gradients and flux signs |
| Shared-overlap arbitrary-field checks | Equal and opposite scalar patch contributions |
| Pressure/velocity coefficient differences | Returned linear coefficients match the flow flux they linearize |
| Affine symmetric stress and rigid rotation | Correct component coupling and zero viscous stress for rigid rotation |
| Recomputed linear residual | The returned vector solves that assembled linear system to its requested tolerance |
| Outer residual history | The steady nonlinear iteration meets its defined stopping criteria |
| Boundary and interface balances | Conservation of the actual stored numerical flux |
| Owned-field comparison across rank counts | Partition-independent results within the observed rounding differences |
| Analytic solution or mesh refinement | Spatial accuracy; rank agreement alone cannot establish it |

The reported `continuity_l1` measures the sum of absolute owned control-volume
mass defects, normalized by inlet mass flow. It uses the **pressure-corrected
stored flux**. The separate `raw_velocity_divergence_l1/inlet_mass` diagnostic
uses element velocity gradients. A tiny first quantity does not imply that the
interpolated velocity is pointwise divergence-free.

The printed interface `momentum_balance` is the cancellation across the internal
interface, including pressure traction. It is not a whole-channel force-balance
measurement including wall reactions.

### Results available on 21 October 2026

| Test | Observed result | Scope |
|---|---|---|
| PI scalar mesh, CUDA on 1/2/4 ranks | 83 iterations; independent relative residual `7.944e-12`; maximum error in `phi=x/2` about `2.217e-10` | User-run GPU results, retrieved and checked |
| New host regression suite | 14/14 pass | Flux, scalar, procedural flow and MPI checks |
| New procedural flow sanitizer check | Pass on 4 MPI ranks | Host AddressSanitizer/UndefinedBehaviorSanitizer; not a GPU race check |
| PI no-slip flow, host on 1/2/4 ranks | 479 outer iterations on each; maximum field difference across rank counts below `1.12e-16` | Same mesh and controls; no mesh-convergence claim |
| PI no-slip flow mass balance | Boundary and interface imbalance below `3.23e-16` relative to inlet mass flow; corrected continuity L1 below `2.05e-15` | Numerical-flux balance |
| PI analytic plug flow, host | Velocity maximum error about `2.93e-9`; pressure maximum error about `5.01e-8` | Uniform-velocity boundary check |
| PI no-slip flow, CUDA on 1/2/4 ranks | 479 iterations on each; velocity/pressure differences across GPU counts and against the host below `1.12e-16` absolute | Retrieved owned fields checked by original node index |
| PI analytic plug flow, CUDA on 1 rank | Velocity error `2.930795e-9`; pressure error `5.005033e-8` | Matches the host check; 495 iterations |
| Partitioned no-slip flow, CUDA on 1/2/4 ranks | Maximum velocity difference `2.78e-17`, pressure `5.56e-17` against the earlier CUDA reference | User-supplied field comparison, all ranks pass |
| Partitioned no-slip flow, HIP on 1/2/4 ranks | Maximum velocity difference `2.78e-17`, pressure `1.67e-16` against the same reference | User-supplied field comparison, all ranks pass |
| Analytic plug flow, HIP on 1 rank | 495 iterations; velocity error `2.930795e-9`, pressure error `5.005033e-8` | User-supplied solver output; corrected continuity `2.469e-15` |

The no-slip case still has a maximum interpolated interface velocity jump of
about `2.394e-3`, roughly 2.4% of the chosen inlet peak speed. The coupling is weak,
so exact trace equality is not enforced. Conservation and rank agreement do not
make this jump disappear or establish its adequacy. The new procedural refinement
study below measures how the discrepancy behaves on a controlled mesh family.

The CUDA no-slip runs report `raw_velocity_divergence_l1/inlet_mass = 0.03835290`
on all three rank counts. This is a volume-weighted velocity-gradient diagnostic
normalized by inlet mass flow, not a 3.8% boundary mass leak: the corrected-flux
boundary balance is at roundoff. It is nevertheless a nonzero discrete
incompressibility measure. Together with the trace jump, it requires accuracy
assessment against refinement and a conforming reference. The new study supports
decreasing errors for the procedural channel; it does not assign an acceptable
error threshold to the PI's original mesh.

The near-roundoff conservation numbers reflect algebraic cancellation and accurate
linear solves. They are not estimates of the physical solution's accuracy.
The CUDA build and run logs, eight owned-field CSV files and 54-file kit manifest
were retrieved for local inspection. The manifest matches the prepared kit; the
CSV checks found finite fields, unique ownership and identical coordinates across
rank counts. The runs validate this kit, not subsequent source changes.

### Controlled refinement and partition checks

The host accuracy driver uses steady flow with density 1, dynamic viscosity 0.1,
inlet-profile scale 0.1 and the same wall/outlet conditions. Each half-channel has
`n` streamwise intervals. Shared-node and matching-interface controls have `n`
intervals in each transverse direction on both blocks. The nonmatching mesh has
`n` transverse intervals on the left and `2n` on the right. The refinement ratios
remain fixed instead of tending toward a matching interface.

The shared-node control has no interface flux. Matching-interface cases use the
same tetrahedra but separate interface node identities. They isolate the coupling
change. Nonmatching cases are compared at coincident nodes with a shared-node
reference at resolution `2n`; that comparison also includes bulk mesh error.
All ten runs converged with the outer tolerance tightened to `1e-9`.

| `n` | Matching versus shared max velocity difference | Nonmatching versus finer shared max velocity difference | Nonmatching max pressure difference | Nonmatching reported trace jump |
|---|---:|---:|---:|---:|
| 2 | 1.7643e-3 | 3.3493e-2 | 5.8468e-2 | 2.063e-2 |
| 4 | 3.0094e-4 | 8.4812e-3 | 1.9458e-2 | 1.136e-2 |
| 8 | 6.7790e-5 | 2.0392e-3 | 9.3446e-3 | 4.032e-3 |

Velocity and pressure differences decrease in both maximum and volume-integrated
norms. The latter integrates the P1 interpolant of nodal differences, not the
entire fine-grid error field. `check_nonconformal_flow_accuracy.py` independently
reconstructs element gradients and integrates absolute/squared raw divergence;
it does not substitute corrected flux continuity for velocity incompressibility.
The nonmatching raw divergence L2 falls from 6.8486e-3 at `n=4` to 4.3310e-3 at
`n=8`; the corresponding shared-node `n=8` value is 4.2574e-3. The coarsest-to-middle
raw-divergence norm is not monotone, so these results do not establish a uniform
order for that diagnostic.

After changing storage, the partitioned PI no-slip runs on 1/2/4 host MPI ranks
still converge in 479 iterations and match the previously validated CUDA fields
within 1.12e-16 absolute. A smaller four-rank partitioned plug channel passes
AddressSanitizer and UndefinedBehaviorSanitizer. The host suite remains 14/14,
and eight Python test groups cover input, partitioning and independent comparisons.
The user-run CUDA and HIP field comparisons above also pass after these storage
changes. The earlier refinement runs
used the replicated layout; no numerical flux or boundary formula changed in the
storage refactor.

Test sources: [scalar checks](https://github.com/dganellari/mars/blob/v0.2.0/tests/nonconformal/scalar.cpp),
[distributed checks](https://github.com/dganellari/mars/blob/v0.2.0/tests/nonconformal/distributed.cpp),
[scalar channel](https://github.com/dganellari/mars/blob/v0.2.0/tests/nonconformal/channel.cpp),
[flow flux checks](https://github.com/dganellari/mars/blob/v0.2.0/tests/nonconformal/flow.cpp) and
[analytic flow fixture](https://github.com/dganellari/mars/blob/v0.2.0/tests/nonconformal/flow_channel.cpp).

## 9. How the PI uses this work

The current entry point is `mars_nonconformal_flow`, with a prepared Tet4 mesh.
The workflow is:

1. Supply two disjoint Tet4 blocks with the required exterior and interface side
   sets. Prepare the mesh with `--flow-boundaries`; preparation validates and
   writes input but does not launch a solver. For distributed input, partition this
   prepared file for the intended rank count; each rank then reads its own partition.
2. Build the example for the chosen backend. The current independent source kit
   has passed CUDA checks without changing the frozen release.
3. Run the analytic plug check, then the stationary-wall channel. The driver
   accepts `--mesh` or `--partitioned-mesh`, plus `--output`, `--check-plug`, `--iterations`, `--inlet-speed`,
   `--viscosity` and `--density`.
4. Check linear and outer convergence, mass and momentum transfer, and the
   interface jump. Compare owned velocity/pressure fields across rank counts.
5. After convergence, inspect the per-rank CSV output. It contains original
   one-based file node indices, coordinates, three velocity components and
   physical pressure; existing output files are not overwritten.

The present executable is tailored to the PI's 2-by-1-by-1 channel and its named
boundary roles. General geometries and boundary models require further driver
integration. Actual build and cluster-launch commands are supplied separately;
this note documents the implementation rather than defining a second run recipe.

Before treating this as a generally available flow capability, we still need
broader accuracy studies and general mesh/domain integration. The partitioned channel
path passes the reported CUDA/HIP checks, while performance scaling remains unmeasured.
Curved or moving interfaces,
mixed Hex8/Tet4 coupling, pyramids and high-Reynolds-number cases remain outside
the validated scope. There is no general stability or energy estimate established
by the current tests, and no claim of OpenAccel field parity from this channel work.

## 10. Relation to Nalu/OpenAccel and the earlier interface discussion

The broad approach is consistent with surface-flux coupling described in the
[Nalu theory manual's non-conformal boundary section](https://nalu.readthedocs.io/en/latest/source/theory/boundaryConditions.html#non-conformal-boundary-condition):
opposing element data, averaged diffusive fluxes, jump penalties and implicit
cross-interface columns. That documentation describes a dual-pass formulation;
our current geometric integration instead uses shared overlap patches to obtain
opposite contributions on nonmatching faces. This is not a claim of identical
implementation or validated equivalence to the PI's STK/OpenAccel example.

The requirement to use pressure-stabilized mass flux consistently in advection
and continuity is also explained in
[Nalu's pressure-stabilization section](https://nalu.readthedocs.io/en/latest/source/theory/pressureStabilization.html#the-role-of-dot-m).
The equations and code mapping in this note are derived from the current MARS
source; the external references provide background, not validation of this code.
