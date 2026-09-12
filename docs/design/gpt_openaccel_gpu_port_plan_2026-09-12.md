# OpenAccel numerical method on the MARS GPU backend

Author: GPT/Codex. Date: 2026-09-12.
Status: architecture and implementation plan requested by the user.
Implementer: Claude, with bounded GPT contributions when useful.
No implementation agent was started by this planning task.

## Decision

Reproduce the complete selected OpenAccel flow algorithm using MARS mesh,
device storage, assembly, communication and solvers. Numerical parity is the
first target; measured GPU improvements are the second. Do not develop another
boundary law to compensate for a partial port.

The current Chorin/BDF2 solver remains a regression and performance baseline.
This plan supersedes the earlier recommendation to keep extending its outlet
experiment. SIMPLE is not yet implemented in production MARS. "OpenAccel
parity" will mean a named reference revision, configuration, fixture and set
of executed comparisons, not a solver name or a generally successful CFD run.

The user's PI uses OpenAccel with an inward-normal velocity inlet, pressure
outlet and other sides no-slip. The clarified fluid has rho=1000 kg/m^3 and
mu=1e-3 Pa s, hence MARS nu=1e-6 m^2/s. Earlier nu=1e-4 runs were an intentional
MARS workaround. Exact OpenAccel algorithm/options and the convergence norm
behind the PI's 1e-6 criterion still need to be identified. Do not infer SIMPLE
versus SIMPLEC, an advection scheme, or a pressure-boundary subtype from that
physical description alone.

The current baseline has completed a user-run one-rank public water-viscosity
test: 200 steps at 38.5 ms/step, final full continuity RMS 5.36e-13 /s. This is
a short startup result, not steady-state convergence or OpenAccel parity.

## Responsibilities and effort

Effort labels are recommendations, not changes to either application's model
setting. Default to one active implementation or review agent per work package.
Claude and GPT normally work sequentially at handoff checkpoints. Parallel work
is useful only for disjoint files/tasks; no delegation is authorized by this
document alone. Ultra is not the default for any stage.

| Stage | Owner | Recommended effort and reason | Exit evidence |
|---|---|---|---|
| 0. Pin reference and write numerical contract | GPT architect | High: signs, boundary elimination and iteration ordering interact | Reviewed source-to-equation map and frozen public fixture/configuration |
| 1. Build reference comparison harness | Claude; GPT checks harness contract | Medium once outputs are specified; High only if extracting the reference needs difficult dependency work | Executed reference outputs with provenance and self-checks |
| 2. Device state, topology and communication | Claude | Medium: reuse established MARS patterns against a settled state contract | GPU ownership/halo/geometry gates; no solver exposed prematurely |
| 3. Momentum, influence coefficients and flux kernels | Claude; GPT reviews mathematical diff | High: discretization, relaxation and boundary terms must agree | Frozen-state matrix/action/RHS/coefficient/flux parity |
| 4. Pressure correction and one complete segregated iteration | Claude; GPT reviews the coupled update | High: pressure, velocity, flux and boundary states interact | One-iteration parity including all updates and residuals |
| 5. Outer iteration, physical-time history and MPI integration | Claude; GPT reviews lifecycle and acceptance | High for integration; Medium for running settled gates | Converged public cases at water properties and 1/2/4-rank agreement |
| 6. Performance optimization | Claude by measured hotspot; GPT for mathematical changes | Medium for storage/batching; High for communication or solver changes | Before/after timings and unchanged acceptance results |
| 7. Private application validation | User/PI on their systems | Their workflow; GPT only handles permitted public work | Owner-reviewed BCs, flow, convergence and reference comparison |

For each package, come back down to Low for status, documentation, commit and
handoff. Escalate to Max only for a named unresolved derivation that actually
needs deeper review. An Ultra step needs separate justification and explicit
approval. Do not launch multiple reviewers by default.

## Stage 0: the architectural contract

The source-derived [version 1 numerical contract](gpt_openaccel_numerical_contract_2026-09-12.md)
now defines the provisional profile, equations, update order, typed reuse and
Claude's first reference-harness ticket. Read [SYNC.md](../../SYNC.md) before
acting on a handoff; it records newer coordination updates. Reference execution
and production GPU parity remain subsequent gates.

Start from the previously inspected public OpenAccel revision
`0d69041ba1afda63e9e4328d9e0d9834bba37756` as a provisional navigation reference,
not as an assertion that the PI ran that revision. Record the PI's version and
nonconfidential algorithm settings when available. If unavailable, select and
label a provisional public configuration; do not block all infrastructure work
or claim parity with the PI's exact run.

The contract must identify:

1. PDE, physical versus kinematic pressure, viscosity convention and units.
2. Unknown locations, control volumes, shape functions and quadrature. Retain
   the reference's sampling and interpolation, including boundary traces.
3. Momentum time/advection/diffusion terms, source terms, linearization and
   under-relaxation, including the exact matrix stage used for influence coefficients.
4. Stored transport flux, newly reconstructed flux, correction flux and history;
   distinguish mass from volume flux and componentwise influence coefficients.
5. Pressure equation, gauge and boundary elimination. Derive the implemented
   SIMPLE/SIMPLEC approximation; do not replace it with our exact Chorin Jacobian.
6. Ordered updates of pressure, velocity, mass flux, reconstructed gradients,
   boundary values and nonlinear residuals. State which quantities remain frozen.
7. Outer-iteration acceptance, linear tolerances, residual normalization and
   physical-time history. Separate steady momentum convergence from continuity
   and changes in a scalar RMS speed.

Required output: a table mapping every term to an exact reference source symbol,
its frozen/current state, units, proposed MARS typed interface and its parity gate.
Check exact signatures before calling an existing MARS component reusable.

The prior review points to these public source families; reread the selected
revision rather than treating old line numbers as the contract:

- `src/equation/flow/segregatedFlowEquations.cpp`: iteration/update order.
- `src/assemble/flow/segregatedFlow/navierStokes/`: momentum and coefficients.
- `src/assemble/flow/segregatedFlow/pressureCorrection/`: pressure and BC terms.
- `src/model/flow/flowModel.cpp`: fluxes, gradients and boundary updates.

Stage 0 must also settle opening/wall corner precedence and inlet normals. MARS
currently uses area-weighted per-node normals and lets opening nodes override
walls at shared edges. The source/units findings are in
[the BC contract](gpt_flow_bc_contract_2026-09-12.md).

## Implementation boundaries and reuse

Create a separate opt-in segregated solver and a public-fixture driver. Proposed
names are `fem/segregated/` and `mars_segregated_flow`; these are architecture
placeholders, not existing targets or runnable commands. Avoid adding another
large algorithm inside `mars_ns_pump_solver.hpp`. Shared fixes should be small,
independently justified commits that preserve the current baseline.

| MARS asset | Reuse rule |
|---|---|
| Device vectors, mesh connectivity, DOF maps, halo topology | Reuse typed interfaces; verify the new fields' ownership and exchange semantics |
| Sparse assembly, reductions and facet ownership | Reuse infrastructure; derive values and quadrature from the reference |
| Opening geometry and flux helpers | Reuse only after sampling/normal/affine and frozen-state parity checks |
| FGMRES and GPU Hypre/AMG | Reuse where matrix properties and achieved residuals permit; record linear-solver differences from the reference |
| Prepared Hypre setup | Reuse only while matrix values, configuration and distributed mapping remain valid |
| Public gates and profiling | Reuse test mechanics; define acceptance for the new method rather than copying Chorin criteria |
| Chorin correction, experimental trace coupling and mass-only response | Not the numerical definition of the new solver |

SIMPLE's momentum and influence coefficients can change on every outer iteration.
The current prepared solver's physical-step lifetime is therefore not a valid
cache-invalidation policy for the new algorithm. Use explicit matrix/coefficient
revisions and mapping/configuration changes; rebuild collectively when required.
Hierarchy lagging across changed matrices is a later, separately tested optimization.

Keep fields, geometry, assembly, flux updates and numerical vector operations
on the device. CPU orchestration and API-required scalar decisions are allowed;
do not promise literally zero host activity. Use persistent buffers, unique
facet ownership, reverse-add then publish where appropriate, and peer exchanges.
Overlap communication only when data dependencies and buffer lifetimes permit.
No host field round-trips in the production iteration. Host oracle calculations
and small fixture exports belong to explicitly separate validation paths.

## Reference harness and acceptance ladder

Use the same public mesh, boundary tags, fields and input properties in both
implementations. Canonical geometric keys identify samples and DOFs across
different orderings/partitions. Record source revisions, full configuration,
precision and output hashes. First make the reference execute; a re-derived
NumPy replica alone is not evidence of OpenAccel implementation parity. If a
full reference build is unavailable, a narrowly extracted actual reference
routine can qualify a local term, with that limitation explicit.

Gate in this order:

1. Geometry/quadrature, constant/affine fields, nonuniform coefficients and
   boundary masks, including bent normals and opening/wall intersections.
2. Frozen momentum action/RHS, relaxed diagonal and influence coefficients.
3. Interior and boundary flux samples, assembled pressure action/RHS and gauge.
   Check finite-difference sensitivities against the *declared frozen-state*
   approximation; SIMPLE need not equal the full nonlinear Jacobian.
4. One complete outer iteration: compare pressure, velocity, stored transport
   flux, reconstructed gradients, boundary updates and reported residuals.
5. Multiple outer iterations and, when selected, BDF startup/history. The
   reference's relaxed intermediate flux need not be divergence-free after one
   subiteration. Require the reference convergence contract at outer acceptance.
6. Public channel at the target water properties, then a public contraction/
   expansion with a downstream chamber to exercise changing area and backflow.
   Match transient trajectories or steady residuals according to the chosen mode.
7. 1/2/4-rank invariants, including an empty-opening rank, before scaling claims.

Compare relative/absolute operator and field errors with declared scales and
linear-solve tolerances; do not demand bitwise equality across reduction orders.
Use the same accepted flux for local cancellation, full continuity and global
boundary balance. Do not substitute visual similarity, a small linear residual,
or a flat scalar speed history for these gates.

## Review workflow and usage control

GPT owns the contract and the final interpretation of mathematical evidence.
Claude owns implementation and its first validation pass. GPT may write a
bounded shared helper, oracle or surgical repair when that is more efficient;
declare file ownership first. Do not have both agents rewrite the same component.

Each Claude handoff should contain only:

```
Base commit / new commit(s):
Contract section implemented:
Files and numerical terms changed:
Tests actually run, with command and result:
Known failures or unverified GPU/MPI scope:
Specific review question, if any:
```

GPT reviews the scoped diff and its necessary call sites, not the entire repo
again. Low-risk scaffolding receives a Low/Medium review. Mathematical kernels,
BCs, relaxation, gauges, ownership and iteration history receive a High review
before integration. "Quick review" means a bounded scope, not skipping the
mathematics. Stop re-running unchanged checks once they pass unless new evidence
or a dependent change warrants another run.

Three mandatory architectural checkpoints are: the completed contract, first
complete reference-matched iteration, and distributed water-property validation.
Use concise committed handoffs rather than forwarding both agents' full histories
and raw logs repeatedly. Pull permitted public logs once and retain local summaries.
This avoids duplicate work; it does not promise a quota saving or metered cost.

The current user preference is to keep GPT/Codex and Claude memory separate.
This repository document is the shared plan. Neither agent should edit the other's
memory. User/PI private mesh and field inspection remains outside GPT's access.

## First Claude task

The Stage 0 contract is linked above. Claude should take its section 9 ticket:
build the public reference export/comparison harness before production kernels.

After the architect closes the contract's numerical questions, Claude's first
coding ticket is Stage 1: build the public
reference harness and produce the specified intermediate outputs. This prevents
a large GPU port from accumulating before numerical equivalence is measurable.
Subsequent implementation tickets should cover a coherent operator or update,
with concrete gates and a small reviewable diff. Commit and push completed
authorized work; the user performs cluster builds and runs.
