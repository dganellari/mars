# Average-pressure outlet: implementation handoff

Author: GPT/Codex, with Codex implementation and review agents.
Date: 2026-09-09.
Status: implemented; host gates pass; CUDA build,
CUDA/MPI execution, domain-halo validation, and manufactured flow validation pending.

## Review follow-up: pressure-anchor guard

GPT/Codex, 2026-09-10. Status: implemented; CUDA/MPI validation pending.
Following Claude's missing-CSR-entry finding, outlet matrix assembly now flags any
missing opposite DOF, invalid owned-row mapping, missing CSR slot, or unresolved
facet caught by its existing connectivity guards. Each owned face row must have its
opposite-node column; non-owned rows are skipped. This avoids the incorrect assumption
that each local matrix facet requires three local writes.

A persistent device failure flag is reset on each refresh. One scalar is copied to
the host for `MPI_Allreduce(MAX)`; every rank participates, including empty-facet ranks.
A failure aborts before publishing `Apre` or entering Hypre. No geometry or patch
counts are printed. The numerical entry and default outlet path are unchanged.
This checks expected entries in the supplied facet list; it does not prove that the
list is complete, the resulting matrix is nonsingular, or the iteration contracts.

The CUDA kernel gate now includes 11 fault/ownership cases (22 assertions), including
a valid ghost column with only one owned face row, missing DOFs/CSR slots, and empty
ranks receiving the failure reduction. A serial host harness extracted the production
kernel, CSR lookup, and these fixture cases: all 22 assertions passed with
`-Wall -Wextra -Werror` and address/undefined-behavior sanitizers. Atomics and MPI were
serial substitutes in that harness; CUDA compilation and real 1/2/4-rank execution
remain pending because this host has no `nvcc`.

## Provenance

The takeover started at `c3cfcc53575c95f84f687429db7cf8036676bcd0` on `cstone`.
Claude had committed D0 (frozen context), D1 (boundary gradient), and D4 (free outlet
pressure rows), and left an uncommitted 117-line facet evaluator/scatter scaffold.
GPT extended that scaffold and implemented D2, D3, the correction loop, and gates.
Claude's earlier implementation remains attributed to Claude. Unrelated working-tree
changes were preserved. No confidential geometry, outputs, or case documents were read.

The mathematical contract remains
[the GPT boundary specification](gpt_outlet_boundary_spec_2026-09-09.md).
[Claude's status](outlet_trace_status.md) records the earlier incomplete state.

## What is now connected

- `mars_outlet_flux.hpp` contains the production facet arithmetic, tetrahedral
  gradient, face-only coefficient, and pressure partial. It has no CUDA runtime
  dependency. Both boundary reporting and continuity scattering call this arithmetic.
- `assemble_outlet_continuity` accumulates interior VMS flux and unique opening
  samples into every continuity row. Inlet samples carry prescribed velocity;
  outlet samples include the frozen trace and RC difference. There is one reverse-add,
  no Dirichlet row removal, no duplicate raw outlet source, and no flux-history update.
- `refresh_outlet_pressure_operator` restores a persistent bare-stiffness snapshot,
  adds compact interior and outlet derivatives using the current frozen coefficients,
  and refreshes the actual `Apre` values. Startup uses BDF1 coefficients; later steps
  use BDF2 coefficients. Additions never accumulate across timesteps.
- The driver builds a separate outlet facet view for matrix-row owners, including
  every facet touching an owned row. Residual facets retain unique ownership.
  Matrix assembly uses all visible incident tetrahedra and writes owned rows;
  residual assembly visits each element/facet once and reverse-adds contributions.
- `mars_outlet_correction.hpp` performs the full-residual correction loop. It reuses
  the existing pressure-increment gradient, including D1's boundary term. It builds
  the VMS context once after diffusion and retains it throughout all corrections.
  Prescribed velocities stay fixed and BDF history advances once per physical step.
- Outlet pressure uses Hypre GMRES. An environment override cannot force PCG onto
  the nonsymmetric outlet matrix. The driver also rejects that conflicting request.

The pressure unknown and increment have physical pressure units. With
`h = dtEff/rho`, the integrated residual and correction are

```
R = interior flux + prescribed inlet flux + stabilized outlet flux
Apre * delta_phi = -R/h
p_trial = p_base + omega * delta_phi
u_trial = u_base - omega*h*Q*G_v(delta_phi)
```

The outlet compact addition is
`-(rho/dtEff)*D_f*(A_f/3 dot grad N_opp)`, in the face-node row and opposite-node
column. It is nonsymmetric. `Apre` approximates the full correction Jacobian
`B Q M^-1 B^T + K_D/h`; its linear residual is never used as the continuity gate.
There is no added pressure pin or mean subtraction in this mode.

## Convergence and failure behavior

The loop measures the volume-weighted RMS of `R_i/V_i`, its maximum over all owned
continuity rows, and the boundary flux imbalance. It checks that `sum_owned R_i`
agrees with the independent boundary reporter. That identity uses the configured
absolute flux tolerance plus a precision-scaled relative allowance for reductions.

A full trial produces `delta_R = R_trial - R_base`. The initial damping is
`min(omega_max, -<R,delta_R>/<delta_R,delta_R>)`, with weights `1/V_i`.
Each accepted update must pass a measured Armijo decrease; rejected trials restart
from the same base state. This is a residual-minimizing correction step, not a
claim that arbitrary meshes admit a convergent stationary iteration.

Momentum failure, invalid geometry/coefficient data, nonfinite residuals, conservation
mismatch, a non-descent direction, failed backtracking, or exhausted corrections
abort the run through MPI. Failure cannot be reported as a completed physical step.
There is no automatic true-J Krylov fallback. If this approximate correction does
not contract on a flow gate, that fallback is a separate numerical extension.

Defaults and units:

| Flag | Default | Meaning |
|---|---:|---|
| `--outlet-max-corrections` | 100 | Maximum pressure corrections per physical step |
| `--outlet-rtol` | 1e-6 | Relative reduction of initial RMS/max and relative flux balance |
| `--outlet-div-tol` | 1e-8 | Absolute RMS/max tolerance, 1/s |
| `--outlet-flux-tol` | 1e-12 | Absolute boundary balance tolerance, volume/s |
| `--outlet-max-damping` | 1 | Upper bound on the measured step, in (0,1] |

RMS and maximum each use `max(absolute_tolerance, relative_tolerance*initial_norm)`.
Boundary balance uses the absolute flux tolerance plus relative tolerance times
the sum of incoming/outgoing magnitudes. The new opt-in `[outlet-continuity]` report
contains full residual norms, balance, correction count, and last damping.
`MARS_SOLVE_TRACE` also prints individual contraction factors.

## Supported first path

The new mode is opt-in:

```
--solver=hypre --vms-stab --rc-implicit --outlet=do-nothing --outlet-beta=0.05
```

It requires constant density, fixed timestep, a velocity-driven inlet, `pumpDp=0`,
`relaxMass=1`, and the full VMS gradient difference. Accepted `relaxU` is preserved.
The first implementation rejects alternative pressure operators, PSPG, RC-only/blend,
opening-flux patches, pressure-driven inlets, open-normal projection, implicit
advection, adaptive timestep, and the legacy `--correctors>1` control. Use the new
correction limit instead. Positive, finite tetrahedral orientation is required.

The ordinary mode remains selected when `outletBeta<0`. Its intended numerical
formula is unchanged; the shared reporting evaluator and necessary device-vector
API fixes also affect code used by that mode. GPU regression validation is pending.

## Executed validation

The C++ host gate executes the production scalar evaluator and line-search helpers:

```sh
clang++ -std=c++20 -O2 -Wall -Wextra -Werror scripts/outlet_boundary_gate.cpp -o /private/tmp/mars_outlet_boundary_gate
/private/tmp/mars_outlet_boundary_gate
```

1,396 host checks pass, including an AddressSanitizer/UndefinedBehaviorSanitizer
build. They cover affine fields, all tetrahedral faces, skewed/scaled geometry,
unequal vertex velocities, face-only coefficients, pressure finite differences,
BDF1/BDF2 scaling, an independently constructed full correction JVP/adjoint,
pressure-level anchoring, bent/opposed areas, and partition algebra. The earlier
unstable compact fixture converges in 39 residual-optimal corrections.

The existing Python trace and boundary-flux gates also pass. The latter remains
an explicitly labeled algebra replica. Host partition tests are not MPI tests.

## Prepared GPU gates and remaining validation

CMake registers and installs:

- `mars_outlet_boundary_gate`: executes the production scalar evaluator on CUDA.
- `mars_outlet_kernel_gate`: executes the actual interior/boundary scatters,
  boundary reporter, compact CSR derivatives, and boundary-aware gradient kernels.
  It compares pressure finite differences and the full pressure/velocity JVP with
  scalar/nodal coefficients and BDF1/BDF2 scaling. Run with 1, 2, and 4 MPI ranks.
  The public one/two-tet fixtures include ranks with no element/facet work.

CTest names are `marsOutletBoundaryEvaluator` and `marsOutletKernels` for single-rank
execution. Use the established GPU binding when launching on the cluster.
The kernel gate uses replicated synthetic arrays and real MPI reductions;
it does not exercise `ElementDomain` halo communication or the full Hypre stepper.

This host has no CUDA compiler. Neither CUDA target nor `mars_pump` has been compiled
or run here. Required before interpreting a changed physical case:

1. Build the solver and both gates on CUDA; execute the scalar gate and 1/2/4-rank kernel gates.
2. Check the actual domain halo, facet ownership, and full stepper on a public synthetic
   mesh, including an empty-opening rank and BDF startup. Compare full residuals and
   physical boundary fluxes across partitions.
3. Run the manufactured pressure-outlet channel specified in the mathematical contract,
   with spatial and temporal refinement and pressure-level checks.

The trace and reconstructed gradient are lagged once per timestep. Convergence certifies
the frozen residual only; refreshed trace/reconstruction defects and BDF2 boundary
accuracy have not been established. No flow improvement or OpenAccel parity is claimed.
