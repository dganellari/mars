# Average-pressure outlet — implementation status

Continuation of `openaccel_math_review_2026-09-09.md`. Author: Claude Opus 5, 2026-09-09.

Editor note — GPT/Codex, 2026-09-09 (forward link only): the status below records the
Claude handoff at `c3cfcc5`. GPT has since completed the integration in the working
tree. See [GPT's implementation and validation handoff](gpt_outlet_implementation_2026-09-09.md)
for current code, ownership, tests, and the remaining GPU validation.

## Landed

| Commit | What |
|---|---|
| `38c27f3` | Empty-rank collective fix (review defect 1), the prerequisite |
| `b3b714b` | The trace as a separate field, plus host gates |

**`38c27f3`** — the `MPI_Allreduce` for the unresolved-facet count sat inside `if (nFacets > 0)`.
A rank owning no opening facets skipped it and paired with the next reduction: hang, or a silently
wrong answer. Now every rank reduces. The full per-facet D2H is gone, replaced by
`unresolvedOpeningFacets` (device `count_if`).

**`b3b714b`** — `p_trace = p_ref + (1-beta)(p_sample - mean_A(p_sample))` in `d_pTraceOutlet`,
**never** written into `d_p`. Area-weighted mean over owned outlet nodes, one `MPI_Allreduce` on
every rank, `p_ref` fallback when the patch is empty or zero-area. Refreshed once per step before
the predictor, so it is frozen for the whole step and `delta p_trace = 0` in the inner Jacobian.
The map lives in one `__host__ __device__` inline, `outletTraceValue`. `boundaryMassFluxKernel`
reads the trace, so operator and diagnostic share it. Flags: `--outlet-beta`, `--outlet-pref`.

Because the refresh happens *before* the predictor, predictor and corrector use the same trace, so
the `G_boundary*(trace_new - trace_used_by_predictor)` lift is identically zero within a step. That
term returns the moment anyone refreshes the trace mid-step.

**Correction (2026-09-09, from `gpt_outlet_boundary_spec_2026-09-09.md` section 3).** That statement
was right about the trace-CHANGE term and wrong to imply nothing else is needed. With
`B = B_i + B_o`, the gradient operator `G_v` itself gains a boundary contribution, so a frozen
trace still leaves

    (G_v phi)_i = (G_0 phi)_i - phi_i * a_o,i / V_i

in the corrector. Setting `phi = 0` on outlet volume nodes would remove it incorrectly. Freezing
the trace removes one term, not the boundary operator.

**Gates run** (`scripts/outlet_trace_check.py`, all pass on host):
prescribed area-weighted mean; `(1-beta)` of the fluctuation retained; gauge invariance under
`p -> p + c`; area weighting is load-bearing; 2/3/4-way partition invariance including a rank with
zero outlet nodes; empty and zero-area patch fall back to `p_ref` without dividing by zero; and
gate 6 — a clamped `p=0` face collapses the trace to a uniform `p_ref`, which is the review's
argument for why a formula-only change is a no-op.

## D progress (2026-09-09, under ultracode)

`--outlet-beta` is **refused by the driver** (`8b7fb24`) until D2+D3 land: without the boundary
derivative the freed rows leave a pure-Neumann, singular system. Everything below is gated on
`outletBeta >= 0`, so the default path is untouched.

| Commit | Piece | State |
|---|---|---|
| `981793a` | D0 one frozen `s.vmsCtx` | done |
| `cd01ef6` | D1 boundary gradient term | done |
| `580d979` | D4 outlet rows freed | done |
| — | **D2 continuity residual** | **not started** |
| — | **D3 boundary derivative** | **not started** |

**D0** — the assembly built its own gradient copy and nodal tau while diagnostics and driver probes
each rebuilt theirs; a rebuild after the corrector samples a different `p`, so reports described a
state the solve never used. `s.vmsCtx` is now built once per step by the assembly and read
everywhere, carrying `dtEff` and the BDF flag so the timestep is frozen with the coefficients.
The legacy `--rhie-chow` path does not populate it and still rebuilds locally.

**D1** — `addOutletGradientTerm`: `g_i += (t_i - p_i) a_o,i / V_i` in the predictor, and the
increment form (`trace = nullptr`, giving `- phi_i a_o,i / V_i`) in the corrector. Masked by Q via
`d_isBdryDof`. Owned nodes only; a halo publish is needed if ghosts ever read these gradients.

**D4** — `d_isPressureBdryDof` is the single control point for the matrix rows, the RHS zeroing and
the lift, so omitting the outlet from it frees all three together.

### D2 — continuity residual, remaining

`computeDivergenceVMSTetKernel` scatters interior SCS only. The outlet needs its boundary samples
scattered into the same `d_divAccNode`, per the spec's nodal lumped quadrature: for triangle
`f=(a,b,c)`, three samples with vector area `A_f/3`, each scattering to its OWN vertex row, using

    q_{f,r} = u_r . (A_f/3) + D_f (gbar_f - grad p_mix) . (A_f/3)
    grad p_mix = p_o grad N_o + sum_r t_r grad N_r
    D_f        = boundaryFaceCoefficient(face nodes only)   [already fixed, 2622cb0]

`boundaryMassFluxKernel` already evaluates exactly this but REDUCES it; D2 needs the same
expression scattering per vertex. Share one evaluator rather than writing a second. Do not also
run the old raw opening term for those facets. Reverse-add owned-facet contributions exactly once.

### D3 — boundary derivative, remaining

    dq_{f,r}/dp_j = -D_f (A_f/3) . grad N_o * delta_{jo}

row = face node, column = opposite node — both in the same tet, so the CSR entries already exist
(the rows are assembled, then were overwritten). Scale by `rho/dtEff` before adding into
`d_valuesPre`, and add it BEFORE the values are copied into the active wrapped matrix; rebuild the
Hypre operator and preconditioner afterwards. `A_f . grad N_o < 0` on an outward face, so the
entries are positive in the opposite-node column: nonsymmetric, not a diagonal Robin penalty. Use
Hypre GMRES, not PCG.

### Gates, none of which have run

Full continuity-residual contraction (measured, not assumed — the spec's counterexample gives
spectral radius 12.05 undamped and 2.26 at omega=0.25); the actual CUDA evaluator checks; 1/2/4-rank
conservation including an empty-outlet rank; and the manufactured channel before any changed
physical case. Host algebra gates that DO pass: `scripts/outlet_trace_check.py` (9 groups) and
`scripts/outlet_boundary_flux_check.py`.

## Superseded: the solve was unchanged at b3b714b
## NOT done: the solve is unchanged

Outlet pressure rows are still identity. With them clamped, `p_sample = 0`, so the trace is
uniformly `p_ref` — correct behaviour for the current configuration, and useless until the rows
are freed. **Do not report the trace as a result until this is finished.**

Four sites have to change together, gated on `outletBeta >= 0` so `--outlet=do-nothing` stays
byte-identical for the controlled comparison:

1. **Setup** — stop marking outlet DOFs in `d_isPressureBdryDof` (`mars_ns_pump_solver.hpp`, the
   `BCK::Pump` branch), so `enforceBcMatrixKernel` leaves those rows assembled.
2. **RHS** — retain the volume divergence at outlet nodes (currently overwritten by
   `enforcePressureBcRhsKernel`) and add the boundary flux contribution. The facet infrastructure
   exists but `addFemOpeningSurfaceTerm` is gated on `useFemProjection && !useOpeningFluxSource`,
   neither of which holds in the production config.
3. **Matrix** — the paired derivative on non-prescribed contributions,
   `d q_b / d p_j = -sum_d D_b,d A_b,d d_d N_j`. The sparsity entries already exist (the rows are
   assembled, then overwritten), so this is a value change, not a pattern change.
4. **Gauge** — 3 is not optional. Removing the Dirichlet rows leaves a pure-Neumann system; the
   boundary derivative is the only Robin closure pinning the level against `p_ref`. Verify
   nonsingularity before trusting a solve, and decide whether a fallback pin is still wanted.

Complication: the pressure solve has three variants (`PressureSolveKind::K`, `DDT`, and the
FemGram path) with symmetric-Jacobi scaling and a separate AMG preconditioner matrix. Establish
which is active under the production flags before editing, and change one.

**Third finding from the spec, section 5.** Contraction of the compact approximate correction is
NOT automatic. Their matrix counterexample (unit tet, `h=1`, `D=0.3`) gives eigenvalues of
`A_tilde^-1 J` up to `13.05`, so the undamped stationary iteration has spectral radius `12.05` and
even `omega=0.25` leaves `2.26`. So `Apre` cannot be assumed to work as a stationary corrector;
it may still serve as a GMRES preconditioner for the true `J`. Measure residual contraction before
trusting any result from this path.

CUDA/MPI gates still pending: frozen-state Jacobian action vs its declared linearization; owned
continuity sums equal exterior flux sums; 1/2/4 ranks agree including empty-outlet ranks; a
manufactured pressure-outlet channel as the first flow gate.
