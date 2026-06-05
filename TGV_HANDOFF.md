# TGV / Periodic NS — Handoff for Next Session

## UPDATE 2026-06-05 (corrector D^T adjoint fix, commit 75b454a) — the verified seam-leak fix

History on this seam residual: collapse (4096) correct; then TWO refuted fixes — advN broadcast
(inert, deleted) and post-corrector velocity broadcast (bea6d38, REFUTED by run: div_max still
3.5@40/12@50, reverted 220375b). Seam VELOCITY drift was NOT the cause.

ROOT CAUSE (mechanically confirmed by reading the code, not theory): the reduced-periodic CORRECTOR
gradient g=M^-1 D^T phi (ns_solver.hpp ~7850) was NOT the exact adjoint of the D the reduced pressure
CG inverted, AT CROSS-RANK SEAM MASTERS. The OPERATOR applyDDTReduced does, at its g=M^-1 D^T phi
stage: applyDivTranspose -> reverseHalo -> **maybePeriodicSum (cross-rank fold)** ~5072 -> normalize
-> **periodicBroadcastSameRank + crossRankPeriodicBroadcast** ~5117/5133 -> exchangeHalo. The corrector
REPLICATED only applyDivTranspose -> reverseHalo -> normalize -> exchangeHalo, DROPPING the fold and
the broadcast (its comment even claimed "matching the reduced operator exactly: NO maybePeriodicSum"
-- which is wrong; the operator DOES maybePeriodicSum). So for a cross-rank pair (master+slave = two
distinct owned DOFs), the operator used the MERGED both-half-CV value while the corrector applied each
half SPLIT and separately mass-normalized -> M^-1 D^T(corrector) != M^-1 D^T(operator) at the seam ->
the projection identity D(u** - dt/rho M^-1 D^T phi)=0 leaks a small div at the seam EVERY step,
accumulates, blows ~step 40. This explains ALL observables: slow monotone growth (not step-0),
seam-localized (div-split 3-8), single-rank perfectly clean (no cross-rank pair exists), and both
velocity fixes refuted (the leak is in the gradient restriction/normalization pairing, not velocity).

FIX (75b454a): the corrector gradient now does maybePeriodicSum x3 (fold) before normalize AND
periodicBroadcastSameRank + crossRankPeriodicBroadcast x3 after normalize -- byte-for-byte mirroring
the operator's g=M^-1 D^T phi stage, so the corrector's D^T is the exact adjoint of the operator's D
at the seam. Single-rank safe (maybePeriodicSum/broadcast are no-ops with no cross-rank pairs).

VERIFY (4-rank cube16 --solver=cg --pressure-solve=DDT --num-steps=110): div(u_n) ENTRY should STOP
growing (was 1.7e-3@10 -> 1.3e-2@30 -> blow@40), div-split ratio relax toward ~1, survive 110 steps,
KE match single-rank. Then --skew=0 (scheme-independent check) and --solver=hypre (same operator, the
assembled DDT path -- this fix is on the matrix-free corrector so Hypre uses the same corrector and
should also benefit). If it STILL leaks: add the [div-global] signed/abs MPI_Allreduce probe (in the
wpruqm49u workflow output) -- signed~0 & abs grows => still a closure dipole; signed drifts => fold
non-telescoping (H3).



## UPDATE 2026-06-05 (advection-seam fix, commit bea6d38) — the last isolated bug

After the true collapse (4096) the only remaining blow-up (~step 40-50, scheme-independent) was the
ADVECTION seam. Root cause (verified in code, not guessed): the post-corrector velocity broadcast
master->slave on s.d_u/v/w (ns_solver.hpp ~8269) was gated OFF by `!routeReducedPeriodicCorr` --
a leftover from the EMULATION era (when the slave owned a row and had its "own projected velocity"
to protect). Under owner-migration the slave is a GHOST with no own projection: u^{n+1}[slave] must
equal u^{n+1}[master]. With the broadcast skipped, the seam slave/ghost velocity drifted from its
master by one step, and the next predictor's advection flux scatter (~6932) integrated the drifted
u^n -> a seam momentum imbalance in mdot=vf.areaVec (common to upwind AND skew -> scheme-independent,
matches the --skew=0==--skew=1 discriminator) that grew div(u^n) every step.

FIX (bea6d38): (A) un-gate the post-corrector velocity broadcast so u^{n+1} syncs master->slave every
step (seam-consistent u^n at the next predictor). (B) deleted the inert advN master->slave broadcast
(predictor reads advN only at owned master DOFs, dof<numOwnedDofs; the slave-ghost slots it wrote
were never read -- advN[master] is already the full-CV sum from reverseExchangeNodeHaloAdd +
maybePeriodicSum). Single-rank safe: gate was numRanks>1; single-rank broadcast is master->same-rank-
slave within rank, harmless. NOTE the prior afc48e9 u^n re-sync regression was WITH the emulation
still present (it fought the reduced-DOF slave projection); that world is gone, so re-enabling the
broadcast is now correct.

VERIFY (run 4-rank cube16 --solver=cg --pressure-solve=DDT --num-steps=110): div(u_n) entering each
step should STOP growing (was 0.49@35 -> 2.5@40 -> 7@50) and stay O(1e-3); div-split ratio should
relax toward ~1 (was stuck 2.5-3.5); survive 110 steps; KE matches single-rank. Then --skew=0 and
--solver=hypre should ALSO survive (same root cause). If it STILL grows, the secondary suspect is
seam areaVec orientation (the [advN-sum] owned total would differ 1-rank vs 4-rank) -- but velocity
drift is the verified primary.



## UPDATE 2026-06-05 (MILESTONE) — TRUE cross-rank collapse landed and is CORRECT (sum numOwnedDofs = 4096)

Implemented owner-migration true cross-rank periodic collapse (commit 612faa2): a cross-rank slave
no longer owns a DOF -- it migrates onto its master's GHOST dof, so the periodic pair is ONE global
DOF. ONE predicate change (drop `d_own[master]==1` guard in the DOF-collapse Pass-1, ns_solver.hpp
~2706) + one Hypre-seed guard (`dof < nOwn`, ~4544) + DELETED all the cross-rank emulation (503
lines: per-matvec crossRankPeriodicPairSumDof, crossRankPeriodicBroadcastDof, dot-mask,
preconditioner-diag override, Path-B XR-slave masks, assembled-DDT slave-col fold). The advection-
seam fixes (advN restore, corrector group-average, seam mass mirror) were KEPT (separate term).

VERIFIED CORRECT on 4-rank cube16 (--solver=cg --pressure-solve=DDT):
- `[owned-node check] sum(numOwnedDofs over ranks) = 4096` -- EXACTLY the single-rank global unique
  DOF count (was 4657 with the emulation double-count). This is the smoking-gun proof: the periodic
  pairs are genuinely one DOF each, no orphan, no double-count.
- Step-0 KE = 0.125 = V0^2/8 (the correct physical TGV value) on 4-rank, matching single-rank.
- Solves healthy: cg_uvw=4, cg_p~115, DDT CG converges to ~1e-9 every step.
The cross-rank DOF machinery is now DONE and provably correct. No more emulation.

### REMAINING (now cleanly ISOLATED): the advection-seam residual, NOT the DOF coupling.
The 4-rank run still blows up ~step 40-50 (div_max ~2.5 @ step40, ~7 @ step50), with div-split
ratio ~2.5-3.5 (boundary RMS / interior RMS) that does NOT fall to ~1, and div(u_n) entering each
step grows (0.003 -> 0.49 by step 35). Because the DOF collapse is now PROVABLY exact (4096), this
residual is NOT cross-rank DOF coupling. It is the ADVECTION-SEAM term -- consistent with the earlier
--skew=0 discriminator (upwind AND skew blow up identically -> scheme-independent seam-flux residual).
NEXT: attack the advection seam directly (the advN/momentum-flux consistency at the cross-rank seam),
now that everything else is eliminated. The single-rank run is clean (div~1e-14), so the seam term is
purely a multi-rank advection-consistency issue at the periodic boundary.

### BONUS to check: --solver=hypre should now ALSO work (or get much further). The half-strength
seam-row problem that made Hypre blow up at step 0 is GONE -- there is no slave row anymore, the
master row is the single complete row. Re-run `--solver=hypre --pressure-solve=DDT` 4-rank; the
step-0 velocity-solve garbage should be fixed. (The advection-seam residual may still cap it ~step
40-50 like CG, but the catastrophic step-0 failure should be resolved.)

---

## UPDATE 2026-06-05 (CORRECTION, SUPERSEDED by the milestone above) — Hypre periodic path WAS broken pre-collapse; CG matrix-free was the only viable route

TESTED `--solver=hypre --pressure-solve=DDT` on 4-rank cube16 periodic: it BLOWS UP FROM STEP 0.
The VELOCITY diffusion solve (not pressure) returns garbage immediately: after DIFF, |u**|=2.34
(should be ~1), div(u**)=5.76, cg_uvw=2/2/2 (Hypre "converged" in 2 iters to a wrong answer) ->
|phi|=9048 -> NaN -> "RHS contains NaN/Inf" abort. (Step-0 KE=0.125 IS correct now = V0^2/8, so the
KE-diagnostic slave-skip fix works; the field starts right and the velocity solve corrupts it.)

WHY (the fundamental reason, do not retry Hypre for periodic without solving this): the Hypre
velocity solve uses the ASSEMBLED nu*K matrix whose cross-rank periodic seam rows are HALF-STRENGTH
(the master row lacks the slave-side stiffness -- the same bug the CG path had). The CG path fixed
it with Strategy-B: per-matvec crossRankPeriodicPairSumDof in the spmvPostCallback, reconstructing
the merged seam operator EVERY iteration. Hypre is a black-box AMG solver -- it CANNOT call a
per-matvec callback, so it never gets the seam merge and diverges on the inconsistent assembled
operator. This matches the pre-existing memory note "Stage 4 (Hypre AMG) reverted pending RHS-NaN
debug" / "DDT+Hypre validation pending" -- Hypre periodic was ALREADY known-broken; we re-confirmed it.

CONCLUSION: for multi-rank periodic, the matrix-free CG path (--solver=cg) is the ONLY viable
architecture -- it is the only one that can apply the per-matvec seam reconstruction the periodic
operator requires. Hypre would need the seam baked into the ASSEMBLED matrix (true global-DOF
collapse: skip the slave row entirely, alias its column to the master's global id), which the
current Hypre IJ wrapper cannot do (it assigns row ids positionally). That is a deeper IJ-wrapper
change, not a config flag. Finish the CG path instead (the advection-seam residual is the last item).

(Below was the earlier optimistic note that the Hypre route was "already correct" -- it is NOT;
kept for context but superseded by the test above.)

---

## UPDATE 2026-06-05 — TWO periodic paths exist; Hypre route is the genuinely-correct one (already wired)

After eliminating geometry (area-vectors), edge/corner cross-rank pairing (probe FRAGMENTED=0),
and the step-0 KE mismatch (an 8/9 DIAGNOSTIC double-count, now fixed — KE skips slaves) as the
cause of the matrix-free residual, I read the Hypre path. Key findings:

- `--solver=cg` (default): `routeReducedPeriodic` matrix-free P^T A P, gated to `solverKind==CG`
  (ns_solver ~7842). This is the path we refined all session — stable ~40-50 steps, div-split ~2.5,
  then a small scheme-independent seam residual blows it up. The cross-rank pairs here are NOT
  collapsed: solved as TWO DOF rows, forced equal post-solve (the emulation; its residual is the
  remaining bug). KEEP THIS PATH (user wants both).
- `--solver=hypre`: `SolverKind::Hypre` -> `routeReducedPeriodic` is FALSE -> takes the ASSEMBLED
  DDT branch (ns_solver ~7848, uses `s.AddT`). For periodic, the assembled cross-rank slave row is
  made a Dirichlet identity AND its column is FOLDED into the master row
  (`foldPeriodicSlaveColIntoMasterKernel`, ns_solver ~4385), so the master equation carries the
  merged physics in an EXACTLY-symmetric assembled operator; `crossRankPeriodicBroadcastDof`
  restores x[slave]=x[master] post-solve. THIS IS the genuinely-correct reduced-DOF treatment
  (MFEM/deal.II style), and it is ALREADY FULLY WIRED — no code change was needed.

NOTE I was wrong earlier that "Hypre uses true global numbering so sidesteps the seam" — the Hypre
IJ wrapper assigns ROW global ids positionally (ilower+i) and only remaps COLUMNS, so it CANNOT
merge two rows by aliasing; that is exactly why the code uses the identity-row + column-fold (Path
B) instead. The fold makes the assembled operator correct; this is more robust than the matrix-free
per-matvec emulation.

### ACTION to get a working multi-rank periodic solver: build with Hypre and run --solver=hypre.
```
# build (Alps): add -DMARS_ENABLE_HYPRE=ON to the cmake config, rebuild mars_tgv.
# run 4-rank periodic via the assembled-DDT + BoomerAMG path:
srun ... mars_tgv --mesh=cube16.mesh --box-lo=0 --box-hi=1 --V0=1 --nu=0.05 --dt=1e-4 \
    --num-steps=110 --pressure-solve=DDT --solver=hypre --skew=1 --adapt-every=0
```
Success = div_max stays ~roundoff, KE decays, survives 110+ steps (the assembled fold is exact, so
no seam emulation residual). If Hypre BoomerAMG struggles on the near-singular pure-Neumann DDT
operator, the DDT eps diagonal shift (MARS_DDT regularization, ns_solver ~3235) and the pin/removeMean
are already in place. The CG matrix-free path stays the no-Hypre fallback (stable ~40 steps).

Both `--solver=cg` and `--solver=hypre` are preserved (user requirement).



## UPDATE 2026-06-04 (LATE) — advection seam fix landed; one secondary seam source remains

Best state now = commit **afc48e9** (advN restore + revert of the u^n re-sync that regressed).
- The ADVECTION seam fix (ea3dd76) was the big one this round: maybePeriodicSum(advN) folds
  sum-only (broadcastBack=false) leaving advN[slave]=0, so the predictor gave the cross-rank
  slave ZERO advection increment -> qStar[slave]!=qStar[master] every step (standing seam jump).
  Restoring advN[slave]:=advN[master] after the fold (mirroring the corrector group-average)
  dropped div-split from ~6-7 to ~2.2-2.5, stabilized |p| (~0.4-0.58, |phi| tiny) cleanly
  through ~step 40, and is divergence-neutral. KEEP.
- TRIED and REVERTED (b9c842f -> afc48e9): re-affirming u^n[slave]=u^n[master] at PREDICTOR
  ENTRY (re-enabling the gated-off start-of-step velocity broadcast on the reduced path). It
  REGRESSED (div-split back to ~5-6, blew up earlier). So u^n input seam-inconsistency was NOT
  the secondary source; do not retry that.
- REMAINING: still blows up ~step 40-50. div-split is ~2.2-2.5 (not ~1) and div(u_n) entering
  still creeps. A smaller secondary seam source remains. Per the advection audit, the leading
  remaining candidate is the skew form's -1/2 q (div u) term: while div(u^n) is nonzero at the
  seam it pumps KE there (positive feedback). NEXT: confirm with the --skew=0 (upwind) run — if
  upwind survives 110 steps where skew blows at ~50, the residual is the skew energy-conservation
  term at the seam, and the fix is either (a) make the seam div(u^n) the corrector sees actually
  ~0 (the projection closes interior but div-split ~2.5 says the seam still carries div), or (b)
  use the pure-divergence advection form at seam faces. Do NOT keep re-broadcasting velocity
  fields (that direction regressed).

Lesson: each seam fix this session pushed the blowup later (immediate -> 30 -> 40 -> ~50) by
removing one seam-consistency defect; the advN fix was real and large. Stopped iterating when a
fix regressed — bank afc48e9, run the --skew=0 discriminator next to target the last term.

### --skew=0 DISCRIMINATOR RESULT (run 2026-06-04, decisive):
UPWIND (--skew=0) blows up at essentially the SAME step (~40-50) as SKEW, with the SAME div-split
(~2.5 early) and the SAME div(u_n) creep (0.001 -> 0.05-0.09 by step 30). So the residual
instability is NOT scheme-specific — it is present in BOTH advection forms. This RULES OUT the
skew -1/2 q(div u) energy term as the cause (which only exists for skew). Per the audit's own
mapping, "upwind also blows up" points at suspect D: a residual SEAM-FLUX inconsistency common to
any advection scheme — area-vector ORIENTATION or periodic-image-FACE double-count/sign at the
cross-rank seam. NOTE: advN-sum being rank-invariant (a global scalar) CANNOT detect a sign-flipped
pair of seam faces (they cancel in the sum), so the audit's earlier desk-refutation of D is not
valid; D needs a DIRECT per-face check.

### NEXT SESSION — target suspect D directly (do NOT add more field broadcasts; the velocity
### re-sync already regressed):
1. The advection/divergence read areaVec[off] per element-face (mars_ns_solver.hpp ~837/988/1081,
   scsLR at 654). The periodic-image halo element (delivered by cstone's periodic Box) sits at its
   REAL shifted coords; its precomputed area vector must have orientation CONSISTENT with the
   on-rank element sharing that physical seam face. Check the area-vector precompute
   (precomputeAreaVectorsGpu / wherever d_areaVec_{x,y,z} is filled) for whether a periodic-image
   element's SCS face normal points the same physical direction as its on-rank twin. A sign flip
   there is a per-face seam error invisible to advN-sum.
2. Cheap probe: instrument max|advN[slave]-advN[master]| right AFTER the advN restore (it should be
   ~0 now); and separately dump, for one known cross-rank seam face, the areaVec on both the owning
   rank and the halo rank — they must be equal-and-opposite (outward from each element) i.e. the
   shared face flux cancels. If they don't, that is D.
3. If area vectors are clean, the remaining candidate is that the seam div(u^n) the corrector
   leaves (~div-split 2.5) is real (not a measurement artifact) and the projection simply is not
   closing the seam to roundoff — revisit whether the reduced operator's D and the corrector's D
   are EXACT transposes at the seam ONE more time, but only after D (area vectors) is cleared.

---

## UPDATE 2026-06-04 — multi-rank periodic: catastrophic blowup fixed, standing seam instability remains (NEW diagnosis)

Branch `cstone`. 4-rank cube16 periodic TGV (`--pressure-solve=DDT --skew=1 --dt=1e-4`).

### Where it stands
- Single-rank periodic TGV: correct (div_max ~3e-6 flat, KE decays, |p|~0.4). This is the reference.
- 4-rank: was an instant catastrophic blowup (|p| 0.13 -> 7447). Now STABLE for ~40-50 steps
  (|p| flat ~0.41, KE decays cleanly, cg_uvw=3-4, cg_p~114), THEN diverges. So all the
  catastrophic-mode bugs are fixed; a slower standing instability remains.

### The 7 fixes that landed this session (all real, all KEEP — they made the pressure path
### + the assembled velocity operator fully cross-rank periodic-consistent):
1. RHS divergence fold order (518662d)
2. predictor slave-velocity broadcast gate (f7d432b)
3. post-solve phi master->slave prolongation + MARS_PROJ_PROBE (e776ad2)
4. **seam mass mirror** — d_massNode[slave] was 0, zeroing the slave gradient; mirror the
   combined pair mass to both slots (bf822ed). THIS killed the catastrophic blowup.
5. velocity operator cross-rank consistency — XR-slave Dirichlet-identity rows + re-enabled
   post-solve broadcast (cd1f0cd); then **Strategy B**: keep half-strength seam rows, merge
   Ap per-matvec via crossRankPeriodicPairSumDof in the CG spmvPostCallback + dot mask
   (31ffb59, f57bc00). Dropped cg_uvw from 19-20 to 4.
6. velocity Jacobi diagonal + RHS seam merge for Strategy-B parity (a41bcb2).
7. corrector seam velocity GROUP-AVERAGE — periodicDivideAtMastersKernel; sets all members
   of a periodic group (master+slaves, edges=2 corners up to 8) to the group mean, which is
   divergence-neutral (eaef201).

### THE KEY NEW FINDING (do not repeat the projection grind):
A dt-halving test REFUTED the "O(dt) projection/transpose accumulation" theory:
- dt=1e-4  -> blows up ~step 50, physical t~0.005
- dt=2.5e-5 (4x smaller, 4x more steps, same T) -> blows up ~step 80-90, physical t~0.002 (EARLIER)
A true projection-consistency error would blow at the SAME or LATER physical time with smaller
dt. It blew EARLIER. So the remaining instability is **NOT** a pressure-projection / D-vs-D^T
transpose error. The pressure operator is correct (symmetric, CG->1e-9, |A phi - b|/|b|~9e-11).

What the data says instead: the `[div-split]` boundary/interior RMS ratio is pinned ~6-7 from
step 0 in BOTH runs (dt-independent), and div(u_n) ENTERING each step grows in physical time
while all solves stay healthy. That signature = a **standing instability at the periodic seam
in the ADVECTION / momentum-flux term** (the SCS face mdot / u.grad u across the seam faces),
or an unstabilized checkerboard mode the periodic faces don't damp — NOT the projection.

### NEXT SESSION SHOULD investigate (in order):
1. The advection / SCS mass-flux (mdot) treatment at cross-rank periodic faces: is the
   convective flux across the seam computed consistently (same mdot both sides, periodic-image
   element seen once)? Compare the advection assembly's cross-rank/periodic handling to how the
   diffusion operator was fixed (Strategy B). advN-sum diagnostic is already rank-invariant, so
   the GLOBAL advection sum is fine — the issue is LOCAL at the seam.
2. Whether the skew-symmetric advection (`--skew=1`) needs the same seam group-consistency the
   velocity got; try `--skew=0` (upwind) as a cheap discriminator — if upwind is stable, it's a
   skew-advection seam energy-conservation issue.
3. Rhie-Chow / Bochev-Dohrmann stabilization behavior at the periodic seam (boundary RMS 6-7x
   interior suggests the seam faces aren't getting the same pressure-velocity coupling damping).

### Do NOT:
- Chase more pressure-projection / D^T-transpose seam fixes (refuted by the dt test).
- Switch to Hypre for this (Hypre is a pressure-operator alternative; the pressure operator is
  not the problem — it would inherit the advection-seam instability).

### Diagnostics available (env-gated, default off):
MARS_PROJ_PROBE (P1 |A phi-b|/|b|, P3 ||Du^{n+1}||/||Du**||, seam|u_s-u_m|),
MARS_PERIODIC_XR_SYMCHECK (assembled velocity seam row symmetry: A_local/A_peer).

### Tutorial: docs/periodic_tgv_tutorial.md (beginner-level, real code refs).

---

Status as of 2026-05-24. Branch: `cstone`.

## What works (verified, do not undo)

### Shared NS solver header
- **`backend/distributed/unstructured/fem/mars_ns_solver.hpp`**:
  full extraction of NSStepper, setup, runNsStep from
  `mars_amr_ns_projection.cu`. Cavity and channel drivers build against it
  with no behavioural change. Cavity regression verified to give bit-identical
  output to pre-refactor (cg_iter_uvw=5, cg_iter_p=92).
- `runNsStep` is split into four composable sub-functions:
  `runPredictorStep`, `runImplicitDiffusionStep`, `runPressureSolveStep`,
  `runCorrectorStep`. Each callable in isolation — gives a hook for unit
  testing the pressure-Poisson sub-step against known velocity fields. This
  was an explicit ask for the DDT-debug effort.

### Periodic BC infrastructure (NEW)
- **`backend/distributed/unstructured/fem/mars_periodic_bc.hpp`**:
  builds the slave→master partner map by geometric matching across opposite
  faces of an axis-aligned box. Verified on `cube16.mesh` to find 817 slaves
  (= 3×17² − 3×17 + 1 corner-flatten = exactly correct).
- `NSStepper::BCKind::Periodic` mode added. Plumbing in setupNSStepper:
  - Skip Dirichlet velocity BC mark, skip pressure pin, skip pressure
    matrix-row enforcement
  - **DOF collapse + compaction**: `nodeToDof[slave] = nodeToDof[master]`,
    then exclusive-scan renumber so active DOFs are dense in `[0, nLive)`.
    s.numOwnedDofs becomes `nLive` (e.g. 4096 on cube16, was 4913). No
    orphan rows in the sparsity. Avoids the "matrix has identity rows
    at slave DOF indices" hack.
  - `removeMean()` after pressure solve to fix null-space gauge
  - **RHS mean projection** before pressure CG (K and DDT paths both): the
    pure-Neumann Laplacian has constant-mode null space and CG needs a
    consistent RHS.
- **Hex kernel diagonal fix** in
  `backend/distributed/unstructured/fem/mars_cvfem_hex_kernel_tensor.hpp`:
  was checking `i == j` (corner indices), now checks `col_dof == row_dof`.
  Identical in cavity (1:1 nodeToDof) but critical under collapse where
  two corners can map to the same DOF.
- **Atomic mass-gather fix** in `mars_ns_solver.hpp`
  `gatherOwnedNodeMassToDofKernel` — was plain assignment, now atomicAdd.
  Latent bug in cavity (idempotent), fatal under collapse.

### TGV driver
- **`examples/distributed/unstructured/mars_tgv.cu`**:
  loads a periodic-cube mesh, builds the periodic map, sets up NSStepper
  in `BCKind::Periodic` mode, applies the TGV initial condition
  `u = V0*sin(kx)*cos(ky)*cos(kz)` etc. with `k = 2π/L`, removes pressure
  mean, runs the time loop, writes per-step KE / |u| / |p| / div_max /
  CG iters, writes VTU/PVTU/PVD every N steps.
- CLI flags: `--mesh=`, `--V0=`, `--nu=`, `--Re=`, `--dt=`, `--num-steps=`,
  `--box-lo=`, `--box-hi=`, `--vtu-output=`, `--vtu-every=`,
  `--pressure-solve=K|DDT`, `--solver=cg|hypre`.
- Diagnostics infrastructure for per-phase max-abs of every relevant
  field (gated by `g_nsDebugStepsLeft`).

### Verified end-to-end behavior
- Run on cube16, single rank: periodic map = 817 slaves (matches theory).
- Setup: 4096 active equations (was 4913 with slaves redundant).
- Step 0 KE = 0.1322530 = V0²/8 (exactly the TGV analytic value for a
  unit cube — confirms IC + lumped mass are mutually consistent).
- Step 1 with DDT: `div_max = 3e-12` — projection is at machine
  precision. Pressure CG converges in ~92 iters, velocity CG in 4-5.
- Step 10 with DDT: still bounded, `div_max = 5e-12`. Pipeline works
  end-to-end.

### Build
- New target `mars_tgv` added to
  `examples/distributed/unstructured/CMakeLists.txt`. Tested green.

## Open problem (do NOT iterate on this with rsync-based debug)

**TGV blows up around step 30** with any of:
- 1st-order upwind + forward-Euler advection (default) — `KE: 0.1582`
- Rotational pressure correction (Timmermans/Guermond)
  `p += phi − ν·div(u*)` — `KE: 0.4370` (worse)
- Central skew-symmetric advection `q_face = ½(q[L] + q[R])` — `KE: 94.7` (much worse)

The diagnosis: **the explicit forward-Euler advection scheme is the actual
stability bottleneck**, not the spatial discretization or the boundary
correction. Cavity stays stable because Dirichlet walls dissipate the
numerical noise. Periodic NS has no wall dissipation, so the noise
compounds exponentially regardless of spatial form.

## What every serious code does instead (the real fix)

NekRS / Nek5000 / MFEM-Navier / deal.II step-35 all use the
Karniadakis-Israeli-Orszag (KIO 1991) split scheme:

- **BDFk** time discretization on the implicit (diffusion) side:
  k=2: `(3 u^(n+1) − 4 u^n + u^(n-1)) / (2 dt) = ...`
  k=3: `(11 u^(n+1) − 18 u^n + 9 u^(n-1) − 2 u^(n-2)) / (6 dt) = ...`
- **EXTk** extrapolation on the explicit (nonlinear) side:
  k=2: `2·(u·∇u)^n − (u·∇u)^(n-1)`
  k=3: `3·(u·∇u)^n − 3·(u·∇u)^(n-1) + (u·∇u)^(n-2)`
- Second/third-order accurate in time, stable up to CFL ~ 1, no
  boundary-splitting error issues at all.

## Concrete implementation plan for BDF2/EXT2

1. **Add history vectors to NSStepper**:
   - `d_u_prev`, `d_v_prev`, `d_w_prev` — velocity at step n-1
   - `d_advU_n`, `d_advV_n`, `d_advW_n` — advection RHS at step n (to use
     at step n+1 via extrapolation)
   - `d_advU_nm1`, `d_advV_nm1`, `d_advW_nm1` — at step n-1
2. **Modify `runImplicitDiffusionStep`**:
   - LHS: `(M · 3/(2 dt) + ν K)` instead of `(M/dt + ν K)`
   - RHS: `M · (4 u^n − u^(n-1)) / (2 dt) + (2·advN − advNm1)`
   - For the first step (no n-1 history) fall back to forward Euler
3. **AMR-aware migration**: when AMR refines/coarsens, the history vectors
   need to be transferred too. The existing `AmrManager::adaptMeshMultiField`
   takes a vector of fields — pass the BDF history alongside the primary
   velocity/pressure.
4. **Rebuild velocity matrix**: when `dt` doesn't change, no rebuild
   needed (existing behaviour). When `dt` changes (AMR-adapted CFL), rebuild
   once. NSStepper already has the machinery — just expose a
   `setupNSStepper(...)`-equivalent that only rebuilds the velocity matrix.
5. **Skew-symmetric advection kernel**: once BDF2/EXT2 is in, the skew-
   symmetric form (which made things worse with forward-Euler) becomes
   the right choice. The combination is what NekRS uses.

Estimated scope: ~800 lines, 2-3 sessions of compile-test iterations.

## Knobs to flip when validating

In `mars_tgv.cu`, after setupNSStepper:
```cpp
// Once BDF2/EXT2 is implemented:
s.timeIntegrator = NSStepper<...>::TimeIntegrator::BDF2_EXT2;
// And the central skew-symmetric form is now safe:
s.useSkewSymmetricAdvection = true;
```

## Useful diagnostic commands

Smoke test (single rank, no AMR, write IC + 10 steps):
```
srun --account=csstaff --time=00:03:00 --nodes=1 --ntasks-per-node=1 \
  ~/affinity/bind_numa.sh \
  /capstor/scratch/cscs/gandanie/git/mars/daint-gpu/examples/distributed/unstructured/mars_tgv \
  --mesh=/capstor/scratch/cscs/gandanie/git/mars/cube16.mesh \
  --box-lo=0 --box-hi=1 --V0=1 --nu=0.05 --dt=1e-4 --num-steps=10 \
  --vtu-output=/tmp/tgv_smoke --vtu-every=5
```

Cavity regression (verify shared header changes don't break Dirichlet):
```
srun --account=csstaff --time=00:03:00 --nodes=1 --ntasks-per-node=1 \
  ~/affinity/bind_numa.sh \
  /capstor/scratch/cscs/gandanie/git/mars/daint-gpu/examples/distributed/unstructured/mars_amr_ns_projection \
  --mesh=/capstor/scratch/cscs/gandanie/git/mars/cube16.mesh \
  --bc=cavity --dt=1e-3 --num-steps=20
```

Expected cavity: `cg_iter_p=92 cg_iter_uvw=5`, `|u| ~ 0.19` at step 20.
If those change, the header refactor broke something.

## Files touched in this work

```
backend/distributed/unstructured/fem/mars_ns_solver.hpp     [NEW, 2900 lines]
backend/distributed/unstructured/fem/mars_periodic_bc.hpp   [NEW, 380 lines]
backend/distributed/unstructured/fem/mars_cvfem_hex_kernel_tensor.hpp [1-line bug fix]
examples/distributed/unstructured/mars_tgv.cu               [NEW, 460 lines]
examples/distributed/unstructured/mars_amr_ns_projection.cu [slimmed: was 2500, now 360 lines]
examples/distributed/unstructured/CMakeLists.txt            [+mars_tgv target]
```

---

## 2026-05-28/29 update — multi-rank periodic deep dive

### Verified facts (do NOT re-litigate)

1. **cstone IS sufficient for periodic halo delivery.** Empirical diagnostic
   on cube16/4-rank shows 1238 halo elements per-rank-sum touch the periodic
   min-faces (xmin/ymin/zmin). cstone's `applyPbc` in `boxoverlap.hpp:270` is
   in the halo-discovery path. `MARS_NODEHALO_HOST=1` is NOT required for
   periodic; the v2 GPU peer discovery works through cstone's periodic Box.
2. **Stale comment in `domain.hpp:542-545`** claiming v2 doesn't work for
   periodic has been deleted.
3. **The lumped-mass / divergence / gradient kernels used to skip halo
   elements**. They looped `[startElem, startElem+numLocal)` (owned-only)
   while the matrix assembler loops `[0, s.elementCount)` (owned+halo with
   an `ownership[node]==1` write filter). Lumped mass was patched to the
   matrix-assembler pattern (write owned-only nodes inside an owned+halo loop,
   no reverse-halo needed). Divergence / gradient / divT kernels still use
   `[0, s.elementCount)` but WITHOUT the ownership filter inside (they over-
   count after reverseExchangeNodeHaloAdd). **Apply the same ownership-filter
   pattern to those 4 kernels** as the next step.

### Remaining bug (multi-rank periodic cube16 4-rank)

Symptoms:
- `KE step-0 = 1.39e-01` (multi-rank) vs `1.58e-01` (single-rank baseline).
  Mass kernel fix improved correctness but didn't close the gap, because
  cross-rank periodic-pair DOFs are still TWO separate own=1 equations.
- `cg_uvw=-2/-2/-2` (pAp<=0 breakdown) inside the velocity Helmholtz CG.
- KE collapses to 0 by step 10.

Root cause (high confidence):
- For a periodic slave on rank A with master on rank D, the slave is own=1
  on rank A (its own DOF) and the master is own=1 on rank D (its own DOF).
  These are two separate equations for the same physical point. The matrix
  IS coupled through cstone's periodic ghost halo (sparsity builder sees
  periodic-image halo elements; emits slave-row x master-ghost-column entries)
  AND the existing periodic-broadcast in the CG callback syncs p[slave]<-
  p[master_ghost] before each SpMV. Yet CG still breaks down — likely because
  the matrix is not SPD on the duplicated equation pair, or there's a
  pre-conditioner / dot-product double-count.

### CG-state instrumentation (NEXT SESSION FIRST ACTION)

Add per-iter print of pAp, rho, alpha to `mars_cg_solver.hpp` solve() under
a verbose-with-state flag (already have `verbose_`; add another field
`stateDump_` to avoid spamming the converging case). Run:

```bash
srun --export=ALL,MARS_AMR_REUSE_DOMAIN=1 \
  --account=csstaff --time=00:03:00 --nodes=1 --ntasks-per-node=4 \
  ~/affinity/bind_numa.sh \
  /capstor/scratch/cscs/gandanie/git/mars/daint-gpu/examples/distributed/unstructured/mars_tgv \
  --mesh=/capstor/scratch/cscs/gandanie/git/mars/cube16.mesh \
  --box-lo=0 --box-hi=1 --V0=1 --nu=0.05 --dt=1e-4 --num-steps=2 \
  --pressure-solve=DDT --skew=1
```

Look at pAp values across the first few CG iters:
- If pAp is small-positive but drifting → preconditioner / dot-product
  inconsistency between ranks
- If pAp is NEGATIVE → matrix is genuinely not SPD; the cross-rank duplicate
  equation pair is the culprit; need to either demote cross-rank slaves
  (and exchange their p value via cstone halo on a halo-range DOF) OR use
  the new `spmvPostCallback` (already wired in `mars_cg_solver.hpp`) to sum
  Ap[slave] and Ap[master] via an MPI Alltoallv pair-sum.

### Quick-diagnostic toggle

Set `MARS_PERIODIC_HALO_DBG=1` in srun --export to enable the one-shot
print of (sum-owned, sum-halo, sum-halo-on-min-face) at setupNSStepper.
Confirms cstone halo behaviour on the run's mesh.

### Files touched in this session

```
backend/distributed/unstructured/domain.hpp                        [stale comment removed]
backend/distributed/unstructured/fem/mars_ns_solver.hpp            [mass kernel ownership filter; all 7 scatter sites loop [0,elementCount); diagnostic added]
backend/distributed/unstructured/solvers/mars_cg_solver.hpp        [spmvPostCallback hook added]
```

### Test commands (canonical)

Single-rank (baseline, must remain 0.158 KE):
```bash
srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=1 \
  --export=ALL,MARS_AMR_REUSE_DOMAIN=1 \
  ~/affinity/bind_numa.sh \
  /capstor/scratch/cscs/gandanie/git/mars/daint-gpu/examples/distributed/unstructured/mars_tgv \
  --mesh=/capstor/scratch/cscs/gandanie/git/mars/cube16.mesh \
  --box-lo=0 --box-hi=1 --V0=1 --nu=0.05 --dt=1e-4 --num-steps=50 \
  --pressure-solve=DDT --skew=1
```

Multi-rank (broken, target):
```bash
srun --export=ALL,MARS_AMR_REUSE_DOMAIN=1 \
  --account=csstaff --time=00:03:00 --nodes=1 --ntasks-per-node=4 \
  ~/affinity/bind_numa.sh \
  /capstor/scratch/cscs/gandanie/git/mars/daint-gpu/examples/distributed/unstructured/mars_tgv \
  --mesh=/capstor/scratch/cscs/gandanie/git/mars/cube16.mesh \
  --box-lo=0 --box-hi=1 --V0=1 --nu=0.05 --dt=1e-4 --num-steps=50 \
  --pressure-solve=DDT --skew=1
```

---

## 2026-05-29 update — cross-rank MPI exchange landed; CG breakdown root cause isolated

Setup works on cube16/4-rank periodic: `[periodic-xr] peers=3 recv=323`, `[DBG mass] uniform 2.44e-04`, `Step 0 KE = 1.582e-01` (matches single-rank exactly).

Velocity CG still diverges at step 1. `MARS_PERIODIC_XR_SYMCHECK=1` proves why:
- `A[slave_row, master_ghost_col] = 0` on slave-owner ranks
- `A[master_row, slave_ghost_col] != 0` on master-owner ranks

Periodic-image halo elements on the slave-owner rank have 8 corners that are ALL ghosts (rank D's owned masters). NONE is rank A's owned slave node. So the sparsity builder never emits the `slave_dof → master_ghost_dof` edge. The matrix is structurally non-symmetric on the seam.

Next session: pick one of
- **(A) Sparsity surgery** — MPI-fetch master rows, inject matching values into slave rows. Architecturally correct, non-trivial.
- **(B) Dirichlet-style slave row + post-solve patch** — `A[slave,:]=0, A[slave,slave]=1`, `b[slave]=0`, then `x[slave] = x[master]` via MPI after solve.
- **(C) spmvPostCallback zeros `Ap[slave]` each iter + final exchange.** Lightest, no matrix surgery.

Infrastructure for (B)/(C) is in this commit: `crossRankPeriodicPairSum[Dof]` in `mars_periodic_bc.hpp`, `setSpmvPostCallback` in `mars_cg_solver.hpp`.

Test command for next session:
```bash
srun --export=ALL,MARS_AMR_REUSE_DOMAIN=1,MARS_PERIODIC_HALO_DBG=1,MARS_CG_TRACE=1,MARS_PERIODIC_XR_SYMCHECK=1 \
  --account=csstaff --time=00:03:00 --nodes=1 --ntasks-per-node=4 ~/affinity/bind_numa.sh \
  /capstor/scratch/cscs/gandanie/git/mars/daint-gpu/examples/distributed/unstructured/mars_tgv \
  --mesh=/capstor/scratch/cscs/gandanie/git/mars/cube16.mesh \
  --box-lo=0 --box-hi=1 --V0=1 --nu=0.05 --dt=1e-4 --num-steps=50 \
  --pressure-solve=DDT --skew=1 2>&1 | grep -E "DBG mass|cg-trace|^Step|periodic-xr|xr-symcheck|cg_uvw"
```

---

## 2026-05-31 update — cross-rank DOF collapse implemented; one final fix left

Committed: ec41ec9. Architectural cross-rank DOF collapse landed. Setup runs clean, no NaN/OOB. Velocity CG still returns -2 because master row's off-diagonal stiffness is ~30% of expected (avel-probe rank 0: nnzOff=5 vs expected 17, sumAbsOff=0.009 vs 0.04).

Root cause: rank D's assembler skips the periodic-image-halo elements from rank A because none of their corners are rank-D-owned. With Pass-1 redirect, slave_ghost_on_D has nodeToDof=master_owned_dof_on_D, but the gate `if (ownership[i] == 0) continue` doesn't see that.

**Fix for next session (one-line change in two files):**

In `backend/distributed/unstructured/fem/mars_cvfem_hex_kernel.hpp:464` and `backend/distributed/unstructured/fem/mars_cvfem_hex_kernel_tensor.hpp:249`, relax the row-corner ownership gate from:
```cpp
if (ownership == 0 || row_dof < 0 || row_dof >= numOwnedRows) continue;
```
to:
```cpp
// Accept ownership==0 corners whose nodeToDof was redirected (via periodic
// DOF collapse) to a LOCAL OWNED DOF. Lets rank D's periodic-image halo
// elements scatter slave-side stiffness into the master row.
if ((ownership == 0 && (row_dof < 0 || row_dof >= numOwnedRows))
    || row_dof < 0
    || row_dof >= numOwnedRows) continue;
// Equivalent shorter form:
// if (row_dof < 0 || row_dof >= numOwnedRows) continue;
```
The shorter form drops the `ownership` check entirely: any corner whose dof is a valid owned row is allowed through. Cross-rank slave_ghosts now scatter into the master row; ghost corners that don't map to a local owned DOF still get filtered by the `row_dof < 0 || row_dof >= numOwnedRows` bound.

Test:
```bash
make -j40 mars_tgv
srun --export=ALL,MARS_AMR_REUSE_DOMAIN=1,MARS_PERIODIC_AVEL_DBG=1,MARS_CG_TRACE=1 \
  --account=csstaff --time=00:03:00 --nodes=1 --ntasks-per-node=4 ~/affinity/bind_numa.sh \
  /capstor/scratch/cscs/gandanie/git/mars/daint-gpu/examples/distributed/unstructured/mars_tgv \
  --mesh=/capstor/scratch/cscs/gandanie/git/mars/cube16.mesh \
  --box-lo=0 --box-hi=1 --V0=1 --nu=0.05 --dt=1e-4 --num-steps=50 \
  --pressure-solve=DDT --skew=1
```
Success: `[avel-probe rank 0]` shows `nnzOff ~ 17, sumAbsOff ~ 0.04`. CG converges with `cg_uvw=3-5/step`. KE decays smoothly from 0.158.

Per-rank periodic DOF collapse counts (these are correct and verified):
- rank 0: 81 → 1296 active
- rank 1: 216 → 1008
- rank 2: 344 → 952
- rank 3: 176 → 840
- total: 4096 (= 16³, truly unique)

Diagnostics available via env vars (all one-shot, gated):
`MARS_PERIODIC_HALO_DBG`, `MARS_PERIODIC_XR_CHECK`, `MARS_PERIODIC_XR_SYMCHECK`, `MARS_PERIODIC_PATHB_DBG`, `MARS_PERIODIC_AVEL_DBG`, `MARS_CG_TRACE`.
