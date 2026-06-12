# TGV Multi-Rank Periodic — Handoff (2026-06-07, BREAKTHROUGH UPDATE)

**Read this before touching code.** The primary feedback loop is FOUND and FIXED. The bug is reframed from "projection won't close" to "two feedback loops amplify the seam residual"; one is fixed, one is precisely localized.

## BREAKTHROUGH (this session)

The reframe that cracked it: stop chasing PROJ-P3 → 0 (the channel never closes it either — see docs/poiseuille_tutorial.md Part 6). The real failure is the **growing div_max feedback loop**, not the nonzero projection residual. Two loops:

**PRIMARY LOOP — FIXED.** The predictor's `grad(p)` used the FOLDED periodic seam reduction (`maybePeriodicSum` at NS:7129-7131) while the corrector/operator under H2 use PER-SLOT (no fold). The mismatch `(G_fold − G_perslot)·p` is LINEAR in accumulated p → dt-independent geometric feedback gain ≈ 1.17 (which is why ~50 projection-side fixes never moved it — the loop is in the predictor↔corrector operator mismatch, not the projection). Fix: `MARS_PRED_PERSLOT_GRADP=1` makes the predictor `grad(p)` per-slot too, matching the corrector. Gain drops 1.17 → ~1.05. (Fable diagnosed this; Opus missed it for two days.)

**SECONDARY LOOP — localized, not yet fixed.** After the primary fix, a weaker loop (gain ~1.05) remains: the skew-symmetric advection's `mdot` at the cross-rank periodic face. Skew (Verstappen) advection conserves discrete KE only when `mdot` telescopes across shared faces; at the periodic seam the master-side and slave-image elements compute `mdot` from slightly-inconsistent `u`, so it does NOT cancel and injects KE each step. PROOF: `--skew=0` (upwind) stays BOUNDED (KE oscillates 200-960, CG healthy through step 200) while `--skew=1` goes to Inf by step 160 — upwind's numerical diffusion damps the injection. CONSTRAINT: the fix is NOT a velocity broadcast — `MARS_TGV_USTART_BCAST=1` (force u[slave]:=u[master] at step start to make mdot consistent) made it WORSE (step 140 vs 160) because broadcasting u re-injects divergence the projection didn't account for. The fix must reconcile `mdot` across the periodic face WITHOUT changing stored `u`.

**Empirical progress:** baseline blew at step ~40 (KE → 1e158). With `MARS_PRED_PERSLOT_GRADP=1 MARS_ROTATIONAL_P=1`: KE decays correctly through step 100, div_max bounded to ~16 at step 100, survives to step ~160. **4× longer survival, clean physics 2.5× longer.** Not stable yet, but a qualitative change.

## Best-known config (env flags, all default OFF, all pump-safe)

```
MARS_PRED_PERSLOT_GRADP=1   # PRIMARY FIX: predictor grad(p) per-slot at seam (de7e366)
MARS_ROTATIONAL_P=1         # damps residual pressure accumulation (Timmermans)
```

Tested combinations on cube16 4-rank (step where KE turns up / blows):
| Config | Blow-up step |
|--------|-------------|
| baseline | ~40 |
| PRED_PERSLOT only | ~110 |
| PRED_PERSLOT + ROTATIONAL | **~160 (best)** |
| PRED_PERSLOT + NONINCREMENTAL | ~110 (non-inc made it worse — loop is NOT pressure-accumulation) |
| PRED_PERSLOT + ROTATIONAL + H2=0 | ~150 (H2 broadcast helps slightly) |
| PRED_PERSLOT + ROTATIONAL + USTART_BCAST | ~140 (worse — velocity broadcast re-injects div) |
| `--skew=0` upwind + PRED_PERSLOT + ROTATIONAL | BOUNDED to step 200 (numerical diffusion damps secondary loop) |

## NEXT SESSION: the secondary-loop fix

The secondary loop is in the skew advection `mdot` at the periodic seam. Read the skew kernel (NS:914-936): per face, `mdot` is computed once and used at both L and R, so per-face KE-conservation is fine WITHIN an element. The leak is that the periodic-image element (slave-owner rank) and the master-side element compute DIFFERENT `mdot` for the SAME physical periodic face (from inconsistent `u`), breaking the cross-element telescoping.

Candidate fixes (untested):
1. **Reconcile mdot across the periodic face** — after the per-element mdot is computed but before the flux scatter, average/fold the mdot at the periodic seam so master and slave-image faces use the SAME value. This is the mdot-level analog of the grad(p) per-slot fix. Does NOT touch u (avoids the USTART re-injection problem).
2. **Use the Verstappen rotational advection form explicitly** — `N_skew = N_div − ½ q (div u)` computed so the `div u` term uses the seam-folded divergence. The pump solver mentions this form (NS:961-984).
3. **Blend upwind dissipation at seam nodes only** — keep skew interior, add first-order upwind within one element of the periodic seam. Pragmatic; `--skew=0` proves it would bound the loop.

The agents (Fable + Opus workflows) hit the session token limit at 4:20pm Zurich; relaunch the secondary-loop diagnosis workflow then if option 1's mdot-reconcile needs the algebra verified.

---
## ORIGINAL HANDOFF (pre-breakthrough, kept for reference)

## Current state (HEAD = `25734c7`)

- Single-rank periodic TGV: **works perfectly** (PROJ-P3 ≈ 1e-14, KE decays cleanly, no blow-up)
- Multi-rank periodic TGV (4-rank cube16): **blows up around step 35-45** in 110-step run
- All other BC modes (Pump, Cavity, Channel) are **untouched** and confirmed working
- Probe suite is intact and informative

## The bug in one line

In multi-rank periodic CG, the pressure-correction projection identity `D(u^{n+1}) = 0` does NOT close at master DOFs. PROJ-P3 ratio stays at 0.4-1.0 (should be ~1e-14). Divergence accumulates step-over-step until catastrophic blow-up.

## Bisection of failure modes

The session has eliminated many wrong directions. **Use this map** to avoid relitigating:

### KNOWN GOOD (verified by probes, don't touch)

- **Foldcount probe** (888045b, 1696617, e1a27f6): on cube16 4-rank, `defect=0, max_master=8, buckets={3375,675,45,1}` matching single-rank exactly. Telescoping partner table is bit-correct.
- **Mass mirror** (NS:3115-3123): `m_M = m_S = m_merged = m_M_geom + m_S_geom` post-mirror. Verified by `XRANK-SLAVE-MASS` probe (0 zero-mass slaves; max-mass on slave = exactly `m_geom_M + m_geom_S`).
- **Area vectors** (mass + AREAVEC probes): owner and halo deliver bit-identical per-element area vectors. `SEAMAX` probe confirms inter-element seam orientation correct.
- **Operator SPD** (SYMPROBE rel = 0): the bare reduced operator `A_red = D M^-1 D^T` with `applyPeriodic=false, reducedPeriodicFold=true` is self-adjoint to roundoff.
- **CG convergence** (PROJ-P1 = 1e-10): linear solve inverts `A_red φ = b` correctly.
- **Phi at periodic pairs** (PHI_SEAM probe): `s.d_phi[slave] = s.d_phi[master]` bit-equal after post-CG broadcast.
- **Pump path untouched**: all session edits gate on `BCKind::Periodic && periodicMap != nullptr`. Pump uses `BCKind::Pump`, `periodicMap = nullptr`, takes the else-branch unchanged. The cross-rank-slave write fix (`e949daa`) is dead code on pump (`dof < numOwnedDofs` for all owned pump nodes).

### KNOWN BROKEN (DO NOT RETRY)

- **`MARS_OP_GBCAST_REDUCED=1`** (operator master→slave g-broadcast inside `applyDDTPerNode`): breaks SYMPROBE (rel = 0.47), CG diverges. The g-broadcast at NS:5395 is the documented adjointness-breaker; the in-file comment at NS:5391-5394 is correct.
- **BDF2 fix at NS:7536**: `invDt = 3/(2*dt)` under BDF2 (instead of `1/dt`) caused immediate catastrophic blow-up (`|Du**|` grew 7e-3 → 1e97 over 30 steps). The workflow's algebra was confident but empirically wrong. Reverted in `25734c7`. Don't repeat without first understanding WHY the algebra failed.
- **Routing multi-rank periodic CG corrector to BARE reduced branch** (without FOLD_G): PROJ-P3 = 0.999 every step. The reduced-CG corrector with bare gradient does nothing measurable to the divergence.
- **`MARS_CORR_NO_GFOLD=1`**: PROJ-P3 went from 0.508 to 0.828 (worse). Per-slot g without fold is wrong.
- **Sign-flip g[S] before broadcast** (workflow H1): never wired but algebra has open questions (sign-mask at triple-image corners is non-trivial).

### THE PRECISE EMPIRICAL SIGNATURE

OP_GRADDUMP across cube8/cube16/cube32/cube64 shows: at periodic-pair nodes, the operator's intermediate `g` has:
- **Tangential components**: bit-identical between master and slave slots
- **Periodic-axis-normal component**: NOT symmetric; ratio g[S]/g[M] varies 0.6 to 2.0 across pairs (sometimes opposite-signed, sometimes same-sign with different magnitude)

8 of 9 probed pairs across cube8/cube32 showed **opposite-sign in periodic-axis-normal**. cube64 had mixed signs. cube16 had 0.65× same-sign in many cases.

This is the structural asymmetry that the corrector can't undo with any simple per-slot manipulation.

## Two open strong diagnoses (both confirmed by code-grep, neither fixable by me)

### Diagnosis 1: BDF2 pressure-RHS mismatch (workflow wtb9updx2, confirmed at NS:7536 + NS:8266-8267)

`runPressureSolveStep` uses `invDt = 1/dt` (BDF1) unconditionally. `runImplicitDiffusionStep` correctly uses `invDt = 3/(2·dt)` under BDF2. `runCorrectorStep` uses `dtEff = 2·dt/3` under BDF2.

Algebra says: corrector subtracts `dt/dtEff = 3/2` too little under BDF2. PROJ-P3 should be `2/3 ≈ 0.67` constant after step 2.

**But the empirical PROJ-P3 oscillates 0.6-1.0, NOT constant 0.67. And patching `invDt = 3/(2·dt)` made things catastrophically WORSE, not better.**

So the algebra is structurally wrong somewhere — possibly the corrector ALSO has a hidden 3/2 elsewhere that already compensates, OR the relationship between `invDt`, `b`, and `dtEff` is not the simple proportionality the workflow assumed.

Reference comparison point: `runImplicitDiffusionStep` at NS:7424-7426 (works under BDF2). Compare its pipeline to `runPressureSolveStep`'s — find what's different in scaling.

### Diagnosis 2: H2 broadcast destroys slot-wise field once eddies advect (workflow wtb9updx2 diag 2)

With `MARS_TGV_H2=1` (currently default ON in `25734c7`), the pre-probe master→slave broadcast at NS:8743-8761 overwrites `u^{n+1}[slave]`. The first 1-2 steps benefit from this (clean periodic anti-symmetry). After eddies advect, the anti-symmetry breaks, broadcast destroys the slot-wise projection-identity field, PROJ-P3 ratio jumps to ~1.0.

Falsifier (untested due to BDF2 attempt): run with `MARS_TGV_H2=0`. If PROJ-P3 collapses to stable ~0.33 every step from step 2 onward, H2 broadcast IS the issue. If it still oscillates, H2 is not the cause.

## What every commit kept (HEAD = `25734c7` carries all of these)

| Commit | Purpose | Status |
|--------|---------|--------|
| 612faa2 | Owner-migration cross-rank DOF collapse (`sum(numOwnedDofs)=4096`) | ✓ Required |
| 1696617 | Race-free single-pass on flattened partner table | ✓ Required |
| 888045b | `periodicPairSumKernel` owned-slave gate (kills double-fold post-halo) | ✓ Required |
| e1a27f6 | Revert `periodicFoldToMasterKernel` gate (pre-halo sites need ungated) | ✓ Required |
| 808140e | Defer u-broadcast past PROJ-P3 probe (probe now reads truthful slot-wise field) | ✓ Required |
| e949daa | `applyCorrectorPerNodeKernel` writes cross-rank slaves (`dof >= numOwnedDofs`) | ✓ Required |
| 1109de9 + ccf3ccd | Falsifier probes: SEAMAX, XRANK-SLAVE-MASS, PHI_SEAM, GRADSLOT, GRADDUMP, GRADTRACE, OP_GRADDUMP, WHEREMAX | ✓ All work |
| 999d4e7 | PROJ-FACTOR probe (A_periodic vs A_bare on converged phi) | ✓ Works |
| b51c5a0 | 12-digit precision + `|dOut-dIn|` explicit print | ✓ Required |
| a79e721 | MARS_TGV_H2 default ON (best-step-1 config) | ⚠️ Partial — see Diag 2 |

## Probe inventory (all env-gated, default OFF; all triple-gated for pump safety)

```
MARS_PROJ_PROBE=1            # PROJ-P1, PROJ-P2, PROJ-P3, dOut-dIn
MARS_PROJ_FACTOR_PROBE=1     # A_periodic vs A_bare on converged phi
MARS_OP_GRADDUMP=1           # operator's internal g at first cross-rank slave
MARS_GRADDUMP_PROBE=1        # corrector's gradPhi at first cross-rank slave
MARS_GRADTRACE_PROBE=1       # gAcc snapshot pre/post-reverseExchange
MARS_GRADSLOT_PROBE=1        # rel |gS-gM| over all cross-rank slaves
MARS_PERIODIC_AREAVEC_PROBE=1   # per-element Stokes closure (sum of 12 SCS)
MARS_PERIODIC_SEAMAX_PROBE=1    # inter-element seam mirror orientation
MARS_XRANK_SLAVE_MASS_PROBE=1   # cross-rank slave mass mirror check
MARS_PHI_SEAM_PROBE=1        # |phi[slave] - phi[master]| split xr/sr
MARS_WHEREMAX_PROBE=1        # locate the global-max |Du^{n+1}| node
MARS_DDT_SYMPROBE=1          # SYMPROBE for A_red self-adjointness
MARS_PERIODIC_FOLDCOUNT=1    # buckets eq{1,2,4,8} for partner table
MARS_TGV_H2=0                # opt out of pre-probe u-broadcast
MARS_USE_REDUCED_CORR=1      # opt into reduced-CG corrector branch (worse)
MARS_USE_ELSE_CORR=1         # opt back to else-branch corrector
MARS_CORR_FOLD_G=0           # opt out of corrector fold+broadcast
MARS_CORR_NO_GFOLD=1         # disable corrector fold (different from CORR_FOLD_G=0)
MARS_OP_GBCAST_REDUCED=1     # operator g-broadcast (BROKEN — breaks SPD)
```

## Reproduce the broken state

```bash
# On Alps (cube16, 4-rank, 110-step blow-up around step 40):
MARS_PROJ_PROBE=1 \
srun --export=ALL,MARS_AMR_REUSE_DOMAIN=1,MARS_PROJ_PROBE=1 \
  --account=csstaff -p debug --nodes=1 --ntasks-per-node=4 \
  ~/affinity/bind_numa.sh \
  /capstor/scratch/cscs/gandanie/git/mars/daint-gpu/examples/distributed/unstructured/mars_tgv \
  --mesh=/capstor/scratch/cscs/gandanie/git/mars/cube16.mesh \
  --box-lo=0 --box-hi=1 --V0=1 --nu=0.05 --dt=1e-4 \
  --num-steps=110 --pressure-solve=DDT --skew=1 --adapt-every=0 \
  2>&1 | grep -E "^Step|PROJ-P3"
```

Mesh sizes available (or generate via `scripts/generate_hex_cube.py`):
- `cube8.mesh` (smallest, fastest)
- `cube16.mesh` (canonical)
- `cube32.mesh` (CG breaks: SYMPROBE rel = 0.092 — A_red non-SPD on this partition geometry; separate bug)
- `cube64.mesh` (CG converges; PROJ-P3 ~ 1.3)

Single-rank reference: same command with `--ntasks-per-node=1`. Gives PROJ-P3 ≈ 1e-14, no blow-up.

## Recommended next-session entry points

In order of likely tractability:

1. **Understand why the BDF2 fix blew up.** The algebra was confident, the code-grep was correct, both `runImplicitDiffusionStep` and `runCorrectorStep` use the BDF2 pattern correctly. Yet patching `runPressureSolveStep`'s `invDt = 3/(2·dt)` caused immediate catastrophic blow-up. There's a hidden assumption about pressure-RHS scaling that I missed. Find it; the answer is probably one line.

2. **Test `MARS_TGV_H2=0` for the full 110 steps.** Diag 2 predicted PROJ-P3 collapses to stable ~0.33 every step. I didn't actually run this empirically. If it works, ship H2=0 as default. If not, eliminate Diag 2.

3. **Investigate the cube32 SYMPROBE failure** (rel = 0.092). This is structural in the operator at that mesh size. The operator should be self-adjoint regardless of mesh — if it's not on cube32 specifically, there's a real bug in cstone's owner-migration partitioning at that geometry. Likely a separate bug from the main TGV failure but might illuminate it.

4. **Re-derive A_red φ = b symbolically on paper** before any more code attempts. The session has burned ~20 fix attempts that all had plausible algebra but wrong empirics. The discrete operator structure has subtleties this session has not modeled. Particularly: the relationship between the operator's `g` (per-slot, no fold) and the projection identity's expected `gradPhi` (which the corrector subtracts). The session's algebra has consistently predicted PROJ-P3 = 2 or 3 or 0.5 when reality was 0.5-1.0 oscillating.

5. **Sanity check**: compare to MFEM/deal.II reference implementation of periodic Chorin. They have working multi-rank periodic NS. Their corrector/pressure pipeline structure might reveal what MARS is doing differently.

6. **If TGV multi-rank periodic is a hard demo requirement**: cavity and channel work multi-rank today. Single-rank periodic TGV also works. Switching demo target is the safest path to a working deliverable.

## What I cannot do from here

I've worked this for 2+ days, ~50 fix attempts. My algebra keeps predicting wrong outcomes. The pattern is consistent: I derive what seems like the algebraically correct fix, ship it, it makes things worse or doesn't move PROJ-P3 enough. The probe data is genuinely rich but my model of the discrete operator at cross-rank periodic seams is incomplete.

A fresh session (different model, or me with full reset) might find the answer faster than I would by iterating more.

---

*Compiled 2026-06-07. HEAD = `25734c7`. All keep-commits + probe suite intact. Pump-safe. Cavity/channel safe. Single-rank periodic safe.*
