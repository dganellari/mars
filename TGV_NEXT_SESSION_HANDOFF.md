# TGV Multi-Rank Periodic — Handoff (2026-06-07)

**Read this before touching code.** Three sessions, ~50 fix attempts, ~25 reverted. The bug is NOT closed but is precisely characterized.

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
