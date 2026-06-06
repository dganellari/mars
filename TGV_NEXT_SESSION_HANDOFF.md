# TGV Multi-Rank Periodic — Handoff for the Next Session

This is a focused handoff for whoever continues the **multi-rank periodic TGV** work. The main
TGV_HANDOFF.md has the full blow-by-blow; this is the actionable state + the one thing left to verify.

## TL;DR — where it stands (2026-06-06)

Multi-rank periodic TGV (4-rank cube16) went from **instant blowup → believed FIXED**, pending one
verification run. The whole problem was making the periodic seam cross-rank-consistent. Four distinct
bugs were found and fixed, each grounded in a confirmed code defect (not guesses):

1. **True cross-rank DOF collapse** (`612faa2`) — the architecture. A cross-rank periodic slave no
   longer owns a DOF; it migrates onto its master's ghost DOF, so the pair is ONE global DOF.
   PROVEN by `sum(numOwnedDofs)=4096` (= single-rank unique count; was 4657). Deleted ~500 lines of
   emulation. This is the MFEM/deal.II-style design the user demanded.
2. **Corrector D^T adjoint fix** (`75b454a`) — the reduced corrector gradient was missing the
   cross-rank fold (`maybePeriodicSum`) + broadcast the OPERATOR does, so M^-1 D^T(corrector) !=
   M^-1 D^T(operator) at the seam -> projection leaked div each step. Dropped div@step40 2.47 -> 0.38.
3. **Telescoping fold fix** (`3eb6dee`) — `maybePeriodicSum` ran on the FLATTENED partner in one
   unordered pass, so an edge/corner node (slave of a local intermediate that is itself cross-rank)
   didn't telescope -> a fixed seam node mis-folded the same way every step (the monotone one-rank
   `[advN-sum]` drift). Fix: keep a DIRECT parent table (`d_periodicPartnerDirect`), fold on it with a
   3-hop loop. Added `MARS_PERIODIC_FOLDCOUNT` probe.
4. **Cross-rank send ultimate-master routing** (`e8b0c3a`) — #3 tripped the corner-chain `MPI_Abort(73)`.
   Root: `crossRankPeriodicPairSum` is a SINGLE MPI round. Fix: GATE the send on the direct parent but
   ROUTE/KEY by the FLATTENED ULTIMATE master (terminal, owned, never a slave by flatten construction),
   so each slave lands directly on its owned ultimate master in ONE round. This also fixed a LATENT
   SILENT-DROP the abort had masked. Independently confirmed correct by workflow (high confidence).

## LATEST RUN RESULT (2026-06-06, e8b0c3a built) — fix is CLOSE but NOT done: a DOUBLE-FOLD remains

The 1-step MARS_PERIODIC_FOLDCOUNT run printed:
```
[periodic-foldcount] global_master_sum=5464 expected(owned_nodes)=4913 defect=551 max_master=13 nonzero_slave_residual=0
```
GOOD NEWS: no MPI_Abort (the e8b0c3a cross-rank routing fix worked), numOwnedDofs=4096 still correct,
nonzero_slave_residual=0 (every slave slot was zeroed -> nothing DROPPED). The single step ran:
[div-global] signed_sum=-2e-15 (rank-invariant, ~0 -- GOOD), but div_max=0.94 after one corrector
(div-split ratio 784 -- still a big seam divergence).

THE REMAINING BUG: defect = +551 (POSITIVE) and max_master=13. A telescoping fold gives each ultimate
master EXACTLY its group_size and max_master <= 8 (a cube corner has 8 images). max_master=13 > 8 and
positive defect => some masters are DOUBLE-COUNTED (a node's contribution lands on a master MORE than
once). residual=0 rules out DROPPED; this is pure OVER-accumulation.

LIKELY CAUSE (task #1 for next session): the LOCAL 3-hop fold (maybePeriodicSum stage-a,
periodicPairSumKernel x3 on d_periodicPartnerDirect) and the CROSS-RANK send (stage-b, gate on direct
parent, route to ultimate master) are NOT mutually exclusive for some edge/corner node -> it is folded
locally AND sent cross-rank, OR the 3-hop loop folds a node onto an intermediate that ALSO receives it
via another hop. Note: the local fold ZEROES the slave after folding, so a node folded locally then
"sent" cross-rank would send 0 (harmless) -- UNLESS the send reads a slot the local fold did NOT zero
(because stage-a's gate own[directParent]==1 skipped it, but stage-b's gate own[directParent]!=1 ALSO
skipped it -> neither, OR a different node double-routes). The +551 and max_master=13 are a FIXED,
reproducible miscount -> use MARS_PERIODIC_FOLDCOUNT to drive it to defect=0.
DEBUG PLAN: instrument WHICH nodes have field>group_size after the fold (dump their mask popcount +
direct vs ultimate partner) to see if it's edge (2-bit), corner (3-bit), or the same-rank vs cross-rank
boundary. The defect 551 ~ the count of multi-bit seam nodes double-counted. Fix the stage-a/stage-b
partition so each owned slave's contribution reaches its ultimate master EXACTLY once.

## THE ONE THING TO DO FIRST: verify the fix (it has NOT been run yet)

The last code change (`e8b0c3a`) is committed + rsync'd to Alps but the user has NOT run it. Rebuild
on Alps and run:

```
# 1-step foldcount check (the decisive topology gate):
srun ... mars_tgv --mesh=.../cube16.mesh --box-lo=0 --box-hi=1 --V0=1 --nu=0.05 --dt=1e-4 \
    --num-steps=1 --pressure-solve=DDT --skew=1 --adapt-every=0
#   with MARS_PERIODIC_FOLDCOUNT=1 in --export
```
PASS = NO `MPI_Abort(73)`, and `[periodic-foldcount] ... defect=0 ... nonzero_slave_residual=0`.
Then the 110-step run (drop FOLDCOUNT): `[advN-sum]` monotone one-rank drift GONE, `div(u_n)` flat
(was growing 1.7e-3 -> 1.3e-2 -> blow), survives 110 steps, KE matches single-rank.

If `defect=0` and it survives 110 steps -> multi-rank periodic TGV is DONE (DOF collapse + projection
adjoint + fold all correct). If `defect != 0`, the foldcount output says how many nodes and which
direction (negative=dropped, positive=double) -> refine the fold from there.

## Diagnostic gates available (all env-gated, off by default)
- `MARS_PERIODIC_FOLDCOUNT=1` — ones-vector through maybePeriodicSum; defect==0 => fold telescopes.
- `MARS_PERIODIC_EDGE_AUDIT=1` — FRAGMENTED==0 => every edge/corner slave reaches its all-min master.
- `MARS_PROJ_PROBE` — P1 |A*phi-b|/|b|, P2/P3 div ratios, seam|u_s-u_m|.
- `[div-global]` (auto in debug window) — signed_sum ~0 & rank-invariant; signed drift => fold issue.
- SYMPROBE (two applyDDTReduced + cross-dot) rel ~1e-12 => A=A^T (the operator stayed symmetric).
- `[advN-sum]` — per-rank owned advection sum; a monotone one-rank drift = the fold-telescoping bug (now fixed).

## Things RULED OUT (do not re-chase)
- Hypre for periodic: blows up (velocity solve can't do per-matvec seam merge; now the collapse may have
  unblocked it but unverified). CG matrix-free is the working path.
- Seam VELOCITY drift (two fixes refuted by runs: bea6d38 reverted, advN broadcast deleted as inert).
- Area-vector geometry (suspect D): refuted — no element straddles the seam, halo area vectors not read
  by the owned-local operator loop.
- Skew energy term: refuted by --skew=0 (upwind blows up identically -> scheme-independent).

## Key files / line anchors (verify against current code before trusting)
- `backend/distributed/unstructured/fem/mars_ns_solver.hpp`: DOF collapse ~2700; maybePeriodicSum
  ~6510 (3-hop direct-parent fold + cross-rank); corrector reduced branch ~7840; [div-global] probe
  after ~8345; FOLDCOUNT probe after the pressure-K assembly lap.
- `backend/distributed/unstructured/fem/mars_periodic_bc.hpp`: PeriodicMap struct (d_periodicPartner +
  d_periodicPartnerDirect) ~104; direct snapshot before flatten ~974; buildCrossRankPeriodicMap STEP-1
  (gate on direct, route by ultimate) ~699; periodicPairSumKernel ~300.
- Memory: `project_multirank_periodic_owner_migration` (the architecture, supersedes the old emulation
  notes).

## Workflow constraints (from CLAUDE.md + memory)
- Commits SHORT, no Co-Authored-By. rsync local->Alps only (I push edits); user runs build + srun.
- Both `--solver=cg` and `--solver=hypre` paths must stay (user requirement).
- cube16 is PUBLIC; the pump mesh is NDA — never read it.
