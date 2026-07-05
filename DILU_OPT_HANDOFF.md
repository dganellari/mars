# DILU Gordon-Bell optimization — handoff (autonomous session 2026-06-22)

> **SUPERSEDED (2026-07-05).** This handoff describes a mid-arc state. The work completed: GPU-native
> setup ~19×, DILU default + hybrid smoother (`MARS_ACM_DILU_LEVELS=2`), factor-once coarsest, and the
> full **distributed ACM** (multi-rank validated at 1/2/4 GPUs, identical physics, +14%/+40% iters).
> Current state + knobs: `internal-notes/coupled_solver_progress.md` §8–9 and
> `docs/reference/06_mars_solvers.md` §11. Notes below are kept for history; several conclusions in
> them (e.g. the launch-bound V-cycle hypothesis) were later corrected by profiling.

Worked autonomously while you were AFK. Goal: optimize the GPU multicolor block-DILU smoother
(`MARS_ACM_SMOOTHER=dilu`) toward Gordon-Bell / beating AMGX. Could **not compile locally** (Mac), so
every step is **committed separately** (build per-commit to isolate any compile issue) and
**adversarially verified** (logic + CUDA-API compile-correctness, since no local build).

## Baseline measured (pump Re=1000, `--acm-rebuild=1`, before this session)
- V-cycle ~26 ms (≈700 kernel launches/cycle → ~22 ms is pure launch overhead, ~4 ms compute floor)
- setup 2.7 s/Picard, host-serial (std::set block-CSR + coloring + 72M-nnz extract + factor)
- ACM iters 35→291 (grows with convection; DILU is the weaker smoother)
- **DILU MUST run `--acm-rebuild=1`** — reuse degrades (a stale hierarchy sends a weak DILU to the 1980 cap)

## Committed this session (oldest → newest)
| commit | step | what | risk |
|---|---|---|---|
| `f162134` | 2 | fused residual `acmResidualKernel = b−Ax` (1 launch, not 2) | bit-identical |
| `e24b7b3` | 0 | stream-thread the cycle + `acmVcycleGpu` recursive→iterative | bit-identical |
| `f15c83c` | 1 | **CUDA-graph the V-cycle** (capture down/up halves, QR eager between) | verify GO; graceful eager fallback |
| `cfd571d` | 4 | **GPU-side setup** (GPU value extract + color-parallel factor) | verify GO `wy6lacp8a`; bit-identical to host |

All four committed + adversarially verified. Build per-commit if any compile issue surfaces (no local compile this session).

## VALIDATION RUN (do this first on return)
Rebuild, then the same probe as before:
```
MARS_ACM_SMOOTHER=dilu MARS_NODEHALO_V2=1 srun --export=ALL --account=csstaff --time=20 --nodes=1 \
  --ntasks-per-node=1 ~/affinity/bind_numa.sh \
  ./examples/distributed/unstructured/mars_coupled_model_operator \
  --mesh=$MESH_BIG --acm-pump --inlet-ss=$INLET_BIG --outlet-ss=$OUTLET_BIG --inlet-speed=0.5 \
  --Re=1000 --supg --relax=0.3 --picard=10 --acm-rebuild=1
```
What to look for:
1. **`[acm-graph] V-cycle captured (down+up graphs)`** on stderr → the graph engaged. If instead you see
   `[acm-graph] V-cycle capture unavailable -> eager path`, the solver is still **correct** (eager fallback)
   but didn't get the speedup — send me the surrounding CUDA error for triage.
2. **`solve=` per Picard** should drop from ~4 s toward **~1.5 s** (V-cycle 26→~9 ms × iters).
3. **`setup=` per Picard** should drop from 2.7 s toward **~1.5 s** (GPU value extract+factor; the host
   std::set + aggregation remain — that's the next win, see below).
4. **Same physical answer**: `|u|max≈610`, mass conserved, ACM iters in the same ballpark (35/88/144…).
   The graph + GPU-setup are convergence-PRESERVING — iters must match the pre-session run. If iters
   changed materially, something's wrong → tell me.

A/B against `MARS_ACM_SMOOTHER=ilu` (unchanged) is still the cross-check for correctness.

## UPDATE 2 (2026-06-22) — profile-driven coarse-solve fix (committed) + GPU-native-setup design

**nsys profile finding (corrected my guesses):** the V-cycle was NOT launch-bound nor smoother-bound. GPU does
only ~2 ms/V-cycle; the ~27 ms wall was sync/host-orchestration, **dominated by the coarsest cusolver QR
re-factoring + ~47 cudaMallocs EVERY V-cycle** (~21% GPU + the malloc/sync) on a fixed 876-DOF op. The DILU
smoother is only ~0.5 ms/V-cycle (so Step 3/coalesce is moot — good thing we profiled).

**FIX committed `4cb94e7` — factor-once dense-LU coarse solve (GPU), verify GO (w3w6nq9ii):** cusolverDn `getrf`
ONCE/Picard at setup (`acmDenseFromCsrKernel` → dense col-major) + `getrs` per V-cycle (in-place RHS,
stream-ordered, NO malloc/factor/sync). Convergence-equivalent (same direct coarse solve). **VALIDATE:** re-run
the probe — expect `[acm-coarse] dense LU n=… -> factored once`, `solve=` dropping sharply, same iters/`|u|max`.

**GPU-NATIVE SETUP design done (w2tfmtgdl)** — to kill the ~1.85 s host setup wall (the `std::set` block-CSR +
greedy coloring + sequential aggregation), "no host work" per your directive. Staged build order (design's own
risk ranking, lowest-risk first, each behind a gate):
- **§1a GPU block-CSR build** (`acmBuildBlockCsrGpu`, shared DILU+ILU): emit (I,J)+transpose+diag keys → sort →
  unique → split → brow/lower_bound/bdiagPtr. Mirrors `buildCoarseOperator`. **Convergence-PRESERVING**,
  bit-identical-testable vs the host `std::set`. Lowest risk — do first.
- **§1b GPU Jones-Plassmann coloring** on the NODE block graph (NOT cuSPARSE csrcolor — that colors the scalar
  4nN matrix and splits nodal blocks across colors, breaking block-DILU). Convergence-preserving if valid
  (every node colored, adjacent differ). ~15-25 colors vs host's 12 (a few more launches, µs).
- **§1c GPU parallel directional aggregation** — **LAST, research-grade.** Strength graph (emit vals², sort,
  reduce_by_key, symmetrize) + parallel maximal-weight matching (each free node proposes its STRONGEST free
  strong neighbor = keeps Mavriplis-directional) + grow-to-kmax + absorb-leftovers. **CRITICAL:** agg ids must
  be contiguous 0..nCoarse-1 (`buildCoarseOperator` indexes agg[] directly — gaps → silent coarse-op
  corruption; thrust::unique+lower_bound remap + assert max(agg)<nCoarse). **CHANGES convergence** (different
  aggregates than host greedy) → **loses the bit-exact oracle**; gate on (i) GPU-vs-GPU determinism, (ii) host
  V-cycle fed the device agg, (iii) FlexGMRES iter-parity vs the host-agg baseline on cavity/channel.
- §2 device-only level loop: delete the per-level D2H (vcycle.hpp:608-612); only ~2 control ints/level survive
  (nCoarse, numColors). Run the build on the default stream (preserves inv(D_k)-before-use). Graph capture is
  post-setup so unaffected.
Full design: workflow w2tfmtgdl output. Risk #2 (agg-id contiguity) and #7 (ILU IKJ factor is NOT a drop-in —
only the §1a structure ports; leave ILU factor host or confirm ILU is off the GB critical path).

## Remaining (deliberately deferred — need YOUR cluster validation, didn't want to do blind)
- **Step 3 — segment-reorder + coalesce** (V-cycle ~9→~6 ms). Reorders each block-row into
  `[lower-color | diag | upper-color]` so the solve kernels are branchless/coalesced. Lowest remaining win
  (the graph already did most of the V-cycle), and it touches the `bblk` layout that Step 4's kernels also
  use — so it wants one careful combined pass, not a blind edit.
- **Structure-reuse (the rest of the setup win, 2.7→~0.3 s)** — reuse the coloring + block-CSR pattern +
  aggregation across Picard iters, refactor only the values. This *freezes the aggregation* (a convergence
  approximation — standard frozen-coarsening AMG), so it needs an iters-vs-`--acm-rebuild=1` comparison on
  the cluster to confirm it doesn't cost too many iters. Safe to add once you can validate.
- **GPU aggregation + GPU block-CSR build** (Jones-Plassmann) — the true 10⁹-scale setup fix; only matters
  at scale where the host std::set is the wall.

## AMGX — PARKED (your call) until you confirm it's installable on Alps. Design done (`mars_amgx_solver.hpp`
wrapper + FGMRES+AGGREGATION+MULTICOLOR_DILU config); needs a from-source build on Alps + `export AMGX_DIR`.
It's the comparison yardstick, not a dependency — the DILU work stands alone.

## The contribution framing
Beating AMGX comes from two levers we have: (1) **directional Mavriplis aggregation** (anisotropy-aware,
fewer iters than AMGX's generic SIZE_8) and (2) these optimized GPU-native DILU kernels. The
`nLev=66,279 → numColors=12` (5,523× fewer serial steps) result is the headline "why naive GPU-ILU fails".
