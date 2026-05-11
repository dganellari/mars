# cstone Domain Reuse: Implementation Status

## Where we are (as of 2026-05-10)

Targeting a ~80% reduction in cstone sync time during AMR rebuilds (focusTree.converge dominates at ~2.9 s of 3.6 s at 977M elements / 16 ranks). Domain-reuse strategy: keep the cstone Domain object alive across AMR refinements so `firstCall_` stays `false` and `converge()` is skipped.

## What's implemented

**MARS-side** (`backend/distributed/unstructured/`):
- `ElementDomain::resyncFromDevice(coords, conn)` — replaces internal mesh state, resets all lazy-cached sub-objects, calls `cstone::Domain::resetForAMRRefinement` then `sync()`.
- `mars_amr_manager.hpp::rebuildDomainFromDevice` — branches on `MARS_AMR_REUSE_DOMAIN=1` env. Default off; reuse opt-in.

**cstone-side** (FetchContent patches via `scripts/patch_cstone_amr_reset.py`):
- `Domain::resetForAMRRefinement(LocalIndex newCount)`: resets `bufDesc_`, `prevBufDesc_`, `layout_`; calls focus tree reset.
- `FocusedOctree::resetFocusRangeForNewDistribution()`: resets `prevFocusStart`, `prevFocusEnd`, sets `rebalanceStatus_ = valid`.

## Iteration history

| Attempt | Symptom |
|---|---|
| Just call `sync()` again with new arrays | Wrong output: stale `bufDesc_.size` ignores new mesh |
| Add `setEndIndex(newCount)` | SEGFAULT inside `focusTree.updateTree` (rebalanceStatus_ check throws) |
| Apply explorer's full reset (`status=invalid`) | Throws "update of criteria required" — `invalid` doesn't pass `updateTree`'s precondition |
| Change `status=valid` instead | Runs to completion. **L2 norm drifts 9% from default (7.67 vs 7.01).** AMR marks 32768/32768 elements vs 32749/32768 in default. |

**Current state:** functionally correct enough to run, but produces a different mesh discretization than the firstCall_ path. Not shippable.

## Suspected remaining cause

Focus tree internal state (centers, geometric centers, leaf counts, possibly `octreeAcc_` itself) is being reused across AMR rebuilds. The next sync's `updateMinMac/updateTree/updateCounts/updateGeoCenters` calls modify these in-place starting from prior state — so the result depends on history. firstCall_'s `converge()` loop drives this to a fixed point (history-independent). Domain reuse skips that fixed-point step → non-deterministic w.r.t. construction history.

## Open question for upstream

Is "AMR-style sync" with discontinuous particle-count changes meaningful in cstone's design? Or is `firstCall_` an architectural assumption that cannot be retrofitted around?

If the latter, cstone needs a new mechanism: probably a `Domain::amrSync()` that does what `sync()` does but explicitly re-converges from scratch on a known stale focus tree. That's a real upstream design discussion, not a 25-line patch.

## Disposition

**Domain reuse remains opt-in via env flag, default off.** Do not flip default until L2 norm matches bit-exact across paths.

Next step: Explorer agent #2 mapping which focus tree fields incrementally use prior state. Awaiting report.
