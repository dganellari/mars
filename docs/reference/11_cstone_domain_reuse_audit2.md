# cstone Domain Reuse — Audit #2 (drift root cause)

Source: parallel agent, 2026-05. Triggered by L2-norm drift (7.67 vs 7.01)
and 133-element count divergence after the first patch landed.

## TL;DR

**The first audit's plan was incomplete.** Resetting `prevFocusStart/End` and
`bufDesc_` is necessary but not sufficient. The drift comes from cstone's
**incremental focus tree logic** (`focusTransfer`, `updateFocus`, `focusTransfer`)
which reads old `leaves_` and old `leafCountsAcc_` to decide *delta* refinement
vs *fresh build*. With AMR-style 8× particle creation, the old counts are
meaningless, but cstone uses them to decide "this leaf already has bucketSize
particles, no refinement needed" → fewer leaves get marked → different mesh
topology → different solution.

## Specific fields needing reset (added to the existing patch)

Ranked by likelihood the agent assigned:

| # | Field | Confidence | Why |
|---|---|---|---|
| 1 | `leaves_` + `leavesAcc_` | 95% | `focusTransfer()` reads them to decide incremental transitions |
| 2 | `leafCountsAcc_` + `countsAcc_` | 90% | Old per-leaf counts misguide refinement decisions |
| 3 | `macsAcc_` | 80% | `updateMinMac(..., true)` accumulates rather than replaces; old bits stick |
| 4 | `centersAcc_` | 75% | `moveCenters` uses old octree structure indirectly |
| 5 | `geoCentersAcc_` + `geoSizesAcc_` | 70% | Computed at end of `updateTree`, but only AFTER focusTransfer made decisions |
| 6 | `octreeAcc_` | 65% | `numLeafNodes` field stale; assertions fire on size mismatch |

## Proposed expanded `resetFocusRangeForNewDistribution()`

```cpp
void resetFocusRangeForNewDistribution()
{
    prevFocusStart   = 0;
    prevFocusEnd     = 0;
    rebalanceStatus_ = valid;

    // Tree structure: force fresh octree construction
    leaves_.clear();
    leaves_.push_back(0);
    leaves_.push_back(nodeRange<KeyType>(0));
    if constexpr (HaveGpu<Accelerator>{}) leavesAcc_ = leaves_;

    // Particle counts: old counts are meaningless after 8× particle creation
    leafCountsAcc_.clear();
    leafCountsAcc_.push_back(bucketSize_ + 1);
    countsAcc_ = leafCountsAcc_;

    // MAC/center state: derived from tree, must be invalidated
    macsAcc_.clear();        macsAcc_.resize(1);
    centersAcc_.clear();     centersAcc_.resize(1);
    geoCentersAcc_.clear();  geoCentersAcc_.resize(1);
    geoSizesAcc_.clear();    geoSizesAcc_.resize(1);
    globalCentersAcc_.clear(); globalCentersAcc_.resize(1);

    // Full octree: reset to single-leaf so updateFocus rebuilds
    octreeAcc_.clear();
    octreeAcc_.resize(1);
}
```

## Two algorithmic concerns the agent flagged

### Concern A: `focusTransfer()` may not handle single-leaf input

The agent's caveat (10% pessimistic case): `focusTransfer()` and `updateFocus()`
might assume a non-trivial prior tree exists. Resetting to single-leaf input
might produce undefined behavior rather than a fresh build.

If that's true, no host-side reset can achieve fresh construction. Only path
forward would be calling `converge()` explicitly for AMR reuse, defeating the
optimization.

### Concern B: `updateMinMac(..., true)` accumulates by design

`Domain::sync()` calls `updateMinMac` with `accumulate=true`. The old MAC bits
are kept and only new bits added. For Domain reuse with refresh focus tree,
this is wrong — we want fresh MACs.

Fix would be either:
- Call `updateMinMac(..., false)` once after `reset`, then do the regular sync
  (requires a flag in cstone to detect "first sync after reset")
- Or zero out `macsAcc_` so `accumulate=true` operates on a clean slate

The proposed reset zeroes `macsAcc_` — which is option (b), the simpler one.

## Test predictions (use to verify each fix is right)

| Fix | Expected | Read as |
|---|---|---|
| Reset `leaves_` + `leafCountsAcc_` only | Element count diff drops from 133 to <20 | Validates incremental-refinement hypothesis |
| Add `macsAcc_` reset | L2 drift drops from 9% to <3% | MAC accumulation was 2nd-order |
| Add `centersAcc_/geoCentersAcc_` reset | L2 drift drops to <1% | Center reuse was small but real |
| All resets | L2 norm matches default within 1e-12 | Domain reuse fully equivalent |
| Even all resets fail to drop drift below 3% | Stop. cstone's `firstCall_` is structural | Pivot to upstream collab |

## What to do

Apply the expanded reset. It's a single-file change (the patch script),
re-runnable, idempotent. Test cube16/4 first — if element count diff drops
below 20 and L2 norm drift below 1%, scale up to cube256/16 to confirm
generality. If element count diff stays large, we've hit Concern A and need
to back off.

## Calibrated honesty

The agent gave 65% / 25% / 10% probabilities for "this works / one more
iteration needed / fundamental dead end." That's **not** "very likely to
work." It's **probably** going to work. Worth one more iteration but if it
fails, that's the signal to pivot.
