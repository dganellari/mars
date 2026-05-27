# cstone Domain Reuse State Map (for AMR-style sync)

Source: parallel agent audit, 2026-05. Triggered by SEGFAULT after `setEndIndex(newCount)` + `sync()` with 8× larger element count.

## TL;DR

`setEndIndex()` alone doesn't work because three subsystems hold persistent state:

1. **Domain-level**: `bufDesc_` has 3 fields; only `.end` updates. `.start`, `.size`, and `prevBufDesc_` stay stale.
2. **FocusedOctree**: `rebalanceStatus_` precondition will throw on next `updateTree()` (it's set to `invalid` after first sync, but the precondition expects `valid`).
3. **Focus tree internal range**: `prevFocusStart` / `prevFocusEnd` are sticky; they assume incremental motion, not 8× particle count jumps.

Minimum patch: **~25 lines across 2 cstone files**, adds `Domain::resetForAMRRefinement(newCount)` + `FocusedOctree::resetFocusRangeForNewDistribution()`.

## Crash trace

User code:
```cpp
domain_->setEndIndex(newCount);   // Only sets bufDesc_.end
domain_->sync(...);                 // 8× larger arrays
```

What happens:
1. `distribute()` runs OK (uses `assignment_`, not `bufDesc_`).
2. `focusTree_.updateTree()` throws (or crashes) at the precondition `if (rebalanceStatus_ != valid)` (line 100 of `octree_focus_mpi.hpp`). After first sync `rebalanceStatus_ == invalid`.
3. If somehow that passes: `updateCounts()` walks 8× particle keys but loop bounds use stale `leavesAcc_` size → GPU memory access violation.

## State fields requiring reset

| Field | File:Line | Why critical | What to set |
|---|---|---|---|
| `bufDesc_.start` | domain.hpp:685 | Used in distribute() offsets | 0 |
| `bufDesc_.end` | domain.hpp:685 | already handled by setEndIndex | newCount |
| `bufDesc_.size` | domain.hpp:685 | Used in updateLayout() / setupHalos() | newCount |
| `prevBufDesc_` | domain.hpp:685 | Validated in reapplySync() | match new bufDesc_ |
| `layout_` | domain.hpp:698 | Stale w.r.t. new focus tree | `{0, newCount}` |
| `focusTree_.rebalanceStatus_` | octree_focus_mpi.hpp:794 | Precondition for updateTree() | `invalid` (0) |
| `focusTree_.prevFocusStart` | octree_focus_mpi.hpp:770 | Focus transfer logic assumes (0,0) for new region | 0 |
| `focusTree_.prevFocusEnd` | octree_focus_mpi.hpp:772 | Same | 0 |

Auto-managed (do not touch): `global_.*` (rebuilt in `assign()`), `halos_.*` (reset in `exchangeRequests()`), `focusTree_.leaves_*/counts_/macs_/centers_` (rebuilt in their respective updates).

## Proposed patch (25 lines)

**`cstone/domain/domain.hpp`** (add public method):
```cpp
void resetForAMRRefinement(LocalIndex newCount)
{
    bufDesc_.start = 0;
    bufDesc_.end   = newCount;
    bufDesc_.size  = newCount;
    prevBufDesc_   = bufDesc_;
    layout_.clear();
    layout_.push_back(0);
    layout_.push_back(newCount);
    focusTree_.resetFocusRangeForNewDistribution();
}
```

**`cstone/focus/octree_focus_mpi.hpp`** (add public method):
```cpp
void resetFocusRangeForNewDistribution()
{
    prevFocusStart   = 0;
    prevFocusEnd     = 0;
    rebalanceStatus_ = invalid;
}
```

**MARS-side caller** (replaces `setEndIndex` in `resyncFromDevice`):
```cpp
domain_->resetForAMRRefinement(elementCount_);
sync(d_conn_local, d_coords_local);
```

## Expected outcome

If patch correct, AMR rebuild at cube256/16:
- Pre-patch (default): ~7.9 s (firstCall_ converge dominates)
- Post-patch (Domain reuse): ~0.6–0.7 s

That's the 82% win the audit predicted, end-to-end measurable in MARS-side timing.

## Risk assessment per agent

> "**Ship the patch**: NOT fighting cstone's design—completing it. The author's `setEndIndex()` comment was incomplete; `resetForAMRRefinement()` is the API that should have been there."

The agent rates this as **upstream-friendly** because it implements the incomplete API the cstone author suggested in their own comments. Worth proposing to Sebastian Keller as a contribution.

## Caveats and unknowns

1. The agent's line numbers should be **verified** against the FetchContent copy before patching — local `cornerstone-octree/` may be a different version than what's fetched.

2. The patch hasn't been tested. The agent identified the state fields by reading the code; the actual fix may need 1-2 iterations of "patch, build, test, find next missed field." Default-off env gate must remain until verified.

3. **GB-scale risk**: at 10⁵ ranks, `firstCall_=true` also triggers a host-replicated global tree Allgather inside `distribute()`. Domain reuse skips converge but the Allgather inside `distribute()` may still scale linearly with rank count. The audit didn't measure this. Worth a follow-up if Domain reuse works.

## Recommendation

Apply patch in next session, gated under `MARS_AMR_REUSE_DOMAIN=1`. Validate cube16/4 → cube256/16 → cube512/32 bit-exact match progression. Ship if all three match host.
