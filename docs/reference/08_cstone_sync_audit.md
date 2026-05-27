# cstone::Domain::sync() AMR Bottleneck Audit

Source: parallel agent audit, 2026-05. Measured cost at cube256/16 AMR Level 1: ~7.87 s per rebuild.

## TL;DR

**MARS triggers `focusTree_.converge()` on every AMR rebuild** because `rebuildDomainFromDevice()` constructs a fresh `Domain`, which has `firstCall_=true`. cstone's converge loop runs synchronously with `MPI_Allreduce` per iteration. The data payload is tiny (1 int per Allreduce), but the loop runs 4–10 iterations, each with local tree rebalance, and **doesn't reuse prior focus tree state across AMR levels**.

The 60 MB Allreduce in the global octree update is **not** the bottleneck — the convergence loop's repeated tree rebuilds + GPU halo discovery on scattered indices are.

## Root cause

[mars_amr_manager.hpp:619](backend/distributed/unstructured/amr/mars_amr_manager.hpp#L619):
```cpp
domain_ = std::make_unique<Domain>(...);   // fresh Domain, firstCall_ = true
```

[cstone/domain/domain.hpp:188-191](cornerstone-octree/include/cstone/domain/domain.hpp#L188):
```cpp
if (firstCall_) {
    focusTree_.converge(box(), keyView, peers, global_.assignment(), ...);
}
```

The previous focus tree (which is mostly identical because only boundary/refined regions changed) is discarded.

## Estimated cost breakdown at 977M global / 16 ranks

| Phase | Estimate | Observed |
|-------|----------|----------|
| `distribute()` global tree + assignment | 250 ms | included in 7.87 s |
| `focusTree_.converge()` (firstCall) | 50–100 ms (4 iters × ~10 ms + Allreduce latency) | dominant unnecessary work |
| `discoverHalos()` GPU kernel | 100–300 ms expected | ~1–2 s actual (scattered MAC walk at scale) |
| `setupHalos()` particle exchange | 200–400 ms | ~500–1000 ms (load imbalance) |
| `updateLayout` reordering | ~50 ms | included |

Realistic total without firstCall: ~600–800 ms. Observed: 7870 ms. **Gap of ~10×.**

## Recommendations (priority order)

### Priority 1 — incremental sync path (30–50% speedup per AMR rebuild)
Reuse the existing Domain object across AMR levels instead of constructing a new one. Either:

- **Option A (simpler)**: Don't reconstruct Domain. Hold it across AMR rebuilds and call `sync()` with new connectivity. Works only if cstone internally accepts changed element count without re-converging from scratch — needs verification.

- **Option B (cleaner, requires cstone patch)**: Add `Domain::incrementalSync(...)` that skips `focusTree_.converge()` when only a subset of elements changed. Pass the set of refined/coarsened elements as a hint.

### Priority 2 — pre-warm Domain to flip firstCall_ false (5–10%)
Add an explicit dummy `sync()` after Domain construction so subsequent AMR rebuilds see `firstCall_=false`. Workaround, not a real fix.

### Priority 3 — optimize discoverHalos GPU kernel (10–30%)
[cstone/halos/exchange_halos_gpu.cuh](cornerstone-octree/include/cstone/halos/exchange_halos_gpu.cuh) — vectorize the MAC octree walk. Out of scope without cstone patch; flag for upstream.

### Priority 4 — bucket size sweep (5–15%)
Try bucketSize=128 and 256. Smaller global tree → faster converge. Trade-off: coarser load balance.

### Priority 5 — batch AMR rebuilds (10–15%)
Accumulate 2–3 mark+refine cycles before triggering a sync. Requires AMR control-flow refactor.

## Confidence

**High confidence**: focusTree_.converge() triggers on every AMR rebuild via firstCall_ — verified in code at the cited file:line.

**Medium confidence**: convergence loop is the dominant cost. Numbers in the table are estimates; we haven't profiled with nsys yet to confirm. Validation step before acting on Priority 1 should be a `cudaProfiler` / `nsys` run on cube256/16 AMR L0→L1.

**Low confidence**: per-iteration count of converge loop. Estimate is 4–10; actual depends on spatial heterogeneity and could be higher post-AMR refinement.

## Action plan if we do Priority 1 Option A first (simplest test)

1. In `mars_amr_manager.hpp`, replace `domain_ = std::make_unique<Domain>(...)` with calling a new `domain_->resync(...)` method that:
   - Replaces internal element / coord vectors
   - Calls `cstoneDomain_->sync()` (the same one — `firstCall_` is now `false`)
2. Validate matrix/rhs norm bit-exact against current Domain-replacing path.
3. Measure rebuild wallclock. Expected drop: 7.87 s → ~5 s if convergence is the issue, ~3 s if more.

If Option A doesn't show meaningful improvement, focusTree_.converge() isn't the dominant cost; pivot to discoverHalos profiling.

## Cross-check

Conflicts with my earlier hand-wave that "cstone's global tree is not distributed" — partially accurate (it IS replicated for partition decisions), but the data volume isn't the issue at this rank count. The real issue is repeated work the algorithm doesn't need to redo on rebuild.
