# MFEM's Node-Ownership and Halo Exchange — Reference for MARS

Sourced from a parallel deep-read of MFEM source by an Explore agent
(2026-05). Saved here as durable context. Citations are to MFEM source
files; verify against current github HEAD before treating as gospel.

## TL;DR

MFEM avoids the runtime claim/merge protocol entirely. **Node ownership
is deterministic and pre-computed at mesh construction**:

1. For each shared vertex, MFEM builds an `IntegerSet` containing the
   ranks of all elements touching it.
2. `IntegerSet::Recreate()` (`general/sets.cpp`) sorts the ranks
   ascending and removes duplicates.
3. `GroupTopology::Create()` (`general/communication.cpp`) calls
   `PickElementInSet(group_id)`, which returns `data[0]` — the **lowest
   rank in the sorted set**.
4. That lowest rank is the master for the entire simulation. No
   per-iteration negotiation, no claim merge, no min-rank-wins recompute.

The **3+-way corner problem doesn't arise** because every rank that
touches the same vertex computes the same sorted IntegerSet → the same
lowest-rank master.

## Key data structures

| Structure | File | Purpose |
|-----------|------|---------|
| `IntegerSet` | `general/sets.cpp` | Sorted, deduped rank list per shared entity |
| `ListOfIntegerSets` | `general/sets.cpp` | One IntegerSet per group |
| `GroupTopology` | `general/communication.cpp` | Master per group, group ⇄ ranks topology |
| `ldof_group[ldof]` | `fem/pfespace.cpp` | Group ID this DOF belongs to (0 = local-only) |
| `ldof_ltdof[ldof]` | `fem/pfespace.cpp` | True-DOF index if I'm master, -2 if not |
| `svert_lvert` | `mesh/pmesh.hpp` | shared-vertex → local-vertex map |
| `group_svert` | `mesh/pmesh.hpp` | group → shared-vertex Table |

## Construction flow (per rank)

```
1. Mesh partitioner assigns rank to each element (input).
2. ParMesh builds vertex→element adjacency (Table).
3. For each vertex, replace element indices with their partition ranks.
4. If this rank's MyRank appears AND others appear → vertex is shared.
5. Build IntegerSet from those ranks; sort ascending, dedupe.
6. GroupTopology stores all distinct IntegerSets.
7. For each group: master rank = IntegerSet[0] (lowest).
8. ParFiniteElementSpace marks each DOF as -1 (mine) or -2 (someone else's).
9. Master ranks number their owned DOFs; broadcast counts so slaves can
   compute global true-DOF numbering.
```

This is O(local boundary), runs once.

## Per-iteration halo exchange

`GroupCommunicator` provides `Bcast()` (master → slaves) and `Reduce()`
(slaves → master, with reduction op).

**Two modes:**
- `byGroup`: per-group point-to-point. One Isend per (group, slave).
- `byNeighbor` (default): aggregate all DOFs destined for neighbor N
  into one buffer, one Isend per neighbor. Minimizes message count.

**Not persistent** (fresh Isend/Irecv per call), **not Alltoallv**. All
MPI is point-to-point organized by neighbor rank.

## How a 3-way corner resolves

Vertex shared by ranks {1, 3, 7}:
- Each rank independently builds `IntegerSet = {1, 3, 7}` (sorted).
- All three call `PickElementInSet → 1`.
- Rank 1 marks DOF as owned (-1 → ltdof = N).
- Ranks 3 and 7 mark DOF as -2.
- Per-iteration: rank 1 broadcasts to ranks 3 and 7. Ranks 3 and 7
  receive from rank 1.

Zero ambiguity. Zero per-iteration logic to determine who sends to whom.

## What MFEM doesn't do (and why it matters)

- **No "I might own it / they might own it" runtime resolution.**
  Ownership is set in stone at mesh construction.
- **No SFC-key key/value claim list exchange.** MFEM uses local element
  adjacency; the rank-set comes from the partitioner output.
- **No yielding of duplicate ownership during runtime.** If a corner is
  shared, the sorted rank list answers "who owns it" deterministically.

## How this differs from our v2

| Aspect | MARS v2 | MFEM |
|--------|---------|------|
| When ownership is set | At first NodeHaloTopology build (per-iteration if rebuilt) | Once at mesh construction |
| Claim mechanism | (sfc_key, my_rank) sent to peers, min-rank-wins on receive | Sorted IntegerSet of ranks touching the entity |
| 3+-way correctness | Requires every rank in the corner-set to be in the peer subgroup | Automatic — every rank touching the vertex sees the same rank set |
| Per-iter cost | Pack/MPI/Unpack (correct shape) | Same: pack/MPI/Unpack |
| Initialization cost | One subgroup p2p exchange + sort_unique merge | Local table build + small Allgatherv on shared-DOF counts |

## Recommendations distilled

The agent gave 10 recommendations; the load-bearing ones for us:

**R1. Replace claim-merge with sorted-rank-set determinism.** Don't ask
"who claims it" — ask "what's the sorted set of ranks that touch it?"
Lowest is master. Same answer on every rank because the set is the same.

**R2. Build ownership once.** GroupTopology in MFEM lives for the
simulation. Our `NodeHaloTopology` already does this — but our v2 attempt
re-derives ownership each construction; it should be cached and reused.

**R3. byNeighbor aggregation.** Already what `exchangeNodeHalo()` does
post-construction. ✓

**R5. Validate consensus with Allreduce(MAX) of owner per key.** Cheap
debug check we should add to `validateAgainstHost`.

**R6. Use SFC key as the unique node identifier.** Already what we do
via `d_localToGlobalSfcMap_`. ✓

**R8. Test 3-way and 8-way corners explicitly.** We're hitting this
on cube16/4 (8-way center), cube256/16 (multiple 8-way corners), etc.

## What this means for our in-flight v2 fix

The bug we hit (rank 0 receiving zero recvs) is not a fundamental
architectural problem. Our approach **is** sorted-rank-set determinism,
just expressed differently:

- MFEM: every rank builds `IntegerSet{ranks_touching_corner}` from local
  element adjacency, then takes lowest.
- MARS v2: every rank publishes its claim, every other rank merges all
  peers' claims, takes lowest.

Mathematically equivalent **if the published claim set on every rank
exactly equals the rank's true element-touching set**. The bug we found
(`d_nodeOwnership` was Step-2-biased, so we under-claimed) is exactly
the failure mode of "claim set ≠ touching set." The fix already in
flight (`markLocallyClaimedKernel` rebuilding "I claim from local
elements") closes that gap.

So the MFEM read is **validating, not redirecting**. Our trajectory is
correct; we just had a subtle bug in *how* we computed "I touch this
corner via a local element."

## Two specific changes to consider after Gate 1 passes

Both are post-correctness optimizations, not blockers:

1. **Cache ownership across NodeHaloTopology rebuilds.** AMR triggers
   a rebuild via `domain_ = make_unique<Domain>(...)`. If the partition
   doesn't change, ownership doesn't change either. Skip the tiebreaker
   exchange.

2. **Add the Allreduce-MAX consensus check (R5) to `validateAgainstHost`.**
   Currently we only diff against the host path's output. Adding "every
   rank's owner array is consistent with every other rank's" catches
   bugs the host comparison can't (e.g., if the host path itself has a
   determinism issue at scale that we've never noticed).

## File index (MFEM source)

```
mfem/general/
  sets.cpp             IntegerSet::Recreate (sort + dedupe)
                       IntegerSet::PickElement (returns data[0])
  communication.hpp    GroupTopology, GroupCommunicator, Mode enum
  communication.cpp    GroupTopology::Create
                       GroupCommunicator::Bcast / Reduce
                       byGroup vs byNeighbor MPI patterns
mfem/mesh/
  pmesh.hpp            ParMesh public surface (svert_lvert, group_svert)
  pmesh.cpp            FindSharedVertices, FindSharedFaces, etc.
mfem/fem/
  pfespace.hpp         ParFiniteElementSpace
  pfespace.cpp         ConstructTrueDofs (-1/-2 marker logic)
                       True DOF numbering broadcast
```

## Cite-with-care note

Agent didn't paste exact line numbers because it browsed via web
fetches; treat the file paths as authoritative but verify line numbers
against current MFEM HEAD if you need to cite them in a paper.
