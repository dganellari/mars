# Interface Search and Ghosting

The two capabilities MARS needs before it can stand in for STK in a code like OpenAccel:
a **geometric search** to pair up the two sides of an interface, and a **custom ghosting** to move
the opposing side's data across ranks. Together they are what non-conformal (sliding) interfaces
and overset grids are built on.

Both are GPU-native. Each header carries two implementations of the same algorithm — the device one
under the plain name, and a host reference suffixed `Host` that exists so the logic can be pinned
down with `mpirun` on a laptop and so the device path has something to be compared against
element-for-element.

| Capability | Header | STK equivalent |
|---|---|---|
| Coarse search | `fem/mars_coarse_search.hpp` | `stk::search::coarse_search` |
| Fine search | `fem/mars_fine_search.hpp` | Nalu's `DgInfo` projection |
| Custom ghosting | `fem/mars_ghost_registry.hpp` | `stk::mesh::Ghosting` |

---

## 1. Coarse search — `coarseSearch`

Given **domain** boxes (the queries, e.g. slave-side face AABBs) and **range** boxes (the targets,
e.g. master-side faces), both distributed, return every overlapping `(domain, range)` pair with the
owning rank of each side.

```cpp
mars::fem::CoarseSearch<double> cs;          // hold one per interface: the scratch persists
thrust::device_vector<mars::fem::SearchPair> pairs;
cs.search(d_slaveBoxes, d_slaveIds, nSlave,
          d_masterBoxes, d_masterIds, nMaster,
          peers, comm, pairs);               // peers: a superset is fine
```

`SearchPair` holds two `IdentProc{id, proc}` — STK's own shape, because the partner usually lives
on another rank and whoever consumes the pair has to know where to fetch the data.

**Communication plan.** Every exchange is point-to-point over `candidatePeers`, the same contract
`GhostRegistry::build` takes. There is deliberately **no `MPI_Alltoall` or `Allgather`** anywhere:
an earlier draft used them for the envelopes and counts, which is O(P) per rank in three places and
does not survive large rank counts. Peers come from cstone's discovery
(`cstone/traversal/peers.hpp`) or from the domain's existing halo peer list.

Routing is one **envelope** (bounding box) per peer, not a distributed tree. Interface and overset
boxes live on a *surface*, so an envelope prunes most of the traffic and the fan-out stays small.
A tree buys nothing until the query set stops being a surface.

**STK parity**, checked against Trilinos `develop` source rather than documentation:

- The **overlap predicate is identical**. STK tests disjointness as
  `amax[d] < bmin[d] || bmax[d] < amin[d]`, which is what `Aabb::overlaps` computes, and it is
  **closed** — at `amax == bmin` the strict `<` is false, so touching boxes *do* intersect. That is
  what a conformal seam needs, where opposing faces share an edge exactly.
- `enforceSymmetry` mirrors STK's `enforceSearchResultSymmetry` and **defaults to true**, as STK's
  does: the pair reaches the range owner as well as the domain owner. The DG-IP assembly needs
  exactly that, since master rows carry opposing-element columns and vice versa.
- We always sort and deduplicate; STK makes that optional (`sortSearchResults`, default false).
  Ours is the stricter guarantee, so the result never depends on message arrival order — which is
  what lets the device path be diffed against the host one.

## 2. Fine search — `projectPointToTriangle` / `selectDonor`

The coarse search narrows a point to a handful of candidate faces. The fine search decides **which**
of them owns it and **where** on it — the barycentric weights the DG-IP flux interpolates with.

```cpp
auto hit = mars::fem::selectDonor(p, faceIds, cornerA, cornerB, cornerC, nCand);
// hit.face, hit.w0/w1/w2, hit.distance, hit.inside
```

Closed form, no iteration: a P1 tet face is a planar triangle. Uses the seven-region closest-point
test (Ericson, *Real-Time Collision Detection*, ch. 5) — three vertex regions, three edge regions,
then the interior.

**Why not "project onto the plane, then clamp the weights":** independently clamping barycentric
coordinates does not give the closest point once the triangle is obtuse. The region test does.

**Points are allowed to land outside every candidate.** On a curved or slightly mismatched
interface the exact projection often falls just outside all of them; the nearest *clamped*
projection still wins, and `inside` records which happened. Without that, such integration points
would silently find no donor.

## 3. Custom ghosting — `GhostRegistry`

A **named** ghosting, matching `bulkData.create_ghosting(name)`: several can be live at once, each
with its own duplicated communicator, and owner ranks are tracked per ghost rather than as a
per-entity array.

```cpp
GhostRegistry<KeyType> reg("fluid-solid-interface");
reg.build(nEntities, owner, key, participates, myRank, peers, comm);
if (!reg.isConsistent()) { /* the two sides disagree about identity or participation */ }

reg.forward(d_field);      // owner -> ghost, ghost slots overwritten
reg.reverseAdd(d_field);   // ghost -> owner, contributions summed (the assembly direction)
```

**`reverseAdd` uses no atomics.** `sendIdx_` is indexed per `(peer, entity)`, so an owned entity
ghosted by several peers is written by several slots. `atomicAdd` would be fast enough — contention
is 1–3 peers per entity — but it sums in completion order, so the result would not be bitwise
reproducible run to run. The collision pattern comes from the ghosting topology and is **fixed at
build time**, so `sendIdx_` is inverted once into `entity -> contributing slots` and one thread
sums each entity in a fixed order. Deterministic, and free per call.

## 4. Building a ghosting from a domain — `createNodeGhosting`

`fem/mars_domain_ghosting.hpp`. The equivalent of `bulkData.create_ghosting(name)`:

```cpp
auto ghosting = createNodeGhosting(domain, "fluid-solid-interface", comm);   // all shared nodes
auto iface    = createNodeGhosting(domain, "iface", comm, myMask);           // a chosen subset
```

A **free function taking the domain**, not a method on `ElementDomain`, because `domain.hpp` is a
lower layer than `fem/` and having it include the registry would invert that.

Everything the registry needs is already in the node halo topology. The domain stores only an
owned/not-owned flag, but the halo's **recv list per peer** says which foreign rank each ghost came
from, so the owner-rank map falls out of it. The peer list is by construction a superset of the
ranks that can own a shared node — a node's owner is always a mesh neighbour.

The default mask marks every node this rank sends *or* receives. That matters: a node must be
marked on the OWNER as well as on every rank that ghosts it, or the two sides disagree and
`isConsistent()` reports it. Taking it from both halo lists satisfies that by construction.

Setup-time only — the registry's build is host-side, so this copies the halo lists and SFC keys
down once. The per-step exchange stays on the device.

---

## Testing

| Gate | What it covers |
|---|---|
| `mars_coarse_search_host_mpi.cpp` | Host reference: fixture topology plus a **randomized brute-force oracle**. The overlap set is *defined*, not chosen, so an O(N²) serial pass over the global boxes is ground truth. |
| `mars_ghost_registry_host_mpi.cpp` | Host reference: send/recv symmetry, forward, reverseAdd, ownerRank, two-registry isolation, inconsistent-participation detection. |
| `mars_search_ghost_mpi.cu` | **Device vs host on identical input**, element-for-element. Also the only translation unit that compiles the device halves at all — both headers are header-only and every other consumer is a `.cpp`. |

Run the device gate at **4 ranks or more**: the `reverseAdd` collision it is really testing cannot
occur on one rank.

## Not done yet

- **Incremental `change_ghosting(add, remove)`.** OpenAccel calls this every step
  (`interface.cpp`, `updateGhostings_()`), because a *sliding* interface re-pairs continuously.
  Our registry is build-once, so sliding would pay a full rebuild per step. Designed below,
  not written.

  ### What STK actually does

  Read from Trilinos `develop`, not from documentation.

  **The remove path communicates.** `BulkData::internal_change_ghosting`
  (`stk_mesh/base/BulkData.cpp`) hands the receiver's remove list to

  ```cpp
  stk::mesh::impl::comm_sync_send_recv(*this, removeRecvGhosts, newSendGhosts, removeSendGhosts);
  ```

  which packs, for each recv-ghost being dropped, its `EntityKey` **to that entity's owner**
  (`GhostCommHelper::pack_and_communicate_buffers`, `stk_mesh/baseImpl/MeshImplUtils.cpp`). The owner
  unpacks them into `removeSendGhosts` as `(key, proc)` and then erases its own obligation in place:

  ```cpp
  entity_comm_map_erase(key, EntityCommInfo(ghosting.ordinal(), proc));
  ```

  The same round carries the mirror case — a rank that wants to ghost an entity it does *not* own
  forwards `(key, destProc)` to the owner, which pushes it onto `newSendGhosts`. So there is one
  symmetrisation step, and it is receiver-driven in both directions. **The delta reasoning above is
  confirmed by the source**: an owner cannot learn to stop sending by local bookkeeping.

  **They genuinely edit, and they never store an ordered list.** `Ghosting` holds no arrays at all —
  `Ghosting::send_list`/`receive_list` (`Ghosting.cpp`) are *derived on demand* by scanning the one
  global `internal_comm_list()` and filtering `ec->ghost_id == m_ordinal`. The persistent state is a
  single `EntityCommDatabase`: entity → `{(ghost_id, proc)}`, edited by `insert`/`erase`, with
  emptied rows marked `key = EntityKey()` and swept later by one
  `remove_if(IsInvalid())` in `delete_unneeded_entries_from_the_comm_list()`.

  **The trick we were missing is the ordering invariant, not the protocol.** `communicate_field_data(const Ghosting&, ...)`
  (`FieldParallel.cpp`) computes its send sizes, its recv sizes and its pack order from a single walk
  of `internal_comm_list()` — and that list is a vector kept **sorted by `EntityKey`**
  (`insert_keep_sorted_and_unique`, `inplace_merge` in `add_comm_list_entries_for_entities`). Both
  sides of a pair walk the same global key order, so slot alignment is a *consequence of the keys*
  and is never negotiated. That is why STK can splice one entity into the middle of a ghosting
  without telling anybody a position.

  Our `build()` pairs slots by **the requester's local index order** instead. Correct, and cheaper to
  build — but it is exactly what would make an edit expensive, because inserting one key in the
  middle renumbers every slot after it. **So the first change is the invariant, not the protocol:
  each peer's send/recv list becomes sorted by key.** Then an edit is a sorted merge on each side
  independently and no position ever crosses the wire.

  **No cheaper trick exists in the source.** There is no always-send-and-discard path, and no
  generation or version counter used for message matching — `synchronized_count()` is a
  cache-invalidation stamp for the *device* map (`volatile_fast_shared_comm_map` regenerates wholesale
  whenever it advances), not part of any protocol. Worth knowing: STK does **not** update its device
  index maps incrementally either.

  ### What not to copy

  **STK's remove makes a round trip.** OpenAccel's owner computes `sendGhostsToRemove` from its own
  `send_list`, ships those keys to the receivers (`communicateToFillRecvGhostsToRemove`), and STK's
  `comm_sync_send_recv` then ships them straight back to the owner — two rounds to move one fact,
  because the public API only accepts a *receiver*-side remove list. Ours should take the remove list
  in whichever direction the caller already has it and communicate only the missing direction.

  **`comm_sync_send_recv` is O(P).** It uses `stk::CommSparse` with no proc hint, which resolves to
  the unknown-pattern exchanger and counts receives with `MPI_Ireduce_scatter_block` over a
  `parallel_size()`-length array (`ReceiveCounter.cpp`). STK knows the cost — `CommSparse.hpp` says
  of the two-argument overload that "*Allowing the send/recv procs to be specified allows the
  CommSparse class to avoid a round of communication*" — but the ghosting path does not take it. We
  already have the peer list, so we pay neither.

  **Closure has no analogue here.** Half of `internal_change_ghosting` is
  `filter_ghosting_remove_receives` and `has_upward_recv_ghost_connectivity`, refusing a removal while
  some surviving ghost still needs the entity in its closure. Our entities are a flat keyed set, so a
  caller that ghosts elements and wants their nodes expands the set itself. The registry stays flat.

  ### The MARS protocol

  `update(addKeysPerPeer, removeKeysPerPeer)`, receiver-driven, point-to-point over
  `candidatePeers` — **not** `peers_`, which `build()` compacts to peers with traffic, so a peer that
  first acquires traffic this step would otherwise have no channel. `candidatePeers` therefore has to
  be retained by `build()`.

  Two messages per peer, both on the registry's own duplicated communicator:

  | | Contents | Tag |
  |---|---|---|
  | 1 | `{nAdd, nRemove}` — symmetric over `candidatePeers`, so no deadlock and no probing | `kTagDeltaCounts` |
  | 2 | `nAdd` added keys then `nRemove` removed keys, one buffer | `kTagDelta` |

  Owner side, per peer, on the `(key, localId)` pair vector that replaces the bare `sendIdx_` run:

  1. sorted-difference out the removed keys; a key that is not currently served is a **divergence**, counted, not silently dropped;
  2. resolve the added keys through the retained `keyToLocal` table — unresolved keys keep the existing `-1` sentinel and the existing counter;
  3. sorted-merge the added keys in, dropping duplicates so a repeated add cannot make the owner send twice.

  Requester side applies the identical edits to its own `(key, localId)` run with no message — it
  computed the delta, so it already knows the answer. Both sides then recompact the CSR
  (`peers_`, `sendOffsets_`, `recvOffsets_`, `sendIdx_`, `recvIdx_`), rebuild `ghostOwner_`, and set
  `idxUploaded_ = false` so the next `forward`/`reverseAdd` re-uploads the index lists and re-inverts
  the `reverseAdd` collision CSR. Regenerating the device side wholesale is the same choice STK makes,
  and it is cheap next to the exchange it feeds.

  **State that has to be kept** beyond what `build()` stores today, all O(interface):

  - the key alongside every send and recv slot — `build()` currently discards `gotKeys`/`reqKeys` once resolved;
  - the sorted `keyToLocal` table, currently a build-local temporary, so an added key resolves without rescanning all entities;
  - `candidatePeers`, for the reason above;
  - `myRank_`, `comm_` — already kept.

  ### Failure modes

  | Mode | Symptom | Guard |
  |---|---|---|
  | Delta reaches one side only | send and recv counts disagree; the *next* exchange misaligns and writes plausible garbage | counts round is symmetric over `candidatePeers`; a per-peer hash of the post-edit key list, compared in the same round, turns the misalignment into an error at the step that caused it |
  | Remove of a key not served | the two sides' idea of the current set has already diverged | count it and `MPI_Allreduce`, exactly as `globalUnresolvedRequests()` does for adds |
  | Duplicate add | owner sends twice, requester expects once | the merge drops duplicates — the same scrub STK does in `verify_and_filter_add_send` |
  | Owner changes under a sliding interface | ownership is an input to `build()`, not to `update()` | out of scope: migration means rebuild. `update()` must assert the ownership of every key it touches is unchanged |
  | Empty delta | a wasted round every step at the far end of a slide | zero-count messages still cost a round trip; skip the key message per peer, keep the count message |

  ### Gates

  | Gate | What it proves |
  |---|---|
  | `IncrementalMatchesRebuild` | the one that matters: drive a random add/remove sequence through `update()`, then `build()` the same final set from scratch, and require `peers_`, both offset arrays and both index arrays to match element-for-element. It makes the new path answerable to the path already trusted. |
  | `SendRecvSymmetryAfterUpdate` | the existing `SendRecvSymmetry` gate, re-run after every step of the sequence rather than once |
  | `ForwardAfterRemoveIsExact` | a removed slot is left untouched — neither zeroed nor written with a stale neighbour value — and every surviving slot is still exact |
  | `RemoveThenAddSameKey` | idempotence, and no double-send |
  | `EmptyDeltaIsNoOp` | `update({}, {})` leaves the CSR bitwise identical |
  | `NewPeerAppears` | an add whose owner carried no traffic before still lands, i.e. the delta really goes over `candidatePeers` |
  | device re-run in `mars_search_ghost_mpi.cu` | `forward`/`reverseAdd` still match host after an update, which is also the only check that the upload invalidation fired |
  | sliding smoke test | N steps of a rotating pairing with `isConsistent()` every step, and `reverseAdd` of all-ones still summing to the ghost count |

  At 4 ranks or more, for the same reason the device gate is.

- **The local overlap test in the coarse search is brute force**, O(nQ × nRange) per peer. Fine for
  interface surfaces; the batching around it is arranged so cstone's `collisions_gpu` traversal can
  be dropped in without touching the communication plan.

## Measured, and rejected

**Persistent MPI requests** (`MPI_Send_init`/`Startall`) for the count exchange. The sizes are fixed
at one int per peer, so it looked like the textbook case — but benchmarked against `Isend`/`Irecv`
on 2–12 ranks it was consistently slower and degraded with peer count: **+1.0 %** at 4 peers,
**+3.3 %** at 8, **+18.6 %** at 12 (median of 2000 iterations). The count exchange is also only
1–5 % of the routine's time to begin with. Measured on shared-memory OpenMPI; worth re-checking on
Alps/cray-mpich before drawing conclusions about the network case.
