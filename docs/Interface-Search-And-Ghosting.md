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

- **`ElementDomain::create_ghosting(name)`** — nothing calls the registry yet.
- **Incremental `change_ghosting(add, remove)`.** OpenAccel calls this every step
  (`interface.cpp`, `updateGhostings_()`), because a *sliding* interface re-pairs continuously.
  Our registry is build-once, so a sliding interface would pay a full rebuild per step.
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
