# NodeHaloTopology v2 — Bug Report (2026-05)

Source: parallel session by Daniel that ran v2 at scale (cube512/32 ranks
and cube1024/256 ranks). Current Mac session was about to add the same
fix from a different angle. Captured here so we don't re-derive it.

## Failure pattern (validated)

Across all mismatching ranks at 32+ ranks, the asymmetry is systematic:

```
v2 send sizes ≤ host send sizes   (often slightly smaller, sometimes 100+ off)
v2 recv sizes ≥ host recv sizes   (often larger, sometimes 2-3× larger)
```

Sample (cube512, 32 ranks, validate mode):
```
[Rank 25] send peer  4  size: host=12      v2=10
[Rank 25] send peer  6  size: host=27758   v2=27711
[Rank 25] recv peer  9  size: host=840     v2=1422   ← 1.7× larger
[Rank 25] recv peer 23  size: host=1841    v2=2676
[Rank 25] recv peer 27  size: host=941     v2=2779   ← 3× larger

[Rank 23] send peer 24  size: host=103962  v2=103744
[Rank 23] recv peer 17  size: host=60144   v2=61201
[Rank 23] recv peer 25  size: host=737     v2=2436
[Rank 23] recv peer 29  size: host=745     v2=2616
```

Cube16/4 ranks reproduces the same pattern (this session, 2026-05-09):
```
[Rank 3] recv peer 1  size: host=61   v2=130   ← 2.1× larger
[Rank 3] recv peer 2  size: host=404  v2=480
[Rank 0] recv peer 2  size: host=137  v2=205
[Rank 0] recv peer 3  size: host=215  v2=266
```

## Root cause (high confidence)

cstone halo bookkeeping doesn't model node ownership. v2's
`buildFromCstoneHalos` walks element corners in cstone halo ranges and
filters by local `d_nodeOwnership_`, but it never exchanges per-key
ownership claims with peers. So:

- A boundary node N may be touched by element halos from peers P1 and P2.
- v2 emits (peer=P1, node=N) and (peer=P2, node=N) into recvNodeIds_,
  because both halo ranges contain elements with N as a corner.
- After sort+unique on (peer, node), both entries survive — they're
  distinct (peer, node) pairs.
- Result: v2 expects to receive N from both P1 and P2.
- Host path resolves this via global `keyOwner[sfc_key]` map (built from
  `MPI_Allgatherv`), which says "N is owned by exactly one rank."
  Host's recv list for N has exactly one entry.

Send side mirror: when this rank's `d_nodeOwnership_` says it owns N
but the global tiebreaker (host) yields ownership to a lower-rank peer,
host removes N from send list while v2 keeps it.

## Numeric impact

Drift is small but **scale-dependent and growing** — more ranks → more
3+-way corners → more wrong-owner decisions → larger drift. CG would
converge to a wrong fixed point.

| Mesh / Ranks       | Host norm                          | v2 norm                            | Drift             |
|--------------------|------------------------------------|------------------------------------|-------------------|
| cube512 / 32       | matrix=1.463879e+01, rhs=2.993200e-03 | matrix=1.463887e+01, rhs=2.995761e-03 | ~5e-6 / ~9e-4 |
| cube1024 / 256     | matrix=2.071203e+01, rhs=1.555497e-03 | matrix=2.071215e+01, rhs=1.558524e-03 | ~6e-6 / ~2e-3 |

Mesh-load speedup confirmed:
- cube512/32: 68 s → 24 ms (2800×)
- cube1024/256: ~27 s → 28 ms (1000×)

Assembly slowdown with v2 halo (suspected: extra atomic contention from
duplicate contributions):
- cube1024/256: 9.2 ms (host) → 32.8 ms (v2)

## Fix design

**Add a peer-subgroup all-to-all tiebreaker** before §4.3 of the v2
build:

1. Each rank emits `(sfc_key, my_rank)` for every node currently flagged
   `ownership ∈ {1, 2}` (locally owned or shared). This is the rank's
   *claim*.
2. Send the claim list to each peer in `findPeersMac` set
   (point-to-point Isend/Irecv on device pointers via cstone's
   mpiSendGpuDirect / Recv).
3. On the receive side, for each incoming `(sfc_key, peer_rank)`:
   binary-search local `d_localToGlobalSfcMap_` for `sfc_key`. If found
   and `peer_rank < my_rank`: yield ownership (set
   `d_nodeOwnership_[local_id] = 0`).
4. Resync `d_nodeOwnership_` on device.
5. Then run the existing emit + filter + sort_unique with the now-
   authoritative ownership.

**Why all-to-all-in-peer-subgroup, not pairwise?** With pairwise
tiebreakers, A and B each compare against each other but neither knows
about C's claim on the same key. 3+-way ties resolve incorrectly. The
all-to-all in the peer subgroup ensures every rank that could touch a
key sees every other rank's claim, so min-rank-wins is computed
identically by all participants.

**Cost.** Boundary nodes per rank ≈ √(elem/rank). For 1M elem/rank ≈
10⁴ keys × 12 bytes × ~30 peers = ~3 MB sent per rank. One round of
peer p2p. O(boundary), no Allgatherv. Keeps the scaling win.

## Verification gate

```bash
# Cube16, 4 ranks (smallest reproducer):
MARS_NODEHALO_VALIDATE=1 srun --account=csstaff --nodes=1 \
  --ntasks-per-node=4 --time=00:10:00 ~/affinity/bind_numa.sh \
  daint-gpu/examples/distributed/unstructured/mars_cvfem_graph \
  --mesh=$SCRATCH/git/mars/cube16.mesh \
  --kernel=tensor --iterations=1
```

Expect: `[Rank N] NodeHaloTopo v2 == host (Gate 1 pass)` on every rank.

Then escalate:
```bash
# cube512/32 — the case the other session hit
srun -l --ntasks-per-node=4 --nodes=8 -A csstaff -t 00:15:00 \
  --export=ALL,MARS_NODEHALO_VALIDATE=1 \
  daint-gpu/examples/distributed/unstructured/mars_cvfem_graph \
  cube512.mesh/ --kernel=tensor --iterations=5
```

Then cube1024/256 once cube512 is clean.
