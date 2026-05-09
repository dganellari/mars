# MARS + Cornerstone Reference Library

Durable, citation-heavy reference docs gathered from deep reads of both repos.
Use these as the source of truth for API surfaces, algorithms, and gotchas.
Each doc was written from a full-source read by an exploration agent and saved
here so we can re-consult across sessions.

Generated 2026-04-30 on branch `cstone`.

## Cornerstone-Octree (the dependency)

- [01_cstone_overview.md](01_cstone_overview.md) — repo layout, core types, Domain
  class, halo subsystem, SFC machinery, MPI utilities, reusable primitives
  table, build flags. The "what is in cstone and how do I call it" doc.
- [04_cstone_focus_tree.md](04_cstone_focus_tree.md) — algorithmic depth on
  Hilbert/Morton encoding, cornerstone tree, full octree linking, GlobalAssignment,
  SfcAssignment, makeSfcAssignment, FocusedOctree lifecycle (converge,
  discoverHalos, computeLayout), exchangeRequestKeys, findPeersMac dual-tree
  traversal, full sync() data flow, primitives reusable for MARS NodeHaloTopology v2.

## MARS (this repo)

- [02_mars_overview.md](02_mars_overview.md) — top-level architecture,
  ElementDomain layout, AdjacencyData CSR, HaloData / NodeHaloTopology, CVFEM
  assembly path, solver flow, AMR pipeline, examples, helpers, connectivity
  conventions, **identified bugs in current state**, build system.
- [03_mars_cvfem_kernels.md](03_mars_cvfem_kernels.md) — CVFEM math
  (SCS formulation, 12 sub-control-surfaces per hex, full vs graph sparsity),
  reference element, 12 kernel variants in detail (tensor / colored / AoS /
  perip / smem-cache / wmma / wgmma / etc.), coloring algorithm, dispatch
  enum, bottleneck analysis, roofline.
- [05_mars_amr_module.md](05_mars_amr_module.md) — AmrManager lifecycle,
  Doerfler vs OctreeNative marking, hex 1→8 refinement (19 new nodes), tet
  Bey's red, solution transfer by parentage, error indicators, rebuild flow,
  InitTimings/AmrStats, MPI patterns, octree integration with cstone, scaling
  gotchas, canonical AMR loop.
- [06_mars_solvers.md](06_mars_solvers.md) — solver inventory, CG / BiCGSTAB /
  GMRES / Hypre PCG, Jacobi & AMG preconditioners, Hypre integration, sparse
  matrix CSR, sparsity builder (graph 7 NNZ/row vs full 27 NNZ/row), DOF
  mapping, BC handling, per-iteration cost breakdown, identified bugs
  (int32 indexing, BC symmetry).

## How these were generated

Four parallel Explore agents, each pointed at a specific subsystem, with
instructions to actually read whole files (not just grep). Reports were
captured verbatim and saved here. Treat the file:line citations as a
snapshot — verify before acting on them, since both repos are evolving.

## How to use

When starting work that touches a subsystem, **read the relevant doc first**.
Don't trust the line numbers blindly — open the cited file and confirm. The
docs are reference material, not gospel; the code is the truth.

When something here turns out wrong or stale, fix the doc in the same commit
as the code change. These docs decay fast otherwise.
