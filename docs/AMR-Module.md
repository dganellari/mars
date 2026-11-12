# AMR Module

*Step 10b of the step-by-step tutorial (see the [index](index.md)). Previous: [step 10a](Matrix-Free-Tutorial.md); next: [the appendix](step-appendix-cmake.md).*

The native AMR module in `backend/distributed/unstructured/amr/` refines Hex8 and
Tet4 meshes on CUDA/HIP devices. `AmrManager` uses `ElementDomain` and Cornerstone
for SFC partitioning, redistribution, tree construction and halo exchange. It
supports nonconforming red refinement and an optional
[conforming tetrahedral LEPP policy](Conforming-Tet-AMR.md). Neither requires Kokkos. The separate older Kokkos bisection module lives in
`backend/kokkos/bisection/`; see [the backend layout](Backend-Layout.md).

## Refinement and field transfer

Each adaptation performs these steps:

1. Mark owned elements from their errors, using a threshold or Cornerstone leaves.
2. For red refinement, exchange halo marks and enforce mesh 2:1 balance until
   every rank reports no further changes. LEPP instead closes shared edge stars.
3. Emit children only for owned parents. Red refinement produces eight children
   per marked parent; LEPP bisects terminal edges and may refine neighbouring tets.
4. Fill old field halos before interpolating new edge, face and interior values.
5. Rebuild through Cornerstone, which redistributes elements and their properties,
   rebuilds its trees and halos, and lets `ElementDomain` reconstruct node storage.
6. Route authoritative field values to the new SFC node owners, then fill needed
   ghost nodes. Reuse one sparse migration plan for all fields in that adaptation.

For nonconforming refinement the original bounding box stays fixed. A node uses a
canonical point inside its SFC cell, so a later emission of the same hanging node
recovers the same key and position. Fresh element and node key buffers are cleared
before encoding: Cornerstone preserves removal markers, which must not be inherited
from reused storage. This policy quantizes coordinates at the key resolution.

An existing node's old owner supplies its value exactly, including when that node
is emitted again as a new midpoint. If a key is entirely new, a deterministic
emitter supplies the interpolated value. Conflicting old owners and missing
sources cause collective failure instead of silently inserting zeros. Ranks with
no local input elements still participate in the rebuild and communication.

## Marking and Cornerstone reuse

`MarkingStrategy::Doerfler` is the existing maximum-relative threshold policy:
`error > refineFraction * globalMaximum`. Despite its historical name, it does
not select a set whose accumulated error reaches a bulk fraction.

`MarkingStrategy::OctreeNative` uses Cornerstone's existing global leaves and
sorted element keys. Each rank computes errors for its owned leaf ranges and calls
`computeNodeOpsGpu` to decide splits. A scalar maximum is reduced across ranks;
there is no additional global array reduction of leaf errors. Assignment boundaries
align with global leaves, and counts forbid sibling merging, so remote leaf errors
are not needed for the owned split decisions.

`AmrOctree` obtains element levels from Cornerstone's current focused-tree layout,
including halos, rather than reconstructing a second element-to-leaf map. Physical
mesh balance still uses element sizes and connectivity: spatial tree balance alone
does not constrain hanging nodes. Cornerstone retains its own global tree; this
change does not eliminate that data structure or establish large-scale performance.

## Using the manager

```cpp
using Manager = mars::amr::AmrManager<mars::TetTag, uint64_t, double>;
Manager::Config config;
config.strategy = mars::amr::MarkingStrategy::OctreeNative;
config.maxLevels = 3;
Manager manager(config);
manager.initialize(meshFile, rank, numRanks);

// Errors cover local elements, including halos; marking counts owned elements.
cstone::DeviceVector<double> next;
auto stats = manager.adaptMesh(errors.data(), solution.data(), next);
solution = std::move(next);
```

`adaptMeshMultiField(errors, oldFields, newFields)` transfers several node fields.
Every rank must supply equally sized field lists and valid pointers for its local
arrays. Use SFC node ownership for distributed transfer. `AmrStats` records owned
element counts, refined parents and phase timings; `shouldContinue(stats)` checks
the configured level/error limits and whether any elements were refined.

## Galerkin equations with hanging nodes

Red refinement leaves hanging edge nodes and, for hexes, face nodes. A continuous
P1/Q1 field expresses their values as weighted combinations of independent nodes.
`HangingNodes` finds these relations, routes them through SFC owners and resolves
dependent constraints. `DofSpace` then supplies prolongation `P`, additive
restriction `P^T` and the constrained operator `P^T A P`. Dirichlet unknowns are
removed symmetrically before PCG.

`mars_amr_galerkin` exercises this path for Poisson on public cube meshes, with
P1 tetrahedra or Q1 hexahedra and BoomerAMG-preconditioned CG. `--problem=linear`
checks an affine manufactured solution, allowing for coordinate quantization;
`--problem=peak` reports an error against the localized exact solution.
The [peak comparison](AMR-Peak-Example.md) uses `--refinement=red|lepp` and
`--uniform` to compare local and uniform resolution at the same equation and solver settings.
`--audit-mesh` and `--audit-algebra` enable the optional conservation and algebra
checks described in [the validation guide](../tests/amr/nonconforming/README.md).

## Validation scope

The source snapshots tested on NVIDIA GH200 and AMD MI300A pass the combined
Cornerstone/red-refinement cube checks at one, four and eight ranks, including
repeated eight-rank tet cases across two nodes. Audited rebuilds preserve owned
cell counts and two cell fingerprints, with no invalid or inverted tets.
Assembled and composed constrained operators agree within the `1e-10` gate.

Separate immediate-transfer checks pass twelve cases per GPU platform: tet/hex,
one/eight ranks, both marking policies, staged AMR messages and empty local input.
Every case performs three adaptations before checking preserved old values and
constant/affine fields, independently of a subsequent PDE solve. CPU/MPI checks
exercise the shared sparse migration and tree-layout logic at one/two/four/eight
ranks, plus actual CPU Cornerstone redistribution and property exchange.
These CPU checks do not turn the native GPU manager into a CPU backend.

See the validation guide for job IDs, source seals, exact tolerances and the
separate status of integration/package testing. These results do not establish
scaling or general-mesh accuracy.

## Limits

- Refinement only: coarsening is not implemented.
- Nonconforming Galerkin acceptance covers the public nonperiodic single-block
  cube fixtures. Periodic constraints, mixed elements, arbitrary geometry and
  nonconforming side-set ancestry are not established by these checks. Use the
  separately validated conforming Tet4 policy for triangle side-set preservation.
- Nonconforming CVFEM and Navier–Stokes AMR drivers do not constrain hanging nodes
  or supply the missing coarse/fine interface flux. A successful run does not
  establish a consistent discretization there.
- Empty local input is covered; a globally undersized mesh with fewer elements
  than ranks remains unsupported.
- `--staged` stages AMR sparse messages. It does not replace all Cornerstone or
  assembled-operator GPU-aware MPI traffic with host transport.
- The gradient-jump indicator is not an error against an exact solution, and its
  current CVFEM implementation can depend on the partition. Old warm-start timing
  examples are not performance evidence for this implementation.

## See also

- [Conforming tetrahedral AMR](Conforming-Tet-AMR.md)
- [Halo Management](Halo-Management.md)
- [Node Halo Topology](Node-Halo-Topology.md)
- [Solution Transfer](Solution-Transfer.md)
- [ElementDomain Overview](ElementDomain-Overview.md)
