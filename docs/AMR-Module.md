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
2. For red refinement, exchange generation requests at shared vertices and enforce mesh 2:1 balance until
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

`MarkingStrategy::OctreeNative` uses Cornerstone's existing global leaves and
sorted element keys. Each rank computes errors for its owned leaf ranges and calls
`computeNodeOpsGpu` to decide splits. A scalar maximum is reduced across ranks;
there is no additional global array reduction of leaf errors. Assignment boundaries
align with global leaves, and counts forbid sibling merging, so remote leaf errors
are not needed for the owned split decisions.

Both policies use `refineFraction` as a **maximum-relative threshold**. The
historical `Doerfler` name does not mean bulk-error marking: an element is marked
when its error is strictly greater than `refineFraction * globalMaximum`.
`OctreeNative` applies that test to each global leaf's largest owned-element
error, then marks the elements in that leaf. It passes `bucketSize + 1` to
Cornerstone for a selected leaf and `bucketSize` otherwise. Thus changing the
bucket size changes spatial grouping, but not the error threshold. Zero error
marks nothing. A finest-resolution SFC leaf can no longer split spatially, but
its selected elements still refine up to their mesh-generation cap. Several
elements can share one key; spatial tree depth is not their refinement limit.
Coarsening remains disabled.

Mesh generations are separate from Cornerstone tree depth. Input cells start at
zero; every red split or LEPP bisection increments its parent's generation.
Cornerstone carries this integer alongside connectivity through its existing
sort, redistribution and halo exchange. No second partition or replicated mesh
is constructed. Cornerstone still retains its global spatial tree.

Red closure compares generations at shared vertices and at the regular masters
of hanging vertices. A reusable sparse communication plan collects the maximum
requested generation at each SFC node owner. Sweeps continue to collective
convergence. Physical edge lengths do not enter this decision, so an initially
graded conforming mesh does not refine merely because neighbouring cells differ
in size. This is a conservative vertex-balance rule, not a minimal face-only closure.

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
the adaptation-round/error limits and whether any elements were refined.
`maxLevels` limits each element's refinement generation; `maxAdaptations` separately
limits manager calls that change the mesh. Neither is a Cornerstone leaf depth.
Red marking cannot exceed the generation cap. If LEPP conformity would require a
bisection beyond that cap, the call fails collectively before replacing the domain;
it never returns a partially conforming mesh.

## Galerkin equations with hanging nodes

Red refinement leaves hanging edge nodes and, for hexes, face nodes. A continuous
P1/Q1 field expresses their values as weighted combinations of independent nodes.
The refiner records each split edge or hex face by its sorted parent corner
keys and the key of the node it emitted. After redistribution, exact entity
matching retains records only where a coarse edge or face still exists.
`HangingNodes` consumes these records, routes them through SFC owners and resolves
dependent constraints. There is no distance tolerance. Equal midpoint keys alone
do not establish incidence: reusing an old key requires the same recorded parent,
and an unrelated collision fails collectively. Hex body-centre collisions are
also rejected. `HangingNodes` requires manager refinement history; a file containing
an already nonconforming mesh cannot reconstruct that history from coordinates. `DofSpace` then supplies prolongation `P`, additive
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

The generation, exact-incidence and explicit-threshold changes pass 23 local
CPU/MPI checks and 46 cases on each of NVIDIA GH200 and AMD MI300A. The native
campaign covers 1/4/8 ranks, including two nodes, with graded Tet4/Hex8 and thin
Tet4 refinement, constraint algebra, field transfer, nonbinary marking and
maximum-generation checks. Conforming LEPP also passes three refinement rounds
on the public non-box MFEM Escher mesh. All checked mesh/constraint counts agree
across ranks and platforms. See the [validation record](../tests/amr/README.md#generation-and-topology-validation)
for jobs, source identity and transport scope.

These tests establish mesh and constraint behavior beyond cubes. They do not
establish PDE accuracy on every distorted mesh. The affine Galerkin solver checks
still use cube domains; the example's boundary selection is based on its bounding box.

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
scaling or arbitrary-geometry PDE accuracy.

## Limits

- Refinement only: coarsening is not implemented. Initial meshes must be conforming,
  single-block and nonperiodic. Keep one refinement policy per mesh history:
  red generations and directional bisection generations are not interchangeable.
  Refinement history is not serialized for restart.
- Nonconforming Galerkin acceptance covers the public nonperiodic single-block
  cube fixtures, while refinement and constraint checks also cover graded/thin
  meshes. Periodic constraints, mixed elements, arbitrary-geometry PDE accuracy
  and nonconforming side-set ancestry remain unvalidated. Use the
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
