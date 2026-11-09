# Conforming tetrahedral refinement

`TetRefinement::ConformingBisection` combines longest-edge propagation (LEPP)
with the existing cornerstone redistribution and domain rebuild. The legacy
Kokkos refiner and the distributed adapter share longest-edge selection,
terminal-edge decisions and child connectivity in `core/mars_simplex_bisection.hpp`.
The distributed implementation uses native device arrays; it does not require
Kokkos. Red refinement remains the default.

```cpp
using Manager = mars::amr::AmrManager<mars::TetTag, uint64_t, double>;
Manager::Config config;
config.tet_refinement = mars::amr::TetRefinement::ConformingBisection;
config.storeSideSets = true; // Preserve triangle side sets when reading Exodus.
Manager manager(config);
```

The candidate passes the maintained CPU, NVIDIA GH200 and AMD MI300A validation
matrix at one, two and four ranks, with the multi-rank cases spanning two nodes.
The GPU cases exercise the actual cornerstone rebuild, boundary matching and
field migration. The tested Cornerstone correction is published and pinned at
`c13b9a9a6c0f5897ec4f23283e4ff53716303c5e`. The policy remains experimental until
the integrated release regressions pass.
See [the results and validation instructions](../tests/amr/README.md).

## Refinement and conformity

Input tetrahedra must already form a conforming, nonperiodic, single-block mesh.
Only mark `1` requests refinement. A marked tetrahedron requests its longest
edge; equal lengths use the canonical pair of global endpoint keys. A requested
edge makes every incident tetrahedron request its own longest edge. Propagation
continues until requests agree across ranks. An edge is terminal only when every
incident tetrahedron selects it. All incident owned tetrahedra then bisect it,
including unmarked neighbours, before the next topology round.

A tetrahedron has one selected longest edge, so it cannot belong to two terminal
edge stars in the same round. Both sides of a shared triangle therefore perform
the same split. Endpoint replacement preserves signed orientation and halves
volume in exact arithmetic. Midpoints created in one round can be parents in a
later round. The two-tet regression distinguishes this algorithm from simply
splitting the initially requested edge: LEPP produces five tetrahedra, whereas
that shortcut produces four.

The distributed edge directory uses the SFC owner of the smallest endpoint key
and sparse peer exchange. Nodes use their own key; faces use their smallest key.
This follows the existing spatial partition instead of scattering records among
unrelated ranks. Self-routed records stay on device and bypass MPI and host staging.
The directory stores the complete local edge incidence in CSR form,
without the old fixed 30-incident limit. It exchanges edge records, not a global
copy of the mesh. One registration per incident rank identifies shared edge stars
exactly. Subsequent propagation exchanges only shared edges; local-only stars stay
local even when their coordinator is another rank. The communication plan is reused
during request propagation; incidence and routing are rebuilt after each topology-changing round. Incremental
updates of only the affected stars are not implemented, and no performance
improvement is claimed without measurements. Shared edges still exchange their
state on each propagation pass; changed-edge-only messages are not implemented.
Initial registration and node/face matching retain conservative owner-routed
records. A local halo alone is not used as proof that an edge star is complete,
and a locally paired face is not discarded if it could carry a remote label.

## Identity, boundaries and fields

A sparse node registry rejects inconsistent shared coordinates or two different
parent edges acquiring the same midpoint key. Before cornerstone rebuilds the
domain, another check rejects distinct coordinates colliding under the new box.
Both that check and the device constructor use the same bounding-box calculation.
SFC keys are identities within their box; they are re-encoded after box changes.

Triangle side sets carry face ancestry through each split. A parent face
containing the split edge produces two tagged child faces; the other faces
retain their labels on one child. Multiple labels on a face are supported.
After redistribution, sparse face matching restores labels and derives local
boundary-node lists on the device. `set_side_set_from_face_keys(name, faces)` accepts
flat triangle triples of current node keys for procedural meshes. Existing
node-only side sets cannot identify face ancestry and are rejected. Exodus
face data is materialized only on the first face-key request. Existing node-only
consumers do not build face-key vectors. File reading and coordinate resolution
still occur on the host; a later face request can cause a second side-set read.

P1 nodal transfer uses `u(midpoint) = 0.5*u(a) + 0.5*u(b)`. Birth batches ensure
that a later midpoint reads values already transferred to its parents. The
manager first refreshes required old ghosts. After rebuild it constructs one
sparse field-migration plan and reuses its routes and owned-node lookup across
all fields. Original-node values come only from old owners; newly created
midpoints come from their incident ranks. Conflicting authoritative copies and
missing owned-node values are rejected collectively. The existing node halo
then supplies active ghosts.

The default exchange passes device buffers to GPU-aware MPI. Setting
`config.stage_refinement_messages = true` stages the new LEPP, boundary and field
messages through host buffers. This does not change cornerstone or node-halo
transport: the integrated GPU path still requires GPU-aware MPI. Peer counts,
side-set names, convergence scalars and allocation sizes are host control data.

## Supported scope and remaining gates

The maintained GPU fixture uses 64-bit keys and double precision. Other type
instantiations passing a syntax check are not native validation. Local midpoint
parent indices must fit in 32 bits; existing domain/side-set index limits also
apply. Input with hanging faces is unsupported. Nonfinite coordinates,
degenerate input tetrahedra, duplicate local keys and detected key collisions
are rejected. These checks do not prove general mesh quality or shape regularity.

Coarsening and refinement-history reconstruction are not implemented. Mark `2`
does not request refinement, although closure may still split that cell.
Multiple blocks, periodic meshes, node-only boundary flags, vote-based ownership
and experimental domain reuse are rejected. `maxLevels` retains the existing
marker's focus-tree-depth meaning; it is not a bisection-generation cap.

Cornerstone decodes coordinates from finite-resolution SFC keys. Direct P1
transfer preserves affine fields to roundoff; post-rebuild checks additionally
allow the documented coordinate-quantization bound. The CPU adapter validates
the shared algorithm and communication; `AmrManager` integration uses the GPU
execution policy. CPU checks alone cannot establish GPU execution, distributed
cornerstone rebuilds, scaling or release readiness.

The three-round manager fixture passes conformity, volume, boundary, identity and
field checks on both GPU platforms. Its final element counts differ: 3315 on
NVIDIA and 3309 on AMD. Each count is consistent across the tested rank counts
and communication modes on that platform. Bitwise identical refinement topology
across GPU vendors is not established, and no scaling or performance claim is made.
