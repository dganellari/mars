# Element-local (cell-wise) DSS on unstructured hex meshes

Status: design, 2026-10-09. Extends the structured cell-wise solver
(`backend/distributed/unstructured/solvers/mars_cellwise_*.hpp`) to any conforming hex
mesh held in a cstone `ElementDomain`, on one or many GPUs. Executable spec:
`marsir-mlir/test/cellwise_unstructured_ref.py`.

## Credit

The method comes from M. Wichrowski's two papers: "Coalesced Matrix-Free Finite
Elements in Cell-Wise Storage" (arXiv:2607.02335) and "Coalesced Matrix-Free Geometric
Multigrid on Persistent Cell-Wise Storage" (arXiv:2607.03413). From them: element-local
storage as the solver state, DSS only inside the preconditioner, the equivalence with
the assembled solve, the face/line/vertex DSS with orientation codes, and the
cell-wise geometric multigrid. MARS adds the nonsymmetric CVFEM operator (BiCGStab),
p = 7, element-granularity tables built on the GPU, canonical-order sums, and the
multi-GPU exchange. The papers are single-GPU.

## What it is

Vectors stay element-local: every element stores its own 8x8x8 copy of the p=7 GLL
nodes it touches (node (a,b,c) at (a*8+b)*8+c, axis a from corner 0 to 1, b from 0 to
3, c from 0 to 4: the VTK corner order `ElementDomain` keeps). The operator works
element by element. Elements communicate only through the DSS: every copy of a node
gets the sum of all its copies.

The structured version finds copies by index arithmetic. Here they come from compact
tables, built once on the device from the elements' corner SFC keys.

## Shared entities

A p=7 hex has 216 interior nodes (never shared), 6 x 36 face-interior nodes, 12 x 6
edge-interior nodes and 8 corners. The shared node sets of faces, edges and vertices
are disjoint, so each can be handled by its own kind of work item, in any order:

- **Face pair**: the two elements of an interior face, with a 3-bit relative
  orientation (the dihedral map between their face frames).
- **Edge star**: every (element, local edge, reversed bit) holding the edge. Any
  valence (3, 4, 5, 6, ...).
- **Vertex star**: every (element, local corner) holding the vertex.

Every copy sums all copies of its node in the **canonical order** and keeps the sum: a
face node adds its two copies (addition is commutative, so their order does not
matter), an edge or vertex node sums its star in canonical order. This is the paper's
face/line/vertex decomposition (Wichrowski, arXiv 2607.02335, Alg 3) at element
granularity, with no macro-blocks. All copies of a node get the same bits and nothing
is atomic.

**Canonical order**: copies sorted by their element's global identity, the sorted
tuple of its 8 corner SFC keys (compared lexicographically, at setup only). It does not
depend on the partition or on local numbering, so the DSS gives the same bits on any
number of ranks, and every copy of a node holds identical bits.

Face frames use the convention of `hexFaceCanonicalPosDev`
(`fem/mars_ho_dof_handler_gpu.hpp`): origin at the smallest corner key, first axis
toward the smaller of its two neighbours, computed from GLOBAL keys. Edge position runs
from the lower to the higher corner key.

## Kernels

One warp per element, elements in grid-stride order, so the elements in flight form a
contiguous window and, on an SFC-ordered mesh, the neighbours' copies are still in L2
when they are read. Per element:

1. Stage the 512 values in shared memory with one coalesced load.
2. Faces: one loop adds the other copy to each of the 216 face nodes, through a
   64 x 36 orientation table (frame codes of both copies, position on the face).
3. Stars: one loop sums the stars of the 72 edge nodes and 8 corners.
4. One coalesced store writes all 512 results through the consumer.

Interior nodes need no work. Each loop handles one kind of node, so the lanes of a warp
follow one path, and the reference-hex index arithmetic comes from tables built once
per thread block. The element's descriptors are loaded once, each on the lane whose
number is the entity's bit (faces 0-5, edges 8-19, corners 20-27), and shuffled to the
lanes that need them; one ballot gives the Dirichlet mask.

Consumers: `dss` (every copy = s); the Jacobi preconditioner (z = 0 on Dirichlet nodes,
else s / diag, plus the dot products of the Krylov step); later the Chebyshev step.

Dot products count each global node once: exactly one copy per node is the **counting
copy** (the first in canonical order). Each element has a 32-bit counting mask in the
same bit layout, with bit 31 (interior nodes) always set. The weights stay exact, which
the structured path got from powers of two.

Measured on one GH200, 64^3 elements (90.5 M unique nodes), ncu:

| kernel | time | DRAM read / write |
|---|---|---|
| this one, plain cube | 1.25 ms | 1.25 / 1.06 GB |
| this one, rotated frames | 1.34 ms | 1.78 / 1.06 GB |
| structured gather (cascade order) | 2.20 ms | 1.67 / 1.06 GB |
| one work item per shared node, writing all its copies (removed) | 2.9 ms | 4.4 / 2.15 GB |
| one warp per element in node order, no staging (removed) | 4.1 ms | 1.46 / 1.06 GB |

The one-writer form writes every memory sector twice, in parts, each time from a
different pass; the node-order warp form was bound by instructions (2.25e9 warp
instructions against 0.27e9 here), because every warp mixed all four kinds of node.

## Building the tables (device, setup only)

Input: per local element, the 8 corner keys (`ElementDomain::indices<I>()`), and dense
local node ids (`getElementToNodeConnectivity()`) for grouping.

1. Faces: emit (sorted 4 local corner ids as two uint64, payload e*6+f), stable radix
   sort as in the HO DOF handler (stages 1-2), adjacent equal keys give pairs.
   Orientation from the global keys of each side's frame.
2. Edges: (packed sorted local corner ids, e*12+k), radix sort, segments give stars.
   Reversed bit = key(corner at local t=0) > key(corner at t=p).
3. Vertices: (local corner id, e*8+c), radix sort, segments give stars.
4. Canonical order inside each star (small segments): insertion sort per star by the
   element identity tuple.
5. Dirichlet flags: faces with one copy on the physical boundary; their edges and
   vertices inherit the flag.
6. Counting-copy bits.

No host loop over elements. The existing `FaceTopology` (`mars_face_topology_gpu.hpp`)
is not used: it groups on the host and drops the 4th face corner.

## Several ranks

Each rank holds only its local elements' copies. Copies on other ranks arrive by a
**single-phase exchange of raw copies**: every rank sends its copies of each shared
node to every other rank that holds the node, and all ranks then sum all copies in the
canonical order. Results are bit-identical to one GPU, as in the structured version.
(The alternative, per-rank partial sums combined in rank order, sends fewer values at
edges and vertices but gives up identity with one GPU. The surface is dominated by
faces, which have one copy per side either way.)

- **Holder discovery** (setup): each rank sends (entity key, its rank, the identities of
  its copies) for every entity on a rank boundary to a rendezvous rank,
  `SfcNodeOwner(min corner key)` (`mars_sfc_ownership.hpp`), over the NBX
  `sparseExchange`. The rendezvous rank returns every holder's full list. This does not
  depend on the node-halo peer list, which can miss co-holders under SFC ownership,
  and needs no node owner (HoHalo needs `MARS_OWNERSHIP=vote`, DofSpace needs SFC
  ownership).
- **Lists**: per peer, the copies to send in canonical order; on receipt they fill
  ghost slots that `face_nbr` and the star tables reference at their canonical
  positions. Counts follow from the shared entity sets, so sends and receives match by
  construction.
- **Exchange**: persistent device buffers, a duplicated communicator, GPU-aware
  Isend/Irecv. Elements with no ghost copies run while the messages travel; the rest
  run after the wait.
- **Counting copy**: the first copy in canonical order may be remote, in which case no
  local copy counts that node.

HoHalo and the `elemDof` scatter/gather stay useful as setup-time oracles.

## Multigrid on unstructured meshes

- **p-multigrid** (7 -> 3 -> 1) on the same mesh, with exactly nested element-local
  transfers, then AMG (hypre BoomerAMG, already in MARS) on the assembled p=1 operator.
  Works for any imported mesh; this is the general path.
- **h-multigrid** at fixed p=7 when the mesh comes from uniform refinement of a coarse
  hex mesh. Needs the parent key and octant of every element carried through the
  cstone resync (the AMR does not keep them today), and siblings kept on one rank or
  one more exchange in the transfers. Adaptive meshes would need the paper's local
  smoothing and refinement-edge masks (arXiv 2607.03413); not in the first version.

## Validation

1. Numpy spec: both kernel shapes and the multi-rank exchange, on structured blocks with
   every element in a random proper rotation of its frame (all 8 face orientations,
   reversed edges) and on extruded meshes with edges shared by 3 and 5 elements:
   copies bit-identical, equal to a canonical-order assembly, and 2/3/5-rank splits
   bit-identical to one domain.
2. GPU, one rank: DSS(1) = valence, DSS of a continuous field = valence x field,
   continuity, and on an unrotated cube the same values as the structured DSS.
3. GPU, several ranks: bit-identical to one rank.
4. Solver: Jacobi BiCGStab on a rotated cube with the same iteration count as the
   structured solver.

## Open risks

- Geometry storage: the p=7 metric is 32 KB per element (8x the values). The 1B-element
  target needs geometry recomputed from the corners inside the operator.
- Coordinates are SFC-quantized unless `storeOriginalCoords=true`.
- Rotated frames read 1.4x the DRAM of the plain cube (1.78 against 1.25 GB): the other
  copy of a face is read through the orientation map. Real meshes have such frames
  everywhere.
