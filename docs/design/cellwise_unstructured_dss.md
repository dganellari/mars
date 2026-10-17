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

Each rank's tables cover its own elements and, after them, **ghost elements**: every
element of another rank that shares a vertex with one of its own (one element layer;
cstone's star-completed halo is the same set). All copies of every shared node are then
in the tables, so finding a node's holders needs no communication, the canonical order
and the counting copy come out the same on every rank, and the DSS is bit-identical to
one rank.

- **Ghost data**: a ghost element needs its corner keys (for the tables) and, at each
  DSS, the values of its shared nodes. These go into per-ghost blocks of 512 values
  (about 10% extra memory at 64^3 elements per rank), of which only the shared nodes
  are filled.
- **Lists** (setup, on the device): a local face whose other copy is a ghost sends its
  36 nodes to the ghost's owner and receives the ghost's; a star holding local and ghost
  copies receives every ghost copy and sends every local copy once to each rank that
  owns a ghost copy. Both sides sort these nodes by (peer, element id, local node), so
  the messages carry values only. At setup one all-to-all checks the counts and one
  exchange of hashed keys checks the order.
- **Exchange**: pack into one device buffer, GPU-aware Isend/Irecv per peer on a
  duplicated communicator, unpack into the ghost blocks. Elements that read no ghost
  copy (none of their vertex stars holds a ghost) run while the messages travel; the
  rest run after the wait.
- **Counting copy**: the first copy in canonical order may be a ghost, and then no local
  copy counts that node; every node is counted exactly once over all ranks.
- **Dirichlet**: a local element's face with no other copy among local and ghost
  elements is on the physical boundary. A ghost element's face that contains a local
  vertex always has its other copy present, so ghost faces never mark a shared edge or
  vertex by mistake.

Code: `UnstructuredHalo` in `mars_cellwise_topology.hpp`; `block_elements` gives a
structured sub-block and its ghost frame for tests. Measured on GH200 (rotated frames):
the DSS on 4 and 8 ranks (32^3 elements) is bit-identical to one domain, and the Jacobi
solve at ne = 8, deform 0.05 takes 151 iterations on one rank and 153 on 4 and 8, with
the same residual history to 4 digits.

## From an ElementDomain

`build_tables` (`mars_cellwise_domain.hpp`) makes one rank's tables from a hexahedral
`ElementDomain`: its local elements (the domain's local range, in its order) and the
ghost layer. The ghost layer is completed through the node owners (`SfcNodeOwner`), so
it does not depend on the reach of the cornerstone halo search, which completes stars
only for owned nodes:

1. Every rank sends each element to the owners of the element's vertices (one NBX
   `sparseExchange`); each owner then holds the complete star of every node it owns.
2. Each owner sends every star to the ranks that have an element in it (a second
   exchange).
3. The ghosts are the received elements, once each, that share a vertex with a local
   element.

A ghost travels as a record of 10 keys: its id, owner rank and 8 corner keys; it needs no
coordinates, because the operator never runs on it. The id is a hash of the sorted
corner keys: the domain's element SFC codes repeat (two elements can share the corner a
code is taken from), and the exchange needs ids that are unique within a rank (checked
at setup). Local corner ids are the ranks of the corner keys among the rank's distinct
keys. Needs SFC node ownership (the default; not `MARS_OWNERSHIP=vote`).

The domain's device constructor keeps only SFC-quantized coordinates, so a
manufactured solution is about 1e-7 instead of 0 on the boundary; `zero_dirichlet` sets
it to 0 on the Dirichlet copies before b = A u is formed.

Measured on GH200 (`mars_marsir_cellwise_solve --domain`, the cube of 8^3 elements with
deform 0.05 and rotated frames, handed to the domain as an unstructured mesh): DSS(1)
equals the copy counts implied by the geometry exactly, DSS(1 + x + 2y + 3z) equals
count x f to 7e-15, and the Jacobi solve takes 154 iterations on 1, 4 and 8 ranks, with
error 2.3e-11.

## Multigrid on unstructured meshes

**p-multigrid** on the same mesh, degrees 7 -> 3 -> 1, works for any imported mesh:
the Knaus CVFEM operator rediscretized at each degree, Chebyshev-Jacobi smoothing
through the DSS (degree 3, range [lmax/15, 1.1 lmax], lmax by 30 power steps, as in the
structured h-multigrid), element-local GLL interpolation between degrees (exactly
nested; restriction is its transpose on the raw residual), one V-cycle per BiCGStab
preconditioner call. The tables, counting bits and Dirichlet flags are the same at every
degree; only the node layout of each entity changes.

Numpy reference `marsir-mlir/test/cellwise_pmg_ref.py` (deform 0.05, (pre, post) = (0, 3),
BiCGStab to 1e-10):

| coarse solve at p = 1 | ne = 2 | 4 | 8 | 16 | 32 |
|---|---|---|---|---|---|
| exact (an ideal AMG) | 7 | 7 | 7 | 6 | - |
| 10 Chebyshev steps | | 7 | 7 | 7 | 9 |
| 3 Chebyshev steps | | | 7 | 9 | |

Without the p = 3 level (7 -> 1) it takes 15 iterations at ne = 8; (2, 2) smoothing saves
one iteration. The coarse level needs only a modest-accuracy solve, but a fixed number
of Chebyshev steps degrades with ne, so large meshes need AMG there. Plan:

1. GPU p-multigrid with k Chebyshev steps at p = 1 (`mars_cellwise_pmg.hpp`): all
   cell-wise, no assembly. The p = 3 and p = 1 operators are `cvfem_apply<P>`
   (`mars_cellwise_operator.hpp`), which takes its 1D matrices as a kernel argument: MARS's
   `ho_cvfem_apply_launch` reads them from one constant bank that holds a single degree,
   while a V-cycle alternates degrees. Same thread layout as that kernel (one thread per
   tangential slot, the element in shared memory), MARSIR's metric layout.
   GH200 (`mars_marsir_cellwise_solve --mg --unstructured`, deform 0.05, 10 coarse
   steps): lambda_max at p = 7, 3, 1 and the iteration counts equal the reference to the
   printed digits (7 at ne = 4 and 8); 7 with rotated frames, on 4 ranks and from an
   ElementDomain. ne = 64 in the table below (Jacobi needs about 1530 iterations).
2. hypre BoomerAMG at p = 1 on the assembled p = 1 operator (`mars_cellwise_amg.hpp`,
   `--amg`, needs a Hypre build). Global vertex numbers: the counting copy of each vertex
   writes its number, one DSS at p = 1 gives it to every copy. The 8x8 element matrices
   come from 8 unit-vector applies of `cvfem_apply<1>`; entries of rows owned by other
   ranks go to their owners (`sparseExchange`), and `hypreFromEntries` builds the ParCSR
   matrix on the GPU. Settings of `HypreAmgPcgSolver` (PMIS, extended+i, l1-Jacobi);
   `--amg-cycles k` V-cycles per coarse solve (default 1). GH200, one V-cycle: 7
   iterations at ne = 8 on one rank, on 4 ranks with rotated frames and on 4 ranks from
   an ElementDomain; AMG setup 640 ms at ne = 64.

GH200, ne = 64 (90.5 M unique nodes), deform 0.05, one rank:

| coarse solve at p = 1 | iterations | V-cycle | solve |
|---|---|---|---|
| 10 Chebyshev steps | 19 | 44.5 ms | 2.11 s |
| 30 Chebyshev steps | 11 | 51.9 ms | 1.41 s |
| 100 Chebyshev steps (`--coarse` default) | 6 | 76.8 ms | 1.13 s |
| one BoomerAMG V-cycle (`--amg`) | 10 | 41.7 ms | 1.10 s |

**h-multigrid** at fixed p = 7 when the mesh comes from uniform refinement of a coarse
hex mesh. Needs the parent key and octant of every element carried through the
cstone resync (the AMR does not keep them today), and siblings kept on one rank or
one more exchange in the transfers. Adaptive meshes would need the paper's local
smoothing and refinement-edge masks (arXiv 2607.03413); not in the first version.

## All-MARSIR form

Every device kernel of the Jacobi-preconditioned solve generated by MLIR passes
(marsir-mlir), next to the hand-written mixed form, which is kept:

- `mir.dss` (DSS on the `UnstructuredTopology` tables) with an optional per-copy
  `epilogue` region and reductions into per-element partials; the Jacobi
  preconditioner with its dot products is a six-line epilogue. `mir.cellwise` is the
  same per-copy map without the gather (the Krylov vector updates).
- `--mir-lower-cellwise` lowers both to per-element parallel phases (the hand
  kernel's staged gather and summation order; reference-hex tables from
  `mars_cellwise_hex.hpp`). `--mir-warp-per-element` maps them onto one warp per
  element: lane-cyclic phases, a barrier only between phases that share scratch,
  reductions as per-lane folds plus an xor-shuffle tree (a 32-slot OR is the
  Dirichlet ballot), element-level stores by lane 0. Partial folding and the Krylov
  scalar steps are structured loops from the front end, through the same mapping.
- `--mir-noalias-descriptors` restores `noalias`/`readonly` on kernel pointers (the
  MLIR 19 NVVM lowering drops them), so read-only loads use `ld.global.nc`.
- CPU gates (warp emulation, 32 threads per element, TSan): DSS and preconditioner
  bit-exact with the hand kernels' formulas on random tables with ghosts, the BiCGStab
  kernels bit-exact with host formulas, at p = 1, 3, 7.

GH200, deform 0.05, rotated frames (`test/dss_gpu_pipeline.py`,
`mars_cellwise_unstructured_dss --dss-ptx/--precond-ptx`,
`mars_marsir_cellwise_solve --precond-ptx --bicgstab-ptx`):

| 64^3, one rank | hand | MARSIR |
|---|---|---|
| DSS, plain / rotated | 1.05 / 1.19 ms | 1.25 / 1.36 ms (same bits) |
| precond + 2 dots | 2.22 ms | 2.13 ms (same z bits) |
| BiCGStab iteration | 20.64 ms | 20.73 ms (every kernel generated) |

Solve iterations: 151 at ne = 8 (153 hand), 151 on 4 ranks, 154 from an
ElementDomain on 4 ranks; 1388 at 64^3 (hand 1291), all converged. The dot products
are summed in another order, and Jacobi's ~1300 iterations amplify that rounding.
Tried for the DSS gap and measured slower: tables in shared memory, fewer blocks
striding over the elements, descriptors prefetched in the staging phase.

**p-multigrid V-cycle generated too** (`--mg --pmg-ptx`): Chebyshev steps (`mir.dss` with a
`minus` operand and an epilogue), masked residual copies (`mir.cellwise` with a Dirichlet
flag), transfers (`mir.interpolate`, sum-factorized), dots; the p = 3 and p = 1 operators
come from the same operator spec as p = 7 with the small-element schedule
(`--mir-cellwise-schedule`: many elements per block, a thread per output point, field and
result staged in shared memory, fills fused into reductions, scratch packed by liveness).
The tensor-core schedule at p = 3 sinks an mma under a lane guard (the convergence gate
rejects it), so it is not used there. Setup (diagonals, power iteration) and BoomerAMG
keep their code. GH200, 64^3, deform 0.05, rotated frames:

| | hand | MARSIR |
|---|---|---|
| V-cycle, 100 Chebyshev steps at p = 1 | 76.6 ms | 53.0 ms |
| V-cycle, BoomerAMG at p = 1 | 41.7 ms | 37.2 ms |
| solve (6 / 10 iterations) | 1130 / 1096 ms | 812 / 995 ms |
| Chebyshev step p = 1 / p = 7 | 0.30 / 4.08 ms | 0.11 / 3.60 ms |
| prolong / restrict 7 <-> 3 | 3.0 / 2.8 ms | 0.56 / 0.47 ms |
| operator p = 1 / p = 3 | 0.060 / 0.465 ms | 0.068 / 0.580 ms |

Same iterations everywhere (7 at ne = 8 on 1 rank and on 4 ranks from an ElementDomain).
The operator numbers are after `--mir-fuse-phases` (phases with the same iteration space
and only same-point dependences merged, also across non-conflicting phases, so the faces'
same steps form one phase; p = 3: 91 -> 21 phases, 77 -> 19 barriers) and with 4 warps per
block at p = 3 (11% faster than 8: 30 instead of 61 KB of shared memory); before fusion they
were about 0.09 / 1.3 ms, and 0.073 / 0.634 ms before register blocking (a short parallel
dimension an input does not depend on runs inside the thread: the sweep's column loaded once
for all faces, as in the hand kernel) and fused sibling reduction loops. Remaining: 1.14x /
1.25x of the hand kernels (BoomerAMG V-cycle not re-measured since fusion: 37.2 ms).

## Mixed precision (numpy, `marsir-mlir/test/cellwise_mixed_ref.py`)

Can the V-cycle run in FP32 inside the FP64 Krylov solver? Not with today's
left-preconditioned BiCGStab: an FP32 V-cycle is not linear at the 1e-7 level, the
solver's recursive residual still reaches 1e-10 in 7 iterations, but the true error stalls
at 1.5e-8 to 3e-8. With right preconditioning (A M y = b, x = M y; the residual is the true
assembled b - A x, the form flexible Krylov methods use for a varying preconditioner) an
FP32 V-cycle gives the same iterations and the same error as FP64, to the printed digits.
deform 0.05, p-multigrid 7-3-1 (0, 3):

| ne, coarse | left FP64 (today) | left FP32 | right FP32, stop 1e-8 | right FP32, stop 1e-10 |
|---|---|---|---|---|
| 4, exact | 7 it, 5.3e-11 | 7 it, 2.9e-8 | | 9 it, 2.3e-13 |
| 8, exact | 7 it, 2.6e-11 | 7 it, 2.3e-8 | 8 it, 4.9e-11 | 9 it, 7.2e-14 |
| 8, 100 Chebyshev | 7 it, 2.6e-11 | 7 it, 1.5e-8 | 8 it, 4.9e-11 | 9 it, 7.2e-14 |

(errors are the true algebraic error against the manufactured solution; right FP64 is
identical to right FP32 in every case.) So the path is a right-preconditioned BiCGStab
(one more DSS per operator application, to assemble the residual) with an FP32 V-cycle:
one more iteration at matched accuracy, for V-cycles that move half the bytes. GPU cost
not measured yet.

## Validation

1. Numpy spec: both kernel shapes and the multi-rank exchange, on structured blocks with
   every element in a random proper rotation of its frame (all 8 face orientations,
   reversed edges) and on extruded meshes with edges shared by 3 and 5 elements:
   copies bit-identical, equal to a canonical-order assembly, and 2/3/5-rank splits
   bit-identical to one domain.
2. GPU, one rank: DSS(1) = valence, DSS of a continuous field = valence x field,
   continuity, and on an unrotated cube the same values as the structured DSS.
3. GPU, several ranks: bit-identical to one rank (done, 4 and 8 ranks).
4. Solver: Jacobi BiCGStab on a rotated cube with the same iteration count as the
   structured solver (done: 62 at ne = 3, 151 at ne = 8 on one rank).

## Open risks

- Geometry storage: the p=7 metric is 32 KB per element (8x the values). The 1B-element
  target needs geometry recomputed from the corners inside the operator.
- Coordinates are SFC-quantized unless `storeOriginalCoords=true`.
- Rotated frames read 1.4x the DRAM of the plain cube (1.78 against 1.25 GB): the other
  copy of a face is read through the orientation map. Real meshes have such frames
  everywhere.
