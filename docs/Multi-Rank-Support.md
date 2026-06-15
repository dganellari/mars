# Multi-Rank Support

MARS runs distributed across multiple GPUs with MPI. Each rank owns a contiguous
space-filling-curve (SFC) chunk of the mesh plus a halo; all mesh data lives on the
device, and MPI is used only to partition the mesh and to exchange ghost data on the
device (CUDA-aware MPI, no host round-trips).

This page is the conceptual overview. The concrete distributed mechanics are documented
in detail elsewhere — this page points to them rather than re-deriving them:

- [Halo Management](Halo-Management.md) — the element halo and node-ownership semantics.
- [Node Halo Topology](Node-Halo-Topology.md) — the per-node halo (`NodeHaloTopology`,
  `exchangeNodeHalo`) used for per-solver-iteration ghost-DOF exchange.
- [AMR Module](AMR-Module.md) — multi-rank adaptive refinement and its validation
  status.

> **Multi-rank status (honest).** Non-periodic distributed assembly and AMR are
> rank-invariant and validated (the cube16/4-rank DOF count matches single-rank
> exactly, and assembled norms match to several digits). Periodic (Taylor–Green)
> multi-rank and some inlet-driven channel multi-rank paths have **known limitations
> under active work** — see [AMR Module → Known limitations](AMR-Module.md#known-limitations)
> and the [periodic TGV tutorial](periodic_tgv_tutorial.md) for the precise current
> state. Single-rank is the validated reference for those cases.

## How partitioning works

There is no separate "MPI manager" class. MPI is passed into the domain as a
`(rank, numRanks)` pair and handled by the embedded cornerstone domain:

```cpp
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
class ElementDomain {
    // cstone owns the SFC partition + element halo discovery across ranks
    std::unique_ptr<DomainType> domain_;   // cstone::Domain<KeyType, RealType, AcceleratorTag>
    int rank_, numRanks_;
    // ...
public:
    int rank()     const { return rank_; }
    int numRanks() const { return numRanks_; }
};
```

The rank and rank count are taken at construction:

```cpp
MPI_Init(&argc, &argv);
int rank, numRanks;
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

// Each rank reads its share, builds SFC keys, and hands them to the cstone domain.
ElementDomain<HexTag, double, uint64_t, cstone::GpuTag> domain(meshFile, rank, numRanks);

if (rank == 0)
    std::cout << "Rank 0 owns " << domain.getElementCount() << " elements (owned + halo)\n";

// ... run assembly / solve ...
MPI_Finalize();
```

The same binary runs on one rank or many; the mesh is split automatically by the SFC.
For a runnable end-to-end example see
[`examples/distributed/unstructured/mars_cvfem_graph.cu`](../examples/distributed/unstructured/mars_cvfem_graph.cu)
and the [Quickstart](Quickstart.md).

## Two halo layers

MARS keeps **two** distributions of the same physical mesh, the standard
dual-decomposition used by other distributed FEM frameworks (PETSc DMPlex, Trilinos
Tpetra, MFEM):

1. **Element halo (cstone-native).** `cstone::Domain` partitions and load-balances
   *elements* by SFC key and discovers the per-element halo. Per-element device arrays
   are exchanged with `exchangeHalos(...)`. This is used during mesh build and for
   per-element fields.

2. **Node halo (`NodeHaloTopology`).** FEM/CVFEM computes on *nodes* (DOFs), not
   elements. A separate node-keyed topology is built lazily and exchanged with
   `exchangeNodeHalo(nodeArray, nodeToDof)` — the hot-path call wired into the solver's
   ghost-exchange callback, running once per Krylov iteration. See
   [Node Halo Topology](Node-Halo-Topology.md) for construction and the pack/exchange/
   scatter path.

```cpp
// Per-element field exchange (cstone element halo)
domain.exchangeHalos(std::tie(d_field), sendBuffer, receiveBuffer);

// Per-node DOF exchange (node halo) — typical solver halo callback:
auto haloExchange = [&domain, dofMap = d_node_to_dof.data()]
    (cstone::DeviceVector<RealType>& p) {
        domain.exchangeNodeHalo(p, dofMap);
    };
```

## Node ownership convention

Multi-rank assembly hinges on a single ownership map, `getNodeOwnershipMap()` (a
`DeviceVector<uint8_t>`, one entry per local node), with three states:

- `0` — **pure ghost**: the node is owned by another rank; this rank holds a copy only.
- `1` — **pure owned**: owned by this rank and interior to it.
- `2` — **shared**: owned by this rank and also sitting on a halo boundary.

For assembly, **owned and shared both get an equation row** — the test is
`(owner == 1 || owner == 2)`; only a pure ghost (`0`) is skipped. This is the same
convention documented in [FEM Assembly](FEM-Assembly.md) and used throughout the
kernels. The `NodeHaloTopology` build also runs a lowest-rank-wins tiebreaker so that
each globally-unique node is owned by exactly one rank (this is what makes the
distributed DOF count match the single-rank count exactly).

## DOF numbering across ranks

`buildDofMappingGpu` (in `mars_cvfem_utils.hpp`) turns the ownership map into a local
DOF numbering entirely on the device via two `exclusive_scan`s over the owned/ghost
masks: owned nodes are numbered `[0, numOwned)`, ghosts `[numOwned, numTotalDofs)`. It
returns `numOwned`, the number of equation rows this rank owns. The **sum of `numOwned`
across ranks equals the global unique node count** — each shared node is counted by
exactly one rank. The solver uses `numOwned` to scope its `MPI_Allreduce` on dot
products so that shared DOFs are not double-counted.

## Characteristic sizes (computed during construction)

Per-node characteristic sizes are computed on the GPU during domain construction by
`calculateCharacteristicSizes(...)` and stored on the device. They are mesh-local; any
global statistic (min/max/sum across ranks) is an ordinary `MPI_Allreduce` over the
rank-local values — there is no MARS-specific collective wrapper.

## Adaptive refinement across ranks

Multi-rank AMR follows mark → refine → rebuild → transfer, rebuilding the cstone domain
(and therefore both halo layers) after refinement, and using `exchangeNodeHalo` to fill
ghost DOFs during solution transfer. See [AMR Module](AMR-Module.md) for the pipeline
and its known limitations.

## See Also

- [ElementDomain Overview](ElementDomain-Overview.md)
- [Halo Management](Halo-Management.md)
- [Node Halo Topology](Node-Halo-Topology.md)
- [FEM Assembly](FEM-Assembly.md)
- [SFC Mapping](SFC-Mapping.md)
- [AMR Module](AMR-Module.md)
