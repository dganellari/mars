# Step 1: the mesh

The first step reads a mesh onto the GPUs and prints how MARS splits it between MPI ranks. There is
no physics yet. The aim is to understand the three ideas everything else builds on: the domain, the
halo, and node ownership.

The program is `examples/distributed/unstructured/mars_step01_mesh.cu`, about 60 lines.

## Run it

Build MARS with the examples (see the [Quickstart](Quickstart.md)), generate a cube of 16×16×16
hexahedra, and run on one and on four ranks:

```bash
python3 scripts/generate_hex_cube.py -n 16 --output cube16
./build/examples/distributed/unstructured/mars_step01_mesh --mesh=cube16
mpirun -np 4 ./build/examples/distributed/unstructured/mars_step01_mesh --mesh=cube16
```

On one rank it prints one row: the rank holds all 4096 elements and all 4913 nodes, with no halo. On
four ranks each row holds about a quarter of the elements plus a halo, and the last line is the same
as on one rank:

```
mesh: 4096 elements, 4913 nodes
```

## The code

### One GPU per rank

```cpp
int deviceCount = 0;
cudaGetDeviceCount(&deviceCount);
if (deviceCount > 0) cudaSetDevice(rank % deviceCount);
```

Each rank picks its GPU before anything touches the device. On a node with four GPUs, ranks 0 to 3
get GPUs 0 to 3. The same code runs on AMD GPUs: MARS maps the CUDA calls to HIP.

### The domain

```cpp
using Domain = ElementDomain<HexTag, double, uint64_t, cstone::execution::Gpu>;
Domain domain(meshDir, rank, numRanks);
```

`ElementDomain` is the mesh, stored on the GPU. Its template arguments are the element type
(`HexTag` or `TetTag`), the floating-point type, the integer type of the space-filling-curve keys,
and where the data lives.

The constructor does three things. It reads the mesh, it sorts the elements along a space-filling
curve (a line through 3D space that keeps nearby elements close together), and it gives each rank one
contiguous piece of that curve. That piece is the rank's *own* elements. Each rank also receives a
*halo*: copies of neighbouring elements that belong to other ranks, so that all the elements around
each node the rank owns (see below) are available locally.

### Own elements and halo elements

```cpp
const auto& conn = domain.getElementToNodeConnectivity();
long held  = std::get<0>(conn).size();
long owned = domain.localElementCount();
long halo  = held - owned;
```

The connectivity has one column per element corner (eight for a hexahedron), and one entry per
element the rank holds. The rank's own elements are the range `[startIndex(), endIndex())`; the halo
elements come before and after it. Loops that must count every element exactly once over all ranks,
like the volume in [step 2](step-02-volume.md), use only the own range.

### Owned nodes and ghost nodes

```cpp
long nodes        = domain.getNodeCount();
const auto& owner = domain.getNodeOwnershipMap();
long ownedNodes   = thrust::reduce(thrust::device, owner.data(), owner.data() + nodes, 0L);
```

The nodes are numbered locally from 0 on every rank: all nodes of all elements the rank holds. A
node on the border between two ranks appears on both, so MARS makes exactly one of them its owner.
The ownership map has one byte per local node: 1 if this rank owns it, 0 if it is a *ghost*, a copy
whose owner is another rank. Solvers compute on owned nodes and refresh the ghosts from their owners
when they need them.

`thrust::reduce` sums the map on the GPU. `thrust::device` says where it runs; the data never goes
through the host.

### Counting over all ranks

```cpp
long global[2] = {owned, ownedNodes};
MPI_Allreduce(MPI_IN_PLACE, global, 2, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
```

Every element is owned by one rank and every node by one rank, so the sums over all ranks are the
mesh's element and node counts, whatever the number of ranks. That is what the last line prints, and
what the release check `marsReleaseStep01_npN` tests.

## Try this

- Run on 2, 3 and 8 ranks: the halo grows with the number of rank borders.
- Generate a longer box, `--nx 64 --ny 16 --nz 16`, and see how the pieces of the curve change shape.

Next: [step 2](step-02-volume.md) runs a first computation on the GPU.
