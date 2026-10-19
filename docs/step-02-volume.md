# Step 2: a first loop on the GPU

Step 2 computes something from the mesh: the volume of every element, then the volume of the whole
domain and the integral of `x` over it. On the unit cube the exact answers are 1 and 1/2, so the
program checks itself. The new ideas are how to reach the mesh data from GPU code and which elements
a loop must cover.

The program is `examples/distributed/unstructured/mars_step02_volume.cu`.

## Run it

```bash
python3 scripts/generate_hex_cube.py -n 16 --output cube16
mpirun -np 4 ./build/examples/distributed/unstructured/mars_step02_volume --mesh=cube16
```

```
volume = 1.000000000000000
integral of x = 0.500000000000000
```

The result is the same on any number of ranks, up to rounding in the last digits.

## The code

### Mesh data on the GPU

```cpp
const auto& conn = domain.getElementToNodeConnectivity();
ElementIntegral volume{{std::get<0>(conn).data(), ..., std::get<7>(conn).data()},
                       domain.getNodeX().data(), domain.getNodeY().data(), domain.getNodeZ().data(),
                       false};
```

The connectivity from [step 1](step-01-mesh.md) has eight columns, one per hexahedron corner; entry
`e` of column `k` is the local node id of corner `k` of element `e`. The node coordinates are three
arrays `x`, `y`, `z` indexed by local node id, built on the first call to `getNodeX()`: MARS creates
such data only when it is first asked for.

By default MARS keeps each node only as its space-filling-curve key and decodes the coordinates
from it, which rounds them to the curve's grid, a few parts in ten million of the box size. That is
fine for solving a PDE, but this step checks the volume to 1e-12, so it asks for the coordinates
exactly as read from the file:

```cpp
Domain domain(meshDir, rank, numRanks, true);   // true: store the original coordinates
```

All these arrays live in GPU memory. `.data()` gives raw device pointers, and the `ElementIntegral`
struct carries them to the GPU.

### The element computation

```cpp
struct ElementIntegral
{
    const uint64_t* corner[8];
    const double *x, *y, *z;
    bool moment;

    __host__ __device__ double operator()(size_t e) const { ... }
};
```

`operator()` runs on the GPU, once per element. It gathers the eight corners of element `e`, splits
the hexahedron into six tetrahedra around its diagonal from corner 0 to corner 6, and adds up their
volumes (the triple product divided by six). For elements with flat faces, like these, the split is
exact. With `moment = true` it weights each tetrahedron by `x` at its centroid, which integrates any
linear function exactly.

Corner order matters here: corners 0 to 3 are the bottom face counterclockwise, 4 to 7 the top face
above them. That is the order `generate_hex_cube.py` writes and MARS keeps.

### The loop: own elements only

```cpp
thrust::counting_iterator<size_t> first(domain.startIndex()), last(domain.endIndex());
local[0] = thrust::transform_reduce(thrust::device, first, last, volume, 0.0, thrust::plus<double>());
```

`thrust::transform_reduce` calls `volume(e)` for every `e` from `startIndex()` to `endIndex()` on
the GPU and sums the results. The range is the rank's own elements. Halo elements also have
connectivity and coordinates on this rank, but they are owned by another rank, which counts them
itself: including them here would add parts of the domain twice. Run with `first(0)` and
`last(held elements)` on four ranks to see the volume come out above 1.

### Summing over ranks

```cpp
MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
```

Each rank has the integral over its own piece; the sum over all ranks is the integral over the
domain. The release checks `marsReleaseStep02Volume_npN` and `marsReleaseStep02Moment_npN` require
both values to within 1e-12.

The pattern of this step, a functor over the element's corners run by Thrust on the rank's own elements,
is how you write your own operators on top of MARS: an element matrix, a residual, an error estimate.

## Try this

- Integrate `x*x` instead: the centroid rule is no longer exact, and the error shrinks as the mesh
  is refined (exact value 1/3).
- Generate the cube with `--scale 2`: the volume becomes 8.

Next: [step 3](step-03-poisson.md) solves a partial differential equation.
