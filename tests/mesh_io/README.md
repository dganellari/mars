# Partitioned Exodus input

The Exodus and binary readers return the same element slices, sorted global node
IDs, local connectivity and coordinates. ElementDomain performs its existing SFC
redistribution afterwards. Exodus input uses independent netCDF hyperslab reads;
a parallel-netCDF build is not required. Coordinate requests coalesce nearby
IDs with at most four values requested per needed node and 8 MiB of scratch.
Physical disk traffic can differ because of filesystem caching and netCDF/HDF5
chunking. Metadata storage scales with block count, not element count.

When file-domain side-set storage is enabled on all ranks, `ElementDomain`
reads disjoint slices of the boundary face lists during construction. Each
reader processes at most 4,096 faces per batch, plus their referenced
connectivity and bounded coordinate ranges. Shared coordinates/connectivity may
be requested by more than one reader; the element/side list entries themselves
are read once globally. Names and small block metadata are read on every rank.

After volume redistribution, the GPU compacts the owned/ghost node keys into
spatial prefixes. Sparse MPI registers these prefixes in a distributed directory,
then routes face batches to candidate ranks. There is no global mesh-key copy,
all-gather or all-to-all directory. Only routing metadata is downloaded from the
GPU; netCDF face records already originate on the host. Incoming scratch scales
with delivered batches, including coarse-prefix candidates, rather than being a
fixed 4,096-face allocation. The GPU performs the exact local/block lookup.
Persistent boundary storage contains local node hits and complete oriented Tet4
faces, with file-entry order and repeated faces preserved.

Coincident nodes are matched by key and block when block identities are enabled;
conforming material blocks share nodes by default. Empty local sets retain their
names. The one-face test exercises readers with no entries. Getters make no MPI
calls: construction and file-domain resync finish distribution first, and all
ranks must agree whether to enable side-set storage. AMR device domains continue
to receive transferred boundary sets. Invalid input aborts distributed setup
rather than stranding ranks in communication.

The explicit `readExodusSideSetsTet4/Hex8` APIs still return complete boundaries
for callers asking for global output, using bounded read scratch. Legacy callers
of these APIs do not become distributed automatically. The separate raw-mesh API
is unchanged. No solver equations, quadrature or boundary values change here.

Build and run on a host with MPI and netCDF:

```sh
cmake -S tests/mesh_io -B build-mesh-io
cmake --build build-mesh-io --parallel 4
ctest --test-dir build-mesh-io --output-on-failure
```

Set `NETCDF_INCLUDE_DIR` and `NETCDF_LIBRARY` if they are not on CMake's search
path. Add `MARS_MESH_IO_CUDA=ON` with `CMAKE_CUDA_ARCHITECTURES=90`, or
`MARS_MESH_IO_HIP=ON` with `CMAKE_HIP_ARCHITECTURES=gfx942`, to compile and run the
same checks through the GPU deduplication path. Use the site's MPI launcher and
GPU binding when running `mars_exodus_partition FIXTURE_DIRECTORY` directly.
The MPI tests include 8 ranks to cover ranks with no elements. Set
`CORNERSTONE_INCLUDE_DIRS` to the matching Cornerstone `include` directory to
enable the actual sparse routing checks; it is required for the native domain gate.

The real netCDF fixtures cover Tet4/Hex8, multiple blocks including a null block,
32/64-bit connectivity variables, split/packed coordinates, sparse node IDs,
and missing optional element totals. Tests compare every returned array against
an independent fixture oracle and the binary reader. Instrumented netCDF calls
check local connectivity reads, coordinate request amplification and file closure
on exceptions. Malformed shapes, IDs and totals must fail. A large duplicate list
exercises GPU deduplication in CUDA/HIP builds. The no-netCDF build must reject
Exodus input explicitly.

`mars_exodus_side_sets FIXTURE_DIRECTORY` exercises the same streamed reader and
lookup on the host or selected GPU backend. It checks Exodus face winding, global
API parity, post-redistribution node/ghost subsets, coincident block identities,
conforming material blocks without duplicate boundary entries,
missing/inactive/locally empty sets, sparse coordinates, schema/ID refusals, and
scratch/read bounds. Its boundary list exceeds two batches. Routing checks use
instrumented netCDF reads and a global entry-coverage oracle to reject replicated
or missing list reads, then compare local nodes and oriented faces after sparse MPI.

After installing the candidate MARS GPU library, add `MARS_MESH_IO_DOMAIN=ON`
and its installation to `CMAKE_PREFIX_PATH`. This builds `mars_exodus_domain`,
run with a fixture directory and `conformal` or `blocks`, a public Hex8/Tet4
ElementDomain test with single and multiple blocks. Each mode uses a
separate process because the identity option is cached. The driver explicitly
disables or enables `MARS_BLOCK_NODE_IDENTITY` before constructing a domain.
Separate block identities use SFC-decoded coordinates: Hex8 does not support
combining them with stored original coordinates. The other cases retain stored
coordinates. It checks local and ghost boundary labels, Tet4 face winding/order,
and a rank-local getter against geometric oracles after actual SFC redistribution. It requires one GPU per rank. The standalone
reader checks alone do not validate ElementDomain or AMR.

The conforming-AMR validation must additionally run its imported Tet4 three-round
refinement/repartition/field-transfer gate against the same installed candidate.
The file-domain loader now prepares its local face ancestry in the same distributed
pass as node labels; the refinement/transfer cycle remains a separate native gate.
