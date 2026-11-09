# Optional structured Kokkos backend

The structured backend provides distributed Cartesian Quad4 and Hex8 meshes with
Morton or Hilbert ordering, nodal DOF numbering and field halo exchange. Enable it
explicitly with `MARS_ENABLE_KOKKOS=ON`. It needs MPI, including for one-rank runs.
It is independent of the unstructured cornerstone backend.

## Build and use

Use a Kokkos installation with the execution space you want to test. The checked
dependency version is 4.4.01. This CPU configuration leaves the other backends and
application examples disabled:

```sh
cmake -S . -B build-structured \
  -DCMAKE_BUILD_TYPE=Release \
  -DKokkos_DIR="$Kokkos_DIR" \
  -DMARS_ENABLE_KOKKOS=ON -DMARS_ENABLE_MPI=ON \
  -DMARS_ENABLE_DISTRIBUTED_BACKEND=ON \
  -DMARS_ENABLE_UNSTRUCTURED=OFF -DMARS_ENABLE_AMR_BACKEND=OFF \
  -DMARS_ENABLE_SERIAL_BACKEND=OFF -DMARS_ENABLE_KOKKOS_KERNELS=OFF \
  -DMARS_ENABLE_EXAMPLES=OFF -DMARS_ENABLE_TESTS=ON \
  -DMARS_STRUCTURED_TEST_RANKS='2;4'
cmake --build build-structured --target mars morton_test hilbert_test \
  mesh_manager_test dof_handler_test mesh_ghost_layer_test --parallel 4
ctest --test-dir build-structured -L structured --output-on-failure -V
```

`Kokkos_DIR` is the directory containing `KokkosConfig.cmake`, not its source tree.
For CUDA or HIP, build Kokkos with `Kokkos_ENABLE_COMPILE_AS_CMAKE_LANGUAGE=ON`
and the matching `Kokkos_ENABLE_CUDA=ON` or `Kokkos_ENABLE_HIP=ON`. MARS uses native
CMake GPU languages and rejects a mismatched or wrapper-only Kokkos package.
For CUDA with GPU-aware MPI, configure Kokkos itself with
`Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=OFF`. Kokkos 4.4 defaults to pooled CUDA
allocations above 40 kB; some MPI transports cannot register that memory.
This changes the allocator, not the device location of Views or the MPI data
path. See the [Kokkos compatibility note](https://kokkos.org/kokkos-core-wiki/known-issues.html).
The validation dependency builder sets this option explicitly. Native allocation
comparisons confirm this failure on the tested GH200/Cray MPI stack: pooled
allocations fail across nodes, while ordinary CUDA allocations and corrected
Views pass.
Keep the C++ host compiler separate from the GPU compiler:

| Setting | NVIDIA CUDA | AMD HIP |
| --- | --- | --- |
| MARS backend | `MARS_ENABLE_CUDA=ON` | `MARS_ENABLE_HIP=ON` |
| GPU compiler | `CMAKE_CUDA_COMPILER=/path/to/nvcc` | `CMAKE_HIP_COMPILER=/path/to/rocm/llvm/bin/clang++` |
| Architecture example | `CMAKE_CUDA_ARCHITECTURES=90` for Hopper | `CMAKE_HIP_ARCHITECTURES=gfx942` for MI300 |

Select the architecture of the actual GPU and use the same compiler and
architecture settings when building Kokkos, MARS and the installed consumer.
The [cluster validation kit](../tests/structured/cluster/README.md) exports pinned
dependency sources and builds all three offline in private directories. It
provides separate CUDA and HIP launch recipes and records source, dependency,
compiler, binary and test evidence.

`MARS_ENABLE_RDMA=ON` passes device buffers directly to MPI and requires GPU-aware
MPI. `MARS_ENABLE_RDMA=OFF` stages both sends and receives through host memory.
Run both modes on each GPU platform; a staged pass does not validate direct device
MPI. The kit disables CUDA UVM so host-accessible allocations cannot hide a device
communication failure.

Every suite and the installed consumer check the execution space, run a device
kernel and an MPI exchange, and print a receipt for each rank. GPU runs require
one distinct physical GPU per rank on each host. Explicit Kokkos device selection
and launcher visibility are respected; otherwise the test assigns visible devices
by local MPI rank. MIG instances sharing a PCI bus ID are outside this gate.
`MPIEXEC_EXECUTABLE`, `MPIEXEC_NUMPROC_FLAG` and `MPIEXEC_PREFLAGS` configure CTest's
launcher. The kit runs from an allocated batch shell so its one-rank tests also
execute on a compute node.

The installed-consumer gate instantiates the public header, creates a mesh and
DOFs, and checks an affine field on owned and ghost nodes:

```sh
cmake --install build-structured --prefix "$PWD/install-structured"
cmake -S tests/structured/install-consumer -B build-structured-consumer \
  -DCMAKE_PREFIX_PATH="$PWD/install-structured"
cmake --build build-structured-consumer --parallel 4
ctest --test-dir build-structured-consumer -L structured-install --output-on-failure -V
```

Link downstream targets to `Mars::mars`; the package exports C++20, MPI and Kokkos
requirements. Kokkos remains OFF by default for existing unstructured users.

## Capability checks

| CTest label | Executables | Checked data |
| --- | --- | --- |
| `structured-sfc` | `morton_test`, `hilbert_test` | Encoding, decoding and curve ordering |
| `structured-mesh` | `mesh_manager_test` | Every Cartesian cell occurs once globally; exact ghost-cell set; explicit coordinates/connectivity; rectangular, thin and full-curve grids; balanced sorted partitions; local-only key generation |
| `structured-dof` | `dof_handler_test` | Q1/Q2 scalar and two-component numbering; one owner per physical lattice node/component; inverse maps and consistent ghost IDs |
| `structured-halo` | `mesh_ghost_layer_test` | Scalar and two-component fields; selected and all-field gathers/additive scatter; owned-data access; periodic cell fields; reusable exchanges, channel isolation and retained device views |
| `structured-install` | `mars_structured_consumer` | Installed public API, dependency resolution, global Q1 ownership count and owned-to-ghost values |

The mesh, DOF and halo suites run at one, two and four ranks by default. The rank
matrix is configurable through `MARS_STRUCTURED_TEST_RANKS`; an empty value keeps
only one-rank tests. CTest timeouts turn a communication deadlock into a failure.
Each field in `DM` stores all local components in interleaved node/component order.
Set a runtime block size before constructing its data manager, and keep it fixed
for that manager's lifetime. Selecting tuple fields for an exchange leaves the
other fields unchanged.

The checks exposed missing dependent-template arguments in the public data
manager, an invalid Hilbert corner-key bound during explicit point generation,
and a DOF count-exchange hang caused by assuming symmetric ownership between
peers. Point generation sorts and deduplicates local cell corners on device.
Boundary construction sorts sparse peer/cell pairs, avoiding a dense
cell-by-rank allocation. Both GPU backends use the same Kokkos scan, sort and
compaction paths as the CPU tests. Memory and layout follow the selected Kokkos
execution space; host-staged MPI covers both directions and collectives.
Morton and Hilbert DOF construction now share the local/ghost candidate algorithm,
without allocating predicates over the global key range. Invalid mesh sizes,
degree-scaled key overflow and failed SFC map construction are rejected before
later communication.

## Validation boundary

The local 2026-11-04 validation uses macOS arm64 and Open MPI 5.0.5 with Kokkos
4.4.01. The final AppleClang 21 Debug builds using its Serial execution space pass
with direct and host-staged MPI; an LLVM 19/libomp Release build passes with two
OpenMP threads and direct MPI. Each configuration passes 11 CTests and 168
GoogleTest case/rank-count combinations, for 504 combinations total. This counts
each MPI case once, not once per process. Installed-consumer tests pass at one,
two and four ranks in both communication modes. Kokkos's Serial execution space
is separate from the excluded MARS serial backend.

The halo suite includes eight two-component cases and checks selected tuple
exchange, field allocation, owned-data access and runtime block size. Explicit
mesh checks verify reference corner order and positive affine Jacobians. Runtime
preflight also checks view/element wrapper lifetimes, nested environments, map
failure detection and short-valued MPI operations.

CUDA compatibility syntax checks cover the nine library sources, five suites
and consumer in both modes using LLVM 19, CUDA 12.4 headers and Kokkos 4.4.01.
These use a private macOS compatibility configuration; they are not native
compilation, linking or GPU execution.

Native results on CPU, AMD and NVIDIA are in [Release status](#release-status).

These tests do not establish old PDE application correctness, staggered
discretizations, periodic DOF identification, arbitrary
polynomial orders, empty partitions, large-scale performance or GPU race freedom.

## Release status

The structured backend is supported. Commit `409a3c87` passed every correctness
gate on CPU, AMD and NVIDIA on 2026-11-08, from one frozen kit (SHA256
`2b97166aa35310f4c1fcedc1cd04093e8b228f8373ee6d4d3523616dcc4aca50`):

| Check | Linux CPU | AMD HIP / MI300A | NVIDIA CUDA / GH200 |
| --- | --- | --- | --- |
| 11 structured and 3 installed-consumer CTests at 1/2/4 ranks; the 2- and 4-rank cases span two nodes | Pass: Serial and OpenMP, direct and staged MPI | Pass: direct and staged MPI | Pass: direct and staged MPI |
| Kokkos-OFF build | Pass | — | — |
| Existing unstructured GPU release suite (48 tests, Kokkos OFF) | — | Pass | Pass |
| Setup and halo timing against the reference snapshot (diagnostic) | 360/360 | 180/180 | 176/180 |

Jobs (preparation/validation): CPU 727102/727103, HIP 727104/727105, CUDA
5295845/5295846. The unstructured suite ran on 2026-11-05 (HIP 721522, CUDA
5262346). Later commits change only structured sources, tests and documentation.

For this optional structured release, the timing comparison is reported as a
diagnostic. The original acceptance rule and its four CUDA rejections are
retained; those rejections are accepted as non-blocking, rather than reclassified
as passing measurements. The reference is snapshot `7c0684ad` of the same backend.
The threshold is 1.10 times the reference median plus 50 microseconds for setup
or 5 microseconds for halo exchange.

The four rejected comparisons all use direct MPI. These are maximum-rank elapsed
medians over five paired rounds:

| Layout | Case | Phase | Reference (ms) | Candidate (ms) | Change |
| --- | --- | --- | ---: | ---: | ---: |
| Two nodes, 2 ranks | Hilbert Quad4 full | DOF setup | 42.425 | 49.321 | +16.3% |
| Two nodes, 4 ranks | Morton Quad4 rectangular | Mesh setup | 41.290 | 51.117 | +23.8% |
| Two nodes, 4 ranks | Hilbert Hex8 rectangular | Mesh setup | 36.171 | 54.682 | +51.2% |
| One node, 4 ranks | Morton Hex8 rectangular | DOF setup | 7.949 | 8.881 | +11.7% |

Candidate and reference ranges overlap, and all four cases pass the timing rule
with host-staged MPI. These observations do not establish noise as the cause or
rule out a performance regression. The rectangular paths are unchanged by the
last full-Hilbert patch, but the comparison includes earlier changes since
`7c0684ad`. The four timing exceptions therefore remain unresolved. All halo
comparisons pass; this release makes no general performance-parity or scaling claim.
The original CUDA aggregate status remains `FAILED_OR_TIMING_REJECTED` because
its direct-MPI timing stage returned `TIMING_REJECTED`; its correctness tests pass.

Earlier candidates, their timing rejections and the pooled CUDA allocation
diagnosis (2026-11-04 to 2026-11-06) are in the history of this file.
For regressions, rerun the [cluster validation kit](../tests/structured/cluster/README.md)
on all three platforms with one candidate source. A pass on one platform does
not cover another. AMR is a separate release scope.

## Device views and layout lifetime

`Mesh`, `DofHandler` and `DM` are host owners. Capture their `device_view()` results
in kernels. Mesh, element, DOF and field device views are nonpolymorphic values;
they contain device allocation handles and geometry/numbering metadata, without
an MPI context, host mirrors or communication state. Managed Kokkos handles retain
the storage after the original owner is destroyed. A view is a layout snapshot,
not a copy of the field values: writes through a field view remain shared.

```cpp
const auto dofs = handler.device_view();
const auto field = data.get_dof_data<0>();
Kokkos::parallel_for("initialize", dofs.get_dof_size(), KOKKOS_LAMBDA(mars::Integer i) {
    double point[3] = {};
    dofs.get_local_dof_coordinates(i, point);
    field(i) = dofs.is_owned(i) ? point[0] : 0.0;
});
```

Finish setup before constructing dependent handlers and managers. Rebuilding a
host mesh replaces its allocations; existing handlers and views retain the old
layout. Construct a new handler and data manager for the new mesh or component
layout. Do not mutate shared topology arrays while a dependent view is in use.
The host mesh/handler still borrows its context; keep that context alive, unmoved
and non-null for host operations. Device snapshots do not borrow the context.

Copy construction of `DM` and `UserData` remains shallow and requires no MPI
collective. Copies share field storage and host communication plans. The plan
owner uses a managed host Kokkos view; device copies retain the handle without
accessing the host plan. The compatibility test exercises legacy manager captures
as well as copies that outlive the original manager. Prefer `device_view()` for
new kernels to keep host control state out of their captures.

## Partition construction

Full Hilbert squares/cubes use the same distributed partition algorithm as
rectangular meshes. The direct-interval shortcut was removed after it triggered
setup regressions in the native comparison. Tests still compare against Cartesian
reference keys and exercise three-rank subcommunicators inside the four-rank gate.

Both curves use the same distributed algorithm. Each rank generates and sorts
only its Cartesian share, with at most `ceil(global_cells / ranks)` keys. Parallel
binary searches and reductions of rank-sized counts select exact SFC split points.
A sparse MPI distributed graph discovers the destinations' incoming peers, and
keys move directly to those ranks. Each final partition differs in size by at most
one cell. No rank gathers the global cell-key array or scans the enclosing curve's
empty coordinate range.

Partition metadata remains proportional to the number of ranks, and selection
uses at most the key width in collective rounds. This is a memory bound, not a
measured strong-scaling result. Sparse boundary construction sorts and deduplicates
`(peer, cell)` and `(peer, owned DOF)` pairs. DOF scratch is proportional to boundary
cell/DOF candidates, rather than owned DOFs times all neighboring ranks.

## Reusable halo exchange

A `DM` constructs an isolated exchange channel collectively on its mesh
communicator. Its plan stores directed peer segments and local boundary/ghost
indices. Each selected field's device buffers, MPI request storage and optional
host staging buffers are allocated once and reused. The plan has one exchange in
flight at a time and posts all selected fields before returning from `start`.

```cpp
data.start_gather<0>();
// Independent work may use other storage while communication is outstanding.
do_independent_work();
data.wait();
```

`gather_ghost_data<...>()` and `scatter_add_ghost_data<...>()` remain synchronous
convenience methods. `start_scatter<...>()` followed by `wait()` returns received
contributions in boundary order. Those returned buffers are borrowed until the
next exchange of the same field; copy them if they must live longer. Gather
updates ghosts, and additive scatter accumulates remote contributions at their
unique owner. Unselected fields are unchanged. Peers must select matching fields
and directions; two independent managers can communicate concurrently because
they use separate channels. Copies of a manager share one plan: starting an
exchange through one copy makes it active through every copy, and any copy may
wait for it. Do not change an in-flight field's layout or access
its receive values before `wait()`.

Only plan/channel construction is collective. `start` and `wait` use sparse
point-to-point requests and add no all-rank collectives. Existing one-shot free
helpers retain their blocking, shared-communicator behavior; for repeated
operations use `DM` or an explicit `HaloExchangePlan<...>`. Callback helpers now
post communication before invoking the callback. The blocking cell-field
`UserData` API also reuses its field buffers and peer requests, on its existing
shared communicator; calls on that communicator must remain consistently ordered.
A `UserData` copy made before the first exchange acquires its ghost-buffer handles
when it exchanges. Copies made after that exchange share the existing ghost views.

Destroy plans and every manager copy before MPI and Kokkos finalize. Destruction
of the last host copy completes outstanding requests before releasing buffers and
the channel. Host plan initialization/destruction is sequential, so MPI cleanup
stays on the calling thread. Actual MPI
progress during independent work depends on the MPI implementation; the start/wait
contract does not establish a performance gain or measured overlap.
