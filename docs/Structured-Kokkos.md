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
Views pass. The corrected MARS correctness gates pass; setup timing gates remain
open as detailed below.
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

The subsequent native builds and tests pass with CUDA 12.9 on NVIDIA GH200
and HIP/ROCm 7.2.3 on AMD MI300A, using Kokkos 4.4.01. Both direct and host-staged
MPI pass all 11 structured and 3 installed-consumer CTests at 1/2/4 ranks on one
node, with a distinct GPU per rank. Each mode covers 168 GoogleTest
case/rank-count combinations, including full-Hilbert direct partitioning and
shallow manager copies. CUDA jobs are 5259527 (direct) and 5259528 (staged);
HIP jobs are 720226 (direct) and 720227 (staged). Their JUnit results, runtime
receipts and retrieved artifact hashes were independently checked against the
same frozen source snapshot, whose manifest SHA256 is
`14eb4f9f7030292f51a392bcd195c6957a7da6a3ff50e421b6c603258ac9976b`.

These tests do not establish old PDE application correctness, staggered
discretizations, periodic DOF identification, arbitrary
polynomial orders, empty partitions, large-scale performance or GPU race freedom.

## Release gate

All applicable checks must pass on **CPU, AMD GPUs and NVIDIA GPUs** before
release. Use the same candidate source on all three platforms and retain its
hashes, build configuration and runtime receipts. A pass on one platform does
not cover another. Required tests that fail, are skipped or have no evidence
keep the release blocked.

The pending DOF optimization sorts identical owned/local candidate lists once,
while retaining independent final key arrays. It passes local Serial/direct and
OpenMP/staged checks, but needs fresh native validation on all three platforms.
The table below records the preceding candidate's evidence; those passes do not
validate this new change. The current suite adds minimum valid partitions and
contains 171 case/rank-count combinations.

| Required check | Linux CPU | AMD HIP / MI300A | NVIDIA CUDA / GH200 |
| --- | --- | --- | --- |
| Native build, 11 structured CTests and 3 installed-consumer CTests at 1/2/4 ranks | Pass, Serial/OpenMP and both MPI modes | Pass, both MPI modes | Pass, both MPI modes |
| MPI across at least two nodes, including ownership, partition and halo checks | Pass, all four configurations | Pass, both MPI modes | Pass, both MPI modes |
| Bounded mesh setup and repeated halo performance checks | Pass, all four configurations | Pass, both MPI modes | Blocked: ten setup timing rejections; MPI abort resolved |
| Kokkos-OFF build and existing applicable release regressions | Core build passes | Pass, 48 tests in each launch layout | Pass, 48 tests in each launch layout |

Results checked on 2026-11-05 use structured source manifest
`1518a67c6efe04f184adf6c1487fb9dd7392f95ded360235745d20b4e8a58270`.
Preparation jobs 721518 (CPU), 721520 (HIP) and 5262344 (CUDA) pass all 11+3
single-node CTests in each configuration. Linux CPU uses GCC 13.4 and Kokkos
4.4.01, with four OpenMP threads. Each configuration covers 168 GoogleTest
case/rank-count combinations. GPU regression jobs 721522 (HIP) and 5262346
(CUDA) pass all 48 existing release tests in each single-node and two-node
launch layout, with Kokkos OFF.

Validation jobs 721898 (CPU), 721899 (HIP) and 5264408 (CUDA) reuse those
prepared binaries and pass the 11+3 CTest matrix in every configuration, with
two/four-rank cases spanning two nodes. The earlier allocation and batch-path
launch failures are resolved. Runtime receipts, driver hashes, JUnit results
and raw benchmark samples were independently checked.

HIP passes all 180 setup/halo performance comparisons. CPU Serial/direct,
Serial/staged and OpenMP/staged pass; OpenMP/direct rejects one of 90 comparisons:
four-rank single-node Morton Quad4 mesh setup has candidate/reference medians
2.084/1.810 ms, above the predefined limit of 2.041 ms.

CUDA/staged rejects eight of 90 comparisons. CUDA/direct completes its
single-node samples, with four of 54 comparisons rejected, then the reference
binary aborts in MPI during its first two-node, two-rank Morton Quad4 warmup.
No candidate multi-node benchmark runs in that mode. The largest observed setup
slowdown is four-rank single-node Hilbert Hex8 rectangular mesh generation:
4.793/2.688 ms candidate/reference with direct MPI, and 4.963/2.889 ms staged.
These earlier failures are retained. Later allocator validation resolves the MPI
abort, while setup timing regressions still block release.

The follow-up diagnostics (CPU 722241, CUDA 5265708) retain all failures and
independently verified 338 log/trace hashes. CPU passes all 18 focused control
comparisons; the subsequent full confirmation is recorded below.
CUDA reproduces 12 focused setup rejections and four failed launches. The MPI
error is `MPI_Irecv` reporting an OFI tagged receive failure (`Invalid argument`).
One failing trace allocates a 64 KiB receive buffer, accepts same-node receives,
and fails on the first remote-node receive. The original Kokkos build enables
`cudaMallocAsync`; this motivated the isolated allocation comparison below.

CPU confirmation 722701 passes all 90 OpenMP/direct comparisons, using the exact
original binaries and limits. Combined with the other three passing CPU
configurations, this closes the current CPU gate without deleting the earlier
rejection.

CUDA preparation 5267897 and validation 5267898 use the same MARS implementation
with Kokkos asynchronous allocation disabled (source/runner manifest
`018c79d35bc4bd2d761f98a40a3147628889b3220771a31271cca22288dad135`).
Both modes pass all 11+3 CTests on one/two nodes. All 100 benchmark invocations
complete with correct results and no MPI failure, but direct/staged reject six
and four of their 90 timing comparisons respectively. All halo timing comparisons
pass. Direct multi-node setup has 35–60 ms candidate medians in several cases;
the reference has similarly large delays in different cases. Their timing cause
is still under investigation.

The independent 60-launch allocation probe passes all ordinary-CUDA and corrected
View cases. All six cross-node raw `cudaMallocAsync` cases fail in OFI receive,
as do four original-View cases above Kokkos's allocation threshold. Smaller
original Views and every single-node case pass. CUDA pointer attributes confirm
the pooled/device allocation types. Payload, padding, synchronization and buffer
lifetime checks separate this MPI allocation issue from MARS mesh code.
The review independently verifies 271 downloaded artifacts and reparses all
150 CPU/CUDA benchmark invocations and allocator outcomes.

The private partition-buffer lifetime trial (5269246/5269247) passes correctness
but still rejects eight setup comparisons. It is not adopted. Its 24 diagnostic
profiles add global GPU fences and do not reproduce the large multi-node delays,
so they cannot establish the cause of the uninstrumented stalls. The unfenced
follow-up below uses the exact original binaries; production timing acceptance
remains blocked. The trial audit independently verifies 338
downloaded artifacts, 100 benchmark invocations and all 24 profiles.

Unfenced diagnostic5269661 passes all126 invocations and3 activation preflights;
609 downloaded artifact hashes and the reconstructed timers were independently
checked. Three of126 focused control comparisons reject. The repeatable
single-rank Morton Quad4 DOF slowdown is concentrated in CUDA temporary allocation
and free calls around the two sorts. The large multi-node stalls do not recur in
these untouched controls on a different node pair; their cause remains unresolved.
The new full validation uses the earlier failing CUDA node pair and retains
unfenced diagnostic traces after the uninstrumented gate.

CPU validation includes Kokkos Serial and OpenMP execution, with direct and
host-staged communication paths. Both use host-accessible buffers on CPU;
they do not establish GPU-aware MPI. The MARS serial backend stays excluded.
Both GPU platforms require direct device MPI and host-staged MPI, including
the multi-node checks. The single-node results above do not cover inter-node
communication.

Measure setup and repeated halo exchanges on each platform with fixed workloads,
rank/thread counts and source versions. Define the comparison and acceptance
criteria before running; include correctness checks with the timings. A
correctness pass alone does not establish acceptable performance or scaling.

For regressions, validate the supported CPU/core configuration on CPU and the
existing unstructured release suite on both GPUs. The unstructured GPU backend
has no CPU release configuration. Keep the backend's optional, OFF-by-default
behavior and verify supported build combinations. The cluster kit prepares CPU/GPU capability and installed-package checks,
multi-node checks and paired setup/halo measurements; a separate runner audits
the existing 48-test unstructured GPU release suite. CPU and NVIDIA performance
gates still need resolution and passing cluster receipts.
AMR remains a separate release scope.

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

Full Hilbert squares/cubes with power-of-two side lengths directly generate each
rank's consecutive key interval and partition offsets on device. This path skips
sorting, splitter searches and key redistribution. Rank sizes differ by at most
one cell, including uneven rank counts. Map construction and collective input/map
validation remain. Tests compare against Cartesian reference keys and exercise
three-rank subcommunicators inside the four-rank gate.

The general case uses the same distributed algorithm for both curves. Each rank generates and sorts
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
