# Conforming tetrahedral AMR checks

These checks cover the optional LEPP policy documented in
[Conforming-Tet-AMR.md](../../docs/Conforming-Tet-AMR.md). Native GPU validation
uses the installed MARS package. Run all three platforms against the same source
snapshot before considering this policy for release.

## CPU

The standalone geometry test needs only a C++20 compiler. The distributed tests
also need MPI and the pinned cornerstone headers. The legacy adapter test uses
an installed CPU MARS library with Kokkos and its AMR backend enabled.
The recipe pins `Kokkos_DIR` for both the library and the installed consumer;
the gate rejects a GPU Kokkos or a version other than 4.4.1.

`prepare.sh` builds that library and the tests from offline dependencies:

```bash
bash tests/amr/prepare.sh cpu "$PWD/build-amr-cpu" /path/to/cornerstone /path/to/kokkos
ctest --test-dir build-amr-cpu/checks --output-on-failure --no-tests=error -V \
  --output-junit "$PWD/build-amr-cpu/results.xml"
```

Set `CXX` to the intended C++ compiler and `MARS_AMR_BUILD_JOBS` to the build
parallelism. The geometry-only shortcut is:

```bash
cmake -S tests/amr -B build-tet-geometry
cmake --build build-tet-geometry
ctest --test-dir build-tet-geometry --output-on-failure
```

The complete CPU suite has nine CTests:

- Geometry: all 144 vertex-permutation/edge choices, 720 equal-length traversal
  orders, complete edge stars, signed volume, paired interior faces and affine transfer.
  All eight red-refinement children are also checked under all 24 parent permutations.
- Actual legacy Kokkos LEPP: the two-tet case produces five active tets and seven
  nodes after splitting the neighbour's longer edge first. The requested-edge
  control produces four tets and six nodes.
- Shared LEPP engine at one, two and four MPI ranks, each with direct and staged
  messages. Fixtures exercise changed local numbering, empty ranks, a 48-cell
  edge star, successive midpoint generations, constant/affine transfer, boundary
  ancestry and repeated partition changes. Canonical topology is compared with
  one rank. Sparse boundary matching tests overlapping labels, duplicate source
  records and lost-face rejection. One migration plan transfers three fields,
  excluding poisoned replicas and rejecting missing or conflicting owner values.
  Malformed connectivity, faces, geometry, midpoint identities and parentage are
  rejected collectively; split communicators are also tested.
  A PMPI send counter checks that disconnected local partitions send no edge,
  node or face records. Sparse-plan checks cover mixed local/remote records,
  reverse ordering and reuse, in direct and staged modes.

- Cornerstone layout on four CPU ranks with 32- and 64-bit keys. Shared tet
  representative nodes force deep SFC leaves. Every element slot must have one
  supplier, poisoned halo connectivity must be replaced, and node owners must
  hold complete stars. The old ten-retry limit fails this check for 64-bit keys.

CPU LEPP partition changes are fixture redistributions. They do not execute the
GPU-only device-data `ElementDomain` constructor or cornerstone GPU kernels.

## CUDA and HIP

Within the intended compiler/runtime environment, build and install the candidate:

```bash
bash tests/amr/prepare.sh cuda "$PWD/build-amr-cuda" /path/to/cornerstone
# On the AMD system instead:
bash tests/amr/prepare.sh hip "$PWD/build-amr-hip" /path/to/cornerstone
```

`MARS_AMR_GPU_ARCH` defaults to `90` for CUDA and `gfx942` for HIP. `CUDACXX`
and `HIPCXX` can select the GPU compiler. Cornerstone is supplied through CMake's
existing source-tree override; no download is needed. Tests link the installed
`Mars::mars` and `Mars::mars_unstructured` targets, including the installed
shared bisection header.

The native checks compile in separate CUDA/HIP files for mesh setup, migration,
direct LEPP, red refinement and manager integration. The CPU oracle and `main`
compile as ordinary C++. Test compilation is serial and keeps going after an
independent file fails, so one build records all failing units.

CUDA compilation uses `compile_cuda.py`. Each invocation keeps the Release flags,
records its command, exit status, elapsed time and peak child RSS under
`checks/compile-receipts/`, and has a ten-minute timeout. On Linux, each compiler
process inherits a 64 GiB virtual-address-space limit. This bounds individual
processes, not the sum of their RSS. Compiler core dumps are disabled so a crash
cannot fill scratch with an image of that address space. `MARS_AMR_COMPILE_MEMORY_GIB` and
`MARS_AMR_COMPILE_TIMEOUT` override these diagnostic bounds. The launcher does not
change the compiler's optimization or target architecture.

There are ten CTests per GPU platform: the host geometry test plus nine native
cases. The native matrix is one, two and four ranks with needed-node halos,
full halos and staged AMR messages. Each native case checks:

- The existing device sparse exchange over three variable-size rounds, before
  mesh setup. The unrefined cube is checked for positive cells, CPU/device SFC
  owner agreement, complete element-halo coverage and node presence on the owner
  before lazy halo construction.
- A two-tet device LEPP fixture with dependent midpoint generations and empty
  ranks; connectivity and boundary labels must match the CPU core.
- A partition-crossing longest-edge request on a 320-tet procedural cube,
  checked against CPU LEPP, positive volumes and paired interior faces.
- Constant and affine fields through actual midpoint birth batches.
- Forced field ownership changes and exclusion of poisoned replicas.
- Default red refinement, its P1 transfer, and a fresh device-data rebuild.
- Three asymmetric `AmrManager` refinement/repartition cycles, with poisoned
  ghost inputs, unique ownership, consistent shared coordinates, conformity,
  transferred fields and seven preserved face labels (one overlaps another).
- Coarsening marks leaving the direct refiner unchanged.

The MPI library must support device buffers even in staged AMR mode: that option
stages only the new LEPP/boundary/migration messages. Cornerstone and the existing
node halo retain their transport. Device selection uses node-local rank modulo
visible GPUs, or device zero when only one is visible. Logs report host, rank,
local rank, PCI identity and GPU name; verify placement before accepting a run.

## Cluster matrix

`cluster.sh` prepares and runs the suite inside an existing **two-node Slurm
batch allocation**. `mpi-launch.sh` uses one node for one rank and two nodes for
two/four ranks. It supplies `--export=ALL`, CPU binding and failure propagation.
It expects the CSCS uenv/default-view workflow; set `MARS_AMR_UENV` to the chosen
image. The user submits the allocation; none of these instructions authorize
an assistant to execute cluster commands.

For a frozen bundle, submit `batch.sh` with two arguments: the backend and the
absolute bundle directory. Slurm copies the script into its spool directory,
so the bundle must not be resolved from the batch script's own path. The local
regression checks this behavior and failure propagation:

```bash
python3 tests/amr/test_batch.py
```

The batch command is:

```bash
bash tests/amr/cluster.sh cpu /absolute/results/cpu /absolute/cornerstone /absolute/kokkos
# Or, in the corresponding GPU environment:
bash tests/amr/cluster.sh cuda /absolute/results/cuda /absolute/cornerstone
bash tests/amr/cluster.sh hip /absolute/results/hip /absolute/cornerstone
```

The frozen recipe uses `--mem=0` with the exclusive batch allocation and build
step. This alone did not solve the CUDA failure: job 5299797 had 849.61 GB on the
build node and reached 787.56 GB RSS before nvcc was killed. The split compilation
and compiler bounds isolated that failure.
Job 5303947 compiled the driver and mesh units but exposed unsupported extended
CUDA lambdas in constructors, private captured types, and a function with deduced
return type. Its two largest units reached the compiler limit. The current code
uses named device-callable operations for those paths. The separately compiled
CPU oracle now explicitly requests C++20. CUDA jobs 5304289 and 5305197 compile
and link all seven units successfully with these corrections.

Use a fresh result directory for each backend. Keep `build.log`, `ctest.log`,
`results.xml`, `ctest-exit.txt`, compiler versions, CMake caches/compile commands,
installed headers and binaries. A missing status file or nonzero process/CTest
exit is not a pass. Do not infer completion from a printed PASS line. Retain
failures; use the exact tested source/dependency manifest for any release claim.

## Numerical bounds and status

The affine field is `1 + 2*x - 3*y + 0.5*z`. Direct refinement allows `1e-12`;
constants allow `1e-13`. Each cornerstone rebuild may move a coordinate within
an SFC interval `h_d = box_length_d / (2^axis_bits_d - 1)`, so the manager test
accumulates `2*h_x + 3*h_y + 0.5*h_z` per rebuild. Midpoint averaging cannot
increase an existing maximum field error. Interior faces must still match
exactly. The small test oracle gathers geometry on the host; production closure,
face matching and migration do not gather a global mesh.

The frozen v6 candidate passed all three native gates on 2026-11-09:

| Platform | Job | Result |
| --- | --- | --- |
| Daint CPU | 5305220 | 9/9 CTests |
| Daint NVIDIA GH200 | 5305197 | 10/10 CTests |
| Beverin AMD MI300A | 729093 | 10/10 CTests |

All batch and CTest exits are zero, with no skipped tests. The source manifest is
`a7c8f8c4db50a830a64504899f97bd4ef58a0c7d9e969666a4e43af036a5b6f4`.
Rank receipts confirm one node for one rank, two nodes for two/four ranks and a
distinct GPU for each rank. Both GPU platforms pass needed/full node halos and
staged AMR messages, including all three manager refinement/repartition rounds.
This validates the one-rank halo accessor, red-child orientation and four-rank
halo-depth corrections on the tested paths. Earlier failures remain retained.

The final manager fixture contains 3315 tetrahedra on NVIDIA and 3309 on AMD.
These counts repeat across each platform's tested rank counts and modes; both
pass the geometric and field invariants. Cross-vendor identical topology is not
claimed. The cause of the count difference has not been isolated.

These gates do not establish coarsening, periodic/multiblock support, mesh-quality
bounds, geometric generation limits, Exodus file-ingestion coverage or scaling.
The tested dependency correction is published and pinned. The integrated release
checks below complete the validation matrix for this scope.

## Integrated release results

The integrated numerical/runtime source is commit
`2ee88f1fc84b2fdd634bbc7aaf01631780d07002`, with the same Cornerstone pin as above.
Commit `a4ece5961307c5f3e81de6a5657c30c12b675239` changes only the validation
runner, its tests and instructions. Results audited on 2026-11-09:

| Platform | Job | Result |
| --- | --- | --- |
| Daint CPU | 5307495 | 9 AMR checks plus red-tet geometry: 10/10 |
| Daint NVIDIA GH200 | 5307505 | Single-node release suite: 50/50 |
| Daint NVIDIA GH200 | 5310365 | Two-node release suite: 50/50, batch exit 0 |
| Beverin AMD MI300A | 729339 | Single-node release suite: 50/50 |
| Beverin AMD MI300A | 729814 | Two-node release suite: 50/50, batch exit 0 |

The release suite includes the existing unstructured solver/assembly checks,
red-tet orientation and single-rank halo access. Four-rank cases in the two-node
layout span both nodes. All expected cases ran, with no skips or failed tests.

Both original GPU jobs exited 1 after their single-node suites passed: the
auditor expected a GoogleTest footer that CTest had truncated from JUnit output.
The matching full logs establish that the named halo test ran and passed. The
corrected auditor rechecked those results, and the continuation jobs ran only
the missing two-node suites using the original executables. No numerical test
or tolerance changed; the original failure records are retained.

The continuation verifies source/dependency manifests and 37 original files,
including all 17 executables per platform, before and after execution. It
preserves test commands and properties while relocating output directories.
Executable hashes were first recorded after the original jobs completed, before
continuation; they were not captured during those original executions. Result
audits also checked exact inventories, artifact hashes and rank placement.

The integrated bundle manifest is
`ca2379b65dcccef25447e72aa0194fae86bcd7adeaaa5fc984404dc5c48fd492`;
its full source manifest is
`5c0b857cea04f017c9a0db16da361a6ac92acca06ae657ada8f1988f13aee526`.
The continuation manifest is
`0a1535a97ef8fcb770eae14098103ed74895bd563565f06429437b928eb6fb77`.
These are correctness results, not scaling measurements.

## Cornerstone halo-depth correction

CUDA job 5304289 compiled all seven native units successfully (peak child RSS
1.58 GiB) and passed seven of ten tests. All three four-rank cases failed before
refinement. The same failure is reproducible with Cornerstone's CPU domain:
the old ten-retry limit leaves a leaf spanning a rank boundary, with 21 allocated
halo elements that no rank supplies. Checking only known node stars missed
these unfilled slots.

The corrected dependency uses `maxTreeLevel<KeyType>` retries in `sync`,
`syncGrav` and `addHalos`. A failed layout is reduced as a boolean so positive and
negative error codes cannot cancel. `sync` and `syncGrav` throw collectively
before publishing an invalid layout. Halo requests are exchanged after the
layout converges. This changes no AMR decision or MPI message protocol.

MARS pins Cornerstone's `mars-extensions` commit
`c13b9a9a6c0f5897ec4f23283e4ff53716303c5e`. Its production header is identical to
the correction used by all three v6 jobs. Those frozen bundles retain their
original provenance: `2a95e63b0ca5454712ab83f27505055c078887dc` plus the recorded
`cornerstone-halo-depth.patch`. The published commit also adds a four-rank
32/64-bit regression; it fails with the old header and passes with the fix.
