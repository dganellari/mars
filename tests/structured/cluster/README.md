# Structured Kokkos release validation

The CUDA dependency recipe explicitly disables `Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC`
for MPI compatibility. This option belongs to the Kokkos build, so changing it
requires rebuilding Kokkos, MARS and the installed consumers. It keeps device
buffers and direct GPU-aware MPI enabled. Results for the current source are in the
[release status](../../../docs/Structured-Kokkos.md#release-status).

This kit builds and checks the optional structured backend using CPU Kokkos
Serial/OpenMP, NVIDIA CUDA or AMD HIP. The MARS serial and bisection backends
remain excluded. A separate regression runner checks the existing unstructured
GPU release configuration with Kokkos disabled.
Each run covers the five suites at their registered one/two/four-rank matrix
(11 CTests, 171 case/rank-count combinations) and the installed consumer at
one/two/four ranks (3 CTests). One GoogleTest case is counted once per rank count,
not once per process.

Run `direct` and `staged` separately for each backend. `direct` passes device
buffers to MPI; `staged` copies communication buffers through host memory. A
staged pass establishes GPU execution with host-staged MPI, not GPU-aware MPI.

The release workflow has two phases. A one-node preparation builds dependencies,
passes capability/install checks and prepares candidate/reference benchmark
binaries. A dependent two-node validation reuses those binaries, audits the
multi-node capability/install tests and measures performance. Neither a prepared
build nor the structured result alone establishes full release readiness; the
existing unstructured GPU release suite must also pass. See the
[release status](../../../docs/Structured-Kokkos.md#release-status).

## CPU and distributed release workflow

`run.py --backend serial` and `--backend openmp --threads 4` run the same 11+3
capability/install matrix as CUDA/HIP. CPU builds omit `--gpu-compiler` and `--arch`.
The CPU receipt proves its execution space and thread count, and distinguishes
the direct/staged code paths without claiming device-buffer MPI. One-rank tests
also use the supplied MPI launcher. `build_deps.py` accepts these CPU backends.

From an allocated batch shell, prepare one platform with `release.py`:

```sh
python3 /SNAPSHOT/source/tests/structured/cluster/release.py \
  --phase prepare --backend cpu --source /SNAPSHOT/source \
  --reference /REFERENCE/source --work /PREPARATION \
  --cxx /HOST_CXX --launcher /PATH/TO/srun --threads 4 --jobs 8 \
  --launcher-arg=--cpu-bind=cores
```

Use `--backend cuda --arch 90 --gpu-compiler /CUDA/bin/nvcc`, or
`--backend hip --arch gfx942 --gpu-compiler /ROCM/llvm/bin/clang++`, for GPUs.
Preparation covers both transport modes; CPU covers both Serial and OpenMP.
It also checks a Kokkos-OFF CPU core build. Add the site's explicit image/view and
binding flags through repeated `--launcher-arg=...` arguments.

After successful preparation, run this inside a two-node allocation with four
ranks, at least four CPU cores per rank and, for GPUs, distinct devices per rank:

```sh
python3 /SNAPSHOT/source/tests/structured/cluster/release.py \
  --phase validate --backend cpu --source /SNAPSHOT/source \
  --reference /REFERENCE/source --prepared /PREPARATION \
  --work /VALIDATION --threads 4
```

`mpi_launch.py` launches one rank on one node and distributes the two/four-rank
checks over two nodes. Receipts must name at least two different hosts; a
single-node launch cannot satisfy this gate. Binary, installed dependency and
source checksums must still match preparation. The launcher uses the same image
and binding arguments recorded at preparation.

Only the test backend, thread count and MPI transport settings are reused from
the preparation environment. Slurm allocation identifiers and GPU visibility
must come from the current job. Validation and performance receipts record
that environment and the driver hashes. A corrected driver directory can run
against an unchanged snapshot and prepared build; keep `run.py`, `release.py`,
`performance.py` and `mpi_launch.py` together. This avoids rebuilding validated
binaries when only orchestration changes.

## Performance comparison

The reference is the separately frozen, previously native-validated `7c0684ad`
snapshot. This measures the subsequent structured changes against that version;
it does not measure the whole redesign against the original legacy code.
The same benchmark source is compiled against both installed libraries.

Six fixed Q1 workloads cover Morton/Hilbert Quad4 and Hex8 rectangles plus full
Hilbert squares/cubes. Each process warms all workloads before measurement.
Every measured halo iteration changes the owned affine field, poisons ghosts,
and verifies all resulting values; global cell and unique-DOF counts must match.
Five paired rounds alternate reference/candidate order at 1/2/4 ranks on one node
and 2/4 ranks over two nodes. Timings are maximum-rank host elapsed intervals with
Kokkos fences, measured separately for mesh setup, DOF setup and repeated gather.
They are not GPU-only kernel timings or a scaling claim.

The predeclared acceptance limit is candidate median <= 1.10 times reference
median plus 50 microseconds for mesh/DOF setup or 5 microseconds per halo exchange.
The absolute allowance prevents tiny timings from dominating the decision.
Missing/invalid correctness receipts fail the run; a timing regression returns
`TIMING_REJECTED`, with raw samples retained. Other configurations still run.
The runner reports PASS only when every comparison passes. For optional
structured support, the four documented CUDA timing rejections are accepted as
non-blocking; their original `TIMING_REJECTED` status and raw samples remain.
Overlapping sample ranges alone do not establish noise or rule out a regression.
A nonzero aggregate result requires checking correctness separately from timing;
a correctness failure still blocks release. See the measured exceptions and
scope of this decision in the [release status](../../../docs/Structured-Kokkos.md#release-status).

## Existing GPU release regressions

`export_release.py` starts from the structured manifest's committed base and
overlays only the structured snapshot. It excludes uncommitted AMR work, private
documents and geometry; the shipped Poiseuille mesh is the sole included data
fixture. Export Cornerstone at `c13b9a9a6c0f5897ec4f23283e4ff53716303c5e` and the same
pinned GoogleTest into its dependencies directory. `regression.py` builds the
existing Kokkos-OFF GPU drivers using GPU-enabled Hypre, netCDF and NumPy. It
requires all 50 release tests: the original 48, including flow/Poiseuille, plus
red-tetrahedron orientation and one-rank halo-access regressions. It rejects
skipped or missing cases and requires the halo GoogleTest to execute its case.
`--multi-node` repeats this suite with
multi-rank runs spanning two nodes. Keep these receipts separate from the
structured capability/performance results.

## Prepare an offline source snapshot

Use clean source checkouts of these exact public dependencies:

| Dependency | Version | Git revision |
| --- | --- | --- |
| Kokkos | 4.4.01 | `15dc143e5f39949eece972a798e175c4b463d4b8` |
| GoogleTest | release-1.12.0 | `15460959cbbfa20e66ef0b5ab497367e47fc0a04` |

The exporter requires those revisions and a clean working tree. It copies only
their Git-tracked files, with no `.git` directory or untracked build outputs.
For MARS it copies the current tracked and untracked implementation under a
small source allowlist. Meshes, internal notes, other backends, worktrees and
unrelated documentation are excluded. Every copied file gets a SHA256 entry.
Do not edit the resulting snapshot: changed, missing or extra files fail the gate.

```sh
python3 tests/structured/cluster/export.py \
  --source "$PWD" \
  --kokkos-source /absolute/path/to/kokkos-4.4.01 \
  --gtest-source /absolute/path/to/googletest-1.12.0 \
  --destination /absolute/path/to/new-structured-snapshot
```

`source-manifest.json` records the MARS base revision and actual copied file
hashes, including local changes. Each dependency has its own revision and hash
manifest. Transfer this directory using the normal reviewed transfer workflow.
No network fetch runs in the build or test job.

## Allocation and environment

Load the site's working C++ compiler, MPI, Python 3, CMake >=3.24
and Ninja environment before submission; GPU jobs also need CUDA or ROCm. The kit does not guess module names,
accounts, partitions, containers or compiler paths. Use a host C++ compiler for
`CXX`, CUDA `nvcc` for the CUDA compiler, and ROCm `llvm/bin/clang++` for the HIP
compiler. Do not use `nvcc_wrapper` or `hipcc` as `CXX`: the build uses native
CMake CUDA/HIP languages. For CUDA the selected CXX is also set explicitly as
`CMAKE_CUDA_HOST_COMPILER` in Kokkos, MARS and the installed consumer.

The direct mode needs GPU-aware MPI. On Cray MPICH this may require
`MPICH_GPU_SUPPORT_ENABLED=1` in the loaded environment. Preserve any required
site transport settings; the result records that variable and
`MPICH_GPU_IPC_ENABLED`. The preflight actually sends a Kokkos device buffer
through an MPI ring and checks it in a device kernel before issuing a receipt.
Staged mode performs the same check through host buffers.

Request four whole GPUs for the four-rank gate. Within a host, each rank must
report a different physical PCI bus ID. A launcher may expose only one assigned
GPU per rank; logical device zero is then valid on every rank. MIG instances
sharing a PCI bus ID are deliberately not accepted. A single node with four
available GPUs is the smallest allocation. If it has fewer, allocate four ranks
across the required nodes and supply the site's GPU binding arguments.

Use the batch shell directly. Do not run `job.sh` inside an `srun` step holding
the resources: its nested two/four-rank CTests need to create their own steps.
The script rejects that nesting. Native GPU CTests launch even the one-rank cases
through MPI, so MPI initialization occurs inside a scheduler step on a compute
node. The batch shell only orchestrates those steps.

The following is an example submission, not a site-specific resource recipe.
Replace all capitalized values and the allocation flags with the site's valid
values. Native dependency and MARS compilation happen inside this bounded job.
The first command tests GH200/H100 architecture 90; the second tests gfx942.

```sh
sbatch --account=ACCOUNT --nodes=1 --ntasks=4 --ntasks-per-node=4 \
  --gpus-per-node=4 --cpus-per-task=8 --time=01:00:00 --export=ALL \
  --output=/RESULT_PARENT/structured-cuda-%j.log \
  /SNAPSHOT/source/tests/structured/cluster/job.sh \
  /SNAPSHOT /NEW_RESULT_DIRECTORY cuda 90 /HOST_CXX /CUDA/bin/nvcc direct

sbatch --account=ACCOUNT --nodes=1 --ntasks=4 --ntasks-per-node=4 \
  --gpus-per-node=4 --cpus-per-task=8 --time=01:00:00 --export=ALL \
  --output=/RESULT_PARENT/structured-hip-%j.log \
  /SNAPSHOT/source/tests/structured/cluster/job.sh \
  /SNAPSHOT /ANOTHER_NEW_RESULT_DIRECTORY hip gfx942 /HOST_CXX /ROCM/llvm/bin/clang++ direct
```

Replace `direct` with `staged` and use another fresh result directory for the
host-staged checks. Supported pinned architecture mappings are CUDA 80, 86, 90
and HIP gfx90a, gfx942. Set `BUILD_JOBS` in the submission environment to change
the default eight build processes. Extra arguments after the transport are
forwarded to `run.py`, for example `--launcher-arg=--cpu-bind=cores`. The kit
always includes `srun --export=ALL`; pass site-specific GPU binding flags using
additional `--launcher-arg=...` arguments.

## Build and run separately

`job.sh` performs these two steps. They can also be used in another scheduler's
batch shell with an explicit launcher, or with an already suitable private
Kokkos installation. The dependency builder enables native CMake GPU language
support, the requested GPU backend, Serial host support and C++20; its compiler
commands and configuration log are retained.

```sh
python3 /SNAPSHOT/source/tests/structured/cluster/build_deps.py \
  --source /SNAPSHOT/dependencies/kokkos --work /NEW_KOKKOS_BUILD \
  --backend cuda --arch 90 --cxx /HOST_CXX --gpu-compiler /CUDA/bin/nvcc

python3 /SNAPSHOT/source/tests/structured/cluster/run.py \
  --source /SNAPSHOT/source --work /NEW_CUDA_DIRECT_RESULTS \
  --backend cuda --transport direct --arch 90 \
  --cxx /HOST_CXX --gpu-compiler /CUDA/bin/nvcc \
  --kokkos-dir /NEW_KOKKOS_BUILD/install/lib/cmake/Kokkos \
  --gtest-source /SNAPSHOT/dependencies/googletest \
  --launcher /ABSOLUTE/PATH/TO/srun --launcher-arg=--export=ALL
```

For HIP set `--backend hip --arch gfx942` (or `gfx90a`) and its native HIP compiler.
For an MPI launcher use `--launcher /path/to/mpiexec --rank-flag=-n`; add its
binding options as separate `--launcher-arg=...` values. The tool never opens a
remote shell or submits a job. `--plan-only` validates offline inputs and writes
commands, but reports `PLANNED_NOT_VALIDATED`, never `PASS`.

## Acceptance and receipts

`result.json` is `PASS` only if configuration, build, installation and every
required CTest succeeded. The runner rejects zero tests, missing/duplicate tests,
skips, filtered GoogleTest suites, the wrong execution space or transport,
missing rank receipts or shared physical GPUs. CPU jobs require the requested
Serial/OpenMP execution space and thread count. Ordered runtime
receipts report backend, MPI ranks, local rank, physical device and Kokkos version.
Ordered test receipts require the complete suites: Morton 4, Hilbert 5, mesh 14,
DOF 17 and halo 23 cases per process, all passing with no skips.
CTest runs one case at a time to avoid competing GPU allocations. Inherited
`GTEST_*` options and diagnostic structured-runtime overrides are cleared.

Each result retains exact command arrays, compiler/version logs, MPI compiler
identification, CMake caches, CTest XML/logs, source/dependency checksums and
built executable/install-tree checksums. Only a short allowlist of environment
variables is recorded. The source and dependencies are checked again at the end.
Failure leaves a nonpassing `result.json` and the logs produced so far. Keep the
Kokkos build result and the `cuda-direct`, `cuda-staged`, `hip-direct` or
`hip-staged` result together.

Local tests of the kit use synthetic CTest output and command fixtures:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 tests/structured/cluster/test_kit.py -v
bash -n tests/structured/cluster/job.sh
```

Those checks establish runner behavior, not native CUDA/HIP compilation or GPU
execution. Hardware results are still required before a release capability claim.

Build/test children unblock only inherited `SIGCHLD` before exec. The parent mask
and all other child signal masks are preserved. Receipts record whether the parent
blocked it; this avoids the confirmed batch-environment CMake startup hang.
