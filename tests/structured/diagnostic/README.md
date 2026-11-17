# Structured setup and MPI diagnostics

This kit investigates structured setup performance with the exact
prepared reference and candidate benchmark binaries. It builds only an observer
library and small MPI activation probes. It does not rebuild MARS, change its
algorithms or relax the release limits. AMD's passing gate is not repeated.

The completed sort-once matrix passes correctness and installed-consumer checks
on CPU, AMD and NVIDIA. AMD passes all timing comparisons; one CPU and eight
CUDA setup comparisons remain rejected. CUDA's pooled-allocation MPI failure is
resolved by the checked Kokkos dependency recipe. The
partition-buffer lifetime trial did not resolve them and is not adopted.

`unfenced.py` targets the remaining CUDA timings. It reuses the original binaries
from preparation5267897 and schedules 70 control invocations plus 56 observations
at all seven affected transport/layout/rank combinations. Controls have five
alternating reference/candidate rounds; MPI-only and CUDA/Kokkos tracing each
have two. All six workloads and correctness checks remain active in every run.

`mpi-quiet` buffers MPI events without printing per-call messages or loading the
Kokkos tool. `host` also loads Kokkos callbacks and CUDA runtime wrappers, with
global fencing disabled. The wrappers forward the original calls and return
codes without adding CUDA operations. Host callbacks measure host time, including
existing waits; they do not measure asynchronous GPU kernel duration. Activation
preflights verify the actual native interposition before interpreting traces.

Timing observers now preserve the original `MPI_Waitall` status argument.
The earlier observer, including job5275809, replaced `MPI_STATUSES_IGNORE` with
a status array even in `mpi-quiet` mode. This can select a different MPI wait
path, so the absence of large stalls in those traces does not establish why
the uninstrumented run stalled. It is a diagnostic confound, not a proven cause.
The original acceptance samples remain valid and their failures are retained.

Both modes record the existing `MPI_Wtime` return values. `setup_trace.py` checks
all 528 timer calls per rank and reconstructs the benchmark's reported max-rank
mesh, DOF and mean-halo times exactly. Per-phase raw calls, sizes and buffer
addresses locate stalls without inferring phases from output order. Nested
CUDA/MPI/Kokkos events overlap and their totals must not be added together.
Failed preflights and invocations are retained while other modes still run.
Neither traced samples nor a passing focused control subset close the full
release gate. The site bundle supplies pinned inputs and submission scripts.

The earlier `diagnose.py` campaign below remains available to reproduce its
historical CPU/CUDA investigation. Its fenced profiles hide the large CUDA
multi-node stalls and therefore cannot establish their uninstrumented cause.

The site bundle pins the source manifests, preparation receipts, installed
dependencies and benchmark hashes. Inputs are checked before and after the run.
The active allocation and device visibility come from the current job. Batch
scripts use an explicit shared-scratch location because Slurm relocates scripts.

CPU runs the failing OpenMP/direct configuration at four ranks on one node:
five alternating reference/candidate control pairs, then one profiled pair.
CUDA runs direct and staged modes. Direct MPI first runs both binaries at two
and four ranks across two nodes, retaining MPI failures and continuing to the
other binary. Controls use five pairs at four ranks on one node for direct MPI,
and on one/two nodes for staged MPI. Each CUDA mode also gets one profiled pair
at four ranks on one/two nodes. There are 12 CPU and 42 CUDA benchmark launches,
plus a two-rank MPI activation probe and one-rank profiled benchmark preflight
per platform. Both observers must activate before the campaign. Each launch has a
five-minute Slurm step limit and a separate 330-second client timeout.

All benchmark launches still run the original six workloads, warmup, changing
affine fields, owned-cell/DOF counts and halo correctness checks. A reference
failure cannot suppress the scheduled candidate diagnostics. Successful control
runs are compared using the original formula and preserved raw samples. These
focused controls do not replace the full release matrix.

The PMPI observer preserves message buffers, counts, status arguments and return
codes by default. `MARS_DIAG_RETAIN_STATUSES=1` explicitly enables replacement
of ignored `MPI_Waitall` statuses for an error investigation. Timing audits
reject this mode and older traces without an explicit status-policy receipt.
When statuses are ignored, aggregate errors remain visible; individual request
errors are reported only when a status array is available. On error it
prints the operation, error class/string, rank, peer, tag and count before the
application aborts. The native preflight intentionally sends to an invalid rank
on an errors-return communicator and requires per-rank observer receipts. This
rejects a missing or inactive interposer before the benchmark campaign starts.

Kokkos profiling records kernel, copy, allocation and fence labels; MPI records
operation durations. Kernel timing requests Kokkos global fences, so it includes
device completion. Raw events retain ordering and allocation sizes for mapping
back to the setup code. Totals cover the original six-case process; they are
not separate per-case or nonoverlapping phase totals. Fences/regions can overlap
other recorded events, so do not add all categories into a whole-run estimate.
CPU affinity receipts are retained for every Linux rank, including controls.

Profiled execution changes synchronization and incurs observer overhead. Its
timings are diagnostic only and are never used as release acceptance samples.
`DIAGNOSTIC_COLLECTED` means all planned launches were attempted, including any
retained failures; it does not mean the backend passed. Missing observer events,
truncated traces or incomplete commands remain explicit errors in the results.

The MPI wrappers can be checked locally with the explicitly linked probe.
Runtime injection on macOS may be suppressed, while Kokkos plugin loading still
works; that is not evidence for Linux PMPI injection. Native jobs therefore
require activation of the ordinary, unlinked probe through `LD_PRELOAD` before
running benchmarks. Local host profiling does not validate GPU execution.

`ctest -R trace_waitall_contract` checks forwarding and return codes against a
controlled PMPI implementation, including ignored/explicit statuses, empty
waits, errors and opt-in retention. It does not test native MPI transport or
establish the cause of the cluster timing failures.
