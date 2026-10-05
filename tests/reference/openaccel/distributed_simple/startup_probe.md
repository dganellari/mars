# Fresh startup comparison

The geometry probes agree on the public cases. The remaining private snapshot
disagreement needs a fresh execution record and the first differing state. This
probe compares the initial nodal velocity/pressure and all 20 completed SIMPLE
iterations, by source node identity. It does not change a numerical operator.

`simple_startup_probe.py prepare` selects the original deck by the SHA-256 in
`case.json`, validates the existing mapping, requires the same mesh and zero
initial fields with no restart, and writes a new private deck. Only the mesh path
spelling, stopping limits and output controls change. Physical parameters,
discretization, boundary conditions and OpenAccel linear-solver settings remain
unchanged. Output uses solver values, not corrected boundary values. The pinned
OpenAccel source writes state zero during initialization when output is enabled.

MARS's optional `--snapshot-iterations 20` writes states 0..20 through its existing
distributed field writer. Default runs do no extra I/O. This is explicitly
sensitive diagnostic output: all fields and numerical logs must remain private
on capstor. Snapshot I/O is counted as output time, but these runs are not suitable
for performance measurements. A run that converges before state 20 is rejected as
incomplete evidence; the option does not override convergence.

The launch wrapper records the command it executes, executable hash, input hashes,
selected runtime environment and resolved `ldd` library hashes before and after
execution, exit status, and output hashes. The comparator verifies these records,
both rank counts, MARS's logged controls, every saved state and node mapping, and
the final snapshot against ordinary final field output. This is a fresh-launch
record, not an attestation of the compute process's loaded libraries or proof of
identical linear solvers. It does not retrofit provenance onto old runs.

Errors use velocity scale U and physical pressure scale rho*U^2, with no pressure
offset removed. The first state whose maximum velocity or pressure error exceeds
1e-5 is reported publicly using only an iteration number and fixed field labels.
Detailed errors, identifiers, counts, paths and fields stay in
`comparison-private.json` and must not be shared. Missing, changed, nonfinite or
incomplete evidence is a rejection, not a parity result. Passing 20 states is not
nonlinear convergence or pump validation. A mismatch locates a state/field, not
yet the responsible assembly or solve stage.

## Daint: OpenAccel terminal

Use the restored OpenAccel uenv/Spack environment and working Python dependencies.
No OpenAccel rebuild is needed. The wrapper launches the visible `srun` command
below; all new files go under the fresh capstor directory. The scratch pointer
passes that directory to the separate MARS terminal.

```bash
(
set -euo pipefail
root=/capstor/scratch/cscs/gandanie/git
scratch=/capstor/scratch/cscs/gandanie
git -C "$root/mars-v010-check" pull --ff-only
python3 -c 'import numpy, netCDF4, yaml'
umask 077
pair=$(mktemp -d "$scratch/simple-startup-XXXXXX")/pair
mkdir -p "$scratch/tmp"
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
probe="$root/mars-v010-check/scripts/simple_startup_probe.py"
python3 "$probe" prepare \
  --case "$scratch/simple-private-observed-3ll0nB/case/case.json" \
  --reference-dir "$root/OpenAccel/prgenv/pi-highres-01" --output-dir "$pair"
printf '%s\n' "$pair" > "$scratch/simple-startup-current.txt"
printf 'Private pair: %s\n' "$pair"
python3 "$probe" run --pair "$pair" --solver openaccel --ranks 4 \
  --executable "$root/OpenAccel-reference-updates-IxgJIp/source/build/openaccel-3D.exe" -- \
  srun --account=csstaff --time=00:15:00 --nodes=1 --ntasks-per-node=4 \
  --cpus-per-task=1 --cpu-bind=cores --export=ALL --kill-on-bad-exit=1
)
```

## Daint: MARS terminal, after the reference capture completes

Use the restored MARS CUDA/Hypre uenv. Build only the changed driver. The short run
normally returns 2 because it has not converged. The capture accepts 0/2, checks
all evidence, and returns 0 on a completed capture; other launch exits fail.
`--kill-on-bad-exit=0` prevents an expected rank exit 2 from killing its peers.

```bash
(
set -euo pipefail
root=/capstor/scratch/cscs/gandanie/git
scratch=/capstor/scratch/cscs/gandanie
cd "$root/mars-v010-check/build-hypre"
python3 -c 'import numpy, netCDF4, yaml'
cmake --build . --parallel 4 --target mars_segregated_simple
pair=$(cat "$scratch/simple-startup-current.txt")
test -f "$pair/reference/launch.json"
umask 077
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
unset CUDA_LAUNCH_BLOCKING MARS_OWNERSHIP MARS_HALO_FACTOR MARS_NODEHALO_ALLOW_INCONSISTENT
export MARS_HYPRE_SPMV_VENDOR=0 MARS_HYPRE_FLEXGMRES=0 MARS_HYPRE_ABSTOL=0
export MARS_HYPRE_VERBOSE=0 MARS_HYPRE_RESIDUAL_AUDIT=0 MARS_SIMPLE_PRESSURE_AUDIT=1
probe=../scripts/simple_startup_probe.py
if ! python3 "$probe" run --pair "$pair" --solver mars --ranks 4 \
  --executable "$PWD/examples/distributed/unstructured/mars_segregated_simple" -- \
  srun --account=csstaff --time=00:15:00 --nodes=1 --ntasks-per-node=4 \
  --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=0 \
  "$HOME/affinity/bind_numa.sh"; then
  if test -f "$pair/mars/run.exit"; then
    python3 ../scripts/simple_public_diagnostics.py \
      --log "$pair/mars/run.log" --exit-file "$pair/mars/run.exit" \
      --output "$pair/failure-public.json"
    cat "$pair/failure-public.json"
  fi
  exit 1
fi
status=0
python3 "$probe" compare --pair "$pair" --output "$pair/comparison-public.json" || status=$?
cat "$pair/comparison-public.json"
printf 'Share only: %s\n' "$pair/comparison-public.json"
exit "$status"
)
```

## Local checks

### Inspect a failed capture without launching again

The original wrapper could hide a nonzero launcher exit behind a second error
when the failed solver had not produced an Exodus file. It now preserves the
failure record before looking for successful output. Library preflight failures
also have separate fixed labels.

For an existing attempt, `inspect` reads local metadata and sanitizes the log.
It launches no solver and does not overwrite the attempt. Its library check
describes the current shell environment; it cannot reconstruct an unrecorded
earlier environment. Missing exit metadata does not prove that no job started.
Only the new public JSON may be shared.

An exit of 134 with no output is not a field-parity failure. The reference may
have aborted before a snapshot was written. OpenAccel's `errorMsg` throws a C++
exception; it does not print the `ERROR:` prefix used by the MARS diagnostic
parser. The inspector therefore also recognizes uncaught-exception messages and
reports fixed startup stages and error categories. Exception text, paths, part
names and numerical values are never exported. Categories are diagnostic hints,
not proven causes; an empty list does not exclude an unrecognized error.
Stages mean that at least one rank printed the marker, not that all ranks passed
that stage. Reinspect the saved attempt below; do not launch another pair yet.
The inspector also reports standard C++ exception classes and a small allowlist
of public source filenames when they appear in the error. Unknown classes become
`other`; private filenames, paths, line numbers and exception text remain local.
The `mesh_ready` marker is at the end of `mesh::read`, before `mesh::setup` registers
geometric fields. Later messages can be buffered, so the marker alone does not
prove the failing call. Field-registration and master-element signatures provide
more specific evidence when available.

For an unclassified exception, `--reference-source` builds a message index from
Git revision `0d69041ba1afda63e9e4328d9e0d9834bba37756`, restricted to public
OpenAccel/Nalu source directories. It never reads untracked files or the current
working-tree contents. Literal error fragments are matched locally; the public
JSON contains only candidate source paths and line numbers from that revision.
These are possible message origins, not a stack trace or a proven cause. Dynamic
messages and diagnostics from separately installed libraries may not match.

```bash
(
set -euo pipefail
root=/capstor/scratch/cscs/gandanie/git
pair=/capstor/scratch/cscs/gandanie/simple-startup-lh05AV/pair
git -C "$root/mars-v010-check" pull --ff-only
summary=$(mktemp -d /capstor/scratch/cscs/gandanie/simple-startup-check-XXXXXX)/public.json
python3 "$root/mars-v010-check/scripts/simple_startup_probe.py" inspect \
  --pair "$pair" --solver openaccel \
  --executable "$root/OpenAccel-reference-updates-IxgJIp/source/build/openaccel-3D.exe" \
  --reference-source "$root/OpenAccel" \
  --output "$summary"
cat "$summary"
printf 'Share only: %s\n' "$summary"
)
```

### Tests

Synthetic tests cover an exact history, an injected pressure mismatch at step 3,
missing snapshots, altered input/run identities, launch failure, nonfinite fields,
missing nodes and preserved controls. The existing field reader tests cover
distributed CSVs, permuted Exodus global IDs and duplicated ghost consistency.
The output gate writes changing snapshots on one and four host MPI ranks under
ASan/UBSan. CUDA compilation and these paired Daint executions remain user-run.

```bash
PYTHONPATH=scripts python3 -m unittest test_simple_startup_probe test_simple_snapshot_compare test_prepare_simple_deck
```
