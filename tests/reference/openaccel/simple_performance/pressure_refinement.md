# Experimental bounded SIMPLE pressure correction

Pressure recovery is disabled by default. `--pressure-refinement 1` explicitly
enables this experiment; the startup log records `pressure_refinement=0` or `1`.
All ranks must select the same mode before iterating. With recovery disabled, a
rejected pressure candidate remains a failure without attempting a correction.

OpenAccel's reference Hypre path does not perform this extra correction loop.
It must stay disabled for controlled OpenAccel comparisons. The snapshot,
startup and pressure-accuracy comparators reject logs with recovery enabled or
attempted, including historical logs that lack the selection line but contain
recovery records. Older logs without either remain eligible for their existing
provenance checks. The same equations do not guarantee the same finite-iteration
history when the linear-solver procedure differs.

When explicitly enabled, the distributed pressure solve attempts a correction
after a finite candidate fails either Hypre's explicit residual check or MARS's
owned-row check. This adds
no new PDE term, quadrature, boundary treatment, pressure gauge, relaxation or
time discretization. Momentum and already accepted pressure solves are unchanged.

For the assembled pressure-increment equation `A phi = b`, form `r = b - A phi`,
solve `A delta = r`, and test `phi + delta` against the original equation.
Both `phi` and `delta` have pressure units; `r` retains the continuity row's
units and scaling. The original matrix and RHS remain fixed throughout.

At most three corrections are allowed, with the original solve and corrections
sharing the original Krylov iteration budget. Each correction starts from zero,
uses relative tolerance 0.1, absolute tolerance zero and minimum iterations zero.
Those inner controls do not change the original acceptance target. All original
Hypre stopping and wrapper acceptance controls are restored before returning.
The existing prepared matrix and preconditioner are reused without another setup.

A correction must do positive work and strictly reduce the compensated residual.
Success then requires all of:

- Hypre's explicit residual for the original equation and target;
- MARS's ordinary owned-row residual after exchanging candidate ghost values;
- an upper bound on the compensated residual at that same target.

The last check includes dot-product evaluation error using the Dot2Err bound
from [Ogita, Rump and Oishi, Algorithm 5.8](https://www.tuhh.de/ti3/paper/rump/OgRuOi05.pdf).
It assumes ordinary IEEE round-to-nearest arithmetic with gradual underflow,
explicit FMA product remainders and no reassociation. Nonfinite/overflowed bounds
cannot pass. The norm uses an additional factor of two in its squared upper
bound to cover positive-reduction rounding within the adapter's checked integer
row-count range. This is an acceptance guard, not a lower bound proving that
a rejected tolerance is unattainable.

No progress, exhausted budget, uncertain arithmetic or failed original-system
checks still reject the solve. The best improving candidate is restored on
failure for diagnostics. Field vectors, compensation and candidate updates stay
on the GPU; MPI exchanges device data. Only scalar control/API results reach the
host. Correction buffers are allocated on the first failure and reused.

`[HypreGMRES] rejected:` records the rejected candidate. A subsequent
`[simple-pressure-refinement] ... accepted=1` records a verified correction.
The public diagnostic exporter keeps the rejection history and distinguishes
recovered candidates from terminal failures. Missing/invalid recovery records,
later rejections, application errors and concatenated runs cannot report success.

## Validation and limits

The local CPU Hypre 2.33.0 wrapper regression passes under ASan/UBSan. Its public
eight-row dyadic problem reproduces GMRES returning above the requested residual
target; correction reaches the exact known solution. GMRES/FlexGMRES control and
resource restoration, failure/budget checks, host MPI integration on 1/2/4 ranks,
and exact-rational dot-bound tests are covered. The dot test includes a case
whose compensated point residual incorrectly rounds to zero while its error
bound correctly prevents acceptance.

These results do not identify the private run's exact Hypre exit branch or prove
that its requested tolerance can be reached. CUDA compilation/execution of this
experimental integration remains pending. Recovery results cannot establish
parity of the unmodified OpenAccel solve procedure.
Earlier component GPU results in `../distributed_matrix/README.md` predate it.

## Alps GPU gate

Use the MARS uenv terminal and the existing `mars-v010-check/build-hypre` build.
The synthetic matrix gate invokes recovery explicitly. The channel parity gate
retains the default without recovery.
Fetch only once, with no other Git operation running in that checkout. The user
executes these commands; all results remain on capstor scratch.

```bash
(
set -euo pipefail
repo=/capstor/scratch/cscs/gandanie/git/mars-v010-check
test "$(git -C "$repo" branch --show-current)" = cstone
git -C "$repo" fetch origin cstone
git -C "$repo" merge --ff-only origin/cstone
cd "$repo/build-hypre"
run=$(mktemp -d /capstor/scratch/cscs/gandanie/simple-pressure-refinement-XXXXXX)
mkdir "$run/tmp"
export TMPDIR="$run/tmp" PYTHONDONTWRITEBYTECODE=1
printf 'Public results: %s\n' "$run"
cmake --build . --parallel 4 --target mars_distributed_matrix_cuda_gate mars_distributed_simple_cuda_gate mars_segregated_simple \
  2>&1 | tee "$run/build.log"
git rev-parse HEAD > "$run/mars-revision.txt"
sha256sum ./mars_distributed_matrix_cuda_gate ./mars_distributed_simple_cuda_gate \
  ./examples/distributed/unstructured/mars_segregated_simple > "$run/executable.sha256"
unset CUDA_LAUNCH_BLOCKING MARS_HYPRE_ABSTOL MARS_HYPRE_MINITER MARS_HYPRE_RESIDUAL_AUDIT MARS_SIMPLE_PRESSURE_AUDIT
export MARS_HYPRE_FLEXGMRES=0 MARS_HYPRE_SPMV_VENDOR=0 MARS_HYPRE_VERBOSE=0
for np in 1 2 4; do
  set +e
  srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node="$np" \
    --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
    "$HOME/affinity/bind_numa.sh" ./mars_distributed_matrix_cuda_gate \
    2>&1 | tee "$run/np$np.log"
  statuses=("${PIPESTATUS[@]}")
  set -e
  printf '%s\n' "${statuses[0]}" > "$run/np$np.exit"
  if (( statuses[0] != 0 || statuses[1] != 0 )); then exit 1; fi
done
srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=1 \
  --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
  "$HOME/affinity/bind_numa.sh" ./mars_distributed_simple_cuda_gate \
  --write-reference "$run/reference.bin" 2>&1 | tee "$run/reference.log"
for np in 1 2 4; do
  srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node="$np" \
    --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
    "$HOME/affinity/bind_numa.sh" ./mars_distributed_simple_cuda_gate \
    --reference "$run/reference.bin" --builder 1 2>&1 | tee "$run/simple-$np.log"
done
printf 'Public results: %s\n' "$run"
)
```

The matrix gate must report both `GPU pressure correction` cases passing: one
with caching, one without. It deliberately forces an inaccurate initial stop
on a synthetic diagonal system, then checks the original target, known solution,
ghost values, setup counts and restored controls through the production refiner.
The SIMPLE runs separately check that successful public-channel solves retain
their entry-by-entry parity. A passing gate is the prerequisite for preparing a
fresh private history capture bound to the rebuilt executable; old capture
provenance must not be edited to substitute that binary.
