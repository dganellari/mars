# First-step comparison with the reference AMG choices

The public `simple-reference-pressure-check-F87prv` report confirms that the
PETSc momentum predictor and both assembled equations match at the first step.
The pressure increments differ, although both satisfy the same saved reference
residual target on the common unscaled system. All twenty printed OpenAccel
pressure targets were met in the original-control history, without a cap or an
ambiguous rounded comparison. Those logged values are not independently
recomputed residuals for every step.

The normal GPU pressure profile deliberately substitutes PMIS for HMIS and
l1-Jacobi for the direct coarse solve. `prepare --reference-amg` removes those
two substitutions while keeping the original pressure and momentum targets,
restart, iteration cap, other supported AMG controls and PETSc momentum profile.
It reuses the saved OpenAccel first step and runs one new MARS step. Recovery
remains disabled. The default device profile is unchanged.

## Algebra and scope

The equation remains `A phi = b`, with the nodal Tet4 CVFEM physical pressure
increment `phi` in Pa. SIMPLE still applies `p += alpha_p phi` and
`u -= d grad(phi)`. Boundary contributions, shifted quadrature, volume scaling,
density, viscosity, pseudo-time and relaxation are unchanged. Only construction
and application of the pressure preconditioner change. Both independent residual
checks still use the original target `max(atol, rtol ||b||)`.

The executable option is `--pressure-amg-mode reference`; it requires an explicit
pressure profile and a one-step audit. The preparation wrapper also requires
`--both-linear-solvers`. It rejects history mode, unsupported settings and changed
baseline evidence before launch. There is no tolerance change or new acceptance
rule. A rebuilt MARS executable is permitted; the captured dependency identities,
rank count, launcher and solver environment remain checked.

**This mode permits CPU work inside Hypre.** In Hypre 2.33,
[HMIS](https://github.com/hypre-space/hypre/blob/v2.33.0/src/parcsr_ls/par_coarsen.c)
clones device matrix structure to host for its Ruge phase and returns the markers
to device; its subsequent PMIS phase follows the execution policy.
[Coarse solver 9](https://github.com/hypre-space/hypre/blob/v2.33.0/src/parcsr_ls/par_gauss_elim.c)
stages the coarse matrix and vectors on host for Gaussian elimination. These are
library-imposed transfers, not a GPU performance optimization. No library rebuild
is needed; this mode requires Hypre 2.33 or newer.

Keeping these choices does not make CPU and GPU solvers identical. Hypre versions,
ordering, ownership and implementation paths can still differ, including the
device PMIS phase inside HMIS. The experiment tests the two explicit substitutions;
it cannot guarantee field agreement or establish nonlinear convergence.

Each rank's post-setup settings must retain HMIS 10 and effective coarse solver 9.
Hypre can replace the coarse solver when coarsening stops early. Such a fallback
fails this comparison, because it would not test the requested direct solve.
`recorded_pressure_settings_verified: true` certifies the selected settings;
`recorded_gpu_pressure_settings_verified: false` in this mode records the explicit
permission for host work, rather than a failed settings check.

## Validation

The host-adapted production wrapper passes its numeric-refresh regression with
real sequential Hypre 2.33.0 and 3.1.0 under AddressSanitizer and
UndefinedBehaviorSanitizer. The added strong-diffusion cases exercise GMRES and
FlexGMRES with a multilevel hierarchy, actual HMIS/direct settings, a prepared
host coarse solve, two numeric updates and a known independent solution.
The existing device-profile and rejection cases remain included.

The Python comparison tests cover the unchanged baseline, original targets,
mixed PETSc/Hypre profile, host-work disclosure before launch, changed saved
arguments and rejected runtime substitutions. C++ option tests reject use outside
the first-step audit. These local checks do not validate CUDA execution or MPI.

The focused CUDA gate below checks the actual controls, more than one level,
host coarse preparation, independently known solutions, exchanged ghosts and
numeric refresh for both Krylov methods. It runs on 1, 2 and 4 GPUs before the
private comparison. A failed gate stops the recipe. CUDA/MPI execution and the
new field comparison remain unvalidated until those user-run results arrive.

## Run in the MARS terminal only

Use the working MARS uenv/venv and existing PETSc-enabled build. Run the Git update
in this terminal only. This rebuilds two MARS targets, not PETSc or Hypre.
The public gates launch directly with `srun`. The subsequent `run` command
executes the saved private launcher once; do not wrap it in another `srun`.
OpenAccel is not relaunched.

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
repo="$scratch/git/mars-v010-check"
test "$(git -C "$repo" branch --show-current)" = cstone
git -C "$repo" fetch origin cstone
git -C "$repo" merge --ff-only refs/remotes/origin/cstone
python3 -c 'import numpy, netCDF4, yaml; print("Python dependencies OK")'
umask 077
mkdir -p "$scratch/tmp"
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
cd "$repo/build-hypre"
cmake --build . --parallel 8 --target mars_segregated_simple mars_distributed_matrix_cuda_gate
run=$(mktemp -d "$scratch/simple-reference-amg-XXXXXX")
printf 'Results: %s\n' "$run"

if ! python3 "$repo/scripts/simple_pressure_profile.py" prepare \
  --baseline-pair "$scratch/simple-first-step-MPQGzD/pair" \
  --defaults-probe "$scratch/simple-hypre-defaults-PrkNXY/probe" \
  --executable "$repo/build-hypre/examples/distributed/unstructured/mars_segregated_simple" \
  --output-dir "$run/pair" --both-linear-solvers --reference-amg \
  --output "$run/preparation-public.json"; then
  cat "$run/preparation-public.json"
  exit 1
fi

for np in 1 2 4; do
  MARS_HYPRE_SPMV_VENDOR=0 srun --account=csstaff --time=00:05:00 \
    --nodes=1 --ntasks-per-node="$np" \
    --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
    "$HOME/affinity/bind_numa.sh" ./mars_distributed_matrix_cuda_gate --reference-amg \
    2>&1 | tee "$run/gate-$np.log"
done

if ! python3 "$repo/scripts/simple_pressure_profile.py" run \
  --pair "$run/pair" --output "$run/launch-public.json"; then
  cat "$run/launch-public.json"
  exit 1
fi
if ! python3 "$repo/scripts/simple_pressure_profile.py" compare \
  --pair "$run/pair" --detail-dir "$run/private-comparison" \
  --output "$run/comparison-public.json"; then
  cat "$run/comparison-public.json"
  exit 1
fi

python3 - "$run/comparison-public.json" <<'PY'
import json, sys
with open(sys.argv[1]) as stream:
    report = json.load(stream)
keys = ("comparison_status", "failed_check", "pressure_amg_mode",
        "hypre_host_staging_permitted", "recorded_pressure_settings_verified",
        "original_pressure_targets_preserved", "original_momentum_targets_preserved",
        "gpu_adaptations", "identical_linear_solvers_verified")
summary = {key: report.get(key) for key in keys}
comparison = report.get("comparison", {})
summary["first_step"] = {key: comparison.get(key) for key in
    ("all_snapshots_match", "first_differing_iteration", "first_differing_fields",
     "first_differing_stage", "first_step_stage_matches", "pressure_common_system_checks")}
print(json.dumps(summary, indent=2))
PY
printf 'Share only: %s\n' "$run/comparison-public.json"
)
```

Agreement would show that retaining these two choices is sufficient for this
first-step comparison. If the pressure increment still differs with verified
settings and passing original residual targets, the substitutions are not the
complete explanation. In either case, keep the comparison bounded; completion
alone does not justify another long nonlinear run.
