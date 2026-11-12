#!/usr/bin/env bash
set -euo pipefail
kit=$1
results="${kit}-results-${SLURM_JOB_ID}"
mkdir "$results"
trap 'code=$?; printf "%s\n" "$code" > "$results/batch-exit.txt"' EXIT
export PYTHONDONTWRITEBYTECODE=1
export MPICH_GPU_SUPPORT_ENABLED=1
export MARS_AMR_COMPILE_MEMORY_GIB=64
export MARS_AMR_COMPILE_TIMEOUT=900
python3 -B "$kit/verify.py" "$kit" > "$results/kit-verify.log"
cp "$kit/manifest.json" "$results/source-manifest.json"
printf 'job=%s nodes=%s\n' "$SLURM_JOB_ID" "$SLURM_JOB_NODELIST" > "$results/allocation.txt"
srun --export=ALL --nodes=1 --ntasks=1 --cpus-per-task=16 \
    bash "$kit/build.sh" "$kit" "$results" > "$results/build.log" 2>&1
status=0
python3 -B "$kit/run.py" "$kit" "$results" > "$results/matrix.log" 2>&1 || status=$?
python3 -B "$kit/verify.py" "$kit" >> "$results/kit-verify.log"
exit "$status"
