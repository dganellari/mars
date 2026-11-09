#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 7 ]]; then
    echo "Usage: job.sh SNAPSHOT RESULT_DIR cuda|hip ARCH CXX GPU_COMPILER direct|staged [run.py launcher options]" >&2
    exit 2
fi
if [[ -z ${SLURM_JOB_ID:-} || ${SLURM_STEP_ID:-batch} != batch ]]; then
    echo "Run this as an sbatch batch script, not from a resource-holding srun step." >&2
    exit 2
fi
snapshot=$1
result_dir=$2
backend=$3
arch=$4
cxx=$5
gpu_compiler=$6
transport=$7
shift 7

if [[ -e $result_dir ]]; then
    echo "Use a fresh result directory for every job." >&2
    exit 2
fi
mkdir -p "$result_dir"
export PYTHONDONTWRITEBYTECODE=1
kit="$snapshot/source/tests/structured/cluster"
python3 "$kit/build_deps.py" --source "$snapshot/dependencies/kokkos" \
    --work "$result_dir/kokkos" --backend "$backend" --arch "$arch" \
    --cxx "$cxx" --gpu-compiler "$gpu_compiler" --jobs "${BUILD_JOBS:-8}"
python3 "$kit/run.py" --source "$snapshot/source" \
    --work "$result_dir/$backend-$transport" --backend "$backend" --transport "$transport" \
    --arch "$arch" --cxx "$cxx" --gpu-compiler "$gpu_compiler" \
    --kokkos-dir "$result_dir/kokkos/install/lib/cmake/Kokkos" \
    --gtest-source "$snapshot/dependencies/googletest" --launcher "$(command -v srun)" \
    --launcher-arg=--export=ALL --jobs "${BUILD_JOBS:-8}" "$@"
