#!/usr/bin/env bash
set -euo pipefail

[[ $# -ge 3 && $# -le 4 && -n ${SLURM_JOB_ID:-} ]] || {
    echo 'Usage inside a two-node batch allocation: cluster.sh cpu|cuda|hip WORK CORNERSTONE_SOURCE [KOKKOS_SOURCE]' >&2
    exit 2
}
backend=$1
work=$2
here=$(cd "$(dirname "$0")" && pwd)
[[ ! -e $work/ctest-exit.txt ]] || { echo "Use a fresh result directory." >&2; exit 2; }
mkdir -p "$work"
work=$(cd "$work" && pwd)
export MARS_AMR_MPI_LAUNCHER="$here/mpi-launch.sh"
export MPICH_GPU_SUPPORT_ENABLED=1 OMP_NUM_THREADS=1
[[ $backend != cpu ]] || export MPICH_GPU_SUPPORT_ENABLED=0
export MPICH_GPU_IPC_ENABLED=0
unset MARS_AMR_REUSE_DOMAIN
cmake --version > "$work/cmake-version.txt"
"${CXX:-mpicxx}" --version > "$work/cxx-version.txt"
if [[ $backend == cuda ]]; then "${CUDACXX:-nvcc}" --version > "$work/gpu-compiler-version.txt"; fi
if [[ $backend == hip ]]; then "${HIPCXX:-hipcc}" --version > "$work/gpu-compiler-version.txt"; fi
uenv=()
[[ -z ${MARS_AMR_UENV:-} ]] || uenv=(--uenv="$MARS_AMR_UENV")
srun "${uenv[@]}" --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=0 --cpu-bind=cores --exclusive --export=ALL \
    --kill-on-bad-exit=1 --uenv-passthrough=use --view=default \
    bash "$here/prepare.sh" "$backend" "$work" "${@:3}" > "$work/build.log" 2>&1
status=0
ctest --test-dir "$work/checks" --output-on-failure --no-tests=error -V \
    --output-junit "$work/results.xml" > "$work/ctest.log" 2>&1 || status=$?
printf '%s\n' "$status" > "$work/ctest-exit.txt"
cat "$work/ctest.log"
exit "$status"
