#!/usr/bin/env bash
set -euo pipefail
mode=$1; destination=$2; hook=$3; binary=$4
shift 4
rank=${SLURM_PROCID:?Run through srun}
mkdir -p "$destination"
if [[ $mode == stock ]]; then
    exec "$binary" "$@"
fi
audit="$destination/rank${rank}.audit.json"
if [[ $mode == hook ]]; then
    exec env BP3_RANGE_AUDIT="$audit" LD_PRELOAD="$hook${LD_PRELOAD:+:$LD_PRELOAD}" "$binary" "$@"
fi
[[ $mode == trace ]]
report="$destination/rank${rank}"
nsys profile --trace=cuda,nvtx,mpi --mpi-impl=mpich --sample=none --cpuctxsw=none \
    --capture-range=cudaProfilerApi --capture-range-end=stop --force-overwrite=false \
    --output="$report" \
    bash "${BASH_SOURCE[0]}" hook "$destination" "$hook" "$binary" "$@" --profile
for kind in nvtx_pushpop_trace nvtx_gpu_proj_sum nvtx_kern_sum cuda_gpu_kern_sum cuda_api_sum mpi_event_sum mpi_event_trace; do
    if ! nsys stats --report "$kind" --format csv --output - "$report.nsys-rep" \
        > "$report.$kind.csv" 2> "$report.$kind.stderr"; then
        # One rank can avoid every traced MPI call. Multi-rank missing data is an error.
        [[ $kind == mpi_event_* && ${SLURM_STEP_NUM_TASKS:-${SLURM_NTASKS:-0}} == 1 ]] && \
            grep -Ei 'does not contain.*MPI|no MPI.*events' "$report.$kind.stderr" "$report.$kind.csv" || exit 1
    fi
done
