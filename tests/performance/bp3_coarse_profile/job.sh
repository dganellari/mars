#!/usr/bin/env bash
#SBATCH --account=csstaff
#SBATCH --time=01:00:00
#SBATCH --output=bp3_coarse_%j.log
set -euo pipefail
kit=${BP3_PROFILE_KIT:-/capstor/scratch/cscs/gandanie/bp3-coarse-profile-kit}
[[ -f $kit/rank.sh ]]
mode=${1:-pilot}
case "$mode" in
    pilot) sizes=(1 8); required_nodes=2 ;;
    scale) sizes=(64); required_nodes=16 ;;
    *) echo 'Use pilot (1/8 ranks) or scale (64 ranks)' >&2; exit 2 ;;
esac
[[ ${SLURM_JOB_NUM_NODES:?Use sbatch} -ge $required_nodes ]]
[[ -z ${LD_PRELOAD:-} ]] || { echo 'Unset LD_PRELOAD; this experiment needs a clean stock control.' >&2; exit 2; }
build=${MARS_BUILD:-/capstor/scratch/cscs/gandanie/git/mars/mlir}
source_dir=$(dirname -- "$build")
binary="$build/examples/distributed/unstructured/mars_marsir_bp"
generated="$source_dir/marsir-mlir/generated"
results=${BP3_PROFILE_RESULTS:-/capstor/scratch/cscs/gandanie/bp3-coarse-profile/$SLURM_JOB_ID}
mkdir -p "$results"
exec > >(tee "$results/job.log") 2>&1
trap 'status=$?; printf "%s\n" "$status" > "$results/job.exit"' EXIT
python3 --version
export MPICH_GPU_SUPPORT_ENABLED=1
exclude=nid005323,nid005617,nid005669,nid005679,nid005742,nid005820,nid005854
git -C "$source_dir" merge-base --is-ancestor 339e6552 HEAD
git -C "$source_dir" rev-parse HEAD > "$results/mars.revision"
git -C "$source_dir" status --short > "$results/mars.status"
git -C "$source_dir" diff --quiet HEAD -- examples/distributed/unstructured/mars_marsir_bp.cu \
    backend/distributed/unstructured/solvers/mars_cellwise_amg.hpp \
    backend/distributed/unstructured/solvers/mars_cellwise_pmg.hpp \
    backend/distributed/unstructured/solvers/mars_cellwise_fcg.hpp
# Rebuild this target once; do not run another build in mlir during this job.
cmake --build "$build" --target mars_marsir_bp --parallel 8
python3 "$kit/build_hook.py" "$binary" "$results/build" > "$results/build.log"
hook="$results/build/libbp3_ranges.so"
sha256sum "$binary" "$hook" > "$results/binary.sha256"
find "$generated" -maxdepth 1 -type f -name '*.ptx' -print0 | sort -z | xargs -0 sha256sum > "$results/ptx.sha256"
nsys stats --help-reports > "$results/nsys-reports.txt"
for report in nvtx_pushpop_trace nvtx_gpu_proj_sum nvtx_kern_sum mpi_event_trace mpi_event_sum; do
    grep -q "$report" "$results/nsys-reports.txt"
done
printenv | grep -E '^(MPICH_|CUDA_VISIBLE_DEVICES|SLURM_JOB_ID|SLURM_JOB_NODELIST|UENV_)' > "$results/environment.txt"

run_case() {
    local ranks=$1 p=$2 ne=$3 kind=$4 label=$5 reference=${6:-}
    local nodes=$(( (ranks + 3) / 4 )) per_node=4
    if ((ranks < 4)); then per_node=$ranks; fi
    local output="$results/$label"
    mkdir "$output"
    srun --nodes="$nodes" --ntasks="$ranks" --ntasks-per-node="$per_node" \
        --kill-on-bad-exit=1 --exclude="$exclude" --export=ALL \
        bash "$kit/rank.sh" "$kind" "$output" "$hook" "$binary" \
        --ptx-dir "$generated" --spec bp3 --p "$p" --ne "$ne" --deform 0.05 \
        --mg --amg --amg-cycles 3 --pre 2 --post 2 --cheb-interval mars --rhs discrete \
        --tol 1e-13 --true-tol 1e-10 --maxit 2000 --solves 3 --reps 1 --amg-print \
        < /dev/null > "$output/solver.log" 2>&1
    local check=(python3 "$kit/analyze.py" "$output/solver.log" --ranks "$ranks" --output "$output/summary.json")
    if [[ -n $reference ]]; then check+=(--reference "$results/$reference/solver.log"); fi
    if [[ $kind != stock ]]; then check+=(--audits "$output"); fi
    if [[ $kind == trace ]]; then check+=(--trace); fi
    "${check[@]}" > "$output/summary.txt"
    echo "PASS $label"
}

# Prove the hooks and NVTX/MPI reports work before allocating large vectors.
run_case 1 3 8 stock smoke_stock
run_case 1 3 8 trace smoke_trace smoke_stock
for ranks in "${sizes[@]}"; do
    case "$ranks" in 1) ne3=64; ne7=32 ;; 8) ne3=128; ne7=64 ;; 64) ne3=256; ne7=128 ;; esac
    for p in 3 7; do
        if [[ $p == 3 ]]; then ne=$ne3; else ne=$ne7; fi
        for round in 1 2 3; do
            stock="p${p}_r${ranks}_stock${round}"
            probe="p${p}_r${ranks}_hook${round}"
            run_case "$ranks" "$p" "$ne" stock "$stock"
            run_case "$ranks" "$p" "$ne" hook "$probe" "$stock"
        done
        run_case "$ranks" "$p" "$ne" trace "p${p}_r${ranks}_trace" "p${p}_r${ranks}_stock3"
    done
done
sha256sum --check "$results/binary.sha256"
sha256sum --check "$results/ptx.sha256" > "$results/ptx-check.log"
python3 "$kit/summarize.py" "$results" > "$results/overview.txt"
echo "All cases passed. Results: $results"
