#!/bin/bash
# Strong scaling study: fixed 32M mesh, varying GPU count
# Usage: bash strong_scaling.sh [BUILD_DIR]
#
# Submits one SLURM job per GPU count; each job logs to LOG_DIR.
# After all jobs finish, run: bash strong_scaling.sh --parse

set -euo pipefail

BINARY=${BUILD_DIR:-$(pwd)}/examples/distributed/unstructured/mars_cvfem_graph
MESH=/users/fabianw/work/asm_test_mesh/block_32M.exo
AFFINITY=~/affinity/bind_numa.sh
ACCOUNT=csstaff
GPUS_PER_NODE=4  # GH200 nodes on santis
ITERATIONS=100
LOG_DIR=$(pwd)/strong_scaling_logs
KERNEL=tensor

# ── Parse mode ───────────────────────────────────────────────────────────────
if [[ "${1:-}" == "--parse" ]]; then
    echo ""
    echo "Strong Scaling — MARS CVFEM Hex (Graph+Lump tensor, 32M mesh)"
    echo ""
    printf "%-6s %-10s %-15s %-11s %-14s %-13s %-14s\n" \
        "GPUs" "Time(ms)" "BW/gpu(GB/s)" "Weak Eff." "Owned imbal." "Matrix norm" "Speedup"
    printf "%-6s %-10s %-15s %-11s %-14s %-13s %-14s\n" \
        "----" "--------" "------------" "---------" "------------" "-----------" "-------"

    baseline_time=""
    for ngpu in 1 2 4 8 16 32 64 128 256 512; do
        logfile="$LOG_DIR/strong_${ngpu}gpu.log"
        [[ ! -f "$logfile" ]] && printf "%-6s  (no log)\n" "$ngpu" && continue

        avg_time=$(grep "Average time:" "$logfile"      | awk '{printf "%.2f", $3}')
        bw=$(grep "Bandwidth:"       "$logfile"         | tail -1 | awk '{printf "%.2f", $2}')
        norm=$(grep "matrix norm:"   "$logfile"         | awk -F'matrix norm: ' '{print $2}' | awk '{printf "%s", $1}' | tr -d ']')
        owned_imbal=$(grep "Owned nodes:" "$logfile"    | awk '{printf "%.2f", $(NF-1)}')

        [[ -z "$avg_time" ]] && printf "%-6s  (parse error)\n" "$ngpu" && continue

        if [[ -z "$baseline_time" ]]; then
            baseline_time=$avg_time
            speedup="1.00x"
            eff="100.00 %"
        else
            speedup=$(awk "BEGIN {printf \"%.2fx\", $baseline_time / $avg_time}")
            eff=$(awk    "BEGIN {printf \"%.2f %%\", ($baseline_time / $avg_time) / $ngpu * 100}")
        fi

        printf "%-6s %-10s %-15s %-11s %-14s %-13s %-14s\n" \
            "$ngpu" "$avg_time" "$bw" "$eff" "${owned_imbal}%" "$norm" "$speedup"
    done
    echo ""
    echo "Strong efficiency = (T_1 / T_p) / p * 100"
    echo "Speedup           = T_1 / T_p"
    exit 0
fi

# ── Submit mode ──────────────────────────────────────────────────────────────
mkdir -p "$LOG_DIR"

submit_job() {
    local ngpu=$1
    local nodes=$(( (ngpu + GPUS_PER_NODE - 1) / GPUS_PER_NODE ))
    local ntasks_per_node=$(( ngpu < GPUS_PER_NODE ? ngpu : GPUS_PER_NODE ))
    local time_limit="00:30:00"
    [[ $ngpu -ge 64  ]] && time_limit="01:00:00"
    [[ $ngpu -ge 256 ]] && time_limit="02:00:00"

    local logfile="$LOG_DIR/strong_${ngpu}gpu.log"

    sbatch --job-name="mars_strong_${ngpu}" \
           --nodes=$nodes \
           --ntasks-per-node=$ntasks_per_node \
           --time=$time_limit \
           --account=$ACCOUNT \
           --output="$logfile" \
           --wrap="srun -l $AFFINITY $BINARY \
               --mesh=$MESH \
               --kernel=$KERNEL \
               --iterations=$ITERATIONS"

    echo "Submitted: ${ngpu} GPU(s) → $nodes node(s), ${ntasks_per_node} tasks/node → $logfile"
}

echo "Submitting strong scaling jobs (32M mesh, $ITERATIONS iterations)..."
echo ""
for ngpu in 1 2 4 8 16 32 64 128 256 512; do
    submit_job $ngpu
done

echo ""
echo "All jobs submitted. When complete, run:"
echo "  bash $0 --parse"
