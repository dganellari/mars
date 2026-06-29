#!/bin/bash
# Drive NodeHaloTopology v2 through Gates 1–7 on Santis.
#
# Usage:
#   bash scripts/run_v2_gates.sh build           # rebuild
#   bash scripts/run_v2_gates.sh gate1           # cube16/64 validate (4 ranks)
#   bash scripts/run_v2_gates.sh gate5           # cube16/64/128 numeric match
#   bash scripts/run_v2_gates.sh gate6           # tet mesh
#   bash scripts/run_v2_gates.sh gate7           # cube256/512 scale test
#   bash scripts/run_v2_gates.sh attempt-1B      # cube1024 -> 1.07B elements
#   bash scripts/run_v2_gates.sh check-prior     # search for prior 1B-element logs
#
# Each "gateN" target is the smallest run that gates the next step.
# Don't skip ahead — Gate 1 failure means v2 CSRs are wrong, fix that first.

set -euo pipefail

# ─── Config ────────────────────────────────────────────────────────────────
ACCOUNT=${ACCOUNT:-csstaff}
GPUS_PER_NODE=4
BUILD_DIR=${BUILD_DIR:-$SCRATCH/git/mars/build}
MESH_DIR=${MESH_DIR:-$SCRATCH/git/meshes}
AFFINITY=${AFFINITY:-~/affinity/bind_numa.sh}
LOG_DIR=${LOG_DIR:-$(pwd)/v2_gate_logs}
ITERS=${ITERS:-50}
KERNEL=${KERNEL:-tensor}

mkdir -p "$LOG_DIR"

CVFEM_ASM=$BUILD_DIR/examples/distributed/unstructured/mars_cvfem_graph
CVFEM_SOLVE=$BUILD_DIR/examples/distributed/unstructured/mars_amr_cvfem_graph
POISSON_TET=$BUILD_DIR/examples/distributed/unstructured/mars_ex1_poisson

# Note on binary choice:
#   mars_cvfem_graph        - assembly-only benchmark, prints "matrix norm"
#   mars_amr_cvfem_graph    - full Poisson solve, prints "||u||" (use for Gate 5 numeric match)
#   mars_ex1_poisson        - tet Poisson solve, prints "L2 norm"

# ─── Helpers ───────────────────────────────────────────────────────────────
submit() {
    local name=$1 nodes=$2 ntasks_per_node=$3 timelimit=$4 logfile=$5
    shift 5
    local cmd="$*"
    sbatch --job-name="$name" \
           --nodes=$nodes \
           --ntasks-per-node=$ntasks_per_node \
           --time=$timelimit \
           --account=$ACCOUNT \
           --output="$logfile" \
           --wrap="$cmd" \
        | tee -a "$LOG_DIR/submission_record.log"
}

require() {
    [[ -x "$1" ]] || { echo "Missing binary: $1"; exit 1; }
    [[ -e "$2" ]] || { echo "Missing mesh:  $2"; exit 1; }
}

# ─── Targets ───────────────────────────────────────────────────────────────

case "${1:-}" in

build)
    set -x
    cd "$BUILD_DIR"
    cmake .. \
        -DMARS_ENABLE_KOKKOS=OFF \
        -DMARS_ENABLE_CUDA=ON \
        -DMARS_ENABLE_TESTS=ON \
        -DMARS_ENABLE_UNSTRUCTURED=ON \
        -DMARS_ENABLE_FEM_EXAMPLES=ON \
        -DCMAKE_CUDA_ARCHITECTURES=90
    make -j mars_cvfem_graph mars_ex1_poisson
    ;;

gate1)
    # Run host + v2 in same process, abort on per-peer node-id mismatch.
    # cube16 first (~4k elem), then cube64 (~262k elem).
    # Use the assembly binary — it triggers NodeHaloTopology construction
    # (which is what we're validating) without needing solver convergence.
    require "$CVFEM_ASM" "$MESH_DIR"
    for sz in 16 64; do
        mesh="$MESH_DIR/cube${sz}.exo"
        [[ -e "$mesh" ]] || mesh="$MESH_DIR/block_${sz}.exo"  # naming fallback
        log="$LOG_DIR/gate1_cube${sz}_4ranks.log"
        submit "v2_g1_c${sz}" 1 4 "00:10:00" "$log" \
            "MARS_NODEHALO_VALIDATE=1 srun -l $AFFINITY $CVFEM_ASM \
                --mesh=$mesh --kernel=$KERNEL --iterations=1"
    done
    echo ""
    echo "After both jobs finish:"
    echo "  grep 'Gate 1 pass\\|MISMATCH\\|abort' $LOG_DIR/gate1_*.log"
    ;;

gate5)
    # Numeric match: solve cube16/64/128 with both modes, diff ||u||.
    # Uses mars_amr_cvfem_graph because it runs the full CG solve and prints ||u||.
    # If your binary is mars_cvfem_graph (assembly only), this catches assembly
    # determinism but not full solver halo correctness.
    bin="$CVFEM_SOLVE"
    [[ -x "$bin" ]] || bin="$CVFEM_ASM"
    require "$bin" "$MESH_DIR"
    for sz in 16 64 128; do
        mesh="$MESH_DIR/cube${sz}.exo"
        [[ -e "$mesh" ]] || mesh="$MESH_DIR/block_${sz}.exo"
        for tag in host v2; do
            envvar=""
            [[ "$tag" == "v2" ]] && envvar="MARS_NODEHALO_V2=1"
            log="$LOG_DIR/gate5_cube${sz}_${tag}.log"
            submit "v2_g5_c${sz}_${tag}" 1 4 "00:30:00" "$log" \
                "$envvar srun -l $AFFINITY $bin \
                    --mesh=$mesh --kernel=$KERNEL --iterations=$ITERS"
        done
    done
    echo ""
    echo "After all jobs finish, compare with:"
    echo "  bash $0 parse-gate5"
    ;;

parse-gate5)
    # Looks for either '||u||=' or 'matrix norm:'  whichever the binary emits.
    extract() {
        local f=$1
        local v=$(grep -ho "||u||=[^ ]*"        "$f" 2>/dev/null | tail -1 | sed 's/||u||=//')
        [[ -z "$v" ]] && v=$(grep -ho "matrix norm: [^,]*"  "$f" 2>/dev/null | tail -1 | sed 's/matrix norm: //')
        echo "$v"
    }
    printf "%-10s %-6s %-26s %s\n" "size" "mode" "norm" "match"
    printf "%-10s %-6s %-26s %s\n" "----" "----" "----" "-----"
    for sz in 16 64 128; do
        host_norm=$(extract "$LOG_DIR/gate5_cube${sz}_host.log")
        v2_norm=$(  extract "$LOG_DIR/gate5_cube${sz}_v2.log")
        match="?"
        if [[ -n "$host_norm" && -n "$v2_norm" ]]; then
            [[ "$host_norm" == "$v2_norm" ]] && match="EXACT" || match="DIFF: $host_norm vs $v2_norm"
        fi
        printf "%-10s %-6s %-26s\n"     "cube${sz}" "host" "$host_norm"
        printf "%-10s %-6s %-26s %s\n"  "cube${sz}" "v2"   "$v2_norm"   "$match"
    done
    ;;

gate6)
    # Tet mesh — exercises NPC=4 path
    require "$POISSON_TET" "$MESH_DIR"
    tet_mesh="$MESH_DIR/tet_unit_cube.exo"
    [[ -e "$tet_mesh" ]] || { echo "No tet mesh at $tet_mesh"; exit 1; }
    for tag in host v2; do
        envvar=""
        [[ "$tag" == "v2" ]] && envvar="MARS_NODEHALO_V2=1"
        log="$LOG_DIR/gate6_tet_${tag}.log"
        submit "v2_g6_${tag}" 1 4 "00:30:00" "$log" \
            "$envvar srun -l $AFFINITY $POISSON_TET --mesh=$tet_mesh --iterations=$ITERS"
    done
    ;;

gate7)
    # Scale test: cube256/512, multiple rank counts. Watch NodeHaloTopo time.
    # Use assembly binary for the throughput measurement (no solver convergence
    # variability). Halo construction is the gate, not solver iteration count.
    require "$CVFEM_ASM" "$MESH_DIR"
    for sz in 256 512; do
        mesh="$MESH_DIR/cube${sz}.exo"
        [[ -e "$mesh" ]] || mesh="$MESH_DIR/block_${sz}.exo"
        for ngpu in 4 16 64; do
            nodes=$(( (ngpu + GPUS_PER_NODE - 1) / GPUS_PER_NODE ))
            ntpn=$(( ngpu < GPUS_PER_NODE ? ngpu : GPUS_PER_NODE ))
            tlim="01:00:00"
            [[ $ngpu -ge 64 ]] && tlim="02:00:00"
            log="$LOG_DIR/gate7_cube${sz}_${ngpu}gpu.log"
            submit "v2_g7_c${sz}_n${ngpu}" $nodes $ntpn $tlim "$log" \
                "MARS_NODEHALO_V2=1 srun -l $AFFINITY $CVFEM_ASM \
                    --mesh=$mesh --kernel=$KERNEL --iterations=$ITERS"
        done
    done
    ;;

parse-gate7)
    # Pull NodeHaloTopo construction time from each log
    printf "%-12s %-6s %-22s %-22s\n" "size" "GPUs" "NodeHaloTopo (ms)" "CG iter avg (ms)"
    printf "%-12s %-6s %-22s %-22s\n" "----" "----" "-----------------" "----------------"
    for sz in 256 512; do
        for ngpu in 4 16 64; do
            log="$LOG_DIR/gate7_cube${sz}_${ngpu}gpu.log"
            [[ -f "$log" ]] || continue
            topo=$(grep -i "NodeHaloTopo" "$log" | grep -oP '\d+\.?\d*\s*ms' | tail -1 | awk '{print $1}')
            cg=$(  grep -i "Average time" "$log" | tail -1 | awk '{print $3}')
            printf "%-12s %-6s %-22s %-22s\n" "cube${sz}" "$ngpu" "${topo:-?}" "${cg:-?}"
        done
    done
    ;;

attempt-1B)
    # cube1024^3 = 1.07B elements. Goal: complete the run, not optimize it.
    # 1024 ranks × 1.05M elem/rank ≈ 256 nodes × 4 GPUs.
    require "$CVFEM_ASM" "$MESH_DIR"
    mesh="$MESH_DIR/cube1024.exo"
    [[ -e "$mesh" ]] || mesh="$MESH_DIR/cube1024_mesh"  # binary format dir
    [[ -e "$mesh" ]] || {
        echo "No cube1024 mesh found. Generate first:"
        echo "  python scripts/generate_hex_cube.py --nx 1024 --ny 1024 --nz 1024 --output $MESH_DIR/cube1024_mesh"
        exit 1
    }
    log="$LOG_DIR/attempt1B_cube1024_1024gpu.log"
    submit "v2_1B" 256 4 "06:00:00" "$log" \
        "MARS_NODEHALO_V2=1 srun -l $AFFINITY $CVFEM_ASM \
            --mesh=$mesh --kernel=$KERNEL --iterations=20"
    echo ""
    echo "10^9 attempt submitted. Check $log when it lands."
    echo "Expect: NodeHaloTopo construction sub-second (was ~70s host)."
    echo "Expect: per-assembly-iter ~few ms (tensor kernel, ~1M elem/rank)."
    echo "If it OOMs or hangs: most likely NodeHaloTopo, then halo exchange."
    echo ""
    echo "Once assembly works, also try the AMR+CVFEM solve:"
    echo "  $CVFEM_SOLVE --mesh=\$mesh --iterations=$ITERS"
    ;;

check-prior)
    # Look for prior 10^9 runs you mentioned. Adjust paths as needed.
    echo "Searching standard scratch locations for prior 1B-element logs..."
    for root in \
        $SCRATCH \
        /capstor/scratch/cscs/$USER \
        /scratch/$USER \
        $HOME/git/mars \
    ; do
        [[ -d "$root" ]] || continue
        echo "--- under $root ---"
        find "$root" -maxdepth 6 \( -name "*1024*" -o -name "*1B*" -o -name "*1e9*" \) 2>/dev/null | head -20
        echo "--- logs mentioning 'NodeHaloTopo' under $root ---"
        grep -rln "NodeHaloTopo" "$root" 2>/dev/null | head -10 || true
    done
    ;;

*)
    cat <<EOF
Usage: $0 {build|gate1|gate5|parse-gate5|gate6|gate7|parse-gate7|attempt-1B|check-prior}

Suggested order:
  1. check-prior        # see if you already ran 10^9 somewhere
  2. build              # rebuild against current source
  3. gate1              # CSR construction matches host
  4. gate5              # full Poisson ||u|| matches
  5. gate6              # tet mesh
  6. gate7              # scale test cube256/512, watch NodeHaloTopo time
  7. attempt-1B         # cube1024 on 1024 GPUs

Don't skip ahead. Gate 1 failure means v2 CSRs are wrong; gate 5 failure
catches subtle issues that gate 1 misses (e.g. epoch/tag handling).
EOF
    ;;
esac
