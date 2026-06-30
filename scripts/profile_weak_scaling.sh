#!/bin/bash
#SBATCH --job-name=mars_profile
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=00:30:00
#SBATCH --account=csstaff
#SBATCH --partition=normal

# Profile weak scaling bottleneck in MARS CVFEM
# Compares 1 rank vs 8 ranks to identify where the 20% slowdown comes from

module load nvhpc cuda nsight-systems

BUILD_DIR=$SCRATCH/git/mars/build
MESH_DIR=$SCRATCH/git/meshes
AFFINITY=${AFFINITY:-~/affinity/bind_numa.sh}
OUTPUT_DIR=$BUILD_DIR/profile_results

mkdir -p $OUTPUT_DIR

echo "============================================"
echo "MARS CVFEM Weak Scaling Profile"
echo "============================================"
echo "Goal: Identify bottleneck causing 20% slowdown from 1→8 ranks"
echo ""

# Profile 1 rank baseline
echo "1. Profiling 1 rank (baseline)..."
nsys profile \
    -o $OUTPUT_DIR/mars_1rank \
    --trace=cuda,nvtx,mpi \
    --cuda-memory-usage=true \
    --stats=true \
    --force-overwrite=true \
    srun -l --nodes=1 --ntasks-per-node=1 $AFFINITY \
        $BUILD_DIR/examples/distributed/unstructured/mars_cvfem_graph \
        --mesh=$MESH_DIR/block_2M.exo \
        --kernel=tensor \
        --iterations=20

echo ""
echo "2. Profiling 8 ranks (weak scaling)..."
nsys profile \
    -o $OUTPUT_DIR/mars_8rank \
    --trace=cuda,nvtx,mpi \
    --cuda-memory-usage=true \
    --stats=true \
    --force-overwrite=true \
    srun -l --nodes=1 --ntasks-per-node=8 $AFFINITY \
        $BUILD_DIR/examples/distributed/unstructured/mars_cvfem_graph \
        --mesh=$MESH_DIR/block_16M.exo \
        --kernel=tensor \
        --iterations=20

echo ""
echo "============================================"
echo "Profile Analysis"
echo "============================================"
echo ""
echo "Generated profile files:"
echo "  - $OUTPUT_DIR/mars_1rank.nsys-rep"
echo "  - $OUTPUT_DIR/mars_8rank.nsys-rep"
echo ""
echo "To analyze:"
echo "  nsys stats --report cuda_gpu_kern_sum $OUTPUT_DIR/mars_1rank.nsys-rep"
echo "  nsys stats --report cuda_gpu_kern_sum $OUTPUT_DIR/mars_8rank.nsys-rep"
echo ""
echo "Key metrics to compare:"
echo "  1. Kernel execution time (should be similar for weak scaling)"
echo "  2. Memory bandwidth utilization (look for saturation)"
echo "  3. L2 cache hit rate (look for thrashing)"
echo "  4. SM occupancy (look for reduced per rank)"
echo "  5. Atomic operations (look for increased contention)"
echo ""

# Generate quick statistical comparison
echo "============================================"
echo "Quick Stats Comparison"
echo "============================================"
echo ""
echo "1 Rank Kernel Stats:"
nsys stats --report cuda_gpu_kern_sum --format table $OUTPUT_DIR/mars_1rank.nsys-rep 2>/dev/null | head -20

echo ""
echo "8 Rank Kernel Stats (per rank):"
nsys stats --report cuda_gpu_kern_sum --format table $OUTPUT_DIR/mars_8rank.nsys-rep 2>/dev/null | head -20

echo ""
echo "Memory bandwidth comparison:"
nsys stats --report cuda_gpu_mem_time_sum --format table $OUTPUT_DIR/mars_1rank.nsys-rep 2>/dev/null | head -15
nsys stats --report cuda_gpu_mem_time_sum --format table $OUTPUT_DIR/mars_8rank.nsys-rep 2>/dev/null | head -15

echo ""
echo "============================================"
echo "For detailed analysis, use Nsight Systems UI:"
echo "  nsys-ui $OUTPUT_DIR/mars_1rank.nsys-rep &"
echo "  nsys-ui $OUTPUT_DIR/mars_8rank.nsys-rep &"
echo "============================================"
