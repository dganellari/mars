#!/bin/bash
# Quick profiling script for weak scaling analysis
# Run directly on allocated node

set -e

BUILD_DIR=${BUILD_DIR:-$SCRATCH/git/mars/build}
MESH_DIR=${MESH_DIR:-$SCRATCH/git/meshes}
AFFINITY=${AFFINITY:-~/affinity/bind_numa.sh}

echo "============================================"
echo "MARS Weak Scaling - Quick Profile"
echo "============================================"
echo ""

# Check if nsys is available
if ! command -v nsys &> /dev/null; then
    echo "ERROR: nsys not found. Load with: module load nsight-systems"
    exit 1
fi

# Create output directory
mkdir -p profile_results
cd profile_results

echo "Step 1: Profile 1 rank (baseline - 4.75ms expected)"
echo "--------------------------------------------"
nsys profile -o mars_1rank \
    --trace=cuda,mpi \
    --stats=true \
    --force-overwrite=true \
    srun -l --nodes=1 --ntasks-per-node=1 $AFFINITY \
        $BUILD_DIR/examples/distributed/unstructured/mars_cvfem_graph \
        --mesh=$MESH_DIR/block_2M.exo \
        --kernel=tensor \
        --iterations=10

echo ""
echo "Step 2: Profile 8 ranks (weak scaling - 5.70ms expected)"
echo "--------------------------------------------"
nsys profile -o mars_8rank \
    --trace=cuda,mpi \
    --stats=true \
    --force-overwrite=true \
    srun -l --nodes=1 --ntasks-per-node=8 $AFFINITY \
        $BUILD_DIR/examples/distributed/unstructured/mars_cvfem_graph \
        --mesh=$MESH_DIR/block_16M.exo \
        --kernel=tensor \
        --iterations=10

echo ""
echo "============================================"
echo "Analysis: Kernel Time Comparison"
echo "============================================"
nsys stats --report cuda_gpu_kern_sum mars_1rank.nsys-rep 2>/dev/null | grep -A 5 "cvfem_hex_assembly"
echo ""
nsys stats --report cuda_gpu_kern_sum mars_8rank.nsys-rep 2>/dev/null | grep -A 5 "cvfem_hex_assembly"

echo ""
echo "============================================"
echo "Analysis: Memory Operations"
echo "============================================"
nsys stats --report cuda_gpu_mem_time_sum mars_1rank.nsys-rep 2>/dev/null | head -10
echo ""
nsys stats --report cuda_gpu_mem_time_sum mars_8rank.nsys-rep 2>/dev/null | head -10

echo ""
echo "============================================"
echo "Files generated:"
echo "  - mars_1rank.nsys-rep"
echo "  - mars_8rank.nsys-rep"
echo ""
echo "Detailed analysis commands:"
echo "  nsys stats --report cuda_gpu_kern_sum mars_1rank.nsys-rep"
echo "  nsys stats --report gpumemsizesum mars_8rank.nsys-rep"
echo "============================================"
