#!/bin/bash
#SBATCH --account=csstaff
#SBATCH --job-name=ho1e12shflmr
#SBATCH --nodes=512
#SBATCH --ntasks-per-node=4
#SBATCH --time=00:30:00
#SBATCH --chdir=/capstor/scratch/cscs/gandanie/git/mars/daint-gpu
#SBATCH --output=/capstor/scratch/cscs/gandanie/git/mars/ho1e12shflmr_%j.out
#SBATCH --uenv-passthrough=use
#
# Step 2 of the isolation ladder: shuffle no-overlap (proven green in step 1) PLUS
# FI_MR_CACHE_MONITOR=userfaultfd. Run ONLY after ho_trillion_512nodes_shfl.sbatch
# completes. If this one MPI-errors while step 1 was green, the MR monitor is the culprit.
unset MARS_HO_RECOMPUTE MARS_HO_OVERLAP

echo "=== 2048 ranks / 512 nodes, ~1.0 TRILLION DOF (cube2500), store-d_G PA SHUFFLE, NO overlap + MR=userfaultfd ==="
srun --account=csstaff --nodes=512 --ntasks-per-node=4 \
  --export=ALL,MARS_HO_SHFL=1,MARS_HO_STAGE_TIMING=1,MARS_GLOBAL_BUCKETSIZE=256,MARS_HO_SKIP_HOST_GATE=1,FI_MR_CACHE_MONITOR=userfaultfd \
  ./examples/distributed/unstructured/mars_ho_dist_apply_test --ncells=2500 --p=4
