#!/bin/bash
#SBATCH --account=csstaff
#SBATCH --job-name=ho1e12shfl
#SBATCH --nodes=512
#SBATCH --ntasks-per-node=4
#SBATCH --time=00:30:00
#SBATCH --chdir=/capstor/scratch/cscs/gandanie/git/mars/daint-gpu
#SBATCH --output=/capstor/scratch/cscs/gandanie/git/mars/ho1e12shfl_%j.out
#SBATCH --uenv-passthrough=use
#
# Trillion-DOF store-d_G PA SHUFFLE on 2048 GPU, NO overlap (comm exposed).
# The only 2048-GPU trillion that completed (10.24 TDOF/s) was no-overlap WITHOUT shuffle. The
# shuffle + FI_MR_CACHE_MONITOR=userfaultfd attempt died with an MPI error (log lost with the
# build dir), so this rerun isolates ONE variable: shuffle alone on the proven no-overlap path,
# NO MR monitor. If green, run ho_trillion_512nodes_shfl_mr.sbatch next to isolate the monitor.
# Output path is absolute so the log survives independent of the build dir.
unset MARS_HO_RECOMPUTE MARS_HO_OVERLAP FI_MR_CACHE_MONITOR

echo "=== 2048 ranks / 512 nodes, ~1.0 TRILLION DOF (cube2500), store-d_G PA SHUFFLE, NO overlap, no MR monitor ==="
srun --account=csstaff --nodes=512 --ntasks-per-node=4 \
  --export=ALL,MARS_HO_SHFL=1,MARS_HO_STAGE_TIMING=1,MARS_GLOBAL_BUCKETSIZE=256,MARS_HO_SKIP_HOST_GATE=1 \
  ./examples/distributed/unstructured/mars_ho_dist_apply_test --ncells=2500 --p=4
