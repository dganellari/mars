#!/bin/bash
#SBATCH --account=csstaff
#SBATCH --job-name=ho1e12PA
#SBATCH --nodes=512
#SBATCH --ntasks-per-node=4
#SBATCH --time=00:20:00
#SBATCH --chdir=/capstor/scratch/cscs/gandanie/git/mars/daint-gpu
#SBATCH --output=ho1e12PA_%j.out
#SBATCH --uenv-passthrough=use
#
# Trillion-DOF store-d_G PARTIAL ASSEMBLY on 512 nodes / 2048 GPU (491M DOF/GPU).
# store-d_G ~68 GB peak < 80 GB usable -> fits (recompute run's own ceiling line:
# store-d_G 138.6 B/DOF -> 5.8e8 DOF/GPU, and 491M < 580M).
# 448 (8.72M elem/rank) crashed in the cstone build (illegal address); 512 (7.63M/rank)
# clears it -- the recompute run already built fine here.
unset MARS_HO_RECOMPUTE              # CRITICAL: force store-d_G PA, kill any leaked shell export
export MARS_GLOBAL_BUCKETSIZE=256   # cstone global-tree Allreduce < 2 GiB (the trillion fix)
export MARS_HO_SKIP_HOST_GATE=1     # device A.1 is the gate at scale
export MARS_HO_RESOLVE_TIMING=1
export MARS_HO_HALO_TIMING=1

echo "=== 2048 ranks / 512 nodes, ~1.0 TRILLION DOF (cube2500), store-d_G PA ==="
srun --account=csstaff --nodes=512 --ntasks-per-node=4 --export=ALL \
  ./examples/distributed/unstructured/mars_ho_dist_apply_test --ncells=2500 --p=4
