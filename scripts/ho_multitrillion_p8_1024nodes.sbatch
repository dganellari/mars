#!/bin/bash
#SBATCH --account=csstaff
#SBATCH --job-name=ho6e12p8
#SBATCH --nodes=1024
#SBATCH --ntasks-per-node=4
#SBATCH --time=00:30:00
#SBATCH --chdir=/capstor/scratch/cscs/gandanie/git/mars/pump
#SBATCH --output=ho6e12p8_1024_%j.out
#SBATCH --uenv-passthrough=use
#
# ~6 TRILLION DOF at p=8 on 1024 nodes / 4096 GPU, RECOMPUTE (matrix-free, no d_G).
# cube2271 @ p=8 -> (2271*8+1)^3 = 18169^3 = 6.00e12 DOF (same target as the 512-node attempt).
# WHY 1024 (the 512-node run OOM'd): the crash was in the NUMBERING, not the apply --
#   elements/rank at 512 nodes = 2271^3/2048 = 5.72M -> elemDof table (5.72M*9^3*4B = 16.7 GB)
#   + the DOF-key sort over ~2.9e9 DOF/rank blew past 80 GB during buildDistributedGpu.
#   At 1024 nodes elements/rank halve to 2271^3/4096 = 2.86M -> elemDof ~8.3 GB and the sort
#   over ~1.46e9 DOF/rank -> numbering peak ~52 GB < 80 GB. The recompute apply (~64 GB at 512)
#   also halves per rank, so both phases fit comfortably.
export MARS_HO_RECOMPUTE=1          # no d_G: the only path that fits ~6T
export MARS_NODEHALO_V2=1
export MARS_GLOBAL_BUCKETSIZE=256   # cstone global-tree Allreduce < 2 GiB (the trillion fix)
export MARS_HO_SKIP_HOST_GATE=1     # device A.1 is the gate
export MARS_HO_RESOLVE_TIMING=1
export MARS_HO_HALO_TIMING=1

echo "=== 4096 ranks / 1024 nodes, ~6.0 TRILLION DOF (cube2271, p=8), RECOMPUTE (no d_G) ==="
srun --account=csstaff --nodes=1024 --ntasks-per-node=4 --export=ALL \
  ./examples/distributed/unstructured/mars_ho_dist_apply_test --ncells=2271 --p=8
