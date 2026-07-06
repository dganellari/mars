#!/bin/bash
#SBATCH --account=csstaff
#SBATCH --job-name=ho6e12p8
#SBATCH --nodes=640
#SBATCH --ntasks-per-node=4
#SBATCH --time=00:30:00
#SBATCH --chdir=/capstor/scratch/cscs/gandanie/git/mars/pump
#SBATCH --output=ho6e12p8_640_%j.out
#SBATCH --uenv-passthrough=use
#
# ~6 TRILLION DOF at p=8 on 640 nodes / 2560 GPU, RECOMPUTE (matrix-free, no d_G).
# cube2271 @ p=8 -> 18169^3 = 6.00e12 DOF (same target as the 512-node attempt that OOM'd).
# WHY 640 (not 512, not 1024): the 512-node crash was in the NUMBERING (buildDistributedGpu),
#   which peaks at ~27 B/DOF. At 512 nodes 6T = 2.93e9 DOF/rank -> ~80+ GB -> OOM. The numbering
#   fits below ~2.6e9 DOF/rank. At 640 nodes 6T = 6e12/2560 = 2.34e9 DOF/rank -> numbering ~64 GB,
#   recompute apply ~51 GB -> both clear 80 GB. 640 is the minimum "a bit more than 512" that lands
#   the FULL 6T; on exactly 512 the safe max is ~4.7-5T (cube ~2100-2150).
export MARS_HO_RECOMPUTE=1          # no d_G: the only path that fits ~6T
export MARS_NODEHALO_V2=1
export MARS_GLOBAL_BUCKETSIZE=256   # cstone global-tree Allreduce < 2 GiB (the trillion fix)
export MARS_HO_SKIP_HOST_GATE=1     # device A.1 is the gate
export MARS_HO_RESOLVE_TIMING=1
export MARS_HO_HALO_TIMING=1

echo "=== 2560 ranks / 640 nodes, ~6.0 TRILLION DOF (cube2271, p=8), RECOMPUTE (no d_G) ==="
srun --account=csstaff --nodes=640 --ntasks-per-node=4 --export=ALL \
  ./examples/distributed/unstructured/mars_ho_dist_apply_test --ncells=2271 --p=8
