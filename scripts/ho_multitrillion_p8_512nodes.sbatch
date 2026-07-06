#!/bin/bash
#SBATCH --account=csstaff
#SBATCH --job-name=ho6e12p8
#SBATCH --nodes=512
#SBATCH --ntasks-per-node=4
#SBATCH --time=00:30:00
#SBATCH --chdir=/capstor/scratch/cscs/gandanie/git/mars/pump
#SBATCH --output=ho6e12p8_%j.out
#SBATCH --uenv-passthrough=use
#
# ~6 TRILLION DOF at p=8 on 512 nodes / 2048 GPU, RECOMPUTE (matrix-free, no d_G).
# cube2271 @ p=8 -> (2271*8+1)^3 = 18169^3 = 6.00e12 DOF.
# WHY this fits AND dodges the cstone 32-bit bug (which capped p=4 at ~1T):
#   - elements/rank = 2271^3 / 2048 = 5.72M  (< the ~7.6M that built fine at the 512-node p=4
#     trillion; the 448-node p=4 crash was at 8.72M). At p=8 each element carries 8^3=512 owned
#     DOF, so far fewer elements per DOF -> we slide well under the cstone element-count limit.
#   - recompute ~22 B/DOF -> 2.93e9 DOF/GPU * 22 = ~64 GB/GPU (fits 80 GB usable; ~80% of the
#     ~3.6e9 DOF/GPU recompute ceiling at p=8). store-d_G would be ~157 GB = OOM, so RECOMPUTE is
#     required. (p=8 memory ceiling caps near ~7.4T; this 6T leaves headroom for partition imbalance.)
# Throughput is recompute-limited (~1.5 GDOF/s/GPU); this is the SCALE headline (~6T on <half of
# Alps), not a throughput number. Gate: device A.1 ~ 1e-17 at ~6 trillion DOF.
export MARS_HO_RECOMPUTE=1          # no d_G: the only path that fits ~6T on 512 nodes
export MARS_NODEHALO_V2=1
export MARS_GLOBAL_BUCKETSIZE=256   # cstone global-tree Allreduce < 2 GiB (the trillion fix)
export MARS_HO_SKIP_HOST_GATE=1     # device A.1 is the gate
export MARS_HO_RESOLVE_TIMING=1
export MARS_HO_HALO_TIMING=1

echo "=== 2048 ranks / 512 nodes, ~6.0 TRILLION DOF (cube2271, p=8), RECOMPUTE (no d_G) ==="
srun --account=csstaff --nodes=512 --ntasks-per-node=4 --export=ALL \
  ./examples/distributed/unstructured/mars_ho_dist_apply_test --ncells=2271 --p=8
