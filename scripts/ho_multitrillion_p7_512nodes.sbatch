#!/bin/bash
#SBATCH --account=csstaff
#SBATCH --job-name=ho5e12p7
#SBATCH --nodes=512
#SBATCH --ntasks-per-node=4
#SBATCH --time=00:30:00
#SBATCH --chdir=/capstor/scratch/cscs/gandanie/git/mars/pump
#SBATCH --output=ho5e12p7_%j.out
#SBATCH --uenv-passthrough=use
#
# ~5 TRILLION DOF at p=7 on 512 nodes / 2048 GPU, RECOMPUTE (matrix-free, no d_G).
# cube2447 @ p=7 -> (2447*7+1)^3 = 17130^3 = 5.03e12 DOF.
# WHY this fits AND dodges the cstone 32-bit bug (which capped p=4 at ~1T):
#   - elements/rank = 2447^3 / 2048 = 7.15M  (< the ~7.6M that built fine at the 512-node p=4
#     trillion; the 448-node p=4 crash was at 8.72M). HIGH p has far fewer elements per DOF, so
#     we slide UNDER the cstone element-count limit while scaling DOF up 5x.
#   - recompute footprint ~22.5 B/DOF -> 2.45e9 DOF/GPU * 22.5 = ~55 GB/GPU (fits the 80 GB usable).
#     store-d_G would be ~169 GB = OOM at p=7 this dense, so RECOMPUTE is required at this scale.
# Throughput here is recompute-limited (~1.5-2 GDOF/s/GPU); this is the SCALE headline (~5T on
# <half of Alps), not a throughput number. Gate: device A.1 ~ 1e-17 at ~5 trillion DOF.
export MARS_HO_RECOMPUTE=1          # no d_G: the only path that fits ~5T on 512 nodes
export MARS_NODEHALO_V2=1           # standard HO node-halo path (no-op env, kept for parity)
export MARS_GLOBAL_BUCKETSIZE=256   # cstone global-tree Allreduce < 2 GiB (the trillion fix)
export MARS_HO_SKIP_HOST_GATE=1     # device A.1 is the gate
export MARS_HO_RESOLVE_TIMING=1     # per-stage resolve breakdown
export MARS_HO_HALO_TIMING=1        # per-stage GPU halo-build breakdown

echo "=== 2048 ranks / 512 nodes, ~5.0 TRILLION DOF (cube2447, p=7), RECOMPUTE (no d_G) ==="
srun --account=csstaff --nodes=512 --ntasks-per-node=4 --export=ALL \
  ./examples/distributed/unstructured/mars_ho_dist_apply_test --ncells=2447 --p=7
