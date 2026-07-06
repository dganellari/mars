#!/bin/bash
#SBATCH --account=csstaff
#SBATCH --job-name=ho_p8mt
#SBATCH --ntasks-per-node=4
#SBATCH --time=00:30:00
#SBATCH --chdir=/capstor/scratch/cscs/gandanie/git/mars/pump
#SBATCH --output=ho_p8mt_%j.out
#SBATCH --uenv-passthrough=use
#
# p=8 MULTI-TRILLION, RECOMPUTE (matrix-free, no d_G). Submit with --nodes + P8NCELLS:
#   sbatch --nodes=640 --export=ALL,P8NCELLS=2271 ho_multitrillion_p8.sbatch   # ~6.0T, 4.57M elem/rank (tight)
#   sbatch --nodes=768 --export=ALL,P8NCELLS=2271 ho_multitrillion_p8.sbatch   # ~6.0T, 3.81M elem/rank (comfortable)
#   sbatch --nodes=512 --export=ALL,P8NCELLS=2050 ho_multitrillion_p8.sbatch   # ~4.4T, 4.20M elem/rank (safe, no extra nodes)
# The 512/cube2271 (6T) attempt OOM'd in the NUMBERING (buildDistributedGpu), NOT the apply:
# elemDof table + DOF-key sort scale as elements/rank * (p+1)^3, and 5.72M elem/rank blew past
# 80 GB. The limiter is elements/rank -> keep it under ~4.5M (more nodes OR smaller cube).
export MARS_HO_RECOMPUTE=1          # no d_G: the only path that fits multi-trillion
export MARS_NODEHALO_V2=1
export MARS_GLOBAL_BUCKETSIZE=256   # cstone global-tree Allreduce < 2 GiB
export MARS_HO_SKIP_HOST_GATE=1     # device A.1 is the gate
export MARS_HO_RESOLVE_TIMING=1
export MARS_HO_HALO_TIMING=1

echo "=== ${SLURM_NNODES} nodes / $((SLURM_NNODES*4)) GPU, p=8 RECOMPUTE, cube${P8NCELLS} ==="
srun --account=csstaff --ntasks-per-node=4 --export=ALL \
  ./examples/distributed/unstructured/mars_ho_dist_apply_test --ncells=${P8NCELLS} --p=8
