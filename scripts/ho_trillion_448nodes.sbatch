#!/bin/bash
#SBATCH --account=csstaff
#SBATCH --job-name=ho1e12fit
#SBATCH --nodes=448
#SBATCH --ntasks-per-node=4
#SBATCH --time=03:00:00
#SBATCH --chdir=/capstor/scratch/cscs/gandanie/git/mars/daint-gpu
#SBATCH --output=ho1e12fit_%j.out
#SBATCH --uenv-passthrough=use
#
# FULL trillion matvec via the RECOMPUTE (matrix-free, no d_G) path. The store-d_G path
# DEADLOCKED at 448 nodes: ~78 GB/rank is at the 96 GB ceiling, so imbalanced ranks OOM'd on
# the 63 GB d_G, returned, and the rest hung forever on the alloc MPI_Barrier. Recompute drops
# d_G -> ~15 GB/rank -> fits with huge margin, no OOM, no deadlock. (It also fits the trillion
# back on 400 nodes / 1600 GPU now, since there is no d_G ceiling.) cube2500 ~= 1.005e12 DOF.
export MARS_NODEHALO_V2=1           # v2 per-node halo -- the standard HO run path (was MISSING)
export MARS_HO_RECOMPUTE=1          # no d_G: recompute metric inline -> fits, no OOM/deadlock
export MARS_GLOBAL_BUCKETSIZE=256   # cstone global-tree Allreduce < 2 GiB (the trillion fix)
export MARS_HO_SKIP_HOST_GATE=1     # device A.1 is the gate
export MARS_HO_RESOLVE_TIMING=1     # per-stage resolve breakdown
export MARS_HO_HALO_TIMING=1        # per-stage GPU halo-build breakdown

echo "=== 1792 ranks / 448 nodes, ~1.0 TRILLION DOF (cube2500), RECOMPUTE (no d_G) ==="
srun --account=csstaff --nodes=448 --ntasks-per-node=4 --export=ALL \
  ./examples/distributed/unstructured/mars_ho_dist_apply_test --ncells=2500 --p=4
