#!/bin/bash
#SBATCH --account=csstaff
#SBATCH --job-name=ho1e12p7
#SBATCH --nodes=448
#SBATCH --ntasks-per-node=4
#SBATCH --time=02:00:00
#SBATCH --chdir=/capstor/scratch/cscs/gandanie/git/mars/daint-gpu
#SBATCH --output=ho1e12p7_%j.out
#SBATCH --uenv-passthrough=use
#
# TRILLION DOF at p=7: (7*1429+1)^3 = 1.001e12 DOF with only 1429^3 = 2.9B elements --
# 5.4x FEWER elements than cube2500/p=4 (15.6B). The cstone domain build (mesh + octree +
# sync) is element-bound and ate the ~1h wall on the p=4 run; at p=7 it drops to ~12 min and
# fits comfortably. d_G is leaner too (~100 B/DOF at p=7 vs ~139 at p=4): 559M DOF/GPU -> ~65 GB
# on 96 GB. Same binary, STORE-d_G default. Gate: device A.1 ~ 1e-18 at a trillion DOF.
export MARS_GLOBAL_BUCKETSIZE=256   # cstone global-tree Allreduce < 2 GiB (the trillion fix)
export MARS_HO_SKIP_HOST_GATE=1     # device A.1 is the gate
export MARS_HO_RESOLVE_TIMING=1     # per-stage resolve breakdown
export MARS_HO_HALO_TIMING=1        # per-stage GPU halo-build breakdown

echo "=== 1792 ranks / 448 nodes, ~1.0 TRILLION DOF at p=7 (cube1429), 5x fewer elements ==="
srun --nodes=448 --ntasks-per-node=4 --export=ALL \
  ./examples/distributed/unstructured/mars_ho_dist_apply_test --ncells=1429 --p=7
