#!/bin/bash
#SBATCH --account=csstaff
#SBATCH --job-name=ho_nt
#SBATCH --nodes=400
#SBATCH --ntasks-per-node=4
#SBATCH --time=01:00:00
#SBATCH --chdir=/capstor/scratch/cscs/gandanie/git/mars/daint-gpu
#SBATCH --output=ho_nt_%j.out
#SBATCH --uenv-passthrough=use
#
# Near-trillion matvec that FITS the per-GPU d_G ceiling on the store-d_G (--PA) binary.
# cube2500 = 628M DOF/GPU needed ~87 GB and OOM'd; cube2450 -> ~591M DOF/GPU
# (~945B total on 1600 GPUs), ~82 GB, fits with margin. Gives A.1 + TDOF/s at
# near-trillion scale with NO rebuild (current validated binary). Pass ncells as arg 1.
# (The full 1.0T matvec needs the --MF/no-d_G path once it validates.)
NCELLS=${1:-2450}
export MARS_GLOBAL_BUCKETSIZE=256   # cstone global-tree Allreduce < 2 GiB (the trillion fix)
export MARS_HO_SKIP_HOST_GATE=1     # device A.1 is the gate
export MARS_HO_RESOLVE_TIMING=1     # per-stage resolve breakdown
export MARS_HO_HALO_TIMING=1        # per-stage GPU halo-build breakdown

echo "=== 1600 ranks / 400 nodes, cube${NCELLS} (fits the d_G ceiling) ==="
srun --nodes=400 --ntasks-per-node=4 --export=ALL \
  ./examples/distributed/unstructured/mars_ho_dist_apply_test --ncells=${NCELLS} --p=4
