#!/bin/bash
#SBATCH --account=csstaff
#SBATCH --job-name=ho1e12
#SBATCH --nodes=400
#SBATCH --ntasks-per-node=4
#SBATCH --time=00:45:00
#SBATCH --chdir=/capstor/scratch/cscs/gandanie/git/mars/daint-gpu
#SBATCH --output=ho1e12_%j.out
#
# HO matrix-free apply, p=4, 1600 ranks (400 nodes x 4 GPU), ~625M DOF/GPU.
# Global: cube2500 -> 10001^3 = ~1.0 TRILLION DOF.
# Expect: device A.1 PASS (~1e-19), apply ~6 GDOF/s/GPU (~9.5 TDOF/s aggregate), comm ~9%.
#
# Insurance: a 1024-rank (640B DOF) run first, so a result lands even if the full
# trillion hits a surprise at 1600 ranks. Both steps run inside this one allocation.

export MARS_NODEHALO_V2=1

echo "=== step 1: 1024 ranks, ~640B DOF (cube2155) ==="
srun --nodes=256 --ntasks-per-node=4 --export=ALL \
  ./examples/distributed/unstructured/mars_ho_dist_apply_test --ncells=2155 --p=4

echo "=== step 2: 1600 ranks, ~1.0 TRILLION DOF (cube2500) ==="
srun --nodes=400 --ntasks-per-node=4 --export=ALL \
  ./examples/distributed/unstructured/mars_ho_dist_apply_test --ncells=2500 --p=4
