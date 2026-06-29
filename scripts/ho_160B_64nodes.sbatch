#!/bin/bash
#SBATCH --account=csstaff
#SBATCH --job-name=ho160B
#SBATCH --nodes=64
#SBATCH --ntasks-per-node=4
#SBATCH --time=01:00:00
#SBATCH --chdir=$SCRATCH/git/mars/daint-gpu
#SBATCH --output=ho160B_%j.out
#SBATCH --uenv-passthrough=use
#
# HO matrix-free apply, p=4, 256 ranks (64 nodes x 4 GPU), ~625M DOF/GPU.
# Global: cube1357 -> 5429^3 = ~160 BILLION DOF.
# Expect: device A.1 PASS (~1e-19), apply ~6 GDOF/s/GPU, comm ~9% (peers maxed at 26).

export MARS_NODEHALO_V2=1
srun --export=ALL ./examples/distributed/unstructured/mars_ho_dist_apply_test --ncells=1357 --p=4
