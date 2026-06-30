#!/bin/bash
#SBATCH --account=csstaff
#SBATCH --job-name=ho1e12
#SBATCH --nodes=400
#SBATCH --ntasks-per-node=4
#SBATCH --time=03:00:00
#SBATCH --output=ho1e12_%j.out
#SBATCH --uenv-passthrough=use
#
# HO matrix-free apply, p=4, 1600 ranks (400 nodes x 4 GPU), ~625M DOF/GPU.
# Global: cube2500 -> 10001^3 = ~1.0 TRILLION DOF.
# Expect: device A.1 PASS (~1e-19), apply ~6 GDOF/s/GPU (~9.5 TDOF/s aggregate), comm ~9%.
#
# Trillion-only: all 400 nodes go straight to cube2500 (no 640B insurance step, which
# previously ate the wall before the trillion could start). Run 640B as its own job if
# you want an insurance result.

export MARS_NODEHALO_V2=1
export MARS_GLOBAL_BUCKETSIZE=256   # keeps cstone global-tree Allreduce < 2 GiB (the trillion fix)
export MARS_HO_SKIP_HOST_GATE=1     # skip the slow host-halo A.1 cross-check; device A.1 is the gate
export MARS_HO_RESOLVE_TIMING=1     # per-stage resolve breakdown (cheap) -> shows P2P weak-scales at 1600 GPU
export MARS_HO_HALO_TIMING=1        # per-stage GPU halo-build breakdown (cheap) -> same, for the halo
# GPU DOF numbering, GPU P2P ownership resolve, AND GPU halo build are ALL the default now (no env).
# MARS_HO_HOST=1 would force the host resolve+halo -- we are NOT setting it, so the full GPU path runs.

# SLURM does not expand vars in #SBATCH directives, so cd to the build dir here.
cd "${BUILD_DIR:-$SCRATCH/git/mars/daint-gpu}"
echo "=== 1600 ranks / 400 nodes, ~1.0 TRILLION DOF (cube2500) ==="
srun --nodes=400 --ntasks-per-node=4 --export=ALL \
  ./examples/distributed/unstructured/mars_ho_dist_apply_test --ncells=2500 --p=4
