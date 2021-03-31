#!/bin/bash -l
#SBATCH --job-name=mars_build_jenkins
#SBATCH --time=01:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=36
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --constraint=mc
#SBATCH --error=/scratch/snx3000/gandanie/build/jenkins-mars-%j.err
#SBATCH --output=/scratch/snx3000/gandanie/build/jenkins-mars-%j.out

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

set -e
set -x

id
ls -lh
pwd -P
env
hostname

DATE=$(date +"%Y%m%d%H%M")
SCRATCH=/scratch/snx3000/anfink

for cpu_gpu in cpu ; do
  mkdir mars_${cpu_gpu}
  pushd mars_${cpu_gpu}
  source $SCRATCH/build/${cpu_gpu}/${cpu_gpu}_jenkins_Release/environment
  cmake -DCMAKE_BUILD_TYPE=release -DCMAKE_VERBOSE_MAKEFILE=true ..
  make -j16
  popd
done
