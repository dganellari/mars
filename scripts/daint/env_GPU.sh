#!/bin/bash
# File: env_GPU.sh

# load env
source $APPS/UES/anfink/gpu/environment

module load daint-gpu
module load nvidia-nsight-systems

export MPICH_RDMA_ENABLED_CUDA=1; export CUDA_LAUNCH_BLOCKING=1
