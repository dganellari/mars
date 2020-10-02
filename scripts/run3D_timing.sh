#!/bin/bash

SP_PATH=./st_example
export MPICH_RDMA_ENABLED_CUDA=1; export CUDA_LAUNCH_BLOCKING=1

for nx in {5,10,20,30,40,50}
do
   srun -Cgpu --pty -n1 -N1 --reservation=eurohack time $SP_PATH $nx 0
done