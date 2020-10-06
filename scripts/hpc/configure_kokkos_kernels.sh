#!/bin/bash

#configure_kokkos_kernels.sh

NVCC_WRAPPER=nvcc_wrapper
KOKKOS_DIR==$HOME/installations/kokkos
cmake ..  \
    -DCMAKE_INSTALL_PREFIX=$HOME/installations/kokkos-kernels \
    -DKokkos_ENABLE_DEBUG=ON \
    -DKokkos_CXX_STANDARD=14 \
    -DCMAKE_CUDA_COMPILERa=${NVCC_WRAPPER} \
    -DCMAKE_CXX_COMPILER=${NVCC_WRAPPER}