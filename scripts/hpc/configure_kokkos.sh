#!/bin/bash

#configure_kokkos.sh

NVCC_WRAPPER=nvcc_wrapper
cmake ..  \
    -DKokkos_ARCH_PASCAL61=ON \
    -DCMAKE_INSTALL_PREFIX=$HOME/installations/kokkos \
    -DKokkos_CUDA_DIR=$CUDATOOLKIT_DIR \
    -DKokkos_ENABLE_CUDA=ON \
    -DKokkos_ENABLE_CUDA_CONSTEXPR=ON \
    -DKokkos_ENABLE_CUDA_LAMBDA=ON \
    -DKokkos_ENABLE_CUDA_LDG_INTRINS=ON \
    -DKokkos_ENABLE_CUDA_RELOCATABLE=ON \
    -DKokkos_ENABLE_DEBUG=ON \
    -DKokkos_CXX_STANDARD=17 \
    -DCMAKE_CUDA_COMPILERa=${NVCC_WRAPPER} \
    -DCMAKE_CXX_COMPILER=${NVCC_WRAPPER}
