#!/bin/bash
mkdir build
cd build
cmake .. -DKokkos_DIR=/Users/dylan/Documents/Summer-Internship/Installations/kokkos/lib/cmake/Kokkos -DKokkosKernels_DIR=/Users/dylan/Documents/Summer-Internship/Installations/kokkos_kernels/lib/cmake/KokkosKernels -DMARS_ENABLE_CXXOPTS=ON
make -j4 && make install