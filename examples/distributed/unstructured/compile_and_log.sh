#!/bin/bash
cd /Users/gandanie/scratch/santis/mars/examples/distributed/unstructured
/usr/local/cuda/bin/nvcc -std=c++17 -O3 \
  -I../../../backend/distributed/unstructured \
  -I../../../backend/distributed/unstructured/utils \
  -I../../../base \
  -I../../../core \
  -I../../../io \
  -I/Users/gandanie/scratch/view_mars_dev/include \
  -I/Users/gandanie/spack/opt/spack/darwin-sonoma-m2/apple-clang-15.0.0/openmpi-4.1.6-pknbrs55rcovztldhn35tyuea6pa2mqj/include \
  -L../../lib \
  -L/Users/gandanie/spack/opt/spack/darwin-sonoma-m2/apple-clang-15.0.0/openmpi-4.1.6-pknbrs55rcovztldhn35tyuea6pa2mqj/lib \
  -lmars \
  -lmpi \
  -lcublas -lcusparse \
  mars_ex_beam_tet_distributed.cu -o mars_ex_beam_tet_distributed > compile.log 2>&1
