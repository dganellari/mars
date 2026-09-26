# Owned-row distributed matrix gates. Shared MARS CMake is not modified.
#
# 1) Host/CPU-MPI gate, standalone:
#      cmake -S tests/reference/openaccel/distributed_matrix -B build-dmatrix && cmake --build build-dmatrix
#      ctest --test-dir build-dmatrix
# 2) CUDA/Hypre gate inside an existing MARS CUDA+Hypre build (same cache, same `mars` target):
#      cmake -S <mars> -B <build> -DCMAKE_PROJECT_mars_INCLUDE=<this dir>/inject.cmake
cmake_minimum_required(VERSION 3.22)
project(mars_distributed_matrix_gates LANGUAGES CXX)
enable_testing()
include(gates.cmake)
