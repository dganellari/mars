# Distributed SIMPLE gates. Shared MARS CMake is not modified.
#   cmake -S tests/reference/openaccel/distributed_simple -B build-dsimple && cmake --build build-dsimple
#   ctest --test-dir build-dsimple
# The CUDA/Hypre gate is added to an existing MARS build by ../distributed_matrix/inject.cmake.
cmake_minimum_required(VERSION 3.22)
project(mars_distributed_simple_gates LANGUAGES CXX)
enable_testing()
include(gates.cmake)
