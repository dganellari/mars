cmake_minimum_required(VERSION 3.18)
project(mfem_bp3_driver LANGUAGES CXX)

find_package(MFEM 4.8 CONFIG REQUIRED)
find_package(MPI REQUIRED COMPONENTS CXX)
find_package(CUDAToolkit REQUIRED)
if(NOT MFEM_USE_MPI OR NOT MFEM_USE_CUDA)
    message(FATAL_ERROR "The BP3 GPU baseline requires MFEM with MPI and CUDA")
endif()

add_executable(mfem_bp3 mfem_bp3.cpp)
target_compile_features(mfem_bp3 PRIVATE cxx_std_17)
# MFEM's exported host target omits CUDA headers needed by its public headers.
target_link_libraries(mfem_bp3 PRIVATE mfem MPI::MPI_CXX CUDA::cudart CUDA::cusparse CUDA::cublas CUDA::curand)
