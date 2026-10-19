# Appendix: your own CMake project

The steps build inside the MARS tree. For your own code, install MARS once and build your project
against it with `find_package(Mars)`. The project `examples/usage_from_external_cmake_project` does
exactly that, and the release checks build and run it on NVIDIA and AMD GPUs.

## Install MARS

From the MARS build directory of [the quickstart](Quickstart.md):

```bash
cmake --install build --prefix $HOME/mars-install
```

The prefix gets the libraries, the headers in the layout of the source tree (so
`#include "backend/distributed/unstructured/domain.hpp"` works as in the steps), the bundled
cornerstone headers, and the CMake package files.

## The project

```cmake
cmake_minimum_required(VERSION 3.24)
project(my_mars_app LANGUAGES C CXX)

find_package(Mars 0 REQUIRED)

# CUDA or HIP, as MARS was built; the GPU architectures default to MARS's
set(CMAKE_${MARS_GPU_LANGUAGE}_ARCHITECTURES "${MARS_GPU_ARCHITECTURES}")
set(CMAKE_${MARS_GPU_LANGUAGE}_STANDARD 20)
enable_language(${MARS_GPU_LANGUAGE})

add_executable(my_app main.cu)
set_source_files_properties(main.cu PROPERTIES LANGUAGE ${MARS_GPU_LANGUAGE})
target_link_libraries(my_app PRIVATE Mars::mars)
if(MARS_GPU_LANGUAGE STREQUAL "CUDA")
    set_target_properties(my_app PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()
```

`find_package(Mars)` sets `MARS_GPU_LANGUAGE` (`CUDA` or `HIP`) and `MARS_GPU_ARCHITECTURES`, so
the same file builds against an NVIDIA or an AMD installation. Enable the GPU language before you
create a target. With CUDA, the MARS libraries contain relocatable device code, so the program needs
a device link step: `CUDA_SEPARABLE_COMPILATION ON` asks CMake for it. Linking `Mars::mars` brings
everything else: the include paths, the definitions the MARS headers need, the device-lambda flags
of nvcc, and the libraries (MPI, the GPU math libraries, Hypre and netCDF when MARS uses them).

`examples/usage_from_external_cmake_project/gpu_example.cu` is a complete `main.cu`: it generates
the unit cube on every rank, lets cornerstone distribute it, and sums the element volumes on the GPU,
as in [step 2](step-02-volume.md).

## Build and run

```bash
cmake -S my_mars_app -B my_mars_app/build -DCMAKE_PREFIX_PATH=$HOME/mars-install \
      -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++
cmake --build my_mars_app/build
mpirun -np 4 my_mars_app/build/my_app
```

For the example project this prints `MARS GPU example: volume = 1.000000000000000 (ok)`.

Use the compilers MARS was built with. With HIP, also the Hypre and netCDF that MARS found: their
paths are recorded in the installed package.
