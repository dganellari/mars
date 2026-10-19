cmake_minimum_required(VERSION 3.24)
project(mars_usage LANGUAGES C CXX)

# -DCMAKE_PREFIX_PATH=<MARS install prefix>
find_package(Mars 0 REQUIRED)

# The GPU language (CUDA or HIP) and architectures MARS was built for, enabled before any target is
# made: the targets that link Mars::mars need them. The definitions and nvcc flags the MARS headers
# need come with Mars::mars.
if(MARS_GPU_LANGUAGE)
    if(NOT DEFINED CMAKE_${MARS_GPU_LANGUAGE}_ARCHITECTURES)
        set(CMAKE_${MARS_GPU_LANGUAGE}_ARCHITECTURES "${MARS_GPU_ARCHITECTURES}")
    endif()
    set(CMAKE_${MARS_GPU_LANGUAGE}_STANDARD 20)
    enable_language(${MARS_GPU_LANGUAGE})
endif()

# host-only smoke: the installed headers compile and Mars::mars links
add_executable(mars_example main.cpp)
target_link_libraries(mars_example PUBLIC Mars::mars)

# a GPU program, the appendix of the step-by-step tutorial
if(MARS_GPU_LANGUAGE)
    add_executable(mars_gpu_example gpu_example.cu)
    set_source_files_properties(gpu_example.cu PROPERTIES LANGUAGE ${MARS_GPU_LANGUAGE})
    target_link_libraries(mars_gpu_example PRIVATE Mars::mars)
    # With CUDA the MARS libraries hold relocatable device code, so the program needs a device link
    # step; CMake adds it on its own only inside the MARS build.
    if(MARS_GPU_LANGUAGE STREQUAL "CUDA")
        set_target_properties(mars_gpu_example PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    endif()
endif()

enable_testing()
add_test(NAME Mars.mars_example COMMAND mars_example)
if(TARGET mars_gpu_example)
    add_test(NAME Mars.mars_gpu_example COMMAND mars_gpu_example)
endif()
