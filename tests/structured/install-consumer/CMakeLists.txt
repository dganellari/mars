cmake_minimum_required(VERSION 3.24)
project(mars_structured_consumer LANGUAGES CXX)

find_package(Mars REQUIRED)
find_package(MPI REQUIRED COMPONENTS CXX)

if(MARS_GPU_LANGUAGE)
    if(NOT DEFINED CMAKE_${MARS_GPU_LANGUAGE}_ARCHITECTURES)
        set(CMAKE_${MARS_GPU_LANGUAGE}_ARCHITECTURES "${MARS_GPU_ARCHITECTURES}")
    endif()
    set(CMAKE_${MARS_GPU_LANGUAGE}_STANDARD 20)
    set(CMAKE_${MARS_GPU_LANGUAGE}_STANDARD_REQUIRED ON)
    enable_language(${MARS_GPU_LANGUAGE})
    set_source_files_properties(main.cpp PROPERTIES LANGUAGE ${MARS_GPU_LANGUAGE})
endif()

add_executable(mars_structured_consumer main.cpp)
if(MARS_GPU_LANGUAGE STREQUAL "CUDA" AND Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE)
    set_target_properties(mars_structured_consumer PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()
target_link_libraries(mars_structured_consumer PRIVATE Mars::mars)
enable_testing()
if(MARS_GPU_LANGUAGE OR MARS_STRUCTURED_FORCE_MPI)
    if(NOT MPIEXEC_EXECUTABLE)
        message(FATAL_ERROR "MPIEXEC_EXECUTABLE is required for the structured GPU consumer, including one-rank cases.")
    endif()
    add_test(NAME structured_install
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
            ${MPIEXEC_PREFLAGS} $<TARGET_FILE:mars_structured_consumer> ${MPIEXEC_POSTFLAGS})
else()
    add_test(NAME structured_install COMMAND mars_structured_consumer)
endif()
set_tests_properties(structured_install PROPERTIES LABELS "structured;kokkos;structured-install" PROCESSORS 1 TIMEOUT 120)

set(MARS_STRUCTURED_TEST_RANKS "2;4" CACHE STRING "MPI rank counts for the installed structured consumer")
foreach(ranks IN LISTS MARS_STRUCTURED_TEST_RANKS)
    if(NOT ranks MATCHES "^[1-9][0-9]*$")
        message(FATAL_ERROR "MARS_STRUCTURED_TEST_RANKS must contain positive integers.")
    endif()
    if(NOT MPIEXEC_EXECUTABLE)
        message(FATAL_ERROR "MPIEXEC_EXECUTABLE is required for structured MPI consumer tests.")
    endif()
    add_test(NAME structured_install_mpi${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks}
            ${MPIEXEC_PREFLAGS} $<TARGET_FILE:mars_structured_consumer> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(structured_install_mpi${ranks} PROPERTIES
        LABELS "structured;kokkos;structured-install;mpi" PROCESSORS ${ranks} TIMEOUT 120)
endforeach()
