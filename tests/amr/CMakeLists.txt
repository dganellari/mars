cmake_minimum_required(VERSION 3.24)
project(mars_tet_bisection_checks LANGUAGES CXX)
enable_testing()

get_filename_component(MARS_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
add_executable(tet_bisection_host tet_bisection_host.cpp)
target_compile_features(tet_bisection_host PRIVATE cxx_std_20)
target_include_directories(tet_bisection_host PRIVATE "${MARS_SOURCE_DIR}")
add_test(NAME tet_bisection_host COMMAND tet_bisection_host)
set_tests_properties(tet_bisection_host PROPERTIES LABELS "amr;tet-bisection;host" TIMEOUT 60)

option(MARS_AMR_MPI_CHECKS "Build CPU checks of the distributed LEPP engine" OFF)
if(MARS_AMR_MPI_CHECKS)
    find_package(MPI REQUIRED COMPONENTS CXX)
    find_path(MARS_CSTONE_INCLUDE_DIR cstone/primitives/mpi_wrappers.hpp REQUIRED)
    add_executable(tet_lepp_mpi tet_lepp_mpi.cpp)
    target_compile_features(tet_lepp_mpi PRIVATE cxx_std_20)
    target_include_directories(tet_lepp_mpi PRIVATE "${MARS_SOURCE_DIR}" "${MARS_CSTONE_INCLUDE_DIR}")
    target_link_libraries(tet_lepp_mpi PRIVATE MPI::MPI_CXX)
    add_executable(tet_domain_layout tet_domain_layout.cpp)
    target_compile_features(tet_domain_layout PRIVATE cxx_std_20)
    target_include_directories(tet_domain_layout PRIVATE "${MARS_SOURCE_DIR}" "${MARS_CSTONE_INCLUDE_DIR}")
    target_link_libraries(tet_domain_layout PRIVATE MPI::MPI_CXX)
    add_test(NAME tet_domain_layout_mpi4
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
            ${MPIEXEC_PREFLAGS} $<TARGET_FILE:tet_domain_layout> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(tet_domain_layout_mpi4 PROPERTIES
        LABELS "amr;tet-lepp;host;mpi;halo-layout" PROCESSORS 4 TIMEOUT 120)
    foreach(ranks IN ITEMS 1 2 4)
        foreach(mode IN ITEMS direct staged)
            set(mode_flag)
            if(mode STREQUAL "staged")
                set(mode_flag --staged)
            endif()
            add_test(NAME tet_lepp_${mode}_mpi${ranks}
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks}
                    ${MPIEXEC_PREFLAGS} $<TARGET_FILE:tet_lepp_mpi> ${MPIEXEC_POSTFLAGS} ${mode_flag})
            set_tests_properties(tet_lepp_${mode}_mpi${ranks} PROPERTIES
                LABELS "amr;tet-lepp;host;mpi;tet-conformity;tet-transfer;tet-boundary" PROCESSORS ${ranks} TIMEOUT 120)
        endforeach()
    endforeach()
endif()

option(MARS_AMR_LEGACY_CHECKS "Check the legacy Kokkos LEPP adapter against its terminal-edge witness" OFF)
if(MARS_AMR_LEGACY_CHECKS)
    find_package(Mars REQUIRED)
    if(NOT Kokkos_VERSION VERSION_EQUAL "4.4.1" OR Kokkos_ENABLE_CUDA OR Kokkos_ENABLE_HIP)
        message(FATAL_ERROR "The legacy CPU gate requires the pinned CPU Kokkos 4.4.1 package; check Kokkos_DIR.")
    endif()
    add_executable(tet_legacy_lepp tet_legacy_lepp.cpp)
    target_link_libraries(tet_legacy_lepp PRIVATE Mars::mars)
    target_compile_features(tet_legacy_lepp PRIVATE cxx_std_20)
    add_test(NAME tet_legacy_lepp COMMAND tet_legacy_lepp)
    set_tests_properties(tet_legacy_lepp PROPERTIES LABELS "amr;tet-lepp;legacy;host" TIMEOUT 60)
endif()

option(MARS_AMR_EXODUS_CHECKS "Check public Tet4 Exodus import, labels and AMR" OFF)
option(MARS_AMR_GPU_CHECKS "Build the installed MARS GPU/MPI conforming tet AMR gate" OFF)
if(MARS_AMR_GPU_CHECKS)
    find_package(Mars REQUIRED)
    find_package(MPI REQUIRED COMPONENTS CXX)
    if(NOT TARGET Mars::mars_unstructured OR NOT MARS_GPU_LANGUAGE MATCHES "^(CUDA|HIP)$")
        message(FATAL_ERROR "MARS_AMR_GPU_CHECKS needs an installed CUDA/HIP unstructured MARS library.")
    endif()
    if(MARS_AMR_EXODUS_CHECKS)
        get_target_property(unstructured_definitions Mars::mars_unstructured INTERFACE_COMPILE_DEFINITIONS)
        if(NOT "MARS_HAVE_NETCDF" IN_LIST unstructured_definitions)
            message(FATAL_ERROR "Exodus GPU checks need netCDF enabled in the installed MARS package.")
        endif()
    endif()
    # Imported HIP libraries can select the linker even for a C++ reader target.
    # Set the architecture before creating any target that links those libraries.
    if(NOT DEFINED CMAKE_${MARS_GPU_LANGUAGE}_ARCHITECTURES)
        set(CMAKE_${MARS_GPU_LANGUAGE}_ARCHITECTURES "${MARS_GPU_ARCHITECTURES}")
    endif()
    set(CMAKE_${MARS_GPU_LANGUAGE}_STANDARD 20)
    enable_language(${MARS_GPU_LANGUAGE})
endif()
set(MARS_AMR_TEST_RANKS "1;2;4" CACHE STRING "MPI rank counts for conforming tetrahedral AMR")
foreach(ranks IN LISTS MARS_AMR_TEST_RANKS)
    if(NOT ranks MATCHES "^[1-9][0-9]*$")
        message(FATAL_ERROR "MARS_AMR_TEST_RANKS must contain positive integers.")
    endif()
endforeach()
if(MARS_AMR_EXODUS_CHECKS)
    find_package(Mars REQUIRED)
    find_package(MPI REQUIRED COMPONENTS CXX)
    find_path(MARS_CSTONE_INCLUDE_DIR cstone/primitives/mpi_wrappers.hpp REQUIRED)
    find_path(NETCDF_INCLUDE_DIR netcdf.h HINTS ENV NETCDF_DIR ENV NETCDF_ROOT PATH_SUFFIXES include REQUIRED)
    find_library(NETCDF_LIBRARY netcdf HINTS ENV NETCDF_DIR ENV NETCDF_ROOT PATH_SUFFIXES lib lib64 REQUIRED)

    add_executable(write_tet_exodus write_tet_exodus.cpp)
    target_compile_features(write_tet_exodus PRIVATE cxx_std_20)
    target_include_directories(write_tet_exodus PRIVATE "${NETCDF_INCLUDE_DIR}")
    target_link_libraries(write_tet_exodus PRIVATE "${NETCDF_LIBRARY}")
    set(exodus_fixture "${CMAKE_CURRENT_BINARY_DIR}/tet-amr.exo")
    add_test(NAME tet_exodus_fixture COMMAND write_tet_exodus "${exodus_fixture}")
    set_tests_properties(tet_exodus_fixture PROPERTIES
        LABELS "amr;tet-exodus;host" FIXTURES_SETUP tet_exodus TIMEOUT 30)

    add_executable(tet_exodus_host tet_exodus_host.cpp tet_reference.cpp)
    target_compile_features(tet_exodus_host PRIVATE cxx_std_20)
    target_compile_definitions(tet_exodus_host PRIVATE MARS_HAVE_NETCDF)
    target_include_directories(tet_exodus_host PRIVATE
        "${MARS_SOURCE_DIR}" "${MARS_CSTONE_INCLUDE_DIR}")
    target_include_directories(tet_exodus_host SYSTEM PRIVATE "${NETCDF_INCLUDE_DIR}")
    if(TARGET Kokkos::kokkos)
        # The netCDF/MPI prefix may also contain an unrelated GPU Kokkos.
        target_include_directories(tet_exodus_host BEFORE PRIVATE
            "$<TARGET_PROPERTY:Kokkos::kokkos,INTERFACE_INCLUDE_DIRECTORIES>")
    endif()
    target_link_libraries(tet_exodus_host PRIVATE Mars::mars MPI::MPI_CXX "${NETCDF_LIBRARY}")
    foreach(ranks IN LISTS MARS_AMR_TEST_RANKS)
        add_test(NAME tet_exodus_host_mpi${ranks}
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks}
                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:tet_exodus_host> ${MPIEXEC_POSTFLAGS} "${exodus_fixture}")
        set_tests_properties(tet_exodus_host_mpi${ranks} PROPERTIES
            LABELS "amr;tet-exodus;host;mpi;tet-boundary" FIXTURES_REQUIRED tet_exodus
            PROCESSORS ${ranks} TIMEOUT 120)
    endforeach()
endif()

if(MARS_AMR_GPU_CHECKS)
    set(gpu_checks tet_gpu_driver.cu tet_gpu_mesh.cu tet_gpu_transfer.cu tet_gpu_dependent.cu
                   tet_gpu_refine.cu tet_gpu_red.cu tet_gpu_manager.cu)
    set_source_files_properties(${gpu_checks} PROPERTIES LANGUAGE ${MARS_GPU_LANGUAGE})
    add_executable(tet_bisection_gpu ${gpu_checks} tet_bisection_gpu.cpp tet_reference.cpp)
    target_compile_features(tet_bisection_gpu PRIVATE cxx_std_20)
    target_link_libraries(tet_bisection_gpu PRIVATE Mars::mars Mars::mars_unstructured MPI::MPI_CXX)
    if(MARS_GPU_LANGUAGE STREQUAL "CUDA")
        find_package(Python3 REQUIRED COMPONENTS Interpreter)
        set_target_properties(tet_bisection_gpu PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
        set_property(TARGET tet_bisection_gpu PROPERTY CUDA_COMPILER_LAUNCHER
            "${Python3_EXECUTABLE};${CMAKE_CURRENT_SOURCE_DIR}/compile_cuda.py;${CMAKE_CURRENT_BINARY_DIR}/compile-receipts")
    endif()

    foreach(ranks IN LISTS MARS_AMR_TEST_RANKS)
        foreach(mode IN ITEMS needed full staged)
            set(mode_flag)
            if(mode STREQUAL "full")
                set(mode_flag --full-halo)
            elseif(mode STREQUAL "staged")
                set(mode_flag --staged)
            endif()
            add_test(NAME tet_lepp_gpu_${mode}_mpi${ranks}
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks}
                    ${MPIEXEC_PREFLAGS} $<TARGET_FILE:tet_bisection_gpu> ${MPIEXEC_POSTFLAGS} ${mode_flag})
            set_tests_properties(tet_lepp_gpu_${mode}_mpi${ranks} PROPERTIES
                LABELS "amr;tet-lepp;gpu;mpi;tet-conformity;tet-transfer;tet-boundary;${mode}"
                PROCESSORS ${ranks} TIMEOUT 300)
            if(MARS_AMR_EXODUS_CHECKS)
                add_test(NAME tet_exodus_gpu_${mode}_mpi${ranks}
                    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks}
                        ${MPIEXEC_PREFLAGS} $<TARGET_FILE:tet_bisection_gpu> ${MPIEXEC_POSTFLAGS}
                        ${mode_flag} --exodus "${exodus_fixture}")
                set_tests_properties(tet_exodus_gpu_${mode}_mpi${ranks} PROPERTIES
                    LABELS "amr;tet-exodus;gpu;mpi;tet-conformity;tet-transfer;tet-boundary;${mode}"
                    FIXTURES_REQUIRED tet_exodus PROCESSORS ${ranks} TIMEOUT 300)
            endif()
        endforeach()
    endforeach()
endif()
