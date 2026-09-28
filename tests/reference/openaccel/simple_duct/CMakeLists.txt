# Rectangular-duct SIMPLE validation: host CPU/MPI tests only. Shared MARS CMake is not modified;
# the CUDA runs use the unchanged production executable mars_segregated_simple (README.md).
#   cmake -S tests/reference/openaccel/simple_duct -B build-duct && cmake --build build-duct -j4
#   ctest --test-dir build-duct -LE long      # analytic, mesh, comparator, coarse 1/2/4-rank duct
#   ctest --test-dir build-duct -L long       # host refinement 4/8/16 on 1/2/4 ranks (hours on 4 cores)
cmake_minimum_required(VERSION 3.22)
project(mars_simple_duct LANGUAGES CXX)
enable_testing()
find_package(MPI REQUIRED COMPONENTS CXX)
find_package(Python3 REQUIRED COMPONENTS Interpreter)
get_filename_component(_duct_root "${CMAKE_CURRENT_SOURCE_DIR}/../../../.." ABSOLUTE)
set(_duct_segregated "${_duct_root}/backend/distributed/unstructured/fem/segregated")
set(_duct "${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(duct_analytic_check analytic_check.cpp)
target_compile_features(duct_analytic_check PRIVATE cxx_std_20)

add_executable(duct_host_run duct_host_run.cpp)
target_include_directories(duct_host_run PRIVATE "${_duct_segregated}" "${_duct}")
target_compile_features(duct_host_run PRIVATE cxx_std_20)
target_compile_definitions(duct_host_run PRIVATE OMPI_SKIP_MPICXX=1 MPICH_SKIP_MPICXX=1)
target_link_libraries(duct_host_run PRIVATE MPI::MPI_CXX)

add_test(NAME ductAnalyticCxx COMMAND duct_analytic_check)
add_test(NAME ductAnalyticPython COMMAND ${Python3_EXECUTABLE} "${_duct}/test_duct_analytic.py" --cxx $<TARGET_FILE:duct_analytic_check>)
add_test(NAME ductMesh COMMAND ${Python3_EXECUTABLE} "${_duct}/test_duct_mesh.py" --cxx $<TARGET_FILE:duct_host_run>)
add_test(NAME ductComparator COMMAND ${Python3_EXECUTABLE} "${_duct}/test_duct_compare.py")
set_tests_properties(ductAnalyticCxx ductAnalyticPython ductMesh ductComparator PROPERTIES TIMEOUT 300)

# The production native reader on generated Exodus files, where netCDF is available.
find_path(MARS_DUCT_NETCDF_INCLUDE_DIR netcdf.h HINTS ${NETCDF_INCLUDE_DIRS})
find_library(MARS_DUCT_NETCDF_LIBRARY netcdf HINTS ${NETCDF_LIBRARY_DIRS})
if(MARS_DUCT_NETCDF_INCLUDE_DIR AND MARS_DUCT_NETCDF_LIBRARY)
    add_executable(duct_exodus_check duct_exodus_check.cpp)
    target_include_directories(duct_exodus_check PRIVATE "${_duct_segregated}" "${_duct}" "${MARS_DUCT_NETCDF_INCLUDE_DIR}")
    target_compile_features(duct_exodus_check PRIVATE cxx_std_20)
    target_compile_definitions(duct_exodus_check PRIVATE MARS_HAVE_NETCDF OMPI_SKIP_MPICXX=1 MPICH_SKIP_MPICXX=1)
    target_link_libraries(duct_exodus_check PRIVATE MPI::MPI_CXX "${MARS_DUCT_NETCDF_LIBRARY}")
    add_test(NAME ductExodusReader COMMAND ${Python3_EXECUTABLE} "${_duct}/test_duct_exodus.py" $<TARGET_FILE:duct_exodus_check>
             ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_PREFLAGS})
    set_tests_properties(ductExodusReader PROPERTIES TIMEOUT 120 PROCESSORS 2)
else()
    message(STATUS "simple_duct: netCDF not found; ductExodusReader is not built")
endif()

add_test(NAME ductHostStudyScript COMMAND ${Python3_EXECUTABLE} "${_duct}/test_run_host_study.py")
set_tests_properties(ductHostStudyScript PROPERTIES TIMEOUT 300)

# Production SIMPLE on the duct, host build: coarse 1/2/4-rank parity (a few minutes), and the
# refinement study (label long). Results persist in the build tree with manifests (executable
# SHA-256, arguments, mesh); they are reused only when the manifest matches and are rerun after
# a rebuild that changes duct_host_run (--rerun-stale).
function(_duct_study name levels ranks)
    add_test(NAME ${name} COMMAND ${Python3_EXECUTABLE} "${_duct}/run_host_study.py" "${CMAKE_CURRENT_BINARY_DIR}/${name}"
             --host-run $<TARGET_FILE:duct_host_run> --levels ${levels} --ranks ${ranks} --rerun-stale
             "--mpiexec=${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}" "--mpi-args=${MPIEXEC_PREFLAGS}" ${ARGN})
endfunction()
_duct_study(ductHostParity 4 1,2,4 --parity-only)
set_tests_properties(ductHostParity PROPERTIES TIMEOUT 1800 PROCESSORS 4)
_duct_study(ductHostRefinement 4,8,16 1,2,4)
set_tests_properties(ductHostRefinement PROPERTIES TIMEOUT 28800 PROCESSORS 4 LABELS long)
