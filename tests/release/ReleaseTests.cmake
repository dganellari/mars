# Release checks: end-to-end runs of the documented v0.1.0 drivers on meshes generated at test
# time, so a fresh clone can verify its GPU build without any mesh download.
#
#   ctest -L release          (needs a GPU, MPI, python3 + numpy; runs in a few minutes)
#
# Included from examples/distributed/unstructured/CMakeLists.txt, so the driver targets exist.
# MARS_RELEASE_TEST_RANKS sets the multi-rank count (ranks share GPUs if there are fewer).

find_package(Python3 COMPONENTS Interpreter)
if(NOT Python3_Interpreter_FOUND)
    message(STATUS "Release tests skipped: python3 not found (the meshes are generated with python3 + numpy)")
    return()
endif()
if(NOT MPIEXEC_EXECUTABLE)
    message(STATUS "Release tests skipped: no MPI launcher (MPIEXEC_EXECUTABLE) found")
    return()
endif()

set(MARS_RELEASE_TEST_RANKS 4 CACHE STRING "Rank count for the multi-rank release tests")

set(_rel_dir  ${CMAKE_BINARY_DIR}/release_meshes)
set(_rel_hex  ${_rel_dir}/hex16)
set(_rel_tet  ${_rel_dir}/tet8)
set(_rel_mpi  ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG})
set(_rel_np   ${MARS_RELEASE_TEST_RANKS})

# --- fixture: generate the meshes once --------------------------------------------------------
add_test(NAME marsReleaseMeshHex
         COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/generate_hex_cube.py
                 --nx 16 --ny 16 --nz 16 --output ${_rel_hex})
add_test(NAME marsReleaseMeshTet
         COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/generate_tet_cube.py
                 --nx 8 --ny 8 --nz 8 --output ${_rel_tet})
set_tests_properties(marsReleaseMeshHex marsReleaseMeshTet PROPERTIES
    FIXTURES_SETUP marsReleaseMeshes LABELS "release" TIMEOUT 120)

# --- assembly: same matrix/RHS norms on 1 rank and on N ranks ---------------------------------
set(_rel_check ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/release/check_rank_invariance.py
               --np ${_rel_np} "--mpiexec=${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
               "--preflags=${MPIEXEC_PREFLAGS}" --)
add_test(NAME marsReleaseHexAssembly
         COMMAND ${_rel_check} $<TARGET_FILE:mars_cvfem_graph> --mesh=${_rel_hex} --iterations=1 --quiet)
add_test(NAME marsReleaseTetAssembly
         COMMAND ${_rel_check} $<TARGET_FILE:mars_cvfem_graph_tet> --mesh=${_rel_tet} --iterations=1 --quiet)

# --- solve: CVFEM Poisson (single-rank driver) must converge ----------------------------------
add_test(NAME marsReleasePoisson
         COMMAND ${_rel_mpi} 1 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:mars_cvfem_poisson> ${MPIEXEC_POSTFLAGS}
                 --mesh=${_rel_hex})

# --- Navier-Stokes projection: cavity and channel, 1 and N ranks --------------------------------
# The driver exits non-zero on any failed linear solve (all ranks stop together).
foreach(_bc cavity channel)
    foreach(_np 1 ${_rel_np})
        add_test(NAME marsReleaseNs_${_bc}_np${_np}
                 COMMAND ${_rel_mpi} ${_np} ${MPIEXEC_PREFLAGS} $<TARGET_FILE:mars_amr_ns_projection>
                         ${MPIEXEC_POSTFLAGS} --mesh=${_rel_hex} --bc=${_bc} --num-steps=10)
        set_tests_properties(marsReleaseNs_${_bc}_np${_np} PROPERTIES FAIL_REGULAR_EXPRESSION "[=: ](-?nan|NaN)[ ,\n]")
    endforeach()
endforeach()

get_property(_rel_tests DIRECTORY PROPERTY TESTS)
list(FILTER _rel_tests INCLUDE REGEX "^marsRelease(Hex|Tet|Poisson|Ns)")
set_tests_properties(${_rel_tests} PROPERTIES
    FIXTURES_REQUIRED marsReleaseMeshes LABELS "release;gpu" TIMEOUT 600)
