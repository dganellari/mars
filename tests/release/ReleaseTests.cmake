# Release checks: end-to-end runs of the documented v0.1.0 drivers on meshes generated at test
# time, so a fresh clone can verify its GPU build without any mesh download.
#
#   ctest -L release          (needs a GPU, MPI, python3 + numpy; runs in a few minutes)
#
# Included from examples/distributed/unstructured/CMakeLists.txt, so the driver targets exist.
# MARS_RELEASE_TEST_RANKS sets the multi-rank count (ranks share GPUs if there are fewer).
# Launcher flags go in MPIEXEC_PREFLAGS (a CMake list), e.g. on a Slurm cluster:
#   -DMPIEXEC_EXECUTABLE=$(which srun) "-DMPIEXEC_PREFLAGS=--account=X;--time=00:10:00;--export=ALL"

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
set(_rel_hexs ${_rel_dir}/hex16_x100)
set(_rel_hexc ${_rel_dir}/hex16_cycled)
set(_rel_mpi  ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG})
set(_rel_np   ${MARS_RELEASE_TEST_RANKS})

# --- fixture: generate the meshes once --------------------------------------------------------
add_test(NAME marsReleaseMeshHex
         COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/generate_hex_cube.py
                 --nx 16 --ny 16 --nz 16 --output ${_rel_hex})
add_test(NAME marsReleaseMeshTet
         COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/generate_tet_cube.py
                 --nx 8 --ny 8 --nz 8 --output ${_rel_tet})
# Same cube with edges of 6.25 instead of 1/16: halo coverage must not depend on length units.
add_test(NAME marsReleaseMeshHexScaled
         COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/generate_hex_cube.py
                 --nx 16 --ny 16 --nz 16 --scale 100 --output ${_rel_hexs})
# Same cube with each element's reference axes along (y, z, x): catches a wrong gradient transform,
# which axis-aligned numbering hides because the Jacobian is then diagonal.
add_test(NAME marsReleaseMeshHexCycled
         COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/generate_hex_cube.py
                 --nx 16 --ny 16 --nz 16 --cycle-axes --output ${_rel_hexc})
set_tests_properties(marsReleaseMeshHex marsReleaseMeshTet marsReleaseMeshHexScaled marsReleaseMeshHexCycled PROPERTIES
    FIXTURES_SETUP marsReleaseMeshes LABELS "release" TIMEOUT 120)

# --- assembly: same norms AND the same owned rows (by node SFC key) on 1 rank and on N ranks ---
# MPIEXEC_PREFLAGS is a CMake list; the checker takes the flags as one space-separated string.
string(REPLACE ";" " " _rel_preflags "${MPIEXEC_PREFLAGS}")
set(_rel_check ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/release/check_rank_invariance.py
               --np ${_rel_np} "--mpiexec=${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
               "--preflags=${_rel_preflags}")
set(_rel_rows ${CMAKE_BINARY_DIR}/release_rows)
add_test(NAME marsReleaseHexAssembly
         COMMAND ${_rel_check} --rows=${_rel_rows}/hex --
                 $<TARGET_FILE:mars_cvfem_graph> --mesh=${_rel_hex} --iterations=1 --quiet)
add_test(NAME marsReleaseHexAssemblyScaled
         COMMAND ${_rel_check} --rows=${_rel_rows}/hex_x100 --
                 $<TARGET_FILE:mars_cvfem_graph> --mesh=${_rel_hexs} --iterations=1 --quiet)
add_test(NAME marsReleaseHexNumbering
         COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/release/check_same_rows.py
                 --rows=${_rel_rows}/hex_numbering --mesh-a=${_rel_hex} --mesh-b=${_rel_hexc} --
                 ${_rel_mpi} 1 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:mars_cvfem_graph> ${MPIEXEC_POSTFLAGS}
                 --mesh={mesh} --iterations=1 --quiet)
add_test(NAME marsReleaseTetAssembly
         COMMAND ${_rel_check} --rows=${_rel_rows}/tet --
                 $<TARGET_FILE:mars_cvfem_graph_tet> --mesh=${_rel_tet} --iterations=1 --quiet)

# --- solve: CVFEM Poisson, -Δu = 1 on the unit cube, u = 0 on the boundary, 1 and N ranks ---------
# Exact centre value 0.05621; the discrete maximum on the 16^3 mesh must be within 5% of it.
add_test(NAME marsReleasePoisson
         COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/release/check_value.py
                 "--regex=Max:\\s*([-+0-9.eE]+)" --lo 0.0534 --hi 0.0590 --
                 ${_rel_mpi} 1 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:mars_cvfem_poisson> ${MPIEXEC_POSTFLAGS}
                 --mesh=${_rel_hex})
add_test(NAME marsReleasePoissonCycled
         COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/release/check_value.py
                 "--regex=Max:\\s*([-+0-9.eE]+)" --lo 0.0534 --hi 0.0590 --
                 ${_rel_mpi} 1 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:mars_cvfem_poisson> ${MPIEXEC_POSTFLAGS}
                 --mesh=${_rel_hexc})
add_test(NAME marsReleasePoisson_np${_rel_np}
         COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/release/check_value.py
                 "--regex=Max:\\s*([-+0-9.eE]+)" --lo 0.0534 --hi 0.0590 --
                 ${_rel_mpi} ${_rel_np} ${MPIEXEC_PREFLAGS} $<TARGET_FILE:mars_cvfem_poisson> ${MPIEXEC_POSTFLAGS}
                 --mesh=${_rel_hex})

# Galerkin P1 Poisson example on the tet cube, 1 and N ranks: same band (exact max 0.05621)
foreach(_np 1 ${_rel_np})
    add_test(NAME marsReleaseEx1Poisson_np${_np}
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/release/check_value.py
                     "--regex=Max:\\s*([-+0-9.eE]+)" --lo 0.0534 --hi 0.0590 --
                     ${_rel_mpi} ${_np} ${MPIEXEC_PREFLAGS} $<TARGET_FILE:mars_ex1_poisson> ${MPIEXEC_POSTFLAGS}
                     --mesh=${_rel_tet})
endforeach()

# --- Navier-Stokes (fem/mars_navier_stokes.hpp, needs Hypre), 1 and N ranks -------------------
# The examples exit non-zero on any failed linear solve (all ranks stop together).
if(MARS_ENABLE_HYPRE)
    foreach(_np 1 ${_rel_np})
        # Closed box: the pressure null space and the unreachable edge nodes.
        add_test(NAME marsReleaseCavity_np${_np}
                 COMMAND ${_rel_mpi} ${_np} ${MPIEXEC_PREFLAGS} $<TARGET_FILE:mars_lid_driven_cavity>
                         ${MPIEXEC_POSTFLAGS} --mesh=${_rel_hex} --num-steps=10 --report-every=5)
        # Planar channel with an inlet and an outlet, generated on every rank.
        add_test(NAME marsReleaseChannel_np${_np}
                 COMMAND ${_rel_mpi} ${_np} ${MPIEXEC_PREFLAGS} $<TARGET_FILE:mars_poiseuille_flow>
                         ${MPIEXEC_POSTFLAGS} --cells=200,40 --num-steps=10 --report-every=5)
        # Taylor-Green vortex in the periodic unit cube: low-Re viscous decay. After 100 steps
        # KE / KE_Stokes is 1.0015; a broken periodic coupling on N ranks moves it far outside the band.
        add_test(NAME marsReleaseTgv_np${_np}
                 COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/release/check_value.py
                         "--regex=TGV final:.*KE/KE_Stokes=([-+0-9.eE]+)" --lo 1.0013 --hi 1.0017 --
                         ${_rel_mpi} ${_np} ${MPIEXEC_PREFLAGS} $<TARGET_FILE:mars_tgv> ${MPIEXEC_POSTFLAGS}
                         --mesh=${_rel_hex} --box-lo=0 --box-hi=1 --nu=0.05 --dt=1e-4 --num-steps=100
                         --report-every=100)
        foreach(_t Cavity Channel Tgv)
            set_tests_properties(marsRelease${_t}_np${_np} PROPERTIES FAIL_REGULAR_EXPRESSION "[=: ](-?nan|NaN)[ ,\n]")
        endforeach()
    endforeach()
endif()

# --- high-order matrix-free operator (experimental); both drivers generate their own cube -------
# One GPU, p = 1..7: GPU metric and element apply against the host reference, A*1 = 0,
# A*linear = 0 at interior DOFs, and a sheared-element gate. Exits 1 if any gate fails.
if(TARGET mars_cvfem_ho_matfree_test)
    add_test(NAME marsReleaseHoMatfree
             COMMAND ${_rel_mpi} 1 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:mars_cvfem_ho_matfree_test> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(marsReleaseHoMatfree PROPERTIES LABELS "release;gpu" TIMEOUT 600)
endif()
# N ranks, p = 3 (edge and face DOFs with more than one node, so their orientation matters):
# one distributed matvec of u = 1 must vanish on every owned DOF. The driver exits 0 even when
# a gate fails, so check_value reads the device gate and the regex catches the host gate.
if(TARGET mars_ho_dist_apply_test)
    add_test(NAME marsReleaseHoDistApply_np${_rel_np}
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/release/check_value.py
                     "--regex=device A\\.1 = ([-+0-9.eE]+) \\[" --lo 0 --hi 1e-8 --
                     ${_rel_mpi} ${_rel_np} ${MPIEXEC_PREFLAGS} $<TARGET_FILE:mars_ho_dist_apply_test>
                     ${MPIEXEC_POSTFLAGS} --ncells=16 --p=3)
    set_tests_properties(marsReleaseHoDistApply_np${_rel_np} PROPERTIES
        LABELS "release;gpu" TIMEOUT 600 FAIL_REGULAR_EXPRESSION "\\[FAIL\\]|CUDA error")
endif()

get_property(_rel_tests DIRECTORY PROPERTY TESTS)
list(FILTER _rel_tests INCLUDE REGEX "^marsRelease(Hex|Tet|Poisson|Ex1Poisson|Cavity|Channel|Tgv)")
set_tests_properties(${_rel_tests} PROPERTIES
    FIXTURES_REQUIRED marsReleaseMeshes LABELS "release;gpu" TIMEOUT 600)
