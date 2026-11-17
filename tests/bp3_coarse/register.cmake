# Launch MPI from CTest so the fixture owns its process-manager connection.
if(NOT MPIEXEC_EXECUTABLE)
    message(FATAL_ERROR "The BP3 coarse test requires an MPI launcher")
endif()
add_test(NAME marsBP3CoarseSingle
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:mars_bp3_coarse_test> ${MPIEXEC_POSTFLAGS} balanced 511)
set_tests_properties(marsBP3CoarseSingle PROPERTIES LABELS "bp3-coarse" TIMEOUT 180
    ENVIRONMENT "BP3_CONSOLIDATE_ROWS=100;BP3_FIXTURE_PROFILE=mars;MPICH_GPU_SUPPORT_ENABLED=1")
