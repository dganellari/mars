set(MARS_STRUCTURED_TEST_RANKS "2;4" CACHE STRING "MPI rank counts for structured Kokkos tests")

function(mars_add_structured_test name source capability distributed)
    if(MARS_GPU_LANGUAGE)
        set_source_files_properties(${source} PROPERTIES LANGUAGE ${MARS_GPU_LANGUAGE})
    endif()
    add_executable(${name} ${source})
    if(MARS_ENABLE_CUDA AND Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE)
        set_target_properties(${name} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    endif()
    target_link_libraries(${name} PRIVATE mars ${MARS_G_TEST_LIBRARIES})
    target_include_directories(${name} PRIVATE
        ${PROJECT_SOURCE_DIR}/backend/distributed/testsuite/common
        ${PROJECT_SOURCE_DIR}/tests/structured)
    if(MARS_GPU_LANGUAGE)
        if(NOT MPIEXEC_EXECUTABLE)
            message(FATAL_ERROR "MPIEXEC_EXECUTABLE is required for structured GPU tests, including one-rank cases.")
        endif()
        add_test(NAME ${name}
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${name}> ${MPIEXEC_POSTFLAGS})
    else()
        add_test(NAME ${name} COMMAND ${name})
    endif()
    set_tests_properties(${name} PROPERTIES LABELS "structured;kokkos;structured-${capability}" PROCESSORS 1 TIMEOUT 120)
    if(distributed)
        foreach(ranks IN LISTS MARS_STRUCTURED_TEST_RANKS)
            if(NOT ranks MATCHES "^[1-9][0-9]*$")
                message(FATAL_ERROR "MARS_STRUCTURED_TEST_RANKS must contain positive integers.")
            endif()
            if(NOT MPIEXEC_EXECUTABLE)
                message(FATAL_ERROR "MPIEXEC_EXECUTABLE is required for structured MPI tests.")
            endif()
            add_test(NAME ${name}_mpi${ranks}
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks}
                    ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${name}> ${MPIEXEC_POSTFLAGS})
            set_tests_properties(${name}_mpi${ranks} PROPERTIES
                LABELS "structured;kokkos;structured-${capability};mpi" PROCESSORS ${ranks} TIMEOUT 120)
        endforeach()
    endif()
    set_property(GLOBAL APPEND PROPERTY MARS_TESTS ${name})
    set_property(GLOBAL APPEND PROPERTY MARS_TEST_TARGETS ${name})
endfunction()
