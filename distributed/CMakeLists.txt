if(MPI_FOUND)

    list(
        APPEND
        MARS_MPI
        .
        sfc
        mars_env
        kokkos
        communication
        communication/util)

    find_project_files(CMAKE_CURRENT_SOURCE_DIR MARS_MPI LOCAL_HEADERS
                       LOCAL_SOURCES)

    add_library(mars_mpi STATIC ${LOCAL_HEADERS} ${LOCAL_SOURCES})
    target_compile_features(mars_mpi PUBLIC cxx_std_14)
    target_link_libraries(mars_mpi mars_kokkos mars_core)

    target_include_directories(mars_mpi PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
    target_include_directories(mars_mpi PUBLIC ${CMAKE_BINARY_DIR})
    target_include_directories(mars_mpi PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
    target_include_directories(mars_mpi PUBLIC ${DISTRIBUTED_INCLUDES})
    target_include_directories(mars_mpi BEFORE PUBLIC ${DISTRIBUTED_INCLUDES})

    foreach(MODULE ${MARS_MPI})
        target_include_directories(mars_mpi BEFORE PUBLIC ${MODULE})
    endforeach(MODULE)

    set(MARS_LIBRARIES
        "${MARS_LIBRARIES};mars_mpi;${DISTRIBUTED_LIBRARIES}"
        PARENT_SCOPE)
    set(MARS_INCLUDES
        "${MARS_INCLUDES};${DISTRIBUTED_INCLUDES}"
        PARENT_SCOPE)

endif()
