if(Kokkos_FOUND)

list(APPEND MARS_EXAMPLES distributed_dm unstructured)

find_project_files(CMAKE_CURRENT_SOURCE_DIR MARS_EXAMPLES LOCAL_HEADERS
                   LOCAL_SOURCES)

add_library(mars_examples STATIC ${LOCAL_HEADERS} ${LOCAL_SOURCES})
target_compile_features(mars_examples PUBLIC cxx_std_14)
target_link_libraries(mars_examples mars_kokkos mars_core mars_mpi)

target_link_libraries(mars_examples mars_io)

target_include_directories(mars_examples PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(mars_examples PUBLIC ${CMAKE_BINARY_DIR})
target_include_directories(mars_examples PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)

target_include_directories(mars_examples
                           PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/unstructured)

if(MPI_FOUND)

target_include_directories(mars_examples
                           PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/distributed_dm)
target_include_directories(mars_examples
                           PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/distributed_dm/staggered)
endif()

foreach(MODULE ${MARS_EXAMPLES})
    target_include_directories(mars_examples BEFORE PUBLIC ${MODULE})
endforeach(MODULE)

set(MARS_LIBRARIES
    "${MARS_LIBRARIES};mars_examples;${DISTRIBUTED_LIBRARIES}"
    PARENT_SCOPE)
set(MARS_INCLUDES
    "${MARS_INCLUDES};${DISTRIBUTED_INCLUDES}"
    PARENT_SCOPE)

endif()
