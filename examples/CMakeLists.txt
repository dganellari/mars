# if(Kokkos_FOUND)

# list(APPEND MARS_EXAMPLES distributed_dm unstructured)

# find_project_files(CMAKE_CURRENT_SOURCE_DIR MARS_EXAMPLES LOCAL_HEADERS
#                    LOCAL_SOURCES)

# add_library(mars_examples STATIC ${LOCAL_HEADERS} ${LOCAL_SOURCES})
# target_compile_features(mars_examples PUBLIC cxx_std_14)
# target_link_libraries(mars_examples mars_kokkos mars_core mars_mpi)

# target_link_libraries(mars_examples mars_io)

# target_include_directories(mars_examples PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
# target_include_directories(mars_examples PUBLIC ${CMAKE_BINARY_DIR})
# target_include_directories(mars_examples PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)

# target_include_directories(mars_examples
#                            PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/unstructured)

# if(MPI_FOUND)

# target_include_directories(mars_examples
#                            PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/distributed_dm)
# target_include_directories(mars_examples
#                            PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/distributed_dm/staggered)
# endif()


# foreach(MODULE ${MARS_EXAMPLES})
#     target_include_directories(mars_examples BEFORE PUBLIC ${MODULE})
# endforeach(MODULE)

# set(MARS_LIBRARIES
#     "${MARS_LIBRARIES};mars_examples;${DISTRIBUTED_LIBRARIES}"
#     PARENT_SCOPE)
# set(MARS_INCLUDES
#     "${MARS_INCLUDES};${DISTRIBUTED_INCLUDES}"
#     PARENT_SCOPE)
# set(MARS_HEADERS
#     "${MARS_HEADERS};${LOCAL_HEADERS}"
#      PARENT_SCOPE)


# endif()
list(APPEND EXAMPLES_MODULES 
    .
    unstructured
    distributed_dm
    distributed_dm/staggered)


scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "${EXAMPLES_MODULES}" MARS_BUILD_INCLUDES MARS_HEADERS MARS_SOURCES)

set(MARS_BUILD_INCLUDES
    ${MARS_BUILD_INCLUDES}
    PARENT_SCOPE)

set(MARS_HEADERS
    ${MARS_HEADERS}
    PARENT_SCOPE)

set(MARS_SOURCES
    ${MARS_SOURCES}
    PARENT_SCOPE)




# if(Kokkos_FOUND)
#     add_executable(st_example mars_st_example.cpp)
#     target_link_libraries(st_example PRIVATE mars mars_io)


#     add_custom_target(mars_examples)
#     # add_dependencies()
# endif()