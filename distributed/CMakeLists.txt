# if(MPI_FOUND)


#     # Create an append list names MARS_MPI set to 0 if it were not defined before.
#     list(
#         APPEND
#         MARS_MPI
#         .
#         sfc
#         mars_env
#         kokkos
#         communication
#         communication/util)

#     find_project_files(CMAKE_CURRENT_SOURCE_DIR MARS_MPI LOCAL_HEADERS
#                        LOCAL_SOURCES)

#     add_library(mars_mpi STATIC ${LOCAL_HEADERS} ${LOCAL_SOURCES})
#     message("Headers:${LOCAL_HEADERS}, Sources: ${LOCAL_SOURCES}")
#     target_compile_features(mars_mpi PUBLIC cxx_std_14)

#     if(Kokkos_FOUND)
#         target_link_libraries(mars_mpi mars_kokkos mars_core)
#     else()
#         target_link_libraries(mars_mpi mars_core)
#     endif()

#     target_include_directories(mars_mpi PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
#     target_include_directories(mars_mpi PUBLIC ${CMAKE_BINARY_DIR})
#     target_include_directories(mars_mpi PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
#     target_include_directories(mars_mpi PUBLIC ${DISTRIBUTED_INCLUDES})
#     target_include_directories(mars_mpi BEFORE PUBLIC ${DISTRIBUTED_INCLUDES})
#     message("Boh:${CMAKE_CURRENT_SOURCE_DIR}")

#     foreach(MODULE ${MARS_MPI})
#         target_include_directories(mars_mpi BEFORE PUBLIC ${MODULE})
#     endforeach(MODULE)

#     set(MARS_LIBRARIES
#         "${MARS_LIBRARIES};mars_mpi;${DISTRIBUTED_LIBRARIES}"
#         PARENT_SCOPE)
#     set(MARS_INCLUDES
#         "${MARS_INCLUDES};${DISTRIBUTED_INCLUDES}"
#         PARENT_SCOPE)
#     set(MARS_HEADERS
#         "${MARS_HEADERS};${LOCAL_HEADERS}"
#         PARENT_SCOPE)

#     install(
#         TARGETS mars_mpi
#         RUNTIME DESTINATION bin
#         LIBRARY DESTINATION lib
#         ARCHIVE DESTINATION lib)
# endif()

list(APPEND DISTRIBUTED_MODULES 
    .
    communication
    communication/util
    kokkos
    mars_env
    sfc)


scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "${DISTRIBUTED_MODULES}" MARS_BUILD_INCLUDES MARS_HEADERS MARS_SOURCES)

set(MARS_BUILD_INCLUDES
    ${MARS_BUILD_INCLUDES}
    PARENT_SCOPE)

set(MARS_HEADERS
    ${MARS_HEADERS}
    PARENT_SCOPE)

set(MARS_SOURCES
    ${MARS_SOURCES}
    PARENT_SCOPE)
