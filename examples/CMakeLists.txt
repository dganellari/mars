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

if(Kokkos_FOUND)
    add_executable(st_example mars_st_example.cpp)
    target_link_libraries(st_example PRIVATE mars)
endif()