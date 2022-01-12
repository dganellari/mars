list(
    APPEND
    MARS_ADIOS2_MODULES
    .)

find_project_files(CMAKE_CURRENT_SOURCE_DIR MARS_ADIOS2_MODULES LOCAL_HEADERS
                   LOCAL_SOURCES)


add_library(mars_adios2  ${LOCAL_HEADERS} ${LOCAL_SOURCES})
target_link_libraries(mars_adios2 PUBLIC adios2::adios2 MPI::MPI_C)
target_link_libraries(mars_adios2 PUBLIC mars_mpi)
target_include_directories(mars_adios2 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(MARS_WITH_ADIOS2 TRUE)
set(MARS_WITH_ADIOS2 TRUE PARENT_SCOPE)


set(MARS_LIBRARIES
    "${MARS_LIBRARIES};mars_adios2"
    PARENT_SCOPE)

set(MARS_HEADERS
    "${MARS_HEADERS};${LOCAL_HEADERS}"
     PARENT_SCOPE)

