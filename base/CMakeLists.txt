list(APPEND MARS_BASE_MODULES .)

find_project_files(CMAKE_CURRENT_SOURCE_DIR MARS_BASE_MODULES LOCAL_HEADERS
                   LOCAL_SOURCES)

add_library(mars_serial INTERFACE)
# target_compile_features(mars_serial INTERFACE cxx_std_14)
target_link_libraries(mars_serial INTERFACE ${BASE_LIBRARIES})

target_include_directories(mars_serial INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(mars_serial INTERFACE ${CMAKE_BINARY_DIR})
target_include_directories(mars_serial INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_include_directories(mars_serial INTERFACE ${BASE_INCLUDES})
target_include_directories(mars_serial BEFORE INTERFACE ${BASE_INCLUDES})

target_include_directories(mars_serial INTERFACE ${Kokkos_TPL_INCLUDE_DIRS}
                                                 ${Kokkos_INCLUDE_DIRS})
set(MARS_LIBRARIES
    "${MARS_LIBRARIES};mars_serial;${BASE_LIBRARIES}"
    PARENT_SCOPE)
set(MARS_INCLUDES
    "${MARS_INCLUDES};${BASE_INCLUDES}"
    PARENT_SCOPE)
