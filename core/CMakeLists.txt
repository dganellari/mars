list(APPEND MARS_CORE_MODULES .)

find_project_files(CMAKE_CURRENT_SOURCE_DIR MARS_CORE_MODULES LOCAL_HEADERS
                   LOCAL_SOURCES)

add_library(mars_core STATIC ${LOCAL_HEADERS} ${LOCAL_SOURCES})
target_compile_features(mars_core PUBLIC cxx_std_14)
target_link_libraries(mars_core mars_serial)

target_include_directories(mars_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(mars_core PUBLIC ${CMAKE_BINARY_DIR})
target_include_directories(mars_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_include_directories(mars_core PUBLIC ${CORE_INCLUDES})
target_include_directories(mars_core BEFORE PUBLIC ${CORE_INCLUDES})

if(Kokkos_FOUND)
    target_include_directories(mars_core PUBLIC ${Kokkos_TPL_INCLUDE_DIRS}
                                                ${Kokkos_INCLUDE_DIRS})
endif()

set(MARS_LIBRARIES
    "${MARS_LIBRARIES};mars_core;${CORE_LIBRARIES}"
    PARENT_SCOPE)
set(MARS_INCLUDES
    "${MARS_INCLUDES};${CORE_INCLUDES}"
    PARENT_SCOPE)
set(MARS_HEADERS
    "${MARS_HEADERS};${LOCAL_HEADERS}"
    PARENT_SCOPE)

install(
    TARGETS mars_core
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)
