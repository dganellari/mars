find_package(ParMoonolith)

if(MOONOLITH_FOUND)

    message(STATUS "MOONOLITH_INCLUDES=${MOONOLITH_INCLUDES}")

    list(APPEND MARS_MOONOLITH_MODULES .)

    find_project_files(CMAKE_CURRENT_SOURCE_DIR MARS_MOONOLITH_MODULES
                       LOCAL_HEADERS LOCAL_SOURCES)

    add_library(mars_moonolith STATIC ${LOCAL_HEADERS} ${LOCAL_SOURCES})
    target_compile_features(mars_moonolith PUBLIC cxx_std_14)
    target_link_libraries(mars_moonolith ${MOONOLITH_LIBRARIES})
    add_dependencies(mars_moonolith par_moonolith)
    target_include_directories(mars_moonolith
                               PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
    target_include_directories(mars_moonolith
                               PUBLIC ${CMAKE_SOURCE_DIR}/generation)

    target_include_directories(mars_moonolith PUBLIC ${CMAKE_BINARY_DIR})
    target_include_directories(mars_moonolith
                               PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
    target_include_directories(mars_moonolith PUBLIC ${MOONOLITH_INCLUDES})
    target_include_directories(mars_moonolith BEFORE
                               PUBLIC ${MOONOLITH_INCLUDES})

    if(Kokkos_FOUND)
        target_include_directories(
            mars_moonolith PUBLIC ${Kokkos_TPL_INCLUDE_DIRS}
                                  ${Kokkos_INCLUDE_DIRS})
    endif()

    set(MARS_LIBRARIES
        "${MARS_LIBRARIES};mars_moonolith;${MOONOLITH_LIBRARIES}"
        PARENT_SCOPE)
    set(MARS_INCLUDES
        "${MARS_INCLUDES};${MOONOLITH_INCLUDES}"
        PARENT_SCOPE)
    set(MOONOLITH_FOUND
        ${MOONOLITH_FOUND}
        PARENT_SCOPE)
    set(WITH_PAR_MOONOLITH
        ON
        PARENT_SCOPE)

    set(CMAKE_CXX_COMPILER ${MOONOLITH_CXX_COMPILER})
    set(CMAKE_C_COMPILER ${MOONOLITH_C_COMPILER})

    set(CMAKE_CXX_COMPILER
        ${MOONOLITH_CXX_COMPILER}
        PARENT_SCOPE)
    set(CMAKE_C_COMPILER
        ${MOONOLITH_C_COMPILER}
        PARENT_SCOPE)

endif()
