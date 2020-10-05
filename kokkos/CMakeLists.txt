if(Kokkos_FOUND)

    list(APPEND MARS_KOKKOS_MODULES .)

    find_project_files(CMAKE_CURRENT_SOURCE_DIR MARS_KOKKOS_MODULES
                       LOCAL_HEADERS LOCAL_SOURCES)

    add_library(mars_kokkos INTERFACE)
    # target_compile_features(mars_kokkos INTERFACE cxx_std_14)
    target_link_libraries(mars_kokkos INTERFACE ${KOKKOS_LIBRARIES})

    target_include_directories(mars_kokkos
                               INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
    target_include_directories(mars_kokkos INTERFACE ${CMAKE_BINARY_DIR})
    target_include_directories(mars_kokkos
                               INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/..)
    target_include_directories(mars_kokkos INTERFACE ${KOKKOS_INCLUDES})
    target_include_directories(mars_kokkos BEFORE INTERFACE ${KOKKOS_INCLUDES})

    target_include_directories(mars_kokkos INTERFACE ${Kokkos_TPL_INCLUDE_DIRS}
                                                     ${Kokkos_INCLUDE_DIRS})
    set(MARS_LIBRARIES
        "${MARS_LIBRARIES};mars_kokkos;${KOKKOS_LIBRARIES}"
        PARENT_SCOPE)
    set(MARS_INCLUDES
        "${MARS_INCLUDES};${KOKKOS_INCLUDES}"
        PARENT_SCOPE)

endif()
