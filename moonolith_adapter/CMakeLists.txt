find_package(ParMoonolith)


if(MOONOLITH_FOUND)

    list(APPEND MARS_MOONOLITH_MODULES 
        .
    )

    find_project_files(CMAKE_CURRENT_SOURCE_DIR MARS_MOONOLITH_MODULES LOCAL_HEADERS LOCAL_SOURCES)

    add_library(mars_moonolith STATIC ${LOCAL_HEADERS} ${LOCAL_SOURCES})
    target_link_libraries(mars_moonolith ${MOONOLITH_LIBRARIES})
    target_include_directories(mars_moonolith PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
    target_include_directories(mars_moonolith PUBLIC ${CMAKE_BINARY_DIR})
    target_include_directories(mars_moonolith PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
    target_include_directories(mars_moonolith PUBLIC ${MOONOLITH_INCLUDES})

    if(Kokkos_FOUND)
        target_include_directories(mars_moonolith PUBLIC ${Kokkos_TPL_INCLUDE_DIRS} ${Kokkos_INCLUDE_DIRS})
    endif()

    set(MARS_LIBRARIES "${MARS_LIBRARIES};mars_moonolith" PARENT_SCOPE)
    set(MOONOLITH_FOUND ${MOONOLITH_FOUND} PARENT_SCOPE)
    set(WITH_PAR_MOONOLITH ON PARENT_SCOPE)
endif()
