# find_package(
#     VTK
#     COMPONENTS vtkCommonCore vtkCommonDataModel vtkFiltersGeneral vtkIOXML vtkIOParallel vtkIOParallelXML
#     REQUIRED)

# list(APPEND MARS_VTK_MODULES .)

# find_project_files(CMAKE_CURRENT_SOURCE_DIR MARS_VTK_MODULES LOCAL_HEADERS
#  LOCAL_SOURCES)

# add_library(mars_vtk STATIC ${LOCAL_HEADERS} ${LOCAL_SOURCES})
# target_link_libraries(mars_vtk PUBLIC mars_core mars_mpi)
# target_include_directories(mars_vtk PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
# target_include_directories(mars_vtk PUBLIC ${CMAKE_BINARY_DIR})
# target_include_directories(mars_vtk PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
# target_include_directories(mars_vtk PUBLIC ${VTK_INCLUDES})
# target_include_directories(mars_vtk BEFORE PUBLIC ${VTK_INCLUDES})

# if(Kokkos_FOUND)
#     target_include_directories(mars_vtk PUBLIC ${Kokkos_TPL_INCLUDE_DIRS}
#       ${Kokkos_INCLUDE_DIRS})
# endif()

# message(STATUS "VTK_VERSION: ${VTK_VERSION}")

# if(VTK_VERSION VERSION_LESS "8.90.0")
#     include(${VTK_USE_FILE})
#     target_link_libraries(mars_vtk PUBLIC ${VTK_LIBRARIES})
# else()
#     target_link_libraries(mars_vtk PUBLIC ${VTK_LIBRARIES})
#     vtk_module_autoinit(TARGETS mars_vtk MODULES ${VTK_LIBRARIES})
# endif()


# set(MARS_WITH_VTK TRUE)
# set(MARS_WITH_VTK TRUE PARENT_SCOPE)


# set(MARS_LIBRARIES
#     "${MARS_LIBRARIES};mars_vtk"
#     PARENT_SCOPE)

# set(MARS_HEADERS
#     "${MARS_HEADERS};${LOCAL_HEADERS}"
#     PARENT_SCOPE)

scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "." MARS_BUILD_INCLUDES MARS_HEADERS MARS_SOURCES)

set(MARS_BUILD_INCLUDES
    ${MARS_BUILD_INCLUDES}
    PARENT_SCOPE)

set(MARS_HEADERS
    ${MARS_HEADERS}
    PARENT_SCOPE)

set(MARS_SOURCES
    ${MARS_SOURCES}
    PARENT_SCOPE)