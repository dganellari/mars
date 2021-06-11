# find_package(
#     VTK
#     COMPONENTS vtkCommonCore vtkCommonDataModel vtkFiltersGeneral vtkIOXML
#                # vtkInteractionStyle vtkRenderingCore vtkRenderingFreeType
#                # vtkRenderingOpenGL2
#     QUIET)
# if(NOT VTK_FOUND)
    # message("Skipping IO: ${VTK_NOT_FOUND_MESSAGE}")
    add_library(mars_io INTERFACE)
    target_include_directories(mars_io INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

# else(VTK_FOUND)

#     list(APPEND MARS_IO_MODULES .)

#     find_project_files(CMAKE_CURRENT_SOURCE_DIR MARS_IO_MODULES LOCAL_HEADERS
#                        LOCAL_SOURCES)

#     add_library(mars_io STATIC ${LOCAL_HEADERS} ${LOCAL_SOURCES})
#     target_include_directories(mars_io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
#     target_include_directories(mars_io PUBLIC ${CMAKE_BINARY_DIR})
#     target_include_directories(mars_io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
#     target_include_directories(mars_io PUBLIC ${VTK_INCLUDES})
#     target_include_directories(mars_io BEFORE PUBLIC ${VTK_INCLUDES})

#     if(Kokkos_FOUND)
#         target_include_directories(mars_io PUBLIC ${Kokkos_TPL_INCLUDE_DIRS}
#                                                   ${Kokkos_INCLUDE_DIRS})
#     endif()

#     message(STATUS "VTK_VERSION: ${VTK_VERSION}")
#     if(VTK_VERSION VERSION_LESS "8.90.0")
#         # old system
#         include(${VTK_USE_FILE})
#         # add_executable(WriteVTU MACOSX_BUNDLE WriteVTU.cxx)
#         target_link_libraries(mars_io PUBLIC ${VTK_LIBRARIES})
#     else()
#         # include all components add_executable(mars_io MACOSX_BUNDLE
#         # mars_io.cxx)
#         target_link_libraries(mars_io PUBLIC ${VTK_LIBRARIES})
#         # vtk_module_autoinit is needed
#         vtk_module_autoinit(TARGETS mars_io MODULES ${VTK_LIBRARIES})
#     endif()

#     set(MARS_WITH_VTK
#         ON
#         PARENT_SCOPE)

# endif()

set(MARS_LIBRARIES
    "${MARS_LIBRARIES};mars_io"
    PARENT_SCOPE)
# set(MARS_INCLUDES "${MARS_INCLUDES};${VTK_INCLUDES}" PARENT_SCOPE)

set(MARS_HEADERS
    "${MARS_HEADERS};${LOCAL_HEADERS}"
     PARENT_SCOPE)


