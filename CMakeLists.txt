cmake_minimum_required(VERSION 3.5...3.19)
cmake_policy(SET CMP0048 NEW)
project(
  mars
  VERSION 0.1.0
  LANGUAGES "CXX"
  HOMEPAGE_URL "https://bitbucket.org/zulianp/mars/src/master/"
  DESCRIPTION
    "Mars is an open-source mesh management library designed to handle N-dimensional elements (N <= 4)"
)

# Std CMake imports
include(CMakePackageConfigHelpers)
include(CMakeDependentOption)
# Cmake imports.
include(cmake/MarsDefaults.cmake)
include(cmake/MarsOptions.cmake)
include(cmake/MarsCMakeFunctions.cmake)
include(cmake/MarsDependencies.cmake)

# ##############################################################################
# ##############################################################################
# Library
# ##############################################################################

list(APPEND MARS_HEADERS mars.hpp mars_env.hpp)
list(APPEND MARS_SOURCES mars.cpp mars_env.cpp)

set(MARS_MODULES base core distributed examples io kokkos)

foreach(MODULE ${MARS_MODULES})
  add_subdirectory("${MODULE}")
endforeach(MODULE)

add_library(mars ${MARS_HEADERS} ${MARS_SOURCES})

# ##############################################################################
# include

foreach(include ${MARS_BUILD_INCLUDES})
  target_include_directories(mars BEFORE PUBLIC $<BUILD_INTERFACE:${include}>)
endforeach()

target_include_directories(
  mars BEFORE
  PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
         $<INSTALL_INTERFACE:include>)

target_include_directories(mars PUBLIC ${MARS_DEP_INCLUDES})

# ##############################################################################
# link

target_link_libraries(mars PUBLIC ${MARS_DEP_LIBRARIES})
# target_link_libraries(mars PUBLIC ${MPI_LIBRARIES})

# ##############################################################################
# properties

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel." FORCE)

  message(STATUS "[Status] CMAKE_BUILD_TYPE=Release")
endif(NOT CMAKE_BUILD_TYPE)

if(NOT DEFINED CMAKE_CXX_STANDARD_REQUIRED)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()
if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 14)
endif()

target_compile_features(mars PUBLIC cxx_std_14)

# ##############################################################################
# ##############################################################################
# Extra targets
# ##############################################################################

if(MARS_ENABLE_BENCHMARK)
  include(bench/benchmarks.cmake)
endif()

if(MARS_ENABLE_TESTING)
  include(tests/tests.cmake)
endif()

if(MARS_ENABLE_VTK)
  find_package(
    VTK
    COMPONENTS vtkCommonCore vtkCommonDataModel vtkFiltersGeneral vtkIOXML
               vtkIOParallel vtkIOParallelXML
    REQUIRED)
  add_subdirectory(backend/vtk)
endif()

if(MARS_ENABLE_ADIOS2)
  find_package(ADIOS2 REQUIRED)
  add_subdirectory(backend/adios2)
endif()


# ##############################################################################
# ##############################################################################
# ##############################################################################

# #################  EXECUTABLES AND INSTALL ##################

add_executable(mars_exec main.cpp)
target_link_libraries(mars_exec mars)

message(STATUS "${MARS_INCLUDES}")
include_directories(${CMAKE_BINARY_DIR})

install(

install(FILES ${MARS_HEADERS} DESTINATION include)
install(FILES ${CMAKE_BINARY_DIR}/mars_config.hpp DESTINATION include)

# ##############################################################################

endif()

message("MARS_DEP_INCLUDES: ${MARS_DEP_INCLUDES}")
message("MARS_DEP_LIBRARIES: ${MARS_DEP_LIBRARIES}")
