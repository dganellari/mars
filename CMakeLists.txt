cmake_minimum_required(VERSION 3.5...3.19)
cmake_policy(SET CMP0048 NEW)

#Check if the compiler is set at command line.
#Needs to be before project(...) where the cxx compiler is set to the default value.
if(CMAKE_CXX_COMPILER)
    set(CMAKE_CXX_COMPILER_SET_EXTERNALLY TRUE)
endif()

project(
  mars
  VERSION 0.2.0
  LANGUAGES "CXX" "C"
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
include(cmake/MarsClangFormat.cmake)

if(MARS_ENABLE_CUDA)
    message(STATUS "${Magenta}For MARS with Cuda make sure to pass the option -DCMAKE_CXX_COMPILER=/path/to/nvcc_wrapper. Otherwise it will be set automatically to the Kokkos CXX Compiler.${ColourReset}")
endif()

# ##############################################################################
# ##############################################################################
# Library
# ##############################################################################

list(APPEND MARS_HEADERS mars.hpp mars_env.hpp)
list(APPEND MARS_SOURCES mars.cpp mars_env.cpp)

# set(MARS_MODULES base core distributed examples io kokkos)
set(MARS_MODULES base core examples io kokkos)

foreach(MODULE ${MARS_MODULES})
  add_subdirectory("${MODULE}")
endforeach(MODULE)

if(MARS_ENABLE_SERIAL_BACKEND)
    add_subdirectory("backend/serial")
endif()

if(MARS_ENABLE_DISTRIBUTED_BACKEND)
    add_subdirectory("backend/distributed")
endif()

if(MARS_ENABLE_AMR_BACKEND)
    add_subdirectory("backend/amr")
endif()


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



# target_compile_features(mars PUBLIC cxx_std_17)

# ##############################################################################
# ##############################################################################
# Extra targets
# ##############################################################################

if(MARS_ENABLE_TESTS)
  include(cmake/tests.cmake)
  add_subdirectory("backend/distributed/testsuite")
  add_subdirectory("backend/adapters/testsuite")
endif()

if(MARS_ENABLE_VTK)
  find_package(
    VTK
    COMPONENTS vtkCommonCore vtkCommonDataModel vtkFiltersGeneral vtkIOXML
               vtkIOParallel vtkIOParallelXML
    REQUIRED)
  add_subdirectory(backend/vtk)
endif()

# ##############################################################################
# ##############################################################################
# ##############################################################################

# #################  EXECUTABLES AND INSTALL ##################

add_executable(mars_exec main.cpp)
target_link_libraries(mars_exec mars)

# Adios2 Executables go here.

include_directories(${CMAKE_BINARY_DIR})

install(
  TARGETS mars mars_exec
  EXPORT MarsTargets
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/base/mars_config.hpp.in
               ${CMAKE_BINARY_DIR}/mars_config.hpp)
install(FILES ${MARS_HEADERS} DESTINATION include)
install(FILES ${CMAKE_BINARY_DIR}/mars_config.hpp DESTINATION include)

install(
    EXPORT MarsTargets
    FILE MarsTargets.cmake
    NAMESPACE Mars::
    DESTINATION lib/cmake)

# generate the config file that is includes the exports
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/MarsConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/MarsConfig.cmake"
  INSTALL_DESTINATION "lib/cmake/")

# generate the version file for the config file
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/MarsConfigVersion.cmake"
  VERSION "${MARS_VERSION}"
  COMPATIBILITY AnyNewerVersion)

# install the configuration file
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/MarsConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/MarsConfigVersion.cmake
        DESTINATION lib/cmake/)

include(cmake/PostInstall.cmake)

# ##############################################################################

# #################  Extra Targets, Complete ##################

add_custom_target(mars_complete)

if(MARS_ENABLE_TESTS)
  add_dependencies(mars_complete mars_test)
endif()


# ##############################################################################
target_format_setup("mars")
# ##############################################################################
