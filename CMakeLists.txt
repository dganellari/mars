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
  HOMEPAGE_URL "https://github.com/dganellari/mars"
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

# ######################### Library Configuration #############################

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

if(MARS_ENABLE_MOONOLITH)
  add_subdirectory(moonolith_adapter)
endif()

if(MARS_ENABLE_CXXOPTS)
  include(cxxopts/cxxopts.cmake)
endif()

add_library(mars ${MARS_HEADERS} ${MARS_SOURCES})

foreach(include ${MARS_BUILD_INCLUDES})
  target_include_directories(mars BEFORE PUBLIC $<BUILD_INTERFACE:${include}>)
endforeach()

target_include_directories(
  mars BEFORE
  PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
         $<INSTALL_INTERFACE:include>)

target_include_directories(mars PUBLIC ${MARS_DEP_INCLUDES})
target_link_libraries(mars PUBLIC ${MARS_DEP_LIBRARIES})

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel." FORCE)

  message(STATUS "[Status] CMAKE_BUILD_TYPE=Release")
endif(NOT CMAKE_BUILD_TYPE)

# ###############################Extra Targets##################################

if(MARS_ENABLE_TESTS)
  include(cmake/Tests.cmake)
  add_subdirectory("backend/distributed/testsuite")
endif()

if(MARS_ENABLE_VTK)
  find_package(
    VTK
    COMPONENTS vtkCommonCore vtkCommonDataModel vtkFiltersGeneral vtkIOXML
               vtkIOParallel vtkIOParallelXML
    REQUIRED)
  add_subdirectory(backend/vtk)
endif()

if(MARS_ENABLE_UNSTRUCTURED)
    find_package(ADIOS2 REQUIRED)

    if(ADIOS2_FOUND)
        message(STATUS "ADIOS2 found: ${ADIOS2_DIR}")
        # Check for CUDA support in ADIOS2 if MARS_ENABLE_CUDA is enabled
        if(MARS_ENABLE_CUDA)
            if(NOT ADIOS2_HAVE_CUDA)
                message(FATAL_ERROR "ADIOS2 was found but it does not have CUDA support. Please install ADIOS2 with CUDA support or disable MARS_ENABLE_CUDA.")
            endif()
        endif()
        add_subdirectory(backend/adios2)
        if(MARS_ENABLE_TESTS)
            add_subdirectory("backend/adios2/testsuite")
        endif()

        # Add the ADIOS2 target to the list of targets to install
        target_include_directories(mars PUBLIC ${ADIOS2_INCLUDE_DIRS})
        target_link_libraries(mars PUBLIC ${ADIOS2_LIBRARIES})
    else()
        message(FATAL_ERROR "ADIOS2 not found. Please install ADIOS2 or disable MARS_ENABLE_UNSTRUCTURED.")
    endif()

    # Now handle Cornerstone dependency for unstructured support
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
    find_package(cornerstone REQUIRED)

    if(cornerstone_FOUND)
        # Add explicit include directory for cornerstone
        message(STATUS "Cornerstone include directory: ${cornerstone_INCLUDE_DIR}")
        # Get the actual source directory from FetchContent
        set(CORNERSTONE_SRC_DIR "${CMAKE_BINARY_DIR}/_deps/cornerstone_fetch-src")
        # The proper include path should be the 'include' directory inside the source
        if(EXISTS "${CORNERSTONE_SRC_DIR}/include/cstone/domain/domain.hpp")
            message(STATUS "Found domain.hpp in expected location")

            # Create an interface library for cornerstone
            add_library(cornerstone_headers INTERFACE)
            target_include_directories(cornerstone_headers INTERFACE
                $<BUILD_INTERFACE:${CORNERSTONE_SRC_DIR}/include>
                $<INSTALL_INTERFACE:include>
            )
            # Link against this interface library
            target_link_libraries(mars PUBLIC cornerstone_headers)

            # Add cornerstone_headers to the list of targets to install
            list(APPEND MARS_TARGETS cornerstone_headers)
        else()
            message(STATUS "Could not find domain.hpp in expected location")
            # Try alternate approach - copy headers to build directory include path
            file(COPY ${CORNERSTONE_SRC_DIR}/include/ DESTINATION ${CMAKE_BINARY_DIR}/include/)
            message(STATUS "Copied cornerstone headers to ${CMAKE_BINARY_DIR}/include/")
            # This path is already in the include directories
        endif()

        # Define a preprocessor macro to indicate cornerstone support
        target_compile_definitions(mars PUBLIC MARS_HAS_CORNERSTONE)
    else()
        message(FATAL_ERROR "Cornerstone not found. Unstructured mesh support will be limited.")
    endif()
endif()

# #################  EXECUTABLES AND INSTALL ##################

include_directories(${CMAKE_BINARY_DIR})

install(
  TARGETS mars ${MARS_TARGETS}
  EXPORT MarsTargets
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION include)

# Export the cxxopts target if it was fetched
if(TARGET cxxopts)
    install(TARGETS cxxopts
        EXPORT MarsTargets
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        INCLUDES DESTINATION include
    )
endif()

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

# #################  Extra Targets, Complete ##################

add_custom_target(mars_complete)

if(MARS_ENABLE_TESTS)
    # Get the list of tests from the global MARS_TESTS property
    get_property(MARS_TESTS GLOBAL PROPERTY MARS_TESTS)

    # Add the tests as dependencies
    add_dependencies(mars_complete ${MARS_TESTS})
endif()

target_format_setup("mars")
