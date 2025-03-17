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

if(MARS_ENABLE_CUDA)
    enable_language(CUDA)
endif()

# Add near the top, after enabling CUDA but before any add_subdirectory calls
if(CMAKE_CUDA_COMPILER)
    add_compile_options($<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda>)
    add_compile_options($<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>)
endif()

# Std CMake imports
include(CMakePackageConfigHelpers)
include(CMakeDependentOption)
# Cmake imports.
include(cmake/MarsDefaults.cmake)
include(cmake/MarsOptions.cmake)
include(cmake/MarsCMakeFunctions.cmake)
include(cmake/MarsDependencies.cmake)
include(cmake/MarsClangFormat.cmake)

# ######################### Library Configuration #############################

# set(MARS_MODULES base core distributed examples io kokkos)
set(MARS_MODULES base core examples io kokkos)

foreach(MODULE ${MARS_MODULES})
  add_subdirectory("${MODULE}")
endforeach(MODULE)

if(MARS_ENABLE_SERIAL_BACKEND)
    add_subdirectory("backend/serial")
endif()

# Now add distributed backend - dependencies handled in unstructured directory
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

if(CMAKE_CUDA_COMPILER AND TARGET mars)
    # Add CUDA flags to the mars target
    target_compile_options(mars PRIVATE 
        $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda>
        $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>)
endif()

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

# Configure cornerstone and ADIOS2 for main target if available from unstructured module
if(MARS_ENABLE_UNSTRUCTURED)
    # Define a preprocessor macro to indicate cornerstone support
    target_compile_definitions(mars PUBLIC MARS_HAS_CORNERSTONE)

     # Make sure ADIOS2 targets are available in this scope
    find_package(ADIOS2 QUIET)
    
    # Link dependencies if variables have been set
    if(ADIOS2_FOUND)
        target_include_directories(mars PUBLIC ${ADIOS2_INCLUDE_DIRS})
        target_link_libraries(mars PUBLIC ${ADIOS2_LIBRARIES})
    endif()
    
    # Link against cornerstone headers if available
    if(TARGET cornerstone_headers)
        target_link_libraries(mars PUBLIC cornerstone_headers)
    endif()
    
    # Add ADIOS2 subdirectories
    add_subdirectory(backend/adios2)
    if(MARS_ENABLE_TESTS)
        add_subdirectory("backend/adios2/testsuite")
    endif()
endif()

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