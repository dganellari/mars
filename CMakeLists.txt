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
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MARS_ENABLE_CUDA)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    include_directories(SYSTEM ${CUDAToolkit_INCLUDE_DIRS})
    message(STATUS "Added CUDA include directories globally: ${CUDAToolkit_INCLUDE_DIRS}")

    add_compile_definitions(USE_CUDA)
    add_compile_definitions(MARS_ENABLE_CUDA)

    set(CMAKE_CUDA_STANDARD 17)
    # Critical: Use the direct flag approach which overrides any CMake standard settings
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --extended-lambda --expt-relaxed-constexpr --Werror=cross-execution-space-call")
endif()

# Add near the top, after enabling CUDA but before any add_subdirectory calls
# if(CMAKE_CUDA_COMPILER)
#     add_compile_options($<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda>)
#     add_compile_options($<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>)
# endif()

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

# Add CUDA include paths directly to mars target for C++ compilation of CUDA headers
if(MARS_ENABLE_CUDA AND TARGET mars)
    target_include_directories(mars SYSTEM PRIVATE ${CUDAToolkit_INCLUDE_DIRS})
    message(STATUS "Added CUDA include directories to mars target")
endif()

target_include_directories(
  mars BEFORE
  PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
         $<INSTALL_INTERFACE:include>)

target_include_directories(mars PUBLIC ${MARS_DEP_INCLUDES})
target_link_libraries(mars PUBLIC ${MARS_DEP_LIBRARIES})

# Find Cornerstone - centralize all Cornerstone finding and setup here
if(MARS_ENABLE_UNSTRUCTURED)
    set(CORNERSTONE_INSTALL_DIR "${CORNERSTONE_INSTALL_DIR}" CACHE PATH "Path to Cornerstone installation directory")
    
    # First try find_package
    find_package(cornerstone QUIET)
    
    # Check if we have valid paths directly from command line
    if(CORNERSTONE_INSTALL_DIR)
        set(CORNERSTONE_SRC_DIR "${CORNERSTONE_INSTALL_DIR}" CACHE PATH "Cornerstone source directory" FORCE)
        set(CORNERSTONE_INCLUDE_DIRS "${CORNERSTONE_INSTALL_DIR}/include" CACHE PATH "Cornerstone include directory" FORCE)
        message(STATUS "Using Cornerstone from specified install dir: ${CORNERSTONE_INSTALL_DIR}")
    endif()
    
    # Fallback: search for headers in standard locations
    if(NOT CORNERSTONE_INCLUDE_DIRS)
        find_path(CORNERSTONE_INCLUDE_PATH cstone/domain/domain.hpp
            HINTS
            ${CMAKE_BINARY_DIR}/_deps/cornerstone_fetch-src/include
            ${cornerstone_INCLUDE_DIR}
            /usr/local/include
            /usr/include
            /opt/local/include
        )
        
        if(CORNERSTONE_INCLUDE_PATH)
            get_filename_component(CORNERSTONE_SRC_DIR "${CORNERSTONE_INCLUDE_PATH}" DIRECTORY)
            set(CORNERSTONE_SRC_DIR "${CORNERSTONE_SRC_DIR}" CACHE PATH "Cornerstone source directory" FORCE)
            set(CORNERSTONE_INCLUDE_DIRS "${CORNERSTONE_INCLUDE_PATH}" CACHE PATH "Cornerstone include directory" FORCE)
            message(STATUS "Found Cornerstone headers at: ${CORNERSTONE_INCLUDE_PATH}")
        endif()
    endif()
    
       # Check if Cornerstone is already found or needs to be fetched
    if(NOT CORNERSTONE_INCLUDE_DIRS)
        message(STATUS "Cornerstone headers not found at specified paths. Trying auto-fetch...")
        
        # Include FetchContent
        include(FetchContent)
        
        # Define the fetch
        FetchContent_Declare(
            cornerstone_fetch
            GIT_REPOSITORY https://github.com/exawind/cornerstone.git
            GIT_TAG master  # Or specific tag/commit
        )
        
        # Set build options for fetched Cornerstone before fetching
        set(CSTONE_ENABLE_CUDA ${MARS_ENABLE_CUDA} CACHE BOOL "Enable CUDA in fetched Cornerstone" FORCE)
        set(CSTONE_ENABLE_MPI ON CACHE BOOL "Enable MPI in fetched Cornerstone" FORCE)
        set(CSTONE_BUILD_EXAMPLES OFF CACHE BOOL "Don't build Cornerstone examples" FORCE)
        
        # Make it available - this will build it
        FetchContent_MakeAvailable(cornerstone_fetch)
        
        # Set paths based on fetched content
        set(CORNERSTONE_SRC_DIR ${cornerstone_fetch_SOURCE_DIR} CACHE PATH "Cornerstone source directory" FORCE)
        set(CORNERSTONE_INCLUDE_DIRS "${cornerstone_fetch_SOURCE_DIR}/include" CACHE PATH "Cornerstone include directory" FORCE)
        
        # If we're using CUDA, look for the GPU library in the build
        if(MARS_ENABLE_CUDA)
            set(CORNERSTONE_GPU_LIBRARY "${cornerstone_fetch_BINARY_DIR}/gpu/libcstone_gpu.a" 
                CACHE FILEPATH "Path to cstone_gpu library" FORCE)
            message(STATUS "Set CORNERSTONE_GPU_LIBRARY to ${CORNERSTONE_GPU_LIBRARY} (will be built)")
        endif()
        
        message(STATUS "Fetched Cornerstone from repository")
    endif()
    
    # Setup the include directory structure for relative includes
    if(CORNERSTONE_INCLUDE_DIRS)
        # Create a directory for symlinks to make relative paths work
        file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/cornerstone_includes_temp)
        
        # Check if the target path already exists and handle it appropriately
        if(EXISTS "${CMAKE_BINARY_DIR}/cornerstone_includes")
            # Check if it's a symbolic link
            get_filename_component(REAL_PATH "${CMAKE_BINARY_DIR}/cornerstone_includes" REALPATH)
            if(NOT "${REAL_PATH}" STREQUAL "${CMAKE_BINARY_DIR}/cornerstone_includes")
                # It's a symlink, remove it
                file(REMOVE "${CMAKE_BINARY_DIR}/cornerstone_includes")
                message(STATUS "Removed existing symlink at ${CMAKE_BINARY_DIR}/cornerstone_includes")
            else()
                # It's a real directory, remove it recursively
                file(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/cornerstone_includes")
                message(STATUS "Removed existing directory at ${CMAKE_BINARY_DIR}/cornerstone_includes")
            endif()
        endif()
        
        # Create a symbolic link to the Cornerstone include directory
        if(EXISTS "${CORNERSTONE_INCLUDE_DIRS}/cstone")
            # If the directory already has a "cstone" subdirectory, link to the parent
            execute_process(
                COMMAND ${CMAKE_COMMAND} -E create_symlink 
                        ${CORNERSTONE_INCLUDE_DIRS}
                        ${CMAKE_BINARY_DIR}/cornerstone_includes
                RESULT_VARIABLE SYMLINK_RESULT
            )
            if(NOT SYMLINK_RESULT EQUAL 0)
                message(WARNING "Failed to create symlink (error ${SYMLINK_RESULT}), using direct include path instead")
            endif()
        else()
            # Otherwise, create cstone subdirectory and link there
            file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/cornerstone_includes)
            execute_process(
                COMMAND ${CMAKE_COMMAND} -E create_symlink 
                        ${CORNERSTONE_INCLUDE_DIRS}
                        ${CMAKE_BINARY_DIR}/cornerstone_includes/cstone
                RESULT_VARIABLE SYMLINK_RESULT
            )
            if(NOT SYMLINK_RESULT EQUAL 0)
                message(WARNING "Failed to create symlink (error ${SYMLINK_RESULT}), using direct include path instead")
            endif()
        endif()
        
        # Make this a global include path with high priority
        include_directories(BEFORE ${CMAKE_BINARY_DIR}/cornerstone_includes)
        message(STATUS "Created Cornerstone include structure at: ${CMAKE_BINARY_DIR}/cornerstone_includes")
        
        # Create the interface target if it doesn't exist
        if(NOT TARGET cornerstone::cornerstone)
            add_library(cornerstone::cornerstone INTERFACE IMPORTED)
            set_target_properties(cornerstone::cornerstone PROPERTIES
                INTERFACE_INCLUDE_DIRECTORIES "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/cornerstone_includes>;$<BUILD_INTERFACE:${CORNERSTONE_INCLUDE_DIRS}>")
            message(STATUS "Created cornerstone::cornerstone target with include directories")
        endif()
        
        # Handle GPU library
        if(MARS_ENABLE_CUDA)
            # If we fetched Cornerstone, the GPU library path is already set
            # Otherwise, search for it
            if(NOT CORNERSTONE_GPU_LIBRARY OR NOT EXISTS "${CORNERSTONE_GPU_LIBRARY}")
                find_library(CORNERSTONE_GPU_LIBRARY
                    NAMES cstone_gpu libcstone_gpu
                    HINTS
                    ${CORNERSTONE_SRC_DIR}
                    ${CORNERSTONE_SRC_DIR}/build
                    ${CORNERSTONE_SRC_DIR}/lib
                    ${CORNERSTONE_SRC_DIR}/gpu
                    ${CORNERSTONE_SRC_DIR}/gpu/include/cstone
                    ${cornerstone_fetch_BINARY_DIR}/gpu
                    ${CMAKE_BINARY_DIR}/cornerstone_build
                    ${CMAKE_BINARY_DIR}/_deps/cornerstone_fetch-build
                    ${CMAKE_BINARY_DIR}/_deps/cornerstone_fetch-build/gpu
                    NO_DEFAULT_PATH
                )
            endif()
            
            if(CORNERSTONE_GPU_LIBRARY AND EXISTS "${CORNERSTONE_GPU_LIBRARY}")
                message(STATUS "Found Cornerstone GPU library: ${CORNERSTONE_GPU_LIBRARY}")
                if(NOT TARGET cstone_gpu)
                    add_library(cstone_gpu STATIC IMPORTED GLOBAL)
                    set_target_properties(cstone_gpu PROPERTIES
                        IMPORTED_LOCATION "${CORNERSTONE_GPU_LIBRARY}"
                        INTERFACE_INCLUDE_DIRECTORIES "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/cornerstone_includes>;$<BUILD_INTERFACE:${CORNERSTONE_INCLUDE_DIRS}>")
                    
                    find_package(CUDAToolkit REQUIRED)
                    set_property(TARGET cstone_gpu PROPERTY
                        INTERFACE_LINK_LIBRARIES "CUDA::cudart")
                    
                    message(STATUS "Created cstone_gpu target")
                endif()
            else()
                message(WARNING "Cornerstone GPU library not found: ${CORNERSTONE_GPU_LIBRARY}")
            endif()
        endif()
                    
        # Mark that we found Cornerstone
        set(CORNERSTONE_FOUND TRUE)
        
        # Define a preprocessor macro to indicate cornerstone support
        target_compile_definitions(mars PUBLIC MARS_HAS_CORNERSTONE)
        
        # Link to the cornerstone target
        target_link_libraries(mars PUBLIC cornerstone::cornerstone)
        
        # Also directly add the include directories for absolute certainty
        target_include_directories(mars PUBLIC 
            $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/cornerstone_includes>
            $<BUILD_INTERFACE:${CORNERSTONE_INCLUDE_DIRS}>
            $<INSTALL_INTERFACE:include>) 
    else()
        set(CORNERSTONE_FOUND FALSE)
        message(WARNING "Cornerstone headers not found, falling back to ADIOS2")
    endif() 
    
    # Make sure ADIOS2 targets are available in this scope
    find_package(ADIOS2 QUIET)
    
    # Link ADIOS2 dependencies if found
    if(ADIOS2_FOUND)
        target_include_directories(mars PUBLIC ${ADIOS2_INCLUDE_DIRS})
        target_link_libraries(mars PUBLIC ${ADIOS2_LIBRARIES})
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

# Final check to see if we created the target
if(TARGET cstone_gpu)
    message(STATUS "cstone_gpu target is available for installation")
    get_target_property(GPU_LOCATION cstone_gpu IMPORTED_LOCATION)
    message(STATUS "cstone_gpu location: ${GPU_LOCATION}")
    # Add to MARS_TARGETS only if it exists in this scope
    list(APPEND MARS_TARGETS cstone_gpu)
    # Install the target (avoid installation for imported targets)
    if(EXISTS "${GPU_LOCATION}")
        message(STATUS "Installing GPU library directly: ${GPU_LOCATION}")
        get_filename_component(GPU_LIB_NAME "${GPU_LOCATION}" NAME)
        install(FILES "${GPU_LOCATION}"
                DESTINATION ${CMAKE_INSTALL_LIBDIR}
                RENAME "${GPU_LIB_NAME}")
    endif()
else()
    message(STATUS "cstone_gpu target not created, skipping installation")
    # Fallback: try using CORNERSTONE_GPU_LIBRARY directly
    if(EXISTS "${CORNERSTONE_GPU_LIBRARY}")
        message(STATUS "Installing GPU library directly from specified path: ${CORNERSTONE_GPU_LIBRARY}")
        get_filename_component(GPU_LIB_NAME "${CORNERSTONE_GPU_LIBRARY}" NAME)
        install(FILES "${CORNERSTONE_GPU_LIBRARY}"
                DESTINATION ${CMAKE_INSTALL_LIBDIR}
                RENAME "${GPU_LIB_NAME}")
    endif()
endif()

# Install the standard targets, but remove cstone_gpu from MARS_TARGETS if it exists
if(TARGET cstone_gpu)
    list(REMOVE_ITEM MARS_TARGETS cstone_gpu)
endif()

# Install the standard targets
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

# Fix for the add_dependencies error
if(MARS_ENABLE_TESTS)
    get_property(MARS_TEST_TARGETS GLOBAL PROPERTY MARS_TEST_TARGETS)
    if(MARS_TEST_TARGETS)
        # Only add dependencies if there are actually targets
        add_dependencies(mars_complete ${MARS_TEST_TARGETS})
    endif()
endif()

target_format_setup("mars")
