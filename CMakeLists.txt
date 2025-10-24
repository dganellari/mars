cmake_minimum_required(VERSION 3.22)

include(CheckLanguage)

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

enable_testing()
if(MARS_ENABLE_TESTS)
    include(cmake/Tests.cmake)
endif()


if(MARS_ENABLE_CUDA)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
        find_package(CUDAToolkit REQUIRED)
        set(CMAKE_CUDA_STANDARD 20)
        # 1. Add global-scope settings to ensure CUDA includes are visible everywhere
        include_directories(SYSTEM BEFORE ${CUDAToolkit_INCLUDE_DIRS})
        # 2. Add via compiler flags to ensure they're used in command-line compilation
        add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-isystem${CUDAToolkit_INCLUDE_DIRS}>)
        add_compile_options($<$<COMPILE_LANGUAGE:CUDA>:-isystem${CUDAToolkit_INCLUDE_DIRS}>)
        add_compile_definitions(USE_CUDA MARS_ENABLE_CUDA)
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --extended-lambda --expt-relaxed-constexpr")

        # Add environment variables for CUDA paths
        set(ENV{CPATH} "${CUDAToolkit_INCLUDE_DIRS}:$ENV{CPATH}")
        set(CUDA_TOOLKIT_ROOT_DIR "${CUDAToolkit_ROOT}" CACHE PATH "Path to CUDA Toolkit" FORCE)
        message(STATUS "Configured CUDA includes globally: ${CUDAToolkit_INCLUDE_DIRS}")
    else()
        message(STATUS "No CUDA support")
        set(MARS_ENABLE_CUDA OFF CACHE BOOL "Disable CUDA support" FORCE)
        set(CSTONE_WITH_CUDA OFF CACHE BOOL "Disable CUDA support" FORCE)
    endif()
endif()

if(MARS_ENABLE_HIP)
    check_language(HIP)
    if(CMAKE_HIP_COMPILER)
        enable_language(HIP)
        include(cmake/rocmlibs_target.cmake)
        set(CMAKE_HIP_STANDARD 20)
        add_compile_definitions(USE_CUDA MARS_ENABLE_HIP)

        include_directories(SYSTEM BEFORE ${HIP_INCLUDE_DIRS})
        add_compile_options($<$<AND:$<COMPILE_LANGUAGE:CXX>,$<NOT:$<STREQUAL:$<TARGET_PROPERTY:TYPE>,EXECUTABLE>>>:-isystem${HIP_INCLUDE_DIRS}>)
        add_compile_options($<$<COMPILE_LANGUAGE:HIP>:-isystem${HIP_INCLUDE_DIRS}>)

        set(CMAKE_HIP_SOURCE_FILE_EXTENSIONS hip;cu;hip.cpp)
        file(GLOB_RECURSE HIP_FILES

            "${CMAKE_SOURCE_DIR}/**/*.hip.cpp"
            "${CMAKE_SOURCE_DIR}/**/*.hip"
            "${CMAKE_SOURCE_DIR}/**/*.cu")
        foreach(HIP_FILE ${HIP_FILES})
            set_source_files_properties(${HIP_FILE}
                PROPERTIES LANGUAGE HIP)
        endforeach()
        # Set the platform flag for all HIP compilations
        add_compile_definitions(__HIP_PLATFORM_AMD__=1)
        # Set correct flags for hipcc
        set(HIP_STANDARD ${CMAKE_CXX_STANDARD})
        set(CMAKE_HIP_STANDARD_REQUIRED ${CMAKE_CXX_STANDARD_REQUIRED})
        message(STATUS "HIP compiler: ${CMAKE_HIP_COMPILER}")
        message(STATUS "HIP source extensions: ${CMAKE_HIP_SOURCE_FILE_EXTENSIONS}")

        set(ENV{CPATH} "${HIP_INCLUDE_DIRS}:$ENV{CPATH}")
        set(HIP_ROOT_DIR "${HIP_ROOT_DIR}" CACHE PATH "Path to HIP Root" FORCE)
        message(STATUS "Configured HIP includes globally: ${HIP_INCLUDE_DIRS}")

    else()
        message(STATUS "No HIP support")
        set(MARS_ENABLE_HIP OFF CACHE BOOL "Disable HIP support" FORCE)
        set(CSTONE_WITH_HIP OFF CACHE BOOL "Disable HIP support" FORCE)
    endif()
endif()

if(MARS_ENABLE_HIP AND MARS_ENABLE_CUDA)
    message(FATAL_ERROR "Cannot enable both HIP and CUDA at the same time.")
endif()

if(MARS_ENABLE_HIP)
    if(NOT DEFINED HIP_GPU_ARCHITECTURES)
        set(HIP_GPU_ARCHITECTURES "gfx942" CACHE STRING "AMD GPU targets (e.g., gfx90a, gfx942)")
    endif()
    set(CMAKE_HIP_ARCHITECTURES ${HIP_GPU_ARCHITECTURES})
    message(STATUS "HIP GPU architectures: ${CMAKE_HIP_ARCHITECTURES}")
    message(STATUS "HIP GPU architectures: ${HIP_GPU_ARCHITECTURES}")
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

# Find Cornerstone - centralize all Cornerstone finding and setup here
if(MARS_ENABLE_UNSTRUCTURED)
    set(CORNERSTONE_INSTALL_DIR "${CORNERSTONE_INSTALL_DIR}" CACHE PATH "Path to Cornerstone installation directory")
    find_package(cornerstone QUIET)
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

    if(NOT CORNERSTONE_INCLUDE_DIRS)
        message(STATUS "Cornerstone headers not found, fetching source...")
        
        include(FetchContent)
        
        # Set options BEFORE declaring the content
        if(MARS_ENABLE_CUDA)
            set(CSTONE_WITH_CUDA ON CACHE BOOL "Enable CUDA for Cornerstone" FORCE)
            set(CSTONE_WITH_GPU_AWARE_MPI ON CACHE BOOL "Enable GPU-aware MPI in Cornerstone" FORCE)
        elseif(MARS_ENABLE_HIP)
            set(CSTONE_WITH_HIP ON CACHE BOOL "Enable HIP for Cornerstone" FORCE)
            set(CSTONE_WITH_GPU_AWARE_MPI ON CACHE BOOL "Enable GPU-aware MPI in Cornerstone" FORCE)
        endif()
        
        # Disable ALL Cornerstone build targets except the library
        set(CSTONE_BUILD_TESTS OFF CACHE BOOL "Disable Cornerstone tests" FORCE)
        set(CSTONE_BUILD_EXAMPLES OFF CACHE BOOL "Disable Cornerstone examples" FORCE)
        set(BUILD_TESTING OFF CACHE BOOL "Disable global testing" FORCE)
        set(CSTONE_BUILD_UNIT_TESTS OFF CACHE BOOL "Disable Cornerstone unit tests" FORCE)
        set(CSTONE_BUILD_PERFORMANCE_TESTS OFF CACHE BOOL "Disable Cornerstone performance tests" FORCE)
        
        # Declare the content
        FetchContent_Declare(
            cornerstone_fetch
            GIT_REPOSITORY https://github.com/sekelle/cornerstone-octree.git
            GIT_TAG master
            OVERRIDE_FIND_PACKAGE
        )
    
        # Use the modern, non-deprecated approach
        FetchContent_MakeAvailable(cornerstone_fetch)
        
        # Get the source directory that was populated
        FetchContent_GetProperties(cornerstone_fetch SOURCE_DIR CORNERSTONE_SRC_DIR)
        
        # Set paths for includes
        set(CORNERSTONE_INCLUDE_DIRS "${CORNERSTONE_SRC_DIR}/include" CACHE PATH "Cornerstone include directory" FORCE)
        
        message(STATUS "Fetched Cornerstone from: ${CORNERSTONE_SRC_DIR}")
        
        # Now only build cstone_gpu if needed
        if(MARS_ENABLE_CUDA OR MARS_ENABLE_HIP)
            if(TARGET cstone_gpu)
                message(STATUS "cstone_gpu target will be built from source")
                set(CORNERSTONE_GPU_LIBRARY "$<TARGET_FILE:cstone_gpu>" CACHE PATH "Path to cstone_gpu library" FORCE)
            else()
                message(WARNING "cstone_gpu target not found in Cornerstone build")
            endif()
        endif()
    endif() 
    
    # Create CUDA/HIP support interface targets
    if(MARS_ENABLE_CUDA)
        if(NOT TARGET Mars::cuda_support)
            add_library(Mars_cuda_support INTERFACE)
            add_library(Mars::cuda_support ALIAS Mars_cuda_support)
            
            target_include_directories(Mars_cuda_support SYSTEM INTERFACE 
                ${CUDAToolkit_INCLUDE_DIRS})
            target_compile_definitions(Mars_cuda_support INTERFACE 
                USE_CUDA MARS_ENABLE_CUDA)
            
            target_link_libraries(Mars_cuda_support INTERFACE CUDA::cudart)
        endif()
        
        # Make cornerstone depend on CUDA support
        if(TARGET cornerstone::cornerstone)
            target_link_libraries(cornerstone::cornerstone INTERFACE Mars::cuda_support)
        endif()
        
        if(TARGET cstone_gpu)
            target_link_libraries(cstone_gpu INTERFACE Mars::cuda_support)
        endif()
    endif()
    
    if(MARS_ENABLE_HIP)
        if(NOT TARGET Mars::hip_support)
            add_library(Mars_hip_support INTERFACE)
            add_library(Mars::hip_support ALIAS Mars_hip_support)
            
            target_include_directories(Mars_hip_support SYSTEM INTERFACE 
                ${HIP_INCLUDE_DIRS})
            target_compile_definitions(Mars_hip_support INTERFACE 
                USE_CUDA MARS_ENABLE_HIP __HIP_PLATFORM_AMD__=1)
            
            if(NOT TARGET mars::rocmlibs)
                include(cmake/rocmlibs_target.cmake)
            endif()
            target_link_libraries(Mars_hip_support INTERFACE mars::rocmlibs)
        endif()
        
        # Make cornerstone depend on HIP support
        if(TARGET cornerstone::cornerstone)
            target_link_libraries(cornerstone::cornerstone INTERFACE Mars::hip_support)
        endif()
        
        if(TARGET cstone_gpu)
            target_link_libraries(cstone_gpu INTERFACE Mars::hip_support)
        endif()
    endif()
    
    # Setup the include directory structure for relative includes
    if(CORNERSTONE_INCLUDE_DIRS)
        # Create a directory for symlinks to make relative paths work
        file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/cornerstone_includes_temp)
        
        if(EXISTS "${CMAKE_BINARY_DIR}/cornerstone_includes")
            get_filename_component(REAL_PATH "${CMAKE_BINARY_DIR}/cornerstone_includes" REALPATH)
            if(NOT "${REAL_PATH}" STREQUAL "${CMAKE_BINARY_DIR}/cornerstone_includes")
                file(REMOVE "${CMAKE_BINARY_DIR}/cornerstone_includes")
                message(STATUS "Removed existing symlink at ${CMAKE_BINARY_DIR}/cornerstone_includes")
            else()
                file(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/cornerstone_includes")
                message(STATUS "Removed existing directory at ${CMAKE_BINARY_DIR}/cornerstone_includes")
            endif()
        endif()
        
        # Create a symbolic link to the Cornerstone include directory
        if(EXISTS "${CORNERSTONE_INCLUDE_DIRS}/cstone")
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
        if(MARS_ENABLE_CUDA OR MARS_ENABLE_HIP)
            if(NOT CORNERSTONE_GPU_LIBRARY OR NOT EXISTS "${CORNERSTONE_GPU_LIBRARY}")
                find_library(CORNERSTONE_GPU_LIBRARY
                    NAMES cstone_gpu libcstone_gpu
                    HINTS
                    ${CORNERSTONE_SRC_DIR}
                    ${CORNERSTONE_SRC_DIR}/build
                    ${CORNERSTONE_SRC_DIR}/lib
                    ${CORNERSTONE_SRC_DIR}/gpu
                    ${CORNERSTONE_SRC_DIR}/include/cstone
                    ${cornerstone_fetch_BINARY_DIR}/gpu
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
                    
                    if(MARS_ENABLE_CUDA)
                        set_property(TARGET cstone_gpu PROPERTY INTERFACE_LINK_LIBRARIES "Mars::cuda_support")
                    elseif(MARS_ENABLE_HIP)
                        set_property(TARGET cstone_gpu PROPERTY INTERFACE_LINK_LIBRARIES "Mars::hip_support")
                    endif()
                    message(STATUS "Created cstone_gpu target")
                endif()
            else()
                if(TARGET cstone_gpu)
                    message(STATUS "Using Cornerstone GPU library built from source")
                else()
                    message(WARNING "No Cornerstone GPU library found or built: ${CORNERSTONE_GPU_LIBRARY}")
                endif()
            endif()
        endif()
        
        set(CORNERSTONE_FOUND TRUE)
    else()
        set(CORNERSTONE_FOUND FALSE)
        message(WARNING "Cornerstone headers not found, falling back to ADIOS2")
    endif()
endif()

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

if(MARS_ENABLE_UNSTRUCTURED AND MARS_ENABLE_ADIOS2)
    # Add ADIOS2 subdirectories
    add_subdirectory(backend/adios2)
    if(MARS_ENABLE_TESTS)
        add_subdirectory("backend/adios2/testsuite")
    endif()
endif()

if(MARS_ENABLE_CXXOPTS)
  include(cxxopts/cxxopts.cmake)
endif()

add_library(mars ${MARS_HEADERS} ${MARS_SOURCES})

if(CORNERSTONE_FOUND)
    # Define a preprocessor macro to indicate cornerstone support
    target_compile_definitions(mars PUBLIC MARS_HAS_CORNERSTONE)
    # Link to the cornerstone target
    target_link_libraries(mars PUBLIC cornerstone::cornerstone)

    target_include_directories(mars PUBLIC
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/cornerstone_includes>
        $<BUILD_INTERFACE:${CORNERSTONE_INCLUDE_DIRS}>
        $<INSTALL_INTERFACE:include>) 
endif()

if(MARS_ENABLE_UNSTRUCTURED AND MARS_ENABLE_ADIOS2)
    # Make sure ADIOS2 targets are available in this scope
    find_package(ADIOS2 QUIET)
        
    # Link ADIOS2 dependencies if found
    if(ADIOS2_FOUND)
        target_include_directories(mars PUBLIC ${ADIOS2_INCLUDE_DIRS})
        target_link_libraries(mars PUBLIC ${ADIOS2_LIBRARIES})
    endif()
endif()

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

# Add HIP include paths directly to mars target for C++ compilation of HIP headers
if(MARS_ENABLE_HIP AND TARGET mars)
    target_include_directories(mars SYSTEM PRIVATE ${HIP_INCLUDE_DIRS})
    message(STATUS "Added HIP include directories to mars target")
endif()

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

# Final check for cstone_gpu and installation handling
if(TARGET cstone_gpu)
    message(STATUS "cstone_gpu target is available")
    get_target_property(GPU_LOCATION cstone_gpu IMPORTED_LOCATION)
    message(STATUS "cstone_gpu location: ${GPU_LOCATION}")

    # Check if it's an imported target
    get_target_property(IS_IMPORTED cstone_gpu IMPORTED)
    
    if(IS_IMPORTED)
        # For imported targets, install the library file directly
        if(EXISTS "${GPU_LOCATION}")
            message(STATUS "Installing GPU library directly")
            get_filename_component(GPU_LIB_NAME "${GPU_LOCATION}" NAME)
            install(FILES "${GPU_LOCATION}"
                    DESTINATION lib
                    RENAME "${GPU_LIB_NAME}")
        endif()
    else()
        # For locally built targets, add to export set
        list(APPEND MARS_TARGETS cstone_gpu)
        set_property(TARGET cstone_gpu PROPERTY EXPORT_NAME Mars::cstone_gpu)
    endif()
else()
    # Fallback for when target doesn't exist but we have a path
    if(EXISTS "${CORNERSTONE_GPU_LIBRARY}")
        message(STATUS "Installing GPU library from path")
        get_filename_component(GPU_LIB_NAME "${CORNERSTONE_GPU_LIBRARY}" NAME)
        install(FILES "${CORNERSTONE_GPU_LIBRARY}"
                DESTINATION lib
                RENAME "${GPU_LIB_NAME}")
    endif()
endif()

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
    get_property(MARS_TEST_TARGETS GLOBAL PROPERTY MARS_TEST_TARGETS)
    if(MARS_TEST_TARGETS)
        # Only add dependencies if there are actually targets
        add_dependencies(mars_complete ${MARS_TEST_TARGETS})
    endif()
endif()

target_format_setup("mars")
