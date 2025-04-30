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


# Modify your CUDA setup section (around line 21):
if(MARS_ENABLE_CUDA)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
        find_package(CUDAToolkit REQUIRED)
        set(CMAKE_CUDA_STANDARD 17)
        
        # 1. Add global-scope settings to ensure CUDA includes are visible everywhere
        include_directories(SYSTEM BEFORE ${CUDAToolkit_INCLUDE_DIRS})
        
        # 2. Add via compiler flags to ensure they're used in command-line compilation
        add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-isystem${CUDAToolkit_INCLUDE_DIRS}>)
        add_compile_options($<$<COMPILE_LANGUAGE:CUDA>:-isystem${CUDAToolkit_INCLUDE_DIRS}>)
        
        # 3. Define macros but only via CMake (not in source files)
        add_compile_definitions(USE_CUDA MARS_ENABLE_CUDA)
        
        # 4. Set CUDA flags (unchanged)
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
        set(CMAKE_HIP_STANDARD 17)
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

# Near the beginning of your CMakeLists.txt, add this if it doesn't exist:
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

        # In the Cornerstone fetch section of your main CMakeLists.txt
    if(NOT CORNERSTONE_INCLUDE_DIRS)
        message(STATUS "Cornerstone headers not found, fetching source...")

        # Include FetchContent
        include(FetchContent)
        
        # Define and populate
        FetchContent_Declare(
            cornerstone_fetch
            GIT_REPOSITORY https://github.com/sekelle/cornerstone-octree.git
            GIT_TAG master  # Or specific tag/commit
        )

        # Set options for Cornerstone before making available
        if(MARS_ENABLE_CUDA)
            set(CSTONE_WITH_CUDA ON CACHE BOOL "Enable CUDA for Cornerstone" FORCE)
            set(CSTONE_WITH_GPU_AWARE_MPI ON CACHE BOOL "Enable GPU-aware MPI in Cornerstone" FORCE)
        elseif(MARS_ENABLE_HIP)
            set(CSTONE_WITH_HIP ON CACHE BOOL "Enable HIP for Cornerstone" FORCE)
            set(CSTONE_WITH_GPU_AWARE_MPI ON CACHE BOOL "Enable GPU-aware MPI in Cornerstone" FORCE)  
        endif()        

        # Any other Cornerstone options you want to control
        set(CSTONE_BUILD_TESTS OFF CACHE BOOL "Disable Cornerstone tests" FORCE)
        set(CSTONE_BUILD_EXAMPLES OFF CACHE BOOL "Disable Cornerstone examples" FORCE)
 
        # Make the content available to use
        FetchContent_MakeAvailable(cornerstone_fetch)
        
        # Set paths for includes - these will be used later
        set(CORNERSTONE_SRC_DIR "${cornerstone_fetch_SOURCE_DIR}" CACHE PATH "Cornerstone source directory" FORCE)
        set(CORNERSTONE_INCLUDE_DIRS "${cornerstone_fetch_SOURCE_DIR}/include" CACHE PATH "Cornerstone include directory" FORCE)
        
        # Add the subdirectory to build it using its own CMakeLists.txt
        add_subdirectory(${cornerstone_fetch_SOURCE_DIR})
        
        # Now cstone_gpu target should be available if GPU was enabled
        if(MARS_ENABLE_CUDA OR MARS_ENABLE_HIP)
            if(TARGET cstone_gpu)
                message(STATUS "Using cstone_gpu target from Cornerstone's own build system")
                set(CORNERSTONE_GPU_LIBRARY "$<TARGET_FILE:cstone_gpu>" CACHE PATH "Path to cstone_gpu library" FORCE)
            endif()
        endif()
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
                if(TARGET cstone_gpu)
                    message(STATUS "Using Cornerstone GPU library built from source")
                else()
                    message(WARNING "No Cornerstone GPU library found or built: ${CORNERSTONE_GPU_LIBRARY}")
                endif()
            endif()
        elseif(MARS_ENABLE_HIP)
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
                message(STATUS "Found Cornerstone GPU library for HIP: ${CORNERSTONE_GPU_LIBRARY}")
                if(NOT TARGET cstone_gpu)
                    add_library(cstone_gpu STATIC IMPORTED GLOBAL)
                    set_target_properties(cstone_gpu PROPERTIES
                        IMPORTED_LOCATION "${CORNERSTONE_GPU_LIBRARY}"
                        INTERFACE_INCLUDE_DIRECTORIES "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/cornerstone_includes>;$<BUILD_INTERFACE:${CORNERSTONE_INCLUDE_DIRS}>")
                    
                    if(NOT TARGET mars::rocmlibs)
                        include(cmake/rocmlibs_target.cmake)
                    endif()

                    set_property(TARGET cstone_gpu PROPERTY
                        INTERFACE_LINK_LIBRARIES "mars::rocmlibs")
                    
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
                    
        # Mark that we found Cornerstone
        set(CORNERSTONE_FOUND TRUE)
        
        if(MARS_ENABLE_CUDA)
            # Make sure CUDA runtime headers are available to cornerstone
            if(TARGET cornerstone::cornerstone)
                target_include_directories(cornerstone::cornerstone SYSTEM INTERFACE ${CUDAToolkit_INCLUDE_DIRS})
                message(STATUS "Added CUDA include directories to cornerstone target: ${CUDAToolkit_INCLUDE_DIRS}")
            endif()
            
            # Also add to cstone_gpu if it exists
            if(TARGET cstone_gpu)
                target_include_directories(cstone_gpu SYSTEM INTERFACE ${CUDAToolkit_INCLUDE_DIRS})
            endif()
        endif()

        if(MARS_ENABLE_HIP)
            # Make sure HIP runtime headers are available to cornerstone
            if(TARGET cornerstone::cornerstone)
                target_include_directories(cornerstone::cornerstone SYSTEM INTERFACE ${HIP_INCLUDE_DIRS})
                message(STATUS "Added HIP include directories to cornerstone target: ${HIP_INCLUDE_DIRS}")
            endif()
            
            # Also add to cstone_gpu if it exists
            if(TARGET cstone_gpu)
                target_include_directories(cstone_gpu SYSTEM INTERFACE ${HIP_INCLUDE_DIRS})
            endif()
        endif()
    else()
        set(CORNERSTONE_FOUND FALSE)
        message(WARNING "Cornerstone headers not found, falling back to ADIOS2")
    endif()     
endif()

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

# Move these lines from the Cornerstone setup to after add_library(mars)
if(CORNERSTONE_FOUND)
    # Define a preprocessor macro to indicate cornerstone support
    target_compile_definitions(mars PUBLIC MARS_HAS_CORNERSTONE)
    
    # Link to the cornerstone target
    target_link_libraries(mars PUBLIC cornerstone::cornerstone)
    
    # Also directly add the include directories for absolute certainty
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

# Fix for the add_dependencies error
if(MARS_ENABLE_TESTS)
    get_property(MARS_TEST_TARGETS GLOBAL PROPERTY MARS_TEST_TARGETS)
    if(MARS_TEST_TARGETS)
        # Only add dependencies if there are actually targets
        add_dependencies(mars_complete ${MARS_TEST_TARGETS})
    endif()
endif()

target_format_setup("mars")
