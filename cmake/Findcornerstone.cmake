# cmake/Findcornerstone.cmake - Find module for cornerstone
# Supports both installed version and automatic Git fetching

# First check if there's an installation directory specified
if(CORNERSTONE_INSTALL_DIR)
    # Check for header files in installed location
    if(EXISTS "${CORNERSTONE_INSTALL_DIR}/include/cstone")
        message(STATUS "Found installed Cornerstone headers at ${CORNERSTONE_INSTALL_DIR}/include")
        set(cornerstone_INCLUDE_DIR "${CORNERSTONE_INSTALL_DIR}/include")
        set(cornerstone_FOUND TRUE)
        
    # Check if this is a build directory (not installed)
    elseif(EXISTS "${CORNERSTONE_INSTALL_DIR}/CMakeCache.txt")
        message(STATUS "Cornerstone build directory detected")
        
        # For build directories, the headers are usually in the source dir
        if(EXISTS "${CORNERSTONE_INSTALL_DIR}/../include/cstone" OR 
           EXISTS "${CORNERSTONE_INSTALL_DIR}/src/include/cstone")
            # Try first possibility - source dir is parent of build dir
            if(EXISTS "${CORNERSTONE_INSTALL_DIR}/../include/cstone")
                set(cornerstone_INCLUDE_DIR "${CORNERSTONE_INSTALL_DIR}/../include")
                message(STATUS "Found Cornerstone headers at ${cornerstone_INCLUDE_DIR}")
            # Try second possibility - source is in src subdir of build dir
            else()
                set(cornerstone_INCLUDE_DIR "${CORNERSTONE_INSTALL_DIR}/src/include")
                message(STATUS "Found Cornerstone headers at ${cornerstone_INCLUDE_DIR}")
            endif()
            
            set(cornerstone_FOUND TRUE)
            
            # Look for any built library in the build directory
            if(MARS_ENABLE_CUDA)
                find_library(CORNERSTONE_GPU_LIBRARY
                    NAMES cornerstone_gpu cstone_gpu libcstone_gpu
                    PATHS "${CORNERSTONE_INSTALL_DIR}/lib64"
                          "${CORNERSTONE_INSTALL_DIR}/lib"
                          "${CORNERSTONE_INSTALL_DIR}"
                          "${CORNERSTONE_INSTALL_DIR}/cstone"
                    NO_DEFAULT_PATH
                )
                set(CORNERSTONE_SRC_DIR "${CORNERSTONE_INSTALL_DIR}" CACHE PATH "Cornerstone source directory" FORCE)
                
                if(CORNERSTONE_GPU_LIBRARY)
                    message(STATUS "Found Cornerstone GPU library: ${CORNERSTONE_GPU_LIBRARY}")
                else()
                    message(STATUS "Cornerstone GPU library not found in build dir, GPU support may be limited")
                endif()
            endif()
            
            # We found headers in the build directory, return early
            return()
        else()
            message(STATUS "Cornerstone build directory detected, but couldn't find headers")
        endif()
    else()
        message(STATUS "Cornerstone installation directory specified but headers not found at ${CORNERSTONE_INSTALL_DIR}/include")
    endif()
    
    # Try one more pattern - maybe it's a source directory?
    if(EXISTS "${CORNERSTONE_INSTALL_DIR}/include/cstone/octree.hpp")
        message(STATUS "Found Cornerstone source directory")
        set(cornerstone_INCLUDE_DIR "${CORNERSTONE_INSTALL_DIR}/include")
        set(CORNERSTONE_SRC_DIR "${CORNERSTONE_INSTALL_DIR}" CACHE PATH "Cornerstone source directory" FORCE)
        set(cornerstone_FOUND TRUE)
        return()
    endif()
endif()

# If we get here, we need to find or fetch the source
find_path(cornerstone_INCLUDE_DIR
    NAMES cstone/octree.hpp cstone/domain/domain.hpp
    PATHS ${CMAKE_INSTALL_PREFIX}/include
          $ENV{CORNERSTONE_DIR}/include
          /usr/include
          /usr/local/include
    DOC "Cornerstone header directory"
)

# If not found and Git fetch is enabled, get it from GitHub
if(NOT cornerstone_INCLUDE_DIR AND NOT CORNERSTONE_NO_GIT_FETCH)
    message(STATUS "Cornerstone headers not found, fetching from GitHub...")
    include(FetchContent)
    
    # Disable testing for cornerstone
    set(BUILD_TESTING OFF CACHE BOOL "Disable cornerstone tests" FORCE)
    
    # Configure CUDA options for fetched cornerstone
    if(DEFINED MARS_ENABLE_CUDA AND MARS_ENABLE_CUDA)
        set(CSTONE_WITH_CUDA ON CACHE BOOL "Enable CUDA for cornerstone" FORCE)
        set(CSTONE_WITH_GPU_AWARE_MPI ON CACHE BOOL "Enable GPU-aware MPI" FORCE)
        if(DEFINED CMAKE_CUDA_ARCHITECTURES)
            set(CSTONE_CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES} CACHE STRING "CUDA architectures" FORCE)
        endif()
        
        # Set CUDA flags
        if(DEFINED MPI_CXX_COMPILER)
            set(CSTONE_CUDA_FLAGS "-std=c++17 -ccbin=${MPI_CXX_COMPILER} --expt-relaxed-constexpr --extended-lambda" 
                CACHE STRING "CUDA flags for Cornerstone" FORCE)
            set(CSTONE_EXTRA_CUDA_FLAGS "-std=c++17 -ccbin=${MPI_CXX_COMPILER} --expt-relaxed-constexpr --extended-lambda" 
                CACHE STRING "Extra CUDA flags for Cornerstone" FORCE)
        endif()
    else()
        set(CSTONE_WITH_CUDA OFF CACHE BOOL "Disable CUDA for cornerstone" FORCE)
        set(CSTONE_WITH_GPU_AWARE_MPI OFF CACHE BOOL "Disable GPU-aware MPI" FORCE)
    endif()
    
    # Declare cornerstone dependency
    FetchContent_Declare(
        cornerstone_fetch
        GIT_REPOSITORY https://github.com/sekelle/cornerstone-octree
        GIT_TAG master
    )
    
    # Make the content available
    FetchContent_MakeAvailable(cornerstone_fetch)
    
    # Set the include directory and source directory
    set(cornerstone_INCLUDE_DIR "${cornerstone_fetch_SOURCE_DIR}/include")
    set(CORNERSTONE_SRC_DIR "${cornerstone_fetch_SOURCE_DIR}" CACHE PATH "Cornerstone source directory" FORCE)
    message(STATUS "Fetched cornerstone from git: ${CORNERSTONE_SRC_DIR}")
endif()

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(cornerstone DEFAULT_MSG
    cornerstone_INCLUDE_DIR
)

# Create an interface library for header-only usage
if(cornerstone_FOUND AND NOT TARGET cornerstone::cornerstone)
    add_library(cornerstone::cornerstone INTERFACE IMPORTED)
    set_target_properties(cornerstone::cornerstone PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${cornerstone_INCLUDE_DIR}"
    )
endif()

# Create GPU target if library was found
if(CORNERSTONE_GPU_LIBRARY AND NOT TARGET cstone_gpu)
    add_library(cstone_gpu UNKNOWN IMPORTED)
    set_target_properties(cstone_gpu PROPERTIES
        IMPORTED_LOCATION "${CORNERSTONE_GPU_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${cornerstone_INCLUDE_DIR}"
    )
endif()

mark_as_advanced(cornerstone_INCLUDE_DIR CORNERSTONE_GPU_LIBRARY)

# Export variables for patch application
if(TARGET cornerstone_fetch)
    set(CORNERSTONE_SRC_DIR "${cornerstone_fetch_SOURCE_DIR}" CACHE PATH "Cornerstone source directory" FORCE)
endif()