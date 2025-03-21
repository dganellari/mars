# cmake/Findcornerstone.cmake - Find module for cornerstone
# Supports both installed version and automatic Git fetching

# First try to find the installed version
find_path(cornerstone_INCLUDE_DIR
    NAMES octree.h cornerstone/octree.h cornerstone-octree.h
    PATHS ${CORNERSTONE_ROOT_DIR}
          ${CORNERSTONE_ROOT_DIR}/include
          ${CMAKE_INSTALL_PREFIX}
          ${CMAKE_INSTALL_PREFIX}/include
          $ENV{CORNERSTONE_DIR}
          $ENV{CORNERSTONE_DIR}/include
          /usr/include
          /usr/local/include
    PATH_SUFFIXES cornerstone cornerstone-octree
    DOC "Cornerstone header directory"
)

# Add debug output
message(STATUS "Looking for cornerstone headers...")
message(STATUS "  CORNERSTONE_ROOT_DIR = ${CORNERSTONE_ROOT_DIR}")
message(STATUS "  cornerstone_INCLUDE_DIR = ${cornerstone_INCLUDE_DIR}")

# If not found and Git fetch is enabled, get it from GitHub
if(NOT cornerstone_INCLUDE_DIR AND NOT CORNERSTONE_NO_GIT_FETCH)
    message(STATUS "Cornerstone headers not found, fetching from GitHub...")
    # Include FetchContent module (available in CMake 3.11+)
    include(FetchContent)
    # Disable testing for cornerstone
    set(BUILD_TESTING OFF CACHE BOOL "Disable cornerstone tests" FORCE)
    # Pass through CUDA/GPU options if they exist in the parent project
    if(DEFINED CSTONE_WITH_CUDA)
        # User has already defined this, so we'll respect it
        message(STATUS "Using parent project's CSTONE_WITH_CUDA=${CSTONE_WITH_CUDA}")
    else()
        # Set default based on whether MARS uses CUDA
        if(DEFINED MARS_ENABLE_CUDA AND MARS_ENABLE_CUDA)
            set(CSTONE_WITH_CUDA ON CACHE BOOL "Enable CUDA for cornerstone" FORCE)
            set(CSTONE_WITH_GPU_AWARE_MPI ON CACHE BOOL "Enable GPU-aware MPI" FORCE)
        else()
            set(CSTONE_WITH_CUDA OFF CACHE BOOL "Disable CUDA for cornerstone" FORCE)
            set(CSTONE_WITH_GPU_AWARE_MPI OFF CACHE BOOL "Disable GPU-aware MPI" FORCE)
        endif()
    endif()
    # Declare cornerstone dependency
    FetchContent_Declare(
        cornerstone_fetch
        GIT_REPOSITORY https://github.com/sekelle/cornerstone-octree
        GIT_TAG master # Change to a specific tag/commit if needed
    )
    # Use MakeAvailable to avoid the warning
    FetchContent_MakeAvailable(cornerstone_fetch)
    # Set the include directory to the fetched source
    set(cornerstone_INCLUDE_DIR "${cornerstone_fetch_SOURCE_DIR}")
    message(STATUS "Fetched cornerstone from git: ${cornerstone_INCLUDE_DIR}")
endif()

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(cornerstone DEFAULT_MSG
    cornerstone_INCLUDE_DIR
)

# Create an INTERFACE library for header-only usage
if(cornerstone_FOUND AND NOT TARGET cornerstone::cornerstone)
    add_library(cornerstone::cornerstone INTERFACE IMPORTED)
    set_target_properties(cornerstone::cornerstone PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${cornerstone_INCLUDE_DIR}"
    )
endif()

mark_as_advanced(cornerstone_INCLUDE_DIR)
