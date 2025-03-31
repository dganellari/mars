# cmake/Findcornerstone.cmake - Find module for cornerstone
# Supports both installed version and automatic Git fetching

# First check if there's an installation directory specified
if(CORNERSTONE_INSTALL_DIR)
    # Check for header files in source or installed location
    if(EXISTS "${CORNERSTONE_INSTALL_DIR}/include/cstone")
        message(STATUS "Found Cornerstone headers at ${CORNERSTONE_INSTALL_DIR}/include")
        set(cornerstone_INCLUDE_DIR "${CORNERSTONE_INSTALL_DIR}/include")
        set(CORNERSTONE_SRC_DIR "${CORNERSTONE_INSTALL_DIR}" CACHE PATH "Cornerstone source directory" FORCE)
        set(cornerstone_FOUND TRUE)
        # Look for GPU library if CUDA is enabled
        if(MARS_ENABLE_CUDA)
            # Allow user to specify GPU library location directly
            set(CORNERSTONE_GPU_LIBRARY_DIR "" CACHE PATH "Path to Cornerstone GPU library directory (optional)")
            if(CORNERSTONE_GPU_LIBRARY_DIR)
                # User specified a directory, search there first
                find_library(CORNERSTONE_GPU_LIBRARY
                    NAMES cornerstone_gpu cstone_gpu libcstone_gpu libcstone_gpu.a
                    PATHS "${CORNERSTONE_GPU_LIBRARY_DIR}"
                    NO_DEFAULT_PATH
                )
            endif()
            # If not found yet, try standard and non-standard locations
            if(NOT CORNERSTONE_GPU_LIBRARY)
                find_library(CORNERSTONE_GPU_LIBRARY
                    NAMES cornerstone_gpu cstone_gpu libcstone_gpu libcstone_gpu.a
                    PATHS "${CORNERSTONE_INSTALL_DIR}/lib64"
                          "${CORNERSTONE_INSTALL_DIR}/lib"
                          "${CORNERSTONE_INSTALL_DIR}/include/cstone"
                    NO_DEFAULT_PATH
                )
            endif()
            if(CORNERSTONE_GPU_LIBRARY)
                if(EXISTS "${CORNERSTONE_GPU_LIBRARY}")
                    message(STATUS "Found Cornerstone GPU library: ${CORNERSTONE_GPU_LIBRARY}")
                else()
                    message(STATUS "Cornerstone GPU library found but path is invalid: ${CORNERSTONE_GPU_LIBRARY}")
                    unset(CORNERSTONE_GPU_LIBRARY CACHE)
                endif()
            else()
                message(STATUS "Cornerstone GPU library not found. If needed, specify using -DCORNERSTONE_GPU_LIBRARY=/path/to/libcstone_gpu.a")
                # Allow manual specification
                set(CORNERSTONE_GPU_LIBRARY "CORNERSTONE_GPU_LIBRARY-NOTFOUND" CACHE FILEPATH "Path to Cornerstone GPU library (set manually if needed)")
            endif()
        endif()
    else()
        message(STATUS "Could not find Cornerstone headers at ${CORNERSTONE_INSTALL_DIR}/include")
    endif()
endif()

# If we get here, we need to find or fetch the source
find_path(cornerstone_INCLUDE_DIR
    NAMES cstone/octree.hpp
    PATHS ${CMAKE_INSTALL_PREFIX}/include
          $ENV{CORNERSTONE_DIR}/include
          /usr/include
          /usr/local/include
          ${CORNERSTONE_INSTALL_DIR}/include
    DOC "Cornerstone header directory"
)

# If found, determine source directory
if(cornerstone_INCLUDE_DIR)
    # Get the parent directory as source directory
    get_filename_component(CORNERSTONE_SRC_DIR "${cornerstone_INCLUDE_DIR}" DIRECTORY)
    set(CORNERSTONE_SRC_DIR "${CORNERSTONE_SRC_DIR}" CACHE PATH "Cornerstone source directory" FORCE)
    message(STATUS "Found Cornerstone headers: ${cornerstone_INCLUDE_DIR}")
    set(cornerstone_FOUND TRUE)
endif()

# If not found and Git fetch is enabled, get it from GitHub
if(NOT cornerstone_FOUND AND NOT CORNERSTONE_NO_GIT_FETCH)
    message(STATUS "Cornerstone headers not found, fetching from GitHub...")
    include(FetchContent)
    # Configure CMake options for fetched cornerstone
    set(BUILD_TESTING OFF CACHE BOOL "Disable cornerstone tests" FORCE)
    # Configure CUDA options if needed
    if(DEFINED MARS_ENABLE_CUDA AND MARS_ENABLE_CUDA)
        set(CSTONE_WITH_CUDA ON CACHE BOOL "Enable CUDA for cornerstone" FORCE)
        set(CSTONE_WITH_GPU_AWARE_MPI ON CACHE BOOL "Enable GPU-aware MPI" FORCE)
        if(DEFINED CMAKE_CUDA_ARCHITECTURES)
            set(CSTONE_CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES} CACHE STRING "CUDA architectures" FORCE)
        endif()
        if(DEFINED MPI_CXX_COMPILER)
            set(CSTONE_CUDA_FLAGS "-std=c++17 -ccbin=${MPI_CXX_COMPILER} --expt-relaxed-constexpr --extended-lambda"
                CACHE STRING "CUDA flags for Cornerstone" FORCE)
            set(CSTONE_EXTRA_CUDA_FLAGS "${CSTONE_CUDA_FLAGS}"
                CACHE STRING "Extra CUDA flags for Cornerstone" FORCE)
        endif()
    else()
        set(CSTONE_WITH_CUDA OFF CACHE BOOL "Disable CUDA for cornerstone" FORCE)
        set(CSTONE_WITH_GPU_AWARE_MPI OFF CACHE BOOL "Disable GPU-aware MPI" FORCE)
    endif()
    # Fetch cornerstone
    FetchContent_Declare(
        cornerstone_fetch
        GIT_REPOSITORY https://github.com/sekelle/cornerstone-octree
        GIT_TAG master
    )
    FetchContent_MakeAvailable(cornerstone_fetch)
    # Set directories
    set(cornerstone_INCLUDE_DIR "${cornerstone_fetch_SOURCE_DIR}/include")
    set(CORNERSTONE_SRC_DIR "${cornerstone_fetch_SOURCE_DIR}" CACHE PATH "Cornerstone source directory" FORCE)
    set(cornerstone_FOUND TRUE)
    message(STATUS "Fetched cornerstone from git: ${CORNERSTONE_SRC_DIR}")
endif()

# If CUDA is enabled, try to find GPU library in additional locations for fetched version
if(MARS_ENABLE_CUDA AND NOT CORNERSTONE_GPU_LIBRARY AND CORNERSTONE_SRC_DIR)
    find_library(CORNERSTONE_GPU_LIBRARY
        NAMES cstone_gpu libcstone_gpu
        PATHS "${CORNERSTONE_SRC_DIR}"
              "${CMAKE_BINARY_DIR}/cornerstone_build"
              "${CMAKE_BINARY_DIR}/_deps/cornerstone_fetch-build"
        NO_DEFAULT_PATH
    )
    if(CORNERSTONE_GPU_LIBRARY)
        message(STATUS "Found Cornerstone GPU library in build directory: ${CORNERSTONE_GPU_LIBRARY}")
    endif()
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
    if(EXISTS "${CORNERSTONE_GPU_LIBRARY}")
        add_library(cstone_gpu STATIC IMPORTED GLOBAL)
        set_target_properties(cstone_gpu PROPERTIES
            IMPORTED_LOCATION "${CORNERSTONE_GPU_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${cornerstone_INCLUDE_DIR}"
        )
        # Add CUDA dependencies
        if(MARS_ENABLE_CUDA)
            find_package(CUDAToolkit REQUIRED)
            set_property(TARGET cstone_gpu PROPERTY
                INTERFACE_LINK_LIBRARIES "CUDA::cudart"
            )
        endif()
        message(STATUS "Created cstone_gpu target with library: ${CORNERSTONE_GPU_LIBRARY}")
    else()
        message(WARNING "GPU library path '${CORNERSTONE_GPU_LIBRARY}' does not exist, cstone_gpu target not created")
        unset(CORNERSTONE_GPU_LIBRARY CACHE)
    endif()
endif()

mark_as_advanced(cornerstone_INCLUDE_DIR CORNERSTONE_GPU_LIBRARY)
