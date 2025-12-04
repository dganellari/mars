# cmake/Findcornerstone.cmake - Find module for Cornerstone
# Supports both installed and fetched versions of Cornerstone

# Check if Cornerstone is installed
if(CORNERSTONE_INSTALL_DIR)
    if(EXISTS "${CORNERSTONE_INSTALL_DIR}/include/cstone")
        message(STATUS "Found Cornerstone headers at ${CORNERSTONE_INSTALL_DIR}/include")
        set(cornerstone_INCLUDE_DIR "${CORNERSTONE_INSTALL_DIR}/include")
        set(CORNERSTONE_SRC_DIR "${CORNERSTONE_INSTALL_DIR}" CACHE PATH "Cornerstone source directory" FORCE)
        set(cornerstone_FOUND TRUE)
    else()
        message(WARNING "Could not find Cornerstone headers at ${CORNERSTONE_INSTALL_DIR}/include")
    endif()
endif()

# If not found, search for Cornerstone headers in standard locations
if(NOT cornerstone_FOUND)
    # Clear cached value if it was NOTFOUND to allow re-searching
    if(cornerstone_INCLUDE_DIR STREQUAL "cornerstone_INCLUDE_DIR-NOTFOUND")
        unset(cornerstone_INCLUDE_DIR CACHE)
    endif()
    
    find_path(cornerstone_INCLUDE_DIR
        NAMES cstone/domain/domain.hpp
        PATHS
            ${CMAKE_BINARY_DIR}/_deps/cornerstone_fetch-src/include
            ${CMAKE_INSTALL_PREFIX}/include
            $ENV{CORNERSTONE_DIR}/include
            /usr/include
            /usr/local/include
            ${CORNERSTONE_INSTALL_DIR}/include
        DOC "Cornerstone header directory"
    )
    if(cornerstone_INCLUDE_DIR)
        get_filename_component(CORNERSTONE_SRC_DIR "${cornerstone_INCLUDE_DIR}" DIRECTORY)
        set(CORNERSTONE_SRC_DIR "${CORNERSTONE_SRC_DIR}" CACHE PATH "Cornerstone source directory" FORCE)
        message(STATUS "Found Cornerstone headers: ${cornerstone_INCLUDE_DIR}")
        set(cornerstone_FOUND TRUE)
    endif()
endif()

# Handle CUDA-specific configurations
if((MARS_ENABLE_CUDA OR MARS_ENABLE_HIP) AND cornerstone_FOUND)
    find_library(CORNERSTONE_GPU_LIBRARY
        NAMES cstone_gpu libcstone_gpu
        PATHS
            "${CORNERSTONE_SRC_DIR}"
            "${CORNERSTONE_SRC_DIR}/lib"
            "${CORNERSTONE_SRC_DIR}/lib64"
            "${CORNERSTONE_SRC_DIR}/build"
            "${CMAKE_BINARY_DIR}/cornerstone_build"
            "${CMAKE_BINARY_DIR}/_deps/cornerstone_fetch-build"
        NO_DEFAULT_PATH
    )
    if(CORNERSTONE_GPU_LIBRARY)
        message(STATUS "Found Cornerstone GPU library: ${CORNERSTONE_GPU_LIBRARY}")
    else()
        message(WARNING "Cornerstone GPU library not found. Specify manually using -DCORNERSTONE_GPU_LIBRARY=/path/to/libcstone_gpu.a")
    endif()
endif()

# Create an interface library for Cornerstone
if(cornerstone_FOUND AND NOT TARGET cornerstone::cornerstone)
    if(cornerstone_INCLUDE_DIR AND NOT cornerstone_INCLUDE_DIR STREQUAL "cornerstone_INCLUDE_DIR-NOTFOUND")
        add_library(cornerstone::cornerstone INTERFACE IMPORTED)
        set_target_properties(cornerstone::cornerstone PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${cornerstone_INCLUDE_DIR}"
        )
    else()
        message(WARNING "Cornerstone headers not found, skipping target creation")
        set(cornerstone_FOUND FALSE)
    endif()
endif()

# Create GPU target if the library was found
if(CORNERSTONE_GPU_LIBRARY AND NOT TARGET cstone_gpu)
    # Validate that include directory is valid
    if(NOT cornerstone_INCLUDE_DIR OR cornerstone_INCLUDE_DIR STREQUAL "cornerstone_INCLUDE_DIR-NOTFOUND")
        message(FATAL_ERROR "Cannot create cstone_gpu target: cornerstone_INCLUDE_DIR is not set properly")
    endif()
    
    add_library(cstone_gpu STATIC IMPORTED GLOBAL)
    set_target_properties(cstone_gpu PROPERTIES
        IMPORTED_LOCATION "${CORNERSTONE_GPU_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${cornerstone_INCLUDE_DIR}"
    )
    if(MARS_ENABLE_CUDA)
        find_package(CUDAToolkit REQUIRED)
        target_link_libraries(cstone_gpu INTERFACE CUDA::cudart)
        set_property(TARGET cstone_gpu PROPERTY 
            INTERFACE_COMPILE_DEFINITIONS "CSTONE_WITH_CUDA;USE_CUDA")
        message(STATUS "Configured cstone_gpu target with CUDA support")
    elseif(MARS_ENABLE_HIP)
        if(NOT TARGET mars::rocmlibs)
            include(${CMAKE_SOURCE_DIR}/cmake/rocmlibs_target.cmake)
        endif()
        target_link_libraries(cstone_gpu INTERFACE mars::rocmlibs)
        set_property(TARGET cstone_gpu PROPERTY 
            INTERFACE_COMPILE_DEFINITIONS "CSTONE_WITH_HIP;USE_CUDA;THRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_HIP")
        message(STATUS "Configured cstone_gpu target with HIP support")
    else()
        message(WARNING "No GPU support enabled. cstone_gpu target will not be linked.")
    endif()
    message(STATUS "Created cstone_gpu target with library: ${CORNERSTONE_GPU_LIBRARY}")
endif()

# Mark variables as advanced
mark_as_advanced(cornerstone_INCLUDE_DIR CORNERSTONE_GPU_LIBRARY)