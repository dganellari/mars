# Set C++20 for host code, C++17 for CUDA (matching Cornerstone's approach)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find MPI first (like Cornerstone does)
find_package(MPI REQUIRED)
message(STATUS "MPI_CXX_COMPILER: ${MPI_CXX_COMPILER}")

# ENABLE CUDA LANGUAGE - SIMPLIFIED TO MATCH CORNERSTONE
if(MARS_ENABLE_CUDA)
    # Check CUDA compiler first
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)

        # Add these specific CUDA flags right after enabling CUDA language
        if(MARS_ENABLE_CUDA)
            # CRITICAL: Force C++17 for CUDA files with direct compiler flag
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -std=c++18 -ccbin=${MPI_CXX_COMPILER}" CACHE STRING "CUDA flags" FORCE)
            
            # These flags help with template metaprogramming in device code
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr --extended-lambda" CACHE STRING "CUDA flags" FORCE)
            
            # Pass these flags to Cornerstone as well - IMPORTANT: Include -ccbin flag!
            set(CSTONE_EXTRA_CUDA_FLAGS "-std=c++17 -ccbin=${MPI_CXX_COMPILER} --expt-relaxed-constexpr --extended-lambda" CACHE STRING "Extra CUDA flags" FORCE)
            
            # Also set CUDA flags directly for Cornerstone - IMPORTANT: Include -ccbin flag!
            set(CSTONE_CUDA_FLAGS "-std=c++17 -ccbin=${MPI_CXX_COMPILER} --expt-relaxed-constexpr --extended-lambda" CACHE STRING "CUDA flags for Cornerstone" FORCE)
            
            # Add preprocessor definitions that might help Cornerstone's implementation
            add_compile_definitions(
                CSTONE_DEVICE_CODE_COMPATIBILITY=1
                _GLIBCXX_USE_CXX11_ABI=1
            )
        endif()
               
        # Find CUDAToolkit (same as Cornerstone)
        find_package(CUDAToolkit REQUIRED)
        include_directories(SYSTEM ${CUDAToolkit_INCLUDE_DIRS})
        
        # Add basic definitions only
        add_compile_definitions(USE_CUDA)
        
        # Set standard exactly like Cornerstone
        set(CMAKE_CUDA_STANDARD 17)
        
        # Set host compiler
        set(CMAKE_CUDA_HOST_COMPILER ${MPI_CXX_COMPILER})
        
        # Apply CRAY workaround if needed (copied directly from Cornerstone)
        if (MPI_CXX_INCLUDE_PATH STREQUAL "" AND DEFINED ENV{CRAY_MPICH2_DIR} OR DEFINED ENV{CRAY_MPICH_DIR})
            message(STATUS "Applying MPI include workaround on CRAY")
            set(MPI_CXX_INCLUDE_PATH $ENV{MPICH_DIR}/include)
        endif()
    else()
        message(STATUS "No CUDA support")
        set(MARS_ENABLE_CUDA OFF)
    endif()
endif()

# Add cornerstone patch to fix template issues
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/utils/cornerstone_patch.hpp
    ${CMAKE_CURRENT_BINARY_DIR}/utils/cornerstone_patch.hpp
    COPYONLY
)
include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
add_compile_definitions(USE_CORNERSTONE_PATCH)

# Find ADIOS2
find_package(ADIOS2 QUIET)

# Try to find pre-built Cornerstone first
set(CORNERSTONE_INSTALL_DIR "" CACHE PATH "Path to Cornerstone installation directory")

if(CORNERSTONE_INSTALL_DIR)
    # Try to use pre-built version first
    set(Cornerstone_DIR "${CORNERSTONE_INSTALL_DIR}/lib/cmake/cornerstone" CACHE PATH "Path to Cornerstone CMake config" FORCE)
    find_package(cornerstone CONFIG QUIET)
    
    if(cornerstone_FOUND)
        message(STATUS "Using pre-built Cornerstone from: ${CORNERSTONE_INSTALL_DIR}")
        message(STATUS "Cornerstone version: ${cornerstone_VERSION}")
        set(USING_PREBUILT_CORNERSTONE TRUE)
    else()
        message(WARNING "Could not find pre-built Cornerstone at ${CORNERSTONE_INSTALL_DIR}, falling back to build from source")
        set(USING_PREBUILT_CORNERSTONE FALSE)
    endif()
else()
    set(USING_PREBUILT_CORNERSTONE FALSE)
    message(STATUS "No CORNERSTONE_INSTALL_DIR specified, building from source")
endif()

# If no pre-built version, use source from FetchContent
if(NOT USING_PREBUILT_CORNERSTONE)
    # Find Cornerstone source
    find_package(cornerstone REQUIRED)
    set(CORNERSTONE_SRC_DIR "${CMAKE_BINARY_DIR}/_deps/cornerstone_fetch-src")
    
    # Manually set include directories since they're not being populated
    set(CORNERSTONE_INCLUDE_DIRS "${CORNERSTONE_SRC_DIR}/include" 
        CACHE PATH "Cornerstone include directories" FORCE)
    
    # Path to the tuple_util.hpp file that needs patching
    set(TUPLE_UTIL_FILE "${CORNERSTONE_SRC_DIR}/include/cstone/util/tuple_util.hpp")

    # Check if the file exists
    if(EXISTS "${TUPLE_UTIL_FILE}")
        # Read the content
        file(READ "${TUPLE_UTIL_FILE}" TUPLE_UTIL_CONTENT)
        
        # Create a backup if needed
        if(NOT EXISTS "${TUPLE_UTIL_FILE}.original")
            file(WRITE "${TUPLE_UTIL_FILE}.original" "${TUPLE_UTIL_CONTENT}")
        endif()

        # Find the zipTuples function with a more robust pattern
        string(REGEX MATCH "template[^{]*zipTuples[^{]*\\{" FOUND_ZIPTUPLE "${TUPLE_UTIL_CONTENT}")
    
        string(REGEX REPLACE 
            "(template[ \t]*<[ \t]*class[ \t]*\\.\\.\\.[ \t]*Tps[ \t]*>[ \t]*\n[ \t]*constexpr[ \t]+auto[ \t]+zipTuples[ \t]*\\()" 
            "#ifndef CSTONE_UTIL_TUPLE_UTIL_HPP_ZIPTUPLE_DEFINED\\n\\1" 
            MODIFIED_CONTENT "${TUPLE_UTIL_CONTENT}")    

        string(REGEX REPLACE 
            "(return std::apply\\(zip, makeIntegralTuple\\(std::make_index_sequence<N>\\{\\}\\)\\);[^\n]*\n[ \t]*)(\\})" 
            "\\1\\2\\n#endif // CSTONE_UTIL_TUPLE_UTIL_HPP_ZIPTUPLE_DEFINED\\n" 
            MODIFIED_CONTENT "${MODIFIED_CONTENT}")

        # Write the modified content back to the file
        file(WRITE "${TUPLE_UTIL_FILE}" "${MODIFIED_CONTENT}")
    else()
        message(WARNING "Could not find ${TUPLE_UTIL_FILE} to apply patch")
    endif()

    # CONFIGURE CORNERSTONE - MATCH CORNERSTONE'S SETTINGS EXACTLY
    if(MARS_ENABLE_CUDA)
        set(CSTONE_WITH_CUDA ON CACHE BOOL "Build with CUDA support" FORCE)
        set(CSTONE_WITH_GPU_AWARE_MPI ON CACHE BOOL "Enable GPU-aware MPI" FORCE)
        set(CSTONE_CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES} CACHE STRING "CUDA architectures" FORCE)
        
        # IMPORTANT: Include -ccbin flag to match your standalone build
        set(CSTONE_CUDA_FLAGS "-std=c++17 -ccbin=${MPI_CXX_COMPILER} --expt-relaxed-constexpr --extended-lambda" 
            CACHE STRING "CUDA flags for Cornerstone" FORCE)
            
        # Set extra flags for both C++ and CUDA
        set(CSTONE_EXTRA_CXX_FLAGS "-D_GLIBCXX_USE_CXX11_ABI=1" CACHE STRING "Extra C++ flags for Cornerstone" FORCE)
        set(CSTONE_EXTRA_CUDA_FLAGS "-std=c++17 -ccbin=${MPI_CXX_COMPILER} --expt-relaxed-constexpr --extended-lambda" 
            CACHE STRING "Extra CUDA flags for Cornerstone" FORCE)
    endif()

    # Include cornerstone headers (before creating targets)
    include_directories(${CORNERSTONE_INCLUDE_DIRS})

    # Add cornerstone as subdirectory (if not already added)
    if(NOT TARGET cornerstone AND NOT TARGET cstone_gpu)
        # Force C++17 for CUDA in all ways possible
        set(CMAKE_CUDA_FLAGS_INIT "-std=c++17 -ccbin=${MPI_CXX_COMPILER}")
        set(CMAKE_CUDA_FLAGS "-std=c++17 -ccbin=${MPI_CXX_COMPILER}" CACHE STRING "CUDA flags" FORCE)
        
        # Add diagnostic message to debug CUDA flags
        message(STATUS "Building Cornerstone with CUDA flags: ${CMAKE_CUDA_FLAGS}")
        message(STATUS "Building Cornerstone with CSTONE_CUDA_FLAGS: ${CSTONE_CUDA_FLAGS}")
        message(STATUS "Building Cornerstone with CSTONE_EXTRA_CUDA_FLAGS: ${CSTONE_EXTRA_CUDA_FLAGS}")
        
        # Set the environment variable to influence the compiler
        set(ENV{NVCC_PREPEND_FLAGS} "-std=c++17 -ccbin=${MPI_CXX_COMPILER}")
        
        # Add the subdirectory with specified options
        add_subdirectory(${CORNERSTONE_SRC_DIR} ${CMAKE_BINARY_DIR}/cornerstone_build EXCLUDE_FROM_ALL)
    endif()
else()
    # Using pre-built, make sure we still have include paths
    include_directories(${cornerstone_INCLUDE_DIRS})
endif()

# Add diagnostic messages
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "CUDA Standard: ${CMAKE_CUDA_STANDARD}")
message(STATUS "CUDA Flags: ${CMAKE_CUDA_FLAGS}")

# CREATE TARGET
add_library(mars_unstructured "")

# Add sources first
target_sources(mars_unstructured PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/domain.cu)
set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/domain.cu PROPERTIES LANGUAGE CUDA)

# Set CUDA properties
if(MARS_ENABLE_CUDA)
    set_target_properties(mars_unstructured PROPERTIES 
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_HOST_COMPILER ${MPI_CXX_COMPILER})
endif()

# Fix: Use generator expressions for include directories
target_include_directories(mars_unstructured PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/mars/backend/distributed/unstructured>
    ${ADIOS2_INCLUDE_DIRS}
    ${MPI_CXX_INCLUDE_DIRS})

# Link dependencies - EXACTLY LIKE CORNERSTONE'S APPROACH
target_link_libraries(mars_unstructured PUBLIC 
    mars
    MPI::MPI_CXX
    ${ADIOS2_LIBRARIES}
    cornerstone::cornerstone)

# Add CUDA runtime link - EXACTLY LIKE CORNERSTONE'S GPU TESTS
if(MARS_ENABLE_CUDA)
    target_link_libraries(mars_unstructured PUBLIC CUDA::cudart)
endif()

# Link with cstone_gpu if available
if(TARGET cstone_gpu)
    target_link_libraries(mars_unstructured PUBLIC cstone_gpu)
endif()

# Export targets
list(APPEND MARS_TARGETS mars_unstructured)
if(TARGET cstone_gpu)
    list(APPEND MARS_TARGETS cstone_gpu)
endif()
set(MARS_TARGETS ${MARS_TARGETS} PARENT_SCOPE)