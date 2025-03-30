# Set C++20 for host code, C++17 for CUDA (matching Cornerstone's approach)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find MPI first (like Cornerstone does)
find_package(MPI REQUIRED)
message(STATUS "MPI_CXX_COMPILER: ${MPI_CXX_COMPILER}")

# ENABLE CUDA LANGUAGE
if(MARS_ENABLE_CUDA)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)

        # Force C++17 for CUDA files with necessary compiler flags
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -std=c++17 -ccbin=${MPI_CXX_COMPILER} --expt-relaxed-constexpr --extended-lambda" 
            CACHE STRING "CUDA flags" FORCE)
            
        # Pass flags to Cornerstone
        set(CSTONE_CUDA_FLAGS "-std=c++17 -ccbin=${MPI_CXX_COMPILER} --expt-relaxed-constexpr --extended-lambda" 
            CACHE STRING "CUDA flags for Cornerstone" FORCE)
        set(CSTONE_EXTRA_CUDA_FLAGS "${CSTONE_CUDA_FLAGS}" 
            CACHE STRING "Extra CUDA flags for Cornerstone" FORCE)
            
        # Add preprocessor definitions
        add_compile_definitions(
            CSTONE_DEVICE_CODE_COMPATIBILITY=1
            _GLIBCXX_USE_CXX11_ABI=1
            USE_CUDA
        )
               
        # Find CUDAToolkit
        find_package(CUDAToolkit REQUIRED)
        include_directories(SYSTEM ${CUDAToolkit_INCLUDE_DIRS})
        
        # Set CUDA standard
        set(CMAKE_CUDA_STANDARD 17)
        set(CMAKE_CUDA_HOST_COMPILER ${MPI_CXX_COMPILER})
        
        # Apply CRAY workaround if needed
        if (MPI_CXX_INCLUDE_PATH STREQUAL "" AND (DEFINED ENV{CRAY_MPICH2_DIR} OR DEFINED ENV{CRAY_MPICH_DIR}))
            message(STATUS "Applying MPI include workaround on CRAY")
            set(MPI_CXX_INCLUDE_PATH $ENV{MPICH_DIR}/include)
        endif()
    else()
        message(STATUS "No CUDA support")
        set(MARS_ENABLE_CUDA OFF)
    endif()
endif()

# Find ADIOS2
find_package(ADIOS2 QUIET)

# Find Cornerstone - either installed or fetched
set(CORNERSTONE_INSTALL_DIR "" CACHE PATH "Path to Cornerstone installation directory")
find_package(cornerstone REQUIRED)

# Check if we're using pre-built or source version
if(CORNERSTONE_SRC_DIR)
    set(USING_PREBUILT_CORNERSTONE FALSE)
    message(STATUS "Using Cornerstone from source")
    
    # Add cornerstone patch only for source version
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/utils/cornerstone_patch.hpp
        ${CMAKE_CURRENT_BINARY_DIR}/utils/cornerstone_patch.hpp
        COPYONLY
    )
    include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
    
    # Set include directories
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

        # Find the zipTuples function
        string(REGEX MATCH "template[^{]*zipTuples[^{]*\\{" FOUND_ZIPTUPLE "${TUPLE_UTIL_CONTENT}")
    
        # Add #ifndef guard before the function
        string(REGEX REPLACE 
            "(template[ \t]*<[ \t]*class[ \t]*\\.\\.\\.[ \t]*Tps[ \t]*>[ \t]*\n[ \t]*constexpr[ \t]+auto[ \t]+zipTuples[ \t]*\\()" 
            "#ifndef CSTONE_UTIL_TUPLE_UTIL_HPP_ZIPTUPLE_DEFINED\\n\\1" 
            MODIFIED_CONTENT "${TUPLE_UTIL_CONTENT}")    

        # Add #endif guard after the function
        string(REGEX REPLACE 
            "(return std::apply\\(zip, makeIntegralTuple\\(std::make_index_sequence<N>\\{\\}\\)\\);[^\n]*\n[ \t]*)(\\})" 
            "\\1\\2\\n#endif // CSTONE_UTIL_TUPLE_UTIL_HPP_ZIPTUPLE_DEFINED\\n" 
            MODIFIED_CONTENT "${MODIFIED_CONTENT}")

        # Write the modified content back to the file
        file(WRITE "${TUPLE_UTIL_FILE}" "${MODIFIED_CONTENT}")
        message(STATUS "Applied patch to ${TUPLE_UTIL_FILE}")
    else()
        set(USING_PREBUILT_CORNERSTONE TRUE)
        message(STATUS "Using pre-built Cornerstone - no patching needed")
        # Using pre-built, make sure we still have include paths
        include_directories(${cornerstone_INCLUDE_DIRS}) 
    endif()

    # Configure Cornerstone CUDA settings
    if(MARS_ENABLE_CUDA)
        set(CSTONE_WITH_CUDA ON CACHE BOOL "Build with CUDA support" FORCE)
        set(CSTONE_WITH_GPU_AWARE_MPI ON CACHE BOOL "Enable GPU-aware MPI" FORCE)
        set(CSTONE_CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES} CACHE STRING "CUDA architectures" FORCE)
        set(CSTONE_EXTRA_CXX_FLAGS "-D_GLIBCXX_USE_CXX11_ABI=1" CACHE STRING "Extra C++ flags for Cornerstone" FORCE)
    endif()

    # Include cornerstone headers
    include_directories(${CORNERSTONE_INCLUDE_DIRS})

    # Add cornerstone as subdirectory if not already added
    if(NOT TARGET cornerstone AND NOT TARGET cstone_gpu)
        # Set environment variable for NVCC
        set(ENV{NVCC_PREPEND_FLAGS} "-std=c++17 -ccbin=${MPI_CXX_COMPILER}")
        
        # Add debug info
        message(STATUS "Building Cornerstone with CUDA flags: ${CMAKE_CUDA_FLAGS}")
        
        # Add the subdirectory with specified options
        add_subdirectory(${CORNERSTONE_SRC_DIR} ${CMAKE_BINARY_DIR}/cornerstone_build EXCLUDE_FROM_ALL)
    endif()
else()
    set(USING_PREBUILT_CORNERSTONE TRUE)
    message(STATUS "Using pre-built Cornerstone - no patching needed")
    
    # Using pre-built, make sure we still have include paths
    include_directories(${cornerstone_INCLUDE_DIRS})
endif()

# Add diagnostic messages
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "CUDA Standard: ${CMAKE_CUDA_STANDARD}")
message(STATUS "CUDA Flags: ${CMAKE_CUDA_FLAGS}")

# Create mars_unstructured target
add_library(mars_unstructured "")
target_sources(mars_unstructured PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/domain.cu)
set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/domain.cu PROPERTIES LANGUAGE CUDA)

# Set CUDA properties
if(MARS_ENABLE_CUDA)
    set_target_properties(mars_unstructured PROPERTIES 
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_HOST_COMPILER ${MPI_CXX_COMPILER})
endif()

# Set include directories with generator expressions
target_include_directories(mars_unstructured PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/mars/backend/distributed/unstructured>
    ${ADIOS2_INCLUDE_DIRS}
    ${MPI_CXX_INCLUDE_DIRS})

# Link dependencies
target_link_libraries(mars_unstructured PUBLIC 
    mars
    MPI::MPI_CXX
    ${ADIOS2_LIBRARIES}
    cornerstone::cornerstone)

# Add CUDA runtime link
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