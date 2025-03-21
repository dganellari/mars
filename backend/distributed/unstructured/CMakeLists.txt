# Set C++20 specifically for unstructured targets
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Check CMake version for Cornerstone compatibility
if(CMAKE_VERSION VERSION_LESS 3.22)
    message(WARNING "Cornerstone recommends CMake 3.22 or higher, but found ${CMAKE_VERSION}. This may cause issues with the unstructured module.")
endif()

# Handle CMP0146 policy for newer CMake versions
if(POLICY CMP0146)
    cmake_policy(SET CMP0146 NEW)
endif()

# 1. ENABLE CUDA LANGUAGE (if needed)
if(MARS_ENABLE_CUDA AND NOT CMAKE_CUDA_COMPILER_ID)
    # Enable CUDA language
    enable_language(CUDA)

    # Find appropriate CUDA toolkit based on CMake version
    if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.18)
        find_package(CUDAToolkit REQUIRED)
    else()
        find_package(CUDA REQUIRED)
    endif()

    # Add this block to make CUDA headers available to C++ files
    add_compile_definitions(USE_CUDA)
    if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.18)
        include_directories(SYSTEM ${CUDAToolkit_INCLUDE_DIRS})
    else()
        include_directories(SYSTEM ${CUDA_INCLUDE_DIRS})
    endif()

    # CUDA-specific settings for unstructured code
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)

    # Start with completely clean CUDA flags - do NOT inherit from CXX flags
    set(CMAKE_CUDA_FLAGS "" CACHE STRING "CUDA flags" FORCE)
    set(CMAKE_CUDA_FLAGS_DEBUG "" CACHE STRING "CUDA Debug flags" FORCE)
    set(CMAKE_CUDA_FLAGS_RELEASE "-O3" CACHE STRING "CUDA Release flags" FORCE)
    set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-O2 -g" CACHE STRING "CUDA RelWithDebInfo flags" FORCE)
endif()

# 2. FIND DEPENDENCIES

# Find ADIOS2
find_package(ADIOS2 REQUIRED)
if(ADIOS2_FOUND)
    message(STATUS "ADIOS2 found: ${ADIOS2_DIR}")
    # Make variables available as CACHE variables for other modules
    set(ADIOS2_FOUND ${ADIOS2_FOUND} CACHE BOOL "ADIOS2 found flag" FORCE)
    set(ADIOS2_INCLUDE_DIRS ${ADIOS2_INCLUDE_DIRS} CACHE PATH "ADIOS2 include directories" FORCE)
    set(ADIOS2_LIBRARIES ${ADIOS2_LIBRARIES} CACHE PATH "ADIOS2 libraries" FORCE)

    # Check for CUDA support in ADIOS2 if MARS_ENABLE_CUDA is enabled
    if(MARS_ENABLE_CUDA)
        if(NOT ADIOS2_HAVE_CUDA)
            message(FATAL_ERROR "ADIOS2 was found but it does not have CUDA support. Please install ADIOS2 with CUDA support or disable MARS_ENABLE_CUDA.")
        endif()
    endif()
else()
    message(FATAL_ERROR "ADIOS2 not found. Please install ADIOS2 or disable MARS_ENABLE_UNSTRUCTURED.")
endif()

# Set CUDA standard for Cornerstone
set(CSTONE_CUDA_STANDARD 17 CACHE STRING "CUDA standard for Cornerstone" FORCE)
set(CSTONE_CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES} CACHE STRING "CUDA architectures for Cornerstone" FORCE)
set(CSTONE_WITH_CUDA ON CACHE BOOL "Build Cornerstone with CUDA support" FORCE)

# Set global CUDA flags (this ensures all CUDA code gets these flags)
if(CMAKE_CUDA_COMPILER)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr --extended-lambda" CACHE STRING "CUDA flags for all targets" FORCE)
endif()


# Find Cornerstone source
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../../cmake")
find_package(cornerstone REQUIRED)

if(cornerstone_FOUND)
    # Get the actual source directory from FetchContent
    set(CORNERSTONE_SRC_DIR "${CMAKE_BINARY_DIR}/_deps/cornerstone_fetch-src")
    
    # Check if headers exist in the expected location
    if(EXISTS "${CORNERSTONE_SRC_DIR}/include/cstone/domain/domain.hpp")
        message(STATUS "Found domain.hpp in expected location")

        # FIX LINKING ISSUE - Explicitly build cornerstone with all components
        if(CMAKE_CUDA_COMPILER AND NOT TARGET cornerstone AND NOT TARGET cstone_gpu)
            # Set options for cornerstone
            set(CSTONE_WITH_CUDA ON CACHE BOOL "Build Cornerstone with CUDA support" FORCE)
            set(CSTONE_WITH_GPU_AWARE_MPI ${MARS_ENABLE_CUDA} CACHE BOOL "Enable GPU-aware MPI" FORCE)
            set(CSTONE_BUILD_EXAMPLES OFF CACHE BOOL "Don't build examples" FORCE)
            set(CSTONE_INSTALL_HEADERS ON CACHE BOOL "Install header files" FORCE)
            set(CSTONE_DISABLE_ERROR_ON_WARNINGS ON CACHE BOOL "Disable warnings as errors" FORCE)

            # Build cornerstone first, separately
            execute_process(
                COMMAND ${CMAKE_COMMAND} -E echo "Building cornerstone in ${CMAKE_BINARY_DIR}/cornerstone_build"
                COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_BINARY_DIR}/cornerstone_build"
                COMMAND ${CMAKE_COMMAND} ${CORNERSTONE_SRC_DIR}
                    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                    -DCSTONE_WITH_CUDA=ON
                    -DCSTONE_WITH_GPU_AWARE_MPI=${MARS_ENABLE_CUDA}
                WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/cornerstone_build"
                RESULT_VARIABLE result
            )

            if(result EQUAL 0)
                execute_process(
                    COMMAND ${CMAKE_COMMAND} --build .
                    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/cornerstone_build"
                    RESULT_VARIABLE build_result
                )

                if(build_result EQUAL 0)
                    message(STATUS "Successfully built cornerstone")
                    # Now add as subdirectory for the targets
                    add_subdirectory(${CORNERSTONE_SRC_DIR} ${CMAKE_BINARY_DIR}/cornerstone_build EXCLUDE_FROM_ALL)
                else()
                    message(FATAL_ERROR "Failed to build cornerstone: ${build_result}")
                endif()
            else()
                message(FATAL_ERROR "Failed to configure cornerstone: ${result}")
            endif()
        endif()
        
    endif()
endif()

# Explicitly add cornerstone include directory to domain.cu
if(cornerstone_FOUND)
    # Make sure the include directory is available
    if(EXISTS "${CORNERSTONE_SRC_DIR}/include")
        message(STATUS "Adding cornerstone include path: ${CORNERSTONE_SRC_DIR}/include")
        include_directories(${CORNERSTONE_SRC_DIR}/include)

        #also add to specific files
        set_property(SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/domain.cu 
                    APPEND PROPERTY COMPILE_OPTIONS "-I${CORNERSTONE_SRC_DIR}/include")
    endif()
endif()

# Find MPI
if(DEFINED ENV{MPICH_DIR})
    # Use environment variable if available
    set(MPI_BASE_DIR "$ENV{MPICH_DIR}")
else()
    # Extract from compiler path
    get_filename_component(MPI_COMPILER_PATH "${MPI_CXX_COMPILER}" DIRECTORY)
    get_filename_component(MPI_BASE_DIR "${MPI_COMPILER_PATH}" DIRECTORY)
endif()

# Set the variables that FindMPIExtended.cmake will look for
set(TPL_MPI_INCLUDE_DIRS "${MPI_BASE_DIR}/include" CACHE PATH "MPI Include path")
if(EXISTS "${MPI_BASE_DIR}/lib/libmpi.so")
    set(TPL_MPI_LIBRARIES "${MPI_BASE_DIR}/lib/libmpi.so" CACHE FILEPATH "MPI library")
else()
    # Try to find libmpi.so in other common locations
    find_library(TPL_MPI_LIBRARIES
        NAMES mpi
        PATHS ${MPI_BASE_DIR}/lib ${MPI_BASE_DIR}/lib64
        NO_DEFAULT_PATH)
endif()

find_package(MPI REQUIRED)
message(STATUS "MPI_CXX_INCLUDE_DIRS: ${MPI_CXX_INCLUDE_DIRS}")
message(STATUS "MPI_CXX_COMPILER: ${MPI_CXX_COMPILER}")

# 3. CREATE TARGET AND CONFIGURATION
if(NOT TARGET mars_unstructured)
    # 1. Target creation
    add_library(mars_unstructured "")
    
    # 2. Target properties - CUDA separable compilation
    if(CMAKE_CUDA_COMPILER)
      set_target_properties(mars_unstructured PROPERTIES 
        CUDA_SEPARABLE_COMPILATION ON)
    endif()
    
    # 3. Include directories
    target_include_directories(mars_unstructured PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include/mars/backend/distributed/unstructured>
    )
    
    # 4. All link libraries
    target_link_libraries(mars_unstructured PUBLIC mars)
    
    # 5. MPI-specific configurations
    if("${MPI_CXX_INCLUDE_DIRS}" STREQUAL "")
        target_include_directories(mars_unstructured PUBLIC ${MPI_BASE_DIR}/include)
        if(CMAKE_CUDA_COMPILER)
            set_property(SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/domain.cu
                         PROPERTY COMPILE_OPTIONS "-I${MPI_BASE_DIR}/include")
        endif()
    else()
        target_include_directories(mars_unstructured PUBLIC ${MPI_CXX_INCLUDE_DIRS})
    endif()
    target_link_libraries(mars_unstructured PUBLIC MPI::MPI_CXX)
    
    # 6. ADIOS2 integration
    if(ADIOS2_FOUND)
        target_include_directories(mars_unstructured PUBLIC ${ADIOS2_INCLUDE_DIRS})
        target_link_libraries(mars_unstructured PUBLIC ${ADIOS2_LIBRARIES})
        if(CMAKE_CUDA_COMPILER)
            foreach(ADIOS_DIR ${ADIOS2_INCLUDE_DIRS})
                set_property(SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/domain.cu
                            APPEND PROPERTY COMPILE_OPTIONS "-I${ADIOS_DIR}")
            endforeach()
        endif()
    endif()
    
    # 7. Link the header-only interface
    target_link_libraries(mars_unstructured PUBLIC cornerstone::cornerstone)
    
    # Also link to the implementation target
    if(TARGET cstone_gpu)
        target_link_libraries(mars_unstructured PUBLIC cstone_gpu)
        # Add cstone_gpu to exported targets list
        list(APPEND MARS_TARGETS cstone_gpu)
    endif() 
    
    # 8. Sources and language properties
    target_sources(mars_unstructured PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/domain.cu)
    set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/domain.cu PROPERTIES LANGUAGE CUDA)
endif()

# These can stay outside the guard
list(APPEND MARS_TARGETS mars_unstructured)
set(MARS_TARGETS ${MARS_TARGETS} PARENT_SCOPE)
