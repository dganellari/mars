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
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr --extended-lambda -std=c++17" CACHE STRING "CUDA flags for all targets" FORCE)
endif()

# Find Cornerstone
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../../cmake")
find_package(cornerstone REQUIRED)

if(cornerstone_FOUND)
    # Add explicit include directory for cornerstone
    message(STATUS "Cornerstone include directory: ${cornerstone_INCLUDE_DIR}")
    # Get the actual source directory from FetchContent
    set(CORNERSTONE_SRC_DIR "${CMAKE_BINARY_DIR}/_deps/cornerstone_fetch-src")
    # Check if headers exist in the expected location
    if(EXISTS "${CORNERSTONE_SRC_DIR}/include/cstone/domain/domain.hpp")
        message(STATUS "Found domain.hpp in expected location")

        # Create an interface library for cornerstone
        add_library(cornerstone_headers INTERFACE)
        target_include_directories(cornerstone_headers INTERFACE
            $<BUILD_INTERFACE:${CORNERSTONE_SRC_DIR}/include>
            $<INSTALL_INTERFACE:include>
        )
        # Add cornerstone_headers to the list of targets to install
        list(APPEND MARS_TARGETS cornerstone_headers)
    else()
        message(STATUS "Could not find domain.hpp in expected location")
        # Try alternate approach - copy headers to build directory include path
        file(COPY ${CORNERSTONE_SRC_DIR}/include/ DESTINATION ${CMAKE_BINARY_DIR}/include/)
        message(STATUS "Copied cornerstone headers to ${CMAKE_BINARY_DIR}/include/")
    endif()
else()
    message(FATAL_ERROR "Cornerstone not found. Unstructured mesh support will be limited.")
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
    
    # 7. Cornerstone integration
    if(TARGET cornerstone_headers)
        target_link_libraries(mars_unstructured PUBLIC cornerstone_headers)
    endif()
    
    # 8. Sources and language properties
    target_sources(mars_unstructured PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/domain.cu)
    set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/domain.cu PROPERTIES LANGUAGE CUDA)
endif()

# These can stay outside the guard
list(APPEND MARS_TARGETS mars_unstructured)
set(MARS_TARGETS ${MARS_TARGETS} PARENT_SCOPE)