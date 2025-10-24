# Set C++20 for this target
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Create main library target
add_library(mars_unstructured "")

# Add the implementation files to the target
target_sources(mars_unstructured PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/domain.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/domain_cuda_impl.cpp
)

# Set language properties based on enabled features
if(MARS_ENABLE_HIP)
    # Default architecture if not defined
    if(NOT DEFINED HIP_GPU_ARCHITECTURES OR "${HIP_GPU_ARCHITECTURES}" STREQUAL "")
        set(HIP_GPU_ARCHITECTURES "gfx942" CACHE STRING "AMD GPU target architecture" FORCE)
    endif()

    # Use HIP for GPU files
    set_source_files_properties(
        ${CMAKE_CURRENT_SOURCE_DIR}/domain.cu 
        PROPERTIES LANGUAGE HIP
    )
    
    # Add HIP-specific instantiations file if it exists
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/cstone_instantiations.hip.cpp")
        target_sources(mars_unstructured PRIVATE 
            ${CMAKE_CURRENT_SOURCE_DIR}/cstone_instantiations.hip.cpp
        )
        set_source_files_properties(
            ${CMAKE_CURRENT_SOURCE_DIR}/cstone_instantiations.hip.cpp 
            PROPERTIES LANGUAGE HIP
        )
    endif()
    
    # HIP-specific properties
    set_target_properties(mars_unstructured PROPERTIES
        HIP_ARCHITECTURES "${HIP_GPU_ARCHITECTURES}"
        HIP_SEPARABLE_COMPILATION ON)

    message(STATUS "HIP_ARCHITECTURES: ${HIP_GPU_ARCHITECTURES}")
        
    # Add HIP-specific defines
    target_compile_definitions(mars_unstructured PRIVATE
        __HIP_PLATFORM_AMD__=1
        $<$<COMPILE_LANGUAGE:HIP>:THRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_HIP>
        $<$<COMPILE_LANGUAGE:HIP>:THRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_CPP>)
elseif(MARS_ENABLE_CUDA)
    # Use CUDA for GPU files
    set_source_files_properties(
        ${CMAKE_CURRENT_SOURCE_DIR}/domain.cu 
        PROPERTIES LANGUAGE CUDA
    )
    # Add CUDA-specific instantiations if needed
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/cstone_instantiations.cu")
        target_sources(mars_unstructured PRIVATE cstone_instantiations.cu)
        set_source_files_properties(cstone_instantiations.cu PROPERTIES LANGUAGE CUDA)
    endif()
    
    # CUDA-specific properties
    set_target_properties(mars_unstructured PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_HOST_COMPILER ${MPI_CXX_COMPILER})
endif()

# For domain_cuda_impl.cpp, set language based on enabled backends
if(MARS_ENABLE_HIP)
    # Use HIP language for this file when HIP is enabled
    set_source_files_properties(
        ${CMAKE_CURRENT_SOURCE_DIR}/domain_cuda_impl.cpp 
        PROPERTIES LANGUAGE HIP
    )
    message(STATUS "Set domain_cuda_impl.cpp to use HIP language")
elseif(MARS_ENABLE_CUDA)
    # For CUDA, keep as CXX but add CUDA include paths
    set_source_files_properties(
        ${CMAKE_CURRENT_SOURCE_DIR}/domain_cuda_impl.cpp 
        PROPERTIES 
            LANGUAGE CXX
            COMPILE_FLAGS "-isystem${CUDAToolkit_INCLUDE_DIRS}"
    )
    message(STATUS "Set domain_cuda_impl.cpp to use CXX language with CUDA includes")
else()
    # Regular C++ for CPU builds
    set_source_files_properties(
        ${CMAKE_CURRENT_SOURCE_DIR}/domain_cuda_impl.cpp 
        PROPERTIES LANGUAGE CXX
    )
endif()

# Set include directories
target_include_directories(mars_unstructured PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/mars/backend/distributed/unstructured>)

# Add cornerstone includes from main setup
target_include_directories(mars_unstructured BEFORE PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/cornerstone_includes>)

# Make sure we can find mars.hpp and other key headers
target_include_directories(mars_unstructured PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/base
    ${CMAKE_BINARY_DIR})

# Link against dependencies
target_link_libraries(mars_unstructured PUBLIC
    mars
    MPI::MPI_CXX)

# Link to Cornerstone if available
if(TARGET cornerstone::cornerstone)
    target_link_libraries(mars_unstructured PUBLIC cornerstone::cornerstone)
endif()

if(TARGET cstone_gpu)
    target_link_libraries(mars_unstructured PUBLIC cstone_gpu)
endif()

# Export targets
list(APPEND MARS_TARGETS mars_unstructured)
set(MARS_TARGETS ${MARS_TARGETS} PARENT_SCOPE)