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
    
    # CUDA-specific properties
    set_target_properties(mars_unstructured PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_HOST_COMPILER ${MPI_CXX_COMPILER})
endif()

# C++ files always use CXX language
set_source_files_properties(
    ${CMAKE_CURRENT_SOURCE_DIR}/domain_cuda_impl.cpp 
    PROPERTIES LANGUAGE CXX
)

# Add after line 47 in your unstructured CMakeLists.txt
# Add include paths directly to the cpp file for GPU header compatibility
if(MARS_ENABLE_HIP)
    set_source_files_properties(
        ${CMAKE_CURRENT_SOURCE_DIR}/domain_cuda_impl.cpp 
        PROPERTIES COMPILE_FLAGS "-isystem${HIP_INCLUDE_DIRS}"
    )
elseif(MARS_ENABLE_CUDA)
    set_source_files_properties(
        ${CMAKE_CURRENT_SOURCE_DIR}/domain_cuda_impl.cpp 
        PROPERTIES COMPILE_FLAGS "-isystem${CUDAToolkit_INCLUDE_DIRS}"
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