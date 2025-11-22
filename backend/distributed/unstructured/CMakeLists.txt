# Set C++20 for this target
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Create main library target
add_library(mars_unstructured STATIC)

# Set include directories
target_include_directories(mars_unstructured PUBLIC
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/base
    ${CMAKE_SOURCE_DIR}/core
    ${CMAKE_SOURCE_DIR}/backend/distributed
    ${CMAKE_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Add all distributed backend subdirectories for header dependencies
target_include_directories(mars_unstructured PRIVATE
    ${CMAKE_SOURCE_DIR}/backend/distributed/mars_env
    ${CMAKE_SOURCE_DIR}/backend/distributed/communication
    ${CMAKE_SOURCE_DIR}/backend/distributed/communication/util
    ${CMAKE_SOURCE_DIR}/backend/distributed/mesh
    ${CMAKE_SOURCE_DIR}/backend/distributed/mesh_manager
    ${CMAKE_SOURCE_DIR}/backend/distributed/octree
    ${CMAKE_SOURCE_DIR}/backend/distributed/sfc
    ${CMAKE_SOURCE_DIR}/backend/distributed/unstructured/utils
)

# System includes (won't propagate warnings)
target_include_directories(mars_unstructured SYSTEM PUBLIC
    ${CMAKE_BINARY_DIR}/cornerstone_includes
    ${CORNERSTONE_INCLUDE_DIRS}
)

# Add the implementation files to the target
target_sources(mars_unstructured PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/domain.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/domain_cuda_impl.cpp
)

# Set language properties based on enabled features
if(MARS_ENABLE_CUDA)
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
        CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}"
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_HOST_COMPILER ${MPI_CXX_COMPILER})
    
    # Suppress future compatibility warning about global function templates
    target_compile_options(mars_unstructured PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:-Xnvlink=--suppress-stack-size-warning>
        $<$<COMPILE_LANGUAGE:CUDA>:--diag-suppress=20281>
    )
endif()

# For domain_cuda_impl.cpp
if(MARS_ENABLE_CUDA)
    set_source_files_properties(
        ${CMAKE_CURRENT_SOURCE_DIR}/domain_cuda_impl.cpp 
        PROPERTIES 
            LANGUAGE CXX
            COMPILE_FLAGS "-isystem${CUDAToolkit_INCLUDE_DIRS}"
    )
endif()

# Link against dependencies
target_link_libraries(mars_unstructured PUBLIC
    ${MARS_DEP_LIBRARIES}
    MPI::MPI_CXX
)

# Link to Cornerstone if available
if(TARGET cornerstone::cornerstone)
    target_link_libraries(mars_unstructured PUBLIC cornerstone::cornerstone)
endif()

# Link to CUDA libraries
if(MARS_ENABLE_CUDA)
    target_link_libraries(mars_unstructured PUBLIC 
        CUDA::cudart 
        CUDA::cusparse 
        CUDA::cublas 
        CUDA::cusolver)
endif()

if(TARGET cstone_gpu)
    target_link_libraries(mars_unstructured PUBLIC cstone_gpu)
endif()

# Export targets
list(APPEND MARS_TARGETS mars_unstructured)
set(MARS_TARGETS ${MARS_TARGETS} PARENT_SCOPE)