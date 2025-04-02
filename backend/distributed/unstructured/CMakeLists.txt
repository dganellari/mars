# Set C++20 for this target
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Create main library target
add_library(mars_unstructured "")

# Add sources and set language properties
target_sources(mars_unstructured PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/domain.cu)
set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/domain.cu PROPERTIES LANGUAGE CUDA)

# Configure CUDA properties if enabled
if(MARS_ENABLE_CUDA)
    set_target_properties(mars_unstructured PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_HOST_COMPILER ${MPI_CXX_COMPILER})
endif()

# Set include directories
target_include_directories(mars_unstructured PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/mars/backend/distributed/unstructured>)

# Add cornerstone includes from main setup
target_include_directories(mars_unstructured BEFORE PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/cornerstone_includes>)

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