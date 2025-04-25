# First, find required packages
if(MARS_ENABLE_CUDA)
    find_package(CUDAToolkit REQUIRED)
elseif(MARS_ENABLE_HIP)
    include(${CMAKE_SOURCE_DIR}/cmake/rocmlibs_target.cmake)
endif()
find_package(MPI REQUIRED)

# Ensure testing is enabled here
enable_testing()

# Define a helper function for adding Mars MPI tests
function(addMarsMpiTest source exename testname)
    add_executable(${exename} ${source})
    
    # Include directories - include all necessary paths here
    target_include_directories(${exename} PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_BINARY_DIR}
        ${CMAKE_BINARY_DIR}/cornerstone_includes
        ${CORNERSTONE_INCLUDE_DIRS}
        ${MPI_CXX_INCLUDE_PATH}
        ${CMAKE_SOURCE_DIR}/backend/distributed/unstructured/utils
    )
    
    target_compile_options(${exename} PRIVATE -Wno-unknown-pragmas)
    
    # Link with all necessary libraries (except GPU-specific ones)
    target_link_libraries(${exename} PRIVATE
        mars_unstructured  # Link to mars_unstructured for all tests
        cstone_gpu
        ${MPI_CXX_LIBRARIES}
        ${MARS_G_TEST_LIBRARIES}
    )
    
    mars_add_test(${testname} ${source})
    
    message(STATUS "Added test: ${testname}")
endfunction()

# Modify the addMarsGpuMpiTest function to ONLY add GPU-specific functionality
function(addMarsGpuMpiTest source exename testname)
    # First, create the basic MPI test
    addMarsMpiTest("${source}" ${exename} ${testname})
    
    # Add GPU-specific libraries based on enabled backend
    if(MARS_ENABLE_CUDA)
        # For CUDA backend
        target_link_libraries(${exename} PRIVATE 
            CUDA::cudart 
            -lcudadevrt)
        
        # Set CUDA language for the source file
        set_source_files_properties(${source} PROPERTIES LANGUAGE CUDA)
        
        # Set CUDA-specific properties
        set_target_properties(${exename} PROPERTIES 
            CUDA_SEPARABLE_COMPILATION ON
            CUDA_HOST_COMPILER ${MPI_CXX_COMPILER})
            
    elseif(MARS_ENABLE_HIP)
        # For HIP backend
        target_link_libraries(${exename} PRIVATE mars::rocmlibs)
        
        # Set HIP language for the source file
        set_source_files_properties(${source} PROPERTIES LANGUAGE HIP)
        
        # Set HIP-specific properties
        set_target_properties(${exename} PROPERTIES 
            HIP_SEPARABLE_COMPILATION ON
            HIP_ARCHITECTURES "${HIP_GPU_ARCHITECTURES}")
    endif()
    
    # Print diagnostic info
    message(STATUS "Test executable: ${exename}")
    message(STATUS "CORNERSTONE_INCLUDE_DIRS: ${CORNERSTONE_INCLUDE_DIRS}")
endfunction()

# # Add domain binary mesh test
# addMarsGpuMpiTest(
#     ${CMAKE_CURRENT_SOURCE_DIR}/mars_domain_binary_mesh.cu
#     mars_domain_binary_test
#     MarsDomainBinaryTest
#     2
# )

addMarsMpiTest(
    ${CMAKE_CURRENT_SOURCE_DIR}/mars_read_mesh_binary.cpp
    mars_read_mesh_binary_test
    marsReadMeshBinaryTest
)

# Print diagnostic messages to help debug linking
if(TARGET cstone_gpu)
    get_target_property(LOCATION cstone_gpu IMPORTED_LOCATION)
    message(STATUS "cstone_gpu location: ${LOCATION}")
    
    get_target_property(INTERFACE_LIBS cstone_gpu INTERFACE_LINK_LIBRARIES)
    message(STATUS "cstone_gpu interface libraries: ${INTERFACE_LIBS}")
endif()

# Print Cornerstone variables for diagnostic purposes
message(STATUS "CORNERSTONE_INCLUDE_DIRS: ${CORNERSTONE_INCLUDE_DIRS}")
message(STATUS "CORNERSTONE_SRC_DIR: ${CORNERSTONE_SRC_DIR}")
