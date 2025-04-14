# First, find required packages
if(MARS_ENABLE_CUDA)
    find_package(CUDAToolkit REQUIRED)
elseif(MARS_ENABLE_HIP)
    include(${CMAKE_SOURCE_DIR}/cmake/rocmlibs_target.cmake)
endif()
find_package(MPI REQUIRED)

# Make sure Tests.cmake is included here if not already included higher up
if(NOT COMMAND mars_add_test)
    include(${CMAKE_SOURCE_DIR}/cmake/Tests.cmake)
endif()

# Define a helper function for adding Mars MPI tests
function(addMarsMpiTest source exename testname ranks)
    add_executable(${exename} ${source})
    
    # Include directories - note how Cornerstone uses PROJECT_SOURCE_DIR/include
    target_include_directories(${exename} PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_BINARY_DIR}
        ${CMAKE_BINARY_DIR}/cornerstone_includes
        ${CORNERSTONE_INCLUDE_DIRS}
        ${MPI_CXX_INCLUDE_PATH}
    )
    
    target_compile_options(${exename} PRIVATE -Wno-unknown-pragmas)
    target_link_libraries(${exename} PRIVATE
        ${MPI_CXX_LIBRARIES}
        ${MARS_G_TEST_LIBRARIES}
    )
    
    # Use the mars_add_mpi_test function if available, otherwise just add_test
    if(COMMAND mars_add_mpi_test)
        mars_add_mpi_test(${testname} EXECUTABLE ${exename} RANKS ${ranks})
    else()
        add_test(NAME ${testname} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} $<TARGET_FILE:${exename}>)
    endif()
endfunction()

# Modify the addMarsGpuMpiTest function to support both CUDA and HIP
function(addMarsGpuMpiTest source exename testname ranks)
    addMarsMpiTest("${source}" ${exename} ${testname} ${ranks})
    
    # Add mars_unstructured which contains domain.cu implementation
    target_link_libraries(${exename} PRIVATE mars_unstructured)

    # Add Cornerstone includes directly
    target_include_directories(${exename} PRIVATE
        ${CMAKE_BINARY_DIR}/cornerstone_includes
        ${CORNERSTONE_INCLUDE_DIRS}
    )
    
    # Link to GPU target
    message(STATUS "Linking to GPU target: cstone_gpu")
    target_link_libraries(${exename} PRIVATE cstone_gpu)
    
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

# Add domain binary mesh test
addMarsGpuMpiTest(
    ${CMAKE_CURRENT_SOURCE_DIR}/mars_domain_binary_mesh.cu
    mars_domain_binary_test
    MarsDomainBinaryTest
    2
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
message(STATUS "CORNERSTONE_MAIN_LIB: ${CORNERSTONE_MAIN_LIB}")