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
endfunction()

# Modify the addMarsGpuMpiTest function to support both CUDA and HIP
# and properly register the test without rank constraints
function(addMarsGpuMpiTest source exename testname)
    # Create the executable directly
    add_executable(${exename} ${source})
    
    # Include all necessary directories
    target_include_directories(${exename} PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_BINARY_DIR}
        ${CMAKE_BINARY_DIR}/cornerstone_includes
        ${CORNERSTONE_INCLUDE_DIRS}
        ${MPI_CXX_INCLUDE_PATH}
    )
    
    # Basic compile options
    target_compile_options(${exename} PRIVATE -Wno-unknown-pragmas)
    
    # Link with standard test libraries first
    target_link_libraries(${exename} PRIVATE
        cstone_gpu
        mars_unstructured
        ${MPI_CXX_LIBRARIES}
        ${MARS_G_TEST_LIBRARIES}
    )

        # ADD THESE LINES - Force-link GPU library to include all symbols
    if(TARGET cstone_gpu)
        get_target_property(GPU_LIB cstone_gpu IMPORTED_LOCATION)
        if(GPU_LIB)
            message(STATUS "Force-linking GPU library: ${GPU_LIB}")
            if(UNIX AND NOT APPLE)
                # Linux syntax
                target_link_options(${exename} PRIVATE "LINKER:--whole-archive,${GPU_LIB},--no-whole-archive")
            elseif(APPLE)
                # macOS syntax
                target_link_options(${exename} PRIVATE "LINKER:-force_load,${GPU_LIB}")
            endif()
        endif()
    endif()

    # Add this inside your function where you set up includes
    if(TARGET cstone_gpu)
        get_target_property(GPU_INCLUDES cstone_gpu INTERFACE_INCLUDE_DIRECTORIES)
        if(GPU_INCLUDES)
            message(STATUS "Found include directories from cstone_gpu: ${GPU_INCLUDES}")
            target_include_directories(${exename} PRIVATE ${GPU_INCLUDES})
        endif()
    endif()

    # Add GPU-specific libraries and settings based on enabled backend
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
            CUDA_RESOLVE_DEVICE_SYMBOLS ON
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

        # Replace with these lines
    add_test(NAME ${testname} COMMAND ${exename})
    set_property(GLOBAL APPEND PROPERTY MARS_TESTS ${testname})
    set_property(GLOBAL APPEND PROPERTY MARS_TEST_TARGETS ${exename})

    # Print diagnostic info
    message(STATUS "Added GPU test: ${testname}")
    message(STATUS "Test executable: ${exename}")
    get_target_property(LINK_LIBS ${exename} LINK_LIBRARIES)
    message(STATUS "GPU Test link libraries: ${LINK_LIBS}")
endfunction()

addMarsMpiTest(
    ${CMAKE_CURRENT_SOURCE_DIR}/mars_read_mesh_binary_basic.cpp
    mars_read_mesh_binary_basic_test
    marsReadMeshBinaryBasicTest
)

addMarsMpiTest(
    ${CMAKE_CURRENT_SOURCE_DIR}/mars_read_mesh_binary_mpi.cpp
    mars_read_mesh_binary_mpi_test
    marsReadMeshBinaryMpiTest
)

addMarsMpiTest(
    ${CMAKE_CURRENT_SOURCE_DIR}/mars_read_mesh_binary_external_mesh.cpp
    mars_read_mesh_binary_external_mesh_test
    marsReadMeshBinaryExternalMeshTest
)

# Add domain binary mesh test
addMarsGpuMpiTest(
    ${CMAKE_CURRENT_SOURCE_DIR}/mars_domain_binary_mesh.cu
    mars_domain_binary_test
    MarsDomainBinaryTest
    MarsDomainBinaryTest
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
