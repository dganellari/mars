# First, find required packages
find_package(CUDAToolkit REQUIRED)
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

# Function for adding GPU-enabled tests
function(addMarsGpuMpiTest source exename testname ranks)
    addMarsMpiTest("${source}" ${exename} ${testname} ${ranks})
    
    # Add mars_unstructured which contains domain.cu implementation
    target_link_libraries(${exename} PRIVATE mars_unstructured)
    
    # Add Cornerstone includes directly
    target_include_directories(${exename} PRIVATE
        ${CMAKE_BINARY_DIR}/cornerstone_includes
        ${CORNERSTONE_INCLUDE_DIRS}
    )
    
    # DIRECTLY link with the specified GPU library from command line
    if(EXISTS "${CORNERSTONE_GPU_LIBRARY}")
        message(STATUS "Directly linking with GPU library: ${CORNERSTONE_GPU_LIBRARY}")
        target_link_libraries(${exename} PRIVATE ${CORNERSTONE_GPU_LIBRARY})
    endif()
    
    # Add CUDA libraries
    target_link_libraries(${exename} PRIVATE 
        CUDA::cudart 
        CUDA::cuda_driver
        -lcudadevrt)
    
    # Print diagnostic messages
    message(STATUS "Test executable: ${exename}")
    message(STATUS "CORNERSTONE_GPU_LIBRARY: ${CORNERSTONE_GPU_LIBRARY}")
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