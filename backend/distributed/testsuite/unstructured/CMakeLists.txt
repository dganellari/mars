# First, find required packages
find_package(CUDAToolkit REQUIRED)
find_package(MPI REQUIRED)

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
    
    # Add GPU-specific links and includes
    target_link_libraries(${exename} PRIVATE
        mars
        mars_unstructured
    )
    
    # Link to the Cornerstone GPU library if available
    if(TARGET cstone_gpu)
        target_link_libraries(${exename} PRIVATE cstone_gpu)
    endif()
    
    target_link_libraries(${exename} PRIVATE CUDA::cudart)
endfunction()

# Add test using our function
addMarsGpuMpiTest(
    ${CMAKE_CURRENT_SOURCE_DIR}/mars_domain_binary_mesh.cu
    mars_domain_binary_test
    MarsDomainBinaryTest
    2
)

# Add extra include directories specifically for this test
target_include_directories(mars_domain_binary_test PRIVATE
    ${CMAKE_BINARY_DIR}/cornerstone_includes
)