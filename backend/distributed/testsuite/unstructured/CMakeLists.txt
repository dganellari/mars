include_directories(${CMAKE_BINARY_DIR}/_deps/cornerstone_fetch-src/include)

# Find cornerstone GPU library FIRST
execute_process(
    COMMAND find ${CMAKE_BINARY_DIR} -name "libcstone_gpu.a" -o -name "libcornerstone_gpu.a"
    OUTPUT_VARIABLE CSTONE_GPU_LIB
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

message(STATUS "Found cornerstone GPU library: ${CSTONE_GPU_LIB}")

# Create tests
mars_add_test(mars_domain_binary_test ${CMAKE_CURRENT_SOURCE_DIR}/mars_domain_binary_mesh.cu)
mars_add_test(mars_read_mesh_binary_test ${CMAKE_CURRENT_SOURCE_DIR}/mars_read_mesh_binary.cpp)
# mars_add_test(mars_read_mesh_adios2_test ${CMAKE_CURRENT_SOURCE_DIR}/mars_read_mesh_adios2.cpp)

# Link mars_unstructured after GPU lib (WITHOUT keywords)
target_link_libraries(mars_domain_binary_test mars_unstructured)
target_link_libraries(mars_read_mesh_binary_test mars_unstructured)

# Link cornerstone GPU first, then mars_unstructured (WITHOUT keywords)
if(CSTONE_GPU_LIB)
    target_link_libraries(mars_domain_binary_test ${CSTONE_GPU_LIB})
    target_link_libraries(mars_read_mesh_binary_test ${CSTONE_GPU_LIB})
endif()
