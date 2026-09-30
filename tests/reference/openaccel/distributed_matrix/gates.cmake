# Gate targets; included by CMakeLists.txt (standalone) or, deferred, by inject.cmake inside
# a MARS configure. Paths are absolute because the including directory differs.
get_filename_component(_dmatrix_root "${CMAKE_CURRENT_LIST_DIR}/../../../.." ABSOLUTE)
set(_dmatrix_gates "${CMAKE_CURRENT_LIST_DIR}")
set(_dmatrix_segregated "${_dmatrix_root}/backend/distributed/unstructured/fem/segregated")
if(NOT TARGET MPI::MPI_CXX)
    find_package(MPI REQUIRED COMPONENTS CXX)
endif()

add_executable(mars_distributed_matrix_host_gate "${_dmatrix_gates}/host_gate.cpp")
target_include_directories(mars_distributed_matrix_host_gate PRIVATE "${_dmatrix_segregated}" "${_dmatrix_gates}")
target_compile_features(mars_distributed_matrix_host_gate PRIVATE cxx_std_20)
target_compile_definitions(mars_distributed_matrix_host_gate PRIVATE OMPI_SKIP_MPICXX=1 MPICH_SKIP_MPICXX=1)
target_link_libraries(mars_distributed_matrix_host_gate PRIVATE MPI::MPI_CXX)
foreach(_ranks 1 2 4)
    add_test(NAME marsDistributedMatrixHost${_ranks}
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${_ranks} ${MPIEXEC_PREFLAGS}
                     $<TARGET_FILE:mars_distributed_matrix_host_gate> ${MPIEXEC_POSTFLAGS})
    # A rank stranded outside a collective hangs; the timeout turns that into a failure.
    set_tests_properties(marsDistributedMatrixHost${_ranks} PROPERTIES TIMEOUT 120 PROCESSORS ${_ranks})
endforeach()

if(TARGET mars AND MARS_ENABLE_CUDA AND MARS_ENABLE_HYPRE)
    # Mirrors mars_segregated_simple_check: link mars, CUDA C++20; CMAKE_CUDA_FLAGS carry --extended-lambda.
    add_executable(mars_distributed_matrix_cuda_gate "${_dmatrix_gates}/cuda_gate.cu")
    target_link_libraries(mars_distributed_matrix_cuda_gate PRIVATE mars)
    target_include_directories(mars_distributed_matrix_cuda_gate PRIVATE "${_dmatrix_segregated}" "${_dmatrix_gates}")
    set_target_properties(mars_distributed_matrix_cuda_gate PROPERTIES CUDA_STANDARD 20 CUDA_STANDARD_REQUIRED ON)
endif()
