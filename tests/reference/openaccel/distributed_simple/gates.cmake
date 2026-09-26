# Distributed SIMPLE gate targets; included by CMakeLists.txt (standalone) or, deferred, by
# ../distributed_matrix/inject.cmake inside a MARS configure. Paths are absolute on purpose.
get_filename_component(_dsimple_root "${CMAKE_CURRENT_LIST_DIR}/../../../.." ABSOLUTE)
set(_dsimple_gates "${CMAKE_CURRENT_LIST_DIR}")
set(_dsimple_segregated "${_dsimple_root}/backend/distributed/unstructured/fem/segregated")
if(NOT TARGET MPI::MPI_CXX)
    find_package(MPI REQUIRED COMPONENTS CXX)
endif()

add_executable(mars_distributed_simple_host_gate "${_dsimple_gates}/simple_gate.cpp")
target_include_directories(mars_distributed_simple_host_gate PRIVATE "${_dsimple_segregated}" "${_dsimple_gates}")
target_compile_features(mars_distributed_simple_host_gate PRIVATE cxx_std_20)
target_compile_definitions(mars_distributed_simple_host_gate PRIVATE OMPI_SKIP_MPICXX=1 MPICH_SKIP_MPICXX=1)
target_link_libraries(mars_distributed_simple_host_gate PRIVATE MPI::MPI_CXX)

# One-rank SimpleRunner references, then distributed comparisons (a stranded rank hangs -> timeout).
function(_dsimple_test name ranks)
    add_test(NAME ${name} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
             $<TARGET_FILE:mars_distributed_simple_host_gate> ${MPIEXEC_POSTFLAGS} ${ARGN})
    set_tests_properties(${name} PROPERTIES TIMEOUT 900 PROCESSORS ${ranks})
endfunction()
set(_ref "${CMAKE_CURRENT_BINARY_DIR}/dsimple")
_dsimple_test(marsDistributedSimpleReference 1 --write-reference ${_ref}-2it.bin)
_dsimple_test(marsDistributedSimpleReferenceShear 1 --write-reference ${_ref}-shear.bin --backflow -1 --iterations 12)
_dsimple_test(marsDistributedSimpleReferenceClosed 1 --write-reference ${_ref}-closed.bin --backflow 1 --iterations 8)
_dsimple_test(marsDistributedSimpleReferenceConverge 1 --write-reference ${_ref}-conv.bin --mesh 8x2x2 --converge 4000)
set_tests_properties(marsDistributedSimpleReference marsDistributedSimpleReferenceShear marsDistributedSimpleReferenceClosed
    marsDistributedSimpleReferenceConverge PROPERTIES FIXTURES_SETUP dsimple_references)
set(_compare)
foreach(_ranks 1 2 4)
    _dsimple_test(marsDistributedSimple2it_${_ranks} ${_ranks} --reference ${_ref}-2it.bin)
    _dsimple_test(marsDistributedSimpleShear_${_ranks} ${_ranks} --reference ${_ref}-shear.bin --backflow -1 --iterations 12)
    _dsimple_test(marsDistributedSimpleClosed_${_ranks} ${_ranks} --reference ${_ref}-closed.bin --backflow 1 --iterations 8)
    _dsimple_test(marsDistributedSimpleConverge_${_ranks} ${_ranks} --reference ${_ref}-conv.bin --mesh 8x2x2 --converge 4000)
    list(APPEND _compare marsDistributedSimple2it_${_ranks} marsDistributedSimpleShear_${_ranks}
                         marsDistributedSimpleClosed_${_ranks} marsDistributedSimpleConverge_${_ranks})
endforeach()
foreach(_ranks 2 4)
    foreach(_fault stale-ghost missing-element duplicate-face)
        _dsimple_test(marsDistributedSimpleFault_${_fault}_${_ranks} ${_ranks} --reference ${_ref}-2it.bin --fault ${_fault})
        list(APPEND _compare marsDistributedSimpleFault_${_fault}_${_ranks})
    endforeach()
endforeach()
set_tests_properties(${_compare} PROPERTIES FIXTURES_REQUIRED dsimple_references)

if(TARGET mars AND MARS_ENABLE_CUDA AND MARS_ENABLE_HYPRE)
    # Mirrors mars_segregated_simple: link mars, MARS_REPLAY_CUDA, CUDA C++20.
    add_executable(mars_distributed_simple_cuda_gate "${_dsimple_gates}/simple_gate.cu")
    target_link_libraries(mars_distributed_simple_cuda_gate PRIVATE mars)
    target_compile_definitions(mars_distributed_simple_cuda_gate PRIVATE MARS_REPLAY_CUDA)
    target_include_directories(mars_distributed_simple_cuda_gate PRIVATE "${_dsimple_segregated}" "${_dsimple_gates}")
    set_target_properties(mars_distributed_simple_cuda_gate PROPERTIES CUDA_STANDARD 20 CUDA_STANDARD_REQUIRED ON)
endif()
