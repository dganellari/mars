cmake_minimum_required(VERSION 3.24)
project(mars_hex_flow_checks LANGUAGES CXX)
enable_testing()
find_package(Mars REQUIRED)
find_package(MPI REQUIRED COMPONENTS CXX)
find_package(Python3 REQUIRED COMPONENTS Interpreter)
if(NOT TARGET Mars::mars_unstructured OR NOT MARS_GPU_LANGUAGE MATCHES "^(CUDA|HIP)$")
    message(FATAL_ERROR "Install the candidate GPU/Hypre unstructured MARS library first")
endif()
if(NOT DEFINED CMAKE_${MARS_GPU_LANGUAGE}_ARCHITECTURES)
    set(CMAKE_${MARS_GPU_LANGUAGE}_ARCHITECTURES "${MARS_GPU_ARCHITECTURES}")
endif()
set(CMAKE_${MARS_GPU_LANGUAGE}_STANDARD 20)
enable_language(${MARS_GPU_LANGUAGE})
set_source_files_properties(gpu.cu PROPERTIES LANGUAGE ${MARS_GPU_LANGUAGE})
add_executable(hex_flow_gpu gpu.cu)
target_link_libraries(hex_flow_gpu PRIVATE Mars::mars Mars::mars_unstructured MPI::MPI_CXX)
if(MARS_GPU_LANGUAGE STREQUAL "CUDA")
    set_target_properties(hex_flow_gpu PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()
foreach(ranks IN ITEMS 1 2 4)
    foreach(mode IN ITEMS flow failed invalid)
        add_test(NAME hex_${mode}_${ranks}
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:hex_flow_gpu> ${MPIEXEC_POSTFLAGS} ${mode})
        set_tests_properties(hex_${mode}_${ranks} PROPERTIES TIMEOUT 300 PROCESSORS ${ranks})
    endforeach()
endforeach()
# The same kinetic energy on 1 and on 2 / 4 ranks (Hypre's AMG differs with the partition: the
# pressure solves agree to their tolerance, hence rtol 1e-6).
string(REPLACE ";" " " _preflags "${MPIEXEC_PREFLAGS}")
foreach(ranks IN ITEMS 2 4)
    add_test(NAME hex_flow_ranks_${ranks}
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../release/check_rank_value.py --np ${ranks}
                "--regex=KE=([-+0-9.eE]+)" "--mpiexec=${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
                "--preflags=${_preflags}" --rtol 1e-6 -- $<TARGET_FILE:hex_flow_gpu> flow)
    set_tests_properties(hex_flow_ranks_${ranks} PROPERTIES TIMEOUT 600 PROCESSORS ${ranks})
endforeach()
