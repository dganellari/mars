# The private callbacks must remain replaceable in the linked Hypre library.
find_path(_mars_bp3_hypre_include HYPRE_config.h PATHS ${HYPRE_INCLUDE_DIRS} NO_DEFAULT_PATH NO_CACHE REQUIRED)
file(READ "${_mars_bp3_hypre_include}/HYPRE_config.h" _mars_bp3_hypre_config)
if(NOT _mars_bp3_hypre_config MATCHES "#define HYPRE_RELEASE_NUMBER[ \t]+23300([ \r\n]|$)"
   OR NOT _mars_bp3_hypre_config MATCHES "#define HYPRE_USING_CUDA([ \r\n]|$)")
    message(FATAL_ERROR "BP3 coarse consolidation requires CUDA-enabled Hypre 2.33.0 headers")
endif()

set(_mars_bp3_hypre_candidates ${HYPRE_LIBRARIES})
list(FILTER _mars_bp3_hypre_candidates INCLUDE REGEX "/lib(HYPRE|hypre)[.-]")
list(LENGTH _mars_bp3_hypre_candidates _mars_bp3_hypre_count)
if(NOT _mars_bp3_hypre_count EQUAL 1)
    message(FATAL_ERROR "BP3 coarse consolidation requires one resolved Hypre library path")
endif()
list(GET _mars_bp3_hypre_candidates 0 _mars_bp3_hypre_library)
if(NOT _mars_bp3_hypre_library MATCHES "\\.so(\\.|$)")
    message(FATAL_ERROR "BP3 coarse consolidation needs shared Hypre with replaceable callbacks")
endif()
find_program(_mars_bp3_readelf NAMES readelf llvm-readelf REQUIRED)
execute_process(COMMAND "${_mars_bp3_readelf}" -rW "${_mars_bp3_hypre_library}"
    RESULT_VARIABLE _mars_bp3_readelf_status OUTPUT_VARIABLE _mars_bp3_relocations
    ERROR_VARIABLE _mars_bp3_readelf_error)
if(NOT _mars_bp3_readelf_status EQUAL 0)
    message(FATAL_ERROR "Cannot inspect Hypre callback relocations: ${_mars_bp3_readelf_error}")
endif()
foreach(_mars_bp3_callback hypre_seqAMGSetup hypre_seqAMGCycle)
    if(NOT _mars_bp3_relocations MATCHES "[^\n]*(JUMP_SLOT|GLOB_DAT)[^\n]*${_mars_bp3_callback}([ \n]|$)")
        message(FATAL_ERROR "Hypre does not expose a replaceable ${_mars_bp3_callback} callback")
    endif()
endforeach()

add_library(mars_bp3_coarse OBJECT
    ${CMAKE_SOURCE_DIR}/backend/distributed/unstructured/solvers/mars_cellwise_amg_consolidate.cu)
target_include_directories(mars_bp3_coarse PRIVATE ${HYPRE_INCLUDE_DIRS})
target_compile_definitions(mars_bp3_coarse PRIVATE BP3_CONSOLIDATE_CUDA)
target_link_libraries(mars_bp3_coarse PRIVATE MPI::MPI_CXX CUDA::cudart)
set_target_properties(mars_bp3_coarse PROPERTIES CUDA_STANDARD 20 CUDA_STANDARD_REQUIRED ON)
target_compile_options(mars_bp3_coarse PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-Wall,-Wextra,-Werror>)

# Only this executable opts into the adapter through the shared P1Amg header.
target_sources(mars_marsir_bp PRIVATE $<TARGET_OBJECTS:mars_bp3_coarse>)
target_compile_definitions(mars_marsir_bp PRIVATE MARS_BP3_COARSE_CONSOLIDATION)
target_link_options(mars_marsir_bp PRIVATE -Wl,--export-dynamic)

if(MARS_ENABLE_TESTS)
    add_executable(mars_bp3_coarse_test ${CMAKE_SOURCE_DIR}/tests/bp3_coarse/fixture.cu)
    target_include_directories(mars_bp3_coarse_test PRIVATE ${CMAKE_SOURCE_DIR} ${HYPRE_INCLUDE_DIRS})
    target_compile_definitions(mars_bp3_coarse_test PRIVATE BP3_CONSOLIDATE_CUDA)
    target_link_libraries(mars_bp3_coarse_test PRIVATE MPI::MPI_CXX "${_mars_bp3_hypre_library}" CUDA::cudart)
    target_link_options(mars_bp3_coarse_test PRIVATE -Wl,--export-dynamic)
    set_target_properties(mars_bp3_coarse_test PROPERTIES CUDA_STANDARD 20 CUDA_STANDARD_REQUIRED ON)
    target_compile_options(mars_bp3_coarse_test PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-Wall,-Wextra,-Werror>)
    include(${CMAKE_SOURCE_DIR}/tests/bp3_coarse/register.cmake)
endif()
