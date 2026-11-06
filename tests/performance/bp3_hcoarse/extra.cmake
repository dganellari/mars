function(bp3_h_targets)
    if(NOT TARGET mars_bp3_coarse)
        message(FATAL_ERROR "The h-level trial needs the accepted CUDA/Hypre adapter")
    endif()
    add_executable(bp3_h_solver "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/private/mars_marsir_bp.cu")
    target_sources(bp3_h_solver PRIVATE $<TARGET_OBJECTS:mars_bp3_coarse>)
    target_compile_definitions(bp3_h_solver PRIVATE MARS_BP3_COARSE_CONSOLIDATION MARS_BP3_HCOARSE_CHECK)
    target_link_libraries(bp3_h_solver PRIVATE mars CUDA::cuda_driver)
    target_link_options(bp3_h_solver PRIVATE -Wl,--export-dynamic)
    set_target_properties(bp3_h_solver PROPERTIES CUDA_STANDARD 20 CUDA_STANDARD_REQUIRED ON
                          RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/hcoarse")
endfunction()

cmake_language(DEFER CALL bp3_h_targets)
