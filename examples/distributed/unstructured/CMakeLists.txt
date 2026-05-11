if(NOT MARS_ENABLE_FEM_EXAMPLES)
    return()
endif()

if(NOT MARS_ENABLE_UNSTRUCTURED OR NOT MARS_ENABLE_CUDA)
    message(STATUS "Skipping FEM examples (requires UNSTRUCTURED + CUDA)")
    return()
endif()

# Poisson example - simple user-facing template
add_executable(mars_ex1_poisson mars_ex1_poisson.cu)

# Link against mars library (public API only)
target_link_libraries(mars_ex1_poisson PRIVATE mars)

# Install
install(TARGETS mars_ex1_poisson RUNTIME DESTINATION bin/examples)

# Beam-tet mesh example - MFEM validation test
add_executable(mars_ex_beam_tet mars_ex_beam_tet.cu)
target_link_libraries(mars_ex_beam_tet PRIVATE mars)
install(TARGETS mars_ex_beam_tet RUNTIME DESTINATION bin/examples)

# Distributed GPU beam-tet example
add_executable(mars_ex_beam_tet_distributed mars_ex_beam_tet_distributed.cu)
target_link_libraries(mars_ex_beam_tet_distributed PRIVATE mars)
install(TARGETS mars_ex_beam_tet_distributed RUNTIME DESTINATION bin/examples)

# CVFEM example - full element-based assembly (27 NNZ/row, mathematically consistent)
add_executable(mars_cvfem_full mars_cvfem_full.cu)
target_link_libraries(mars_cvfem_full PRIVATE mars)
install(TARGETS mars_cvfem_full RUNTIME DESTINATION bin/examples)

# CVFEM example - graph-based assembly matching STK (7 NNZ/row, optimized sparsity)
add_executable(mars_cvfem_graph mars_cvfem_graph.cu)
if(TARGET CUDA::nvToolsExt)
    target_link_libraries(mars_cvfem_graph PRIVATE mars CUDA::nvToolsExt)
else()
    target_link_libraries(mars_cvfem_graph PRIVATE mars)
    target_compile_definitions(mars_cvfem_graph PRIVATE MARS_DISABLE_NVTX)
endif()
install(TARGETS mars_cvfem_graph RUNTIME DESTINATION bin/examples)

# CVFEM tet driver: graph-based assembly for linear tetrahedra
add_executable(mars_cvfem_graph_tet mars_cvfem_graph_tet.cu)
if(TARGET CUDA::nvToolsExt)
    target_link_libraries(mars_cvfem_graph_tet PRIVATE mars CUDA::nvToolsExt)
else()
    target_link_libraries(mars_cvfem_graph_tet PRIVATE mars)
    target_compile_definitions(mars_cvfem_graph_tet PRIVATE MARS_DISABLE_NVTX)
endif()
install(TARGETS mars_cvfem_graph_tet RUNTIME DESTINATION bin/examples)

# CVFEM Poisson solver
add_executable(mars_cvfem_poisson mars_cvfem_poisson.cu)
target_link_libraries(mars_cvfem_poisson PRIVATE mars)
install(TARGETS mars_cvfem_poisson RUNTIME DESTINATION bin/examples)

# AMR Poisson solver - adaptive mesh refinement with CVFEM
add_executable(mars_amr_poisson mars_amr_poisson.cu)
target_link_libraries(mars_amr_poisson PRIVATE mars)
target_include_directories(mars_amr_poisson PRIVATE
    ${CMAKE_SOURCE_DIR}/backend/distributed/unstructured/amr
)
install(TARGETS mars_amr_poisson RUNTIME DESTINATION bin/examples)

# AMR CVFEM Graph solver - adaptive mesh refinement with graph-based assembly
add_executable(mars_amr_cvfem_graph mars_amr_cvfem_graph.cu)
target_link_libraries(mars_amr_cvfem_graph PRIVATE mars)
target_include_directories(mars_amr_cvfem_graph PRIVATE
    ${CMAKE_SOURCE_DIR}/backend/distributed/unstructured/amr
)
install(TARGETS mars_amr_cvfem_graph RUNTIME DESTINATION bin/examples)

# AMR pressure-Poisson solver: same CVFEM operator + AMR loop as
# mars_amr_cvfem_graph, but with f = -div(v . grad v) for a prescribed
# free-stream + 2D tip-vortex velocity field. NASA wing-tip demo target.
add_executable(mars_amr_pressure_poisson mars_amr_pressure_poisson.cu)
target_link_libraries(mars_amr_pressure_poisson PRIVATE mars)
target_include_directories(mars_amr_pressure_poisson PRIVATE
    ${CMAKE_SOURCE_DIR}/backend/distributed/unstructured/amr
)
install(TARGETS mars_amr_pressure_poisson RUNTIME DESTINATION bin/examples)