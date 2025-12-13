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

# CVFEM example - face-based control volume finite element method
add_executable(mars_cvfem_example mars_cvfem_example.cu)
target_link_libraries(mars_cvfem_example PRIVATE mars)
install(TARGETS mars_cvfem_example RUNTIME DESTINATION bin/examples)