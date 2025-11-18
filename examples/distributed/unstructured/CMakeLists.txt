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