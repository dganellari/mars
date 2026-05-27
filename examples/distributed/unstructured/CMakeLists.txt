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

# AMR + CVFEM Poisson for linear tetrahedra (full 4x4 assembly + CG)
add_executable(mars_amr_cvfem_graph_tet mars_amr_cvfem_graph_tet.cu)
target_link_libraries(mars_amr_cvfem_graph_tet PRIVATE mars)
target_include_directories(mars_amr_cvfem_graph_tet PRIVATE
    ${CMAKE_SOURCE_DIR}/backend/distributed/unstructured/amr
)
install(TARGETS mars_amr_cvfem_graph_tet RUNTIME DESTINATION bin/examples)

# AMR pressure-Poisson solver: same CVFEM operator + AMR loop as
# mars_amr_cvfem_graph, but with f = -div(v . grad v) for a prescribed
# free-stream + 2D tip-vortex velocity field. NASA wing-tip demo target.
# When MARS_ENABLE_HYPRE=ON, --solver=hypre selects PCG+BoomerAMG; the
# MARS_ENABLE_HYPRE define and Hypre link flags arrive transitively through
# the mars target (mars_config.hpp + MARS_DEP_LIBRARIES), same as
# mars_ex_beam_tet_distributed -- no extra CMake plumbing needed here.
add_executable(mars_amr_pressure_poisson mars_amr_pressure_poisson.cu)
target_link_libraries(mars_amr_pressure_poisson PRIVATE mars)
target_include_directories(mars_amr_pressure_poisson PRIVATE
    ${CMAKE_SOURCE_DIR}/backend/distributed/unstructured/amr
)
install(TARGETS mars_amr_pressure_poisson RUNTIME DESTINATION bin/examples)

# Advection-diffusion validator: fixed cube mesh, backward Euler, Gaussian IC
# transported by a prescribed solid-body rotation. Phase B.2 of the
# NS-projection plan -- reuses CvfemHexAssembler::assembleFull with gamma=nu
# and nonzero mdot to verify the advection branch before building a full NS
# driver. omega=0 falls back to the Phase B.1 heat-equation regression.
# Picks up Hypre transitively through the mars target.
add_executable(mars_amr_advdiff mars_amr_advdiff.cu)
target_link_libraries(mars_amr_advdiff PRIVATE mars)
target_include_directories(mars_amr_advdiff PRIVATE
    ${CMAKE_SOURCE_DIR}/backend/distributed/unstructured/amr
)
install(TARGETS mars_amr_advdiff RUNTIME DESTINATION bin/examples)

# Divergence operator validator: fixed cube mesh, prescribed velocity field,
# compares CVFEM SCS-face discrete div(u) against analytical references
# (rotation -> 0, sinx -> pi*cos(pi*x)). Phase B.3 of the NS-projection plan
# -- validates the RHS of the pressure Poisson equation before wiring it into
# a full projection driver.
add_executable(mars_amr_div mars_amr_div.cu)
target_link_libraries(mars_amr_div PRIVATE mars)
target_include_directories(mars_amr_div PRIVATE
    ${CMAKE_SOURCE_DIR}/backend/distributed/unstructured/amr
)
install(TARGETS mars_amr_div RUNTIME DESTINATION bin/examples)

# Gradient operator validator: fixed cube mesh, prescribed scalar p field,
# compares CVFEM SCS-face discrete grad(p) against analytical references
# (linear -> (a,b,c), sin -> pi-products, constant -> 0). Phase B.4 of the
# NS-projection plan -- dual of B.3; validates the pressure-correction
# operator u^{n+1} = u* - dt * grad p.
add_executable(mars_amr_grad mars_amr_grad.cu)
target_link_libraries(mars_amr_grad PRIVATE mars)
target_include_directories(mars_amr_grad PRIVATE
    ${CMAKE_SOURCE_DIR}/backend/distributed/unstructured/amr
)
install(TARGETS mars_amr_grad RUNTIME DESTINATION bin/examples)

# Galerkin Laplacian A = D M^{-1} D^T standalone validator (Phase E, Stage 1).
# Matrix-free SpMV with the same per-node + reverse-halo D / D^T scatter as the
# NS driver. Runs five correctness tests (linear, symmetry, nullspace, psd, cg)
# in isolation so the operator math can be locked down before any NS wiring.
add_executable(mars_amr_ddt mars_amr_ddt.cu)
target_link_libraries(mars_amr_ddt PRIVATE mars)
target_include_directories(mars_amr_ddt PRIVATE
    ${CMAKE_SOURCE_DIR}/backend/distributed/unstructured/amr
)
install(TARGETS mars_amr_ddt RUNTIME DESTINATION bin/examples)

# Chorin projection NS driver (Phase C): ties B.1-B.4 together into an
# unsteady incompressible Navier-Stokes solver. Per step: explicit predictor
# (advection + grad p^n), implicit diffusion (3 CG/Hypre solves), pressure
# Poisson (1 solve), corrector (subtract dt/rho * grad phi), pressure update.
# Test case: lid-driven cavity / channel flow. The solver kernels and
# NSStepper struct live in backend/.../fem/mars_ns_solver.hpp so other
# drivers (TGV, future flow-past-sphere) can reuse them. Picks up Hypre
# transitively through mars.
add_executable(mars_amr_ns_projection mars_amr_ns_projection.cu)
target_link_libraries(mars_amr_ns_projection PRIVATE mars)
target_include_directories(mars_amr_ns_projection PRIVATE
    ${CMAKE_SOURCE_DIR}/backend/distributed/unstructured/amr
)
install(TARGETS mars_amr_ns_projection RUNTIME DESTINATION bin/examples)

# Taylor-Green vortex driver. Periodic-cube NS solver + AMR. Currently a
# skeleton that loads a mesh, builds the periodic node-pair map, applies the
# TGV initial condition, and writes the IC as VTU. The full projection time
# loop will be wired against mars_ns_solver.hpp in a follow-up commit once
# BCKind::Periodic is added to NSStepper.
add_executable(mars_tgv mars_tgv.cu)
target_link_libraries(mars_tgv PRIVATE mars)
target_include_directories(mars_tgv PRIVATE
    ${CMAKE_SOURCE_DIR}/backend/distributed/unstructured/amr
)
install(TARGETS mars_tgv RUNTIME DESTINATION bin/examples)