add_subdirectory(communication)

if(MARS_ENABLE_KOKKOS)
    add_subdirectory(sfc)
    add_subdirectory(discretization)
    add_subdirectory(mesh_manager)
endif()

if(MARS_ENABLE_UNSTRUCTURED)
    add_subdirectory(unstructured)
endif()
