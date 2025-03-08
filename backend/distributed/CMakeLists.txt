list(APPEND DISTRIBUTED_MODULES
    .
    communication
    communication/util
    mars_env
    testsuite/communication
    testsuite/unstructured
    )

if(MARS_ENABLE_KOKKOS)
    list(APPEND DISTRIBUTED_MODULES
        sfc
        mesh
        mesh_manager
        testsuite/sfc
        testsuite/common
        testsuite/discretization
        testsuite/mesh_manager
        )
endif()

# Add unstructured if enabled
if(MARS_ENABLE_UNSTRUCTURED)
    # Add directories to search paths
    list(APPEND DISTRIBUTED_MODULES
        octree
        unstructured/utils
        )
    
    # Add the unstructured subdirectory
    add_subdirectory(unstructured)
    
    # Collect information to be used by parent CMake
    if(TARGET mars_unstructured)
        list(APPEND MARS_DEP_LIBRARIES mars_unstructured)
    endif()
endif()

# Export variables to parent scope
set(MARS_DEP_INCLUDES ${MARS_DEP_INCLUDES} PARENT_SCOPE)
set(MARS_DEP_LIBRARIES ${MARS_DEP_LIBRARIES} PARENT_SCOPE)
set(MARS_TARGETS ${MARS_TARGETS} PARENT_SCOPE)

# Scan directories for source files
scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "${DISTRIBUTED_MODULES}" MARS_BUILD_INCLUDES MARS_HEADERS MARS_SOURCES)

# Export build variables to parent scope
set(MARS_BUILD_INCLUDES
    ${MARS_BUILD_INCLUDES}
    PARENT_SCOPE)

set(MARS_HEADERS
    ${MARS_HEADERS}
    PARENT_SCOPE)

set(MARS_SOURCES
    ${MARS_SOURCES}
    PARENT_SCOPE)
