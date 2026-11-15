# First, define regular modules (without test directories)
list(APPEND DISTRIBUTED_MODULES
    .
    communication
    communication/util
    mars_env
    # Removed testsuite directories from here
)

# Define test modules separately
list(APPEND DISTRIBUTED_TEST_MODULES
    testsuite/communication
)

if(MARS_ENABLE_KOKKOS)
    if(NOT MARS_ENABLE_MPI)
        message(FATAL_ERROR "The structured Kokkos backend requires MARS_ENABLE_MPI=ON.")
    endif()
    list(APPEND DISTRIBUTED_MODULES
        structured/sfc
        structured/mesh
        structured/mesh_manager
        # Removed testsuite directories from here
    )
    
    # Add Kokkos test directories to test modules
    list(APPEND DISTRIBUTED_TEST_MODULES
        testsuite/sfc
        testsuite/discretization
        testsuite/mesh_manager
    )
endif()

# Add unstructured if enabled
if(MARS_ENABLE_UNSTRUCTURED)
    # Add directories to search paths
    list(APPEND DISTRIBUTED_MODULES
        unstructured
        unstructured/utils
    )
    list(APPEND DISTRIBUTED_TEST_MODULES testsuite/unstructured)
    
    # Add the unstructured subdirectory
    add_subdirectory(unstructured)
endif()

# Export variables to parent scope
set(MARS_DEP_INCLUDES ${MARS_DEP_INCLUDES} PARENT_SCOPE)
set(MARS_DEP_LIBRARIES ${MARS_DEP_LIBRARIES} PARENT_SCOPE)
set(MARS_TARGETS ${MARS_TARGETS} PARENT_SCOPE)

# Scan directories for library source files (excluding test directories)
scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "${DISTRIBUTED_MODULES}" MARS_BUILD_INCLUDES MARS_HEADERS MARS_SOURCES)

# Export build variables to parent scope
set(MARS_BUILD_INCLUDES ${MARS_BUILD_INCLUDES} PARENT_SCOPE)
set(MARS_HEADERS ${MARS_HEADERS} PARENT_SCOPE)
set(MARS_SOURCES ${MARS_SOURCES} PARENT_SCOPE)

# Handle test directories separately, only if tests are enabled
if(MARS_ENABLE_TESTS)
    if(MARS_ENABLE_KOKKOS)
        include(testsuite/StructuredTests.cmake)
    endif()
    # Process test modules but don't add them to MARS_SOURCES
    foreach(TEST_DIR ${DISTRIBUTED_TEST_MODULES})
        # Add the directory to the build
        add_subdirectory(${TEST_DIR})
    endforeach()
endif()
