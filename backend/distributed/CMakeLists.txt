list(APPEND DISTRIBUTED_MODULES
    .
    communication
    communication/util
    mesh
    mesh_manager
    mars_env
    sfc
    testsuite
    testsuite/common
    )

if(MARS_ENABLE_UNSTRUCTURED)
    list(APPEND DISTRIBUTED_MODULES
        octree
        unstructured
        unstructured/utils
        )
endif()

scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "${DISTRIBUTED_MODULES}" MARS_BUILD_INCLUDES MARS_HEADERS MARS_SOURCES)

set(MARS_BUILD_INCLUDES
    ${MARS_BUILD_INCLUDES}
    PARENT_SCOPE)

set(MARS_HEADERS
    ${MARS_HEADERS}
    PARENT_SCOPE)

set(MARS_SOURCES
    ${MARS_SOURCES}
    PARENT_SCOPE)
