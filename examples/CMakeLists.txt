if(MARS_ENABLE_EXAMPLES)
    list(APPEND EXAMPLES_MODULES 
        .
        util
        unstructured
        distributed
        distributed/staggered)


    scan_directories(${CMAKE_CURRENT_SOURCE_DIR} "${EXAMPLES_MODULES}" MARS_BUILD_INCLUDES MARS_HEADERS MARS_SOURCES)

    set(MARS_BUILD_INCLUDES
        ${MARS_BUILD_INCLUDES}
        PARENT_SCOPE)

    set(MARS_HEADERS
        ${MARS_HEADERS}
        PARENT_SCOPE)

    set(MARS_SOURCES
        ${MARS_SOURCES}
        PARENT_SCOPE)


    set(MARS_TARGETS mars_examples PARENT_SCOPE)

    add_executable(mars_examples mars_discretization.cpp)
    target_link_libraries(mars_examples PRIVATE mars)
    
    # Build MFEM to MARS binary converter tool
    add_executable(mfem_to_mars_binary 
        distributed/unstructured/mfem_to_mars_binary.cpp)
    target_include_directories(mfem_to_mars_binary PRIVATE
        ${CMAKE_SOURCE_DIR}/backend/distributed/unstructured/utils)
endif()
