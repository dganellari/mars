if(MARS_ENABLE_CXXOPTS)
    # ##########################################################################
    find_package(cxxopts QUIET)

    if(NOT cxxopts_FOUND)
        message("cxxopts not found, now installing it.")
        include(FetchContent)

        FetchContent_Declare(
            cxxopts
            GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
            GIT_TAG v3.0.0)
        FetchContent_MakeAvailable(cxxopts)
        set(MARS_DEP_LIBRARIES "${MARS_DEP_LIBRARIES};cxxopts::cxxopts")
        # set(MARS_DEP_INCLUDES "${MARS_DEP_INCLUDES};cxxopts::cxxopts")

        # message(${MARS_DEP_INCLUDES})

        # FetchContent_GetProperties(cxxopts)

        # if(NOT cxxopts_POPULATED)
        #     FetchContent_Populate(cxxopts)
        #     add_subdirectory(${cxxopts_SOURCE_DIR} ${cxxopts_BINARY_DIR})
        # endif()

        # set_target_properties(cxxopts PROPERTIES FOLDER extern)

        # set(MARS_CXXOPTS_LIBRARIES cxxopts::cxxopts)

    else()
        message("Found cxxopts")
        set(MARS_CXXOPTS_LIBRARIES cxxopts::cxxopts)
    endif()
    set(WITH_CXXOPTS ON)

    # ##########################################################################

    # message("MARS_CXXOPTS_LIBRARIES: ${MARS_CXXOPTS_LIBRARIES}")

    # if(TARGET ${MARS_CXXOPTS_LIBRARIES})
    #     message("MARS_CXXOPTS_LIBRARIES: ${MARS_CXXOPTS_LIBRARIES}")
    #     # get_target_property(MARS_CXXOPTS_INCLUDES ${MARS_CXXOPTS_LIBRARIES}
    #                   # INTERFACE_INCLUDE_DIRECTORIES)
    #     set(MARS_DEP_INCLUDES "${MARS_DEP_LIBRARIES};${MARS_CXXOPTS_LIBRARIES}")
    # endif()

    # set(MARS_EXEC_DEP_LIBRARIES "${MARS_EXEC_DEP_LIBRARIES};${MARS_CXXOPTS_LIBRARIES}")


endif(MARS_ENABLE_CXXOPTS)