if(MARS_ENABLE_CXXOPTS)
    # ##########################################################################
    find_package(cxxopts QUIET)

    if(NOT cxxopts_FOUND)
        include(FetchContent)

        FetchContent_Declare(
            cxxopts
            GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
            GIT_TAG v3.0.0)

        FetchContent_GetProperties(cxxopts)

        if(NOT cxxopts_POPULATED)
            FetchContent_Populate(cxxopts)
            add_subdirectory(${cxxopts_SOURCE_DIR} ${cxxopts_BINARY_DIR})
        endif()

        set_target_properties(cxxopts PROPERTIES FOLDER extern)

        set(MARS_CXXOPTS_LIBRARIES cxxopts)

        message ("Setup xcxops")
    else()
        message("Found cxxopts")
        set(MARS_CXXOPTS_LIBRARIES cxxopts::cxxopts)
    endif()
    set(WITH_CXXOPTS ON)

    # ##########################################################################

    message(${cxxopts_SOURCE_DIR})
    # target_include_directories(mars PUBLIC ${cxxopts_SOURCE_DIR}/include)
endif(MARS_ENABLE_CXXOPTS)