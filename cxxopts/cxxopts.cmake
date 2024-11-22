if(MARS_ENABLE_CXXOPTS)
    find_package(cxxopts QUIET)

    if(NOT cxxopts_FOUND)
        message(STATUS "cxxopts not found, fetching from GitHub.")
        include(FetchContent)
        FetchContent_Declare(
            cxxopts
            GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
            GIT_TAG v3.0.0
        )
        FetchContent_MakeAvailable(cxxopts)
    else()
        message(STATUS "Found cxxopts")
    endif()

    # Add cxxopts to the list of dependencies
    list(APPEND MARS_DEP_LIBRARIES cxxopts::cxxopts)
    list(APPEND MARS_DEP_INCLUDES ${cxxopts_INCLUDE_DIRS})

    set(MARS_ENABLE_CXXOPTS ON)
endif()