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

    # cxxopts is header-only CLI tooling used only by in-tree example drivers, not by the
    # mars library or its public API. Link it BUILD-only so examples get it transitively,
    # while it never leaks into the installed Mars::mars interface (cxxopts ships no
    # library, so a leaked reference becomes a broken -lcxxopts for consumers).
    list(APPEND MARS_DEP_LIBRARIES $<BUILD_INTERFACE:cxxopts::cxxopts>)
    list(APPEND MARS_DEP_INCLUDES ${cxxopts_INCLUDE_DIRS})

    set(MARS_ENABLE_CXXOPTS ON)
endif()