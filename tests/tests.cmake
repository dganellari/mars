if(MARS_ENABLE_TESTING)

    set(MARS_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests)

    list(APPEND TEST_MODULES . /adios2)

    find_project_files(${MARS_TEST_DIR} "." UNIT_TESTS_HEADERS
                       UNIT_TESTS_SOURCES)

    # ##########################################################################

    enable_testing()

    find_package(GTest QUIET)

    if(NOT GTest_FOUND)
        include(FetchContent)

        FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG release-1.10.0)

        FetchContent_GetProperties(googletest)

        if(NOT googletest_POPULATED)
            FetchContent_Populate(googletest)
            add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})
        endif()

        set_target_properties(gtest PROPERTIES FOLDER extern)
        set_target_properties(gtest_main PROPERTIES FOLDER extern)
        set_target_properties(gmock PROPERTIES FOLDER extern)
        set_target_properties(gmock_main PROPERTIES FOLDER extern)

        set(MARS_G_TEST_LIBRARIES gtest gmock)
        # add_library(googletest ALIAS gtest)
    else()
        set(MARS_G_TEST_LIBRARIES GTest::GTest)
        # add_library(googletest PUBLIC GTest::GTest)
    endif()

    include(GoogleTest)

    macro(mars_add_test TESTNAME TEST_FILE TEST_SOURCES)
        add_executable(${TESTNAME} ${TEST_FILE} ${ARGN})
        target_link_libraries(${TESTNAME} mars
                              ${MARS_G_TEST_LIBRARIES})

        gtest_add_tests(TARGET ${TESTNAME} SOURCES ${TEST_SOURCES})
        set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)

        # Does not work on cluster env where you cannot run MPI on master node
        # So for the moment it will only work on apple laptops (the test can be
        # run with ./MARS_test anyway but not with make test)
        if(APPLE
           OR WIN32
           OR MARS_ALLOW_DISCOVER_TESTS)
            gtest_discover_tests(
                ${TESTNAME}
                WORKING_DIRECTORY ${PROJECT_DIR}
                PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_DIR}")
        endif()

    endmacro()

    # ##########################################################################

    mars_add_test(
        mars_test ${CMAKE_CURRENT_SOURCE_DIR}/tests/test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/mars_MeshGenerationTest.cpp)
endif(MARS_ENABLE_TESTING)