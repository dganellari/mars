if(MARS_ENABLE_TESTS)

    enable_testing()

    find_package(GTest QUIET)

    if(NOT GTest_FOUND)
        include(FetchContent)

        FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG release-1.12.0)

        FetchContent_GetProperties(googletest)

        if(NOT googletest_POPULATED)
            FetchContent_MakeAvailable(googletest)
        endif()

        # Add include directories for downloaded GTest
        include_directories(${googletest_SOURCE_DIR}/googletest/include)
        include_directories(${googletest_SOURCE_DIR}/googlemock/include)

        set_target_properties(gtest PROPERTIES FOLDER extern)
        set_target_properties(gtest_main PROPERTIES FOLDER extern)
        set_target_properties(gmock PROPERTIES FOLDER extern)
        set_target_properties(gmock_main PROPERTIES FOLDER extern)

        set(MARS_G_TEST_LIBRARIES gtest gmock)
        # add_library(googletest ALIAS gtest)
    else()
        # Add include directories for system GTest
        include_directories(${GTEST_INCLUDE_DIRS})
        set(MARS_G_TEST_LIBRARIES GTest::GTest)
        # add_library(googletest PUBLIC GTest::GTest)
    endif()

    include(GoogleTest)

    macro(mars_add_test TESTNAME TEST_FILE)
        add_executable(${TESTNAME} ${TEST_FILE} ${ARGN})
        target_link_libraries(${TESTNAME} mars ${MARS_G_TEST_LIBRARIES})

        # Automatically include TEST_FILE in the TEST_SOURCES
        gtest_add_tests(TARGET ${TESTNAME} SOURCES ${TEST_FILE} ${ARGN})
        get_filename_component(CURRENT_DIRECTORY_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
        set_target_properties(${TESTNAME} PROPERTIES FOLDER ${CURRENT_DIRECTORY_NAME})

        # Add this test to a global property
        set_property(GLOBAL APPEND PROPERTY MARS_TESTS ${TESTNAME})
        # Does not work on cluster env where you cannot run MPI on master node
        # So for the moment it will only work on apple laptops (the test can be
        # run with ./mars_test anyway but not with make test)
        if(APPLE
           OR WIN32
           OR MARS_ALLOW_DISCOVER_TESTS)
            gtest_discover_tests(
                ${TESTNAME}
                WORKING_DIRECTORY ${PROJECT_DIR}
                PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_DIR}")
        endif()

    endmacro()

endif(MARS_ENABLE_TESTS)
