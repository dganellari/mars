# Always enable testing when this file is included
enable_testing()

# Set a flag to indicate testing is enabled
set(MARS_TESTING_ENABLED TRUE CACHE INTERNAL "Flag to indicate testing is enabled")

# Define testing functionality if MARS_ENABLE_TESTS is set
if(MARS_ENABLE_TESTS)
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
    else()
        # Add include directories for system GTest
        include_directories(${GTEST_INCLUDE_DIRS})
        set(MARS_G_TEST_LIBRARIES GTest::GTest)
    endif()

    include(GoogleTest)

    # Macro for adding standard tests
    macro(mars_add_test TESTNAME TEST_FILE)
        # Create the test executable
        add_executable(${TESTNAME} ${TEST_FILE} ${ARGN})
        target_link_libraries(${TESTNAME} mars ${MARS_G_TEST_LIBRARIES})

        # Add to test system - always use add_test directly first
        add_test(NAME ${TESTNAME} COMMAND ${TESTNAME})

        # Then use gtest_add_tests for additional discovery
        gtest_add_tests(TARGET ${TESTNAME} SOURCES ${TEST_FILE} ${ARGN})
        
        get_filename_component(CURRENT_DIRECTORY_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
        set_target_properties(${TESTNAME} PROPERTIES FOLDER ${CURRENT_DIRECTORY_NAME})

        # Add this test to a global property
        set_property(GLOBAL APPEND PROPERTY MARS_TESTS ${TESTNAME})
        set_property(GLOBAL APPEND PROPERTY MARS_TEST_TARGETS ${TESTNAME})
        
        # Handle platform-specific test discovery
        if(APPLE OR WIN32 OR MARS_ALLOW_DISCOVER_TESTS)
            gtest_discover_tests(
                ${TESTNAME}
                WORKING_DIRECTORY ${PROJECT_DIR}
                PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_DIR}")
        endif()
        
        message(STATUS "Added test: ${TESTNAME}")
    endmacro()
else()
    # Stub macros if tests are disabled
    macro(mars_add_test TESTNAME TEST_FILE)
        # Do nothing
    endmacro()
endif()