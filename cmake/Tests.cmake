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
    else()
        # Add include directories for system GTest
        include_directories(${GTEST_INCLUDE_DIRS})
        set(MARS_G_TEST_LIBRARIES GTest::GTest)
    endif()

    include(GoogleTest)

    macro(mars_add_test TESTNAME TEST_FILE)
        # Create the test executable
        add_executable(${TESTNAME} ${TEST_FILE} ${ARGN})
        target_link_libraries(${TESTNAME} mars ${MARS_G_TEST_LIBRARIES})

        # Link against mars_unstructured for tests that need it
        if(TARGET mars_unstructured)
            # Check if the file path contains "unstructured"
            if("${TEST_FILE}" MATCHES "unstructured" OR
               "${TEST_FILE}" MATCHES "domain.hpp")
                # Link to mars_unstructured
                target_link_libraries(${TESTNAME} mars_unstructured)

                # Add CUDA-specific settings
                if(CMAKE_CUDA_COMPILER)
                    # Check for .cu files in sources
                    set(HAS_CUDA_FILES FALSE)
                    foreach(SOURCE_FILE ${TEST_FILE} ${ARGN})
                        get_filename_component(FILE_EXT ${SOURCE_FILE} EXT)
                        if(FILE_EXT STREQUAL ".cu" OR FILE_EXT STREQUAL ".cuh")
                            set(HAS_CUDA_FILES TRUE)
                        endif()
                    endforeach()

                    # Add CUDA flags for CUDA files
                    if(HAS_CUDA_FILES)
                        target_compile_options(${TESTNAME} PRIVATE
                            $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda>
                            $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>)
                        
                        set_target_properties(${TESTNAME} PROPERTIES
                            CUDA_SEPARABLE_COMPILATION ON)
                    endif()
                endif()
            endif()
        endif()

        # Add to test system
        gtest_add_tests(TARGET ${TESTNAME} SOURCES ${TEST_FILE} ${ARGN})
        get_filename_component(CURRENT_DIRECTORY_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
        set_target_properties(${TESTNAME} PROPERTIES FOLDER ${CURRENT_DIRECTORY_NAME})

        # Add this test to a global property
        set_property(GLOBAL APPEND PROPERTY MARS_TESTS ${TESTNAME})
        
        # Handle platform-specific test discovery
        if(APPLE OR WIN32 OR MARS_ALLOW_DISCOVER_TESTS)
            gtest_discover_tests(
                ${TESTNAME}
                WORKING_DIRECTORY ${PROJECT_DIR}
                PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_DIR}")
        endif()
    endmacro()

endif(MARS_ENABLE_TESTS)