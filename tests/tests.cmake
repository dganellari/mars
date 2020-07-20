set(MARS_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests)

list(APPEND TEST_MODULES
    .
)


# Now simply link against gtest or gtest_main as needed. Eg
add_executable(mars_test ${MARS_TEST_DIR}/test.cpp)
target_link_libraries(mars_test gtest_main)
add_test(NAME mars_test COMMAND mars_test)

set(LOCAL_HEADERS "")
set(LOCAL_SOURCES "")
find_project_files(MARS_TEST_DIR TEST_MODULES LOCAL_HEADERS LOCAL_SOURCES)
target_sources(mars_test PRIVATE ${LOCAL_SOURCES})
target_link_libraries(mars_test mars)

target_include_directories(mars_test PRIVATE ${MARS_TEST_DIR})
target_include_directories(mars_test PRIVATE .)
target_include_directories(mars_test PRIVATE ${TEST_MODULES})


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MARS_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g ${MARS_DEV_FLAGS}")
