set(MARS_BENCH_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bench)

list(APPEND BENCH_MODULES
    .
)


# Now simply link against gtest or gtest_main as needed. Eg
add_executable(mars_bench ${MARS_BENCH_DIR}/mars_benchmarks.cpp)
target_link_libraries(mars_bench benchmark::benchmark)
add_test(NAME mars_bench COMMAND mars_bench)

set(LOCAL_HEADERS "")
set(LOCAL_SOURCES "")
find_project_files(MARS_BENCH_DIR BENCH_MODULES LOCAL_HEADERS LOCAL_SOURCES)
target_sources(mars_bench PRIVATE ${LOCAL_SOURCES})
target_link_libraries(mars_bench mars)

target_include_directories(mars_bench PRIVATE ${MARS_BENCH_DIR})
target_include_directories(mars_bench PRIVATE .)
target_include_directories(mars_bench PRIVATE ${BENCH_MODULES})


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MARS_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g ${MARS_DEV_FLAGS}")
