if(MARS_ENABLE_BENCHMARK)
    # ##########################################################################
    set(BENCHMARK_ENABLE_TESTING OFF)
    # message(STATUS "GOOGLETEST_PATH=${GOOGLETEST_PATH}")

    find_package(benchmark QUIET)

    if(NOT benchmark_FOUND)
        include(FetchContent)

        FetchContent_Declare(
            benchmark
            GIT_REPOSITORY https://github.com/google/benchmark.git
            GIT_TAG v1.5.2)

        FetchContent_GetProperties(benchmark)

        if(NOT benchmark_POPULATED)
            FetchContent_Populate(benchmark)
            add_subdirectory(${benchmark_SOURCE_DIR} ${benchmark_BINARY_DIR})
        endif()

        set_target_properties(benchmark PROPERTIES FOLDER extern)
        set_target_properties(benchmark_main PROPERTIES FOLDER extern)

        set(MARS_BENCH_LIBRARIES benchmark)
    else()
        set(MARS_BENCH_LIBRARIES benchmark::benchmark)
    endif()

    # ##########################################################################

    set(MARS_BENCH_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bench)

    list(APPEND BENCH_MODULES .)

    set(LOCAL_HEADERS "")
    set(LOCAL_SOURCES "")
    find_project_files(MARS_BENCH_DIR BENCH_MODULES LOCAL_HEADERS
                       LOCAL_SOURCES)

    if(NOT LOCAL_SOURCES)
        # message(WARNING "For some reason benchmark was not found
        # automatically")

        list(APPEND LOCAL_SOURCES
             ${MARS_BENCH_DIR}/mars_MeshGenerationBenchmark.cpp)
    endif()

    add_executable(
        mars_bench ${CMAKE_CURRENT_SOURCE_DIR}/bench/mars_benchmarks.cpp
                        ${LOCAL_SOURCES})

    target_link_libraries(mars_bench ${MARS_BENCH_LIBRARIES})

    target_link_libraries(mars_bench mars)

    target_include_directories(
        mars_bench PRIVATE "$<BUILD_INTERFACE:${MARS_BENCH_DIR}>")
    target_include_directories(mars_bench PRIVATE " $<BUILD_INTERFACE:.>")
    target_include_directories(mars_bench
                               PRIVATE " $<BUILD_INTERFACE:${BENCH_MODULES}>")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MARS_DEV_FLAGS}")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g ${MARS_DEV_FLAGS}")


endif(MARS_ENABLE_BENCHMARK)
