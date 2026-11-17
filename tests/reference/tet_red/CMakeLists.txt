# A standalone build keeps the geometry regression available without a GPU toolchain.
cmake_minimum_required(VERSION 3.24)
project(mars_tet_red_check LANGUAGES CXX)
enable_testing()

add_executable(tet_red_host tet_red_host.cpp)
target_compile_features(tet_red_host PRIVATE cxx_std_20)
target_include_directories(tet_red_host PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../../..")
add_test(NAME tet_red_host COMMAND tet_red_host)
set_tests_properties(tet_red_host PROPERTIES LABELS "amr;host" TIMEOUT 60)
