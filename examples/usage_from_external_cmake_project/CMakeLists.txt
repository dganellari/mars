cmake_minimum_required(VERSION 3.5...3.19)
project(mars_usage)

include(CMakeFindDependencyMacro)

find_package(Mars 0 REQUIRED)
find_package(Kokkos REQUIRED)
find_package(KokkosKernels REQUIRED)

# we use one of our examples here
add_executable(evt_exe ${CMAKE_SOURCE_DIR}/../mars_example_post_install.cpp)

# import all the dependencies and targets
target_link_libraries(evt_exe Mars::mars)

# If you want to test it
enable_testing()
add_test(NAME Mars.evt_exe COMMAND evt_exe ${CMAKE_SOURCE_DIR}/../data)