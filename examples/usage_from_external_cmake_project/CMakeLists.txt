cmake_minimum_required(VERSION 3.5...3.19)
project(mars_usage)

include(CMakeFindDependencyMacro)

find_package(Mars 0 REQUIRED)


# we use one of our examples here
add_executable(mars_example_1 ${CMAKE_SOURCE_DIR}/../mars_example_post_install.cpp)	

# import all the dependencies and targets
target_link_libraries(mars_example_1 Mars::mars)

# If you want to test it
enable_testing()
add_test(NAME Mars.mars_example_1 COMMAND mars_example_1 ${CMAKE_SOURCE_DIR}/../data)