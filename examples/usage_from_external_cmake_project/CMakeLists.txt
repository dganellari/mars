cmake_minimum_required(VERSION 3.5...3.19)
project(mars_usage)

include(CMakeFindDependencyMacro)

find_package(Mars 0 REQUIRED)


# self-contained smoke: compiles against installed headers + links Mars::mars
add_executable(mars_example main.cpp)

# import all the dependencies and targets
target_link_libraries(mars_example PUBLIC Mars::mars)

# If you want to test it
enable_testing()
add_test(NAME Mars.mars_example COMMAND mars_example)
