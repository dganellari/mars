cmake_minimum_required(VERSION 3.8.1)

project(mars)
include_directories(.)

list(APPEND MARS_HEADERS 
	simplex.hpp
	)

list(APPEND MARS_SOURCES
	simplex.cpp
	)

set(CMAKE_CXX_DEBUG "-g")

add_executable(mars_exec main.cpp ${MARS_HEADERS} ${MARS_SOURCES})
target_compile_features(mars_exec PUBLIC cxx_std_11)