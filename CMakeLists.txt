cmake_minimum_required(VERSION 3.5...3.19)
cmake_policy(SET CMP0048 NEW)
project(mars 
    VERSION 0.1.0
    LANGUAGES "CXX"
    HOMEPAGE_URL "https://bitbucket.org/zulianp/mars/src/master/"
    DESCRIPTION "Mars is an open-source mesh management library designed to handle N-dimensional elements (N <= 4)")



# FIXME make them usable from avaialble installation
# Don't know what this means.
add_subdirectory(external/googletest)
add_subdirectory(external/benchmark)


# Std CMake imports
include(CMakePackageConfigHelpers)

# Cmake imports.
include(cmake/MarsDefaults.cmake)
include(cmake/MarsOptions.cmake)
include(cmake/MarsCMakeFunctions.cmake)
include(cmake/MarsDependencies.cmake)

# ##############################################################################
# ##############################################################################
# Library
# ##############################################################################

list(APPEND MARS_HEADERS mars.hpp mars_env.hpp)
list(APPEND MARS_SOURCES mars.cpp mars_env.cpp)

set(MARS_MODULES
    base
    core
    distributed
    examples
    io
    kokkos)

foreach(MODULE ${MARS_MODULES})
    add_subdirectory("${MODULE}")
endforeach(MODULE)

add_library(mars STATIC ${MARS_HEADERS} ${MARS_SOURCES})

# ##############################################################################
# include

foreach(include ${MARS_BUILD_INCLUDES})
    target_include_directories(mars BEFORE
                               PUBLIC $<BUILD_INTERFACE:${include}>)
endforeach()


target_include_directories(
    mars BEFORE
    PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
           $<INSTALL_INTERFACE:include> $<INSTALL_INTERFACE:include/kernels>)


target_include_directories(mars PUBLIC ${MARS_DEP_INCLUDES})



# #############################################################################
# link 

target_link_libraries(mars PUBLIC ${MARS_DEP_LIBRARIES})
# target_link_libraries(mars PUBLIC ${MPI_LIBRARIES})



# #############################################################################
# properties

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE
        "Release"
        CACHE STRING "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel." FORCE)

    message(STATUS "[Status] CMAKE_BUILD_TYPE=Release")
endif(NOT CMAKE_BUILD_TYPE)


if(NOT DEFINED CMAKE_CXX_STANDARD_REQUIRED)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()
if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 14)
endif()


# add_library(mars STATIC ${MARS_HEADERS} ${MARS_SOURCES})

# target_link_libraries(mars ${MARS_LIBRARIES})

target_link_libraries(
    mars # $<$<CONFIG:RelWithDebInfo>:asan> # $<$<CONFIG:RelWithDebInfo>:ubsan>
)

target_compile_options(
    mars
    PUBLIC -Wall -Wextra
           # $<$<CONFIG:RelWithDebInfo>:-fsanitize=address -fsanitize=leak
           # -fsanitize=undefined -fsanitize-address-use-after-scope>
    # $<$<CONFIG:RelWithDebInfo>:-fsanitize=memory
    # -fsanitize-memory-track-origins>
)

# target_link_libraries(mars ${MPI_LIBRARIES})
target_compile_features(mars PUBLIC cxx_std_14)
target_include_directories(mars PUBLIC .)

# ##############################################################################
# ##############################################################################
# Extra targets
# ##############################################################################

if(MARS_ENABLE_BENCHMARK)
    include(bench/benchmarks.cmake)
endif()

if(MARS_ENABLE_TESTING)
    include(tests/tests.cmake)
endif()

if(MARS_ENABLE_VTK)
    add_subdirectory(backend/vtk)
endif()


if(MARS_ENABLE_ADIOS2)
    find_package(ADIOS2 REQUIRED)
    add_subdirectory(backend/adios2)
endif()

if(TRY_WITH_MOONOLITH)
    add_subdirectory(moonolith_adapter)

    if(WITH_PAR_MOONOLITH)
        target_link_libraries(mars mars_moonolith)
    endif()
endif()


# ##############################################################################
# ##############################################################################
# ##############################################################################

                ##################  EXECUTABLES AND INSTALL ################## 

add_executable(mars_exec main.cpp)
target_link_libraries(mars_exec mars)

target_include_directories(
    mars PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/external/cxxopts/include)

if(TRY_WITH_MOONOLITH)
    add_subdirectory(moonolith_adapter)

    if(WITH_PAR_MOONOLITH)
        target_link_libraries(mars mars_moonolith)
    endif()
endif()

target_include_directories(mars PUBLIC ${CMAKE_BINARY_DIR})
target_include_directories(mars PUBLIC ${MARS_INCLUDES})

message(STATUS "${MARS_INCLUDES}")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/base/mars_config.hpp.in
               ${CMAKE_BINARY_DIR}/mars_config.hpp)
include_directories(${CMAKE_BINARY_DIR})

# if(MARS_WITH_VTK)
#     install(
#         TARGETS mars_vtk
#         RUNTIME DESTINATION bin
#         LIBRARY DESTINATION lib
#         ARCHIVE DESTINATION lib)
# endif()



# if(MARS_WITH_ADIOS2)
#     install(
#         TARGETS mars_adios2
#         RUNTIME DESTINATION bin
#         LIBRARY DESTINATION lib
#         ARCHIVE DESTINATION lib)
# endif()

# if(Kokkos_FOUND)
# install(
#     TARGETS mars_examples
#     RUNTIME DESTINATION bin
#     LIBRARY DESTINATION lib
#     ARCHIVE DESTINATION lib)
# endif()

install(
    TARGETS mars mars_exec
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)

install(FILES ${MARS_HEADERS} DESTINATION include)
install(FILES ${CMAKE_BINARY_DIR}/mars_config.hpp DESTINATION include)

# ##############################################################################

