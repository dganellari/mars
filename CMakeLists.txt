cmake_minimum_required(VERSION 2.8.12)

project(mars)

# FIXME make them usable from avaialble installation
add_subdirectory(external/googletest)
add_subdirectory(external/benchmark)

option(TRY_WITH_MOONOLITH
       "Use -DTRY_WITH_MOONOLITH=ON for enabling mesh transfer functions." OFF)
option(TRY_WITH_KOKKOS
       "Use -DTRY_WITH_KOKKOS=ON for enabling mesh transfer functions." OFF)

option(MARS_USE_CUDA "Build mars with cuda support" OFF)
option(MARS_NO_RDMA "Build mars with NO RDMA support" OFF)
option(MARS_USE_CUDA "Build mars with cuda support" OFF)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE
        "Release"
        CACHE STRING "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel." FORCE)

    message(STATUS "[Status] CMAKE_BUILD_TYPE=Release")
endif(NOT CMAKE_BUILD_TYPE)

if(UNIX AND NOT APPLE)
    set(LINUX TRUE)
endif()

# ##############################################################################
# ##############################################################################
# ##############################################################################

function(find_project_files rootPath dirPaths headers sources)
    set(verbose TRUE)

    set(theaders ${${headers}})
    set(tsources ${${sources}})

    set(ROOT_PATH ${${rootPath}})
    set(DIR_PATHS ${${dirPaths}})

    foreach(INCLUDE_PATH ${DIR_PATHS})
        include_directories(${ROOT_PATH}/${INCLUDE_PATH})

        file(GLOB TEMP_HPPSRC "${ROOT_PATH}/${INCLUDE_PATH}/*.cpp")
        file(GLOB TEMP_SRC "${ROOT_PATH}/${INCLUDE_PATH}/*.c")
        file(GLOB TEMP_HPPDR "${ROOT_PATH}/${INCLUDE_PATH}/*.hpp")
        file(GLOB TEMP_HDR "${ROOT_PATH}/${INCLUDE_PATH}/*.h")

        source_group(
            ${INCLUDE_PATH} FILES ${TEMP_HPPDR}; ${TEMP_HDR}; ${TEMP_HPPSRC};
                                  ${TEMP_SRC}; ${TEMP_UI})

        set(tsources ${tsources}; ${TEMP_SRC}; ${TEMP_HPPSRC})
        set(theaders ${theaders}; ${TEMP_HDR}; ${TEMP_HPPDR})
    endforeach(INCLUDE_PATH)

    set(${headers}
        ${theaders}
        PARENT_SCOPE)
    set(${sources}
        ${tsources}
        PARENT_SCOPE)
endfunction()

# ##############################################################################
# ##############################################################################
# ##############################################################################

find_package(MPIExtended)

if(MPI_FOUND)
    set(WITH_MPI ON)
endif()

if(NOT DEFINED CMAKE_CXX_STANDARD_REQUIRED)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()
if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 14)
endif()

if(TRY_WITH_KOKKOS)
    if(NOT TRILINOS_DIR)
        message(STATUS "Setting TRILINOS_DIR to $ENV{TRILINOS_DIR}")
        set(TRILINOS_DIR
            $ENV{TRILINOS_DIR}
            CACHE PATH "Directory where Kokkos is installed")
    endif()

    if(NOT KOKKOS_DIR)
        message(STATUS "Setting KOKKOS_DIR to $ENV{KOKKOS_DIR}")
        set(KOKKOS_DIR
            $ENV{KOKKOS_DIR}
            CACHE PATH "Directory where Kokkos is installed")
    endif()

    # FIND_PACKAGE(Trilinos PATHS ${TRILINOS_DIR}/lib/cmake/Trilinos QUIET)

    find_package(
        Kokkos
        HINTS
        ${KOKKOS_DIR}
        ${KOKKOS_DIR}/lib/CMake/Kokkos
        ${KOKKOS_DIR}/lib64/CMake/Kokkos
        ${TRILINOS_DIR}
        ${TRILINOS_DIR}/lib/cmake/Kokkos
        ${TRILINOS_DIR}/lib64/cmake/Kokkos
        QUIET)

    if(Kokkos_FOUND)
        if (TARGET Kokkos::kokkos)
            message(STATUS "Kokkos::kokkos is a target, get include directories from the target")
            get_target_property(Kokkos_INCLUDE_DIRS Kokkos::kokkos INTERFACE_INCLUDE_DIRECTORIES)
            get_target_property(Kokkos_LIBRARIES Kokkos::kokkos INTERFACE_LINK_LIBRARIES)
            get_target_property(Kokkos_LIBRARY_DIRS Kokkos::kokkos INTERFACE_LINK_DIRECTORIES)
        endif()
        message("\nFound Kokkos!  Here are the details: ")
        message(" Kokkos_CXX_COMPILER = ${Kokkos_CXX_COMPILER}")
        message(" Kokkos_INCLUDE_DIRS = ${Kokkos_INCLUDE_DIRS}")
        message(" Kokkos_LIBRARIES = ${Kokkos_LIBRARIES}")
        message(" Kokkos_TPL_LIBRARIES = ${Kokkos_TPL_LIBRARIES}")
        message(" Kokkos_LIBRARY_DIRS = ${Kokkos_LIBRARY_DIRS}")
        message(" Kokkos CUDA Enabled = ${Kokkos_ENABLE_CUDA}")

        # message("   Kokkos_CXX_COMPILER = ${Kokkos_CXX_COMPILER}") message("
        # Kokkos_C_COMPILER = ${Kokkos_C_COMPILER}")

        if(Kokkos_CXX_COMPILER)
            set(CMAKE_C_COMPILER ${Kokkos_C_COMPILER})
            set(CMAKE_CXX_COMPILER ${Kokkos_CXX_COMPILER})
            # message( "Setting CMAKE_CXX_COMPILER to
            # Kokkos_CXX_COMPILER=${Kokkos_CXX_COMPILER}" )
        endif()

        set(WITH_KOKKOS ON)

        find_package(
            KokkosKernels
            HINTS
            ${KOKKOS_DIR}
            ${KOKKOS_DIR}/lib/CMake/KokkosKernels
            ${KOKKOS_DIR}/lib64/cmake/KokkosKernels
            ${TRILINOS_DIR}
            ${TRILINOS_DIR}/lib/cmake/KokkosKernels
            ${TRILINOS_DIR}/lib64/cmake/KokkosKernels
            REQUIRED)

        find_package(OpenMP)

        if(OPENMP_FOUND)
            set(CMAKE_Fortran_FLAGS
                "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
            set(CMAKE_CXX_LINK_FLAGS
                "${CMAKE_CXX_LINK_FLAGS} ${OpenMP_CXX_FLAGS}")
        endif()

    else()
        message(WARNING "Could not find Kokkos!")
        if(MARS_USE_CUDA)
            message(FATAL_ERROR "Could not use CUDA without Kokkos!")
        endif()
    endif()

endif()

# ##############################################################################

list(APPEND MARS_HEADERS mars.hpp mars_env.hpp)
list(APPEND MARS_SOURCES mars.cpp mars_env.cpp)

add_subdirectory(base)
add_subdirectory(core)
add_subdirectory(distributed)
add_subdirectory(examples)
add_subdirectory(io)
add_subdirectory(kokkos)

add_library(mars STATIC ${MARS_HEADERS} ${MARS_SOURCES})

target_link_libraries(mars ${MARS_LIBRARIES})

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

target_link_libraries(mars ${MPI_LIBRARIES})
target_compile_features(mars PUBLIC cxx_std_14)
target_include_directories(mars PUBLIC .)

if(Kokkos_FOUND)

    message(" Kokkos_INCLUDE_DIRS = ${Kokkos_INCLUDE_DIRS}")
    message(" Kokkos_LIBRARIES = ${Kokkos_LIBRARIES}")
    message(" Kokkos_TPL_LIBRARIES = ${Kokkos_TPL_LIBRARIES}")
    message(" Kokkos_LIBRARY_DIRS = ${Kokkos_LIBRARY_DIRS}")

    if(MARS_USE_CUDA AND (NOT DEFINED Kokkos_ENABLE_CUDA OR NOT ${Kokkos_ENABLE_CUDA}))
            message(
                FATAL_ERROR
                    "Enable Kokkos Cuda or unset MARS_USE_CUDA to continue with OpenMP!"
            )
    endif()

    if(Kokkos_ENABLE_CUDA)
        kokkos_check(OPTIONS CUDA_LAMBDA)
    endif()

    target_include_directories(mars SYSTEM PUBLIC ${Kokkos_TPL_INCLUDE_DIRS} ${Kokkos_INCLUDE_DIRS})

    message(
        STATUS
            "Kokkos_INCLUDE_DIRS=${Kokkos_INCLUDE_DIRS}, ${Kokkos_TPL_INCLUDE_DIRS}"
    )

    if(Kokkos_CXX_COMPILER)
        target_link_libraries(mars ${Kokkos_LIBRARIES} ${Kokkos_TPL_LIBRARIES})
    else()
        target_link_libraries(mars ${Kokkos_LIBRARIES} ${Kokkos_TPL_LIBRARIES}
                              -L${Kokkos_LIBRARY_DIRS})
    endif()

    if(KokkosKernels_FOUND)
        if (TARGET Kokkos::kokkoskernels)
            target_link_libraries(mars Kokkos::kokkoskernels)
        else()
            target_include_directories(mars PUBLIC ${KokkosKernels_TPL_INCLUDE_DIRS}
                                                           ${KokkosKernels_INCLUDE_DIRS})

            if(Kokkos_CXX_COMPILER)
                target_link_libraries(mars ${KokkosKernels_LIBRARIES}
                                      ${KokkosKernels_TPL_LIBRARIES})
            else()
                target_link_libraries(
                    mars ${KokkosKernels_LIBRARIES} ${KokkosKernels_TPL_LIBRARIES}
                    -L${KokkosKernels_LIBRARY_DIRS})
            endif()
        endif()
    endif()

endif()

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

target_link_libraries(mars ${MARS_LIBRARIES})
target_include_directories(mars PUBLIC ${CMAKE_BINARY_DIR})
target_include_directories(mars PUBLIC ${MARS_INCLUDES})

message(STATUS "${MARS_INCLUDES}")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/base/mars_config.hpp.in
               ${CMAKE_BINARY_DIR}/mars_config.hpp)
include_directories(${CMAKE_BINARY_DIR})

install(
    TARGETS mars mars_exec
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)

install(FILES ${MARS_HEADERS} DESTINATION include)
install(FILES ${CMAKE_BINARY_DIR}/mars_config.hpp DESTINATION include)

# ##############################################################################

if(Kokkos_FOUND)
    add_executable(st_example examples/mars_st_example.cpp)
    target_link_libraries(st_example mars mars_io)
endif()

# ##############################################################################
#ADIOS STUFF

find_package(ADIOS2 REQUIRED)
target_link_libraries(mars adios2::cxx11_mpi MPI::MPI_C)
if(ADIOS2_FOUND)
    add_executable(adios_example examples/distributed_dm/mars_example_adios2.cpp)
    target_link_libraries(adios_example mars mars_io)
endif()

include(tests/tests.cmake)
include(bench/benchmarks.cmake)
