cmake_minimum_required(VERSION 3.23.0 FATAL_ERROR)

project(MeshTest VERSION 0.1.0 LANGUAGES C CXX)

# project(MeshTest VERSION 0.1.0 LANGUAGES C CXX)

# For Trilinos
# set(Trilinos_DIR /home/roos/apps/trilinos/develop_cuda/lib64/cmake/Trilinos)
# set(Trilinos_DIR /home/roos/apps/trilinos/develop_gcc/lib64/cmake/Trilinos)
# set(Trilinos_DIR /home/roos/apps/trilinos/develop_omp/lib64/cmake/Trilinos)
set(CMAKE_PREFIX_PATH ${Trilinos_DIR} ${CMAKE_PREFIX_PATH})

set(CMAKE_VERBOSE_MAKEFILE on)

# set( Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE true)
add_definitions(-DUSE_STK_SIMD_NONE)

set(CMAKE_CXX_STANDARD 20)

# find_package(MPI REQUIRED)
message(STATUS "OMPI_CXX environment variable: $ENV{OMPI_CXX}")

# add_compile_options("-DOMPI_SKIP_MPICXX")
find_package(Trilinos REQUIRED)

# Echo trilinos build info just for fun
MESSAGE("\nFound Trilinos!  Here are the details: ")
MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
MESSAGE("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
MESSAGE("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
MESSAGE("   CMAKE_CXX_COMPILER = ${Trilinos_CXX_COMPILER}")
MESSAGE("End of Trilinos details\n")

# COMBAK: [fab4100@posteo.net; 2025-09-11] forcing compiler below to trilinos
# compiler is not necessary when build environment is setup properly.  Results
# in infinite loop during cmake phase otherwise.
#
# Make sure to use same compilers and flags as Trilinos
# SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER})
# SET(CMAKE_C_COMPILER ${Trilinos_C_COMPILER})
# SET(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER})

SET(CMAKE_CXX_FLAGS "${Trilinos_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
SET(CMAKE_C_FLAGS "${Trilinos_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS}")
SET(CMAKE_Fortran_FLAGS "${Trilinos_Fortran_COMPILER_FLAGS} ${CMAKE_Fortran_FLAGS}")

message(STATUS "OMPI_CXX environment variable: $ENV{OMPI_CXX}")

INCLUDE_DIRECTORIES(${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS} . ./external ./external_old)
LINK_DIRECTORIES(${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})

find_package(yaml-cpp REQUIRED)

# Also controls kokkos setting of linear solver library here
set(USE_KOKKOS true)

if(USE_KOKKOS)
    message(STATUS "Enabling Kokkos support for sandbox")
    add_compile_options("-DUSE_KOKKOS")
endif()

add_compile_definitions(NO_LINSOLVE_CONF)

# add_compile_definitions(SOLVER_SINGLE_PRECISION) TODO: not currently supported
add_subdirectory("src/solver")

# Add all source files in the source directory to a variable
# file(GLOB_RECURSE SOURCES "src/*.cpp" "src/*.F" "src/*.C" "external/*.cpp" "external/*.C" "external/*.F")
file(GLOB_RECURSE SOURCES "src/*.cpp" "src/*.F" "src/*.C" "external/*.cpp" "external/*.C" "external/*.F" "external_old/*.cpp" "external_old/*.C" "external_old/*.F")

# Find all header files recursively in subdirectories
# file(GLOB_RECURSE HEADER_FILES "src/*.h" "src/*.hpp" "external/*.h" "external/*.hpp")
file(GLOB_RECURSE HEADER_FILES "src/*.h" "src/*.hpp" "external/*.h" "external/*.hpp" "external_old/*.h" "external_old/*.hpp")

list(FILTER SOURCES EXCLUDE REGEX ".*/src/solver/.*")
list(FILTER HEADER_FILES EXCLUDE REGEX ".*/src/solver/.*")

list(FILTER SOURCES EXCLUDE REGEX ".*/src/solver/.*")
list(FILTER HEADER_FILES EXCLUDE REGEX ".*/src/solver/.*")

# Enable Fortran language support (if not already done)
enable_language(Fortran)

# Add the executable target
add_executable(AssemblyTest.exe ${SOURCES})

# Include the directories containing the header files
foreach(header_file ${HEADER_FILES})
    get_filename_component(header_dir ${header_file} DIRECTORY)

    MESSAGE("   header_dir = ${header_dir}")
    target_include_directories(AssemblyTest.exe PRIVATE ${header_dir})
endforeach()

target_include_directories(AssemblyTest.exe PRIVATE ${MPI_INCLUDE_PATH} ${libLinSolve_INC})
target_include_directories(AssemblyTest.exe PRIVATE ./eigen)

TARGET_LINK_LIBRARIES(AssemblyTest.exe ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${libLinSolve} ${MPI_CXX_LIBRARIES} yaml-cpp::yaml-cpp)

# target_compile_options(AssemblyTest.exe PRIVATE -fvisibility=hidden)

# target_link_libraries(${libLinSolve} ${MPI_CXX_LIBRARIES})
