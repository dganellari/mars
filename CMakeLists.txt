cmake_minimum_required(VERSION 2.8.12)

project(mars)

option(TRY_WITH_MOONOLITH
       "Use -DTRY_WITH_MOONOLITH=ON for enabling mesh transfer functions." OFF)
option(TRY_WITH_KOKKOS
       "Use -DTRY_WITH_KOKKOS=ON for enabling mesh transfer functions." OFF)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE
        "Release"
        CACHE STRING "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel." FORCE)

    message(STATUS "[Status] CMAKE_BUILD_TYPE=Release")
endif(NOT CMAKE_BUILD_TYPE)

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
    # set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER}) set(CMAKE_CXX_COMPILER_DEBUG
    # ${MPI_CXX_COMPILER})
    set(WITH_MPI ON)
endif()

option(MARS_USE_CUDA "Build mars with cuda support" OFF)

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
    find_package(Kokkos HINTS ${KOKKOS_DIR}/lib/CMake/Kokkos
                 ${TRILINOS_DIR}/lib/cmake/Kokkos QUIET)

    if(Kokkos_FOUND)
        message("\nFound Kokkos!  Here are the details: ")
        message("   Kokkos_CXX_COMPILER = ${Kokkos_CXX_COMPILER}")
        message("   Kokkos_C_COMPILER = ${Kokkos_C_COMPILER}")

        if(Kokkos_CXX_COMPILER)
            set(CMAKE_C_COMPILER ${Kokkos_C_COMPILER})
            set(CMAKE_CXX_COMPILER ${Kokkos_CXX_COMPILER})
            message("Setting CMAKE_CXX_COMPILER to Kokkos_CXX_COMPILER")
        endif()

        set(WITH_KOKKOS ON)

        # find_package( KokkosContainers HINTS
        # ${KOKKOS_DIR}/lib/CMake/KokkosContainers
        # ${TRILINOS_DIR}/lib/cmake/KokkosContainers QUIET)

        # if(KokkosContainers_FOUND) message( STATUS
        # "KokkosContainers_TPL_LIBRARIES=${KokkosContainers_TPL_LIBRARIES}" )
        # link_libraries(${KokkosContainers_TPL_LIBRARIES}) endif()

        find_package(KokkosKernels HINTS ${KOKKOS_DIR}/lib/CMake/KokkosKernels
                     ${TRILINOS_DIR}/lib/cmake/KokkosKernels QUIET)

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

list(
    APPEND
    MARS_HEADERS
    base/mars.hpp
    base/mars_base.hpp
    base/mars_benchmark.hpp
    base/mars_bisection.hpp
    base/mars_communicator.hpp
    base/mars_config.hpp.in
    base/mars_dof_map.hpp
    base/mars_dual_graph.hpp
    base/mars_edge.hpp
    base/mars_edge_element_map.hpp
    base/mars_edge_node_map.hpp
    base/mars_edge_select.hpp
    base/mars_edge_split.hpp
    base/mars_edge_split_pool.hpp
    base/mars_extract_surface.hpp
    base/mars_fwd.hpp
    base/mars_globals.hpp
    base/mars_imesh.hpp
    base/mars_lagrange_element.hpp
    base/mars_longest_edge.hpp
    base/mars_matrix.hpp
    base/mars_mesh.hpp
    base/mars_mesh_partition.hpp
    base/mars_mesh_reader.hpp
    base/mars_mesh_writer.hpp
    base/mars_newest_vertex.hpp
    base/mars_node_rank.hpp
    base/mars_oldest_edge.hpp
    base/mars_par_bisection.hpp
    base/mars_par_edge_split_pool.hpp
    base/mars_par_mesh.hpp
    base/mars_partitioned_bisection.hpp
    base/mars_quality.hpp
    base/mars_ranked_edge.hpp
    base/mars_red_green_refinement.hpp
    base/mars_simplex.hpp
    base/mars_static_math.hpp
    base/mars_stream.hpp
    base/mars_test.hpp
    base/mars_tracker.hpp
    base/mars_utils.hpp
    base/mars_vector.hpp
    base/mars_visualization.hpp
    base/mars_vtk_writer.hpp
    base/mars_non_simplex.hpp
    base/mars_lepp_bisection.hpp
    base/mars_lepp_benchmark.hpp
    base/mars_prelepp_benchmark.hpp
    base/mars_prelepp_bisection.hpp
    base/mars_mesh_reader.hpp
    base/mars_mesh_writer.hpp
    kokkos/mars_sub_view.hpp
    kokkos/mars_point.hpp
    kokkos/mars_mesh_generation.hpp
    kokkos/mars_device_vector.hpp
    kokkos/mars_memory.hpp)

if(Kokkos_FOUND)
    list(
        APPEND
        MARS_HEADERS
        kokkos/mars_imesh_kokkos.hpp
        kokkos/mars_mesh_kokkos.hpp
        kokkos/mars_test_kokkos.hpp
        kokkos/mars_utils_kokkos.hpp
        kokkos/mars_edge_select_kokkos.hpp
        kokkos/mars_static_math_kokkos.hpp
        kokkos/mars_simplex_kokkos.hpp
        kokkos/mars_non_simplex_kokkos.hpp
        kokkos/mars_edge_kokkos.hpp
        kokkos/mars_lepp_benchmark_kokkos.hpp
        kokkos/mars_bisection_kokkos.hpp
        kokkos/mars_edge_node_map_kokkos.hpp
        kokkos/mars_edge_element_map_kokkos.hpp
        kokkos/mars_mark_kokkos.hpp
        distributed/sfc/mars_sfc_code.hpp
        distributed/sfc/mars_sfc_generation.hpp)
endif()

if(WITH_MPI)
    list(
        APPEND
        MARS_HEADERS
        distributed/mars_distributed_octant.hpp
        distributed/mars_distributed_dof.hpp
        distributed/mars_distributed_data_management.hpp
        distributed/mars_distributed_user_data.hpp
        distributed/mars_distributed_utils.hpp
        distributed/kokkos/mars_distributed_mesh_kokkos.hpp
        distributed/kokkos/mars_distributed_non_simplex_kokkos.hpp
        distributed/kokkos/mars_distributed_simplex_kokkos.hpp
        distributed/mars_distributed_mesh_generation.hpp
        distributed/mars_test_mpi.hpp
        distributed/communication/util/mars_pp_util.hpp
        distributed/communication/mars_context.hpp
        distributed/communication/mars_gathered_vector.hpp
        distributed/communication/mars_mpi_error.hpp
        distributed/communication/mars_mpi.hpp
        distributed/communication/mars_distributed_context.hpp
        distributed/communication/mars_execution_context.hpp
        distributed/mars_env/mars_mpi_guard.hpp
        examples/mars_advection.hpp)
endif()

list(APPEND MARS_SOURCES base/mars_communicator.cpp base/mars_simplex.cpp
     base/mars_test.cpp)

if(WITH_MPI)
    list(
        APPEND
        MARS_SOURCES
        distributed/communication/mars_mpi.cpp
        distributed/communication/mars_mpi_error.cpp
        distributed/communication/mars_mpi_context.cpp
        distributed/communication/mars_execution_context.cpp
        distributed/communication/mars_dry_run_context.cpp
        examples/mars_poisson_operator.cpp)
endif()

set(CMAKE_CXX_DEBUG "-g")

add_library(mars STATIC ${MARS_HEADERS} ${MARS_SOURCES})
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
# target_compile_features(mars PUBLIC cxx_std_11)

if(WITH_MPI)
    target_compile_definitions(mars PUBLIC MARS_HAVE_MPI)
endif()

add_definitions("-std=c++11")
target_include_directories(mars PUBLIC base)
target_include_directories(mars PUBLIC kokkos)
target_include_directories(mars PUBLIC distributed)
target_include_directories(mars PUBLIC distributed/kokkos)
target_include_directories(mars PUBLIC distributed/sfc)
target_include_directories(mars PUBLIC distributed/mars_env)
target_include_directories(mars PUBLIC distributed/communication)
target_include_directories(mars PUBLIC distributed/communication/util)
target_include_directories(mars PUBLIC examples)

if(Kokkos_FOUND)
    message(" Kokkos_INCLUDE_DIRS = ${Kokkos_INCLUDE_DIRS}")
    message(" Kokkos_LIBRARIES = ${Kokkos_LIBRARIES}")
    message(" Kokkos_TPL_LIBRARIES = ${Kokkos_TPL_LIBRARIES}")
    message(" Kokkos_LIBRARY_DIRS = ${Kokkos_LIBRARY_DIRS}")

    if(MARS_USE_CUDA AND NOT Kokkos_TPL_LIBRARIES MATCHES "cuda")
        message(
            FATAL_ERROR
                "Enable Kokkos Cuda or unset MARS_USE_CUDA to continue with OpenMP!"
        )
    endif()

    target_include_directories(mars PUBLIC ${Kokkos_TPL_INCLUDE_DIRS}
                                           ${Kokkos_INCLUDE_DIRS})

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
        message(STATUS "KokkosKernels_FOUND=${KokkosKernels_FOUND}")
        # target_link_libraries(mars PUBLIC Kokkos::Kokkos)

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

add_executable(mars_exec main.cpp)
target_link_libraries(mars_exec mars)

if(TRY_WITH_MOONOLITH)
    add_subdirectory(moonolith_adapter)
endif()

add_subdirectory(io)

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

# add_definitions("-DWITH_MOONOLITH") set(MOONOLITH_DIR
# "$ENV{INSTALL_DIR}/moon") target_include_directories(mars_exec PUBLIC
# ${MOONOLITH_DIR}/include) target_link_libraries(mars_exec
# -L${MOONOLITH_DIR}/lib/ -leps -lexternal -lintersection -lio -lmeshing
# -lmoonolith -lquadrature -lshape -lutility )
