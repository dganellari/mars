cmake_minimum_required(VERSION 2.8.12)

project(mars)

option(TRY_WITH_MOONOLITH "Use -DTRY_WITH_MOONOLITH=ON for enabling mesh transfer functions." ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel."
      FORCE)

  	message(STATUS "[Status] CMAKE_BUILD_TYPE=Release")
endif(NOT CMAKE_BUILD_TYPE)

#################################################################################
#################################################################################
#################################################################################


function(find_project_files rootPath dirPaths headers sources)
	SET(verbose TRUE)

	SET(theaders ${${headers}})
	SET(tsources ${${sources}})

	SET(ROOT_PATH ${${rootPath}})
	SET(DIR_PATHS ${${dirPaths}})

	FOREACH(INCLUDE_PATH ${DIR_PATHS})
		INCLUDE_DIRECTORIES(${ROOT_PATH}/${INCLUDE_PATH})

		file(GLOB TEMP_HPPSRC "${ROOT_PATH}/${INCLUDE_PATH}/*.cpp" )
		file(GLOB TEMP_SRC "${ROOT_PATH}/${INCLUDE_PATH}/*.c" )
		file(GLOB TEMP_HPPDR "${ROOT_PATH}/${INCLUDE_PATH}/*.hpp" )
		file(GLOB TEMP_HDR "${ROOT_PATH}/${INCLUDE_PATH}/*.h" )

		SOURCE_GROUP(${INCLUDE_PATH} FILES ${TEMP_HPPDR}; ${TEMP_HDR}; ${TEMP_HPPSRC}; ${TEMP_SRC}; ${TEMP_UI})

		SET(tsources ${tsources}; ${TEMP_SRC}; ${TEMP_HPPSRC})
		SET(theaders ${theaders}; ${TEMP_HDR}; ${TEMP_HPPDR})
	ENDFOREACH(INCLUDE_PATH)

	SET(${headers} ${theaders} PARENT_SCOPE)
	SET(${sources} ${tsources} PARENT_SCOPE)
endfunction()

#################################################################################
#################################################################################
#################################################################################

find_package(MPIExtended)

# if(MPI_FOUND)
# 	set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
# 	set(CMAKE_CXX_COMPILER_DEBUG ${MPI_CXX_COMPILER})
# 	set(WITH_MPI ON)
# endif()

option(MARS_USE_CUDA "Build mars with cuda support" OFF)


if (NOT TRILINOS_DIR)
  message(STATUS "Setting TRILINOS_DIR to $ENV{TRILINOS_DIR}")
  set(TRILINOS_DIR $ENV{TRILINOS_DIR} CACHE PATH "Directory where Kokkos is installed")
endif()

if (NOT KOKKOS_DIR)
  message(STATUS "Setting KOKKOS_DIR to $ENV{KOKKOS_DIR}")
  set(KOKKOS_DIR $ENV{KOKKOS_DIR} CACHE PATH "Directory where Kokkos is installed")
endif()

#FIND_PACKAGE(Trilinos PATHS ${TRILINOS_DIR}/lib/cmake/Trilinos QUIET)
FIND_PACKAGE(Kokkos HINTS ${KOKKOS_DIR}/lib/CMake/Kokkos ${TRILINOS_DIR}/lib/cmake/Kokkos QUIET)

IF(Kokkos_FOUND)
   	MESSAGE("\nFound Kokkos!  Here are the details: ")
   	MESSAGE("   Kokkos_CXX_COMPILER = ${Kokkos_CXX_COMPILER}")
   	MESSAGE("   Kokkos_C_COMPILER = ${Kokkos_C_COMPILER}")
   
   	IF(Kokkos_CXX_COMPILER)
		SET(CMAKE_C_COMPILER ${Kokkos_C_COMPILER})
		SET(CMAKE_CXX_COMPILER ${Kokkos_CXX_COMPILER})
		MESSAGE("Setting CMAKE_CXX_COMPILER to Kokkos_CXX_COMPILER")
   	ENDIF()
   
   	set(WITH_KOKKOS ON)
   
	FIND_PACKAGE(OpenMP)

	IF(OPENMP_FOUND)
	 	SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
	 	SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	 	SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	 	set (CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${OpenMP_CXX_FLAGS}")
	ENDIF()
   
ELSE()
  	MESSAGE(WARNING "Could not find Kokkos!")
  	IF(MARS_USE_CUDA)
    	MESSAGE(FATAL_ERROR "Could not use CUDA without Kokkos!")
  	ENDIF()
ENDIF()

list(APPEND MARS_HEADERS
	mars.hpp
	mars_base.hpp
	mars_benchmark.hpp
	mars_bisection.hpp
	mars_communicator.hpp
	mars_config.hpp.in
	mars_dof_map.hpp
	mars_dual_graph.hpp
	mars_edge.hpp
	mars_edge_element_map.hpp
	mars_edge_node_map.hpp
	mars_edge_select.hpp
	mars_edge_split.hpp
	mars_edge_split_pool.hpp
	mars_extract_surface.hpp
	mars_fwd.hpp
	mars_globals.hpp
	mars_imesh.hpp
	mars_lagrange_element.hpp
	mars_lepp_benchmark.hpp
	mars_lepp_bisection.hpp
	mars_longest_edge.hpp
	mars_matrix.hpp
	mars_mesh.hpp
	mars_mesh_partition.hpp
	mars_mesh_reader.hpp
	mars_mesh_writer.hpp
	mars_newest_vertex.hpp
	mars_node_rank.hpp
	mars_oldest_edge.hpp
	mars_par_bisection.hpp
	mars_par_edge_split_pool.hpp
	mars_par_mesh.hpp
	mars_partitioned_bisection.hpp
	mars_quality.hpp
	mars_ranked_edge.hpp
	mars_red_green_refinement.hpp
	mars_simplex.hpp
	mars_static_math.hpp
	mars_stream.hpp
	mars_test.hpp
	mars_tracker.hpp
	mars_utils.hpp
	mars_vector.hpp
	mars_visualization.hpp
	mars_vtk_writer.hpp
	generation/mars_point.hpp
	generation/mars_mesh_generation.hpp
	generation/mars_memory.hpp
	)
	
IF(Kokkos_FOUND)
	list(APPEND MARS_HEADERS
		generation/mars_imesh_kokkos.hpp	
	    generation/mars_mesh_kokkos.hpp
	    generation/mars_test_kokkos.hpp
	    generation/mars_utils_kokkos.hpp
	)
ENDIF()

list(APPEND MARS_SOURCES
	mars_communicator.cpp
	mars_simplex.cpp
	mars_test.cpp
	)

set(CMAKE_CXX_DEBUG "-g")

add_library(mars STATIC ${MARS_HEADERS} ${MARS_SOURCES})
target_link_libraries(mars ${MPI_LIBRARIES})
#target_compile_features(mars PUBLIC cxx_std_11)
add_definitions("-std=c++11")
target_include_directories(mars PUBLIC .)

add_executable(mars_exec main.cpp ${MARS_HEADERS} ${MARS_SOURCES})
target_link_libraries(mars_exec mars)

IF(Kokkos_FOUND)
	MESSAGE(" Kokkos_INCLUDE_DIRS = ${Kokkos_INCLUDE_DIRS}")
	MESSAGE(" Kokkos_LIBRARIES = ${Kokkos_LIBRARIES}")
	MESSAGE(" Kokkos_TPL_LIBRARIES = ${Kokkos_TPL_LIBRARIES}")
	MESSAGE(" Kokkos_LIBRARY_DIRS = ${Kokkos_LIBRARY_DIRS}")
	
	IF(MARS_USE_CUDA AND NOT Kokkos_TPL_LIBRARIES MATCHES "cuda")
		MESSAGE(FATAL_ERROR "Enable Kokkos Cuda or unset MARS_USE_CUDA to continue with OpenMP!")
	ENDIF()
	
	target_include_directories(mars PUBLIC ${Kokkos_TPL_INCLUDE_DIRS} ${Kokkos_INCLUDE_DIRS})
	# target_include_directories(mars_exec PUBLIC ${Kokkos_TPL_INCLUDE_DIRS} ${Kokkos_INCLUDE_DIRS})

	message(STATUS "Kokkos_INCLUDE_DIRS=${Kokkos_INCLUDE_DIRS}, ${Kokkos_TPL_INCLUDE_DIRS}")
	
	IF(Kokkos_CXX_COMPILER)
		target_link_libraries(mars_exec mars ${Kokkos_LIBRARIES} ${Kokkos_TPL_LIBRARIES})
	ELSE()
		target_link_libraries(mars_exec mars ${Kokkos_LIBRARIES} ${Kokkos_TPL_LIBRARIES} -L${Kokkos_LIBRARY_DIRS})
	ENDIF()
ENDIF()

if(TRY_WITH_MOONOLITH)
	add_subdirectory(moonolith_adapter)
endif()

target_link_libraries(mars ${MARS_LIBRARIES})
target_include_directories(mars PUBLIC ${CMAKE_BINARY_DIR})
target_include_directories(mars PUBLIC ${MARS_INCLUDES})

message(STATUS "${MARS_INCLUDES}")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mars_config.hpp.in ${CMAKE_BINARY_DIR}/mars_config.hpp)
include_directories(${CMAKE_BINARY_DIR})

install(TARGETS mars mars_exec
	RUNTIME DESTINATION bin 
	LIBRARY DESTINATION lib 
	ARCHIVE DESTINATION lib 
	)

install(FILES ${MARS_HEADERS} DESTINATION include)
install(FILES ${CMAKE_BINARY_DIR}/mars_config.hpp DESTINATION include)



# add_definitions("-DWITH_MOONOLITH")
# set(MOONOLITH_DIR "$ENV{INSTALL_DIR}/moon")
# target_include_directories(mars_exec PUBLIC ${MOONOLITH_DIR}/include)
# target_link_libraries(mars_exec 
# 	-L${MOONOLITH_DIR}/lib/
# 	-leps
# 	-lexternal
# 	-lintersection
# 	-lio
# 	-lmeshing
# 	-lmoonolith
# 	-lquadrature
# 	-lshape
# 	-lutility
# )
