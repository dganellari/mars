cmake_minimum_required(VERSION 3.8.1)

project(mars)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(MPIExtended)

if(MPI_FOUND)
	set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
	set(CMAKE_CXX_COMPILER_DEBUG ${MPI_CXX_COMPILER})
	set(WITH_MPI ON)
endif()

list(APPEND MARS_HEADERS
	mars.hpp 
	mars_base.hpp
	mars_benchmark.hpp
	mars_bisection.hpp
	mars_communicator.hpp
	mars_dof_map.hpp
	mars_dual_graph.hpp
	mars_edge.hpp
	mars_edge_element_map.hpp
	mars_edge_node_map.hpp
	mars_edge_select.hpp
	mars_edge_split.hpp
	mars_edge_split_pool.hpp
	mars_lagrange_element.hpp
	mars_longest_edge.hpp
	mars_matrix.hpp
	mars_mesh.hpp
	mars_mesh_partition.hpp
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
	)

list(APPEND MARS_SOURCES
	mars_communicator.cpp
	mars_simplex.cpp
	mars_test.cpp
	)

set(CMAKE_CXX_DEBUG "-g")

add_library(mars STATIC ${MARS_HEADERS} ${MARS_SOURCES})
target_link_libraries(mars ${MPI_LIBRARIES})
target_compile_features(mars PUBLIC cxx_std_11)
target_include_directories(mars PUBLIC .)

add_executable(mars_exec main.cpp ${MARS_HEADERS} ${MARS_SOURCES})
target_link_libraries(mars_exec mars)

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