cmake_minimum_required(VERSION 3.8.1)

project(mars)


# find_package(MPI REQUIRED)
# set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
# set(CMAKE_CXX_COMPILER_DEBUG ${MPI_CXX_COMPILER})

add_definitions("-DWITH_MPI")
set(CMAKE_CXX_COMPILER mpicc)
set(CMAKE_CXX_COMPILER_DEBUG mpicxx)
set(MPI_LIBRARIES "-L/opt/local/lib/mpich-mp/ -lmpi -lmpicxx")



list(APPEND MARS_HEADERS 
	base.hpp
	bisection.hpp
	dual_graph.hpp
	edge.hpp
	edge_element_map.hpp
	edge_node_map.hpp
	lagrange_element.hpp
	matrix.hpp
	mesh.hpp
	red_green_refinement.hpp
	simplex.hpp
	static_math.hpp
	test.hpp
	vector.hpp
	)

list(APPEND MARS_SOURCES
	simplex.cpp
	test.cpp
	communicator.cpp
	)

set(CMAKE_CXX_DEBUG "-g")

add_library(mars STATIC ${MARS_HEADERS} ${MARS_SOURCES})
target_link_libraries(mars ${MPI_LIBRARIES})
target_compile_features(mars PUBLIC cxx_std_11)
target_include_directories(mars PUBLIC .)

add_executable(mars_exec main.cpp ${MARS_HEADERS} ${MARS_SOURCES})
target_link_libraries(mars_exec mars)
target_compile_features(mars_exec PUBLIC cxx_std_11)


install(TARGETS mars mars_exec
	RUNTIME DESTINATION bin 
	LIBRARY DESTINATION lib 
	ARCHIVE DESTINATION lib 
)

install(FILES ${MARS_HEADERS} DESTINATION include)

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