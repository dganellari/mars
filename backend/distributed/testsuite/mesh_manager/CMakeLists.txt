# backend/distributed/testsuite/comm/CMakeLists.txt

mars_add_test(mesh_manager_test ${CMAKE_CURRENT_SOURCE_DIR}/mars_mesh_generation.cpp)
mars_add_test(mesh_finite_elements_test ${CMAKE_CURRENT_SOURCE_DIR}/mars_finite_elements.cpp)
# mars_add_test(mesh_discretization_test ${CMAKE_CURRENT_SOURCE_DIR}/mars_discretization.cpp)
