# backend/distributed/testsuite/comm/CMakeLists.txt

mars_add_test(mesh_manager_test ${CMAKE_CURRENT_SOURCE_DIR}/test_mesh_generation.hpp ${CMAKE_CURRENT_SOURCE_DIR}/test_mesh_generation.hpp)

mars_add_test(mesh_discretization ${CMAKE_CURRENT_SOURCE_DIR}/test_discretization.hpp ${CMAKE_CURRENT_SOURCE_DIR}/test_discretization.hpp)

