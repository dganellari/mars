# backend/distributed/testsuite/comm/CMakeLists.txt

mars_add_test(morton_test ${CMAKE_CURRENT_SOURCE_DIR}/mars_morton.hpp ${CMAKE_CURRENT_SOURCE_DIR}/mars_morton.hpp)
mars_add_test(hilbert_test ${CMAKE_CURRENT_SOURCE_DIR}/mars_hilbert.hpp ${CMAKE_CURRENT_SOURCE_DIR}/mars_hilbert.hpp)

