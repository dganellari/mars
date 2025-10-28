# CMake generated Testfile for 
# Source directory: /Users/gandanie/scratch/santis/mars/backend/distributed/testsuite/communication
# Build directory: /Users/gandanie/scratch/santis/mars/mac/backend/distributed/testsuite/communication
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
include("/Users/gandanie/scratch/santis/mars/mac/backend/distributed/testsuite/communication/mpi_test[1]_include.cmake")
add_test([=[MPITest.TestGuard]=] "/Users/gandanie/scratch/santis/mars/mac/backend/distributed/testsuite/communication/mpi_test" "--gtest_filter=MPITest.TestGuard")
set_tests_properties([=[MPITest.TestGuard]=] PROPERTIES  SKIP_REGULAR_EXPRESSION "\\[  SKIPPED \\]" _BACKTRACE_TRIPLES "/opt/homebrew/Cellar/cmake/3.28.3/share/cmake/Modules/GoogleTest.cmake;402;add_test;/Users/gandanie/scratch/santis/mars/cmake/Tests.cmake;41;gtest_add_tests;/Users/gandanie/scratch/santis/mars/backend/distributed/testsuite/communication/CMakeLists.txt;3;mars_add_test;/Users/gandanie/scratch/santis/mars/backend/distributed/testsuite/communication/CMakeLists.txt;0;")
add_test([=[MPITest.TestMPIContext]=] "/Users/gandanie/scratch/santis/mars/mac/backend/distributed/testsuite/communication/mpi_test" "--gtest_filter=MPITest.TestMPIContext")
set_tests_properties([=[MPITest.TestMPIContext]=] PROPERTIES  SKIP_REGULAR_EXPRESSION "\\[  SKIPPED \\]" _BACKTRACE_TRIPLES "/opt/homebrew/Cellar/cmake/3.28.3/share/cmake/Modules/GoogleTest.cmake;402;add_test;/Users/gandanie/scratch/santis/mars/cmake/Tests.cmake;41;gtest_add_tests;/Users/gandanie/scratch/santis/mars/backend/distributed/testsuite/communication/CMakeLists.txt;3;mars_add_test;/Users/gandanie/scratch/santis/mars/backend/distributed/testsuite/communication/CMakeLists.txt;0;")
