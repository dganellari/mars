if(EXISTS "/Users/gandanie/scratch/santis/mars/mac/backend/distributed/testsuite/communication/mpi_test[1]_tests.cmake")
  include("/Users/gandanie/scratch/santis/mars/mac/backend/distributed/testsuite/communication/mpi_test[1]_tests.cmake")
else()
  add_test(mpi_test_NOT_BUILT mpi_test_NOT_BUILT)
endif()
