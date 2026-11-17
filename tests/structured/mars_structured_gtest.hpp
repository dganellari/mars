#pragma once

#include "mars_structured_runtime.hpp"
#include <gtest/gtest.h>

namespace mars::test
{

inline int run_structured_tests(int& argc, char**& argv, int expected_cases)
{
    StructuredRuntime runtime(argc, argv);
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    const auto* tests = testing::UnitTest::GetInstance();
    int rank = 0, ranks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    int valid = result == 0 && tests->test_to_run_count() == expected_cases &&
                tests->successful_test_count() == expected_cases && tests->skipped_test_count() == 0;
    MPI_Allreduce(MPI_IN_PLACE, &valid, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    std::ostringstream receipt;
    receipt << "MARS_STRUCTURED_TEST_RECEIPT {\"rank\":" << rank << ",\"ranks\":" << ranks
            << ",\"cases\":" << tests->test_to_run_count()
            << ",\"passed\":" << tests->successful_test_count()
            << ",\"skipped\":" << tests->skipped_test_count() << "}";
    for (int writer = 0; writer < ranks; ++writer)
    {
        if (writer == rank) std::cout << receipt.str() << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
    }
    if (!valid && rank == 0)
        std::cerr << "The structured gate requires every case on every rank; filtering or skipping fails the gate.\n";
    return valid ? 0 : 1;
}

} // namespace mars::test
