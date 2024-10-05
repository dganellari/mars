// File: testsuite/comm/test_mars_mpi.cpp

#include <gtest/gtest.h>
#include "mars.hpp"

TEST(MarsTest, TestMPI) {
    int argc = 1;
    char *argv[] = { (char*)"test_program", NULL };
    testing::internal::CaptureStdout();
    try {
        mars::marsenv::mpi_guard guard(argc, argv, false);
        if (mars::rank(context) == 0) {
            printf("MPI initialized\n");
        }
        std::string output = testing::internal::GetCapturedStdout();
        if (mars::rank(context) == 0) {
            EXPECT_EQ(output, "MPI initialized\n");
        }
    } catch (std::exception &e) {
        FAIL() << "Exception caught: " << e.what();
    }
}

TEST(MarsTest, TestMPIContext) {
    int argc = 1;
    char *argv[] = { (char*)"test_program", NULL };
    testing::internal::CaptureStdout();
    try {
        mars::proc_allocation resources;
        mars::marsenv::mpi_guard guard(argc, argv, false);
        auto context = mars::make_context(resources, MPI_COMM_WORLD);
        if (mars::rank(context) == 0) {
            std::cout << "mpi:      " << (mars::has_mpi(context) ? "yes" : "no") << "\n";
            std::cout << "ranks:    " << mars::num_ranks(context) << "\n";
            std::cout << "rank:    " << mars::rank(context) << "\n" << std::endl;
        }
        std::string output = testing::internal::GetCapturedStdout();
        if (mars::rank(context) == 0) {
            std::string expected_output = "mpi:      yes\nranks:    " + std::to_string(mars::num_ranks(context)) + "\nrank:    0\n\n";
            EXPECT_EQ(output, expected_output);
        }
    } catch (std::exception &e) {
        FAIL() << "Exception caught: " << e.what();
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
