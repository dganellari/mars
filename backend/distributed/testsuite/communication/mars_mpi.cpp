// File: testsuite/comm/test_mars_mpi.cpp

#include <gtest/gtest.h>
#include "mars_template_mpi.hpp"



TEST(MPI, TestGuard) {
    int argc = 1;
    char *argv[] = { (char*)"test_program", NULL };

    char**argv_ptr = argv;
    auto result = mars::test_mpi(argc, argv_ptr);
    EXPECT_TRUE(result);
}

TEST(MPI, TestMPIContext) {
    int argc = 1;
    char *argv[] = {(char *)"test_program", NULL};
    char **argv_ptr = argv;
    testing::internal::CaptureStdout();
    auto context = mars::test_mpi_context(argc, argv_ptr);
    std::string output = testing::internal::GetCapturedStdout();
    if (mars::rank(context) == 0) {
        std::string expected_output =
            "mpi:      yes\nranks:    " + std::to_string(mars::num_ranks(context)) + "\nrank:    0\n\n";
        EXPECT_EQ(output, expected_output);
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
