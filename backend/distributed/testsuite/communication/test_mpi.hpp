// File: testsuite/comm/test_mars_mpi.cpp

#include <gtest/gtest.h>
#include "mars_test_mpi.hpp"



TEST(MarsTest, TestMPI) {
    int argc = 1;
    char *argv[] = { (char*)"test_program", NULL };

    auto result = mars::test_mpi(argc, argv);
    EXPECT_TRUE(result);
}

TEST(MarsTest, TestMPIContext) {
    int argc = 1;
    char *argv[] = {(char *)"test_program", NULL};
    testing::internal::CaptureStdout();
    mars::test_mpi_context(argc, argv);
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
