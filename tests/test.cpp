#include "gtest/gtest.h"

#include "mars_instance.hpp"

int main(int argc, char **argv) {
    using namespace mars;

    MARS::init(argc, argv);

    ::testing::InitGoogleTest(&argc, argv);
    int ok = RUN_ALL_TESTS();

    MARS::finalize();
    return ok;
}
