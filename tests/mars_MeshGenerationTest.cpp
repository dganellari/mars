
#include "gtest/gtest.h"
#include "mars.hpp"

namespace mars {
    TEST(Mesh4, generate) {
        Mesh4 m;
        EXPECT_TRUE(m.n_elements() == 0);
    }
}  // namespace mars
