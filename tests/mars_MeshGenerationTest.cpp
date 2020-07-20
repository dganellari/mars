#include "gtest/gtest.h"

#include "mars.hpp"
using namespace mars;

template <typename MeshType>
class MeshGenerationTest : public testing::Test {
public:
    void SetUp() override {
        constexpr Integer n_x_dim = 10;
        for (Integer d = 0; d < MeshType::Dim; ++d) {
            dims[d] = n_x_dim;
        }
    }

protected:
    std::array<Integer, 4> dims{0, 0, 0, 0};
};

using MeshGenerationTestTypes = ::testing::Types<Mesh2, Mesh3>;

TYPED_TEST_SUITE(MeshGenerationTest, MeshGenerationTestTypes);

TYPED_TEST(MeshGenerationTest, Generate) {
    TypeParam mesh;
    generate_cube(mesh, this->dims[0], this->dims[1], this->dims[2]);
    mesh.build_dual_graph();
    EXPECT_TRUE(mesh.check_side_ordering());
}
