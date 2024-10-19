#include <gtest/gtest.h>
#include "mars_template_mesh_manager.hpp"

using namespace mars;

// Test the mesh generation for 2D meshes for the morton key
TEST(MarsTest, TestMeshGeneration2D_Small) {
    const int x = 10;
    const int y = 10;
    auto num_elements = test_mars_distributed_nonsimplex_mesh_generation_kokkos_2D(x, y);
    auto expected_num_elements = x * y;  // Assuming the mesh is a grid of x by y elements
                                             //
    std::cout << "num_elements: " << num_elements << std::endl;
    ASSERT_EQ(num_elements, expected_num_elements);
}

TEST(MarsTest, TestMeshGeneration2D_Medium) {
    const int x = 100;
    const int y = 100;
    auto num_elements = test_mars_distributed_nonsimplex_mesh_generation_kokkos_2D(x, y);
    auto expected_num_elements = x * y;
    ASSERT_EQ(num_elements, expected_num_elements);
}

TEST(MarsTest, TestMeshGeneration2D_Large) {
    const int x = 1000;
    const int y = 1000;
    auto num_elements = test_mars_distributed_nonsimplex_mesh_generation_kokkos_2D(x, y);
    auto expected_num_elements = x * y;
    ASSERT_EQ(num_elements, expected_num_elements);
}

TEST(MarsTest, TestMeshGeneration2D_NonSquare) {
    const int x = 100;
    const int y = 200;
    auto num_elements = test_mars_distributed_nonsimplex_mesh_generation_kokkos_2D(x, y);
    auto expected_num_elements = x * y;
    ASSERT_EQ(num_elements, expected_num_elements);
}

TEST(MarsTest, TestMeshGeneration3D_Small) {
    const int x = 10;
    const int y = 10;
    const int z = 10;
    auto num_elements = test_mars_distributed_nonsimplex_mesh_generation_kokkos_3D(x, y, z);
    auto expected_num_elements = x * y * z;
    ASSERT_EQ(num_elements, expected_num_elements);
}

TEST(MarsTest, TestMeshGeneration3D_Medium) {
    const int x = 100;
    const int y = 100;
    const int z = 100;
    auto num_elements = test_mars_distributed_nonsimplex_mesh_generation_kokkos_3D(x, y, z);
    auto expected_num_elements = x * y * z;
    ASSERT_EQ(num_elements, expected_num_elements);
}

TEST(MarsTest, TestMeshGeneration3D_Large) {
    const int x = 1000;
    const int y = 1000;
    const int z = 1000;
    auto num_elements = test_mars_distributed_nonsimplex_mesh_generation_kokkos_3D(x, y, z);
    auto expected_num_elements = x * y * z;
    ASSERT_EQ(num_elements, expected_num_elements);
}

TEST(MarsTest, TestMeshGeneration3D_NonCube) {
    const int x = 100;
    const int y = 200;
    const int z = 300;
    auto num_elements = test_mars_distributed_nonsimplex_mesh_generation_kokkos_3D(x, y, z);
    auto expected_num_elements = x * y * z;
    ASSERT_EQ(num_elements, expected_num_elements);
}

// Test the mesh generation for 2D meshes for the hilbert key
TEST(MarsTest, TestMeshGeneration2D_HilbertKey) {
    const int x = 10;
    const int y = 10;
    auto num_elements = test_mars_distributed_nonsimplex_mesh_generation_kokkos_2D<HilbertKey<Unsigned>>(x, y);
    auto expected_num_elements = x * y;  // Assuming the mesh is a grid of x by y elements
    ASSERT_EQ(num_elements, expected_num_elements);
}

TEST(MarsTest, TestMeshGeneration2D_HilbertKey_Medium) {
    const int x = 100;
    const int y = 100;
    auto num_elements = test_mars_distributed_nonsimplex_mesh_generation_kokkos_2D<HilbertKey<Unsigned>>(x, y);
    auto expected_num_elements = x * y;
    ASSERT_EQ(num_elements, expected_num_elements);
}

TEST(MarsTest, TestMeshGeneration2D_HilbertKey_Large) {
    const int x = 1000;
    const int y = 1000;
    auto num_elements = test_mars_distributed_nonsimplex_mesh_generation_kokkos_2D<HilbertKey<Unsigned>>(x, y);
    auto expected_num_elements = x * y;
    ASSERT_EQ(num_elements, expected_num_elements);
}

TEST(MarsTest, TestMeshGeneration2D_HilbertKey_NonSquare) {
    const int x = 100;
    const int y = 200;
    auto num_elements = test_mars_distributed_nonsimplex_mesh_generation_kokkos_2D<HilbertKey<Unsigned>>(x, y);
    auto expected_num_elements = x * y;
    ASSERT_EQ(num_elements, expected_num_elements);
}

// Test the mesh generation for 3D meshes for the hilbert key
TEST(MarsTest, TestMeshGeneration3D_HilbertKey_Small) {
    const int x = 10;
    const int y = 10;
    const int z = 10;
    auto num_elements = test_mars_distributed_nonsimplex_mesh_generation_kokkos_3D<HilbertKey<Unsigned>>(x, y, z);
    auto expected_num_elements = x * y * z;
    ASSERT_EQ(num_elements, expected_num_elements);
}

TEST(MarsTest, TestMeshGeneration3D_HilbertKey_Medium) {
    const int x = 100;
    const int y = 100;
    const int z = 100;
    auto num_elements = test_mars_distributed_nonsimplex_mesh_generation_kokkos_3D<HilbertKey<Unsigned>>(x, y, z);
    auto expected_num_elements = x * y * z;
    ASSERT_EQ(num_elements, expected_num_elements);
}

TEST(MarsTest, TestMeshGeneration3D_HilbertKey_Large) {
    const int x = 1000;
    const int y = 1000;
    const int z = 1000;
    auto num_elements = test_mars_distributed_nonsimplex_mesh_generation_kokkos_3D<HilbertKey<Unsigned>>(x, y, z);
    auto expected_num_elements = x * y * z;
    ASSERT_EQ(num_elements, expected_num_elements);
}

TEST(MarsTest, TestMeshGeneration3D_HilbertKey_NonCube) {
    const int x = 100;
    const int y = 200;
    const int z = 300;
    auto num_elements = test_mars_distributed_nonsimplex_mesh_generation_kokkos_3D<HilbertKey<Unsigned>>(x, y, z);
    auto expected_num_elements = x * y * z;
    ASSERT_EQ(num_elements, expected_num_elements);
}

int main(int argc, char **argv) {
    mars::Env env(argc, argv);
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

