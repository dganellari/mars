#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include "mars_template_mesh_manager.hpp"

using namespace mars;

TEST(DofHandlerTest, TestDOFHandlerEnumerationMorton2D) {
    // Define the dimensions of the mesh
    const int xDim = 10;
    const int yDim = 10;
    const int zDim = 0;
    const int block = 1;

    auto unique = test_mars_distributed_dof_handler<ElementType::Quad4>(xDim, yDim, zDim, block);
    // Use a Kokkos parallel_reduce to check for unique identifiers
    ASSERT_TRUE(unique);
}

TEST(DofHandlerTest, TestDOFHandlerEnumerationMorton2DDegree2) {
    // Define the dimensions of the mesh
    const int xDim = 10;
    const int yDim = 10;
    const int zDim = 0;
    const int block = 1;

    auto unique = test_mars_distributed_dof_handler<ElementType::Quad4, 2>(xDim, yDim, zDim, block);
    // Use a Kokkos parallel_reduce to check for unique identifiers
    ASSERT_TRUE(unique);
}

TEST(DofHandlerTest, TestDOFHandlerEnumerationMorton3D) {
    // Define the dimensions of the mesh
    const int xDim = 10;
    const int yDim = 10;
    const int zDim = 10;
    const int block = 1;

    auto unique = test_mars_distributed_dof_handler<ElementType::Hex8>(xDim, yDim, zDim, block);
    // Use a Kokkos parallel_reduce to check for unique identifiers
    ASSERT_TRUE(unique);
}

TEST(DofHandlerTest, TestDOFHandlerEnumerationMorton3DDegree2) {
    // Define the dimensions of the mesh
    const int xDim = 10;
    const int yDim = 10;
    const int zDim = 10;
    const int block = 1;

    auto unique = test_mars_distributed_dof_handler<ElementType::Hex8, 2>(xDim, yDim, zDim, block);
    // Use a Kokkos parallel_reduce to check for unique identifiers
    ASSERT_TRUE(unique);
}

TEST(DofHandlerTest, TestDOFHandlerEnumerationHilbert2D) {
    // Define the dimensions of the mesh
    const int xDim = 10;
    const int yDim = 10;
    const int zDim = 0;
    const int block = 1;

    auto unique = test_mars_distributed_dof_handler<ElementType::Quad4, 1, true, HilbertKey<Unsigned>>(xDim, yDim, zDim, block);
    // Use a Kokkos parallel_reduce to check for unique identifiers
    ASSERT_TRUE(unique);
}

TEST(DofHandlerTest, TestDOFHandlerEnumerationHilbert2DDegree2) {
    // Define the dimensions of the mesh
    const int xDim = 10;
    const int yDim = 10;
    const int zDim = 0;
    const int block = 1;

    auto unique = test_mars_distributed_dof_handler<ElementType::Quad4, 2, true, HilbertKey<Unsigned>>(xDim, yDim, zDim, block);
    // Use a Kokkos parallel_reduce to check for unique identifiers
    ASSERT_TRUE(unique);
}

TEST(DofHandlerTest, TestDOFHandlerEnumerationHilbert3D) {
    // Define the dimensions of the mesh
    const int xDim = 10;
    const int yDim = 10;
    const int zDim = 10;
    const int block = 1;

    auto unique = test_mars_distributed_dof_handler<ElementType::Hex8, 1, true, HilbertKey<Unsigned>>(xDim, yDim, zDim, block);
    // Use a Kokkos parallel_reduce to check for unique identifiers
    ASSERT_TRUE(unique);
}


TEST(DofHandlerTest, TestDOFHandlerEnumerationHilbert3DDegree2) {
    // Define the dimensions of the mesh
    const int xDim = 10;
    const int yDim = 10;
    const int zDim = 10;
    const int block = 1;

    auto unique = test_mars_distributed_dof_handler<ElementType::Hex8, 2, true, HilbertKey<Unsigned>>(xDim, yDim, zDim, block);
    // Use a Kokkos parallel_reduce to check for unique identifiers
    ASSERT_TRUE(unique);
}

int main(int argc, char **argv) {
    mars::Env env(argc, argv);
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

