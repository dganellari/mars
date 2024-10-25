#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include "mars_template_mesh_manager.hpp"

using namespace mars;

TEST(DofHandler, EnumerationMorton2D) {
    // Define the dimensions of the mesh
    const int xDim = 64;
    const int yDim = 64;
    const int zDim = 0;
    const int block = 1;

    auto unique = test_mars_distributed_dof_handler<ElementType::Quad4>(xDim, yDim, zDim, block);
    // Use a Kokkos parallel_reduce to check for unique identifiers
    ASSERT_TRUE(unique);
}

TEST(DofHandler, EnumerationMorton2DDegree2) {
    // Define the dimensions of the mesh
    const int xDim = 128;
    const int yDim = 256;
    const int zDim = 0;
    const int block = 1;

    auto unique = test_mars_distributed_dof_handler<ElementType::Quad4, 2>(xDim, yDim, zDim, block);
    // Use a Kokkos parallel_reduce to check for unique identifiers
    ASSERT_TRUE(unique);
}

TEST(DofHandler, EnumerationMorton3D) {
    // Define the dimensions of the mesh
    const int xDim = 290;
    const int yDim = 420;
    const int zDim = 210;
    const int block = 1;

    auto unique = test_mars_distributed_dof_handler<ElementType::Hex8>(xDim, yDim, zDim, block);
    // Use a Kokkos parallel_reduce to check for unique identifiers
    ASSERT_TRUE(unique);
}

TEST(DofHandler, EnumerationMorton3DDegree2) {
    // Define the dimensions of the mesh
    const int xDim = 290;
    const int yDim = 420;
    const int zDim = 210;
    const int block = 1;

    auto unique = test_mars_distributed_dof_handler<ElementType::Hex8, 2>(xDim, yDim, zDim, block);
    // Use a Kokkos parallel_reduce to check for unique identifiers
    ASSERT_TRUE(unique);
}

TEST(DofHandler, EnumerationHilbert2D) {
    // Define the dimensions of the mesh
    const int xDim = 4;
    const int yDim = 4;
    const int zDim = 0;
    const int block = 1;

    auto unique = test_mars_distributed_dof_handler<ElementType::Quad4, 1, true, HilbertKey<Unsigned>>(xDim, yDim, zDim, block);
    // Use a Kokkos parallel_reduce to check for unique identifiers
    ASSERT_TRUE(unique);
}

TEST(DofHandler, EnumerationHilbert2DDegree2) {
    // Define the dimensions of the mesh
    const int xDim = 128;
    const int yDim = 128;
    const int zDim = 0;
    const int block = 1;

    auto unique = test_mars_distributed_dof_handler<ElementType::Quad4, 2, true, HilbertKey<Unsigned>>(xDim, yDim, zDim, block);
    // Use a Kokkos parallel_reduce to check for unique identifiers
    ASSERT_TRUE(unique);
}

TEST(DofHandler, EnumerationHilbert3D) {
    // Define the dimensions of the mesh
    const int xDim = 64;
    const int yDim = 64;
    const int zDim = 64;
    const int block = 1;

    auto unique = test_mars_distributed_dof_handler<ElementType::Hex8, 1, true, HilbertKey<Unsigned>>(xDim, yDim, zDim, block);
    // Use a Kokkos parallel_reduce to check for unique identifiers
    ASSERT_TRUE(unique);
}


TEST(DofHandler, EnumerationHilbert3DDegree2) {
    // Define the dimensions of the mesh
    const int xDim = 32;
    const int yDim = 32;
    const int zDim = 32;
    const int block = 1;

    auto unique = test_mars_distributed_dof_handler<ElementType::Hex8, 2, true, HilbertKey<Unsigned>>(xDim, yDim, zDim, block);
    // Use a Kokkos parallel_reduce to check for unique identifiers
    ASSERT_TRUE(unique);
}

// Vector Valued Tests

TEST(DofHandler, BlockEnumerationMorton2D) {
   // Define the dimensions of the mesh
   const int xDim = 64;
   const int yDim = 64;
   const int zDim = 0;
   const int block = 2;
   auto unique = test_mars_distributed_dof_handler<ElementType::Quad4>(xDim, yDim, zDim, block);
   // Use a Kokkos parallel_reduce to check for unique identifiers
   ASSERT_TRUE(unique);
   }

TEST(DofHandler, BlockEnumerationMorton2DDegree2) {
    // Define the dimensions of the mesh
    const int xDim = 128;
    const int yDim = 256;
    const int zDim = 0;
    const int block = 2;
    auto unique = test_mars_distributed_dof_handler<ElementType::Quad4, 2>(xDim, yDim, zDim, block);
    // Use a Kokkos parallel_reduce to check for unique identifiers
    ASSERT_TRUE(unique);
    }

TEST(DofHandler, BlockEnumerationMorton3D) {
    // Define the dimensions of the mesh
    const int xDim = 290;
    const int yDim = 420;
    const int zDim = 210;
    const int block = 2;
    auto unique = test_mars_distributed_dof_handler<ElementType::Hex8>(xDim, yDim, zDim, block);
    // Use a Kokkos parallel_reduce to check for unique identifiers
    ASSERT_TRUE(unique);
    }

TEST(DofHandler, BlockEnumerationMorton3DDegree2) {
    // Define the dimensions of the mesh
    const int xDim = 290;
    const int yDim = 420;
    const int zDim = 210;
    const int block = 2;
    auto unique = test_mars_distributed_dof_handler<ElementType::Hex8, 2>(xDim, yDim, zDim, block);
    // Use a Kokkos parallel_reduce to check for unique identifiers
    ASSERT_TRUE(unique);
    }

TEST(DofHandler, BlockEnumerationHilbert2D) {
    // Define the dimensions of the mesh
    const int xDim = 256;
    const int yDim = 256;
    const int zDim = 0;
    const int block = 2;
    auto unique = test_mars_distributed_dof_handler<ElementType::Quad4, 1, true, HilbertKey<Unsigned>>(xDim, yDim, zDim, block);
    // Use a Kokkos parallel_reduce to check for unique identifiers
    ASSERT_TRUE(unique);
    }

TEST(DofHandler, BlockEnumerationHilbert2DDegree2) {
    // Define the dimensions of the mesh
    const int xDim = 128;
    const int yDim = 128;
    const int zDim = 0;
    const int block = 2;
    auto unique = test_mars_distributed_dof_handler<ElementType::Quad4, 2, true, HilbertKey<Unsigned>>(xDim, yDim, zDim, block);
    // Use a Kokkos parallel_reduce to check for unique identifiers
    ASSERT_TRUE(unique);
    }

TEST(DofHandler, BlockEnumerationHilbert3D) {
    // Define the dimensions of the mesh
    const int xDim = 64;
    const int yDim = 64;
    const int zDim = 64;
    const int block = 2;
    auto unique = test_mars_distributed_dof_handler<ElementType::Hex8, 1, true, HilbertKey<Unsigned>>(xDim, yDim, zDim, block);
    // Use a Kokkos parallel_reduce to check for unique identifiers
    ASSERT_TRUE(unique);
    }


TEST(DofHandler, BlockEnumerationHilbert3DDegree2) {
    // Define the dimensions of the mesh
    const int xDim = 32;
    const int yDim = 32;
    const int zDim = 32;
    const int block = 2;
    auto unique = test_mars_distributed_dof_handler<ElementType::Hex8, 2, true, HilbertKey<Unsigned>>(xDim, yDim, zDim, block);
    // Use a Kokkos parallel_reduce to check for unique identifiers
    ASSERT_TRUE(unique);
    }

int main(int argc, char **argv) {
    mars::Env env(argc, argv);
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

