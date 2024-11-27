#include <gtest/gtest.h>
#include "mars_template_mesh_manager.hpp"

using namespace mars;

#ifdef KOKKOS_ENABLE_CUDA
// Test the mesh generation for 2D meshes for the morton key
TEST(GhostLayer, MeshGhostLayer2D_Small) {
    const int x = 10;
    const int y = 10;
    ASSERT_TRUE(test_mars_distributed_nonsimplex_mesh_ghost_layer_2D(x, y));
}

TEST(GhostLayer, MeshGhostLayer2D_Medium) {
    const int x = 100;
    const int y = 100;
    ASSERT_TRUE(test_mars_distributed_nonsimplex_mesh_ghost_layer_2D(x, y));
}

TEST(GhostLayer, MeshGhostLayer2D_Large) {
    const int x = 1000;
    const int y = 1000;
    ASSERT_TRUE(test_mars_distributed_nonsimplex_mesh_ghost_layer_2D(x, y));
}

TEST(GhostLayer, MeshGhostLayer2D_NonSquare) {
    const int x = 100;
    const int y = 200;
    ASSERT_TRUE(test_mars_distributed_nonsimplex_mesh_ghost_layer_2D(x, y));
}

TEST(GhostLayer, MeshGhostLayer3D_Small) {
    const int x = 10;
    const int y = 10;
    const int z = 10;
    ASSERT_TRUE(test_mars_distributed_nonsimplex_mesh_ghost_layer_3D(x, y, z));
}

TEST(GhostLayer, MeshGhostLayer3D_Medium) {
    const int x = 100;
    const int y = 100;
    const int z = 100;
    ASSERT_TRUE(test_mars_distributed_nonsimplex_mesh_ghost_layer_3D(x, y, z));
}

TEST(GhostLayer, MeshGhostLayer3D_Large) {
    const int x = 1000;
    const int y = 1000;
    const int z = 1000;
    ASSERT_TRUE(test_mars_distributed_nonsimplex_mesh_ghost_layer_3D(x, y, z));
}


TEST(GhostLayer, MeshGhostLayer3D_NonCube) {
    const int x = 100;
    const int y = 200;
    const int z = 300;
    ASSERT_TRUE(test_mars_distributed_nonsimplex_mesh_ghost_layer_3D(x, y, z));
}

// Test the mesh generation for 2D meshes for the hilbert key
TEST(GhostLayer, MeshGhostLayer2D_HilbertKey) {
    const int x = 10;
    const int y = 10;
    ASSERT_TRUE(test_mars_distributed_nonsimplex_mesh_ghost_layer_2D<HilbertKey<Unsigned>>(x, y));
}

TEST(GhostLayer, MeshGhostLayer2D_HilbertKey_Medium) {
    const int x = 100;
    const int y = 100;
    ASSERT_TRUE(test_mars_distributed_nonsimplex_mesh_ghost_layer_2D<HilbertKey<Unsigned>>(x, y));
}

TEST(GhostLayer, MeshGhostLayer2D_HilbertKey_Large) {
    const int x = 1000;
    const int y = 1000;
    ASSERT_TRUE(test_mars_distributed_nonsimplex_mesh_ghost_layer_2D<HilbertKey<Unsigned>>(x, y));
}

TEST(GhostLayer, MeshGhostLayer2D_HilbertKey_NonSquare) {
    const int x = 100;
    const int y = 200;
    ASSERT_TRUE(test_mars_distributed_nonsimplex_mesh_ghost_layer_2D<HilbertKey<Unsigned>>(x, y));
}

// Test the mesh generation for 3D meshes for the hilbert key
TEST(GhostLayer, MeshGhostLayer3D_HilbertKey_Small) {
    const int x = 10;
    const int y = 10;
    const int z = 10;
    ASSERT_TRUE(test_mars_distributed_nonsimplex_mesh_ghost_layer_3D<HilbertKey<Unsigned>>(x, y, z));
}

TEST(GhostLayer, MeshGhostLayer3D_HilbertKey_Medium) {
    const int x = 100;
    const int y = 100;
    const int z = 100;
    ASSERT_TRUE(test_mars_distributed_nonsimplex_mesh_ghost_layer_3D<HilbertKey<Unsigned>>(x, y, z));
}

TEST(GhostLayer, MeshGhostLayer3D_HilbertKey_Large) {
    const int x = 1000;
    const int y = 1000;
    const int z = 1000;
    ASSERT_TRUE(test_mars_distributed_nonsimplex_mesh_ghost_layer_3D<HilbertKey<Unsigned>>(x, y, z));
}


TEST(GhostLayer, MeshGhostLayer3D_HilbertKey_NonCube) {
    const int x = 100;
    const int y = 200;
    const int z = 300;
    ASSERT_TRUE(test_mars_distributed_nonsimplex_mesh_ghost_layer_3D<HilbertKey<Unsigned>>(x, y, z));
}

TEST(GhostLayer, MeshGhostLayer3D_HilbertKey_Large_power_of_2) {
    const int x = 1024;
    const int y = 256;
    const int z = 512;
    ASSERT_TRUE(test_mars_distributed_nonsimplex_mesh_ghost_layer_3D<HilbertKey<Unsigned>>(x, y, z));
}

#endif

int main(int argc, char **argv) {
    mars::Env env(argc, argv);
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

