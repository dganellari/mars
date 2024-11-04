#include <gtest/gtest.h>
#include "mars_template_mesh_manager.hpp"

using namespace mars;

TEST(FEDofMapTest, TopologicalEnumeration2DMorton) {
    bool success = test_mars_distributed_dof_map<ElementType::Quad4, 1, true, MortonKey<Unsigned>>(4, 4, 0, 1);
    ASSERT_TRUE(success);
}

TEST(FEDofMapTest, TopologicalEnumeration3DMorton) {
    bool success = test_mars_distributed_dof_map<ElementType::Hex8, 1, true, MortonKey<Unsigned>>(4, 4, 4, 1);
    ASSERT_TRUE(success);
}

TEST(FEDofMapTest, TopologicalEnumeration2DVectorValuedMorton) {
    bool success = test_mars_distributed_dof_map<ElementType::Quad4, 1, true, MortonKey<Unsigned>>(4, 4, 0, 2);
    ASSERT_TRUE(success);
}

TEST(FEDofMapTest, TopologicalEnumeration3DVectorValuedMorton) {
    bool success = test_mars_distributed_dof_map<ElementType::Hex8, 1, true, MortonKey<Unsigned>>(4, 4, 4, 2);
    ASSERT_TRUE(success);
}

TEST(FEDofMapTest, TopologicalEnumeration2DLargeMorton) {
    bool success = test_mars_distributed_dof_map<ElementType::Quad4, 1, true, MortonKey<Unsigned>>(64, 64, 0, 1);
    ASSERT_TRUE(success);
}

TEST(FEDofMapTest, TopologicalEnumeration3DLargeMorton) {
    bool success = test_mars_distributed_dof_map<ElementType::Hex8, 1, true, MortonKey<Unsigned>>(32, 32, 32, 1);
    ASSERT_TRUE(success);
}

TEST(FEDofMapTest, EdgeCasesMorton) {
    bool success = test_mars_distributed_dof_map<ElementType::Quad4, 1, true, MortonKey<Unsigned>>(1, 1, 0, 1);
    ASSERT_TRUE(success);
}

TEST(FEDofMapTest, TopologicalEnumeration2DHilbert) {
    bool success = test_mars_distributed_dof_map<ElementType::Quad4, 1, true, HilbertKey<Unsigned>>(4, 4, 0, 1);
    ASSERT_TRUE(success);
}

TEST(FEDofMapTest, TopologicalEnumeration3DHilbert) {
    bool success = test_mars_distributed_dof_map<ElementType::Hex8, 1, true, HilbertKey<Unsigned>>(4, 4, 4, 1);
    ASSERT_TRUE(success);
}

TEST(FEDofMapTest, TopologicalEnumeration2DVectorValuedHilbert) {
    bool success = test_mars_distributed_dof_map<ElementType::Quad4, 1, true, HilbertKey<Unsigned>>(4, 4, 0, 2);
    ASSERT_TRUE(success);
}

TEST(FEDofMapTest, TopologicalEnumeration3DVectorValuedHilbert) {
    bool success = test_mars_distributed_dof_map<ElementType::Hex8, 1, true, HilbertKey<Unsigned>>(4, 4, 4, 2);
    ASSERT_TRUE(success);
}

TEST(FEDofMapTest, TopologicalEnumeration2DLargeHilbert) {
    bool success = test_mars_distributed_dof_map<ElementType::Quad4, 1, true, HilbertKey<Unsigned>>(64, 64, 0, 1);
    ASSERT_TRUE(success);
}

TEST(FEDofMapTest, TopologicalEnumeration3DLargeHilbert) {
    bool success = test_mars_distributed_dof_map<ElementType::Hex8, 1, true, HilbertKey<Unsigned>>(32, 32, 32, 1);
    ASSERT_TRUE(success);
}

TEST(FEDofMapTest, EdgeCasesHilbert) {
    bool success = test_mars_distributed_dof_map<ElementType::Quad4, 1, true, HilbertKey<Unsigned>>(1, 1, 0, 1);
    ASSERT_TRUE(success);
}

// Test sparsity pattern

TEST(SparsityPatternTest, SparsityPattern2DMorton) {
    bool success = test_mars_distributed_sparsity_pattern<ElementType::Quad4, 1, true, MortonKey<Unsigned>>(4, 4, 0, 1);
    ASSERT_TRUE(success);
}

TEST(SparsityPatternTest, SparsityPattern3DMorton) {
    bool success = test_mars_distributed_sparsity_pattern<ElementType::Hex8, 1, true, MortonKey<Unsigned>>(4, 4, 4, 1);
    ASSERT_TRUE(success);
}

TEST(SparsityPatternTest, SparsityPattern2DVectorValuedMorton) {
    bool success = test_mars_distributed_sparsity_pattern<ElementType::Quad4, 1, true, MortonKey<Unsigned>>(4, 4, 0, 2);
    ASSERT_TRUE(success);
}

TEST(SparsityPatternTest, SparsityPattern3DVectorValuedMorton) {
    bool success = test_mars_distributed_sparsity_pattern<ElementType::Hex8, 1, true, MortonKey<Unsigned>>(4, 4, 4, 2);
    ASSERT_TRUE(success);
}

TEST(SparsityPatternTest, SparsityPattern2DLargeMorton) {
    bool success = test_mars_distributed_sparsity_pattern<ElementType::Quad4, 1, true, MortonKey<Unsigned>>(64, 64, 0, 1);
    ASSERT_TRUE(success);
}

TEST(SparsityPatternTest, SparsityPattern3DLargeMorton) {
    bool success = test_mars_distributed_sparsity_pattern<ElementType::Hex8, 1, true, MortonKey<Unsigned>>(32, 32, 32, 1);
    ASSERT_TRUE(success);
}

TEST(SparsityPatternTest, SparsityPattern2DHilbert) {
    bool success = test_mars_distributed_sparsity_pattern<ElementType::Quad4, 1, true, HilbertKey<Unsigned>>(4, 4, 0, 1);
    ASSERT_TRUE(success);
}

TEST(SparsityPatternTest, SparsityPattern3DHilbert) {
    bool success = test_mars_distributed_sparsity_pattern<ElementType::Hex8, 1, true, HilbertKey<Unsigned>>(4, 4, 4, 1);
    ASSERT_TRUE(success);
}

TEST(SparsityPatternTest, SparsityPattern2DVectorValuedHilbert) {
    bool success = test_mars_distributed_sparsity_pattern<ElementType::Quad4, 1, true, HilbertKey<Unsigned>>(4, 4, 0, 2);
    ASSERT_TRUE(success);
}

TEST(SparsityPatternTest, SparsityPattern3DVectorValuedHilbert) {
    bool success = test_mars_distributed_sparsity_pattern<ElementType::Hex8, 1, true, HilbertKey<Unsigned>>(4, 4, 4, 2);
    ASSERT_TRUE(success);
}

TEST(SparsityPatternTest, SparsityPattern2DLargeHilbert) {
    bool success = test_mars_distributed_sparsity_pattern<ElementType::Quad4, 1, true, HilbertKey<Unsigned>>(64, 64, 0, 1);
    ASSERT_TRUE(success);
}

TEST(SparsityPatternTest, SparsityPattern3DLargeHilbert) {
    bool success = test_mars_distributed_sparsity_pattern<ElementType::Hex8, 1, true, HilbertKey<Unsigned>>(32, 32, 32, 1);
    ASSERT_TRUE(success);
}

int main(int argc, char **argv) {
    mars::Env env(argc, argv);
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}