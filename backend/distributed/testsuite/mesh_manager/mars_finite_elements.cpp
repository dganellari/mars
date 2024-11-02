#include <gtest/gtest.h>
#include "mars_template_mesh_manager.hpp"

using namespace mars;

template <Integer Type, Integer Degree = 1, bool Overlap = true, class KeyType = MortonKey<Unsigned>>
void run_topological_enumeration_test(const int xDim, const int yDim, const int zDim, const int block) {
    test_mars_distributed_dof_map<Type, Degree, Overlap, KeyType>(xDim, yDim, zDim, block);
}

TEST(FEDofMapTest, TopologicalEnumeration2DMorton) {
    run_topological_enumeration_test<ElementType::Quad4, 1, true, MortonKey<Unsigned>>(4, 4, 0, 1);
}

TEST(FEDofMapTest, TopologicalEnumeration3DMorton) {
    run_topological_enumeration_test<ElementType::Hex8, 1, true, MortonKey<Unsigned>>(4, 4, 4, 1);
}

TEST(FEDofMapTest, TopologicalEnumeration2DVectorValuedMorton) {
    run_topological_enumeration_test<ElementType::Quad4, 1, true, MortonKey<Unsigned>>(4, 4, 0, 2);
}

TEST(FEDofMapTest, TopologicalEnumeration3DVectorValuedMorton) {
    run_topological_enumeration_test<ElementType::Hex8, 1, true, MortonKey<Unsigned>>(4, 4, 4, 2);
}

TEST(FEDofMapTest, TopologicalEnumeration2DLargeMorton) {
    run_topological_enumeration_test<ElementType::Quad4, 1, true, MortonKey<Unsigned>>(64, 64, 0, 1);
}

TEST(FEDofMapTest, TopologicalEnumeration3DLargeMorton) {
    run_topological_enumeration_test<ElementType::Hex8, 1, true, MortonKey<Unsigned>>(32, 32, 32, 1);
}

TEST(FEDofMapTest, EdgeCasesMorton) {
    run_topological_enumeration_test<ElementType::Quad4, 1, true, MortonKey<Unsigned>>(1, 1, 0, 1);
}

TEST(FEDofMapTest, TopologicalEnumeration2DHilbert) {
    run_topological_enumeration_test<ElementType::Quad4, 1, true, HilbertKey<Unsigned>>(4, 4, 0, 1);
}

TEST(FEDofMapTest, TopologicalEnumeration3DHilbert) {
    run_topological_enumeration_test<ElementType::Hex8, 1, true, HilbertKey<Unsigned>>(4, 4, 4, 1);
}

TEST(FEDofMapTest, TopologicalEnumeration2DVectorValuedHilbert) {
    run_topological_enumeration_test<ElementType::Quad4, 1, true, HilbertKey<Unsigned>>(4, 4, 0, 2);
}

TEST(FEDofMapTest, TopologicalEnumeration3DVectorValuedHilbert) {
    run_topological_enumeration_test<ElementType::Hex8, 1, true, HilbertKey<Unsigned>>(4, 4, 4, 2);
}

TEST(FEDofMapTest, TopologicalEnumeration2DLargeHilbert) {
    run_topological_enumeration_test<ElementType::Quad4, 1, true, HilbertKey<Unsigned>>(64, 64, 0, 1);
}

TEST(FEDofMapTest, TopologicalEnumeration3DLargeHilbert) {
    run_topological_enumeration_test<ElementType::Hex8, 1, true, HilbertKey<Unsigned>>(32, 32, 32, 1);
}

TEST(FEDofMapTest, EdgeCasesHilbert) {
    run_topological_enumeration_test<ElementType::Quad4, 1, true, HilbertKey<Unsigned>>(1, 1, 0, 1);
}

//Test sparsity pattern

TEST(SparsityPatternTest, SparsityPattern2DMorton) {
    test_mars_distributed_sparsity_pattern<ElementType::Quad4, 1, true, MortonKey<Unsigned>>(4, 4, 0, 1);
}

TEST(SparsityPatternTest, SparsityPattern3DMorton) {
    test_mars_distributed_sparsity_pattern<ElementType::Hex8, 1, true, MortonKey<Unsigned>>(4, 4, 4, 1);
}

TEST(SparsityPatternTest, SparsityPattern2DVectorValuedMorton) {
    test_mars_distributed_sparsity_pattern<ElementType::Quad4, 1, true, MortonKey<Unsigned>>(4, 4, 0, 2);
}

TEST(SparsityPatternTest, SparsityPattern3DVectorValuedMorton) {
    test_mars_distributed_sparsity_pattern<ElementType::Hex8, 1, true, MortonKey<Unsigned>>(4, 4, 4, 2);
}

TEST(SparsityPatternTest, SparsityPattern2DLargeMorton) {
    test_mars_distributed_sparsity_pattern<ElementType::Quad4, 1, true, MortonKey<Unsigned>>(64, 64, 0, 1);
}

TEST(SparsityPatternTest, SparsityPattern3DLargeMorton) {
    test_mars_distributed_sparsity_pattern<ElementType::Hex8, 1, true, MortonKey<Unsigned>>(32, 32, 32, 1);
}

TEST(SparsityPatternTest, SparsityPattern2DHilbert) {
    test_mars_distributed_sparsity_pattern<ElementType::Quad4, 1, true, HilbertKey<Unsigned>>(4, 4, 0, 1);
}

TEST(SparsityPatternTest, SparsityPattern3DHilbert) {
    test_mars_distributed_sparsity_pattern<ElementType::Hex8, 1, true, HilbertKey<Unsigned>>(4, 4, 4, 1);
}

TEST(SparsityPatternTest, SparsityPattern2DVectorValuedHilbert) {
    test_mars_distributed_sparsity_pattern<ElementType::Quad4, 1, true, HilbertKey<Unsigned>>(4, 4, 0, 2);
}

TEST(SparsityPatternTest, SparsityPattern3DVectorValuedHilbert) {
    test_mars_distributed_sparsity_pattern<ElementType::Hex8, 1, true, HilbertKey<Unsigned>>(4, 4, 4, 2);
}

TEST(SparsityPatternTest, SparsityPattern2DLargeHilbert) {
    test_mars_distributed_sparsity_pattern<ElementType::Quad4, 1, true, HilbertKey<Unsigned>>(64, 64, 0, 1);
}

TEST(SparsityPatternTest, SparsityPattern3DLargeHilbert) {
    test_mars_distributed_sparsity_pattern<ElementType::Hex8, 1, true, HilbertKey<Unsigned>>(32, 32, 32, 1);
}

int main(int argc, char **argv) {
    mars::Env env(argc, argv);
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
