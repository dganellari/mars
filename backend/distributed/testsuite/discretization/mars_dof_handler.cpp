#include "mars_structured_test.hpp"

using namespace mars;

TEST(StructuredDofs, MortonQuad4Q1) { test::check_dofs<ElementType::Quad4, 1, MortonKey<Unsigned>>(5, 3, 0); }

TEST(StructuredDofs, MortonQuad4Q1Vector) { test::check_dofs<ElementType::Quad4, 1, MortonKey<Unsigned>, 2>(5, 3, 0); }

TEST(StructuredDofs, MortonQuad4Q2) { test::check_dofs<ElementType::Quad4, 2, MortonKey<Unsigned>>(5, 3, 0); }

TEST(StructuredDofs, MortonQuad4Q2Vector) { test::check_dofs<ElementType::Quad4, 2, MortonKey<Unsigned>, 2>(5, 3, 0); }

TEST(StructuredDofs, MortonHex8Q1) { test::check_dofs<ElementType::Hex8, 1, MortonKey<Unsigned>>(5, 3, 2); }

TEST(StructuredDofs, MortonHex8Q1Vector) { test::check_dofs<ElementType::Hex8, 1, MortonKey<Unsigned>, 2>(5, 3, 2); }

TEST(StructuredDofs, MortonHex8Q2) { test::check_dofs<ElementType::Hex8, 2, MortonKey<Unsigned>>(5, 3, 2); }

TEST(StructuredDofs, MortonHex8Q2Vector) { test::check_dofs<ElementType::Hex8, 2, MortonKey<Unsigned>, 2>(5, 3, 2); }

TEST(StructuredDofs, HilbertQuad4Q1)
{
    // At four ranks this partition needs a receive-only DOF peer.
    // Count exchange must not assume that DOF communication is symmetric.
    test::check_dofs<ElementType::Quad4, 1, HilbertKey<Unsigned>>(5, 3, 0);
}

TEST(StructuredDofs, HilbertQuad4Q1Vector)
{
    test::check_dofs<ElementType::Quad4, 1, HilbertKey<Unsigned>, 2>(5, 3, 0);
}

TEST(StructuredDofs, HilbertQuad4Q2) { test::check_dofs<ElementType::Quad4, 2, HilbertKey<Unsigned>>(5, 3, 0); }

TEST(StructuredDofs, HilbertQuad4Q2Vector)
{
    test::check_dofs<ElementType::Quad4, 2, HilbertKey<Unsigned>, 2>(5, 3, 0);
}

TEST(StructuredDofs, HilbertHex8Q1) { test::check_dofs<ElementType::Hex8, 1, HilbertKey<Unsigned>>(5, 3, 2); }

TEST(StructuredDofs, HilbertHex8Q1Vector) { test::check_dofs<ElementType::Hex8, 1, HilbertKey<Unsigned>, 2>(5, 3, 2); }

TEST(StructuredDofs, HilbertHex8Q2) { test::check_dofs<ElementType::Hex8, 2, HilbertKey<Unsigned>>(5, 3, 2); }

TEST(StructuredDofs, HilbertHex8Q2Vector) { test::check_dofs<ElementType::Hex8, 2, HilbertKey<Unsigned>, 2>(5, 3, 2); }

int main(int argc, char** argv)
{
    return mars::test::run_structured_tests(argc, argv, 16);
}
