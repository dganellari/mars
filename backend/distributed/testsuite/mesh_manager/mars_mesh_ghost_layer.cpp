#include "mars_structured_test.hpp"
#include "mars_structured_design_check.hpp"

using namespace mars;

TEST(StructuredHalo, MortonQuad4Q1) { test::check_dofs<ElementType::Quad4, 1, MortonKey<Unsigned>>(5, 3, 0, true); }

TEST(StructuredHalo, MortonQuad4Q2) { test::check_dofs<ElementType::Quad4, 2, MortonKey<Unsigned>>(5, 3, 0, true); }

TEST(StructuredHalo, MortonHex8Q1) { test::check_dofs<ElementType::Hex8, 1, MortonKey<Unsigned>>(5, 3, 2, true); }

TEST(StructuredHalo, MortonHex8Q2) { test::check_dofs<ElementType::Hex8, 2, MortonKey<Unsigned>>(5, 3, 2, true); }

TEST(StructuredHalo, MortonQuad4PeriodicCells)
{
    test::check_mesh<ElementType::Quad4, MortonKey<Unsigned>>(5, 3, 0, true);
}

TEST(StructuredHalo, MortonHex8PeriodicCells)
{
    test::check_mesh<ElementType::Hex8, MortonKey<Unsigned>>(5, 3, 2, true);
}

TEST(StructuredHalo, HilbertQuad4Q1) { test::check_dofs<ElementType::Quad4, 1, HilbertKey<Unsigned>>(5, 3, 0, true); }

TEST(StructuredHalo, HilbertQuad4Q2) { test::check_dofs<ElementType::Quad4, 2, HilbertKey<Unsigned>>(5, 3, 0, true); }

TEST(StructuredHalo, HilbertHex8Q1) { test::check_dofs<ElementType::Hex8, 1, HilbertKey<Unsigned>>(5, 3, 2, true); }

TEST(StructuredHalo, HilbertHex8Q2) { test::check_dofs<ElementType::Hex8, 2, HilbertKey<Unsigned>>(5, 3, 2, true); }

TEST(StructuredHalo, MortonQuad4Q1Vector)
{
    test::check_dofs<ElementType::Quad4, 1, MortonKey<Unsigned>, 2>(5, 3, 0, true);
    test::check_dofs<ElementType::Quad4, 1, MortonKey<Unsigned>, 0>(5, 3, 0, true, 2);
}

TEST(StructuredHalo, MortonQuad4Q2Vector)
{
    test::check_dofs<ElementType::Quad4, 2, MortonKey<Unsigned>, 2>(5, 3, 0, true);
}

TEST(StructuredHalo, MortonHex8Q1Vector)
{
    test::check_dofs<ElementType::Hex8, 1, MortonKey<Unsigned>, 2>(5, 3, 2, true);
}

TEST(StructuredHalo, MortonHex8Q2Vector)
{
    test::check_dofs<ElementType::Hex8, 2, MortonKey<Unsigned>, 2>(5, 3, 2, true);
}

TEST(StructuredHalo, HilbertQuad4Q1Vector)
{
    test::check_dofs<ElementType::Quad4, 1, HilbertKey<Unsigned>, 2>(5, 3, 0, true);
}

TEST(StructuredHalo, HilbertQuad4Q2Vector)
{
    test::check_dofs<ElementType::Quad4, 2, HilbertKey<Unsigned>, 2>(5, 3, 0, true);
}

TEST(StructuredHalo, HilbertHex8Q1Vector)
{
    test::check_dofs<ElementType::Hex8, 1, HilbertKey<Unsigned>, 2>(5, 3, 2, true);
}

TEST(StructuredHalo, HilbertHex8Q2Vector)
{
    test::check_dofs<ElementType::Hex8, 2, HilbertKey<Unsigned>, 2>(5, 3, 2, true);
}

TEST(StructuredHalo, HilbertQuad4PeriodicCells)
{
    test::check_mesh<ElementType::Quad4, HilbertKey<Unsigned>>(5, 3, 0, true);
}

TEST(StructuredHalo, HilbertHex8PeriodicCells)
{
    test::check_mesh<ElementType::Hex8, HilbertKey<Unsigned>>(5, 3, 2, true);
}

TEST(StructuredHalo, DeviceSnapshotsAndIsolatedPlans) { test::check_device_snapshots_and_plans(); }
TEST(StructuredHalo, SparseSubsetPlan) { test::check_sparse_subset_plan(); }
TEST(StructuredHalo, SharedManagerCopies) { test::check_manager_copies(); }

int main(int argc, char** argv)
{
    return mars::test::run_structured_tests(argc, argv, 23);
}
