#include "mars_structured_test.hpp"

using namespace mars;

template <class SfcKey>
void check_invalid_mesh_inputs()
{
    auto context = make_context(proc_allocation{}, MPI_COMM_WORLD);
    DistributedMesh<SfcKey, ElementType::Quad4> mesh(context);
    EXPECT_THROW(generate_distributed_cube(mesh, 0, 3, 0), std::invalid_argument);
    EXPECT_THROW(generate_distributed_cube(mesh, -1, 3, 0), std::invalid_argument);
    EXPECT_THROW(generate_distributed_cube(mesh, 3, 3, 1), std::invalid_argument);
    EXPECT_THROW(generate_distributed_cube(mesh, std::numeric_limits<Integer>::max(), 3, 0), std::invalid_argument);
    DofHandler<decltype(mesh), 1> unbuilt(mesh);
    EXPECT_THROW(unbuilt.enumerate_dofs(), std::invalid_argument);
    if (num_ranks(context) > 1)
    {
        EXPECT_THROW(generate_distributed_cube(mesh, 1, 1, 0), std::invalid_argument);
        const Integer local_width = rank(context) == 0 ? 4 : 0;
        EXPECT_THROW(generate_distributed_cube(mesh, local_width, 4, 0), std::invalid_argument);
        const Integer different_width = rank(context) == 0 ? 4 : 5;
        EXPECT_THROW(generate_distributed_cube(mesh, different_width, 4, 0), std::invalid_argument);
    }

    generate_distributed_cube(mesh, 4, 4, 0);
    mesh.set_XDim(static_cast<Integer>(structured_coordinate_limit<ElementType::Quad4, SfcKey>() / 2));
    DofHandler<decltype(mesh), 2> too_fine(mesh);
    EXPECT_THROW(too_fine.enumerate_dofs(), std::invalid_argument);
}

TEST(StructuredMesh, MortonQuad4Rectangular)
{
    test::check_mesh<ElementType::Quad4, MortonKey<Unsigned>>(5, 3, 0);
    check_invalid_mesh_inputs<MortonKey<Unsigned>>();
    check_invalid_mesh_inputs<HilbertKey<Unsigned>>();
}

TEST(StructuredMesh, MortonHex8Rectangular) { test::check_mesh<ElementType::Hex8, MortonKey<Unsigned>>(5, 3, 2); }

TEST(StructuredMesh, MortonQuad4FullCurve) { test::check_mesh<ElementType::Quad4, MortonKey<Unsigned>>(4, 4, 0); }

TEST(StructuredMesh, MortonHex8FullCurve) { test::check_mesh<ElementType::Hex8, MortonKey<Unsigned>>(4, 4, 4); }

TEST(StructuredMesh, HilbertQuad4Rectangular) { test::check_mesh<ElementType::Quad4, HilbertKey<Unsigned>>(5, 3, 0); }

TEST(StructuredMesh, HilbertHex8Rectangular) { test::check_mesh<ElementType::Hex8, HilbertKey<Unsigned>>(5, 3, 2); }

TEST(StructuredMesh, HilbertQuad4FullCurve) { test::check_mesh<ElementType::Quad4, HilbertKey<Unsigned>>(4, 4, 0); }

TEST(StructuredMesh, HilbertHex8FullCurve) { test::check_mesh<ElementType::Hex8, HilbertKey<Unsigned>>(4, 4, 4); }

TEST(StructuredMesh, MortonQuad4Thin) { test::check_mesh<ElementType::Quad4, MortonKey<Unsigned>>(257, 1, 0); }
TEST(StructuredMesh, HilbertQuad4Thin) { test::check_mesh<ElementType::Quad4, HilbertKey<Unsigned>>(1, 257, 0); }
TEST(StructuredMesh, MortonHex8Thin) { test::check_mesh<ElementType::Hex8, MortonKey<Unsigned>>(129, 1, 2); }
TEST(StructuredMesh, HilbertHex8Thin) { test::check_mesh<ElementType::Hex8, HilbertKey<Unsigned>>(1, 129, 2); }

int main(int argc, char** argv)
{
    return mars::test::run_structured_tests(argc, argv, 12);
}
