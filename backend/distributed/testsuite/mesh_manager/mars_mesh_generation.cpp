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

TEST(StructuredMesh, HilbertQuad4FullCurve) {
    for (Integer width : {2, 4, 8})
        test::check_mesh<ElementType::Quad4, HilbertKey<Unsigned>>(width, width, 0, false, true);
}

TEST(StructuredMesh, HilbertHex8FullCurve) {
    for (Integer width : {2, 4, 8})
        test::check_mesh<ElementType::Hex8, HilbertKey<Unsigned>>(width, width, width, false, true);
}

template<Integer Type>
void check_full_hilbert_subcommunicators() {
    MPI_Comm comm = MPI_COMM_NULL;
    // The four-rank gate creates a three-rank group and an independent singleton.
    ASSERT_EQ(MPI_Comm_split(MPI_COMM_WORLD, test::rank() < 3 ? 0 : 1, test::rank(), &comm), MPI_SUCCESS);
    {
        auto context = make_context(proc_allocation{}, comm);
        DistributedMesh<HilbertKey<Unsigned>, Type> mesh(context);
        generate_distributed_cube(mesh, 4, 4, Type == ElementType::Quad4 ? 0 : 4);
        const Integer total = Type == ElementType::Quad4 ? 16 : 64;
        const int ranks = num_ranks(context), rank_id = rank(context);
        const auto partition = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, mesh.get_view_gp());
        const auto keys = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, mesh.get_view_sfc());
        Integer offset = 0;
        for (int peer = 0; peer <= ranks; ++peer) {
            EXPECT_EQ(partition(2 * peer), offset);
            EXPECT_EQ(partition(2 * peer + 1), offset);
            if (peer == rank_id)
                for (std::size_t i = 0; i < keys.extent(0); ++i) EXPECT_EQ(keys(i), offset + i);
            offset += total / ranks + (peer < total % ranks);
        }
        EXPECT_EQ(keys.extent(0), total / ranks + (rank_id < total % ranks));
        EXPECT_EQ(mesh.partition_received_keys(), 0);
        std::vector<int> counts(total, 0), global(total, 0);
        for (std::size_t i = 0; i < keys.extent(0); ++i) {
            EXPECT_LT(keys(i), total);
            if (keys(i) < total) ++counts[keys(i)];
        }
        EXPECT_EQ(MPI_Allreduce(counts.data(), global.data(), total, MPI_INT, MPI_SUM, comm), MPI_SUCCESS);
        for (const int count : global) EXPECT_EQ(count, 1);
    }
    EXPECT_EQ(MPI_Comm_free(&comm), MPI_SUCCESS);
}

TEST(StructuredMesh, HilbertQuad4UnevenGroups) { check_full_hilbert_subcommunicators<ElementType::Quad4>(); }
TEST(StructuredMesh, HilbertHex8UnevenGroups) { check_full_hilbert_subcommunicators<ElementType::Hex8>(); }

TEST(StructuredMesh, MortonQuad4Thin) { test::check_mesh<ElementType::Quad4, MortonKey<Unsigned>>(257, 1, 0); }
TEST(StructuredMesh, HilbertQuad4Thin) { test::check_mesh<ElementType::Quad4, HilbertKey<Unsigned>>(1, 257, 0); }
TEST(StructuredMesh, MortonHex8Thin) { test::check_mesh<ElementType::Hex8, MortonKey<Unsigned>>(129, 1, 2); }
TEST(StructuredMesh, HilbertHex8Thin) { test::check_mesh<ElementType::Hex8, HilbertKey<Unsigned>>(1, 129, 2); }

int main(int argc, char** argv)
{
    return mars::test::run_structured_tests(argc, argv, 14);
}
