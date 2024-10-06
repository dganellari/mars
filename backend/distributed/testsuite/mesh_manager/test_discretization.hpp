#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include "template_mesh_manager.hpp"

TEST(MarsTest, TestDOFHandlerEnumeration) {
    // Define the dimensions of the mesh
    const int xDim = 10;
    const int yDim = 10;
    const int zDim = 10;
    const int block = 1;

    // Create a context
    mars::proc_allocation resources;
    auto context = mars::make_context(resources, MPI_COMM_WORLD);

    // Create a mesh
    using DMesh = DistributedMesh<MortonKey<Unsigned>, ElementType::Quad4>;
    DMesh mesh(context);
    generate_distributed_cube(mesh, xDim, yDim, zDim);

    // Create a DOFHandler and enumerate DOFs
    using DOFHandler = DofHandler<DMesh, 1, 0>;
    DOFHandler dof_handler(mesh);
    dof_handler.enumerate_dofs();
    dof_handler.set_block(block);

    // Check that the DOFHandler has the correct number of DOFs
    Unsigned expected_num_dofs = xDim * yDim * zDim;  // Replace with your expected value
    ASSERT_EQ(dof_handler.get_dof_size(), expected_num_dofs);

    // Check that the block size is set correctly
    ASSERT_EQ(dof_handler.get_block_size(), block);

    // Check that each DOF has a unique identifier
    // Create a Kokkos View to store the DOF identifiers
    Kokkos::View<Unsigned*> dof_ids("dof_ids", dof_handler.get_dof_size());

    // Use a Kokkos parallel_for to fill the View with DOF identifiers
    Kokkos::parallel_for("fill_dof_ids", dof_handler.get_dof_size(), KOKKOS_LAMBDA(const int i) {
        dof_ids(i) = dof_handler.get_dof_id(i);
    });

    // Use a Kokkos parallel_reduce to check for unique identifiers
    bool unique_ids = true;
    Kokkos::parallel_reduce("check_unique_ids", dof_handler.get_dof_size(), KOKKOS_LAMBDA(const int i, bool& update) {
        for (int j = i + 1; j < dof_handler.get_dof_size(); ++j) {
            if (dof_ids(i) == dof_ids(j)) {
                update = false;
            }
        }
    }, unique_ids);

    ASSERT_TRUE(unique_ids);
}

TEST(DofMapTest, TestDOFMAPGeneration) {
    // Define the dimensions of the mesh
    const int xDim = 10;
    const int yDim = 10;
    const int zDim = 10;
    const int block = 1;

    // Create a context
    mars::proc_allocation resources;
    auto context = mars::make_context(resources, MPI_COMM_WORLD);

    // Create a mesh
    using DMesh = DistributedMesh<MortonKey<Unsigned>, ElementType::Quad4>;
    DMesh mesh(context);
    generate_distributed_cube(mesh, xDim, yDim, zDim);

    // Create a DOFHandler and enumerate DOFs
    using DOFHandler = DofHandler<DMesh, 1, 0>;
    DOFHandler dof_handler(mesh);
    dof_handler.enumerate_dofs();
    dof_handler.set_block(block);

    // Check that the DOFHandler has the correct number of DOFs
    Unsigned expected_num_dofs = xDim * yDim * zDim;  // Replace with your expected value
    ASSERT_EQ(dof_handler.get_dof_size(), expected_num_dofs);

    //overlap true
    auto fe = build_fe_dof_map<DOFHandler, 1>(dof_handler);

    // Create an array to record errors
    Kokkos::View<int*> errors("errors", mesh.elements().extent(0));

    // Check the results on the device
    Kokkos::parallel_for("TestDOFMap", mesh.elements().extent(0), KOKKOS_LAMBDA(const int i) {
        auto element = mesh.elements(i);
        auto dofs = dofMap.get_DOFs(element);

        if (dofs.size() != element.num_DOFs()) {
            errors(i) = 1;
        }

        for (const auto& dof : dofs) {
            if (dof.element() != element) {
                errors(i) = 1;
            }
        }
    });

    // Copy the errors back to the host
    auto h_errors = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), errors);

    // Check if any errors were recorded
    for (int i = 0; i < h_errors.extent(0); ++i) {
        if (h_errors(i) != 0) {
            printf("Error in element %d\n", i);
        }
    }

}
// test the DOFMAP generation

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

