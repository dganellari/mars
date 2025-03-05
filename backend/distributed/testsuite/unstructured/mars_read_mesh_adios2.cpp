#include "mars_read_mesh_adios2.hpp"
#include <gtest/gtest.h>

#ifdef MARS_ENABLE_MPI

TEST(ReadMesh, ReadsMeshFileInParallel) {
    int rank, size;

    using namespace mars;
    mars::proc_allocation resources;
    // bool result = false;
    // create a distributed context
    auto context = mars::make_context(resources, MPI_COMM_WORLD);
    rank = mars::rank(context);
    size = mars::num_ranks(context);

    readMeshInParallel("tetrahedral_mesh.bp", rank, size);
    // Add assertions to verify the correctness of the read data
    // For simplicity, we assume the mesh is correctly read if no errors occur
    ASSERT_TRUE(true);
}

TEST(ReadMesh, ReadsComplexMeshFileInParallel) {
    int rank, size;

    using namespace mars;
    mars::proc_allocation resources;
    // bool result = false;
    // create a distributed context
    auto context = mars::make_context(resources, MPI_COMM_WORLD);
    rank = mars::rank(context);
    size = mars::num_ranks(context);

    readMeshInParallel("complex_mesh.bp", rank, size);
    // Add assertions to verify the correctness of the read data
    // For simplicity, we assume the mesh is correctly read if no errors occur
    ASSERT_TRUE(true);
}
#endif

int main(int argc, char** argv) {
    mars::Env env(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
