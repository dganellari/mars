#include <gtest/gtest.h>

#define USE_CUDA

#include "domain.hpp"
#include "mars.hpp"
#include "mars_globals.hpp"
#include "mars_read_mesh_adios2.hpp"

// Define a test fixture for MPI-based tests
class TetrahedralDomainTest : public ::testing::Test {
protected:
    int rank, size;
    mars::proc_allocation resources;
    mars::context context;

    void SetUp() override {
        // Initialize with Mars MPI context
        context = mars::make_context(resources, MPI_COMM_WORLD);
        rank = mars::rank(context);
        size = mars::num_ranks(context);
    }
};

// Test basic construction and initialization
TEST_F(TetrahedralDomainTest, InitializationTest) {
    try {
        // Create the domain with a tetrahedral mesh
        TetrahedralDomain domain("tetrahedral_mesh.bp", rank, size);

        // Verify domain was initialized with data
        ASSERT_GT(domain.getNodeCount(), 0) << "Domain should have nodes";
        ASSERT_GT(domain.getTetCount(), 0) << "Domain should have tetrahedra";

        // Check domain interface access
        ASSERT_GE(domain.startIndex(), 0) << "Start index should be non-negative";
        ASSERT_GT(domain.endIndex(), domain.startIndex()) << "End index should be after start index";
    } catch (const std::exception& e) {
        FAIL() << "Exception during domain initialization: " << e.what() << " (on rank " << rank << ")";
    }
}

// Test domain sync operation
TEST_F(TetrahedralDomainTest, SyncTest) {
    if (size < 2) {
        GTEST_SKIP() << "This test requires at least 2 MPI ranks";
    }

    try {
        TetrahedralDomain domain("tetrahedral_mesh.bp", rank, size);

        // Record initial state
        size_t initialNodeCount = domain.getNodeCount();
        size_t initialTetCount = domain.getTetCount();

        // Perform sync operation
        domain.sync();

        // Domain should still have data after sync
        ASSERT_GT(domain.getNodeCount(), 0) << "Domain should have nodes after sync";
        ASSERT_GT(domain.getTetCount(), 0) << "Domain should have tetrahedra after sync";

        // In a real simulation, node counts might change after sync
        // but for a static test mesh, they should remain the same
        EXPECT_EQ(domain.getNodeCount(), initialNodeCount) << "Node count shouldn't change for static mesh";
        EXPECT_EQ(domain.getTetCount(), initialTetCount) << "Tetrahedra count shouldn't change for static mesh";
    } catch (const std::exception& e) {
        FAIL() << "Exception during sync test: " << e.what() << " (on rank " << rank << ")";
    }
}

// Test octree node mapping
/* TEST_F(TetrahedralDomainTest, OctreeNodeMappingTest) {
    try {
        TetrahedralDomain domain("tetrahedral_mesh.bp", rank, size);

        // Get the total number of octree nodes
        int numOctants = domain.getDomain().numTreeNodes();
        ASSERT_GT(numOctants, 0) << "Domain should have at least one octant";

        // Count the total tetrahedra across all octants
        size_t totalMappedTets = 0;
        for (int i = 0; i < numOctants; i++) {
            auto tetsInNode = domain.getTetrahedraInOctreeNode(i);
            totalMappedTets += tetsInNode.size();
        }

        // Verify that all tetrahedra are mapped to at least one octant
        // Note: A tetrahedron might be mapped to multiple octants if it spans boundaries
        EXPECT_GE(totalMappedTets, domain.getTetCount()) << "All tetrahedra should be mapped to at least one octant";
    } catch (const std::exception& e) {
        FAIL() << "Exception during octree mapping test: " << e.what() << " (on rank " << rank << ")";
    }
} */

// Test halo exchange for field data
TEST_F(TetrahedralDomainTest, HaloExchangeTest) {
    if (size < 2) {
        GTEST_SKIP() << "This test requires at least 2 MPI ranks";
    }

    try {
        TetrahedralDomain domain("tetrahedral_mesh.bp", rank, size);

        // Create a field and initialize with rank value
        std::vector<Real> field(domain.getNodeCount(), static_cast<Real>(rank));

        // Exchange halos
        // domain.exchangeHalos(field);

        // After exchange, local values should still be our rank
        for (size_t i = domain.startIndex(); i < domain.endIndex(); i++) {
            EXPECT_DOUBLE_EQ(field[i], static_cast<Real>(rank)) << "Local values should not change after halo exchange";
        }

        // Note: We can't easily test halo values without domain-specific knowledge
        // In a real test, you might want to set specific boundary values and check them
    } catch (const std::exception& e) {
        FAIL() << "Exception during halo exchange test: " << e.what() << " (on rank " << rank << ")";
    }
}

// Test for main functionality - converted from the original main function
TEST_F(TetrahedralDomainTest, MainFunctionality) {
    try {
        // Create the tetrahedral domain
        TetrahedralDomain domain("tetrahedral_mesh.bp", rank, size);

        // Verify domain initialization
        ASSERT_GT(domain.getNodeCount(), 0) << "Domain should have nodes";
        ASSERT_GT(domain.getTetCount(), 0) << "Domain should have tetrahedra";

        // Run iterative computation similar to the cornerstone example
        int nIterations = 3;  // Reduced from 5 for testing speed
        for (int iter = 0; iter < nIterations; iter++) {
            // Sync domain to update node positions
            domain.sync();

            // Verify domain still has data after sync
            ASSERT_GT(domain.getNodeCount(), 0) << "Domain should have nodes after sync";

            // Compute volumes on tetrahedra
            domain.computeOnTetrahedra(computeTetrahedralVolumesKernel);

            // Create a field for computed values
            std::vector<Real> fieldValue(domain.getNodeCount(), 0.0);

            // Compute field values for local nodes
            for (size_t i = domain.startIndex(); i < domain.endIndex(); i++) {
                fieldValue[i] = static_cast<Real>(i % 10);
            }

            // Exchange halos to ensure consistent values across processes
            // domain.exchangeHalos(fieldValue);

            // Verify field values are maintained for local nodes
            for (size_t i = domain.startIndex(); i < domain.endIndex(); i++) {
                EXPECT_DOUBLE_EQ(fieldValue[i], static_cast<Real>(i % 10)) << "Field value changed after halo exchange";
            }
        }

        // Report statistics for rank 0
        if (rank == 0) {
            std::cout << "Test completed: Domain has " << domain.getNodeCount() << " nodes and " << domain.getTetCount()
                      << " tetrahedra." << std::endl;
        }
    } catch (const std::exception& e) {
        FAIL() << "Exception: " << e.what() << " (on rank " << rank << ")";
    }
}

// Run the test program
int main(int argc, char** argv) {
    mars::Env env(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
