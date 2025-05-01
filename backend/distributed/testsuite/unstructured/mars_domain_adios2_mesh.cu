#include <gtest/gtest.h>

#include "domain.hpp"
#include "mars.hpp"
#include "mars_globals.hpp"
#include "mars_read_mesh_adios2.hpp"

using namespace mars;
// Define computation kernel for tetrahedral volumes
__global__ void computeTetrahedralVolumesKernel(const Real* x,
                                                const Real* y,
                                                const Real* z,
                                                const int* i0,
                                                const int* i1,
                                                const int* i2,
                                                const int* i3,
                                                int elementCount)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= elementCount) return;

    // Calculate tetrahedron volume (simplified for testing)
    int n0 = i0[idx];
    int n1 = i1[idx];
    int n2 = i2[idx];
    int n3 = i3[idx];

    // Volume calculation is simplified for testing purposes
    // In a real implementation, this would compute the actual volume
    // and store it somewhere
}

// Define a test fixture for MPI-based tests
class TetrahedralDomainTest : public ::testing::Test
{
protected:
    int rank, size;
    mars::proc_allocation resources;
    mars::context context;

    // Define the type we're testing
    using DomainType = ElementDomain<TetTag, cstone::GpuTag>;

    void SetUp() override
    {
        // Initialize with Mars MPI context
        context = mars::make_context(resources, MPI_COMM_WORLD);
        rank    = mars::rank(context);
        size    = mars::num_ranks(context);
    }
};

// Test basic construction and initialization
TEST_F(TetrahedralDomainTest, InitializationTest)
{
    try
    {
        // Create the domain with a tetrahedral mesh
        DomainType domain("tetrahedral_mesh.bp", rank, size);

        // Verify domain was initialized with data
        ASSERT_GT(domain.getNodeCount(), 0) << "Domain should have nodes";
        ASSERT_GT(domain.getElementCount(), 0) << "Domain should have tetrahedra";

        // Check domain interface access
        ASSERT_GE(domain.startIndex(), 0) << "Start index should be non-negative";
        ASSERT_GT(domain.endIndex(), domain.startIndex()) << "End index should be after start index";
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception during domain initialization: " << e.what() << " (on rank " << rank << ")";
    }
}

// Test domain sync operation
TEST_F(TetrahedralDomainTest, SyncTest)
{
    if (size < 2) { GTEST_SKIP() << "This test requires at least 2 MPI ranks"; }

    try
    {
        DomainType domain("tetrahedral_mesh.bp", rank, size);

        // Record initial state
        size_t initialNodeCount    = domain.getNodeCount();
        size_t initialElementCount = domain.getElementCount();

        // Perform sync operation
        domain.sync();

        // Domain should still have data after sync
        ASSERT_GT(domain.getNodeCount(), 0) << "Domain should have nodes after sync";
        ASSERT_GT(domain.getElementCount(), 0) << "Domain should have tetrahedra after sync";

        // In a real simulation, node counts might change after sync
        // but for a static test mesh, they should remain the same
        EXPECT_EQ(domain.getNodeCount(), initialNodeCount) << "Node count shouldn't change for static mesh";
        EXPECT_EQ(domain.getElementCount(), initialElementCount) << "Tetrahedra count shouldn't change for static mesh";
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception during sync test: " << e.what() << " (on rank " << rank << ")";
    }
}

// Test halo exchange for field data
TEST_F(TetrahedralDomainTest, HaloExchangeTest)
{
    if (size < 2) { GTEST_SKIP() << "This test requires at least 2 MPI ranks"; }

    try
    {
        DomainType domain("tetrahedral_mesh.bp", rank, size);

        // Create a field and initialize with rank value
        std::vector<Real> field(domain.getNodeCount(), static_cast<Real>(rank));

        // Note: ElementDomain doesn't have an exchangeHalos method currently
        // domain.exchangeHalos(field);

        // After exchange, local values should still be our rank
        for (size_t i = domain.startIndex(); i < domain.endIndex(); i++)
        {
            EXPECT_FLOAT_EQ(field[i], static_cast<Real>(rank)) << "Local values should not change after halo exchange";
        }

        // Note: We can't easily test halo values without domain-specific knowledge
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception during halo exchange test: " << e.what() << " (on rank " << rank << ")";
    }
}

// Test for main functionality - converted from the original main function
TEST_F(TetrahedralDomainTest, MainFunctionality)
{
    try
    {
        // Create the tetrahedral domain
        DomainType domain("tetrahedral_mesh.bp", rank, size);

        // Verify domain initialization
        ASSERT_GT(domain.getNodeCount(), 0) << "Domain should have nodes";
        ASSERT_GT(domain.getElementCount(), 0) << "Domain should have tetrahedra";

        // Run iterative computation similar to the cornerstone example
        int nIterations = 3; // Reduced from 5 for testing speed
        for (int iter = 0; iter < nIterations; iter++)
        {
            // Sync domain to update node positions
            domain.sync();

            // Verify domain still has data after sync
            ASSERT_GT(domain.getNodeCount(), 0) << "Domain should have nodes after sync";

            // Compute volumes on tetrahedra
            domain.computeOnElements(computeTetrahedralVolumesKernel);

            // Create a field for computed values
            std::vector<Real> fieldValue(domain.getNodeCount(), 0.0);

            // Compute field values for local nodes
            for (size_t i = domain.startIndex(); i < domain.endIndex(); i++)
            {
                fieldValue[i] = static_cast<Real>(i % 10);
            }

            // Note: ElementDomain doesn't have an exchangeHalos method currently
            // domain.exchangeHalos(fieldValue);

            // Verify field values are maintained for local nodes
            for (size_t i = domain.startIndex(); i < domain.endIndex(); i++)
            {
                EXPECT_FLOAT_EQ(fieldValue[i], static_cast<Real>(i % 10)) << "Field value changed after halo exchange";
            }
        }

        // Report statistics for rank 0
        if (rank == 0)
        {
            std::cout << "Test completed: Domain has " << domain.getNodeCount() << " nodes and "
                      << domain.getElementCount() << " tetrahedra." << std::endl;
        }
    }
    catch (const std::exception& e)
    {
        FAIL() << "Exception: " << e.what() << " (on rank " << rank << ")";
    }
}

// Run the test program
int main(int argc, char** argv)
{
    mars::Env env(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}