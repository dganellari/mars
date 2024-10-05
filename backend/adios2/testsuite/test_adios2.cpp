// For some reason this does not compile with c++14

// #include "gtest/gtest.h"

// #include "mars.hpp"
// using namespace mars;

// template <typename MeshType>
// class Adios2IOTest : public testing::Test {
// public:
//     void SetUp() override {
//         constexpr Integer n_x_dim = 10;
//         for (Integer d = 0; d < MeshType::Dim; ++d) {
//             dims[d] = n_x_dim;
//         }
//     }

// protected:
//     std::array<Integer, 4> dims{0, 0, 0, 0};
// };

// using Adios2IOTestTypes = ::testing::Types<Mesh2, Mesh3>;

// TYPED_TEST_SUITE(Adios2IOTest, Adios2IOTestTypes);

// TYPED_TEST(Adios2IOTest, Generate) {
//     TypeParam mesh;
//     generate_cube(mesh, this->dims[0], this->dims[1], this->dims[2]);
//     mesh.build_dual_graph();
//     EXPECT_TRUE(mesh.check_side_ordering());
// }

#include "gtest/gtest.h"
#include "mars.hpp"
#include "mars_adios2_IO.hpp"

using namespace mars;

using DMesh2 = ::mars::DistributedMesh<::mars::ElementType::Quad4>;
using DMesh3 = ::mars::DistributedMesh<::mars::ElementType::Hex8>;

template <typename MeshType>
class Adios2IOTest {
public:
    static constexpr Integer Block = 0;
    static constexpr Integer Degree = 1;
    static constexpr int Dim = MeshType::Dim;
    // static constexpr bool Overlap = true;
    using DOFHandler_t = mars::DofHandler<MeshType, Degree, Block>;
    using FEDofMap_t = mars::FEDofMap<DOFHandler_t>;
    using Adios2IO_t = mars::adios2::IO<FEDofMap_t>;

    void set_up() {
        constexpr Integer n_x_dim = 200;
        for (int d = 0; d < Dim; ++d) {
            dims[d] = n_x_dim;
        }
    }

    void run_IO_test() {
        set_up();

        using namespace mars;
        mars::proc_allocation resources;

#ifdef MARS_ENABLE_MPI
        // create a distributed context
        auto context = mars::make_context(resources, MPI_COMM_WORLD);
        int proc_num = mars::rank(context);
#else
        // resources.gpu_id = marsenv::default_gpu();
        // // create a local context
        // auto context = mars::make_context(resources);
#endif

#ifdef MARS_ENABLE_KOKKOS

        Kokkos::Timer timer;

        MeshType mesh(context);
        generate_distributed_cube(mesh, dims[0], dims[1], dims[2]);

        double time_gen = timer.seconds();
        std::cout << "Mesh Generation took: " << time_gen << std::endl;
        /*
                Kokkos::Timer timer_dof;

                DOFHandler_t dof_handler(&mesh);
                dof_handler.enumerate_dofs();
                dof_handler.set_block(1);

                double time_dh = timer_dof.seconds();
                std::cout << "DOFHandler enum took: " << time_dh << std::endl;
                auto fe = build_fe_dof_map<DOFHandler_t>(dof_handler);
         */
        Kokkos::Timer timer_map;

        // builds FEDOFMAP from the mesh and block 1.
        auto fe = build_fe_dof_map<DOFHandler_t>(mesh, 1);
        auto dof_handler = fe.get_dof_handler();

        ViewVectorType<Real> data("data", dof_handler.get_owned_dof_size());

        Kokkos::deep_copy(data, proc_num);

        Adios2IO_t io(fe);
        io.write("out_" + std::to_string(Dim) + ".bp", data);
#endif
    }

protected:
    std::array<Integer, 4> dims{0, 0, 0, 0};
};

TEST(Adios2IOTest2, Generation) { Adios2IOTest<DMesh2>().run_IO_test(); }
TEST(Adios2IOTest3, Generation) { Adios2IOTest<DMesh3>().run_IO_test(); }
