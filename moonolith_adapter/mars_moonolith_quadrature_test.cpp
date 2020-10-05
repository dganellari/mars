#include "mars_moonolith_quadrature_test.hpp"

#include "mars_base.hpp"
#include "mars_mesh.hpp"
#include "mars_bisection.hpp"

#include "moonolith_mesh_adapter.hpp"
#include "moonolith_polygon.hpp"
#include "moonolith_intersect_polygons.hpp"
#include "moonolith_intersect_polyhedra.hpp"
#include "moonolith_intersect_lines.hpp"
#include "mars_moonolith_mesh_adapter.hpp"
#include "moonolith_intersect_polygon_with_h_polyhedron.hpp"
#include "moonolith_intersect_polyline_with_convex_polygon.hpp"
#include "moonolith_duplicate_intersection_avoidance.hpp"
#include "moonolith_build_quadrature.hpp"


namespace mars {

    template<Integer Dim>
    static double mars_compute_quad(
        moonolith::Communicator &comm,
        const moonolith::Quadrature<double, 2> &q_ref,
        const mars::Mesh<Dim, 2> &mesh_master,
        const mars::Mesh<Dim, 2> &mesh_slave)
    {
        using MeshT          = mars::Mesh<Dim, 2>;
        using AlgorithmT     = moonolith::OneMasterOneSlaveAlgorithm<Dim, MeshT>;
        using Adapter = typename AlgorithmT::Adapter;

        double elapsed = MPI_Wtime();

        AlgorithmT algo(comm, moonolith::make_unique<moonolith::CollectionManager<MeshT>>());

        algo.init_simple(
            mesh_master,
            mesh_slave,
            1e-10
        );

        double vol = 0.;

        elapsed = MPI_Wtime() - elapsed;
        moonolith::logger() << "init: " << elapsed  << std::endl;

        ////////////////////////////////////////////////////
        /////////////////// pair-wise method ///////////////////////
        elapsed = MPI_Wtime();


        moonolith::Polygon<double, 2> poly_m, poly_s;
        moonolith::Quadrature2<double> q_out;

        moonolith::BuildQuadrature<moonolith::Polygon<double, 2>> q_builder;

        algo.compute([&](const Adapter &master, const Adapter &slave) -> bool {
            const auto &mm = master.collection();
            const auto &me = master.elem();

            const auto &sm = slave.collection();
            const auto &se = slave.elem();

            make(mm, me, poly_m);
            make(sm, se, poly_s);

            if(q_builder.apply(q_ref, poly_m, poly_s, q_out)) {
                vol += std::accumulate(std::begin(q_out.weights), std::end(q_out.weights), 0.);
                return true;
            }

            return false;
        });

        comm.all_reduce(&vol, 1, moonolith::MPISum());

        elapsed = MPI_Wtime() - elapsed;
        moonolith::logger() << "time(master, slave): " << elapsed  << std::endl;
        moonolith::logger() << "vol " << vol << std::endl;
        return vol;

    }

    static void run_mars_quad_2()
    {
        moonolith::logger() << "--------------------------------\n";
        moonolith::logger() << "run_mars_quad_2" << std::endl;

        moonolith::Communicator comm;

        Integer res = 2;

        mars::Mesh2 mesh_master(true);

        if(comm.rank() == 0 || comm.is_alone()) {
            read_mesh("../data/square_2.MFEM", mesh_master);
            mars::Bisection<mars::Mesh2> b(mesh_master);
            b.uniform_refine(res * 2);
            mesh_master.clean_up();
        }

        mars::Mesh2 mesh_slave(true);

        if(comm.rank() == 1 || comm.is_alone()) {
            read_mesh("../data/square_2.MFEM", mesh_slave);
            mars::Bisection<mars::Mesh2> b(mesh_slave);
            b.uniform_refine(res * 3);
            mesh_slave.clean_up();
        }

        moonolith::logger() << "n_master " << mesh_master.n_elements() << " n_slave " << mesh_slave.n_elements() << " " << std::endl;

        moonolith::Quadrature2<double> q_ref;
        q_ref.points = {{std::sqrt(2.)/2., std::sqrt(2.)/2.}};
        q_ref.weights = {1.};
        auto area = mars_compute_quad(comm, q_ref, mesh_master, mesh_slave);

        if(comm.rank() == 1 || comm.is_alone()) {
            assert(moonolith::approxeq(1., area, 1e-8));
        }
    }


    void run_mars_quadrature_test()
    {
        run_mars_quad_2();
    }

}
