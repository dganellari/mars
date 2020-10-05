#include "mars_moonolith_test.hpp"
#include "mars_base.hpp"
#include "mars_mesh.hpp"
#include "mars_bisection.hpp"
#include "mars_utils.hpp"

#include "moonolith_mesh_adapter.hpp"
#include "moonolith_polygon.hpp"
#include "moonolith_intersect_polygons.hpp"
#include "moonolith_intersect_polyhedra.hpp"
#include "moonolith_intersect_lines.hpp"
#include "moonolith_intersect_polygon_with_h_polyhedron.hpp"
#include "moonolith_intersect_polyline_with_convex_polygon.hpp"
#include "moonolith_duplicate_intersection_avoidance.hpp"

#include "mars_extract_surface.hpp"
#include "mars_moonolith_mesh_adapter.hpp"
#include "mars_moonolith_quadrature_test.hpp"
#include "mars_moonolith_project_test.hpp"
#include "mars_moonolith_l2_transfer_test.hpp"

namespace mars {

    template<Integer Dim, Integer ManifoldDim>
    class IntersectData {};

    template<Integer Dim>
    class IntersectData<Dim, 1> {
    public:
        moonolith::Line<double, Dim> poly_m, poly_s, poly_isect;
        moonolith::IntersectLines<double, Dim> isect;

        template<class Adapter>
        inline double compute(const Adapter &master, const Adapter &slave)
        {
            auto &e_m = master.elem();
            auto &e_s = slave.elem();

            auto &m_m = master.collection();
            auto &m_s = slave.collection();

            make(m_m, e_m, poly_m);
            make(m_s, e_s, poly_s);

            if(isect.apply(poly_m, poly_s, poly_isect)) {
                return moonolith::measure(poly_isect);
            } else {
                return 0.;
            }
        }
    };


    template<Integer Dim>
    class IntersectData<Dim, 2> {
    public:
        moonolith::Polygon<double, Dim> poly_m, poly_s, poly_isect;
        moonolith::IntersectConvexPolygons<double, Dim> isect;

        template<class Adapter>
        inline double compute(const Adapter &master, const Adapter &slave)
        {
            auto &e_m = master.elem();
            auto &e_s = slave.elem();

            auto &m_m = master.collection();
            auto &m_s = slave.collection();

            make(m_m, e_m, poly_m);
            make(m_s, e_s, poly_s);

            if(isect.apply(poly_m, poly_s, poly_isect)) {
                return moonolith::measure(poly_isect);
            } else {
                return 0.;
            }
        }
    };

    template<>
    class IntersectData<3, 3> {
    public:
        moonolith::Polyhedron<double> poly_m, poly_s, poly_isect;
        moonolith::IntersectConvexPolyhedra<double> isect;

        template<class Adapter>
        inline double compute(const Adapter &master, const Adapter &slave)
        {
            auto &e_m = master.elem();
            auto &e_s = slave.elem();

            auto &m_m = master.collection();
            auto &m_s = slave.collection();

            make(m_m, e_m, poly_m);
            make(m_s, e_s, poly_s);

            if(isect.apply(poly_m, poly_s, poly_isect)) {
                return moonolith::measure(poly_isect);
            } else {
                return 0.;
            }
        }
    };

    template<>
    class IntersectData<4, 3> {
    public:
        moonolith::Polyhedron4<double> poly_m, poly_s, poly_isect;
        moonolith::IntersectConvexSurfPolyhedra4<double> isect;

        template<class Adapter>
        inline double compute(const Adapter &master, const Adapter &slave)
        {
            auto &e_m = master.elem();
            auto &e_s = slave.elem();

            auto &m_m = master.collection();
            auto &m_s = slave.collection();

            auto n_m = mars::normal(e_m, m_m.points());
            auto n_s = mars::normal(e_s, m_s.points());

            if(!moonolith::approxeq(mars::dot(n_m, n_s), 1., 1e-10)) {
                return 0.;
            }

            auto v = m_m.point(e_m.nodes[0]) - m_s.point(e_s.nodes[0]);

            //detect elements that are not in the same plane
            if(mars::dot(v, n_s) > 1e-8) {
                return false;
            }

            make(m_m, e_m, poly_m);
            make(m_s, e_s, poly_s);

            if(isect.apply(poly_m, poly_s, poly_isect)) {
                auto ret = isect.isect_measure();
                moonolith::logger() << master.handle() << " -> " << slave.handle() << " " << ret << std::endl;
                return ret;
            } else {
                return 0.;
            }
        }
    };

    template<Integer Dim>
    class Vol2SurfIntersectData {};

    template<class VolElement, class SurfElement>
    class Vol2SurfIntersectDataAux {
    public:

        template<class Adapter>
        inline double compute_aux(const moonolith::Storage<Adapter> &masters, const Adapter &slave)
        {
            auto n = masters.size();
            algo.master_poly.resize(n);

            for(std::size_t i = 0; i < n; ++i) {
                const auto &master = masters[i];
                const auto &m_m = master.collection();
                const auto &e_m = master.elem();
                make(m_m, e_m, algo.master_poly[i]);
            }

            auto &m_s = slave.collection();
            auto &e_s = slave.elem();
            make(m_s, e_s, algo.slave_poly);

            // logger() << "---------------------------------\n" << std::endl;

            if(!algo.apply()) {
                return 0.;
            }

            double ret = 0.;
            for(std::size_t i = 0; i < n; ++i) {
                if(algo.intersected[i]) {
                    const auto len_i = moonolith::measure(algo.isect[i]);
                    ret += len_i;
                    // logger() << masters[i].handle() << "/" << i << ", " << slave.handle() << " -> " << len_i << std::endl;
                }
            }

            auto expected = moonolith::measure(algo.slave_poly);
            assert(moonolith::approxeq(ret, expected, 1e-6));
            return ret;
        }

        moonolith::DuplicateIntersectionAvoidance<VolElement,SurfElement> algo;
    };

    template<>
    class Vol2SurfIntersectData<3> : private Vol2SurfIntersectDataAux<moonolith::Polyhedron<double>, moonolith::Polygon<double, 3>> {
    public:

        template<class Adapter>
        inline double compute(const moonolith::Storage<Adapter> &masters, const Adapter &slave)
        {
            return this->compute_aux(masters, slave);
        }

        template<class Adapter>
        inline double compute(const Adapter &master, const Adapter &slave)
        {
            auto &e_m = master.elem();
            auto &e_s = slave.elem();

            auto &m_m = master.collection();
            auto &m_s = slave.collection();

            make(m_m, e_m, poly_m);
            make(m_s, e_s, poly_s);

            h_poly_m.make(poly_m);

            if(!isect.apply(poly_s, h_poly_m, poly_isect)) {
                return 0.;
            }

            return moonolith::measure(poly_isect);
        }

        moonolith::Polyhedron<double> poly_m;
        moonolith::HPolytope<double, 3> h_poly_m;
        moonolith::Polygon<double, 3> poly_s, poly_isect;
        moonolith::IntersectPolygonWithHPolyhedron<double> isect;
    };



    template<>
    class Vol2SurfIntersectData<2> : private Vol2SurfIntersectDataAux< moonolith::Polygon<double, 2>, moonolith::Line<double, 2> >{
    public:

        template<class Adapter>
        inline double compute(const moonolith::Storage<Adapter> &masters, const Adapter &slave)
        {
            return this->compute_aux(masters, slave);
        }

        template<class Adapter>
        inline double compute(const Adapter &master, const Adapter &slave)
        {
            auto &e_m = master.elem();
            auto &e_s = slave.elem();

            auto &m_m = master.collection();
            auto &m_s = slave.collection();

            make(m_m, e_m, poly_m);
            make(m_s, e_s, poly_s);

            if(!isect.apply(poly_s, poly_m, poly_isect)) {
                return 0.;
            }

            auto ret = moonolith::measure(poly_isect);
            return ret;
        }

        moonolith::Polygon<double, 2> poly_m;
        moonolith::Line<double, 2> poly_s, poly_isect;
        moonolith::IntersectPolylineWithConvexPolygon2<double> isect;
    };

    template<Integer Dim, Integer ManifoldDim>
    static double mars_compute_isect_measure(
        moonolith::Communicator &comm,
        const mars::Mesh<Dim, ManifoldDim> &mesh_master,
        const mars::Mesh<Dim, ManifoldDim> &mesh_slave)
    {
        using MeshT          = mars::Mesh<Dim, ManifoldDim>;
        using AlgorithmT     = moonolith::OneMasterOneSlaveAlgorithm<Dim, MeshT>;
        using IntersectDataT = mars::IntersectData<Dim, ManifoldDim>;
        using Adapter = typename AlgorithmT::Adapter;

        double elapsed = MPI_Wtime();

        AlgorithmT algo(comm, moonolith::make_unique<moonolith::CollectionManager<MeshT>>());

        algo.init_simple(
            mesh_master,
            mesh_slave,
            1e-10
        );

        IntersectDataT isect;

        double vol = 0.;

        elapsed = MPI_Wtime() - elapsed;
        moonolith::logger() << "init: " << elapsed  << std::endl;

        ////////////////////////////////////////////////////
        /////////////////// pair-wise method ///////////////////////
        elapsed = MPI_Wtime();

        algo.compute([&](const Adapter &master, const Adapter &slave) -> bool {
            auto v = isect.compute(master, slave);

            if(v == 0.) return false;

            vol += v;
            return true;
        });

        comm.all_reduce(&vol, 1, moonolith::MPISum());

        elapsed = MPI_Wtime() - elapsed;
        moonolith::logger() << "time(master, slave): " << elapsed  << std::endl;
        moonolith::logger() << "vol " << vol << std::endl;
        return vol;
    }


    template<Integer Dim>
    static double mars_vol2surf_isect_measure(
        moonolith::Communicator &comm,
        const mars::IMesh<Dim> &mesh_master,
        const mars::IMesh<Dim> &mesh_slave)
    {
        using MeshT   = mars::IMesh<Dim>;
        using Adapter = typename moonolith::ManyMastersOneSlaveAlgorithm<Dim, MeshT>::Adapter;

        double elapsed = MPI_Wtime();

        moonolith::ManyMastersOneSlaveAlgorithm<Dim, MeshT> algo(comm, moonolith::make_unique<moonolith::CollectionManager<MeshT>>());

        algo.init_simple(
            mesh_master,
            mesh_slave,
            1e-12
        );

        Vol2SurfIntersectData<Dim> isect;

        double vol = 0.;

        elapsed = MPI_Wtime() - elapsed;
        moonolith::logger() << "init: " << elapsed  << std::endl;

        ////////////////////////////////////////////////////
        /////////////////// many-to-one method ///////////////////////

        elapsed = MPI_Wtime();

        vol = 0.;
        algo.compute([&](const moonolith::Storage<Adapter> &masters, const Adapter &slave) -> bool {
            auto v = isect.compute(masters, slave);

            if(v == 0.) return false;

            vol += v;
            return true;
        });

        comm.all_reduce(&vol, 1, moonolith::MPISum());

        elapsed = MPI_Wtime() - elapsed;
        moonolith::logger() << "time(vector<master>, slave): " << elapsed  << std::endl;
        moonolith::logger() << "vol " << vol << std::endl;

        return vol;
    }

    static void run_mars_intersect_1()
    {
        moonolith::logger() << "--------------------------------\n";
        moonolith::logger() << "run_mars_intersect_1" << std::endl;

        moonolith::Communicator comm;

        Integer res = 2;

        mars::Mesh<1, 1> mesh_master(true);

        if(comm.rank() == 0 || comm.is_alone()) {
            mesh_master.reserve(1, 2);
            mesh_master.add_point({0.});
            mesh_master.add_point({1.});
            mesh_master.add_elem({0, 1});
        }

        mars::Mesh<1, 1> mesh_slave(true);

        if(comm.rank() == 1 || comm.is_alone()) {

            mesh_slave.reserve(1, 3);
            mesh_slave.add_point({0.});
            mesh_slave.add_point({0.5});
            mesh_slave.add_point({1.});
            mesh_slave.add_elem({0, 1});
            mesh_slave.add_elem({1, 2});
        }

        moonolith::logger() << "n_master " << mesh_master.n_elements() << " n_slave " << mesh_slave.n_elements() << " " << std::endl;
        auto area = mars_compute_isect_measure(comm, mesh_master, mesh_slave);

        if(comm.rank() == 1 || comm.is_alone()) {
            assert(moonolith::approxeq(1., area, 1e-8));
        }
    }


    static void run_mars_intersect_2()
    {
        moonolith::logger() << "--------------------------------\n";
        moonolith::logger() << "run_mars_intersect_2" << std::endl;

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
        auto area = mars_compute_isect_measure(comm, mesh_master, mesh_slave);

        if(comm.rank() == 1 || comm.is_alone()) {
            assert(moonolith::approxeq(1., area, 1e-8));
        }
    }

    static void run_mars_intersect_3()
    {
        moonolith::logger() << "--------------------------------\n";
        moonolith::logger() << "run_mars_intersect_3" << std::endl;

        moonolith::Communicator comm;

        Integer res = 1;

        mars::Mesh3 mesh_master(true);

        if(comm.rank() == 0 || comm.is_alone()) {
            read_mesh("../data/cube_6.MFEM", mesh_master);
            mars::Bisection<mars::Mesh3> b(mesh_master);
            b.uniform_refine(res * 2);
            mesh_master.clean_up();

            mars::write_mesh("mesh_m", mesh_master);
        }

        mars::Mesh3 mesh_slave(true);

        if(comm.rank() == 1 || comm.is_alone()) {
            read_mesh("../data/cube_6.MFEM", mesh_slave);
            mars::Bisection<mars::Mesh3> b(mesh_slave);
            b.uniform_refine(res * 3);
            mesh_slave.clean_up();

            mars::write_mesh("mesh_s", mesh_slave);
        }

        moonolith::logger() << "n_master " << mesh_master.n_elements() << " n_slave " << mesh_slave.n_elements() << " " << std::endl;
        auto vol = mars_compute_isect_measure(comm, mesh_master, mesh_slave);

        if(comm.rank() == 1 || comm.is_alone()) {
            assert(moonolith::approxeq(1., vol, 1e-8));
        }
    }

    static void run_mars_surf_intersect_2()
    {
        moonolith::logger() << "--------------------------------\n";
        moonolith::logger() << "run_mars_surf_intersect_2" << std::endl;

        moonolith::Communicator comm;

        Integer res = 2;

        mars::Mesh2 mesh_master(true);
        mars::Mesh<2, 1> surf_mesh_master(true);

        if(comm.rank() == 0 || comm.is_alone()) {
            read_mesh("../data/square_2.MFEM", mesh_master);
            mars::mark_boundary(mesh_master);

            mars::Bisection<mars::Mesh2> b(mesh_master);
            b.uniform_refine(res);
            mesh_master.clean_up();

            mars::write_mesh("mesh_m", mesh_master);

            mesh_master.update_dual_graph(true);
            mars::extract_surface(mesh_master, surf_mesh_master);
        }

        mars::Mesh2 mesh_slave(true);
        mars::Mesh<2, 1> surf_mesh_slave(true);

        if(comm.rank() == 1 || comm.is_alone()) {
            read_mesh("../data/square_2.MFEM", mesh_slave);
            mars::mark_boundary(mesh_slave);

            mars::Bisection<mars::Mesh2> b(mesh_slave);
            b.uniform_refine(res * 3);
            mesh_slave.clean_up();

            mars::write_mesh("mesh_s", mesh_slave);

            mesh_slave.update_dual_graph(true);

            mars::extract_surface(mesh_slave, surf_mesh_slave);
        }


        moonolith::logger() << "n_master " << surf_mesh_master.n_elements() << " n_slave " << surf_mesh_slave.n_elements() << " " << std::endl;
        double arclen = mars_compute_isect_measure(comm, surf_mesh_master, surf_mesh_slave);


        if(comm.rank() == 1 || comm.is_alone()) {
            assert(moonolith::approxeq(4., arclen, 1e-8));
        }
    }


    static void run_mars_vol_2_surf_intersect_2()
    {
        moonolith::logger() << "--------------------------------\n";
        moonolith::logger() << "run_mars_vol_2_surf_intersect_2" << std::endl;

        moonolith::Communicator comm;

        Integer res = 5;

        mars::Mesh2 mesh_master(true);

        if(comm.rank() == 0 || comm.is_alone()) {
            read_mesh("../data/square_2.MFEM", mesh_master);
            mars::mark_boundary(mesh_master);

            mars::Bisection<mars::Mesh2> b(mesh_master);
            b.uniform_refine(res);
            mesh_master.clean_up();

            mars::write_mesh("mesh_m", mesh_master);

            mesh_master.update_dual_graph(true);
        }

        mars::Mesh<2, 1> surf_mesh_slave(true);

        if(comm.rank() == 1 || comm.is_alone()) {
            surf_mesh_slave.reserve(1, 3);
            surf_mesh_slave.add_point({0.4999999, 0.});
            surf_mesh_slave.add_point({0.500001, 0.5});
            surf_mesh_slave.add_point({0.5000001, 1.});
            surf_mesh_slave.add_elem({0, 1});
            surf_mesh_slave.add_elem({1, 2});
        }

        moonolith::logger() << "n_master " <<  mesh_master.n_active_elements() << "/" << mesh_master.n_elements() << " n_slave " << surf_mesh_slave.n_elements() << " " << std::endl;
        double area = mars_vol2surf_isect_measure(comm, mesh_master, surf_mesh_slave);

        if(comm.rank() == 1 || comm.is_alone()) {
            assert(moonolith::approxeq(1., area, 1e-4));
        }
    }

    static void run_mars_surf_intersect_3()
    {
        moonolith::logger() << "--------------------------------\n";
        moonolith::logger() << "run_mars_surf_intersect_3" << std::endl;

        moonolith::Communicator comm;

        Integer res = 2;

        mars::Mesh3 mesh_master(true);
        mars::Mesh<3, 2> surf_mesh_master(true);

        if(comm.rank() == 0 || comm.is_alone()) {
            read_mesh("../data/cube_6.MFEM", mesh_master);
            mars::mark_boundary(mesh_master);

            mars::Bisection<mars::Mesh3> b(mesh_master);
            b.uniform_refine(res);
            mesh_master.clean_up();

            mars::write_mesh("mesh_m", mesh_master);

            mesh_master.update_dual_graph(true);
            extract_surface(mesh_master, surf_mesh_master);
        }

        mars::Mesh3 mesh_slave(true);
        mars::Mesh<3, 2> surf_mesh_slave(true);

        if(comm.rank() == 1 || comm.is_alone()) {
            read_mesh("../data/cube_6.MFEM", mesh_slave);
            mars::mark_boundary(mesh_slave);

            mars::Bisection<mars::Mesh3> b(mesh_slave);
            b.uniform_refine(res * 3);
            mesh_slave.clean_up();

            mars::write_mesh("mesh_s", mesh_slave);

            mesh_slave.update_dual_graph(true);

            extract_surface(mesh_slave, surf_mesh_slave);

            // surf_mesh_slave.describe(logger());
        }


        moonolith::logger() << "n_master " << surf_mesh_master.n_elements() << " n_slave " << surf_mesh_slave.n_elements() << " " << std::endl;
        double area = mars_compute_isect_measure(comm, surf_mesh_master, surf_mesh_slave);

        if(comm.rank() == 1 || comm.is_alone()) {
            assert(moonolith::approxeq(6., area, 1e-8));
        }
    }


    static void run_mars_vol_2_surf_intersect_3()
    {
        moonolith::logger() << "--------------------------------\n";
        moonolith::logger() << "run_mars_vol_2_surf_intersect_3" << std::endl;

        moonolith::Communicator comm;

        Integer res = 6;

        mars::Mesh3 mesh_master(true);

        if(comm.rank() == 0 || comm.is_alone()) {
            read_mesh("../data/cube_6.MFEM", mesh_master);
            mars::mark_boundary(mesh_master);

            mars::Bisection<mars::Mesh3> b(mesh_master);
            b.uniform_refine(res);
            mesh_master.clean_up();

            mars::write_mesh("mesh_m", mesh_master);

            mesh_master.update_dual_graph(true);
        }

        mars::Mesh<3, 2> surf_mesh_slave(true);
        if(comm.rank() == 1 || comm.is_alone()) {
            surf_mesh_slave.reserve(2, 4);
            surf_mesh_slave.add_point({0.0, 0.0, 0.49999999});
            surf_mesh_slave.add_point({1.0, 0.0, 0.5});
            surf_mesh_slave.add_point({1.0, 1.0, 0.5});
            surf_mesh_slave.add_point({0.0, 1.0, 0.500000001});

            surf_mesh_slave.add_elem({0, 1, 2});
            surf_mesh_slave.add_elem({0, 2, 3});

            mars::Bisection<mars::Mesh<3, 2>> b(surf_mesh_slave);
            b.uniform_refine(res*2);
            surf_mesh_slave.clean_up();
        }

        moonolith::logger() << "n_master " << mesh_master.n_elements() << " n_slave " << surf_mesh_slave.n_elements() << " " << std::endl;
        double area = mars_vol2surf_isect_measure(comm, mesh_master, surf_mesh_slave);
        if(comm.rank() == 1 || comm.is_alone()) {
            assert(moonolith::approxeq(1., area, 1e-6));
        }
    }


    static void run_mars_surf_intersect_4()
    {
        moonolith::logger() << "--------------------------------\n";
        moonolith::logger() << "run_mars_surf_intersect_4" << std::endl;

        moonolith::Communicator comm;

        Integer res = 1;

        mars::Mesh4 mesh_master(true);
        mars::Mesh<4, 3> surf_mesh_master(true);

        if(comm.rank() == 0 || comm.is_alone()) {
            read_mesh("../data/pentatope_1.MFEM", mesh_master);
            mars::mark_boundary(mesh_master);

            mars::Bisection<mars::Mesh4> b(mesh_master);
            b.uniform_refine(res);
            mesh_master.clean_up();

            mesh_master.update_dual_graph(true);
            extract_surface(mesh_master, surf_mesh_master);
        }

        mars::Mesh4 mesh_slave(true);
        mars::Mesh<4, 3> surf_mesh_slave(true);

        if(comm.rank() == 1 || comm.is_alone()) {
            read_mesh("../data/pentatope_1.MFEM", mesh_slave);
            mars::mark_boundary(mesh_slave);

            mars::Bisection<mars::Mesh4> b(mesh_slave);
            b.uniform_refine(res*3);
            mesh_slave.clean_up();

            mesh_slave.update_dual_graph(true);

            extract_surface(mesh_slave, surf_mesh_slave);

            // mesh_slave.describe(logger());
            // surf_mesh_slave.describe(logger());
        }

        moonolith::logger() << "n_master " << surf_mesh_master.n_elements() << " n_slave " << surf_mesh_slave.n_elements() << " " << std::endl;
        double area = mars_compute_isect_measure(comm, surf_mesh_master, surf_mesh_slave);


        if(comm.rank() == 1 || comm.is_alone()) {
            assert(moonolith::approxeq(1., area, 1e-8));
        }
    }

    void single_intersect_2()
    {
        static const Integer Dim = 2;
        static const Integer ManifoldDim = 2;
        using MeshT          = mars::Mesh<Dim, ManifoldDim>;
        using AlgorithmT     = moonolith::SingleCollectionOneMasterOneSlaveAlgorithm<Dim, MeshT>;
        using IntersectDataT = mars::IntersectData<Dim, ManifoldDim>;
        using Adapter = typename AlgorithmT::Adapter;

        moonolith::logger() << "--------------------------------\n";
        moonolith::logger() << "single_intersect_2" << std::endl;

        moonolith::Communicator comm;

        auto res = 1;
        MeshT mesh(true);

        if(comm.rank() == 0 || comm.is_alone()) {
            read_mesh("../data/square_2.MFEM", mesh);
            mars::Bisection<mars::Mesh2> b(mesh);
            b.uniform_refine(res * 2);
            mesh.clean_up();

            //duplicate elements
            auto n_elems = mesh.n_elements();
            for(Integer i = 0; i < n_elems; ++i) {
                mesh.elem(i).block = 1;
                auto id = mesh.add_elem(mesh.elem(i));
                mesh.elem(id).block = 2;
            }
        }

        double elapsed = MPI_Wtime();

        AlgorithmT algo(comm, moonolith::make_unique<moonolith::CollectionManager<MeshT>>());

        algo.init(
            mesh,
            {{1, 2}},
            0.
        );

        IntersectDataT isect;

        double vol = 0.;

        elapsed = MPI_Wtime() - elapsed;
        moonolith::logger() << "init: " << elapsed  << std::endl;

        ////////////////////////////////////////////////////
        /////////////////// pair-wise method ///////////////////////
        elapsed = MPI_Wtime();

        algo.compute([&](const Adapter &master, const Adapter &slave) -> bool {
            auto v = isect.compute(master, slave);

            if(v == 0.) return false;

            vol += v;
            return true;
        });

        comm.all_reduce(&vol, 1, moonolith::MPISum());

        elapsed = MPI_Wtime() - elapsed;
        moonolith::logger() << "time(master, slave): " << elapsed  << std::endl;
        moonolith::logger() << "vol " << vol << std::endl;

        assert(moonolith::approxeq(1., vol, 1e-8));
    }

    void run_mars_moonolith_test()
    {
        // //volume-to-volume
        run_mars_intersect_1();
        run_mars_intersect_2();
        run_mars_intersect_3();
        single_intersect_2();

        //surface-to-surface
        run_mars_surf_intersect_2();
        run_mars_surf_intersect_3();

        // if(Communicator().is_alone()) {
            run_mars_surf_intersect_4();
        // }

        //volume-to-surface
        run_mars_vol_2_surf_intersect_2();
        run_mars_vol_2_surf_intersect_3();


        //others
        run_mars_quadrature_test();
        run_mars_project_test();
        run_mars_l2_transfer_test();
    }

}
