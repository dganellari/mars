#include "mars_moonolith_project_test.hpp"

#include "mars.hpp"

#include "moonolith_mesh_adapter.hpp"
#include "moonolith_polygon.hpp"
#include "moonolith_intersect_polygons.hpp"
#include "moonolith_intersect_polyhedra.hpp"
#include "moonolith_intersect_lines.hpp"
#include "mars_moonolith_mesh_adapter.hpp"
#include "moonolith_intersect_polygon_with_h_polyhedron.hpp"
#include "moonolith_intersect_polyline_with_convex_polygon.hpp"
#include "moonolith_project_convex_polygons.hpp"
#include "moonolith_project_polygons.hpp"
#include "moonolith_contact.hpp"
#include "moonolith_matlab_scripter.hpp"
#include "mars_extract_surface.hpp"
#include "mars_moonolith_function_space_adapter.hpp"

#include <cassert>

namespace mars {
	template<int Dim, int ManifoldDim>
	class ProjectData {
	public:
        using Trafo = moonolith::AffineTransform<double, Dim-1, Dim>;
        
        moonolith::AffineContact<double, Dim> contact;
        std::shared_ptr<Trafo> trafo_m, trafo_s;

        ProjectData()
        {
            trafo_m = std::make_shared<Trafo>();
            trafo_s = std::make_shared<Trafo>(); 

            contact.trafo_master = trafo_m;
            contact.trafo_slave  = trafo_s;

            contact.q_rule.points  = { { 0., 0.} };
            contact.q_rule.weights = { 1.};
        }

		template<class Adapter>
		double compute(Adapter &master, Adapter &slave)
		{
            auto &e_m = master.elem();
            auto &e_s = slave.elem();

            auto &m_m = master.collection();
            auto &m_s = slave.collection();

            make(m_m, e_m, contact.master);
            make(m_s, e_s, contact.slave);

            make_transform(m_m, e_m, *trafo_m);
            make_transform(m_s, e_s, *trafo_s);

            if(contact.compute()) {
                return contact.q_slave.weights[0] * moonolith::measure(contact.slave);
            }

			return 0.;
		}
	};
    

	template<Integer Dim, Integer ManifoldDim>
	static double mars_compute_project_measure(
		moonolith::Communicator &comm,
		const mars::Mesh<Dim, ManifoldDim> &mesh_master,
		const mars::Mesh<Dim, ManifoldDim> &mesh_slave,
		const double blow_up)
	{
		using MeshT = mars::Mesh<Dim, ManifoldDim>;
		using ProjectDataT = mars::ProjectData<Dim, ManifoldDim>;
		using Adapter = typename moonolith::ManyMastersOneSlaveAlgorithm<Dim, MeshT>::Adapter;

		double elapsed = MPI_Wtime();

		auto cm = std::make_shared<moonolith::CollectionManager<MeshT>>();
		moonolith::ManyMastersOneSlaveAlgorithm<Dim, MeshT> algo(comm, cm);
		

		algo.init_simple(
			mesh_master,
			mesh_slave,
			blow_up
		);

		ProjectDataT isect;

		double vol = 0.;


		elapsed = MPI_Wtime() - elapsed;
		moonolith::logger() << "init: " << elapsed  << std::endl;
		
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


    static void run_mars_inclined_elem_test_2()
    {
        using Trafo = moonolith::AffineTransform<double, 1, 2>;

        moonolith::logger() << "run_mars_inclined_elem_test_2" << std::endl;
        mars::Mesh<2, 1> mesh;
        mesh.reserve(2, 4);

        std::vector<Integer> elem_master{{
            mesh.add_point({-0.25, -0.566987}),
            mesh.add_point({-0.566987, -0.678606 })
        }};

        std::vector<Integer> elem_slave{{
            mesh.add_point({-0.272727,   -0.5}),
            mesh.add_point({-0.0909091 , -0.5 })
        }};

        mesh.add_elem(elem_master);
        mesh.add_elem(elem_slave);

        // mesh.describe(std::cout);

        moonolith::AffineContact<double, 2> contact;

        std::shared_ptr<Trafo> trafo_m, trafo_s;

        trafo_m = std::make_shared<Trafo>();
        trafo_s = std::make_shared<Trafo>(); 

        contact.q_rule.points  = { { 0.}, {1.} };
        contact.q_rule.weights = { 0.5, 0.5};

        make(mesh, mesh.elem(0), contact.master);
        make(mesh, mesh.elem(1), contact.slave);

        make_transform(mesh, mesh.elem(0), *trafo_m);
        make_transform(mesh, mesh.elem(1), *trafo_s);

        contact.trafo_master = trafo_m;
        contact.trafo_slave  = trafo_s;

        if(contact.compute()) {
            // moonolith::MatlabScripter script;
            // script.close_all();
            // script.hold_on();
            // script.plot(contact.master, "b.-");
            // script.plot(contact.slave, "r.-");
            // script.plot(contact.q_master_physical.points, "b*");
            // script.plot(contact.q_slave_physical.points,  "r*");
            // script.axis_equal();
            // script.save("bug.m");
        }
    }

	static void run_mars_surf_intersect_3()
	{
        moonolith::logger() << " project test " << std::endl;
		moonolith::Communicator comm;

		Integer res = 1;

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
			mars::extract_surface(mesh_master, surf_mesh_master);
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
			
			mars::extract_surface(mesh_slave, surf_mesh_slave);

			// surf_mesh_slave.describe(logger());
		}

		//shift mesh in x direction
		for(Integer i = 0; i < surf_mesh_slave.n_nodes(); ++i) {
			surf_mesh_slave.point(i)(0) += 1.1;
		}

		moonolith::logger() << "n_master " << surf_mesh_master.n_elements() << " n_slave " << surf_mesh_slave.n_elements() << " " << std::endl;
		double area = mars_compute_project_measure(comm, surf_mesh_master, surf_mesh_slave, 0.1001);

		if(comm.rank() == 1 || comm.is_alone()) {
			assert(moonolith::approxeq(1., area, 1e-8));
		}
	}

    void mars_contact_test()
    {
        using MasterElem = moonolith::Triangle<double, 1, 3>;
        using SlaveElem  = moonolith::Quad<double, 1, 3>;

        //dummy quad rule
        moonolith::Quadrature2<double> q_rule;
        q_rule.points  = { { 0., 0.}, {1.0, 0.0}, {0.0, 1.0} };
        q_rule.weights = { 1./3., 1./3., 1./3.};

        MasterElem master;
        SlaveElem slave;

        master.make_reference();
        slave.make_reference();
        
        moonolith::ContactMortar<MasterElem, SlaveElem> algo;

        //always initialize both rules
        algo.set_quadrature_affine(q_rule);
        algo.set_quadrature_warped(q_rule);

        if(algo.assemble(master, slave)) {
            assert(false);
        } 

        slave.set_affine(false);

        if(algo.assemble(master, slave)) {
            assert(false);
        } 

        invert_orientation(slave);

        slave.set_affine(true);

        if(!algo.assemble(master, slave)) {
            assert(false);
        } 
        assert(moonolith::approxeq(0.5, algo.projection_measure(), 1e-5));

        slave.set_affine(false);
        
        algo.clear();
        if(!algo.assemble(master, slave)) {
            assert(false);
        } 

        assert(moonolith::approxeq(0.5, algo.projection_measure(), 1e-5));
        // algo.describe(logger());
    }

    void mars_contact_polymorphic_test()
    {
        using MasterElem = moonolith::Triangle<double, 1, 3>;
        using SlaveElem  = moonolith::Quad<double, 1, 3>;

        using PolymorphicElem = moonolith::Elem<double, 2, 3>;

        //dummy quad rule
        moonolith::Quadrature2<double> q_rule;
        q_rule.points  = { { 0., 0.}, {1.0, 0.0}, {0.0, 1.0} };
        q_rule.weights = { 1./3., 1./3., 1./3.};

        MasterElem master;
        SlaveElem slave;

        master.make_reference();
        slave.make_reference();
        
        moonolith::ContactMortar<PolymorphicElem, PolymorphicElem> algo;

        //always initialize both rules
        algo.set_quadrature_affine(q_rule);
        algo.set_quadrature_warped(q_rule);

        if(algo.assemble(master, slave)) {
            assert(false);
        } 

        slave.set_affine(false);

        if(algo.assemble(master, slave)) {
            assert(false);
        } 

        invert_orientation(slave);

        slave.set_affine(true);

        if(!algo.assemble(master, slave)) {
            assert(false);
        } 
        assert(moonolith::approxeq(0.5, algo.projection_measure(), 1e-5));

        slave.set_affine(false);
        
        algo.clear();
        if(!algo.assemble(master, slave)) {
            assert(false);
        } 

        assert(moonolith::approxeq(0.5, algo.projection_measure(), 1e-5));
        // algo.describe(logger());
    }

	void run_mars_project_test()
	{
        run_mars_inclined_elem_test_2();
		run_mars_surf_intersect_3();
        mars_contact_test();
        mars_contact_polymorphic_test();
	}
}
