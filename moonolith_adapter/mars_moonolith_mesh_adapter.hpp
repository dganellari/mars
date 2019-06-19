#ifndef MOONOLITH_MARS_MESH_ADAPTER_HPP
#define MOONOLITH_MARS_MESH_ADAPTER_HPP 

#include "mars.hpp"
#include "moonolith_mesh_adapter.hpp"
#include "moonolith_polygon.hpp"
#include "moonolith_line.hpp"
#include "moonolith_check_stream.hpp"
#include "moonolith_affine_transform.hpp"

#include "moonolith_elem_segment.hpp"
#include "moonolith_elem_triangle.hpp"
#include "moonolith_elem_tetrahedron.hpp"

#include <numeric>
#include <vector>

namespace mars {
	
	template<class Mesh>
	class MeshAllocator {};

	template<Integer Dim, Integer ManifoldDim>
	class MeshAllocator< mars::Mesh<Dim, ManifoldDim> > {
	public:
		inline static std::unique_ptr< mars::Mesh<Dim, ManifoldDim> > make(const Integer type)
		{
			(void)type;
			return moonolith::make_unique< mars::Mesh<Dim, ManifoldDim> >();
		}
	};

	template<Integer Dim>
	class MeshAllocator< mars::IMesh<Dim> > {
	public:
		inline static std::unique_ptr< mars::IMesh<Dim> > make(const Integer type)
		{
			switch(type)
			{
				case 1: { return moonolith::make_unique< mars::Mesh<Dim, 1> >(); }
				case 2: { return moonolith::make_unique< mars::Mesh<Dim, 2> >(); }
				case 3: { return moonolith::make_unique< mars::Mesh<Dim, 3> >(); }
				case 4: { return moonolith::make_unique< mars::Mesh<Dim, 4> >(false); }
				default: {assert(false); return nullptr; }
			}
		}
	};

	template<class Mesh>
	class BlowUpDir {
	public:
		template<class Elem, class Vector>
		static void build(const Mesh &mesh, const Elem&e, Vector &normal)
		{
			normal.set(1.);
		}
	};


	template<Integer Dim>
	class BlowUpDir< mars::Mesh<Dim, Dim-1> > {
	public:
		static void build(const mars::Mesh<Dim, Dim-1> &mesh, const mars::Simplex<Dim, Dim-1> &e, mars::Vector<double, Dim> &normal)
		{
			normal = mars::normal(e, mesh.points());
		}
	};

	template<Integer Dim>
	class BlowUpDir<mars::IMesh<Dim>> {
	public:
		static void build(const mars::IMesh<Dim> &mesh, const mars::IElem &e, mars::Vector<double, Dim> &normal)
		{
			auto m_ptr = dynamic_cast<const mars::Mesh<Dim, Dim-1> *>(&mesh);
			auto e_ptr = dynamic_cast<const mars::Simplex<Dim, Dim-1> *>(&e);

			if(m_ptr && e_ptr) {
				normal = mars::normal(*e_ptr, m_ptr->points());
			} else {
				normal.set(1.);
			}
		}
	};


	template<class Mesh>
	class MarsCollectionManager {
	public:
		using Elem  = typename Mesh::Elem;
		using Point = typename Mesh::Point;
		static const Integer Dim = Mesh::Dim;

		template<class Adapter>
		static void data_added(const Integer local_idx, const Adapter &a)
		{
			(void)local_idx;
			(void)a;
		}

		static const Elem &elem(const Mesh &mesh, const Integer handle) 
		{
			return mesh.elem(handle);
		}

		static Integer tag(const Mesh &mesh, const Integer handle)
		{
			return mesh.elem(handle).block;
		}

		static Integer n_elements(const Mesh &mesh)
		{
			return mesh.n_elements();
		}

		static Integer elements_begin(const Mesh &mesh)
		{
			return 0;
		}

		static Integer elements_end(const Mesh &mesh)
		{
			return mesh.n_elements();
		}

		static Integer handle(const Mesh &mesh, const Integer element_index)
		{
			return element_index;
		}

		static bool skip(const Mesh &mesh, const Integer &element_index)
		{
			return !mesh.is_active(element_index);
		}

		template<class Bound>
		static void fill_bound(const Mesh &mesh, const Integer handle, Bound &bound, const double blow_up)
		{
			const auto &e = mesh.elem(handle);

			if(e.type() < Dim) {
				Point p, q;
				
				mars::Vector<double, Dim> nn;
				BlowUpDir<Mesh>::build(mesh, e, nn);

				for(Integer i = 0; i < n_nodes(e); ++i) {
					
					auto n = e.node(i);
					q = mesh.point(n);

					for(Integer d = 0; d < Dim; ++d) {
						p[d] = q[d] + blow_up * nn[d];
					}

					bound += p.values;

					for(Integer d = 0; d < Dim; ++d) {
						p[d] = q[d] - blow_up * nn[d];
					}

					bound += p.values;
				}

			} else {

				for(Integer i = 0; i < n_nodes(e); ++i) {
					auto n = e.node(i);
					bound += mesh.point(n).values;
				}
			}
		}

		template<class Iterator>
		static void serialize(const Mesh &mesh, const Iterator &begin, const Iterator &end, moonolith::OutputStream &os)
		{
			CHECK_STREAM_WRITE_BEGIN("serialize", os);


			Integer n_elements = std::distance(begin, end);

			os << mesh.type();
			os << n_elements;

			std::vector<Integer> node_index(mesh.n_nodes(), 0);

			for(auto it = begin; it != end; ++it) {
				const auto &e = mesh.elem(*it);

				for(Integer i = 0; i < n_nodes(e); ++i) {
					node_index[e.node(i)] = 1; 
				}
			}

			Integer n_nodes = std::accumulate(std::begin(node_index), std::end(node_index), 0);

			os << n_nodes;

			Integer cum_sum = 0;
			for(auto &n : node_index) {
				if(n == 1) {
					cum_sum += n;
					n = cum_sum;
				}
			}

			// logger() << "outgoing: " << n_elements << ", " << n_nodes << std::endl;

			std::for_each(std::begin(node_index), std::end(node_index), [](Integer &v) {
				v -= 1;
			});

			//write points
			for(Integer i = 0; i < mesh.n_nodes(); ++i) {
				if(node_index[i] >= 0) {
					const auto &p = mesh.point(i);
					
					for(Integer d = 0; d < Dim; ++d) {
						os << p[d];
					}
				}
			}

			//write element connectivity
			for(auto it = begin; it != end; ++it) {
				const auto &e = mesh.elem(*it);
				Integer block = e.get_block();
				os << block;

				for(Integer i = 0; i < mars::n_nodes(e); ++i) {
					auto n = e.node(i);
					assert(node_index[n] >= 0);
					os << node_index[n];
				}
			}

			CHECK_STREAM_WRITE_END("serialize", os);
		}

		static std::unique_ptr<Mesh> build(moonolith::InputStream &is)
		{

			CHECK_STREAM_READ_BEGIN("serialize", is);

			Integer type, n_elements, n_nodes;
			is >> type;
			is >> n_elements;
			is >> n_nodes;

			auto m = MeshAllocator<Mesh>::make(type);

			// logger() << "incoming: " << n_elements << ", " << n_nodes << std::endl;

			m->reserve(n_elements, n_nodes);

			Point p;
			for(Integer i = 0; i < n_nodes; ++i) {
				for(Integer d = 0; d < Dim; ++d) {
					is >> p[d];
				}

				m->add_point(p);
			}

			//FIXME make it so that this is surely checked in case we something else that simplices
			std::vector<Integer> nodes(type + 1);

			for(Integer i = 0; i < n_elements; ++i) {
				Integer block = -1;
				is >> block;

				for(auto &n : nodes) {
					is >> n;
				}

				auto idx = m->add_elem(nodes);
				m->elem(idx).set_block(block);
			}

			CHECK_STREAM_READ_END("serialize", is);
			return m;
		}

	};
	}

	namespace moonolith {
		template<Integer Dim, Integer ManifoldDim>
		class CollectionManager< mars::Mesh<Dim, ManifoldDim> > : public mars::MarsCollectionManager< mars::Mesh<Dim, ManifoldDim> > {};


		template<Integer Dim>
		class CollectionManager< mars::IMesh<Dim> > : public mars::MarsCollectionManager< mars::IMesh<Dim> > {};
	}

	namespace mars {
	inline void make_aux(const mars::Mesh<2, 1> &mesh, const mars::Mesh<2, 1>::Elem &elem, moonolith::Line<double, 2> &poly)
	{	
		const auto &p0 = mesh.point(elem.nodes[0]);
		const auto &p1 = mesh.point(elem.nodes[1]);

		poly.p0.x = p0[0];
		poly.p0.y = p0[1];

		poly.p1.x = p1[0];
		poly.p1.y = p1[1];
	}

	inline void make(const mars::Mesh<2, 1> &mesh, const mars::Mesh<2, 1>::Elem &elem, moonolith::Line<double, 2> &poly)
	{
		make_aux(mesh, elem, poly);
	}

	inline void make(const mars::Mesh<1, 1> &mesh, const mars::Mesh<1, 1>::Elem &elem, moonolith::Line<double, 1> &poly)
	{
		const auto &p0 = mesh.point(elem.nodes[0]);
		const auto &p1 = mesh.point(elem.nodes[1]);

		poly.p0.x = p0[0];
		poly.p1.x = p1[0];
	}

	inline void make_aux(const mars::Mesh2 &mesh, const mars::Mesh2::Elem &elem, moonolith::Polygon<double, 2> &poly)
	{
		poly.resize(mars::n_nodes(elem));

		for(Integer i = 0; i < mars::n_nodes(elem); ++i) {
			const  auto &p = mesh.point(elem.nodes[i]);
			poly[i].x = p[0];
			poly[i].y = p[1];
		}
	}

	inline void make(const mars::Mesh2 &mesh, const mars::Mesh2::Elem &elem, moonolith::Polygon<double, 2> &poly)
	{
		make_aux(mesh, elem, poly);
	}

	inline void make(const mars::IMesh<2> &mesh, const mars::IMesh<2>::Elem &elem, moonolith::Polygon<double, 2> &poly)
	{
		auto mesh_ptr = dynamic_cast<const mars::Mesh2 *>(&mesh);
		auto elem_ptr = dynamic_cast<const mars::Mesh2::Elem *>(&elem);

		if(!mesh_ptr || !elem_ptr) {
			assert(false);
			std::cerr << "[Error] unknown mesh type" << std::endl;
			return;
		}

		make_aux(*mesh_ptr, *elem_ptr, poly);
	}

	inline void make(const mars::IMesh<2> &mesh, const mars::IMesh<2>::Elem &elem, moonolith::Line<double, 2> &poly)
	{
		auto mesh_ptr = dynamic_cast<const mars::Mesh<2, 1> *>(&mesh);
		auto elem_ptr = dynamic_cast<const mars::Mesh<2, 1>::Elem *>(&elem);

		if(!mesh_ptr || !elem_ptr) {
			assert(false);
			std::cerr << "[Error] unknown mesh type" << std::endl;
			return;
		}

		make_aux(*mesh_ptr, *elem_ptr, poly);
	}

	inline void make(const mars::Mesh<3, 2> &mesh, const mars::Mesh<3, 2>::Elem &elem, moonolith::Polygon<double, 3> &poly)
	{
		poly.resize(mars::n_nodes(elem));

		for(Integer i = 0; i < mars::n_nodes(elem); ++i) {
			const auto &p = mesh.point(elem.node(i));

			poly.points[i].x = p[0];
			poly.points[i].y = p[1];
			poly.points[i].z = p[2];
		}
	}

	inline void make(const mars::IMesh<3> &mesh, const mars::IMesh<3>::Elem &elem, moonolith::Polygon<double, 3> &poly)
	{
		auto mesh_ptr = dynamic_cast<const mars::Mesh<3, 2> *>(&mesh);
		auto elem_ptr = dynamic_cast<const mars::Mesh<3, 2>::Elem *>(&elem);

		if(!mesh_ptr || !elem_ptr) {
			assert(false);
			std::cerr << "[Error] unknown mesh type" << std::endl;
			return;
		}

		poly.resize(mars::n_nodes(elem));

		for(Integer i = 0; i < mars::n_nodes(elem); ++i) {
			const auto &p = mesh.point(elem.node(i));

			poly.points[i].x = p[0];
			poly.points[i].y = p[1];
			poly.points[i].z = p[2];
		}
	}

	inline void make_aux(const mars::Mesh3 &mesh, const mars::Mesh3::Elem &elem, moonolith::Polyhedron<double> &poly)
	{
		poly.el_ptr.resize(4 + 1);
		poly.el_index.resize(12);
		poly.points.resize(4);

		poly.type = moonolith::Polyhedron<double>::TET;

		poly.el_ptr[0] = 0;
		poly.el_ptr[1] = 3;
		poly.el_ptr[2] = 6;
		poly.el_ptr[3] = 9;
		poly.el_ptr[4] = 12;

		for(Integer i = 0; i < mars::n_nodes(elem); ++i) {
			const  auto &p = mesh.point(elem.nodes[i]);
			poly.points[i].x = p[0];
			poly.points[i].y = p[1];
			poly.points[i].z = p[2];
		}

        poly.el_index[0] = 0;
        poly.el_index[1] = 2;
        poly.el_index[2] = 1;
   
        poly.el_index[3] = 0;
        poly.el_index[4] = 3;
        poly.el_index[5] = 2;
   
        poly.el_index[6] = 0;
        poly.el_index[7] = 1;
        poly.el_index[8] = 3;
   
        poly.el_index[9] = 1;
        poly.el_index[10] = 2;
        poly.el_index[11] = 3;
		       
	}

	
	inline void make(const mars::Mesh3 &mesh, const mars::Mesh3::Elem &elem, moonolith::Polyhedron<double> &poly)
	{
		make_aux(mesh, elem, poly);
	}

	inline void make(const mars::IMesh<3> &mesh, const mars::IMesh<3>::Elem &elem, moonolith::Polyhedron<double> &poly)
	{
		auto mesh_ptr = dynamic_cast<const mars::Mesh3 *>(&mesh);
		auto elem_ptr = dynamic_cast<const mars::Mesh3::Elem *>(&elem);

		if(!mesh_ptr || !elem_ptr) {
			assert(false);
			std::cerr << "[Error] unknown mesh type" << std::endl;
			return;
		}

		make_aux(*mesh_ptr, *elem_ptr, poly);
	}

	inline void make(const mars::Mesh<4, 3> &mesh, const mars::Mesh<4, 3>::Elem &elem, moonolith::Polyhedron4<double> &poly)
	{
		poly.el_ptr.resize(4 + 1);
		poly.el_index.resize(12);
		poly.points.resize(4);

		poly.type = moonolith::Polyhedron4<double>::TET;

		poly.el_ptr[0] = 0;
		poly.el_ptr[1] = 3;
		poly.el_ptr[2] = 6;
		poly.el_ptr[3] = 9;
		poly.el_ptr[4] = 12;

		for(Integer i = 0; i < mars::n_nodes(elem); ++i) {
			const  auto &p = mesh.point(elem.nodes[i]);
			poly.points[i].x = p[0];
			poly.points[i].y = p[1];
			poly.points[i].z = p[2];
			poly.points[i].w = p[3];
		}

        poly.el_index[0] = 0;
        poly.el_index[1] = 2;
        poly.el_index[2] = 1;
   
        poly.el_index[3] = 0;
        poly.el_index[4] = 3;
        poly.el_index[5] = 2;
   
        poly.el_index[6] = 0;
        poly.el_index[7] = 1;
        poly.el_index[8] = 3;
   
        poly.el_index[9] = 1;
        poly.el_index[10] = 2;
        poly.el_index[11] = 3;
		       
	}

    inline void make_transform(
        const mars::Mesh<3, 2> &mesh,
        const mars::Mesh<3, 2>::Elem &elem,
        moonolith::AffineTransform<double, 2, 3> &trafo)
    {
        const auto &q0 = mesh.point(elem.nodes[0]);
        const auto &q1 = mesh.point(elem.nodes[1]);
        const auto &q2 = mesh.point(elem.nodes[2]);

        moonolith::Vector<double, 3> p0, p1, p2;

        p0.x = q0[0];
        p0.y = q0[1];
        p0.z = q0[2];

        p1.x = q1[0];
        p1.y = q1[1];
        p1.z = q1[2];

        p2.x = q2[0];
        p2.y = q2[1];
        p2.z = q2[2];

        make(p0, p1, p2, trafo);
    }

    inline void make_transform(
        const mars::Mesh<2, 1> &mesh,
        const mars::Mesh<2, 1>::Elem &elem,
        moonolith::AffineTransform<double, 1, 2> &trafo)
    {
        const auto &q0 = mesh.point(elem.nodes[0]);
        const auto &q1 = mesh.point(elem.nodes[1]);

        moonolith::Vector<double, 2> p0, p1;

        p0.x = q0[0];
        p0.y = q0[1];

        p1.x = q1[0];
        p1.y = q1[1];

        make(p0, p1, trafo);
    }

    inline void make(const mars::Mesh1 &mesh, const mars::Mesh1::Elem &elem, moonolith::Segment<double, 1, 1> &poly)
    {
    	for(Integer i = 0; i < mars::n_nodes(elem); ++i) {
    		const  auto &p = mesh.point(elem.nodes[i]);
    		poly.node(i).x = p[0];
    	}
    }

    inline void make(const mars::Mesh2 &mesh, const mars::Mesh2::Elem &elem, moonolith::Triangle<double, 1, 2> &poly)
    {
    	for(Integer i = 0; i < mars::n_nodes(elem); ++i) {
    		const  auto &p = mesh.point(elem.nodes[i]);
    		poly.node(i).x = p[0];
    		poly.node(i).y = p[1];
    	}
    }

    inline void make(const mars::Mesh3 &mesh, const mars::Mesh3::Elem &elem, moonolith::Tetrahedron<double, 1, 3> &poly)
    {
    	for(Integer i = 0; i < mars::n_nodes(elem); ++i) {
    		const  auto &p = mesh.point(elem.nodes[i]);
    		poly.node(i).x = p[0];
    		poly.node(i).y = p[1];
    		poly.node(i).z = p[2];
    	}
    }
}

#endif //MOONOLITH_MARS_MESH_ADAPTER_HPP
