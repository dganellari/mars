#ifndef MARS_BISECTION_HPP
#define MARS_BISECTION_HPP

namespace mars {
	template<Integer Dim, Integer ManifoldDim>
	class Mesh;

	template<Integer Dim>
	void bisect_element(
		Mesh<Dim, 2> &mesh,
		const Integer element_id,
		EdgeNodeMap &enm,
		const Edge &edge
		// const std::vector<Edge> &edges)
	{
		if(edges.empty()) return;

		mesh.elem(element_id).children.clear();

		Simplex<Dim, 2> s;
		s.parent_id = element_id;

		if(edges.size() == 1) {
			// const Integer v1 = edges[0].nodes[0];
			// const Integer v2 = edges[0].nodes[1];

			const Integer v1 = edge.nodes[0];
			const Integer v2 = edges.nodes[1];

			auto midpoint = enm.get(v1, v2);
			auto opposite = mesh.elem(element_id).opposite(v1, v2);

			if(midpoint == INVALID_INDEX) {
				midpoint = mesh.add_point(0.5 * (mesh.point(v1) + mesh.point(v2)));
				enm.update(v1, v2, midpoint);
			}

			s.nodes[0] = opposite;
			s.nodes[1] = midpoint;
			s.nodes[2] = v1;

			mesh.elem(element_id).children.push_back( mesh.add_elem(s) );

			s.nodes[0] = opposite;
			s.nodes[1] = v2;
			s.nodes[2] = midpoint;

			mesh.elem(element_id).children.push_back( mesh.add_elem(s) );
			return;
		}

		// if(edges.size() == 2) {
		// 	const Integer va1 = edges[0].nodes[0];
		// 	const Integer va2 = edges[0].nodes[1];

		// 	const Integer vb1 = edges[1].nodes[0];
		// 	const Integer vb2 = edges[1].nodes[1];

		// 	auto midpoint_a = enm.get(va1, va2);
		// 	auto opposite_a = mesh.elem(element_id).opposite(va1, va2);
		// 	if(midpoint_a == INVALID_INDEX) {
		// 		midpoint_a = mesh.add_point(0.5 * (mesh.point(va1) + mesh.point(va2)));
		// 		enm.update(va1, va2, midpoint_a);
		// 	}

		// 	auto midpoint_b = enm.get(vb1, vb2);
		// 	auto opposite_b = mesh.elem(element_id).opposite(vb1, vb2);
		// 	if(midpoint_b == INVALID_INDEX) {
		// 		midpoint_b = mesh.add_point(0.5 * (mesh.point(vb1) + mesh.point(vb2)));
		// 		enm.update(vb1, vb2, midpoint_b);
		// 	}

		// 	s.nodes[0] = midpoint_a;
		// 	s.nodes[1] = opposite_b;
		// 	s.nodes[2] = opposite_a;

		// 	mesh.elem(element_id).children.push_back( mesh.add_elem(s) );

		// 	s.nodes[0] = midpoint_a;
		// 	s.nodes[1] = midpoint_b;
		// 	s.nodes[2] = va1;

		// 	mesh.elem(element_id).children.push_back( mesh.add_elem(s) );

		// 	child.nodes[0] = midpoint_a;
		// 	child.nodes[1] = vb1;
		// 	child.nodes[2] = midpoint_b;

		// 	mesh.elem(element_id).children.push_back( add_elem(s) );
		// 	return; 
		// }
	}

	template<Integer Dim, Integer ManifoldDim>
	class Bisection {
	public:
		Bisection(Mesh<Dim, ManifoldDim> &mesh)
		: mesh(mesh)
		{}

		void uniformly_refine(const n_levels)
		{

		}

		void refine(const std::vector<Integer> &elems)
		{
			if(flags.empty()) {
				flags.resize(mesh.n_elements(), NONE);
			}

			for(auto i : elements) {
				const auto &e = mesh.elem(i);

				//FIXME
				Integer edge_num = n_edges(e) * double(rand())/RAND_MAX;
				Integer n1, n2;
				e.edge(edge_num, n1, n2);

			}
		}
		
	private:
		Mesh<Dim, ManifoldDim> &mesh;
		std::vector<Integer> flags;
		std::vector<std::array<Integer, ManifoldDim+1> > side_flags;
		EdgeNodeMap edge_node_map_;
	};
}

#endif //MARS_BISECTION_HPP
