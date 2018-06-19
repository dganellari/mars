#ifndef MARS_EDGE_ELEMENT_MAP_HPP
#define MARS_EDGE_ELEMENT_MAP_HPP

#include <map>
#include <vector>
#include <memory>
#include <cassert>

namespace mars {
	template<Integer Dim, Integer ManifoldDim>
	class Mesh;

	class EdgeElementMap {
	public:
		using ElementVector = std::vector<Integer>;

		class Edge {
		public:
			std::array<Integer, 2> nodes;

			Edge(const Integer a_node, const Integer another_node)
			{
				if(a_node < another_node) {
					nodes[0] = a_node;
					nodes[1] = another_node;
				} else {
					nodes[0] = another_node;
					nodes[1] = a_node;
				}
			}

			inline bool operator<(const Edge &other) const
			{
				if(nodes[0] < other.nodes[0]) {
					return true;
				}

				if(nodes[0] > other.nodes[0]) {
					return false;
				}

				return nodes[1] < other.nodes[1];
			}
		};

		template<Integer Dim, Integer ManifoldDim>
		void update(const Simplex<Dim, ManifoldDim> &e)
		{
			for(Integer ki = 0; ki < n_nodes(e); ++ki) {
				for(Integer kj = ki+1; kj < n_nodes(e); ++kj) {
					auto &vec = mapping_[Edge(e.nodes[ki], e.nodes[kj])];
					if(!vec) {
						vec = std::make_shared<ElementVector>();
					}

					vec->push_back(e.id);
				}
			}
		}

		template<Integer Dim, Integer ManifoldDim>
		void build(const Mesh<Dim, ManifoldDim> &mesh)
		{	
			auto ne = mesh.n_elements();

			for(Integer i = 0; i < ne; ++i) {
				update(mesh.elem(i));
			}
		}

		const ElementVector &elements(const Integer a_node, const Integer another_node) const
		{
			static const ElementVector null_vec;

			auto it = mapping_.find(Edge(a_node, another_node));
			if(it == mapping_.end()) {
				return null_vec;
			}

			assert(it->second);
			return *it->second;
		}

		std::map<Edge, std::shared_ptr<ElementVector> > mapping_;
	};
}

#endif //MARS_EDGE_ELEMENT_MAP_HPP
