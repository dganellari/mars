#ifndef MARS_EDGE_ELEMENT_MAP_HPP
#define MARS_EDGE_ELEMENT_MAP_HPP

#include "edge.hpp"

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
		void update_active(const Mesh<Dim, ManifoldDim> &mesh)
		{
			auto ne = mesh.n_elements();
			for(Integer i = 0; i < ne; ++i) {
				if(mesh.is_active(i)) {
					update(mesh.elem(i));
				}
			}
		}


		template<Integer Dim, Integer ManifoldDim>
		void build(const Mesh<Dim, ManifoldDim> &mesh)
		{	
			mapping_.clear();
			
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
