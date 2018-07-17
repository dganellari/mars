#ifndef MARS_PAR_EDGE_SPLIT_POOL_HPP
#define MARS_PAR_EDGE_SPLIT_POOL_HPP

#include "communicator.hpp"

namespace mars {
	template<class Mesh>
	class Bisection;

	template<Integer Dim, Integer ManifoldDim>
	class ParMesh;

	class ParEdgeSplitPool {
	public:

		template<class ParMesh>
		void build_edge_interface(
			ParMesh &mesh,
			Bisection<typename ParMesh::Mesh> &bisection)
		{

		}

		void synchronize() {}

		template<Integer Dim, Integer ManifoldDim>
		void update_midpoint_parts(
			const EdgeNodeMap &edge_node_map,
			ParMesh<Dim, ManifoldDim> &mesh)
		{

		}

		template<Integer Dim, Integer ManifoldDim>
		void update(
			const ParMesh<Dim, ManifoldDim> &mesh,
			const EdgeElementMap &edge_elem_map,
			const EdgeNodeMap &edge_node_map,
			const std::vector<Edge> &bisected_edges)
		{

		}

		template<Integer Dim, Integer ManifoldDim>
		void write_to_mesh(
			const EdgeNodeMap &edge_node_map,
			ParMesh<Dim, ManifoldDim> &mesh)
		{

		}

		template<Integer Dim, Integer ManifoldDim>
		void read_from_mesh(
			const EdgeNodeMap &edge_node_map,
			const ParMesh<Dim, ManifoldDim> &mesh)
		{

		}

		template<Integer Dim, Integer ManifoldDim>
		void collect_splits_to_local_edges(
			const ParMesh<Dim, ManifoldDim> &mesh,
			std::vector<Edge> &splits)
		{

		}

		bool empty() const
		{
			return false;
		}

		ParEdgeSplitPool(const Communicator &comm)
		: comm_(comm)
		{}
		
		Communicator comm_;
	};
}

#endif //MARS_PAR_EDGE_SPLIT_POOL_HPP
