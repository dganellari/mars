#ifndef MARS_PAR_EDGE_SPLIT_POOL_HPP
#define MARS_PAR_EDGE_SPLIT_POOL_HPP

#include "communicator.hpp"

namespace mars {
	template<Integer Dim, Integer ManifoldDim>
	class Bisection;

	template<Integer Dim, Integer ManifoldDim>
	class ParMesh;

	class ParEdgeSplitPool {
	public:

		template<Integer Dim, Integer ManifoldDim>
		void build_edge_interface(
			ParMesh<Dim, ManifoldDim> &mesh,
			Bisection<Dim, ManifoldDim> &bisection)
		{

		}

		void synchronize() {}

		template<Integer Dim, Integer ManifoldDim>
		void update_midpoint_parts(EdgeNodeMap &edge_node_map, ParMesh<Dim, ManifoldDim> &mesh)
		{

		}

		template<Integer Dim, Integer ManifoldDim>
		void update(
			ParMesh<Dim, ManifoldDim> &mesh,
			EdgeElementMap &edge_elem_map,
			EdgeNodeMap &edge_node_map,
			std::vector<Edge> &bisected_edges)
		{

		}

		template<Integer Dim, Integer ManifoldDim>
		void write_to_mesh(
			EdgeNodeMap &edge_node_map,
			ParMesh<Dim, ManifoldDim> &mesh)
		{

		}

		template<Integer Dim, Integer ManifoldDim>
		void read_from_mesh(
			EdgeNodeMap &edge_node_map,
			ParMesh<Dim, ManifoldDim> &mesh)
		{

		}

		template<Integer Dim, Integer ManifoldDim>
		void collect_splits_to_local_edges(
			ParMesh<Dim, ManifoldDim> &mesh,
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
