#ifndef MARS_EDGE_SELECT_HPP
#define MARS_EDGE_SELECT_HPP

#include "mars_base.hpp"
#include "mars_edge.hpp"
#include "mars_edge_element_map.hpp"
#include "mars_fwd.hpp"

#include <ostream>

namespace mars {

	template<class Mesh>
	class EdgeSelect {
	public:
		static const Integer Dim = Mesh::Dim;
		static const Integer ManifoldDim = Mesh::ManifoldDim;

		virtual ~EdgeSelect() {}
		virtual Integer select(
			const Mesh &mesh,
			const Integer element_id) const
		{
			//first edge selected
			return 0;
		}

		virtual Integer select(
			const Mesh &mesh,
			const Edge &neighbor_edge,
			const Integer element_id) const
		{
			//first edge selected
			return 0;
		}

		virtual bool repair_element()
		{
			return true;
		}

		virtual bool can_refine(
			const Mesh &mesh,
			const Integer element_id) const
		{
			return true;
		}

		virtual void reorder_edge(
			const Mesh &mesh,
			const Integer element_id,
			Integer &v1,
			Integer &v2) const
		{
		}

		virtual void edge_refined(
			const Mesh &mesh,
			const EdgeElementMap &eem,
			const Edge &edge)
		{
			//
		}

		virtual void update(const Mesh &mesh) {}

		virtual void element_refined(
			const Mesh &mesh,
			const Integer element_id,
			const Edge &edge,
			const Integer local_midpoint_id)
		{
			//do smth
		}

		virtual bool is_recursive() const
		{
			return false;
		}

		virtual void describe(std::ostream &os) const {}

		virtual std::string name() const = 0;
	};

}


#endif //MARS_EDGE_SELECT_HPP
