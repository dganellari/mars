#ifndef MARS_LEPP_BISECTION_HPP_
#define MARS_LEPP_BISECTION_HPP_
#include <stack>
#include "mars_bisection.hpp"
#include <err.h>

#define select_err 1

namespace mars {

template<class Mesh_>
class LeppBisection: public Bisection<Mesh_> {

public:
	using Mesh = Mesh_;

	LeppBisection(Mesh &mesh) :
			Bisection<Mesh>(mesh) {
	}

	Integer longest_edge_neighbor(Integer element_id, Edge edge) {

		auto incidents = Bisection<Mesh>::edge_element_map().elements(edge);

		Integer index = INVALID_INDEX; //in case a boundary longest edge

		for (auto i : incidents) {
			if (!Bisection<Mesh>::get_mesh().is_active(i) || i == element_id)
				continue;

			assert(has_edge(Bisection<Mesh>::get_mesh().elem(i), edge.nodes[0], edge.nodes[1]));

			index = i;
		}

		return index;
	}

	void refine_element(const Integer element_id) {

		std::stack<Integer> s;
		s.push(element_id);

		while (!s.empty()){

			Integer top = s.top();

			const Integer edge_num = Bisection<Mesh>::edge_select()->select(Bisection<Mesh>::get_mesh(), top);

			Edge edge;
			Bisection<Mesh>::get_mesh().elem(element_id).edge(edge_num, edge.nodes[0], edge.nodes[1]);
			edge.fix_ordering(); //todo: check if needed at all

			Integer le_nhb = longest_edge_neighbor(top,edge);

			edge_num = edge_select_->select(Bisection<Mesh>::get_mesh(), edge, le_nhb);
			Edge new_edge;
			mesh.elem(element_id).edge(edge_num, new_edge.nodes[0],	new_edge.nodes[1]);
			new_edge.fix_ordering();

			if(le_nhb == INVALID_INDEX || edge == new_edge) { // if top is terminal triangle

				//bisect top;
				bisect_element(top, edge);

				if(le_nhb != INVALID_INDEX){
					//bisect le_nhb
					bisect_element(le_nhb, new_edge);
				}
				s.pop();
			}else
				s.push(le_nhb);

			/*if(terminal_triangle(top)){

				const Integer edge_num = edge_select_->select(mesh, top);
				Edge edge;
				mesh.elem(element_id).edge(edge_num, edge.nodes[0],
						edge.nodes[1]);
				edge.fix_ordering();
				bisect_element(top,edge);

				if(le_nhb != INVALID_INDEX){

					//the edge is the same. TODO: may ommit recalc of edge.
					const Integer edge_num = edge_select_->select(mesh, le_nhb);
					Edge edge;
					mesh.elem(element_id).edge(edge_num, edge.nodes[0],
							edge.nodes[1]);
					edge.fix_ordering();
					bisect_element(le_nhb, edge);
				}
				s.pop();

			}else*/
			//s.push(le_nhb);
		}

	}

private:
	void set_edge_select(
			const std::shared_ptr<EdgeSelect<Mesh_>> &edge_select) final {

		if(edge_select->name() != "LongestEdge")
			errx(select_err,"%s %d %s Error: Calling set_edge_select on LeppBisection with %s is not supported!",
					__FILE__,__LINE__, __func__, edge_select->name().c_str());
	}

};
}

#endif /* MARS_LEPP_BISECTION_HPP_ */
