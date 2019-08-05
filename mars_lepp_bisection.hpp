#ifndef MARS_LEPP_BISECTION_HPP_
#define MARS_LEPP_BISECTION_HPP_

#include "mars_bisection.hpp"

#include <stack>
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

			assert(has_edge(Bisection<Mesh>::get_mesh().elem(i), edge.nodes[0], edge.nodes[1]));

			if (!Bisection<Mesh>::get_mesh().is_active(i) || i == element_id)
				continue;

			if(Bisection<Mesh>::edge_select()->can_refine(Bisection<Mesh>::get_mesh(), i))
				index = i;
		}

		return index;
	}

	/*void refine_element(const Integer element_id) {

		if(!Bisection<Mesh>::edge_select()->can_refine(Bisection<Mesh>::get_mesh(), element_id)) {
			Bisection<Mesh>::get_incomplete_elements().push_back(element_id);
			assert(!Bisection<Mesh>::get_fail_if_not_refine());
			return;
		}

		std::stack<Integer> s;
		s.push(element_id);

		while (!s.empty()){

			Integer top = s.top();

			Integer edge_num = Bisection<Mesh>::edge_select()->select(Bisection<Mesh>::get_mesh(), top);

			Edge edge;
			Bisection<Mesh>::get_mesh().elem(top).edge(edge_num, edge.nodes[0], edge.nodes[1]);
			edge.fix_ordering();

			Integer le_nhb = longest_edge_neighbor(top,edge);

			Edge new_edge;

			//If the longest edge of top is not boundary edge (if there exists a longest edge neighbor).
			if(le_nhb != INVALID_INDEX){

				//get the longest edge of that neighbor.
				edge_num = Bisection<Mesh>::edge_select()->select(Bisection<Mesh>::get_mesh(), le_nhb);
//				edge_num = Bisection<Mesh>::edge_select()->select(Bisection<Mesh>::get_mesh(), edge, le_nhb);

				Bisection<Mesh>::get_mesh().elem(le_nhb).edge(edge_num, new_edge.nodes[0],	new_edge.nodes[1]);
				new_edge.fix_ordering();
			}

			// if top is terminal triangle (if the longest edge is shared or the longest edge is boundary).
			if(le_nhb == INVALID_INDEX || edge == new_edge) {

				//bisect top;
				Bisection<Mesh>::bisect_element(top, edge);

				if(le_nhb != INVALID_INDEX){
					//bisect le_nhb
					Bisection<Mesh>::bisect_element(le_nhb, new_edge);
				}
				s.pop();
			}else
				s.push(le_nhb);
		}

	}*/


	void compute_lepp(std::map<Edge, std::vector<Integer>>& lepp, const Integer element_id){

		Integer edge_num = Bisection<Mesh>::edge_select()->stable_select(Bisection<Mesh>::get_mesh(), element_id);
		Edge edge;
		Bisection<Mesh>::get_mesh().elem(element_id).edge(edge_num, edge.nodes[0], edge.nodes[1]);
		edge.fix_ordering();

		auto incidents = Bisection<Mesh>::edge_element_map().elements(edge);

		 std::vector<Integer> lepp_incidents;
		 std::vector<Integer> lepp_eq;


		if (is_terminal(edge, incidents,lepp_incidents,lepp_eq)) {

			lepp[edge] = lepp_eq;

		} else {

			for (auto i : lepp_incidents) {
				compute_lepp(lepp,i);
			}
		}
	}

	bool is_terminal(const Edge edge, const std::vector<Integer> incidents,
			std::vector<Integer>& lepp_inc, std::vector<Integer>& lepp_eq) {

		bool terminal = true; //in case the elements share the longest edge or there is only one incident (itself)
							 // - meaning the longest edge is on the boundary.

		for (auto i : incidents) {

			if(Bisection<Mesh>::get_mesh().is_active(i) && Bisection<Mesh>::edge_select()->can_refine(Bisection<Mesh>::get_mesh(), i)){

				Edge new_edge;
				const Integer edge_num = Bisection<Mesh>::edge_select()->select(Bisection<Mesh>::get_mesh(), edge, i);
				Bisection<Mesh>::get_mesh().elem(i).edge(edge_num, new_edge.nodes[0], new_edge.nodes[1]);
				new_edge.fix_ordering();

				if(edge != new_edge){
					lepp_inc.push_back(i);
					terminal= false;
				}else
					lepp_eq.push_back(i);
			}
		}

		return terminal;
	}

	void refine_element(const Integer element_id) {

		if (!Bisection<Mesh>::edge_select()->can_refine(Bisection<Mesh>::get_mesh(), element_id)) {
			Bisection<Mesh>::get_incomplete_elements().push_back(element_id);
			assert(!Bisection<Mesh>::get_fail_if_not_refine());
			return;
		}

		while (Bisection<Mesh>::get_mesh().is_active(element_id)) {

			std::map<Edge, std::vector<Integer>> terminal_edges;
			compute_lepp(terminal_edges, element_id);

			for (auto const& star : terminal_edges) {

				for (const Integer element : star.second) {

					if (Bisection<Mesh>::get_mesh().is_active(element)
							&& Bisection<Mesh>::edge_select()->can_refine(Bisection<Mesh>::get_mesh(), element))
						Bisection<Mesh>::bisect_element(element, star.first);
				}
			}

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
