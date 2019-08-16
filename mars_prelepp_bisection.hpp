#ifndef MARS_PRELEPP_BISECTION_HPP_
#define MARS_PRELEPP_BISECTION_HPP_

#include "mars_bisection.hpp"

#include <stack>
#include <err.h>
#include <omp.h>

#define select_err 1

namespace mars {

template<typename T>
class UnitElementTree {

public:

	std::vector<Edge> edges;
	//std::vector<T> incidents;

	std::vector<std::vector<T>> incidents;

	/*UnitElementTree() =default;

	UnitElementTree(Integer size){

		edges = std::vector<Edge>(size);
		incidents = std::vector<T>(size * 10, -1);
	}*/

};

//LEPP re-computation becomes the bottleneck of this algorithm.
template<class Mesh_>
class PrecomputeLeppBisection: public Bisection<Mesh_> {

public:
	using Mesh = Mesh_;

	static UnitElementTree<Integer> tree;

	int max =0;

	PrecomputeLeppBisection(Mesh &mesh) :
			Bisection<Mesh>(mesh) {
	}


	void refine_element(const Integer element_id) {

		if (!Bisection<Mesh>::edge_select()->can_refine(Bisection<Mesh>::get_mesh(), element_id)) {
			Bisection<Mesh>::get_incomplete_elements().push_back(element_id);
			assert(!Bisection<Mesh>::get_fail_if_not_refine());
			return;
		}

		while (Bisection<Mesh>::get_mesh().is_active(element_id)) {

			//UnitElementTree<Integer> tree;//(Bisection<Mesh>::get_mesh().n_active_elements());

			precompute_lepp_incidents(Bisection<Mesh>::get_mesh());

			compute_lepp(element_id);
		}

		//std::cout<<"max: "<<max<<std::endl;
	}

private:

	void precompute_lepp_incidents(Mesh &mesh){

		tree.edges= std::vector<Edge>(mesh.n_elements());
		tree.incidents = std::vector<std::vector<Integer>>(mesh.n_elements());

		for(Integer k=0;k<mesh.n_elements();++k){

			if(!mesh.is_active(k)) continue;


			Integer edge_num = Bisection<Mesh>::edge_select()->stable_select(Bisection<Mesh>::get_mesh(), k);
			//Integer edge_num = Bisection<Mesh>::edge_select()->select(Bisection<Mesh>::get_mesh(), element_id);
			//Integer edge_num = Bisection<Mesh>::edge_select()->select(Bisection<Mesh>::get_mesh(), element_id, Bisection<Mesh>::edge_element_map());

			Edge edge;
			Bisection<Mesh>::get_mesh().elem(k).edge(edge_num, edge.nodes[0], edge.nodes[1]);
			edge.fix_ordering();

			tree.edges[k] = edge;

			auto incidents = Bisection<Mesh>::edge_element_map().elements(edge);

			tree.incidents[k] = std::move(incidents);

			/*int j=0;
			for (auto i : tree.incidents[k]) {

				if (Bisection<Mesh>::get_mesh().is_active(i)
						&& Bisection<Mesh>::edge_select()->can_refine(
								Bisection<Mesh>::get_mesh(), i)) {
					if(j<10){
						tree.incidents[k][j] = i;
						++j;
					}
					else{
						std::cout<<"more than 10 el: "<<k<<std::endl;
					}
					++j;
				}
			}

			if(j>max)
				max=j;*/
		}


	}

	void set_edge_select(
			const std::shared_ptr<EdgeSelect<Mesh_>> &edge_select) final {

		if(edge_select->name() != "LongestEdge")
			errx(select_err,"%s %d %s Error: Calling set_edge_select on LeppBisection with %s is not supported!",
					__FILE__,__LINE__, __func__, edge_select->name().c_str());
	}

	bool is_leaf(const Integer element_id){

		for (const Integer element : tree.incidents[element_id]) {
			if(Bisection<Mesh>::get_mesh().is_active(element) && tree.edges[element_id] != tree.edges[element]){
				return false;
			}
		}

		return true;
	}

	void compute_lepp(const Integer element_id){

	/*	Integer edge_num = Bisection<Mesh>::edge_select()->stable_select(Bisection<Mesh>::get_mesh(), element_id);
		//Integer edge_num = Bisection<Mesh>::edge_select()->select(Bisection<Mesh>::get_mesh(), element_id);
		//Integer edge_num = Bisection<Mesh>::edge_select()->select(Bisection<Mesh>::get_mesh(), element_id, Bisection<Mesh>::edge_element_map());

		Edge edge;
		Bisection<Mesh>::get_mesh().elem(element_id).edge(edge_num, edge.nodes[0], edge.nodes[1]);
		edge.fix_ordering();

		auto incidents = Bisection<Mesh>::edge_element_map().elements(edge);*/

		//auto edge = tree.edges[element_id]; auto incidents =  tree.incidents[element_id];
		//if(is_leaf(element_id)){

			if (is_terminal(tree.edges[element_id], tree.incidents[element_id])) {

				for (const Integer element : tree.incidents[element_id]) {

					if (Bisection<Mesh>::get_mesh().is_active(element)
							&& Bisection<Mesh>::edge_select()->can_refine(Bisection<Mesh>::get_mesh(), element))
						Bisection<Mesh>::bisect_element(element, tree.edges[element_id]);
				}
			}
		//}
	}

	bool is_terminal(const Edge& edge, const std::vector<Integer>& incidents) {

		bool terminal = true; //in case the elements share the longest edge or there is only one incident (itself)
							 // - meaning the longest edge is on the boundary.

		for (auto i : incidents) {

			if(Bisection<Mesh>::get_mesh().is_active(i) && Bisection<Mesh>::edge_select()->can_refine(Bisection<Mesh>::get_mesh(), i)){

			/*	Edge new_edge;
				const Integer edge_num = Bisection<Mesh>::edge_select()->stable_select(Bisection<Mesh>::get_mesh(), i);
				Bisection<Mesh>::get_mesh().elem(i).edge(edge_num, new_edge.nodes[0], new_edge.nodes[1]);
				new_edge.fix_ordering();*/

				if(edge != tree.edges[i]){
		//		if(edge != new_edge){

					terminal= false;
					compute_lepp(i);
				}
			}
		}

		return terminal;
	}

};

template<class Mesh_>
UnitElementTree<Integer> PrecomputeLeppBisection<Mesh_>::tree;

}

#endif /* MARS_LEPP_BISECTION_HPP_ */
