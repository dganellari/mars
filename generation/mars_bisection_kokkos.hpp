#ifndef MARS_BISECTION_KOKKOS_HPP
#define MARS_BISECTION_KOKKOS_HPP

#include "mars_dof_map.hpp"
#include "mars_tracker.hpp"
#include "mars_edge_select.hpp"
#include "mars_edge_select_kokkos.hpp"
#include "mars_longest_edge.hpp"
#include <iostream>

namespace mars {
	template<class Mesh_, class EdgeSelect_ = LongestEdgeSelect<Mesh_>>
	class ParallelBisection {
	public:
		using Mesh     = Mesh_;
		using Elem     = typename Mesh::Elem;
		using SideElem = typename Mesh::SideElem;
		using EdgeElementMap1 = SubManifoldElementMap<2,KokkosImplementation>;

		using ElementVector = TempArray<Integer,8>;

		static const Integer Dim 		 = Mesh_::Dim;
		static const Integer ManifoldDim = Mesh_::ManifoldDim;

		virtual ~ParallelBisection() {}

		ParallelBisection(Mesh *mesh)
		: mesh(mesh),
		  verbose(false),
		  fail_if_not_refine(false)
		{}

		void set_fail_if_not_refine(const bool val) {
			fail_if_not_refine = val;
		}

		bool get_fail_if_not_refine() {
			return fail_if_not_refine;
		}

		/*std::vector<Integer>& get_incomplete_elements(){
			return incomplete_elements_;
		}

		virtual void set_edge_select(const std::shared_ptr<EdgeSelect_> &edge_select)
		{
			edge_select_ = edge_select;
		}

		std::shared_ptr<EdgeSelect_> edge_select()
		{
			return edge_select_;
		}



		Integer add_elem(const Elem &e)
		{
			flags.push_back(NONE);
			Integer id = mesh.add_elem(e);
			
			if(edge_select_->repair_element()) {
				mesh.repair_element(id);
			}

			level.push_back((mesh.elem(id).parent_id == INVALID_INDEX)? 0 : (level[mesh.elem(id).parent_id] + 1)); 
			edge_element_map_.update(mesh.elem(id));
			return id;
		}

		void other_nodes(
			const std::array<Integer, ManifoldDim+1> &nodes,
			const Integer v1, 
			const Integer v2,
			std::array<Integer, ManifoldDim-1>  &opposite_nodes) const
		{
			Integer i = 0;
			for(auto n : nodes) {
				if(n != v1 && n != v2) {
					opposite_nodes[i++] = n;
				}
			}
		}

		void bisect_element(
			const Integer element_id,
			const Edge &edge)
		{
			Integer v1 = edge.nodes[0];
			Integer v2 = edge.nodes[1];

			edge_select_->reorder_edge(mesh, element_id, v1, v2);
			bisect_element(element_id, v1, v2);
		}

		// void bisect_element(
		// 	const Integer element_id,
		// 	const Integer v1,
		// 	const Integer v2)
		// {
		// 	mesh.elem(element_id).children.clear();
		// 	mesh.set_active(element_id, false);

		// 	Elem s;
		// 	s.parent_id = element_id;
			
		// 	if(verbose) {
		// 		std::cout << "bisect(" << v1 << ", " << v2 << ") for " << element_id << std::endl;
		// 	}

		// 	auto midpoint = edge_node_map_.get(v1, v2);

		// 	if(midpoint == INVALID_INDEX) {
		// 		// midpoint = mesh.add_point(0.5 * (mesh.point(v1) + mesh.point(v2)));
		// 		midpoint = mesh.add_point((mesh.point(v1) + mesh.point(v2))/2.);
		// 		edge_node_map_.update(v1, v2, midpoint);
		// 	}

		// 	std::array<Integer, ManifoldDim-1> opposite_nodes;
		// 	other_nodes(mesh.elem(element_id).nodes, v1, v2, opposite_nodes);

		// 	for(Integer i = 0; i < ManifoldDim-1; ++i) {
		// 		s.nodes[2+i] = opposite_nodes[i];
		// 	}

		// 	s.nodes[0] = v1;
		// 	s.nodes[1] = midpoint;

		// 	Integer new_id = add_elem(s);
		// 	mesh.elem(element_id).children.push_back(new_id);

		// 	s.nodes[0] = v2;
		// 	s.nodes[1] = midpoint;

		// 	new_id = add_elem(s);
		// 	mesh.elem(element_id).children.push_back(new_id);

		// 	bisect_side_tags(element_id, Edge(v1, v2), midpoint);
		// 	return;
		// }

		void bisect_element(
			const Integer element_id,
			const Integer v1,
			const Integer v2)
		{
			mesh.elem(element_id).children.clear();
			mesh.set_active(element_id, false);

			Elem s;
			s.parent_id = element_id;
			
			if(verbose) {
				std::cout << "bisect(" << v1 << ", " << v2 << ") for " << element_id << std::endl;
			}

			auto midpoint = edge_node_map_.get(v1, v2);

			if(midpoint == INVALID_INDEX) {
				midpoint = mesh.add_point((mesh.point(v1) + mesh.point(v2))/2.);
				edge_node_map_.update(v1, v2, midpoint);
			}

			std::array<Integer, ManifoldDim-1> opposite_nodes;
			other_nodes(mesh.elem(element_id).nodes, v1, v2, opposite_nodes);

			for(Integer i = 0; i < ManifoldDim-1; ++i) {
				s.nodes[2 + i] = opposite_nodes[i];
			}

			s.nodes[0] = midpoint;
			s.nodes[1] = v1;

			Integer new_id = add_elem(s);
			mesh.elem(element_id).children.push_back(new_id);
			mesh.elem(new_id).block = mesh.elem(element_id).block;

			s.nodes[0] = midpoint;
			s.nodes[1] = v2;
			
			new_id = add_elem(s);
			mesh.elem(element_id).children.push_back(new_id);
			mesh.elem(new_id).block = mesh.elem(element_id).block;

			bisect_side_tags(element_id, Edge(v1, v2), midpoint);
			
			tracker_.element_refined(element_id);

			edge_select_->element_refined(
				mesh,
				element_id,
				Edge(v1, v2),
				midpoint
			);



			return;
		}

		inline Integer side_num(
			const Integer element_id,
			const SideElem &side) const
		{
			auto nodes = side.nodes;
			std::sort(nodes.begin(), nodes.end());

			const auto &e = mesh.elem(element_id);

			SideElem e_side;
			
			for(Integer i = 0; i < n_sides(e); ++i) {
				e.side(i, e_side);
				std::sort(std::begin(e_side.nodes), std::end(e_side.nodes));

				bool same_side = true;
				for(Integer k = 0; k < nodes.size(); ++k) {
					assert(nodes[k] != INVALID_INDEX);
					assert(e_side.nodes[k] != INVALID_INDEX);

					if(nodes[k] != e_side.nodes[k]) {
						same_side = false;
						break;
					} 
				}

				if(same_side) {
					return i;
				}
			}

			return INVALID_INDEX;
		}


		void bisect_side_tags(
			const Integer element_id,
			const Edge &edge,
			const Integer midpoint_id)
		{
			const auto &e = mesh.elem(element_id);

			for(auto c : e.children) {
				std::fill(mesh.elem(c).side_tags.begin(),
					mesh.elem(c).side_tags.end(),
					INVALID_INDEX);
			}	

			SideElem side;
			SideElem child_side;

			for(Integer i = 0; i < n_sides(e); ++i) {
				e.side(i, side);
				const Integer tag = e.side_tags[i];
				if(tag == INVALID_INDEX) continue;

				if(has_edge(side, edge[0], edge[1])) {
					//tag of split sides
					child_side.nodes[0] = midpoint_id;
					
					for(Integer j = 0; j < 2; ++j) {
						const Integer vj = edge[j];
						child_side.nodes[1] = vj;

						Integer local_ind = 2;
						for(auto s : side.nodes) {
							if(s == edge[0] || s == edge[1]) {
								continue;
							}
							child_side.nodes[local_ind++] = s;
						}

						bool found_side = false;
						for(auto c : e.children) {
							auto sn = side_num(c, child_side);

							if(INVALID_INDEX != sn) {
								assert(mesh.elem(c).side_tags[sn] == INVALID_INDEX || 
									mesh.elem(c).side_tags[sn] == tag);

								mesh.elem(c).side_tags[sn] = tag;
								found_side = true;
							}
						}

						assert(found_side);
					}
				} else {
					bool found_side = false;

					for(auto c : e.children) {
						auto sn = side_num(c, side);

						if(INVALID_INDEX != sn) {
							assert(mesh.elem(c).side_tags[sn] == INVALID_INDEX || 
								mesh.elem(c).side_tags[sn] == tag);

							mesh.elem(c).side_tags[sn] = tag;
							found_side = true;
						}
					}

					assert(found_side);
				}
			}
		}

		bool refine_element_recursive(
			const Integer element_id,
			const Edge &edge,
			const Integer max_level)
		{
			assert(has_edge(mesh.elem(element_id), edge.nodes[0], edge.nodes[1]));

			const Integer edge_num = edge_select_->select(mesh, edge, element_id);

			Edge new_edge;
			mesh.elem(element_id).edge(edge_num, new_edge.nodes[0], new_edge.nodes[1]);
			new_edge.fix_ordering();

			bool success = true;
			if(edge == new_edge) {
				bisect_element(element_id, edge);
			} else if(!edge_select_->is_recursive()) {
				bisect_element(element_id, edge);
			} else {
				success = refine_edge(new_edge);
				assert(!mesh.is_active(element_id));
			}

			return success;
		}

		virtual void refine_element(const Integer element_id)
		{
			if(!edge_select_->can_refine(mesh, element_id)) {
				incomplete_elements_.push_back(element_id);
				assert(!fail_if_not_refine);
				return;
			}

			//const Integer edge_num = edge_select_->stable_select(mesh, element_id);
			const Integer edge_num = edge_select_->select(mesh, element_id);
			Edge edge;
			mesh.elem(element_id).edge(edge_num, edge.nodes[0], edge.nodes[1]);
			edge.fix_ordering();
			refine_edge(edge);
		}

		bool refine_edge(const Edge &edge)
		{
			bool visited_all_incidents = false;
			bool has_refined = false;
			bool success = true;
			bool all_done = true;
			auto incidents = edge_element_map_.elements(edge);

			while(!visited_all_incidents) {
				visited_all_incidents = true;

				for(auto i : incidents) {
					if(!mesh.is_active(i)) continue;

					assert(has_edge(mesh.elem(i), edge.nodes[0], edge.nodes[1]));
					all_done = false;

					const bool can_refine = edge_select_->can_refine(mesh, i);
					if(can_refine) {
						if(refine_element_recursive(i, edge, level[i])) {
							visited_all_incidents = false;
							has_refined = true;
							continue;
						} 

						assert(!fail_if_not_refine);
					}

					assert(!fail_if_not_refine);
					success = false;
				}

				const auto &next_incidents = edge_element_map_.elements(edge);

				if(next_incidents.size() > incidents.size()) {
					incidents = next_incidents;
					visited_all_incidents = false;
				}
			}

			if(!success) {
				incomplete_edges_.push_back(edge);
			}

			if(has_refined) {
				bisected_edges_.push_back(edge);
			}

			return success || all_done;
		}

		void uniform_refine(const Integer n_levels)
		{
			if(flags.empty()) {
				flags.resize(mesh.n_elements(), NONE);
				level.resize(mesh.n_elements(), 0);
				edge_element_map_.update(mesh);
				mesh.update_dual_graph();
			}

			for(Integer l = 0; l < n_levels; ++l) {
				const Integer n_elements = mesh.n_elements();
				for(Integer i = 0; i < n_elements; ++i) {
					if(!mesh.is_valid(i)) {
						std::cerr << "tried to refine non-existing element " << i << std::endl;
						continue;
					}

					if(!mesh.is_active(i)) continue;
					refine_element(i);
				}
			}

			mesh.update_dual_graph();
			mesh.tags() = level;
		}

		void refine(const std::vector<Integer> &elements)
		{
			std::vector<Integer> sorted_elements = elements;
			// std::sort(sorted_elements.begin(), sorted_elements.end());

			if(flags.empty()) {
				flags.resize(mesh.n_elements(), NONE);
				level.resize(mesh.n_elements(), 0);
				edge_element_map_.update(mesh);
				mesh.update_dual_graph();
			}

			for(auto i : sorted_elements) {
				if(!mesh.is_valid(i)) {
					std::cerr << "tried to refine non-existing element " << i << std::endl;
					continue;
				}

				if(!mesh.is_active(i)) {
					continue;
				}

				refine_element(i);
			}

			mesh.update_dual_graph();
			mesh.tags() = level;
		}

		const EdgeNodeMap &edge_node_map() const
		{
			return edge_node_map_;
		}

		const EdgeElementMap &edge_element_map() const
		{
			return edge_element_map_;
		}

		EdgeElementMap &edge_element_map()
		{
			return edge_element_map_;
		}

		void set_verbose(const bool val)
		{
			verbose = val;
		}

		const Mesh &get_mesh() const
		{
			return mesh;
		}

		Mesh &get_mesh()
		{
			return mesh;
		}

		const std::vector<Edge> &bisected_edges() const
		{
			return bisected_edges_;
		}

		void clear_bisected_edges()
		{
			bisected_edges_.clear();
		}

		void refine_edges(const std::vector<Edge> &edges)
		{
			if(flags.empty()) {
				flags.resize(mesh.n_elements(), NONE);
				level.resize(mesh.n_elements(), 0);
				edge_element_map_.update(mesh);
				mesh.update_dual_graph();
			}

			for(auto e : edges) {
				refine_edge(e);
			}

			mesh.update_dual_graph();
			mesh.tags() = level;
		}

		bool refine_incomplete()
		{
			std::cout << "incomplete elems: " << incomplete_elements_.size();
			std::cout << " edges: " << incomplete_edges_.size() << std::endl;

			auto elements_to_refine = std::move(incomplete_elements_);

			incomplete_elements_.clear();
			refine(elements_to_refine);

			auto edges_to_refine = std::move(incomplete_edges_);
			incomplete_edges_.clear();

			refine_edges(edges_to_refine);
			return incomplete_elements_.empty() && incomplete_edges_.empty();
		}

		void clear()
		{
			flags.clear();
			level.clear();
			side_flags.clear();
			edge_node_map_.clear();
			edge_element_map_.clear();
		}

		void tracking_begin()
		{
			tracker_.begin_iterate();
		}

		void tracking_end()
		{
			tracker_.end_iterate();
		}

		void undo()
		{
			tracker_.undo_last_iterate(mesh);
		}*/

		void init(Mesh** m)
		{

			Mesh* tmp = (Mesh*) Kokkos::kokkos_malloc(sizeof(Mesh));
			Mesh mCopy = **m;
			Kokkos::parallel_for("CreateMeshObject", 1, KOKKOS_LAMBDA (const int&)
			{
				new ((Mesh*)tmp) Mesh(mCopy); // it only works on a copy since **m is still a host pointer and fails on the device.
			});

			*m = tmp; //make the mesh pointer a device one so that this init func is not neccessary anymore
		}


		void uniform_refine(const Integer n_levels)
		{

			/*Mesh2 sMesh3;
			read_mesh("../data/square_2_def.MFEM", sMesh3);*/

			Mesh3 sMesh3;
			convert_parallel_mesh_to_serial<Dim,ManifoldDim>(sMesh3, *mesh);
			sMesh3.update_dual_graph();
			Kokkos::Timer timer;

			map_.update(sMesh3);
			double time = timer.seconds();
			std::cout << "serial kokkos took: " << time << " seconds." << std::endl;

			map_.describe(std::cout);
			std::cout<<"size: "<<map_.mapping_.size()<<std::endl;

		    Integer bound = Combinations<ManifoldDim + 1, ManifoldDim>::value * sMesh3.n_elements();

			//Euler's formula for graphs
			Integer faces_nr = bound - sMesh3.n_boundary_sides();
			//Integer edges_nr = sMesh3.n_nodes() + faces_nr - 2;
			Integer edges_nr = sMesh3.n_nodes() + faces_nr -1 - sMesh3.n_active_elements();

            bound = 2 * mesh->n_elements();


            const Integer nr_elements = mesh->n_elements();

            init(&mesh);
            //	careful: now the mesh is a device pointer. Not possible to call on the host for ex. mesh.nr_elements()


			edge_element_map_.reserve_map(bound);
			reserve_tree(nr_elements);
			/*if(flags.empty()) {
				flags.resize(mesh.n_elements(), NONE);
				level.resize(mesh.n_elements(), 0);*/

			std::cout<<"updating"<< std::endl;
            Kokkos::Timer timer1;

			edge_element_map_.update(mesh, nr_elements);
				/*mesh.update_dual_graph();
			}*/
			Kokkos::fence();

			double time1 = timer1.seconds();
			std::cout << "paralel kokkos took: " << time1 << " seconds." << std::endl;


			edge_element_map_.describe(std::cout);

			precompute_lepp_incidents(mesh,nr_elements);
/*

			std::cout<<"cap: "<<edge_element_map_.mapping_.capacity()<<std::endl;
			std::cout<<"size: "<<edge_element_map_.mapping_.size()<<std::endl;
			std::cout<<"E: "<<edges_nr<< " F: "<< faces_nr <<std::endl;
			std::cout<<"C: "<<sMesh3.n_elements()<< " V: "<< sMesh3.n_nodes() <<std::endl;
			std::cout<<"sMesh3.n_boundary_sides(): "<<sMesh3.n_boundary_sides()<< " bound: "<< bound <<std::endl;
			std::cout<<"ManifoldDim: "<<ManifoldDim<< " bound: "<< bound <<std::endl;
*/


			refine_mesh(n_levels,nr_elements);

		/*	mesh.update_dual_graph();
			mesh.tags() = level;*/
		}


		struct BuildTree
		{
			Mesh_* mesh;
			UnorderedMap<Side<2,KokkosImplementation>,ElementVector> mapping;
			ViewVectorType<uint32_t> tree;


			BuildTree(UnorderedMap<Side<2, KokkosImplementation>, ElementVector> mp,
					Mesh_* ms, ViewVectorType<uint32_t> tr) :
				mapping(mp), mesh(ms), tree(tr)
			{
			}

			MARS_INLINE_FUNCTION
			void operator()(int element_id) const
			{
				EdgeSelect_ es;
				Integer edge_num = es.select(*mesh, element_id);
				//Integer edge_num = Bisection<Mesh>::edge_select()->select(Bisection<Mesh>::get_mesh(), element_id);
				//Integer edge_num = Bisection<Mesh>::edge_select()->select(Bisection<Mesh>::get_mesh(), element_id, Bisection<Mesh>::edge_element_map());

				Edge edge;
				mesh->elem(element_id).edge(edge_num, edge.nodes[0], edge.nodes[1]);
				edge.fix_ordering();

			//	tree(element_id) = mapping_.find(edge);

			}
		};

		void reserve_tree(Integer size)
		{
			tree_ = ViewVectorType<uint32_t>("tree", size);
		}

		void precompute_lepp_incidents(
				Mesh_ *mesh,
				const Integer nr_elements)
		{

			Kokkos::parallel_for(nr_elements,
					BuildTree(edge_element_map_.mapping_, mesh, tree_));

		}

	struct RefineMesh
	{
		Mesh* mesh;

		/*void init(const Mesh& m)
		{

			Mesh* tmp = (Mesh*) Kokkos::kokkos_malloc(sizeof(Mesh));

			Kokkos::parallel_for("CreateMeshObject", 1, KOKKOS_LAMBDA (const int&)
					{
						new ((Mesh*)tmp) Mesh(m); //same as this.mesh if mesh instead of tmp and this is a host pointer.
					});

			mesh = tmp; //otherwise "this" pointer is still a host pointer within the parallel_for.
		}*/

		RefineMesh(Mesh *m) : mesh(m)
		{
			//init(m);
		}


		KOKKOS_INLINE_FUNCTION
		void operator()(int element_id) const
		{

			/*while(mesh->is_active(element_id)){


			}*/

			if (mesh->is_active(element_id))
			{
				EdgeSelect_ es;
				es.select(*mesh, element_id);
			}
		}
	};

	inline bool refine_mesh(const Integer levels, const Integer nr_elements)
	{
		using namespace Kokkos;

		RefineMesh ref= RefineMesh(mesh);

		//for(Integer l = 0; l < levels; ++l)
		parallel_for(nr_elements, ref);

		/*const Mesh* tmp = ref.mesh;

		parallel_for("DestroyMeshObject",1, KOKKOS_LAMBDA (const int&) {
			tmp->~Mesh();
		});*/

		return true;
	}


	private:
		Mesh *mesh;
		EdgeElementMap1 edge_element_map_;
		EdgeElementMap map_;

		ViewVectorType<uint32_t> tree_;

		/*std::vector<Integer> flags;
		std::vector<Integer> level;
		std::vector<std::array<Integer, ManifoldDim+1> > side_flags;
		EdgeNodeMap edge_node_map_;




		//tracking the refinement 
		std::vector<Edge> bisected_edges_;

		std::vector<Edge>    incomplete_edges_;
		std::vector<Integer> incomplete_elements_;
		Tracker tracker_;*/
		bool verbose;
		bool fail_if_not_refine;

	};
}

#endif //MARS_BISECTION_HPP