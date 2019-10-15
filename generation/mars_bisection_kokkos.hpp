#ifndef MARS_BISECTION_KOKKOS_HPP
#define MARS_BISECTION_KOKKOS_HPP

#include "mars_dof_map.hpp"
#include "mars_tracker.hpp"
#include "mars_edge_select_kokkos.hpp"
#include "mars_longest_edge.hpp"
#include "mars_edge_node_map_kokkos.hpp"
#include "mars_edge_element_map_kokkos.hpp"
#include <iostream>
#include <typeinfo>

namespace mars {
template<class Mesh_, class EdgeSelect_ = LongestEdgeSelect<Mesh_>>
class ParallelBisection {
public:
	using Mesh     = Mesh_;
	using Elem     = typename Mesh::Elem;
	using SideElem = typename Mesh::SideElem;
	using ParallelEdgeElementMap = SubManifoldElementMap<2,KokkosImplementation>;
	using ParallelEdgeNodeMap = DeviceEdgeNodeMap<2>;
	using Edge = Side<2, KokkosImplementation>;


	using ElementVector = TempArray<Integer,20>;

	static const Integer Dim 		 = Mesh_::Dim;
	static const Integer ManifoldDim = Mesh_::ManifoldDim;

	virtual ~ParallelBisection() {}

	ParallelBisection(Mesh *mesh)
	: mesh(nullptr),
	  host_mesh(mesh),
	  verbose(false),
	  fail_if_not_refine(false)
	{}

	void set_fail_if_not_refine(const bool val) {
		fail_if_not_refine = val;
	}

	bool get_fail_if_not_refine() {
		return fail_if_not_refine;
	}

	/*

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

	//does ont perform a deep copy of the view containted in the mesh. Just the mesh object.
	void copy_mesh_to_device()
	{
		Mesh* tmp = (Mesh*) Kokkos::kokkos_malloc(sizeof(Mesh));
		Mesh mCopy = *host_mesh;
		Mesh* oldDeviceMesh = mesh;
		Kokkos::parallel_for("CreateMeshObject", 1, KOKKOS_LAMBDA (const int&)
		{
			new ((Mesh*)tmp) Mesh(mCopy); // two local copies for m and tmp since this->m, this->mesh host pointers
			if (oldDeviceMesh) oldDeviceMesh->~Mesh();
			// it only works on a copy since **m is still a host pointer and fails on the device.
		});

		mesh = tmp; //make the mesh pointer a device one so that this init func is not neccessary anymore
	}


	Integer euler_graph_formula(Mesh *mesh){

		/*Mesh3 sMesh3;
		read_mesh("../data/write/euler.MFEM", sMesh3);*/

		typename Mesh::SerialMesh sMesh;
		convert_parallel_mesh_to_serial<Dim,ManifoldDim>(sMesh, *mesh);

		Kokkos::Timer timer;

		map_.update(sMesh);

		double time = timer.seconds();
		std::cout << "serial kokkos took: " << time << " seconds." << std::endl;

		//map_.describe(std::cout);

		sMesh.update_dual_graph();


		//nr_faces per eleme counted twicee * nr elements
		Integer bound = Combinations<ManifoldDim + 1, ManifoldDim>::value * sMesh.n_elements();

		//Euler's formula for graphs
		Integer faces_nr = (bound + sMesh.n_boundary_sides())/2; //boundary sides counted only once
		//Integer edges_nr = sMesh3.n_nodes() + faces_nr - 2;

		Integer edges_nr = 0;

		if(ManifoldDim == 2)
			edges_nr = sMesh.n_nodes() + sMesh.n_elements() + 1 -2; // +1 because also the outer face is counted
		else
			edges_nr = sMesh.n_nodes() + faces_nr -1 - sMesh.n_elements();

		// edges_nr = 2 * mesh->n_elements();

		return edges_nr;
	}

	void uniform_refine(const Integer n_levels)
	{

		const Integer nr_elements = host_mesh->n_elements();

		const Integer capacity = euler_graph_formula(host_mesh);
		edge_node_map_.reserve_map(capacity);
		edge_element_map_.reserve_map(capacity);
		reserve_tree(nr_elements);
		/*if(flags.empty()) {
		 flags.resize(mesh.n_elements(), NONE);
		 level.resize(mesh.n_elements(), 0);*/


		copy_mesh_to_device();

		Kokkos::Timer timer1;

		edge_element_map_.update(mesh, nr_elements);

		double time1 = timer1.seconds();
		std::cout << "paralel kokkos took: " << time1 << " seconds." << std::endl;

		//edge_element_map_.describe(std::cout);

		Kokkos::Timer timer2;

		precompute_lepp_incidents(mesh,nr_elements);

		double time2 = timer2.seconds();
		std::cout << "paralel precompute_lepp_incidents took: " << time2 << " seconds." << std::endl;


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
			const Integer edge_num = es.stable_select(*mesh, element_id);
			//Integer edge_num = Bisection<Mesh>::edge_select()->select(Bisection<Mesh>::get_mesh(), element_id);
			//Integer edge_num = Bisection<Mesh>::edge_select()->select(Bisection<Mesh>::get_mesh(), element_id, Bisection<Mesh>::edge_element_map());

			Edge edge;
			mesh->elem(element_id).edge(edge_num, edge.nodes[0], edge.nodes[1]);
			edge.fix_ordering();

			tree(element_id) = mapping.find(edge);
			//ElementVector incidents =mapping.value_at(tree(element_id));

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


	MARS_INLINE_FUNCTION
	static bool is_terminal(const Edge& edge, const ElementVector& incidents, Mesh *mesh)
	{

		bool terminal = true; //in case the elements share the longest edge or there is only one incident (itself)
							  // - meaning the longest edge is on the boundary.

		for (auto i = 0; i < incidents.index; ++i)
		{

			if (mesh->is_active(incidents[i]))
			{

				Edge new_edge;

				EdgeSelect_ es;
				const Integer edge_num = es.stable_select(*mesh,
						incidents[i]);
				//const Integer edge_num = Bisection<Mesh>::edge_select()->select(Bisection<Mesh>::get_mesh(), edge, i);

				mesh->elem(incidents[i]).edge(edge_num, new_edge.nodes[0],
						new_edge.nodes[1]);
				new_edge.fix_ordering();

				//	if(edge != tree.edges[i]){
				if (edge != new_edge)
				{
					terminal = false;
				}
			}
		}

		return terminal;
	}

	MARS_INLINE_FUNCTION
	static bool is_leaf(const Integer element_id, Mesh_ *mesh,
			ViewVectorType<uint32_t> tree,
			const UnorderedMap<Side<2, KokkosImplementation>, ElementVector>& mapping)
	{

		Integer index = tree(element_id);

		for (auto i = 0; i < mapping.value_at(index).index; ++i)
		{
			const auto &element = mapping.value_at(index)[i];
			if (mesh->is_active(element)
					&& mapping.key_at(index) != mapping.key_at(tree(element)))
			{
				return false;
			}
		}

		return true;
	}

	struct RefineMesh
	{
		Mesh_* mesh;
		UnorderedMap<Side<2,KokkosImplementation>,ElementVector> mapping;
		ViewVectorType<uint32_t> tree;
		ViewVectorType<Integer> index_count;


		RefineMesh(UnorderedMap<Side<2, KokkosImplementation>, ElementVector> mp,
				Mesh_* ms, ViewVectorType<uint32_t> tr, ViewVectorType<Integer> ic) :
			mapping(mp), mesh(ms), tree(tr), index_count(ic)
		{
		}


		MARS_INLINE_FUNCTION
		void compute_lepp(const Integer element_id, Integer& count) const
		{

			Integer index = tree(element_id);

			if (is_terminal(mapping.key_at(index), mapping.value_at(index), mesh))
			{
				for (auto i = 0; i < mapping.value_at(index).index; ++i)
				{
					auto &element = mapping.value_at(index)[i];
					if (mesh->is_active(element)){
						++count;
					}
				}
			}
		}


		MARS_INLINE_FUNCTION
		void depth_first(const Integer node, Integer& count) const
		{

			if (is_leaf(node, mesh, tree, mapping))
			{
				compute_lepp(node,count);
			}

			Integer index = tree(node);
			for (auto i = 0; i < mapping.value_at(index).index; ++i)
			{
				const auto &child = mapping.value_at(index)[i];
				if (mesh->is_active(child)
						&& mapping.key_at(index) != mapping.key_at(tree(child)))
					depth_first(child, count);
			}
		}

		MARS_INLINE_FUNCTION
		void operator()(const int element_id) const
		{

			if (mesh->is_active(element_id))
			{
				//Integer terminal_elem_count=0;
				Integer terminal_incident_count=0;

				depth_first(element_id, terminal_incident_count);

				//index_count(element_id,0) = terminal_elem_count;
				index_count(element_id+1) = terminal_incident_count;
				//+1 for leaving the first cell 0 and performing an inclusive scan on the rest
				//to have both exclusive and inclusive and the total on the last cell.

			}
		}
	};


	struct ScatterElem : RefineMesh
	{
		//ViewMatrixType<Integer> lepp_elem_index;

		//ViewMatrixType<Integer> lepp_incidents_index;

		UnorderedMap<Integer, Edge> lepp_incidents_map;

		ScatterElem(
				UnorderedMap<Side<2, KokkosImplementation>, ElementVector> mp,
				Mesh_* ms, ViewVectorType<uint32_t> tr,
				ViewVectorType<Integer> ic,// ViewMatrixType<Integer> lei,
				UnorderedMap<Integer, Edge> lim) :
				RefineMesh(mp, ms, tr, ic),// lepp_elem_index(lei),
				lepp_incidents_map(lim)
		{
		}

		MARS_INLINE_FUNCTION
		void compute_lepp(const Integer element_id) const
		{
			Integer index = this->tree(element_id);

			if (is_terminal(this->mapping.key_at(index), this->mapping.value_at(index),  this->mesh))
			{
				Integer count=0;
			//	lepp_elem_index(index_count(element_id,0)++) = element_id;
				for (auto i = 0; i < this->mapping.value_at(index).index; ++i)
				{
					auto &element = this->mapping.value_at(index)[i];
					if (this->mesh->is_active(element))
					{
						/*lepp_incidents_index(this->index_count(element_id),0) = this->mapping.key_at(index);
						lepp_incidents_index(this->index_count(element_id),0) = index;
						lepp_incidents_index(this->index_count(element_id)++,1) = element;
						printf(" (%li - %li - %li) ",lepp_incidents_index(this->index_count(element_id),0),count, element_id);*/

						auto result = lepp_incidents_map.insert(element, this->mapping.key_at(index));
						if (result.failed())
							printf("Exceeded UnorderedMap: lepp_incidents_map capacity\n");
					}
				}
			}
		}

		MARS_INLINE_FUNCTION
		void depth_first(const Integer node) const
		{
			if (is_leaf(node, this->mesh, this->tree, this->mapping))
			{
				compute_lepp(node);
			}

			Integer index = this->tree(node);
			for (auto i = 0; i < this->mapping.value_at(index).index; ++i)
			{
				const auto &child = this->mapping.value_at(index)[i];
				if (this->mesh->is_active(child)
						&& this->mapping.key_at(index) != this->mapping.key_at(this->tree(child)))
					depth_first(child);
			}
		}

		MARS_INLINE_FUNCTION
		void operator()(const int element_id) const
		{
			if (this->mesh->is_active(element_id))
			{
				depth_first(element_id);
			}
		}
	};


	struct BisectElem
	{
		Mesh_* mesh;
		ViewMatrixType<Integer> lepp_incident_index;
		bool verbose;
		using UEMap = UnorderedMap<Side<2,KokkosImplementation>,ElementVector>;
		using UNMap = UnorderedMap<Edge,Integer>;


		UEMap edge_element_map;
		UNMap edge_node_map;

		Integer elem_start_index;
		Integer child_start_index;
		//Integer node_start_index;
		ViewObject<Integer> node_start_index;

		BisectElem(Mesh_* ms, ViewMatrixType<Integer> lii,
				UEMap mp, UNMap nmp, bool v, Integer esi, Integer csi, ViewObject<Integer> nsi) :
				mesh(ms), lepp_incident_index(lii), edge_element_map(mp), edge_node_map(nmp), verbose(
						v), elem_start_index(esi), child_start_index(csi), node_start_index(nsi)
		{
		}

		MARS_INLINE_FUNCTION
		void other_nodes(const SubView<Integer, ManifoldDim + 1> &nodes,
			const Integer v1, const Integer v2,
			TempArray<Integer, ManifoldDim - 1> &opposite_nodes) const
		{
			Integer i = 0;
			for (Integer j=0;j<   ManifoldDim + 1; ++j)
			{
				Integer n = nodes[j];
				if (n != v1 && n != v2)
				{
					opposite_nodes[i++] = n;
				}
			}
		}

		MARS_INLINE_FUNCTION
		void update_elem(Integer elem_new_id) const
		{
/*
			flags.push_back(NONE);
			level.push_back(
				(mesh.elem(id).parent_id == INVALID_INDEX) ?
					0 : (level[mesh.elem(id).parent_id] + 1));*/
		}

		MARS_INLINE_FUNCTION
		void add_new_elem(const Integer elem_new_id, Elem& old_el,
				const Integer midpoint,const Integer v,
				const TempArray<Integer, ManifoldDim - 1>& opposite_nodes) const
		{

			Elem new_el = mesh->elem(elem_new_id);
			new_el.parent_id = old_el.id;

			for (Integer i = 0; i < ManifoldDim - 1; ++i)
			{
				new_el.nodes[2 + i] = opposite_nodes[i];
			}

			new_el.nodes[0] = midpoint;
			new_el.nodes[1] = v;

			mesh->set_active(elem_new_id);
			ParallelEdgeElementMap::update_elem(edge_element_map, new_el);
			new_el.block = old_el.block;
		}

		MARS_INLINE_FUNCTION
		Integer add_node(const Integer v1,	const Integer v2) const
		{
			TempArray<Integer, 2> nodes;
			nodes[0] = v1;
			nodes[1] = v2;

			Side<2, KokkosImplementation> edge(nodes);

			Integer midpoint=0;

			if (ParallelEdgeNodeMap::update(edge, edge_node_map, node_start_index, midpoint))
			{
				assert(midpoint!=INVALID_INDEX);

				mesh->add_point((mesh->point(v1) + mesh->point(v2)) / 2.,
						midpoint);
			}

			return midpoint;
		}

		MARS_INLINE_FUNCTION
		void operator()(const int i) const
		{
			Integer element_id = lepp_incident_index(i, 0);
			Integer v1 = lepp_incident_index(i, 1);
			Integer v2 = lepp_incident_index(i, 2);

			if (verbose)
			{
				printf("bisect(%li , %li) for %li\n", v1, v2, element_id);
			}

			mesh->set_active(element_id, false);
			//printf("elem: %li - %d\n", element_id, mesh->is_active(element_id));

			Integer nodes_new_id = add_node(v1, v2);

			Integer elem_new_id = elem_start_index + 2 * i;
			Integer childrens_id = child_start_index + i;

			Elem old_el = mesh->elem(element_id, childrens_id);

			TempArray<Integer, ManifoldDim - 1> opposite_nodes;
			other_nodes(old_el.nodes, v1, v2, opposite_nodes);

			add_new_elem(elem_new_id, old_el, nodes_new_id, v1, opposite_nodes);
			old_el.children(0) = elem_new_id;

			add_new_elem(++elem_new_id, old_el, nodes_new_id, v2, opposite_nodes);
			old_el.children(1) = elem_new_id;
		}
	};

	void scan(const Integer start, const Integer end, ViewVectorType<Integer> index_count_)
	{
		using namespace Kokkos;

		//paralel scan on the index_count view for both columns at the same time misusing the complex instead.
	/*	parallel_scan (RangePolicy<>(1 , nr_elements +1 ),	KOKKOS_LAMBDA (const int& i,
					complex<Integer>& upd, const bool& final)
		{*/
		parallel_scan (RangePolicy<>(start , end ),	KOKKOS_LAMBDA (const int& i,
					Integer& upd, const bool& final)
		{
			// Load old value in case we update it before accumulating
			const float val_i = index_count_(i);

			upd += val_i;
			/*upd.real() += val_i0;
			upd.imag() += val_i1*/;

			if (final)
			{
				/*index_count_(i,0) = upd.real(); // only update array on final pass
				index_count_(i,1) = upd.imag();// only update array on final pass1*/
				index_count_(i) = upd;
			}
			// For exclusive scan, change the update value after
			// updating array, like we do here. For inclusive scan,
			// change the update value before updating array.
			/*upd.real() += val_i0;
			upd.imag() += val_i1;*/

		});
	}

	void compact_map_to_view(const UnorderedMap<Integer, Edge>& lepp_incidents_map,
			ViewMatrixType<Integer> lepp_incident_index)
	{
		using namespace Kokkos;

		ViewObject<Integer> global_index("global_index");

		parallel_for(lepp_incidents_map.capacity(), KOKKOS_LAMBDA(uint32_t i)
		{
		  if( lepp_incidents_map.valid_at(i) )
		  {
			Integer k = Kokkos::atomic_fetch_add(&global_index(0), 1);

			lepp_incident_index(k,0) = lepp_incidents_map.key_at(i);
			lepp_incident_index(k,1) = lepp_incidents_map.value_at(i).nodes(0);
			lepp_incident_index(k,2) = lepp_incidents_map.value_at(i).nodes(1);

		  }
		});
	}

	inline bool refine_mesh(const Integer levels, const Integer nr_elements)
	{
		using namespace Kokkos;

		ViewVectorType<Integer> index_count_ = ViewVectorType<Integer>(
				"index_count", nr_elements + 1); //1 more for the total sum

		Kokkos::Timer timer2;

		//for(Integer l = 0; l < levels; ++l)
		parallel_for(nr_elements,
			RefineMesh(edge_element_map_.mapping_, mesh, tree_,
					index_count_));

		double time2 = timer2.seconds();
		std::cout << "paralel RefineMesh took: " << time2 << " seconds." << std::endl;

		Kokkos::Timer timer;

		scan(1, nr_elements +1, index_count_);

		double time = timer.seconds();
		std::cout << "paralel scan took: " << time << " seconds." << std::endl;

		auto sum_subview = subview (index_count_, nr_elements);
		//printf("%s\n", typeid(sum_subview).name());

		auto h_ac = Kokkos::create_mirror_view(sum_subview);

		// Deep copy device view to host view.
		deep_copy(h_ac, sum_subview);

		Integer lepp_incidents_count = h_ac(0);

		UnorderedMap<Integer, Edge> lepp_incidents_map(lepp_incidents_count);


		parallel_for(nr_elements, ScatterElem(edge_element_map_.mapping_, mesh, tree_,
				index_count_, lepp_incidents_map));


		Integer lip_size = lepp_incidents_map.size();
		std::cout<<"lip_size: "<<lip_size<<std::endl;

		Integer elem_start_index = host_mesh->n_elements();
		Integer child_start_index = host_mesh->n_childrens();
	//	Integer node_start_index = host_mesh->n_nodes();

		ViewObject<Integer> node_start_index= ViewObject<Integer>("node_start_index");
		//printf("%s\n", typeid(sum_subview).name());

		auto h_aci = Kokkos::create_mirror_view(node_start_index);
		h_aci(0) =  host_mesh->n_nodes();
		deep_copy(node_start_index , h_aci);


		host_mesh->resize_children(lip_size);
		host_mesh->resize_elements(2 * lip_size);
		host_mesh->resize_active(2 * lip_size);
		host_mesh->resize_points(lip_size);

		host_mesh->set_n_elements(host_mesh->n_elements() + 2 * lip_size);
		host_mesh->set_n_nodes(host_mesh->n_nodes() + lip_size);
		host_mesh->set_n_childrens(host_mesh->n_childrens() + lip_size);


		copy_mesh_to_device(); //only the object. No deep copy of the member views.

		ViewMatrixType<Integer> lepp_incident_index = ViewMatrixType<Integer>(
						"lepp_incidents_index", lip_size, 3);

		compact_map_to_view(lepp_incidents_map, lepp_incident_index);

		//start the bisection kernel starting from host_mesh->n_elements() and host_mesh->get_view_children().extent(0) as initial index

		std::cout<<"Starting bisection"<<std::endl;

		edge_node_map_.rehash_map(2*host_mesh->n_elements()); //TODO: maybe a modified euler for capacity?
		//edge_element_map_.describe(std::cout);
		edge_element_map_.rehash_map(2*host_mesh->n_elements()); //TODO: maybe a modified euler for capacity?
		//edge_element_map_.describe(std::cout);

		Kokkos::Timer timer1;

		parallel_for(lip_size,
				BisectElem(mesh, lepp_incident_index, edge_element_map_.mapping_,
						edge_node_map_.mapping_, verbose, elem_start_index,
						child_start_index, node_start_index));

		double time1 = timer1.seconds();
		std::cout << "paralel bisection took: " << time1 << " seconds." << std::endl;

		/*const Mesh* tmp = ref.mesh;

		 parallel_for("DestroyMeshObject",1, KOKKOS_LAMBDA (const int&) {
		 mesh->~Mesh();
		 });

		 fence();

		 kokkos_free(mesh);
		 */

		return true;
	}


private:
	Mesh *mesh;
	Mesh *host_mesh;
	ParallelEdgeElementMap edge_element_map_;
	EdgeElementMap map_;
	ParallelEdgeNodeMap edge_node_map_;

	ViewVectorType<uint32_t> tree_;

	/*std::vector<Integer> flags;
	std::vector<Integer> level;
	std::vector<std::array<Integer, ManifoldDim+1> > side_flags;




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
