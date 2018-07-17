#ifndef PAR_MESH_HPP
#define PAR_MESH_HPP

#include "mesh.hpp"
#include "communicator.hpp"

namespace mars {

	template<Integer Dim_, Integer ManifoldDim_ = Dim_>
	class ParMesh {
	public:
		static const Integer Dim 		 = Dim_;
		static const Integer ManifoldDim = ManifoldDim_;

		using Mesh     = mars::Mesh<Dim, ManifoldDim>;
		using Point    = mars::Vector<Real, Dim>;
		using Elem     = mars::Simplex<Dim, ManifoldDim>;
		using SideElem = mars::Simplex<Dim, ManifoldDim-1>; 

		inline Integer add_point(
			const Point &point,
			const Integer global_id = INVALID_INDEX,
			const Integer owner     = INVALID_INDEX)
		{
			Integer local_id = mesh_.add_point(point);
			node_map_.add(local_id, global_id, owner);
			return local_id;
		}

		inline Integer add_elem(
			const Elem &elem,
			const Integer global_id = INVALID_INDEX,
			const Integer owner     = INVALID_INDEX)
		{
			Integer local_id = mesh_.add_elem(elem);
			elem_map_.add(local_id, global_id, owner);
			return local_id;
		}

		Integer n_owned_nodes() const
		{
			return node_map_.n_owned();
		}

		inline Integer max_global_node_id() const
		{
			return node_map_.max_id();
		}

		inline Integer max_global_elem_id() const
		{
			return elem_map_.max_id();
		}

		bool is_elem_interface(const Integer local_element_id) const
		{
			const auto &e = mesh_.elem(local_element_id);
			for(auto st : e.side_tags) {
				if(is_partition_tag(st)) return true;
			}

			return false;
		}

		void mark_partition_boundary(
			const Integer local_element_id,
			const Integer side_num,
			const Integer adj_partition)
		{
			auto &e = mesh_.elem(local_element_id);
			e.side_tags[side_num] = partition_id_to_tag(adj_partition);
		}

		void reserve_elements(const Integer n_elems) {}

		void reserve_nodes(const Integer n_nodes) {}

		inline Integer n_active_elements() const
		{
			return n_active_elements_;
		}	

		inline Integer n_nodes() const
		{
			return n_nodes_;
		}	

		inline Integer n_local_active_elements() const
		{
			return mesh_.n_active_elements();
		}

		inline Integer n_local_elements() const
		{
			return mesh_.n_elements();
		}

		inline Integer n_local_nodes() const
		{
			return mesh_.n_nodes();
		}	

		ParMesh(const Communicator &comm)
		: comm_(comm),
		  node_map_(comm.rank(), comm.size()),
		  elem_map_(comm.rank(), comm.size()),
		  is_interfaced_(comm.size(), false),
		  n_active_elements_(0), n_nodes_(0)
		{}

		void init(
			Mesh &serial_mesh,
			const std::vector<Integer> &partitioning)
		{
			serial_mesh.update_dual_graph();

			std::vector<Integer> node_partitioning(serial_mesh.n_nodes(), INVALID_INDEX);
			std::vector<std::set<Integer>> shared_nodes(serial_mesh.n_nodes());

			for(Integer i = 0; i < serial_mesh.n_elements(); ++i) {
				if(!serial_mesh.is_active(i)) continue;
				const Integer partition_id = partitioning[i];
				const auto &e = serial_mesh.elem(i);

				if(partition_id == comm_.rank()) {
					const Integer local_element_id = this->add_elem(e, e.id, comm_.rank());

					const auto &adj = serial_mesh.dual_graph().adj(i);
					for(Integer f = 0; f < adj.size(); ++f) {
						if(adj[f] == INVALID_INDEX) continue;

						const Integer adj_partition_id = partitioning[adj[f]];
						if(adj_partition_id != partition_id) {

							this->mark_partition_boundary(
								local_element_id, f, adj_partition_id
							);

							is_interfaced_[adj_partition_id] = true;
						}
					}
				}

				for(Integer k = 0; k < mars::n_nodes(e); ++k) {
					if(node_partitioning[e.nodes[k]] == INVALID_INDEX) {
						node_partitioning[e.nodes[k]] = partition_id;
					}

					shared_nodes[e.nodes[k]].insert(partition_id);
				}
			}

			this->add_and_index_nodes(serial_mesh, node_partitioning);

			for(Integer i = 0; i < this->n_local_nodes(); ++i) {
				const auto global_n = this->node_map().global(i);
				
				for(auto pid : shared_nodes[global_n]) {
					this->node_map().add_partition(i, pid);
				}
			}

			this->n_active_elements_ = serial_mesh.n_elements();
			this->n_nodes_           = serial_mesh.n_nodes();
		}

		void synchronize()
		{
			std::vector<Integer> offsets;
			node_map().pack_for_global(offsets, true);
			comm_.all_reduce(&offsets[0], offsets.size(), MPIMax());
			node_map().unpack_for_global(offsets);

			elem_map().pack_for_global(offsets, true);
			comm_.all_reduce(&offsets[0], offsets.size(), MPIMax());
			elem_map().unpack_for_global(offsets);

			std::array<Integer, 2> buff = { 
				mesh_.n_active_elements(), 
				node_map().n_owned()
			};

			comm_.all_reduce(&buff[0], buff.size(), MPISum());

			n_active_elements_ = buff[0];
			n_nodes_           = buff[1];
		}

		void clean_up()
		{
			//IMPLEMENT ME
		}

		inline const Map &node_map() const
		{
			return node_map_;
		}
		
		inline const Map &elem_map() const
		{
			return elem_map_;
		}

		inline Map &node_map()
		{
			return node_map_;
		}
		
		inline Map &elem_map()
		{
			return elem_map_;
		}

		void describe(std::ostream &os) const
		{
			serial_apply(comm_, [this, &os]() {
				os << "===================\n";
				os << "mesh " << comm_ << "\n";
				mesh_.describe(os);
				node_map_.describe(os);
				elem_map_.describe(os);
				os << "===================" << std::endl;
			});
		}

		Mesh &get_serial_mesh()
		{
			return mesh_;
		}

		const Mesh &get_serial_mesh() const
		{
			return mesh_;
		}

		inline Communicator &comm()
		{
			return comm_;
		}

		inline const Communicator &comm() const
		{
			return comm_;
		}

	private:

		static bool is_partition_tag(const Integer tag)
		{
			return tag < INVALID_INDEX;
		}

		static Integer partition_id_to_tag(const Integer adj_partition)
		{
			return -(adj_partition + 2);
		}

		static Integer tag_to_partition_id(const Integer tag)
		{
			return -(tag + 2);
		}

		Communicator comm_;
		Mesh mesh_;
		Map node_map_;
		Map elem_map_;
		std::vector<bool> is_interfaced_;

		//global size descriptors
		Integer n_active_elements_;
		Integer n_nodes_;



	//aux stuff	
	private:

		void add_and_index_nodes(
			const Mesh &serial_mesh, 
			const std::vector<Integer> &node_partitioning)
		{

			std::vector<Integer> visited(node_partitioning.size(), INVALID_INDEX);

			const auto n_nodes = node_partitioning.size();

			for(Integer i = 0; i < mesh_.n_elements(); ++i) {
				assert(mesh_.is_active(i));
				auto &e = mesh_.elem(i);
				const auto &parent_e = serial_mesh.elem(elem_map_.global(i));

				for(Integer k = 0; k < parent_e.nodes.size(); ++k) {
					const Integer n = parent_e.nodes[k];
					Integer local_n = INVALID_INDEX;

					if(visited[n] == INVALID_INDEX) {
						local_n = mesh_.add_point(serial_mesh.point(n));
						// global_node_id_.push_back(n);
						// node_partition_id_.push_back(node_partitioning[n]);
						// global_to_local_node_[n] = local_n;

						node_map_.add(
							local_n,
							n,
							node_partitioning[n]
						);
						
						visited[n] = local_n;
					} else {
						local_n = visited[n];
					}

					assert(e.nodes[k] == parent_e.nodes[k]);
					e.nodes[k] = local_n;
				}
			}

			for(Integer i = 0; i < mesh_.n_nodes(); ++i) {
				visited[node_map_.global(i)] = INVALID_INDEX;
			}

			mesh_.update_dual_graph();
		}
	};
}

#endif //PAR_MESH_HPP
