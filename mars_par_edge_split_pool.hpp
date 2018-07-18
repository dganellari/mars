#ifndef MARS_PAR_EDGE_SPLIT_POOL_HPP
#define MARS_PAR_EDGE_SPLIT_POOL_HPP

#include "mars_communicator.hpp"
#include <sstream>
#include <string>
#include <iomanip>

namespace mars {
	template<class Mesh>
	class Bisection;

	template<Integer Dim, Integer ManifoldDim>
	class ParMesh;

	class ParEdgeSplitPool {
	public:
		using OutputStream = std::ostringstream;
		using InputStream  = std::istringstream;
		using BufferObject = std::string;

		template<class ParMesh>
		void build_edge_interface(
			ParMesh &mesh,
			Bisection<typename ParMesh::Mesh> &bisection)
		{
			using SideElem = typename ParMesh::SideElem;
		
			std::vector<OutputStream> output(comm_.size());
			std::vector<BufferObject> send_buff(comm_.size());
			std::vector<BufferObject> recv_buff(comm_.size());

			std::vector<std::vector<SideElem>> sides(comm_.size());

			std::vector<SideElem> local_sides;
			for(Integer r = 0; r < comm_.size(); ++r) {
				if(mesh.is_interfaced(r)) {
					mesh.collect_interface_sides(
						r,
						local_sides,
						true,
						true
					);

					Integer n_sides = local_sides.size();
					output[r] << n_sides;

					for(auto &s : local_sides) {
						write(s, output[r]);
					}
				}
			}

			//TODO write additional info

			//create output buffers
			for(Integer r = 0; r < comm_.size(); ++r) {
				if(mesh.is_interfaced(r)) {
					//FIXME
					send_buff[r] = output[r].str();
					comm_.i_send(&send_buff[r][0], send_buff[r].size(), r, r);
				}
			}

			Integer n_interfaced = mesh.n_interfaced();

			for(Integer i = 0; i < n_interfaced; ++i) {
				Integer rank, size;
				while ( !comm_.i_probe_any<byte>( &rank, &size ) ) {}
				recv_buff[rank].resize(size);
				comm_.i_recv(&recv_buff[rank][0], size, rank, comm_.rank());
			}

			for(Integer i = 0; i < n_interfaced; ++i) {
				Integer rank, index;
				while ( !comm_.test_recv_any( &rank, &index ) ) {}

				Integer n_sides = 0;
				InputStream is(recv_buff[rank]);
				is >> n_sides;

				sides[rank].resize(n_sides);

				for(Integer k = 0; k < n_sides; ++k) {
					read(sides[rank][k], is);
				}
			}

			//TODO read additional info

			//wait that all sends have been accomplished
			comm_.wait_all();
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
