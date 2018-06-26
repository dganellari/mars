#ifndef MARS_EDGE_ELEMENT_MAP_HPP
#define MARS_EDGE_ELEMENT_MAP_HPP

#include "edge.hpp"

#include <map>
#include <vector>
#include <memory>
#include <cassert>
#include <set>

namespace mars {
	template<Integer Dim, Integer ManifoldDim>
	class Mesh;

	template<Integer N>
	class SubManifoldElementMap {
	public:
		using ElementVector = std::vector<Integer>;

		virtual ~SubManifoldElementMap() {}

		template<Integer Dim, Integer ManifoldDim>
		void update(const Simplex<Dim, ManifoldDim> &e)
		{
			std::array<Integer, N> s;
			std::array<Integer, N> nodes;
			for(Integer i = 0; i < Combinations<ManifoldDim+1, N>::value; ++i) {
				Combinations<ManifoldDim+1, N>::generate(i, &s[0]);

				for(Integer j = 0; j < N; ++j) {
					nodes[j] = e.nodes[s[j]];
				}

				auto &vec = mapping_[Side<N>(nodes)];
				if(!vec) {
					vec = std::make_shared<ElementVector>();
				}

				vec->push_back(e.id);
			}
		}

		template<Integer Dim, Integer ManifoldDim>
		void update(const Mesh<Dim, ManifoldDim> &mesh)
		{
			auto ne = mesh.n_elements();
			for(Integer i = 0; i < ne; ++i) {
				// if(mesh.is_active(i)) {
					update(mesh.elem(i));
				// }
			}
		}

		template<Integer Dim, Integer ManifoldDim>
		void update_active(const Mesh<Dim, ManifoldDim> &mesh)
		{
			auto ne = mesh.n_elements();
			for(Integer i = 0; i < ne; ++i) {
				if(mesh.is_active(i)) {
					update(mesh.elem(i));
				}
			}
		}

		template<Integer Dim, Integer ManifoldDim>
		void build(const Mesh<Dim, ManifoldDim> &mesh)
		{	
			mapping_.clear();
			
			auto ne = mesh.n_elements();
			for(Integer i = 0; i < ne; ++i) {
				update(mesh.elem(i));
			}
		}

		const ElementVector &elements(const Side<N> &side) const
		{
			static const ElementVector null_vec;

			auto it = mapping_.find(side);
			if(it == mapping_.end()) {
				return null_vec;
			}

			assert(it->second);
			return *it->second;
		}

		template<Integer Dim, Integer ManifoldDim>
		void adj(const Simplex<Dim, ManifoldDim> &e, std::vector<Integer> &res) const
		{
			res.clear();
			std::set<Integer> temp;

			std::array<Integer, N> s;
			std::array<Integer, N> nodes;
			for(Integer i = 0; i < Combinations<ManifoldDim+1, N>::value; ++i) {
				Combinations<ManifoldDim+1, N>::generate(i, &s[0]);

				for(Integer j = 0; j < N; ++j) {
					nodes[j] = e.nodes[s[j]];
				}

				const ElementVector &els = elements(Side<N>(nodes));
					
				for(auto el : els) {
					temp.insert(el);
				}
			}

			res.insert(res.end(), temp.begin(), temp.end());
		}

		void describe(std::ostream &os) const
		{
			os << "-----------------------------------\n";
			os << "SubManifoldElementMap<" << N  << ">\n";

			for(const auto &m : mapping_) {
				for(const auto n : m.first.nodes) {
					os << n << " ";
				}

				os << "-> ";

				for(const auto n : *m.second) {
					os << n << " ";
				}

				os << "\n";
			}

			os << "-----------------------------------\n";
		}

		std::map<Side<N>, std::shared_ptr<ElementVector> > mapping_;
	};

	class EdgeElementMap : public SubManifoldElementMap<2> {
	public:

		using SubManifoldElementMap<2>::elements;

		const ElementVector &elements(const Integer a_node, const Integer another_node) const
		{
			return SubManifoldElementMap<2>::elements({a_node, another_node});
		}
	};

	template<Integer N, Integer LowestLevel = 1>
	class MultilevelElementMap {
	public:
		using ElementVector = std::vector<Integer>;

		template<Integer Dim, Integer ManifoldDim>
		void update(const Simplex<Dim, ManifoldDim> &e)
		{
			map_.update(e);
			sub_elem_map_.update(e);
		}

		template<Integer Dim, Integer ManifoldDim>
		void update(const Mesh<Dim, ManifoldDim> &m)
		{
			map_.update(m);
			sub_elem_map_.update(m);
		}

		template<Integer M>
		const ElementVector &elements(const Side<M> &side) const
		{
			static_assert(M < N, "M < N");
			return sub_elem_map_.elements(side);
		}

		const ElementVector &elements(const Side<N> &side) const
		{
			return map_.elements(side); 
		}

		const ElementVector &elements(const std::vector<Integer> &side) const
		{
			assert(side.size() <= N);
			if(side.size() == N) {
				return elements(Side<N>(side));
			} else {
				return sub_elem_map_.elements(side);
			}
		}

		template<Integer SubManifoldDim, Integer Dim, Integer ManifoldDim>
		void adj(
			const Simplex<Dim, ManifoldDim> &e,
			std::vector<Integer> &res) const
		{
			if(SubManifoldDim == N) {
				return map_.adj(e, res);
			} else {
				return sub_elem_map_.template adj<SubManifoldDim>(e, res);
			}
		}

		template<Integer Dim, Integer ManifoldDim>
		void adj(
			const Integer sub_manifold_dim,
			const Simplex<Dim, ManifoldDim> &e,
			std::vector<Integer> &res) const
		{
			if(sub_manifold_dim == N) {
				return map_.adj(e, res);
			} else {
				return sub_elem_map_.adj(sub_manifold_dim, e, res);
			}
		}

		void describe(std::ostream &os) const
		{
			map_.describe(os);
			sub_elem_map_.describe(os);
		}

	private:
		SubManifoldElementMap<N>  map_;
		MultilevelElementMap<N-1, LowestLevel> sub_elem_map_;
	};

	template<Integer N>
	class MultilevelElementMap<N, N> {
	public:
		using ElementVector = std::vector<Integer>;

		template<Integer Dim, Integer ManifoldDim>
		void update(const Simplex<Dim, ManifoldDim> &e)
		{
			map_.update(e);
		}

		template<Integer Dim, Integer ManifoldDim>
		void update(const Mesh<Dim, ManifoldDim> &m)
		{
			map_.update(m);
		}

		const ElementVector &elements(const Side<N> &side) const
		{
			return map_.elements(side); 
		}

		const ElementVector &elements(const std::vector<Integer> &side) const
		{
			return elements(Side<N>(side));
		}

		void describe(std::ostream &os) const
		{
			map_.describe(os);
		}

		template<Integer SubManifoldDim, Integer Dim, Integer ManifoldDim>
		void adj(
			const Simplex<Dim, ManifoldDim> &e,
			std::vector<Integer> &res) const
		{
			assert(SubManifoldDim == N);
			return map_.adj(e, res);
		}

		template<Integer Dim, Integer ManifoldDim>
		void adj(
			const Integer sub_manifold_dim,
			const Simplex<Dim, ManifoldDim> &e,
			std::vector<Integer> &res) const
		{
			assert(sub_manifold_dim == N);
			return map_.adj(e, res);
		}

	private:
		SubManifoldElementMap<N>  map_;
	};

}

#endif //MARS_EDGE_ELEMENT_MAP_HPP
