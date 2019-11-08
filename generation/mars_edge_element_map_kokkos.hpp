#ifndef MARS_EDGE_ELEMENT_MAP_KOKKOS_HPP
#define MARS_EDGE_ELEMENT_MAP_KOKKOS_HPP

#include "mars_static_math_kokkos.hpp"

#include <map>
#include <vector>
#include <memory>
#include <cassert>
#include <set>

#include "mars_edge_kokkos.hpp"
#include "mars_edge_select_kokkos.hpp"
#include "mars_longest_edge.hpp"
#include "mars_fwd.hpp"

namespace mars {

	template<Integer N>
	class SubManifoldElementMap<N,KokkosImplementation> {
	public:

/*
		Based on Rivara the smallest angle remaining is 27.89 degrees
		which means 360 / 27.89 = 12.9. The max number of incidents is then 13.
		Adding here un upper limit considering that incative elements might be
		still be present in the list then it will be 2* 13= 26.
		The allocation is done using 32 as a grace alloc. for special cases
		or bad quality meshes (very small and very large angles).
*/
		using ElementVector = TempArray<Integer,32>;

		virtual ~SubManifoldElementMap() {}

		template<typename Elem>
		MARS_INLINE_FUNCTION
		static void update_elem(const UnorderedMap<Side<N,KokkosImplementation>,ElementVector>& mapping, const Elem &e)
		{
			using Comb = Combinations<Elem::ManifoldDim_ + 1, N, KokkosImplementation>;
			Comb combinations;

			TempArray<Integer, N> nodes;
			for (Integer i = 0;
					i<Comb::value;
					++i)
			{

				for (Integer j = 0; j < N; ++j)
				{
					nodes[j] = e.nodes[combinations.combs[i][j]];
				}

				auto result= mapping.insert(Side<N, KokkosImplementation>(nodes));

				if(!result.failed()){
					auto &vec = mapping.value_at(result.index());
					vec.insert(e.id);
				}else{
					 printf("Edge Element Map: Exceeded UnorderedMap capacity\n");
					//TODO: handle the case of failure insert. GO to the host for rehash.
				}

			}
		}

		template<Integer Dim, Integer ManifoldDim>
		struct UMapUpdate
		{
			using PMesh = Mesh<Dim, ManifoldDim, KokkosImplementation>;

			PMesh* mesh;
			UnorderedMap<Side<N,KokkosImplementation>,ElementVector> mapping;


			/*void init(const PMesh *m)
			{

				PMesh* tmp = (PMesh*) Kokkos::kokkos_malloc(sizeof(PMesh));

				Kokkos::parallel_for("CreateMeshObject", 1, KOKKOS_LAMBDA (const int&)
				{
					new ((PMesh*)tmp) PMesh(*m); //same as this.mesh if mesh instead of tmp and this is a host pointer.
				});

				mesh = tmp; //otherwise "this" pointer is still a host pointer within the parallel_for.
			}*/

			UMapUpdate(UnorderedMap<Side<N,KokkosImplementation>, ElementVector> mp, PMesh *ms) :
					mapping(mp), mesh(ms)
			{
				//init(ms);
			}

			MARS_INLINE_FUNCTION
			void operator()(int element_id) const
			{
				if(mesh->is_active(element_id))
				{
					update_elem(mapping, mesh->elem(element_id));
				}

			}
		};

		void reserve_map(Integer capacity)
		{
			mapping_ = UnorderedMap<Side<N,KokkosImplementation>,ElementVector>(capacity);
		}

		void rehash_map(Integer capacity)
		{
			mapping_.rehash(capacity);
		}

		template<Integer Dim, Integer ManifoldDim>
		void update(Mesh<Dim, ManifoldDim, KokkosImplementation> *mesh,
				const Integer nr_elements)
		{
			Kokkos::parallel_for(nr_elements,
					UMapUpdate<Dim, ManifoldDim>(mapping_, mesh));
		}

		/*template<Integer Dim, Integer ManifoldDim>
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
		}*/

		MARS_INLINE_FUNCTION
		const ElementVector &elements(const Side<N,KokkosImplementation> &side) const
		{
			static const ElementVector null_vec;

			uint32_t it = mapping_.find(side);
			if(it == mapping_.invalid_index) {
				return null_vec;
			}

			assert(mapping_.value_at(it));
			return mapping_.value_at(it);
		}


	/*	template<Integer Dim, Integer ManifoldDim>
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
		}*/




		struct UMapDescribe
		{
			UnorderedMap<Side<N,KokkosImplementation>,ElementVector> mapping;

			UMapDescribe(UnorderedMap<Side<N,KokkosImplementation>, ElementVector> m) :
					mapping(m)
			{
			}

			MARS_INLINE_FUNCTION
			void operator()(int element_id) const
			{
				int count = 0;
				printf("mapping.size(): %i ", mapping.capacity());
				for (int element_id = 0; element_id < mapping.capacity();
						++element_id)
				{
					if (mapping.valid_at(element_id))
					{
						auto key = mapping.key_at(element_id);
						auto value = mapping.value_at(element_id);

						for (Integer i = 0; i < N; ++i)
						{
							printf("%li ", key.nodes[i]);
						}
						printf("-> ");

						for (Integer i = 0; i < value.index; ++i)
						{
							printf(" %li ", value[i]);
						}
						printf("\n");
						++count;
					}

				}

				printf("count: %i \n", count);

			}
		};

		void describe(std::ostream &os) const
		{
			os << "-----------------------------------\n";
			os << "SubManifoldElementMap<" << N  << ">\n";

			Kokkos::parallel_for(1, UMapDescribe(mapping_));

			os << "-----------------------------------\n";
		}

		inline bool empty() const
		{
			return (mapping_.size() == 0);
		}

		void clear()
		{
			mapping_.clear();
		}

		UnorderedMap<Side<N,KokkosImplementation>,ElementVector> mapping_;
	};

	/*class EdgeElementMap : public SubManifoldElementMap<2> {
	public:

		using SubManifoldElementMap<2>::elements;

		const ElementVector &elements(const Integer a_node, const Integer another_node) const
		{
			return SubManifoldElementMap<2>::elements({a_node, another_node});
		}
	};*/

}

#endif //MARS_EDGE_ELEMENT_MAP_HPP
