#ifndef MARS_STATIC_MATH_KOKKOS_HPP
#define MARS_STATIC_MATH_KOKKOS_HPP

#include "mars_base.hpp"
#include <iostream>
#include <vector>
#include "mars_static_math.hpp"
#include "mars_sub_view.hpp"
#include "mars_utils_kokkos.hpp"

namespace mars {

	template<Integer N, Integer K>
	class CombinationsAux<N, K, KokkosImplementation>
	{
	public:
		static constexpr Integer NChooseK = Factorial<N>::value / (Factorial<K>::value * Factorial<N - K>::value);

		static MARS_INLINE_FUNCTION void apply(Integer combs[NChooseK][K])
		{
			TempArray<Integer, K> data;

			Integer comb_index = 0;
			apply(data, 0, 0, combs, comb_index);
		}

	private:
		static MARS_INLINE_FUNCTION void apply(
				TempArray<Integer, K>& data, //call this with the shared memory view shared pointer
				const Integer index,
				const Integer i,
				Integer combs[NChooseK][K],
				Integer &comb_index)
		{
			if(index == K)
			{
				//std::copy(data, data+K, combs[comb_index++]);
				for(Integer j=0; j<K; ++j)
				{
					combs[comb_index][j] = data[j];
				}

				++comb_index;
				return;
			}

			if(i >= N)
			{
				return;
			}

			data[index] = i;

			apply(data, index+1, i+1, combs, comb_index);

			// current is excluded, replace it with next (Note that
			// i+1 is passed, but index is not changed)
			apply(data, index, i+1, combs, comb_index);
		}
	};


	template<Integer N, Integer ChooseM>
	class Combinations<N, ChooseM, KokkosImplementation>
	{

	public:
		static const Integer value = Factorial<N>::value / (Factorial<ChooseM>::value * Factorial<N - ChooseM>::value);
		Integer combs[value][ChooseM];

		MARS_INLINE_FUNCTION Combinations()
		{
			CombinationsAux<N, ChooseM, KokkosImplementation>::apply(combs);
		}

		template<typename T>
		MARS_INLINE_FUNCTION void choose(const Integer k, const SubView<Integer,N>& in, T* out)
		{

			assert(k < value);

			const auto &comb = combs[k];

			for(Integer i = 0; i < ChooseM; ++i)
			{
				assert(comb[i] < N);
				assert(comb[i] >= 0);

				out[i] = in[comb[i]];
			}
		}

	};

	/*template<Integer N, Integer ChooseM>
	class Combinations<N,ChooseM, KokkosImplementation> {
	public:
		static const Integer value = Factorial<N>::value/(Factorial<ChooseM>::value * Factorial<N-ChooseM>::value);
		Integer combs[value][ChooseM];

		static void print_all()
		{
			for(auto const &c : instance().combs) {
				for(auto n : c) {
					std::cout << n << " ";
				}

				std::cout << std::endl;
			}
		}

		template<typename T>
		 MARS_INLINE_FUNCTION static void choose(const Integer k, const SubView<Integer,N>& in, T* out) {
			assert(k < value);

			const Combinations* a= instance();
		//	Integer *comb = a.combs[k];

			for (Integer i = 0; i < ChooseM; ++i) {
				assert(a->combs[k][i] < N);
				assert(a->combs[k][i] >= 0);
				out[i] = in[a->combs[k][i]];
			//	std::cout<<"p: "<<comb[i]<<" - "<<in[comb[i]]<<std::endl;
			}


		}

		static void generate(const Integer k, Integer comb[ChooseM])
		{
			std::copy(instance().combs[k].begin(), instance().combs[k].end(), comb);
		}

	private:
		MARS_INLINE_FUNCTION Combinations()
		{
			CombinationsAux<N, ChooseM, KokkosImplementation>::apply(combs);
		}

		MARS_INLINE_FUNCTION  static const Combinations* instance()
		{	
			static const Combinations* instance_ = new Combinations();
			return instance_;
		}
	};*/

/*	template<Integer N>
	class Combinations<N, 1,KokkosImplementation> {
	public:
		static const Integer value = N;
		static void generate(const Integer k, Integer comb[1])
		{
			comb[0] = k;
		}

		template<typename T>
		static MARS_INLINE_FUNCTION void choose(const Integer k, const T* in, T* out) {

			assert(k < N);
			out(0) = in(k);
		}
	};*/

}

#endif //MARS_STATIC_MATH_HPP
