#ifndef MARS_STATIC_MATH_KOKKOS_HPP
#define MARS_STATIC_MATH_KOKKOS_HPP

#include "mars_base.hpp"
#include <iostream>
#include <vector>
#include "mars_static_math.hpp"
#include "mars_utils_kokkos.hpp"

namespace mars {

	template<Integer N, Integer K>
	class CombinationsAux<N,K,KokkosImplementation> {
	public:
		static const Integer NChooseK = Factorial<N>::value/(Factorial<K>::value * Factorial<N-K>::value);
		
		static void apply(ViewMatrixTextureC<Integer,NChooseK,K> combs)
		{
			//std::array<Integer, K> data;
			Integer comb_index = 0;
//			apply(data, 0, 0, combs, comb_index);
			apply(0, 0, combs, comb_index);

		}

	private:
		static void apply(
			//Integer* data, //call this with the shared memory view shared pointer
			const Integer index, 
			const Integer i,
			ViewMatrixTextureC<Integer,NChooseK,K> combs,
			Integer &comb_index)
		{
			if(index == K) {
				//std::copy(std::begin(data), std::end(data), std::begin(combs(comb_index++)));
				return;
			}

			if(i >= N) {
				return;
			}

			//data[index] = i;
			combs(comb_index++)[index] = i;

			apply(index+1, i+1, combs, comb_index);
			
			// current is excluded, replace it with next (Note that
			// i+1 is passed, but index is not changed)
			apply(index, i+1, combs, comb_index);
		}
	};

	template<Integer N, Integer ChooseM>
	class Combinations<N,ChooseM, KokkosImplementation> {
	public:
		static const Integer value = Factorial<N>::value/(Factorial<ChooseM>::value * Factorial<N-ChooseM>::value);
		/*std::array<std::array<Integer, ChooseM>, value> combs;*/
		ViewMatrixTextureC<Integer,value,ChooseM> combs;

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
		static void choose(const Integer k, const T* in, T* out) {
			assert(k < value);
			const auto &comb = instance().combs(k);

			for (Integer i = 0; i < ChooseM; ++i) {
				assert(comb(i) < N);
				assert(comb(i) >= 0);

				out(i) = in(comb(i));
			}
		}

		/*static void generate(const Integer k, Integer comb[ChooseM])
		{
			std::copy(instance().combs[k].begin(), instance().combs[k].end(), comb);
		}*/

	private:
		Combinations()
		{
			CombinationsAux<N, ChooseM>::apply(combs);
		}

		inline static const Combinations &instance()
		{	
			static const Combinations instance_;
			return instance_;
		}
	};

	template<Integer N>
	class Combinations<N, 1,KokkosImplementation> {
	public:
		static const Integer value = N;
		/*static void generate(const Integer k, Integer comb[1])
		{
			comb[0] = k;
		}*/

		template<typename T>
		static void choose(const Integer k, const T* in, T* out) {

			assert(k < N);
			out(0) = in(k);
		}
	};

}

#endif //MARS_STATIC_MATH_HPP
