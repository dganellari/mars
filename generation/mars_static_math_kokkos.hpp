#ifndef MARS_STATIC_MATH_KOKKOS_HPP
#define MARS_STATIC_MATH_KOKKOS_HPP

#include "mars_base.hpp"
#include <array>
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

		static MARS_INLINE_FUNCTION void apply(ViewMatrixTextureC<Integer, NChooseK, K> combs)
		{
			TempArray<Integer, K> data;
			//std::cout << "Applying to combs" << std::endl;
			Integer comb_index = 0;
			apply(data, 0, 0, combs, comb_index);
		}

	private:
		static MARS_INLINE_FUNCTION void apply(
				TempArray<Integer, K>& data, //call this with the shared memory view shared pointer
				const Integer index,
				const Integer i,
				ViewMatrixTextureC<Integer, NChooseK, K> combs,
				Integer &comb_index)
		{
			if(index == K)
			{
				//std::copy(data, data+K, combs[comb_index++]);
				for(Integer j=0; j<K; ++j)
				{
					combs(comb_index,j) = data[j];
				}
				//std::cout << "Applying to combs with idx=" << comb_index << std::endl;

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
		static constexpr Integer value = Factorial<N>::value / (Factorial<ChooseM>::value * Factorial<N - ChooseM>::value);
		ViewMatrixTextureC<Integer, value, ChooseM>* combs;

		struct GenerateComb
		{
			ViewMatrixTextureC<Integer, value, ChooseM> combs;

			GenerateComb(ViewMatrixTextureC<Integer, value, ChooseM> cmb) :
					combs(cmb)
			{
			}

			MARS_INLINE_FUNCTION
			void operator()(int i) const
			{
				CombinationsAux<N, ChooseM, KokkosImplementation>::apply(combs);
			}
		};


		MARS_INLINE_FUNCTION Combinations()
		  : combs(nullptr)
		{}

		MARS_INLINE_FUNCTION
		void initialize() {
			combs = new ViewMatrixTextureC<Integer, value, ChooseM>("combinations");
			Kokkos::parallel_for(1, GenerateComb(*combs));
		}

		template<typename T>
		MARS_INLINE_FUNCTION void choose(const Integer k, const SubView<Integer,N>& in, T* out) const
		{

			assert(k < value);

			for(Integer i = 0; i < ChooseM; ++i)
			{
				assert(combs->operator()(k,i) < N);
				assert(combs->operator()(k,i) >= 0);

				out[i] = in[combs->operator()(k,i)];
			}
		}

	};

}

#endif //MARS_STATIC_MATH_HPP
