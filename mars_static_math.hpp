#ifndef MARS_STATIC_MATH_HPP
#define MARS_STATIC_MATH_HPP

#include "mars_base.hpp" 
#include <iostream>
#include <vector>

namespace mars {



	constexpr Integer factorial(const Integer N)
	{
		assert(N>=-1 && " in factorial(N), N>=-1 ");
		switch(N)
		{
			case -1: return 1;
			case  0: return 1;
			case  1: return 1;
			default: return N*factorial(N-1);
		}
	}

	constexpr Integer binomial_coefficient(const Integer N,Integer K)
	{
		return factorial(N)/( factorial(K) * factorial(N-K));
	}



	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////// Number<N>
	//////// Class used to use Numbers as types and avoid constexpr arrays
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	template<Integer N>
	class Number
	{
	public:
		static constexpr Integer value=N;
	};
	using Zero=Number<0>;
	using One=Number<1>;


	template<Integer N>
	class Factorial {
	public:
		static const Integer value = Factorial<N-1>::value * N;
	};

	template<>
	class Factorial<1> {
	public:
		static const Integer value = 1;
	};

	template<>
	class Factorial<0> {
	public:
		static const Integer value = 1;
	};


	template<>
	class Factorial<-1> {
	public:
		static const Integer value = 1;
	};

	template<Integer Base, Integer Expon>
	class Power {
	public:
		static const Integer value = Base * Power<Base, Expon-1>::value;
	};

	template<Integer Base>
	class Power<Base, 1> {
	public:
		static const Integer value = Base;
	};

	template<Integer Base>
	class Power<Base, 0> {
	public:
		static const Integer value = 1;
	};

	template<Integer N, Integer K>
	class CombinationsAux {
	public:
		static const Integer NChooseK = Factorial<N>::value/(Factorial<K>::value * Factorial<N-K>::value);
		
		static void apply(std::array<std::array<Integer, K>, NChooseK> &combs)
		{
			std::array<Integer, K> data;
			Integer comb_index = 0;
			apply(data, 0, 0, combs, comb_index);
		}

	private:
		static void apply(
			std::array<Integer, K> &data,
			const Integer index, 
			const Integer i,
			std::array<std::array<Integer, K>, NChooseK> &combs,
			Integer &comb_index)
		{
			if(index == K) {
				std::copy(std::begin(data), std::end(data), std::begin(combs[comb_index++]));
				return;
			}

			if(i >= N) {
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
	class Combinations {
	public:
		static const Integer value = Factorial<N>::value/(Factorial<ChooseM>::value * Factorial<N-ChooseM>::value);
		std::array<std::array<Integer, ChooseM>, value> combs;

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
		static void choose(
			const Integer k,
			const std::array<T, N> &in,
			std::array<T, ChooseM> &out)
		{
			assert(k < value);
			const auto &comb = instance().combs[k];

			for(Integer i = 0; i < ChooseM; ++i) {
				assert(comb[i] < N);
				assert(comb[i] >= 0);

				out[i] = in[comb[i]];
			}
		}

		static void generate(const Integer k, Integer comb[ChooseM])
		{
			std::copy(instance().combs[k].begin(), instance().combs[k].end(), comb);
		}

		static constexpr std::array<Integer,ChooseM> generate(const Integer k)
		{
   //          std::array<Integer,ChooseM> comb;
			// // std::copy(instance().combs[k].begin(), instance().combs[k].end(), comb);
			// const auto combbs_tmp=instance().combs[k];
			// for(Integer kk=0;kk<combbs_tmp.size();kk++)
			// 	comb[kk]=combbs_tmp[kk];
			return instance().combs[k];
		}

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
	class Combinations<N, 1> {
	public:
		static const Integer value = N;
		static void generate(const Integer k, Integer comb[1])
		{
			comb[0] = k;
		}

		template<typename T>
		static void choose(
			const Integer k,
			const std::array<T, N> &in,
			std::array<T, 1> &out)
		{	

			assert(k < N);
			out[0] = in[k];
		}
	};


	template<Integer N,Integer M>
	class Modulo {
	public:
		static const Integer value = N - M* (N/M);
	};


	// given an index_sequence<std::size_t Is...>, create the corresponding variadic T...
	template <typename T, std::size_t>
	using getTypeSequence = T;

	template<typename T,typename S>
	class IsSame 
	{
	public:
	  static constexpr bool value=std::is_same<T,S>::value;
	};


	template<typename T,typename S>
	class IsDifferent
	{
	public:
	  static constexpr bool value=std::conditional<std::is_same<T,S>::value,Number<0>,Number<1>>::type ::value;
	};
	// template<>
	// class Combinations<3, 2> {
	// public:
	// 	static const Integer value = 3;
	// 	static void generate(const Integer k, Integer comb[2])
	// 	{
	// 		switch(k) {
	// 			case 0:
	// 			{
	// 				comb[0] = 0;
	// 				comb[1] = 1;
	// 				return;
	// 			}

	// 			case 1:
	// 			{
	// 				comb[0] = 0;
	// 				comb[1] = 2;
	// 				return;
	// 			}

	// 			case 2:
	// 			{
	// 				comb[0] = 1;
	// 				comb[1] = 2;
	// 				return;
	// 			}

	// 			default:
	// 			{
	// 				assert(false);
	// 				return;
	// 			}
	// 		}
	// 	}
	// };


	// template<>
	// class Combinations<4, 2> {
	// public:
	// 	static const Integer value = 6;
	// 	static void generate(const Integer k, Integer comb[2])
	// 	{
	// 		switch(k) {
	// 			case 0:
	// 			{
	// 				comb[0] = 0;
	// 				comb[1] = 1;
	// 				return;
	// 			}

	// 			case 1:
	// 			{
	// 				comb[0] = 0;
	// 				comb[1] = 2;
	// 				return;
	// 			}

	// 			case 2:
	// 			{
	// 				comb[0] = 0;
	// 				comb[1] = 3;
	// 				return;
	// 			}

	// 			case 3:
	// 			{
	// 				comb[0] = 1;
	// 				comb[1] = 2;
	// 				return;
	// 			}

	// 			case 4:
	// 			{
	// 				comb[0] = 1;
	// 				comb[1] = 3;
	// 				return;
	// 			}

	// 			case 5:
	// 			{
	// 				comb[0] = 2;
	// 				comb[1] = 3;
	// 				return;
	// 			}

	// 			default:
	// 			{
	// 				assert(false);
	// 				return;
	// 			}
	// 		}
	// 	}
	// };

	// template<>
	// class Combinations<4, 3> {
	// public:
	// 	static const Integer value = 4;
	// 	static void generate(const Integer k, Integer comb[3])
	// 	{
	// 		switch(k) {
	// 			case 0:
	// 			{
	// 				comb[0] = 0;
	// 				comb[1] = 1;
	// 				comb[2] = 2;
	// 				return;
	// 			}

	// 			case 1:
	// 			{
	// 				comb[0] = 0;
	// 				comb[1] = 1;
	// 				comb[2] = 3;
	// 				return;
	// 			}

	// 			case 2:
	// 			{
	// 				comb[0] = 0;
	// 				comb[1] = 2;
	// 				comb[2] = 3;
	// 				return;
	// 			}

	// 			case 3:
	// 			{
	// 				comb[0] = 1;
	// 				comb[1] = 2;
	// 				comb[2] = 3;
	// 				return;
	// 			}
	// 		}
	// 	}
	// };


	// template<>
	// class Combinations<5, 2> {
	// public:
	// 	static const Integer value = 10;
	// 	static void generate(const Integer k, Integer comb[2])
	// 	{
	// 		switch(k) {
	// 			case 0:
	// 			{
	// 				comb[0] = 0;
	// 				comb[1] = 1;
	// 				return;
	// 			}

	// 			case 1:
	// 			{
	// 				comb[0] = 0;
	// 				comb[1] = 2;
	// 				return;
	// 			}

	// 			case 2:
	// 			{
	// 				comb[0] = 0;
	// 				comb[1] = 3;
	// 				return;
	// 			}

	// 			case 3:
	// 			{
	// 				comb[0] = 0;
	// 				comb[1] = 4;
	// 				return;
	// 			}

	// 			case 4:
	// 			{
	// 				comb[0] = 1;
	// 				comb[1] = 2;
	// 				return;
	// 			}

	// 			case 5:
	// 			{
	// 				comb[0] = 1;
	// 				comb[1] = 3;
	// 				return;
	// 			}

	// 			case 6:
	// 			{
	// 				comb[0] = 1;
	// 				comb[1] = 4;
	// 				return;
	// 			}

	// 			case 7:
	// 			{
	// 				comb[0] = 2;
	// 				comb[1] = 3;
	// 				return;
	// 			}

	// 			case 8:
	// 			{
	// 				comb[0] = 2;
	// 				comb[1] = 4;
	// 				return;
	// 			}

	// 			case 9:
	// 			{
	// 				comb[0] = 3;
	// 				comb[1] = 4;
	// 				return;
	// 			}

	// 			default:
	// 			{
	// 				assert(false);
	// 				return;
	// 			}
	// 		}
	// 	}
	// };

	// template<>
	// class Combinations<5, 3> {
	// public:
	// 	static const Integer value = 10;
	// 	static void generate(const Integer k, Integer comb[3])
	// 	{
	// 		switch(k) {
	// 			case 0:
	// 			{
	// 				comb[0] = 0;
	// 				comb[1] = 1;
	// 				comb[2] = 2;
	// 				return;
	// 			}

	// 			case 1:
	// 			{
	// 				comb[0] = 0;
	// 				comb[1] = 1;
	// 				comb[2] = 3;
	// 				return;
	// 			}

	// 			case 2:
	// 			{
	// 				comb[0] = 0;
	// 				comb[1] = 1;
	// 				comb[2] = 4;
	// 				return;
	// 			}

	// 			case 3:
	// 			{
	// 				comb[0] = 0;
	// 				comb[1] = 2;
	// 				comb[2] = 3;
	// 				return;
	// 			}

	// 			case 4:
	// 			{
	// 				comb[0] = 0;
	// 				comb[1] = 2;
	// 				comb[1] = 4;
	// 				return;
	// 			}

	// 			case 5:
	// 			{
	// 				comb[0] = 0;
	// 				comb[1] = 3;
	// 				comb[2] = 4;
	// 				return;
	// 			}

	// 			case 6:
	// 			{
	// 				comb[0] = 1;
	// 				comb[1] = 2;
	// 				comb[2] = 3;
	// 				return;
	// 			}

	// 			case 7:
	// 			{
	// 				comb[0] = 1;
	// 				comb[1] = 2;
	// 				comb[2] = 4;
	// 				return;
	// 			}

	// 			case 8:
	// 			{
	// 				comb[0] = 1;
	// 				comb[1] = 3;
	// 				comb[2] = 4;
	// 				return;
	// 			}

	// 			case 9:
	// 			{
	// 				comb[0] = 2;
	// 				comb[1] = 3;
	// 				comb[2] = 4;
	// 				return;
	// 			}

	// 			default:
	// 			{
	// 				assert(false);
	// 				return;
	// 			}
	// 		}
	// 	}
	// };


	// template<>
	// class Combinations<5, 4> {
	// public:
	// 	static const Integer value = 5;
	// 	static void generate(const Integer k, Integer comb[4])
	// 	{
	// 		switch(k) {
	// 			case 0:
	// 			{
	// 				comb[0] = 0;
	// 				comb[1] = 1;
	// 				comb[2] = 2;
	// 				comb[3] = 3;
	// 				return;
	// 			}

	// 			case 1:
	// 			{
	// 				comb[0] = 0;
	// 				comb[1] = 1;
	// 				comb[2] = 2;
	// 				comb[3] = 4;
	// 				return;
	// 			}

	// 			case 2:
	// 			{
	// 				comb[0] = 0;
	// 				comb[1] = 1;
	// 				comb[2] = 3;
	// 				comb[3] = 4;
	// 				return;
	// 			}

	// 			case 3:
	// 			{
	// 				comb[0] = 0;
	// 				comb[1] = 2;
	// 				comb[2] = 3;
	// 				comb[3] = 4;
	// 				return;
	// 			}

	// 			case 4:
	// 			{
	// 				comb[0] = 1;
	// 				comb[1] = 2;
	// 				comb[2] = 3;
	// 				comb[3] = 4;
	// 				return;
	// 			}

	// 			default:
	// 			{
	// 				assert(false);
	// 				return;
	// 			}
	// 		}
	// 	}
	// };


}

#endif //MARS_STATIC_MATH_HPP
