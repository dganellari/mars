#ifndef MARS_STATIC_MATH_HPP
#define MARS_STATIC_MATH_HPP

#include "base.hpp"
#include <iostream>

namespace mars {
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

	template<Integer N, Integer ChooseM>
	class Combinations {
	public:
		static const Integer value = Factorial<N>::value/(Factorial<ChooseM>::value * Factorial<N-ChooseM>::value);


		//https://www4.uwsp.edu/math/nwodarz/Math209Files/209-0809F-L10-Section06_03-AlgorithmsForGeneratingPermutationsAndCombinations-Notes.pdf
		// static void generate(const Integer k, Integer comb[ChooseM])
		// {
		// 	assert(k < value);

		// 	for(Integer i = 0; i < ChooseM; ++i) {
		// 		comb[i] = i;
		// 	}

		// 	if(k == 0) {
		// 		return;
		// 	}

		// 	for(Integer i = 1; i < value; ++i) {
		// 		Integer m       = ChooseM - 1;
		// 		Integer max_val = value   - 1;

		// 		while(comb[m] == max_val) {
		// 			--m;
		// 			--max_val;
		// 		}

		// 		++comb[ChooseM-1];

		// 		for(Integer j = m + 1; j < ChooseM; ++j) {
		// 			comb[j] = comb[j-1] + 1;
		// 		}

		// 		if(k == i) {
		// 			return;
		// 		}
		// 	}
		// }

		// static void print_all(std::ostream &os = std::cout) 
		// {
		// 	Integer comb[ChooseM];
		// 	for(Integer i = 0; i < value; ++i) {
		// 		generate(i, comb);

		// 		for(Integer k = 0; k < ChooseM; ++k) {
		// 			os << comb[k] << " ";
		// 		}

		// 		os << std::endl;
		// 	}
		// }

	};

	template<Integer N>
	class Combinations<N, 1> {
	public:
		static const Integer value = N;
		static void generate(const Integer k, Integer comb[1])
		{
			comb[0] = k;
		}
	};

	template<>
	class Combinations<3, 2> {
	public:
		static const Integer value = 3;
		static void generate(const Integer k, Integer comb[2])
		{
			switch(k) {
				case 0:
				{
					comb[0] = 0;
					comb[1] = 1;
					return;
				}

				case 1:
				{
					comb[0] = 0;
					comb[1] = 2;
					return;
				}

				case 2:
				{
					comb[0] = 1;
					comb[1] = 2;
					return;
				}

				default:
				{
					assert(false);
					return;
				}
			}
		}
	};


	template<>
	class Combinations<4, 2> {
	public:
		static const Integer value = 6;
		static void generate(const Integer k, Integer comb[2])
		{
			switch(k) {
				case 0:
				{
					comb[0] = 0;
					comb[1] = 1;
					return;
				}

				case 1:
				{
					comb[0] = 0;
					comb[1] = 2;
					return;
				}

				case 2:
				{
					comb[0] = 0;
					comb[1] = 3;
					return;
				}

				case 3:
				{
					comb[0] = 1;
					comb[1] = 2;
					return;
				}

				case 4:
				{
					comb[0] = 1;
					comb[1] = 3;
					return;
				}

				case 5:
				{
					comb[0] = 2;
					comb[1] = 3;
					return;
				}

				default:
				{
					assert(false);
					return;
				}
			}
		}
	};

	template<>
	class Combinations<4, 3> {
	public:
		static const Integer value = 4;
		static void generate(const Integer k, Integer comb[3])
		{
			switch(k) {
				case 0:
				{
					comb[0] = 0;
					comb[1] = 1;
					comb[2] = 2;
					return;
				}

				case 1:
				{
					comb[0] = 0;
					comb[1] = 1;
					comb[2] = 3;
					return;
				}

				case 2:
				{
					comb[0] = 0;
					comb[1] = 2;
					comb[2] = 3;
					return;
				}

				case 3:
				{
					comb[0] = 1;
					comb[1] = 2;
					comb[2] = 3;
					return;
				}
			}
		}
	};


	template<>
	class Combinations<5, 2> {
	public:
		static const Integer value = 10;
		static void generate(const Integer k, Integer comb[2])
		{
			switch(k) {
				case 0:
				{
					comb[0] = 0;
					comb[1] = 1;
					return;
				}

				case 1:
				{
					comb[0] = 0;
					comb[1] = 2;
					return;
				}

				case 2:
				{
					comb[0] = 0;
					comb[1] = 3;
					return;
				}

				case 3:
				{
					comb[0] = 0;
					comb[1] = 4;
					return;
				}

				case 4:
				{
					comb[0] = 1;
					comb[1] = 2;
					return;
				}

				case 5:
				{
					comb[0] = 1;
					comb[1] = 3;
					return;
				}

				case 6:
				{
					comb[0] = 1;
					comb[1] = 4;
					return;
				}

				case 7:
				{
					comb[0] = 2;
					comb[1] = 3;
					return;
				}

				case 8:
				{
					comb[0] = 2;
					comb[1] = 4;
					return;
				}

				case 9:
				{
					comb[0] = 3;
					comb[1] = 4;
					return;
				}

				default:
				{
					assert(false);
					return;
				}
			}
		}
	};

	template<>
	class Combinations<5, 3> {
	public:
		static const Integer value = 10;
		static void generate(const Integer k, Integer comb[3])
		{
			switch(k) {
				case 0:
				{
					comb[0] = 0;
					comb[1] = 1;
					comb[2] = 2;
					return;
				}

				case 1:
				{
					comb[0] = 0;
					comb[1] = 1;
					comb[2] = 3;
					return;
				}

				case 2:
				{
					comb[0] = 0;
					comb[1] = 1;
					comb[2] = 4;
					return;
				}

				case 3:
				{
					comb[0] = 0;
					comb[1] = 2;
					comb[2] = 3;
					return;
				}

				case 4:
				{
					comb[0] = 0;
					comb[1] = 2;
					comb[1] = 4;
					return;
				}

				case 5:
				{
					comb[0] = 0;
					comb[1] = 3;
					comb[2] = 4;
					return;
				}

				case 6:
				{
					comb[0] = 1;
					comb[1] = 2;
					comb[2] = 3;
					return;
				}

				case 7:
				{
					comb[0] = 1;
					comb[1] = 2;
					comb[2] = 4;
					return;
				}

				case 8:
				{
					comb[0] = 1;
					comb[1] = 3;
					comb[2] = 4;
					return;
				}

				case 9:
				{
					comb[0] = 2;
					comb[1] = 3;
					comb[2] = 4;
					return;
				}

				default:
				{
					assert(false);
					return;
				}
			}
		}
	};


	template<>
	class Combinations<5, 4> {
	public:
		static const Integer value = 5;
		static void generate(const Integer k, Integer comb[4])
		{
			switch(k) {
				case 0:
				{
					comb[0] = 0;
					comb[1] = 1;
					comb[2] = 2;
					comb[3] = 3;
					return;
				}

				case 1:
				{
					comb[0] = 0;
					comb[1] = 1;
					comb[2] = 2;
					comb[3] = 4;
					return;
				}

				case 2:
				{
					comb[0] = 0;
					comb[1] = 1;
					comb[2] = 3;
					comb[3] = 4;
					return;
				}

				case 3:
				{
					comb[0] = 0;
					comb[1] = 2;
					comb[2] = 3;
					comb[3] = 4;
					return;
				}

				case 4:
				{
					comb[0] = 1;
					comb[1] = 2;
					comb[2] = 3;
					comb[3] = 4;
					return;
				}

				default:
				{
					assert(false);
					return;
				}
			}
		}
	};

}

#endif //MARS_STATIC_MATH_HPP
