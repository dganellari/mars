#ifndef MARS_STATIC_MATH_HPP
#define MARS_STATIC_MATH_HPP

#include <iostream>
#include <vector>
#include "mars_base.hpp"
#include "mars_fwd.hpp"

namespace mars {
    template <Integer N>
    class Factorial {
    public:
        static const Integer value = Factorial<N - 1>::value * N;
    };

    template <>
    class Factorial<1> {
    public:
        static const Integer value = 1;
    };

    template <>
    class Factorial<0> {
    public:
        static const Integer value = 1;
    };

    template <>
    class Factorial<-1> {
    public:
        static const Integer value = 1;
    };

    template <Integer Base, Integer Expon>
    class Power {
    public:
        static const Integer value = Base * Power<Base, Expon - 1>::value;
    };

    template <Integer Base>
    class Power<Base, 1> {
    public:
        static const Integer value = Base;
    };

    template <Integer Base>
    class Power<Base, 0> {
    public:
        static const Integer value = 1;
    };

    template <Integer N, Integer K, class Implementation_>
    class CombinationsAux {
    public:
        static const Integer NChooseK = Factorial<N>::value / (Factorial<K>::value * Factorial<N - K>::value);

        static void apply(std::array<std::array<Integer, K>, NChooseK> &combs) {
            std::array<Integer, K> data;
            Integer comb_index = 0;
            apply(data, 0, 0, combs, comb_index);
        }

    private:
        static void apply(std::array<Integer, K> &data,
                          const Integer index,
                          const Integer i,
                          std::array<std::array<Integer, K>, NChooseK> &combs,
                          Integer &comb_index) {
            if (index == K) {
                std::copy(std::begin(data), std::end(data), std::begin(combs[comb_index++]));
                return;
            }

            if (i >= N) {
                return;
            }

            data[index] = i;

            apply(data, index + 1, i + 1, combs, comb_index);

            // current is excluded, replace it with next (Note that
            // i+1 is passed, but index is not changed)
            apply(data, index, i + 1, combs, comb_index);
        }
    };

    template <Integer N, Integer ChooseM, class Implementation_>
    class Combinations {
    public:
        static const Integer value = Factorial<N>::value / (Factorial<ChooseM>::value * Factorial<N - ChooseM>::value);
        std::array<std::array<Integer, ChooseM>, value> combs;

        static void print_all() {
            for (auto const &c : instance().combs) {
                for (auto n : c) {
                    std::cout << n << " ";
                }

                std::cout << std::endl;
            }
        }

        template <typename T>
        static void choose(const Integer k, const std::array<T, N> &in, std::array<T, ChooseM> &out) {
            assert(k < value);
            const auto &comb = instance().combs[k];

            for (Integer i = 0; i < ChooseM; ++i) {
                assert(comb[i] < N);
                assert(comb[i] >= 0);

                out[i] = in[comb[i]];
            }
        }

        static void generate(const Integer k, Integer comb[ChooseM]) {
            std::copy(instance().combs[k].begin(), instance().combs[k].end(), comb);
        }

    private:
        Combinations() { CombinationsAux<N, ChooseM>::apply(combs); }

        inline static const Combinations &instance() {
            static const Combinations instance_;
            return instance_;
        }
    };

    template <Integer N>
    class Combinations<N, 1, class Implementation_> {
    public:
        static const Integer value = N;
        static void generate(const Integer k, Integer comb[1]) { comb[0] = k; }

        template <typename T>
        static void choose(const Integer k, const std::array<T, N> &in, std::array<T, 1> &out) {
            assert(k < N);
            out[0] = in[k];
        }
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

}  // namespace mars

#endif  // MARS_STATIC_MATH_HPP
