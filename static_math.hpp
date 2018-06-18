#ifndef MARS_STATIC_MATH_HPP
#define MARS_STATIC_MATH_HPP

#include "base.hpp"

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
	    static const Integer value = 0;
	};
}

#endif //MARS_STATIC_MATH_HPP
