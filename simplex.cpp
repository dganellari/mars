#include "simplex.hpp"

namespace mars {
	template class Simplex<1, 0>;
	template class Simplex<1, 1>;

	template class Simplex<2, 0>;
	template class Simplex<2, 1>;
	template class Simplex<2, 2>;

	template class Simplex<3, 0>;
	template class Simplex<3, 1>;
	template class Simplex<3, 2>;
	template class Simplex<3, 3>;

	template class Simplex<4, 0>;
	template class Simplex<4, 1>;
	template class Simplex<4, 2>;
	template class Simplex<4, 3>;
	template class Simplex<4, 4>;
}