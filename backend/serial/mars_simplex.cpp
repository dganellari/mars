#include "mars_simplex.hpp"

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

    template class Simplex<5, 0>;
    template class Simplex<5, 1>;
    template class Simplex<5, 2>;
    template class Simplex<5, 3>;
    template class Simplex<5, 4>;
    template class Simplex<5, 5>;

    template class Simplex<6, 0>;
    template class Simplex<6, 1>;
    template class Simplex<6, 2>;
    template class Simplex<6, 3>;
    template class Simplex<6, 4>;
    template class Simplex<6, 5>;
    template class Simplex<6, 6>;
}  // namespace mars
