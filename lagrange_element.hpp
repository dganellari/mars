#ifndef LAGRANGE_ELEMENT_HPP
#define LAGRANGE_ELEMENT_HPP

#include <vector>
#include "vector.hpp"

namespace mars {
	template<class Element>
	class LagrangeElement {};

	template<Integer Dim, Integer ManifoldDim>
	class LagrangeElement<Simplex<Dim, ManifoldDim>> {
	public:

		inline void fun(
			const Vector<Real, Dim> &x_ref,
			std::vector<Real> &values) const
		{
			values.resize(ManifoldDim+1, 0.);
			values[0] = 1.;
			
			for(Integer i = 0; i < Dim; ++i) {
				values[0]    -= x_ref[i];
				values[i + 1] = x_ref[i];
			}
		}

		inline void grad(
			const Vector<Real, Dim> &x_ref,
			std::vector<Vector<Real, Dim>> &values) const
		{
			for(Integer i = 0; i < Dim; ++i) {
				values[0](i)   = -1;
				values[i+1]    = Vector<Real, Dim>().zero();
				values[i+1](i) = 1;
			}
		}
	};
}

#endif //LAGRANGE_ELEMENT_HPP
