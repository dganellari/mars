#ifndef GENERATION_MARS_UNIT_MESH_GENERATION_HPP_
#define GENERATION_MARS_UNIT_MESH_GENERATION_HPP_

#include "mars_mesh.hpp"

using namespace std;

namespace mars {
namespace unit_generation {

template<Integer Dim, Integer ManifoldDim>
bool generate_cube(Mesh<Dim, ManifoldDim>& mesh, const Integer xDim,
		const Integer yDim, const Integer zDim) {

	assert(ManifoldDim <= Dim);
	assert(xDim != 0);
	assert(yDim != 0);
	assert(zDim != 0);

	switch (ManifoldDim) {

	case 2: {
		const int n_elements = 2 * xDim * yDim;
		const int n_nodes = (xDim + 1) * (yDim + 1);
		mesh.reserve(n_elements, n_nodes);

		for (Integer i = 0; i <= xDim; ++i) {
			for (Integer j = 0; j <= yDim; ++j) {
				Vector<Real, Dim> p(
						{ static_cast<Real>(i) / static_cast<Real>(xDim),
								static_cast<Real>(j) / static_cast<Real>(yDim),
								0.0 });

				mesh.add_point(p);
			}
		}

		const int offset = yDim + 1;

		for (Integer i = 0; i < xDim; ++i) {
			for (Integer j = 0; j < yDim; ++j) {

				std::array<Integer, ManifoldDim + 1> nodes;

				nodes[0] = i * offset + j;
				nodes[1] = (i + 1) * offset + j;
				nodes[2] = (i + 1) * offset + (j + 1); //just to write it more clear

				mesh.add_elem(nodes);

				nodes[0] = i * offset + j;
				nodes[1] = (i + 1) * offset + (j + 1);
				nodes[2] = i * offset + (j + 1);

				mesh.add_elem(nodes);

			}
		}

		return true;
	}
	case 3: {

		const int n_elements = xDim * yDim * zDim;
		const int n_nodes = (2 * xDim + 1) * (2 * yDim + 1) * (2 * zDim + 1);
		mesh.reserve(n_elements, n_nodes);

		for (Integer i = 0; i <= 2 * xDim; ++i) {
			for (Integer j = 0; j <= 2 * yDim; ++j) {
				for (Integer k = 0; k <= 2 * zDim; ++k) {
					Vector<Real, Dim> p(
							{ static_cast<Real>(i)
									/ static_cast<Real>(2 * xDim),
									static_cast<Real>(j)
											/ static_cast<Real>(2 * yDim),
									static_cast<Real>(k)
											/ static_cast<Real>(2 * zDim), });

					mesh.add_point(p);

				}
			}
		}

		return true;
	}
	default:
		return false;
	}
}

}
}

#endif /* GENERATION_MARS_UNIT_MESH_GENERATION_HPP_ */
