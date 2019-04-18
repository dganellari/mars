#ifndef GENERATION_MARS_UNIT_MESH_GENERATION_HPP_
#define GENERATION_MARS_UNIT_MESH_GENERATION_HPP_

#include "mars_mesh.hpp"

using namespace std;

namespace mars {
namespace unit_generation {

Integer index(const Integer xDim, const Integer yDim, const Integer i,const Integer j,const Integer k){
	//return k+ (2*zDim +1) * (j + i* (2*yDim + 1));
	return i + (2*xDim +1) * (j + k* (2*yDim + 1));
}

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
	case 3: { //building the tetra mesh from the hex27 elem

		const int n_elements = xDim * yDim * zDim * 24; //24 tetrahedrons on one hex27
		const int n_nodes = (2 * xDim + 1) * (2 * yDim + 1) * (2 * zDim + 1);
		mesh.reserve(n_elements, n_nodes);

		for (Integer k = 0; k <= 2 * zDim; ++k) {
			for (Integer j = 0; j <= 2 * yDim; ++j) {
			for (Integer i = 0; i <= 2 * xDim; ++i) {
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

		 int one = 1;
		 int two = 2;
		//add 4 tetrahedra for each face using the center element connectivity which is part of all connectivities.
			for (Integer j = 2; j < (4 * yDim); j += 2) {
				if (j == 2) {
					one = -1 * one;
					two = -1 * two;
				}
		 for (Integer k = 0; k < (4 * zDim); k += 2) {

				for (Integer i = 0; i < (2 * xDim); i += 2) {


					//if(k ==1 && j == 1 && i == 1){
					std::array<Integer, ManifoldDim + 1> nodes;

					nodes[0] = index(yDim, zDim, i, j, k);
					nodes[1] = index(yDim, zDim, i + 2, j, k);
					nodes[2] = index(yDim, zDim, i + 1, j + one, k);
					nodes[3] = index(yDim, zDim, 1, 1, 1);
					mesh.add_elem(nodes);

					nodes[0] = index(yDim, zDim, i, j, k);
					nodes[1] = index(yDim, zDim, i, j + two, k);
					nodes[2] = index(yDim, zDim, i + 1, j + one, k);
					nodes[3] = index(yDim, zDim, 1, 1, 1);
					mesh.add_elem(nodes);

					nodes[0] = index(yDim, zDim, i + 2, j + two, k);
					nodes[1] = index(yDim, zDim, i, j + two, k);
					nodes[2] = index(yDim, zDim, i + 1, j + one, k);
					nodes[3] = index(yDim, zDim, 1, 1, 1);
					mesh.add_elem(nodes);

					nodes[0] = index(yDim, zDim, i + 2, j + two, k);
					nodes[1] = index(yDim, zDim, i + 2, j, k);
					nodes[2] = index(yDim, zDim, i + 1, j + one, k);
					nodes[3] = index(yDim, zDim, 1, 1, 1);
					mesh.add_elem(nodes);
					//}

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
