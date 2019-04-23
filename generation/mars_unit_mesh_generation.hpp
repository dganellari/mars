#ifndef GENERATION_MARS_UNIT_MESH_GENERATION_HPP_
#define GENERATION_MARS_UNIT_MESH_GENERATION_HPP_

#include "mars_mesh.hpp"

using namespace std;

namespace mars {
namespace unit_generation {

constexpr int side_size=6;

Integer index(const Integer xDim, const Integer yDim, const Integer i,
		const Integer j, const Integer k) {
	//return k+ (2*zDim +1) * (j + i* (2*yDim + 1));
	return i + (2 * xDim + 1) * (j + k * (2 * yDim + 1));
}

void add_side(std::vector<Integer>& side, const Integer a,
		const Integer b, const Integer index) {

	if (a != 1 && b != 1) { //add only nodes which are not mid faces or mid edges
		side.push_back(index);
	} else if (a == 1 && b == 1) { // then add only mid faces
		side.push_back(index);
	}
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

		Vector<std::vector<Integer>, side_size> sides;

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

					Integer ind = index(1, 1, i, j, k);

					if (k == 0)
						add_side(sides(0), i, j, ind);

					if (k == 2)
						add_side(sides(1), i, j, ind);

					//if (k == 2 * (zDim - 1))
					if (j == 0)
						add_side(sides(2), i, k, ind);

					//if (j == 2 * (yDim - 1))
					if (j == 2)
						add_side(sides(3), i, k, ind);

					if (i == 0)
						add_side(sides(4), j, k, ind);

					//if (i == 2 * (xDim - 1))
					if (i == 2)
						add_side(sides(5), j, k, ind);

				}
			}
		}

		//build tetrahedra elements from the hex27 faces.
		for (unsigned int i = 0; i < side_size; ++i) {

			Integer tmp = sides(i)[2]; //swap places between the last point and the midface point.
			sides(i)[2] = sides(i)[4];
			sides(i)[4] = tmp;

			for (unsigned int k = 0; k < 4; k++) { //the swap allows for this loop.

				std::array<Integer, ManifoldDim + 1> nodes;

				nodes[0] = sides(i)[k];
				nodes[1] = sides(i)[4]; // midface point always the last element.
				nodes[2] = (k == 3 ? sides(i)[0] : sides(i)[k + 1]); // rotation to catch all combinations.
				nodes[3] = index(yDim, zDim, 1, 1, 1); // the center of the cube.

				mesh.add_elem(nodes);
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
