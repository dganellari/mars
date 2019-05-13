#ifndef GENERATION_MARS_MESH_GENERATION_KOKKOS_HPP_
#define GENERATION_MARS_MESH_GENERATION_KOKKOS_HPP_

#include "mars_mesh_kokkos.hpp"
#include "mars_fwd_kokkos.hpp"
#include <Kokkos_Core.hpp>
#include <functional>

using namespace std;

namespace mars {
namespace generation {

namespace kokkos{

namespace Private {

constexpr int hex_n_sides = 6; // 6 faces in total for the hex27.
constexpr int hex_n_nodes = 27; // 27 nodes for the hex27.
constexpr int hex_side_n_nodes = 9; // 9 nodes per face for the hex27.

//libmesh method to map the sides and nodes.
const unsigned int hex_side_nodes[hex_n_sides][hex_side_n_nodes] = { { 0, 3, 2,
		1, 11, 10, 9, 8, 20 }, // Side 0
		{ 0, 1, 5, 4, 8, 13, 16, 12, 21 }, // Side 1
		{ 1, 2, 6, 5, 9, 14, 17, 13, 22 }, // Side 2
		{ 2, 3, 7, 6, 10, 15, 18, 14, 23 }, // Side 3
		{ 3, 0, 4, 7, 11, 12, 19, 15, 24 }, // Side 4
		{ 4, 5, 6, 7, 16, 17, 18, 19, 25 }  // Side 5
};

/*
 template<Integer Dim, Integer ManifoldDim, class Point_>
 void remove_extra_nodes(Mesh<Dim, ManifoldDim, Point_>& mesh) {

 for (auto it = mesh.points().begin(); it != mesh.points().end();
 no it++){

 if (!(*it).isActive()) {
 mesh.remove_point(it);
 //cout << "removed: " << (it - mesh.points().begin()) << endl;
 } else
 ++it;
 }

 }
 */

template<Integer Dim, Integer ManifoldDim, class Point_>
void remove_extra_nodes(Mesh<Dim, ManifoldDim, Point_>& mesh,
		std::vector<Vector<Real, Dim> >& np, const std::vector<bool>& active) {

	int count = 0;
	for (unsigned int i = 0; i < active.size(); ++i) {
		if (active[i]) {
			np[count] = mesh.point(i);
			++count;
		}

	}

	mesh.setPoints(move(np));

}

Integer index(const Integer xDim, const Integer yDim, const Integer i,
		const Integer j, const Integer k) {
	//return k+ (2*zDim +1) * (j + i* (2*yDim + 1));
	return i + (2 * xDim + 1) * (j + k * (2 * yDim + 1));
}

void add_side(std::vector<Integer>& side, const Integer a, const Integer b,
		const Integer index) {

	if (a != 1 && b != 1) { //add only nodes which are not mid faces or mid edges
		side.push_back(index);
	} else if (a == 1 && b == 1) { // then add only mid faces
		side.push_back(index);
	}
}

void build_hex27(array<Integer, hex_n_nodes>& nodes, const Integer xDim,
		const Integer yDim, const int i, const int j, const int k) {

	nodes[0] = index(xDim, yDim, i, j, k);
	nodes[1] = index(xDim, yDim, i + 2, j, k);
	nodes[2] = index(xDim, yDim, i + 2, j + 2, k);
	nodes[3] = index(xDim, yDim, i, j + 2, k);
	nodes[4] = index(xDim, yDim, i, j, k + 2);
	nodes[5] = index(xDim, yDim, i + 2, j, k + 2);
	nodes[6] = index(xDim, yDim, i + 2, j + 2, k + 2);
	nodes[7] = index(xDim, yDim, i, j + 2, k + 2);
	nodes[8] = index(xDim, yDim, i + 1, j, k);
	nodes[9] = index(xDim, yDim, i + 2, j + 1, k);
	nodes[10] = index(xDim, yDim, i + 1, j + 2, k);
	nodes[11] = index(xDim, yDim, i, j + 1, k);
	nodes[12] = index(xDim, yDim, i, j, k + 1);
	nodes[13] = index(xDim, yDim, i + 2, j, k + 1);
	nodes[14] = index(xDim, yDim, i + 2, j + 2, k + 1);
	nodes[15] = index(xDim, yDim, i, j + 2, k + 1);
	nodes[16] = index(xDim, yDim, i + 1, j, k + 2);
	nodes[17] = index(xDim, yDim, i + 2, j + 1, k + 2);
	nodes[18] = index(xDim, yDim, i + 1, j + 2, k + 2);
	nodes[19] = index(xDim, yDim, i, j + 1, k + 2);
	nodes[20] = index(xDim, yDim, i + 1, j + 1, k);
	nodes[21] = index(xDim, yDim, i + 1, j, k + 1);
	nodes[22] = index(xDim, yDim, i + 2, j + 1, k + 1);
	nodes[23] = index(xDim, yDim, i + 1, j + 2, k + 1);
	nodes[24] = index(xDim, yDim, i, j + 1, k + 1);
	nodes[25] = index(xDim, yDim, i + 1, j + 1, k + 2);
	nodes[26] = index(xDim, yDim, i + 1, j + 1, k + 1);
}
}

template<Integer Dim, Integer ManifoldDim>
bool generate_cube(Parallel_Mesh<Dim, ManifoldDim>& mesh, const Integer xDim,
		const Integer yDim, const Integer zDim) {

	using namespace mars::generation::Private;
	using Elem = mars::Simplex<Dim, ManifoldDim>;

	using namespace Kokkos;


	assert(ManifoldDim <= Dim);
	assert(Dim <= 3);

	switch (ManifoldDim) {

	/*case 0: {

	 assert(xDim == 0);
	 assert(yDim == 0);
	 assert(zDim == 0);

	 mesh.reserve(1, 1);

	 Vector<Real, Dim> p( { 0.0, 0.0, 0.0 });

	 int index = mesh.add_point(p);

	 Simplex<Dim, ManifoldDim> e;
	 e.node(index);

	 mesh.add_elem(e);

	 return true;
	 }*/

	case 1: {

		assert(xDim != 0);
		assert(yDim == 0);
		assert(zDim == 0);


		const int n_elements = xDim;
		const int n_nodes = xDim + 1;
		mesh.reserve(n_elements, n_nodes);


		//todo: try lambda equivalent

		/*parallel_for(n_nodes, [=,&mesh] __device__ __host__(const size_t index){
			mesh.add_point1(index,xDim);
		});*/

		//:todo try to avoid it by calling mesh.AddPoints and remove the view as a parameter but instead get it from structure directly

		ViewMatrixType<Real> points = mesh.get_view_points();
		ViewMatrixType<Integer> elems = mesh.get_view_elems();
		ViewVectorType<bool> active = mesh.get_view_active();

		parallel_for(n_nodes,
				typename Parallel_Mesh<Dim, ManifoldDim>::AddPoint(points,
						xDim));

		Cuda::fence();

		parallel_for(n_elements,
				typename Parallel_Mesh<Dim, ManifoldDim>::AddElem(elems,
						active));

		Cuda::fence();

		return true;
	}

	/*case 2: {

		assert(xDim != 0);
		assert(yDim != 0);
		assert(zDim == 0);

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

		assert(xDim != 0);
		assert(yDim != 0);
		assert(zDim != 0);

		const int n_elements = xDim * yDim * zDim * 24; //24 tetrahedrons on one hex27
		const int n_nodes = (2 * xDim + 1) * (2 * yDim + 1) * (2 * zDim + 1);

		const int n_tetra_nodes = 5 * (xDim * yDim * zDim)
				+ 2 * (xDim * yDim + xDim * zDim + yDim * zDim)
				+ (xDim + yDim + zDim) + 1;

		std::vector<bool> active_nodes(n_nodes);

		// mesh points need resize because the point vector is getting accessed using the [] operator.
		//First the center node is added somewhere at the end of the vector and then all the other nodes are added to it.
		mesh.reserve_elements(n_elements);
		mesh.resize_points(n_tetra_nodes);

		std::array<Integer, hex_side_n_nodes> side;

		int el_id = 1;
		for (Integer k = 0; k < 2 * zDim; k += 2) {
			for (Integer j = 0; j < 2 * yDim; j += 2) {
				for (Integer i = 0; i < 2 * xDim; i += 2) {

					//build the hex27 element which serves to generalize the idea to many hex27.
					//without generating the hex27 element first there is no easy way to create the sides.
					//Locally, using local element indexing it is possible.
					std::array<Integer, hex_n_nodes> hex;
					build_hex27(hex, xDim, yDim, i, j, k);

					//add center of the hex to the new points.
					int centerHex = n_nodes / 2 + el_id;

					Vector<Real, Dim> p(
							{ static_cast<Real>(i+1)
									/ static_cast<Real>(2 * xDim),
									static_cast<Real>(j+1)
											/ static_cast<Real>(2 * yDim),
									static_cast<Real>(k+1)
											/ static_cast<Real>(2 * zDim), });

					mesh.point(centerHex) = p;

					//build tetrahedra elements from the hex27 faces.
					for (unsigned int i = 0; i < hex_n_sides; ++i) {

						//use the hex27 element local connectivity to build the tetrahedra connectivity
						for (unsigned int j = 0; j < hex_side_n_nodes; ++j) {
							side[j] = hex[hex_side_nodes[i][j]];
						}

						//build 4 tetrahedra out of one face.
						for (unsigned int k = 0; k < 4; k++) {

							auto& e = mesh.add_elem();

							e.nodes[0] = side[k] / 2;
							e.nodes[1] = side[8] / 2; // midface point always the last element.
							e.nodes[2] = (k == 3 ? side[0] : side[k + 1]) / 2; // rotation to catch all combinations.
							e.nodes[3] = centerHex; // the center of the cube.

							active_nodes[side[k]] = true;
							active_nodes[side[8]] = true;
							active_nodes[(k == 3 ? side[0] : side[k + 1])] =
									true;
							//active_nodes[hex[26]]=true;
						}

					}
					++el_id;
				}
			}
		}

		//first the element indices as above
		//and then at this moment add the only the needed points to avoid extra node removal.
		for (Integer k = 0; k <= 2 * zDim; ++k) {
			for (Integer j = 0; j <= 2 * yDim; ++j) {
				for (Integer i = 0; i <= 2 * xDim; ++i) {

					Integer in = index(xDim, yDim, i, j, k);

					if (active_nodes[in]) {

						mesh.point(in / 2)[0] = static_cast<Real>(i)
												/ static_cast<Real>(2 * xDim);
						mesh.point(in / 2)[1] = static_cast<Real>(j)
												/ static_cast<Real>(2 * xDim);
						mesh.point(in / 2)[2] = static_cast<Real>(k)
												/ static_cast<Real>(2 * xDim);
					}

				}
			}
		}

		return true;
	}*/
	default: {

		std::cerr << "Not implemented for other dimensions yet" << std::endl;
		return false;
	}
	}
}

/*bool generate_square(Mesh2& mesh, const Integer xDim, const Integer yDim) {
	return generate_cube(mesh, xDim, yDim, 0);
}

bool generate_line(Mesh<1, 1>& mesh, const Integer xDim) {
	return generate_cube(mesh, xDim, 0, 0);
}*/

/*bool generate_point(Mesh<1, 0>& mesh) {
 return generate_cube(mesh,0,0,0);
 }*/

//different approach which works fine for the unit_cube and would have been faster
//but unfortunately it does not generalize.
Mesh3 generate_unit_cube() {

	using namespace mars::generation::Private;

	Mesh3 mesh;

	constexpr Integer xDim = 1, yDim = 1, zDim = 1, Dim = 3;

	const int n_elements = xDim * yDim * zDim * 24; //24 tetrahedrons on one hex27
	const int n_nodes = (2 * xDim + 1) * (2 * yDim + 1) * (2 * zDim + 1);
	mesh.reserve(n_elements, n_nodes);

	Vector<std::vector<Integer>, hex_n_sides> sides;

	for (Integer k = 0; k <= 2 * zDim; ++k) {
		for (Integer j = 0; j <= 2 * yDim; ++j) {
			for (Integer i = 0; i <= 2 * xDim; ++i) {
				Vector<Real, Dim> p(
						{ static_cast<Real>(i) / static_cast<Real>(2 * xDim),
								static_cast<Real>(j)
										/ static_cast<Real>(2 * yDim),
								static_cast<Real>(k)
										/ static_cast<Real>(2 * zDim), });

				mesh.add_point(p);

				Integer ind = index(1, 1, i, j, k);

				//build the faces on the fly
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
	for (unsigned int i = 0; i < hex_n_sides; ++i) {

		Integer tmp = sides(i)[2]; //swap places between the last point and the midface point.
		sides(i)[2] = sides(i)[4];
		sides(i)[4] = tmp;

		for (unsigned int k = 0; k < 4; k++) {

			std::array<Integer, 4> nodes;

			nodes[0] = sides(i)[k];
			nodes[1] = sides(i)[4]; // midface point always the last element.
			nodes[2] = (k == 3 ? sides(i)[0] : sides(i)[k + 1]); // rotation to catch all combinations.
			nodes[3] = index(yDim, zDim, 1, 1, 1); // the center of the cube.

			mesh.add_elem(nodes);
		}

	}

	return mesh;

}

}
}
}

#endif /* GENERATION_MARS_MESH_GENERATION_HPP_ */
