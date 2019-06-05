#ifndef GENERATION_MARS_MESH_GENERATION_KOKKOS_HPP_
#define GENERATION_MARS_MESH_GENERATION_KOKKOS_HPP_
#include <functional>
#include "mars_mesh_kokkos.hpp"


#include "mars_fwd_kokkos.hpp"
#include <Kokkos_Core.hpp>

namespace mars {
namespace generation {
namespace kokkos{
namespace private_ {

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

void build_hex27(std::array<Integer, hex_n_nodes>& nodes, const Integer xDim,
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
bool generate_cube(ParallelMesh<Dim, ManifoldDim>& mesh, const Integer xDim,
		const Integer yDim, const Integer zDim) {

	using namespace mars::generation::private_;
	using Elem = mars::Simplex<Dim, ManifoldDim>;

	using namespace Kokkos;

	assert(Dim <= 3);
	assert(ManifoldDim <= Dim);

	mesh.generate_points(xDim,yDim);

	mesh.generate_elements(xDim,yDim);
}

bool generate_line(ParallelMesh<1, 1>& mesh, const Integer xDim) {
	return generate_cube(mesh, xDim, 0, 0);
}

bool generate_square(ParallelMesh<2, 2>& mesh, const Integer xDim, const Integer yDim) {
	return generate_cube(mesh, xDim, yDim, 0);
}

/*
bool generate_point(Mesh<1, 0>& mesh) {
 return generate_cube(mesh,0,0,0);
 }
*/


}
}
}

#endif /* GENERATION_MARS_MESH_GENERATION_HPP_ */
