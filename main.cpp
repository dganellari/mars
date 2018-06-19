
#include <cstdlib>
#include <iostream>
#include <cassert>
#include "simplex.hpp"
#include "lagrange_element.hpp"
#include "mesh.hpp"

void test_midpoint_index()
{
	using namespace mars;

	//tri
	{
		const auto i01 = midpoint_index<2>(0, 1); assert(i01 == 3);
		const auto i02 = midpoint_index<2>(0, 2); assert(i02 == 4);
		const auto i03 = midpoint_index<2>(1, 2); assert(i03 == 5);
	}

	//tet
	{
		const auto i01 = midpoint_index<3>(0, 1); assert(i01 == 4);
		const auto i02 = midpoint_index<3>(0, 2); assert(i02 == 5);
		const auto i03 = midpoint_index<3>(0, 3); assert(i03 == 6);

		const auto i04 = midpoint_index<3>(1, 2); assert(i04 == 7);
		const auto i13 = midpoint_index<3>(1, 3); assert(i13 == 8);

		const auto i23 = midpoint_index<3>(2, 3); assert(i23 == 9);
	}

	//pentatope
	{
		const auto i01 = midpoint_index<4>(0, 1); assert(i01 == 5);
		const auto i02 = midpoint_index<4>(0, 2); assert(i02 == 6);
		const auto i03 = midpoint_index<4>(0, 3); assert(i03 == 7);
		const auto i04 = midpoint_index<4>(0, 4); assert(i04 == 8);

		const auto i12 = midpoint_index<4>(1, 2); assert(i12 == 9);
		const auto i13 = midpoint_index<4>(1, 3); assert(i13 == 10);
		const auto i14 = midpoint_index<4>(1, 4); assert(i14 == 11);

		const auto i23 = midpoint_index<4>(2, 3); assert(i23 == 12);
		const auto i24 = midpoint_index<4>(2, 4); assert(i24 == 13);

		const auto i34 = midpoint_index<4>(3, 4); assert(i34 == 14);
	}
}

void test_red_refinement_interpolator()
{
	using namespace mars;

	SimplexInterpolator<2> mat2;
	red_refinement_interpolator<2>(mat2);
	mat2.describe(std::cout);

	SimplexInterpolator<3> mat3;
	red_refinement_interpolator<3>(mat3);
	mat3.describe(std::cout);

	SimplexInterpolator<4> mat4;
	red_refinement_interpolator<4>(mat4);
	mat4.describe(std::cout);
}

void test_red_refinement()
{
	using namespace mars;

	Pentatope4 ptope;
	ptope.nodes[0] = 0;
	ptope.nodes[1] = 1;
	ptope.nodes[2] = 2;
	ptope.nodes[3] = 3;
	ptope.nodes[4] = 4;

	std::vector<Vector4r> points(
	{
		{ 0., 0., 0., 0. },
		{ 2., 0., 0., 0. },
		{ 0., 2., 0., 0. },
		{ 0., 0., 2., 0. },
		{ 0., 0., 0., 2. }
	});

	SimplexInterpolator<4> interp;
	std::array<Pentatope4, 16> sub_simplices;
	std::vector<Vector4r> sub_points;

	red_refinement<4, 4, 16>(ptope, points, sub_simplices, sub_points, interp);

}

void test_mfem_mesh()
{	
	using namespace mars;

	Mesh<4, 4> mesh;
	read_mesh("../data/cube4d_24.MFEM", mesh);
	// mesh.describe(std::cout);

	mesh.build_dual_graph();
	// mesh.describe_dual_graph(std::cout);

	mesh.refine_element(0);
	mesh.describe(std::cout, true);
	// mesh.describe(std::cout, false);
}

int main(const int argc, const char *argv[])
{
	using namespace mars;

	// Tetrahedron4 tet; 
	// tet.nodes[0] = 0;
	// tet.nodes[1] = 1;
	// tet.nodes[2] = 2;
	// tet.nodes[3] = 3;

	// Triangle4 tri4;
	// tri4.nodes[0] = 1;
	// tri4.nodes[1] = 2;
	// tri4.nodes[2] = 3;

	// Line2 line;
	// line.nodes[0] = 0;
	// line.nodes[1] = 1;

	// Triangle2 tri2;
	// tri2.nodes[0] = 0;
	// tri2.nodes[1] = 1;
	// tri2.nodes[2] = 2;

	// Pentatope4 ptope;
	// ptope.nodes[0] = 0;
	// ptope.nodes[1] = 1;
	// ptope.nodes[2] = 2;
	// ptope.nodes[3] = 3;
	// ptope.nodes[4] = 4;

	// std::vector<Vector4r> points4(
	// {
	// 	{ 0., 0., 0., 0. },
	// 	{ 2., 0., 0., 0. },
	// 	{ 0., 2., 0., 0. },
	// 	{ 0., 0., 2., 0. },
	// 	{ 0., 0., 0., 2. }
	// });

	// std::vector<Vector2r> points2(
	// {
	// 	{ 0., 0. },
	// 	{ 1., 0. },
	// 	{ 0., 1. }
	// });

	// std::cout << "===================" << std::endl;
	// std::cout << "vol(line2): " << volume(line, points2) << std::endl;
	// std::cout << "vol(tri2):  " << volume(tri2, points2) << std::endl;
	// std::cout << "vol(tri4):  " << volume(tri4, points4) << std::endl;
	// std::cout << "vol(tet4):  " << volume(tet, points4) << std::endl;
	// std::cout << "vol(ptope): " << volume(ptope, points4) << std::endl;
	// std::cout << "===================" << std::endl;

	// Matrix<Real, 4, 4> J44;
	// jacobian(ptope, points4, J44); 
	// std::cout << "J(ptope):\n" << J44 << std::endl;

	// Matrix<Real, 4, 2> Jtri4;
	// jacobian(tri4, points4, Jtri4);
	// check_and_fix_jac(Jtri4);
	// std::cout << "J(tri4):\n" << Jtri4 << std::endl;

	// auto n1 = normal(line, points2);
	// auto n3 = normal(tri4, points4);
	// auto n4 = normal(tet, points4);

	// std::cout << n1 << std::endl;
	// std::cout << n3 << std::endl;
	// std::cout << n4 << std::endl;

	// test_midpoint_index();
	// test_red_refinement_interpolator();
	// test_red_refinement();
	test_mfem_mesh();
	return EXIT_SUCCESS;
}


