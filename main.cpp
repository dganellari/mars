
#include <cstdlib>
#include <iostream>
#include "simplex.hpp"

int main(const int argc, const char *argv[])
{
	using namespace mars;

	Tetrahedron4 tet; 
	tet.nodes[0] = 0;
	tet.nodes[1] = 1;
	tet.nodes[2] = 2;
	tet.nodes[3] = 3;

	Triangle4 tri4;
	tri4.nodes[0] = 1;
	tri4.nodes[1] = 2;
	tri4.nodes[2] = 3;

	Line2 line;
	line.nodes[0] = 0;
	line.nodes[1] = 1;

	Triangle2 tri2;
	tri2.nodes[0] = 0;
	tri2.nodes[1] = 1;
	tri2.nodes[2] = 2;

	std::vector<Vector4r> points4(
	{
		{ 2., 0., 0., 0. },
		{ 0., 0., 0., 0. },
		{ 0., 2., 0., 0. },
		{ 0., 0., 2., 0. }
	});

	std::vector<Vector2r> points2(
	{
		{ 0., 0. },
		{ 1., 0. },
		{ 0., 1. }
	});

	std::cout << "===================" << std::endl;
	// std::cout << "vol(line2): " << volume(line, points2) << std::endl;
	// std::cout << "vol(tri2):  " << volume(tri2, points2) << std::endl;
	std::cout << "vol(tri4):  " << volume(tri4, points4) << std::endl;
	// std::cout << "vol(tet4):  " << volume(tet, points4) << std::endl;
	std::cout << "===================" << std::endl;


	Matrix<Real, 4, 2> Jtri4;
	jacobian(tri4, points4, Jtri4);
	check_and_fix_jac(Jtri4);
	std::cout << "J(tri4):\n" << Jtri4 << std::endl;

	auto n1 = normal(line, points2);
	auto n3 = normal(tri4, points4);
	auto n4 = normal(tet, points4);

	std::cout << n1 << std::endl;
	std::cout << n3 << std::endl;
	std::cout << n4 << std::endl;
	return EXIT_SUCCESS;
}