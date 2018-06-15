
#include <cstdlib>
#include <iostream>
#include "simplex.hpp"

int main(const int argc, const char *argv[])
{
	using namespace mars;

	Matrix<double, 3, 3> m({1, 0, 0,
	                        0, 2, 0,
	                        0, 0, 3});


	Matrix<double, 3, 2> m_sub({1, 0,
	                            0, 2,
	                            0, 0});

	std::cout << det(m)     << std::endl;
	std::cout << det(m_sub) << std::endl;

	auto n = normal(
		{ 2., 0., 0., 0. },
		{ 0., 0., 0., 0. },
		{ 0., 2., 0., 0. },
		{ 0., 0., 2., 0. }
	);

	std::cout << n << std::endl;

	return EXIT_SUCCESS;
}