#ifndef MARS_BISECTION_HPP
#define MARS_BISECTION_HPP

namespace mars {
	template<Integer Dim, Integer ManifoldDim>
	class Mesh;

	template<Integer Dim, Integer ManifoldDim>
	class Bisection {
	public:
		Bisection(Mesh<Dim, ManifoldDim> &mesh)
		: mesh(mesh)
		{}

		void uniformly_refine(const n_levels)
		{

		}

		void refine(const std::vector<Integer> &elems)
		{

		}
		
	private:
		Mesh<Dim, ManifoldDim> &mesh;
	};
}

#endif //MARS_BISECTION_HPP
