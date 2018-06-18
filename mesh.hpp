#ifndef MARS_MESH_HPP
#define MARS_MESH_HPP

#include "simplex.hpp"

#include <vector>
#include <array>
#include <fstream>
#include <sstream>

namespace mars {
	template<Integer Dim, Integer ManifoldDim = Dim>
	class Mesh {
	public:
		void reserve(
			const std::size_t n_elements,
			const std::size_t n_points)
		{
			elements_.reserve(n_elements);
			points_.reserve(n_points);
		}

		inline Simplex<Dim, ManifoldDim> &elem(const Integer id)
		{
			return elements_[id];
		}

		inline Integer add_point(const Vector<Real, Dim> &point)
		{
			points_.push_back(point);
			return points_.size() - 1;
		}

		inline Integer add_elem(const Simplex<Dim, ManifoldDim> &elem)
		{
			elem.id = elements_.size();
			elements_.push_back(elem);
			return elem.id;
		}

		template<std::size_t NNodes>
		Integer add_elem(const std::array<Integer, NNodes> &nodes)
		{
			static_assert(NNodes == std::size_t(ManifoldDim + 1), "does not have the correct number of nodes");
			elements_.emplace_back();
			auto &e = elements_.back();
			e.id = elements_.size() - 1;
			e.nodes = nodes;
			return e.id;
		}
	
		inline void set_refinement_flag(const Integer &element_id, const Integer flag)
		{
			refinement_flag_[element_id] = flag;
		}

		inline void refine(const Integer n_levels)
		{

		}

		void repair()
		{
			for(std::size_t i = 0; i < elements_.size(); ++i) {
				auto &e = elements_[i];
				const Real vol = volume(e, points_);

				if(vol < 0.) {
					if(Dim == 4) {
						std::swap(e.nodes[3], e.nodes[4]);
						const Real vol_after = volume(e, points_);
						assert(vol_after > 0.);
					}
				}
			}
		}

		void describe(std::ostream &os) const
		{
			for(std::size_t i = 0; i < elements_.size(); ++i) {
				const auto &e = elements_[i];
				const Real vol = volume(e, points_);
				os << "---------------------------------\n";
				os << "[" << i << "]: vol: " << vol << ", ";
				for(auto v : e.nodes) {
					os << " " << v;
				}

				os << "\n";

				Simplex<Dim, ManifoldDim-1> side;
				Matrix<Real, Dim, Dim-1> J;

				os << "sides:";
				for(Integer k = 0; k < n_sides(e); ++k) {
					e.side(k, side);
					os << "==============\n";
					jacobian(side, points_, J);
					const Real area   = volume(side, points_);
					const Real u_area = unsigned_volume(side, points_);
					
					J.describe(os);
					os << area << " == " << u_area << std::endl;
				}

				os << "---------------------------------\n";
				os << "\n";
			}

			for(std::size_t i = 0; i < points_.size(); ++i) {
				points_[i].describe(os);
			}
		}

	private:
		std::vector< Simplex<Dim, ManifoldDim> > elements_;
		std::vector< Vector<Real, Dim> > points_;
		
		//refinement
		std::vector<Integer> refinement_flag_;
		std::vector<std::array<Integer, ManifoldDim+1>> dual_graph_;
	};

	bool read_mesh(const std::string &path, Mesh<4, 4> &mesh)
	{
		std::ifstream is(path);
		if(!is.good()) {
			return false;
		}

		int dim = -1;
		int n_elements = -1;
		int n_nodes = -1;
		int n_coords = -1;

		std::string line;
		while(is.good()) {
			std::getline(is, line);

			if(line == "dimension") {
				std::getline(is, line);
				dim = atoi(line.c_str());
				assert(dim == 4);
			} else if(line == "elements") {
				std::getline(is, line);
				n_elements = atoi(line.c_str());

				for(Integer i = 0; i < n_elements; ++i) {
					assert(is.good());
					std::getline(is, line);
					std::stringstream ss(line);
					int attr, type;

					std::array<Integer, 5> nodes;
					ss >> attr >> type;

					for(Integer k = 0; k < 5; ++k) {
						ss >> nodes[k];
					}

					mesh.add_elem(nodes);
				}
			} else if(line == "vertices") {
				std::getline(is, line);
				n_nodes = atoi(line.c_str());
				std::getline(is, line);
				n_coords = atoi(line.c_str());
				assert(n_coords == 4);

				Vector<Real, 4> p;
				for(Integer i = 0; i < n_nodes; ++i) {
					assert(is.good());

					for(Integer k = 0; k < 4; ++k) {
						is >> p(k);
					}

					mesh.add_point(p);
				}

			}
		}

		is.close();

		mesh.repair();
		return true;
	}

	bool mesh_hyper_cube(
		const std::array<Integer, 4> &dims,
		const Vector<Real, 4> &lobo,
		const Vector<Real, 4> &upbo,
		const Mesh<4, 4> &mesh)
	{

		return false;
	}
}

#endif //MARS_MESH_HPP
