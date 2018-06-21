#ifndef MARS_MESH_HPP
#define MARS_MESH_HPP

#include "simplex.hpp"
#include "edge_element_map.hpp"
#include "edge_node_map.hpp"
#include "dual_graph.hpp"
#include "red_green_refinement.hpp"

namespace moonolith {
	using Integer = mars::Integer;
}

#include "moonolith_config.hpp"
#include "moonolith_mesh.hpp"
// #include "moonolith_svg_canvas.hpp"
#include "moonolith_eps_canvas.hpp"
#include "moonolith_plotter.hpp"
#include "moonolith_func_to_color.hpp"

#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <memory>
#include <algorithm>  

namespace mars {

	template<Integer Dim, Integer ManifoldDim = Dim>
	class Mesh {
	public:
		void reserve(
			const std::size_t n_elements,
			const std::size_t n_points)
		{
			elements_.reserve(n_elements);
			active_.reserve(n_elements);
			points_.reserve(n_points);
		}

		inline Simplex<Dim, ManifoldDim> &elem(const Integer id)
		{
			return elements_[id];
		}

		inline const Simplex<Dim, ManifoldDim> &elem(const Integer id) const
		{
			return elements_[id];
		}

		inline bool is_active(const Integer id) const
		{
			return active_[id];
		}

		inline void set_active(const Integer id, const bool val)
		{
			active_[id] = val;
		}

		inline Integer add_point(const Vector<Real, Dim> &point)
		{
			points_.push_back(point);
			return points_.size() - 1;
		}

		inline const Vector<Real, Dim> &point(const Integer i) const
		{
			return points_[i];
		}

		const std::vector<Vector<Real, Dim>> &points() const
		{
			return points_;
		}

		inline Integer add_elem(const Simplex<Dim, ManifoldDim> &elem)
		{
			auto id = elements_.size();
			elements_.push_back(elem);
			elements_.back().id = id;
			active_.push_back(true);
			assert(elements_.back().id == id);
			return elements_.back().id;
		}

		template<std::size_t NNodes>
		Integer add_elem(const std::array<Integer, NNodes> &nodes)
		{
			static_assert(NNodes == std::size_t(ManifoldDim + 1), "does not have the correct number of nodes");
			elements_.emplace_back();
			auto &e = elements_.back();
			e.id = elements_.size() - 1;
			e.nodes = nodes;
			active_.push_back(true);
			assert(e.id == elements_.size() - 1);
			return e.id;
		}

		inline void points(const Integer element_id, std::vector<Vector<Real, Dim>> &pts)
		{
			auto &e = elements_[element_id];
			pts.resize(ManifoldDim + 1);
			
			for(Integer i = 0; i < ManifoldDim + 1; ++i) {
				pts[i] = points_[e.nodes[i]];
			}
		}

		inline void deactivate_children(const Integer id)
		{
			for(auto c : elem(id).children) {
				active_[c] = false;
			}
		}

		void repair_element(const Integer element_id, const bool verbose = false)
		{
			assert(element_id >= 0);

			auto &e = elements_[element_id];
			const Real vol = volume(e, points_);

			if(vol < 0.) {
				if(verbose) {
					std::cout << element_id << " has negative volume" << std::endl;
				}

				switch(Dim) {
					case 1:
					{
						std::swap(e.nodes[0], e.nodes[1]);
						const Real vol_after = volume(e, points_);
						assert(vol_after > 0.);
						break;
					}
					case 2:
					{
						std::swap(e.nodes[1], e.nodes[2]);
						const Real vol_after = volume(e, points_);
						assert(vol_after > 0.);
						break;
					}
					case 3:
					{
						std::swap(e.nodes[2], e.nodes[3]);
						const Real vol_after = volume(e, points_);
						assert(vol_after > 0.);
						break;
					}
					case 4:
					{
						std::swap(e.nodes[3], e.nodes[4]);
						const Real vol_after = volume(e, points_);
						assert(vol_after > 0.);
						break;
					}
					default: 
					{
						assert(false && "implement me");
						break;
					}
				}
			}
		}

		void repair(const bool verbose = false)
		{
			for(std::size_t i = 0; i < elements_.size(); ++i) {
				repair_element(i, verbose);
			}
		}

		void describe(std::ostream &os, const bool print_sides = false) const
		{
			for(std::size_t i = 0; i < elements_.size(); ++i) {
				if(!active_[i]) continue;

				const auto &e = elements_[i];
				const Real vol = volume(e, points_);
				const auto b   = barycenter(e, points_);

				os << "---------------------------------\n";
				os << "[" << i << "]: vol: " << vol << ", ";
				for(auto v : e.nodes) {
					os << " " << v;
				}

				os << "\n";

				if(print_sides) {
					Simplex<Dim, ManifoldDim-1> side;
					Matrix<Real, Dim, Dim-1> J;

					os << "sides:\n";
					for(Integer k = 0; k < n_sides(e); ++k) {
						e.side(k, side);
						os << "==============\n";
						jacobian(side, points_, J);

						const auto n = normal(side, points_);
						const auto sign = dot(points_[side.nodes[0]] - b, n) > 0? 1 : -1;
						const Real u_area = unsigned_volume(side, points_);
						const Real area   = sign * u_area;

						J.describe(os);
						os << area << " == " << u_area << std::endl;
					}
				}

				os << "---------------------------------\n";
				os << "\n";
			}

			for(std::size_t i = 0; i < points_.size(); ++i) {
				os << i << ") ";
				points_[i].describe(os);
			}
		}

		inline Integer n_nodes() const
		{
			return points_.size();
		}

		inline Integer n_elements() const
		{
			return elements_.size();
		}

		inline Integer n_active_elements() const
		{
			Integer ret = 0;
			for(auto a : active_) {
				ret += a;
			}

			return ret;
		}

		bool have_common_side(const Integer e_index_1, const Integer e_index_2) const
		{
			const auto &e1 = elem(e_index_1);
			const auto &e2 = elem(e_index_2);

			Integer n_common_nodes = 0;
			for(Integer i = 0; i < ManifoldDim + 1; ++i) {
				for(Integer j = 0; j < ManifoldDim + 1; ++j) {
					n_common_nodes += e1.nodes[i] == e2.nodes[j];
				}
			}

			assert(n_common_nodes <= ManifoldDim);
			return n_common_nodes == ManifoldDim;
		}

		Integer common_side_num(const Integer e_index_1, const Integer e_index_2) const
		{
			const auto &e1 = elem(e_index_1);
			const auto &e2 = elem(e_index_2);

			Simplex<Dim, ManifoldDim-1> side;
			for(Integer k = 0; k < n_sides(e1); ++k) {
				e1.side(k, side);				
				
				Integer nn = 0;

				for(Integer i = 0; i < ManifoldDim; ++i) {
					const auto side_node = side.nodes[i];

					for(Integer j = 0; j < ManifoldDim + 1; ++j) {
						if(side_node == e2.nodes[j]) {
							nn++;
							break;
						}
					}
				}

				if(nn == ManifoldDim) {
					return k;
				}
			}

			assert(false);
			return INVALID_INDEX;
		}

		void describe_dual_graph(std::ostream &os) const
		{
			dual_graph_.describe(os);
		}

		Integer n_boundary_sides() const
		{
			assert( !dual_graph_.empty() && "requires that build_dual_graph is called first");

			Integer ret = 0;
			for(Integer i = 0; i < n_elements(); ++i) {
				if(!active_[i]) continue;

				const auto &e = elem(i);
				const auto &e_adj = dual_graph_.adj(i);
				for(Integer k = 0; k < e_adj.size(); ++k) {
					const Integer j = e_adj[k];
					if(j == INVALID_INDEX) {
						ret++;
					}
				}
			}

			return ret;
		}

		bool check_side_ordering() const
		{
			assert( !dual_graph_.empty() && "requires that build_dual_graph is called first");
			
			if(ManifoldDim == 4) {
				std::cerr << "not implemented for 4d yet" << std::endl;
				return false;
			}

			Simplex<Dim, ManifoldDim-1> side, other_side;

			for(Integer i = 0; i < n_elements(); ++i) {
				if(!active_[i]) continue;
				
				const auto &e = elem(i);
				const auto &e_adj = dual_graph_.adj(i);
				for(Integer k = 0; k < e_adj.size(); ++k) {
					const Integer j = e_adj[k];
					if(j == INVALID_INDEX) continue;
					e.side(k, side);

					const auto &other = elem(j);
					const auto &other_adj = dual_graph_.adj(j);


					Integer other_side_index = 0;
					{
						auto it = std::find(other_adj.begin(), other_adj.end(), i);
						

						if(it == other_adj.end()) {
							std::cerr << "Bad dual graph for " <<  i << " <-> " << j << std::endl;
							assert(it != other_adj.end());
							return false;
						}

						other_side_index = std::distance(other_adj.begin(), it);
						other.side(other_side_index, other_side);
					}

					auto it = std::find(other_side.nodes.begin(), other_side.nodes.end(), side.nodes[0]);
					assert(it != other_side.nodes.end());

					Integer other_offset = std::distance(other_side.nodes.begin(), it);

					for(Integer q = 0; q < ManifoldDim; ++q) {
						Integer other_q = other_offset - q;
						
						if(other_q < 0) {
							other_q += ManifoldDim;
						}

						if(side.nodes[q] != other_side.nodes[other_q]) {
							std::cerr << "common face not matching for (" << i << ", " << k << ") and (" << j << ", " << other_side_index << ")" << std::endl;
							std::cerr << "[ ";
							for(auto s : side.nodes) {
								std::cerr << s << " ";
							}
							std::cerr << " ]\n";

							std::cerr << "[ ";
							for(auto s : other_side.nodes) {
								std::cerr << s << " ";
							}
							std::cerr << " ]\n";
							break;
						}
					}
				}
			}

			return true;
		}

		DualGraph<ManifoldDim> &dual_graph() { return dual_graph_; }

		void update_dual_graph(const bool force = false)
		{
			dual_graph_.update(*this, force);
		}

		void build_dual_graph()
		{
			update_dual_graph();
		}

		inline Integer root(const Integer id) const
		{	
			Integer current_id = id;
			while(elem(current_id).parent_id != INVALID_INDEX) {
				current_id = elem(current_id).parent_id;
			}

			return current_id;
		}

		std::vector<Integer> &tags() { return tags_; }
		const std::vector<Integer> &tags() const { return tags_; }

	private:
		std::vector< Simplex<Dim, ManifoldDim> > elements_;
		std::vector< Vector<Real, Dim> > points_;
		std::vector<Integer> tags_;
		DualGraph<ManifoldDim> dual_graph_;
		std::vector<bool> active_;
	};

	enum PlotFun
	{
		PLOT_ROOT = 0,
		PLOT_TAG = 1,
		PLOT_ID = 2,
		PLOT_PARENT = 3,
		PLOT_PARENT_TAG = 4
	};

	inline void flag_color(const Integer flag, Real &r, Real &g, Real &b)
	{
		r = 0.; g = 0.; b = 0;
		switch(flag)
		{
			case RED: {
				r = 1.;
				break;
			}
			case GREEN_1:
			case GREEN_2:
			case GREEN_3:
			{
				g = 1.;
				break;
			}
			default: {
				r = 1.; 
				g = 1.;
				b = 1.;
				break;
			}
		}
	}


	template<Integer Dim>
	bool write_mesh(
		const std::string &path,
		const Mesh<Dim, 2> &mesh,
		const Real scale_factor = 1.,
		const PlotFun plot_fun = PLOT_ROOT)
	{

		moonolith::Mesh m;
		m.dim = Dim;

		m.points.resize(mesh.n_nodes() * Dim);
		m.el_index.resize(mesh.n_active_elements() * 3);

		for(Integer i = 0; i < mesh.n_nodes(); ++i) {
			for(Integer d = 0; d < Dim; ++d) {
				m.points[i * Dim + d] = mesh.point(i)(d) * scale_factor;
			}
		}

		m.elem_type.resize(mesh.n_active_elements());
		std::fill(m.elem_type.begin(), m.elem_type.end(), moonolith::ElemType::TRI3);
		m.uniform_elem_type = moonolith::ElemType::TRI3;
		m.has_uniform_elem_type = true;

		m.el_ptr.resize(m.elem_type.size() + 1);

		m.el_ptr[0] = 0;
		for(std::size_t i = 1; i < m.el_ptr.size(); ++i) {
			m.el_ptr[i] = m.el_ptr[i - 1] + 3;
		}

		Integer k = 0;
		for(std::size_t i = 0; i < mesh.n_elements(); ++i) {
			if(mesh.is_active(i)) {
				for(Integer j = 0; j < 3; ++j) {
					m.el_index[k * 3 + j] = mesh.elem(i).nodes[j];
				}

				k++;
			}
		}

		// moonolith::SVGCanvas canvas;
		moonolith::EPSCanvas canvas;
		canvas.set_line_width(0.1/mesh.n_active_elements());

		std::vector<Real> f, hsv;

		if(plot_fun == PLOT_TAG || plot_fun == PLOT_PARENT_TAG) {
			f.resize(mesh.n_active_elements() * 3);
		} else {
			f.resize(mesh.n_active_elements());
		}
		for(std::size_t i = 0, k = 0; i < mesh.n_elements(); ++i) {
			if(mesh.is_active(i)) {
				switch(plot_fun) {
					case PLOT_ROOT: {
						f[k++] = mesh.root(i);
						break;
					}

					case PLOT_TAG: {
						Real r = 0., g = 0., b = 0;
						flag_color(mesh.tags()[i], r, g, b);
						f[k * 3]     = r;
						f[k * 3 + 1] = g;
						f[k * 3 + 2] = b; 

						k++;
						break;
					}

					case PLOT_PARENT_TAG: {
						Real r = 0., g = 0., b = 0;
						if(mesh.elem(i).parent_id != INVALID_INDEX) {
							flag_color(mesh.tags()[mesh.elem(i).parent_id], r, g, b);
						} else {
							flag_color(NONE, r, g, b);
						}

						f[k * 3]     = r;
						f[k * 3 + 1] = g;
						f[k * 3 + 2] = b; 

						k++;
						break;
					}

					case PLOT_PARENT: {
						if(mesh.elem(i).parent_id != INVALID_INDEX) {
							f[k++] = mesh.elem(i).parent_id;
						} else {
							f[k++] = mesh.elem(i).id;
						}
						break;
					}

					default:
					{
						f[k++] = mesh.elem(i).id;
						break;
					}
				}
			}
		}

		if(plot_fun == PLOT_TAG || plot_fun == PLOT_PARENT_TAG) {
			hsv = f;
		} else {
			moonolith::func_to_hsv(f, hsv);
		}

		canvas.fill_mesh(m, hsv);
		canvas.set_color(0,0,0);
		canvas.stroke_mesh(m);

		for(std::size_t i = 0, k = 0; i < mesh.n_elements(); ++i) {
			if(mesh.is_active(i)) {
				auto b = barycenter(mesh.elem(i), mesh.points());
				canvas.draw_text(b(0)*scale_factor, b(1)*scale_factor, 4./mesh.n_active_elements(), "Arial", std::to_string(i), true);
			}
		}

		return canvas.write(path);
	}

	template<Integer Dim, Integer ManifoldDim>
	bool read_mesh(const std::string &path, Mesh<Dim, ManifoldDim> &mesh, const bool verbose = false)
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
				assert(dim == ManifoldDim);
			} else if(line == "elements") {
				std::getline(is, line);
				n_elements = atoi(line.c_str());

				for(Integer i = 0; i < n_elements; ++i) {
					assert(is.good());
					std::getline(is, line);
					std::stringstream ss(line);
					int attr, type;

					std::array<Integer, ManifoldDim+1> nodes;
					ss >> attr >> type;

					for(Integer k = 0; k < ManifoldDim+1; ++k) {
						ss >> nodes[k];
					}

					mesh.add_elem(nodes);
				}
			} else if(line == "vertices") {
				std::getline(is, line);
				n_nodes = atoi(line.c_str());
				std::getline(is, line);
				n_coords = atoi(line.c_str());
				assert(n_coords == Dim);

				Vector<Real, Dim> p;
				p.zero();
				for(Integer i = 0; i < n_nodes; ++i) {
					assert(is.good());

					for(Integer k = 0; k < n_coords; ++k) {
						is >> p(k);
					}

					mesh.add_point(p);
				}

			}
		}

		is.close();

		mesh.repair(verbose);
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
