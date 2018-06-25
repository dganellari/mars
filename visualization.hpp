#ifndef MARS_VISUALIZATION_HPP
#define MARS_VISUALIZATION_HPP

#include "base.hpp"


namespace moonolith {
	using Integer = mars::Integer;
}

#include "moonolith_config.hpp"
#include "moonolith_mesh.hpp"
// #include "moonolith_svg_canvas.hpp"
#include "moonolith_eps_canvas.hpp"
#include "moonolith_plotter.hpp"
#include "moonolith_func_to_color.hpp"

namespace mars {
	template<Integer Dim, Integer ManifoldDim>
	class Mesh;

	template<Integer Dim, Integer ManifoldDim>
	class RedGreenRefinement;

	class RGB {
	public:
		Real r, g, b;
		RGB(const Real r = 0., const Real g = 0., const Real b = 0.) : r(r), g(g), b(b) {}
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
				for(auto n : mesh.elem(i).nodes) {
					auto p = mesh.point(n);
					canvas.set_color(1., 1., 1.);
					canvas.fill_circle(p(0)*scale_factor, p(1)*scale_factor, 1./std::sqrt(mesh.n_active_elements()));
					canvas.set_color(0., 0., 0.);
					canvas.stroke_circle(p(0)*scale_factor, p(1)*scale_factor, 1./std::sqrt(mesh.n_active_elements()));

					canvas.update_box(p(0)*scale_factor + 4./mesh.n_active_elements(), p(1)*scale_factor + scale_factor/10.);
					canvas.update_box(p(0)*scale_factor - 4./mesh.n_active_elements(), p(1)*scale_factor - scale_factor/10.);
				}
			}
		}


		for(std::size_t i = 0, k = 0; i < mesh.n_elements(); ++i) {
			if(mesh.is_active(i)) {
				auto b = barycenter(mesh.elem(i), mesh.points());
				canvas.draw_text(b(0)*scale_factor, b(1)*scale_factor, 1./std::sqrt(mesh.n_active_elements()), "Courier", std::to_string(i), true);

				for(auto n : mesh.elem(i).nodes) {
					auto p = mesh.point(n);
					canvas.draw_text(p(0)*scale_factor, p(1)*scale_factor, 1./std::sqrt(mesh.n_active_elements()), "Courier", std::to_string(n), true);
				}
			}
		}

		return canvas.write(path);
	}

	template<class Canvas, Integer Dim, Integer ManifoldDim>
	void draw_element(
		const Mesh<Dim, ManifoldDim> &m,
		const Integer element_id,
		const Real cx,
		const Real cy,
		const Real scale_factor,
		const RGB &rgb,
		Canvas &canvas)
	{
		const Real node_size = 0.5;
		const Real line_width_e = 0.04;

		const auto &e = m.elem(element_id);
		std::array<Real, 2*(ManifoldDim+1)> xy;

		const Real dAngle = 2.*M_PI/n_nodes(e);

		for(Integer k = 0; k < n_nodes(e); ++k) {
			xy[k * 2]     = cx + std::cos(dAngle * k) * scale_factor;
			xy[k * 2 + 1] = cy + std::sin(dAngle * k) * scale_factor;
		}

		canvas.set_line_width(line_width_e);
		canvas.set_color(rgb.r, rgb.g, rgb.b);
		canvas.stroke_polygon(&xy[0], n_nodes(e));

		for(Integer k = 0; k < n_nodes(e); ++k) {
			canvas.set_color(0., 0., 0.);
			canvas.stroke_circle(xy[k*2], xy[k*2+1], node_size*2.);
			canvas.set_color(1., 1., 1.);
			canvas.fill_circle(xy[k*2], xy[k*2+1], node_size*1.94);
		}

		canvas.set_color(0., 0., 0.);
		for(Integer k = 0; k < n_nodes(e); ++k) {
			canvas.draw_text(xy[k*2], xy[k*2+1], 1, "Courier", std::to_string(e.nodes[k]), true);
		}

		// auto b = barycenter(e);
		canvas.set_color(rgb.r, rgb.g, rgb.b);
		canvas.draw_text(cx, cy - scale_factor * 1.2, 1.5, "Courier", std::to_string(element_id), true);
	}

	template<Integer Dim, Integer ManifoldDim, class Canvas>
	void draw_element_as_child(
		const Mesh<Dim, ManifoldDim> &m,
		const EdgeNodeMap &en_map,
		const Integer element_id,
		const Real cx,
		const Real cy,
		const Real scale_factor,
		const RGB &rgb,
		Canvas &canvas)
	{
		const Real node_size = 0.5;
		const Real line_width_e = 0.04;

		std::array<Real, 2*(ManifoldDim+1)> c_xy;
		std::array<Integer, (ManifoldDim+1)> c_index;
		std::array<Integer, (ManifoldDim+1)> is_midpoint;

		const auto &c = m.elem(element_id);
		const auto &e = m.elem(c.parent_id);

		const auto &e_nodes = e.nodes;
		const auto &c_nodes = c.nodes;

		const Real dAngle = 2.*M_PI/n_nodes(e);

		Integer index = 0;
		for(Integer k = 0; k < ManifoldDim + 1; ++k) {
			auto it = std::find(c_nodes.begin(), c_nodes.end(), e_nodes[k]);

			if(it != c_nodes.end()) {
				c_xy[index * 2]     = cx + std::cos(dAngle * k) * scale_factor;
				c_xy[index * 2 + 1] = cy + std::sin(dAngle * k) * scale_factor;
				c_index[index] = e_nodes[k];
				is_midpoint[index] = false;
				++index;
			}

			for(Integer j = k + 1; j < ManifoldDim + 1; ++j) {
				auto v = en_map.get(e_nodes[k], e_nodes[j]);
				if(v == INVALID_INDEX) continue;
				auto it = std::find(c_nodes.begin(), c_nodes.end(), v);

				if(it != c_nodes.end()) {
					c_xy[index * 2]     = cx + (std::cos(dAngle * k) + std::cos(dAngle * j)) * scale_factor/2.;
					c_xy[index * 2 + 1] = cy + (std::sin(dAngle * k) + std::sin(dAngle * j)) * scale_factor/2.;
					c_index[index] = v;
					is_midpoint[index] = true;
					++index;
				}
			}

			assert(index <= ManifoldDim+1);
		}

		canvas.set_line_width(0.1);
		canvas.set_color(
			rgb.r, rgb.g, rgb.b
		);

		for(Integer i = 0; i < ManifoldDim + 1; ++i) {
			for(Integer j = i + 1; j < ManifoldDim + 1; ++j) {
				canvas.stroke_line(c_xy[i*2], c_xy[i*2+1], c_xy[j*2], c_xy[j*2+1]);
			}
		}

		canvas.set_dashing(0.5);
		
		for(Integer k = 0; k < n_nodes(c); ++k) {
			if(!is_midpoint[k]) continue;
			canvas.set_color(0.5, 0.5, 0.5);
			canvas.stroke_circle(c_xy[k*2], c_xy[k*2+1], node_size*2.);
			canvas.set_color(1., 1., 1.);
			canvas.fill_circle(c_xy[k*2], c_xy[k*2+1], node_size*1.94);
		}

		canvas.set_line_width(line_width_e);

		canvas.set_color(0,0,0);
		for(Integer k = 0; k < n_nodes(c); ++k) {
			canvas.draw_text(c_xy[k*2], c_xy[k*2+1], 1, "Courier", std::to_string(c_index[k]), true);
		}

		canvas.clear_dashing();
	}

	template<Integer Dim, Integer ManifoldDim>
	bool write_element(
		const std::string &path,
		const RedGreenRefinement<Dim, ManifoldDim> &rgr,
		const Integer element_id,
		const Real scale_factor,
		const Integer child_num)
	{
		const auto &m      = rgr.get_mesh();
		const auto &e      = m.elem(element_id);
		const auto &en_map = rgr.edge_node_map();

		const Real margin = 1.3;
		const Real node_size = 0.5;
		const Real line_width_e = 0.04;

		const Real canvas_w = 3 * 2. * (margin * scale_factor);
		const Real canvas_h = 3 * 2. * (margin * scale_factor);

		moonolith::EPSCanvas canvas;
		
		canvas.update_box(-canvas_w/2, -canvas_h/2);
		canvas.update_box(canvas_w/2, canvas_h/2);

		std::array<Real, 2*(ManifoldDim+1)> xy;

		const Real dAngle = 2.*M_PI/n_nodes(e);
		const Real cx = 0.;
		const Real cy = 0.;

		RGB rgb;
		flag_color(rgr.refinement_flag(element_id), rgb.r, rgb.g, rgb.b);

		if(child_num != INVALID_INDEX && !e.children.empty()) {
			draw_element_as_child(
					m,
					en_map,
					e.children[child_num],
					cx,
					cy,
					scale_factor,
					RGB(0., 0., 1.),
					canvas
			);
		}

		draw_element(
			m,
			element_id,
			cx,
			cy,
			scale_factor,
			rgb,
			canvas);



		auto &adj = m.dual_graph().adj(element_id);

		for(Integer i = 0; i < adj.size(); ++i) {
			if(adj[i] == INVALID_INDEX) continue;
			const auto &e_adj = m.elem(adj[i]);

			
			Integer k = 0;
			for(; k < ManifoldDim + 1; ++k) {
				auto en = e.nodes[k];
				if(std::find(e_adj.nodes.begin(), e_adj.nodes.end(), en) == e_adj.nodes.end()) {
					break;
				}
			}

			draw_element(
				m,
				adj[i],
				cx + cos(k * dAngle + M_PI) * 2.5 * scale_factor,
				cy + sin(k * dAngle + M_PI) * 2.5 * scale_factor ,
				scale_factor,
				rgb,
				canvas);

		}

		return canvas.write(path);
	}
}

#endif //MARS_VISUALIZATION_HPP
