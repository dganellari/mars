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
	class MeshPartition;

	template<class Mesh>
	class VTKMeshWriter;

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
		PLOT_PARENT_TAG = 4,
		PLOT_NUMERIC_TAG = 5,
		PLOT_UNIFORM = 6
	};

	class PlotOpts {
	public:
		PlotFun plot_fun;
		Real scale_factor;
		Real node_size;
		Real uniform;
		bool active_only;
		bool show_id;

		std::vector<Integer> node_id;
		std::vector<Integer> element_id;
		std::vector<Real>    fun;

		PlotOpts()
		: plot_fun(PLOT_ID),
		  scale_factor(10.),
		  node_size(2.),
		  uniform(0.),
		  active_only(true),
		  show_id(true)
		{}
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
			case BISECTION:
			{
				b = 1.;
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
	template<Integer Dim, Integer ManifoldDim>
	void mesh_color(
		const Mesh<Dim, ManifoldDim> &mesh,
		const PlotOpts &opts,
		std::vector<Real> &hsv
		)
	{

		const bool active_only = opts.active_only;
		const Integer plot_fun = opts.plot_fun;

		std::vector<Real> f;

		if(plot_fun == PLOT_TAG || plot_fun == PLOT_PARENT_TAG) {
			if(active_only) {
				f.resize(mesh.n_active_elements() * 3);
			} else {
				f.resize(mesh.n_elements() * 3);
			}
		} else {
			if(active_only) {
				f.resize(mesh.n_active_elements());
			} else {
				f.resize(mesh.n_elements());
			}
		} 

	
		for(std::size_t i = 0, k = 0; i < mesh.n_elements(); ++i) {
			if(mesh.is_active(i) || !active_only) {
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

					case PLOT_NUMERIC_TAG:
					{
						f[k++] = mesh.tags()[i];
						break;
					}

					default:
					{
						if(!opts.element_id.empty()) {
							f[k++] = opts.element_id[mesh.elem(i).id];
						} else {
							f[k++] = mesh.elem(i).id;
						}
						break;
					}
				}
			}
		}

		if(plot_fun == PLOT_UNIFORM) {
			std::fill(std::begin(f), std::end(f), opts.uniform);
		} 

		if(plot_fun == PLOT_TAG || plot_fun == PLOT_PARENT_TAG) {
			hsv = f;
		} else if(plot_fun == PLOT_UNIFORM) {
			moonolith::func_to_hsv(f, 0., 1., hsv);
		} else {

			Real max_el = *std::max_element(std::begin(f), std::end(f));
			Real min_el = *std::min_element(std::begin(f), std::end(f));

			if(max_el != min_el) {
				moonolith::func_to_hsv(f, hsv);
			}
		}
	}

	template<class Canvas, Integer Dim>
	void draw_mesh(
		Canvas &canvas,
		const Mesh<Dim, 2> &mesh,
		const PlotOpts &opts)
	{

		moonolith::Mesh m;
		m.dim = Dim;

		m.points.resize(mesh.n_nodes() * Dim);
		m.el_index.resize(mesh.n_active_elements() * 3);

		for(Integer i = 0; i < mesh.n_nodes(); ++i) {
			for(Integer d = 0; d < Dim; ++d) {
				m.points[i * Dim + d] = mesh.point(i)(d) * opts.scale_factor;
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

		canvas.set_line_width(0.1/mesh.n_active_elements());

		std::vector<Real> hsv;
		mesh_color(mesh, opts, hsv);

		canvas.fill_mesh(m, hsv);
		canvas.set_color(0,0,0);
		canvas.stroke_mesh(m);

		for(std::size_t i = 0, k = 0; i < mesh.n_elements(); ++i) {
			if(mesh.is_active(i)) {
				for(auto n : mesh.elem(i).nodes) {
					auto p = mesh.point(n);

					auto n_id = n;

					if(!opts.node_id.empty()) {
						n_id = opts.node_id[n]; 
					}

					if(n_id == INVALID_INDEX) {
						canvas.set_color(1., 0., 0.);
						canvas.stroke_circle(p(0)*opts.scale_factor, p(1)*opts.scale_factor, opts.node_size);
					} else {
						if(opts.show_id) {
							canvas.set_color(0., 0., 0.);
							canvas.stroke_circle(p(0)*opts.scale_factor, p(1)*opts.scale_factor, opts.node_size);

							canvas.set_color(1., 1., 1.);
							canvas.fill_circle(p(0)*opts.scale_factor, p(1)*opts.scale_factor, opts.node_size);
						}
					}
					

					canvas.update_box(p(0)*opts.scale_factor + opts.scale_factor/10., p(1)*opts.scale_factor + opts.scale_factor/10.);
					canvas.update_box(p(0)*opts.scale_factor - opts.scale_factor/10., p(1)*opts.scale_factor - opts.scale_factor/10.);
				}
			}
		}

		if(opts.show_id) {
			canvas.set_color(0., 0., 0.);
			
			for(std::size_t i = 0, k = 0; i < mesh.n_elements(); ++i) {
				if(mesh.is_active(i)) {
					auto b = barycenter(mesh.elem(i), mesh.points());
					auto e_id = i;

					if(!opts.element_id.empty()) {
						e_id = opts.element_id[i]; 
					}

					canvas.draw_text(b(0)*opts.scale_factor, b(1)*opts.scale_factor, opts.node_size, "Courier", std::to_string(e_id), true);

					for(auto n : mesh.elem(i).nodes) {
						auto p = mesh.point(n);

						auto n_id = n;

						if(!opts.node_id.empty()) {
							n_id = opts.node_id[n]; 
						}

						if(n_id != INVALID_INDEX) {
							canvas.draw_text(p(0)*opts.scale_factor, p(1)*opts.scale_factor, opts.node_size, "Courier", std::to_string(n_id), true);
						}
					}
				}
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
		moonolith::EPSCanvas canvas;
		PlotOpts opts;
		opts.scale_factor = scale_factor;
		opts.plot_fun = plot_fun;
		opts.node_size = 1./std::sqrt(mesh.n_active_elements());

		draw_mesh(canvas, mesh, opts);
		return canvas.write(path);
	}

	template<Integer Dim, Integer ManifoldDim>
	bool export_parts(
		const std::string &path,
		const std::vector<std::shared_ptr<MeshPartition<Dim, ManifoldDim>>> &parts)
	{
		bool ok = true;
		for(const auto &p : parts) {
			VTKMeshWriter<Mesh<Dim, ManifoldDim>> w;
			ok &= w.write(path + std::to_string(p->partition_id()) + ".vtu", p->get_mesh());
		}

		return ok;
	}

	template<Integer Dim>
	bool write_mesh_partitions(
		const std::string &path,
		const std::vector<std::shared_ptr<MeshPartition<Dim, 3>>> &parts,
		const PlotFun plot_fun)
	{
		bool ok = true;
		for(const auto &p : parts) {
			VTKMeshWriter<Mesh<Dim, 3>> w;
			ok &= w.write(path + std::to_string(p->partition_id()) + ".vtu", p->get_mesh());
		}

		return ok;
	}


	template<Integer Dim>
	bool write_mesh_partitions(
		const std::string &path,
		const std::vector<std::shared_ptr<MeshPartition<Dim, 2>>> &parts,
		const PlotFun plot_fun)
	{
		moonolith::EPSCanvas canvas;

		Integer n_elements = 0;
		for(const auto &p : parts) {
			n_elements += p->get_mesh().n_elements();
		}

		PlotOpts opts;
		// opts.scale_factor = scale_factor;
		opts.plot_fun = plot_fun;
		opts.node_size = 1./std::sqrt(n_elements);

		if(n_elements > 400) {
			opts.show_id = false;
		}

		for(const auto &p : parts) {
			opts.node_id    = p->node_map().global_id();
			opts.element_id = p->elem_map().global_id();
			opts.uniform = Real(p->partition_id() + 1)/parts.size();
			draw_mesh(canvas, p->get_mesh(), opts);
		}
		return canvas.write(path);
	}


	template<class Canvas, Integer Dim, Integer ManifoldDim, class NodeVector>
	void draw_element_subsurface(
		const Mesh<Dim, ManifoldDim> &m,
		const Integer element_id,
		const NodeVector &sub_surface_nodes,
		const Real cx,
		const Real cy,
		const Real scale_factor,
		const RGB &rgb,
		Canvas &canvas,
		const bool sort_nodes = true)
	{

		std::cout << "rgb: " << rgb.r << " " << rgb.g << " " << rgb.b << std::endl;

		auto e = m.elem(element_id);
		if(sort_nodes) {
			std::sort(e.nodes.begin(), e.nodes.end());
		}

		const Real node_size = 0.5;
		const Real line_width_e = 0.2;
		std::vector<Real> xy(sub_surface_nodes.size()*2);
		std::vector<Integer> node_ids(sub_surface_nodes.size(), INVALID_INDEX);

		const Real dAngle = 2.*M_PI/n_nodes(e);

		Integer index = 0;
		for(Integer k = 0; k < n_nodes(e); ++k) {
			if(std::find(
				sub_surface_nodes.begin(),
				sub_surface_nodes.end(),
				e.nodes[k]) != sub_surface_nodes.end()
				) {
				xy[index * 2]     = cx + std::cos(dAngle * k) * scale_factor;
			xy[index * 2 + 1] = cy + std::sin(dAngle * k) * scale_factor;
			node_ids[index] = e.nodes[k];
			++index;
		}
	}

	canvas.set_line_width(line_width_e);
	canvas.set_color(rgb.r, rgb.g, rgb.b);
	canvas.set_dashing(1./(2. + 2*rand()/double(RAND_MAX)));
	canvas.stroke_polygon(&xy[0], sub_surface_nodes.size());
	canvas.clear_dashing();

	for(Integer k = 0; k < sub_surface_nodes.size(); ++k) {
		canvas.set_color(0., 0., 0.);
		canvas.stroke_circle(xy[k*2], xy[k*2+1], node_size*2.);
		canvas.set_color(1., 1., 1.);
		canvas.fill_circle(xy[k*2], xy[k*2+1], node_size*1.94);
	}

	canvas.set_color(0., 0., 0.);
	for(Integer k = 0; k < sub_surface_nodes.size(); ++k) {
		canvas.draw_text(xy[k*2], xy[k*2+1], 1, "Courier", std::to_string(node_ids[k]), true);
	}

		// auto b = barycenter(e);

	if(e.id != INVALID_INDEX) {
		canvas.set_color(rgb.r, rgb.g, rgb.b);
		canvas.draw_text(cx, cy - scale_factor * 1.2, 1.5, "Courier", std::to_string(e.id), true);
	}
}


	template<class Canvas, Integer Dim, Integer ManifoldDim>
void draw_element_side(
	const Mesh<Dim, ManifoldDim> &m,
	const Integer element_id,
	const Integer side_num,
	const Real cx,
	const Real cy,
	const Real scale_factor,
	const RGB &rgb,
	Canvas &canvas,
	const bool sort_nodes = true)
{
	Simplex<Dim, ManifoldDim-1> side;
	m.elem(element_id).side(side_num, side);
	draw_element_subsurface(
		m,
		element_id,
		side.nodes,
		cx, cy, scale_factor, rgb,
		canvas,
		sort_nodes);
}

	template<class Canvas, Integer Dim, Integer ManifoldDim>
void draw_element(
	const Mesh<Dim, ManifoldDim> &m,
	const Integer element_id,
	const Real cx,
	const Real cy,
	const Real scale_factor,
	const RGB &rgb,
	Canvas &canvas,
	const bool sort_nodes = true)
{
	const Real node_size = 0.5;
	const Real line_width_e = 0.04;

	auto e = m.elem(element_id);
	if(sort_nodes) {
		std::sort(e.nodes.begin(), e.nodes.end());
	}

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
	const Real node_size = 0.6;
	const Real line_width_e = 0.02;

	std::array<Real, 2*(ManifoldDim+1)> c_xy;
	std::array<Integer, (ManifoldDim+1)> c_index;
	std::array<Integer, (ManifoldDim+1)> is_midpoint;
	std::array<Edge, (ManifoldDim+1)> edge;

	auto c = m.elem(element_id);
	std::sort(c.nodes.begin(), c.nodes.end());
	auto e = m.elem(c.parent_id);
	std::sort(e.nodes.begin(), e.nodes.end());

	const auto &e_nodes = e.nodes;
	const auto &c_nodes = c.nodes;

	const Real dAngle    = 2.*M_PI/n_nodes(e);
	const Real dAngleMid = 2.*M_PI/Combinations<ManifoldDim + 1, 2>::value;

	Integer index = 0;
	Integer mid_index = 0;
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
					// c_xy[index * 2]     = cx + (std::cos(dAngle * k) + std::cos(dAngle * j)) * scale_factor/2.;
					// c_xy[index * 2 + 1] = cy + (std::sin(dAngle * k) + std::sin(dAngle * j)) * scale_factor/2.;

				c_xy[index * 2]     = cx + std::cos(dAngleMid/2 + dAngleMid * midpoint_index<ManifoldDim>(k, j)) * scale_factor/1.5;
				c_xy[index * 2 + 1] = cy + std::sin(dAngleMid/2 + dAngleMid * midpoint_index<ManifoldDim>(k, j)) * scale_factor/1.5;
				c_index[index] = v;
				is_midpoint[index] = true;
				edge[index] = Edge(e_nodes[k], e_nodes[j]);
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
		if(!is_midpoint[k]) continue;
		canvas.draw_text(
			c_xy[k*2],
			c_xy[k*2+1],
			0.5,
			"Courier",
			std::to_string(c_index[k]) + "=(" + std::to_string(edge[k].nodes[0]) + "," + std::to_string(edge[k].nodes[1]) + ")", true);
	}

	canvas.clear_dashing();
}



	template<Integer Dim, Integer ManifoldDim>
bool write_element(
	const std::string &path,
	const Mesh<Dim, ManifoldDim> &m,
	const EdgeNodeMap &en_map,
	const Integer element_id,
	const Real scale_factor,
	const Integer child_num)
{

	auto e = m.elem(element_id);
	std::sort(e.nodes.begin(), e.nodes.end());


	std::vector<Real> hsv;

	PlotOpts opts;
	opts.plot_fun = PLOT_ID;
	opts.active_only = false;

	mesh_color(m, opts, hsv);

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
	flag_color(m.tags()[element_id], rgb.r, rgb.g, rgb.b);

	std::set<Integer> neigs;

	if(!e.children.empty()) {
		if(child_num != INVALID_INDEX) {
			const auto c = e.children[child_num];
			draw_element_as_child(
				m,
				en_map,
				c,
				cx,
				cy,
				scale_factor,
				RGB(hsv[c*3], hsv[c*3+1], hsv[c*3+2]),
				canvas
				);

			for(auto a : m.dual_graph().adj(c)) {
				if(a == INVALID_INDEX || m.elem(a).parent_id == element_id) continue;
				neigs.insert(a);
			}

		} else {
			for(auto c : e.children) {
				draw_element_as_child(
					m,
					en_map,
					c,
					cx,
					cy,
					scale_factor,
					RGB(hsv[c*3], hsv[c*3+1], hsv[c*3+2]),
					canvas
					);

				for(auto a : m.dual_graph().adj(c)) {
					if(a == INVALID_INDEX || m.elem(a).parent_id == element_id) continue;
					neigs.insert(a);
				}
			}
		}
	}

	auto e_rgb = RGB(hsv[element_id*3], hsv[element_id*3+1], hsv[element_id*3+2]);

	draw_element(
		m,
		element_id,
		cx,
		cy,
		scale_factor,
		e_rgb,
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

		const auto a_cx = cx + cos(k * dAngle + M_PI) * 2.5 * scale_factor;
		const auto a_cy = cy + sin(k * dAngle + M_PI) * 2.5 * scale_factor;

		draw_element(
			m,
			adj[i],
			a_cx,
			a_cy,
			scale_factor,
			rgb,
			canvas);

		for(auto e_c : neigs) {
			if(m.is_child(adj[i], e_c)) {
				draw_element_as_child(
					m,
					en_map,
					e_c,
					a_cx,
					a_cy,
					scale_factor,
					RGB(hsv[e_c*3], hsv[e_c*3+1], hsv[e_c*3+2]),
					canvas
					);
			}
		}

	}

	return canvas.write(path);
}

	template<Integer Dim, Integer ManifoldDim>
bool write_element(
	const std::string &path,
	const RedGreenRefinement<Dim, ManifoldDim> &rgr,
	const Integer element_id,
	const Real scale_factor,
	const Integer child_num)
	{

		return write_element(
			path,
			rgr.get_mesh(),
			rgr.edge_node_map(),
			element_id,
			scale_factor,
			child_num);
	}



	template<Integer Dim, Integer ManifoldDim>
bool write_element_with_sides(
	const std::string &path,
	const Mesh<Dim, ManifoldDim> &m,
	const Integer element_id,
	const Real scale_factor,
	const Integer side_num,
	const bool skip_invalid_adj_sides = true)
{
	const auto &e = m.elem(element_id);
	auto sorted_nodes = e.nodes;
		// std::sort(sorted_nodes.begin(), sorted_nodes.end());

		// MultilevelElementMap<ManifoldDim, 2> mlem;
		// mlem.update(rgr.get_mesh());

	std::vector<Real> hsv;

	PlotOpts opts;
	opts.plot_fun = PLOT_ID;
	opts.active_only = false;

	mesh_color(m, opts, hsv);

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

	auto e_rgb = RGB(hsv[element_id*3], hsv[element_id*3+1], hsv[element_id*3+2]);

	std::set<Integer> neigs;

	Simplex<Dim, ManifoldDim-1> side;

		// for(Integer sub_manifold = ManifoldDim; sub_manifold >= 2; --sub_manifold) {
	if(side_num != INVALID_INDEX) {
		const auto &adj = m.dual_graph().adj(element_id);
		const auto a = adj[side_num];

		bool skip = a == INVALID_INDEX && skip_invalid_adj_sides;
		if(!skip) {
			auto rgb = RGB(0, 0, 0);

			if(a != INVALID_INDEX) {
				rgb.r = hsv[a*3];
				rgb.g = hsv[a*3 + 1];
				rgb.b = hsv[a*3 + 2];
			}

			draw_element_side(
				m,
				element_id,
				side_num,
				cx,
				cy,
				scale_factor,
				rgb,
				canvas,
				false);

			if(a != INVALID_INDEX) {
				e.side(side_num, side);
				Integer k = 0;
				for(; k < ManifoldDim + 1; ++k) {
					auto en = sorted_nodes[k];
					if(std::find(side.nodes.begin(), side.nodes.end(), en) == side.nodes.end()) {
						break;
					}
				}


				const auto a_cx = cx + cos(k * dAngle + M_PI) * 2.5 * scale_factor;
				const auto a_cy = cy + sin(k * dAngle + M_PI) * 2.5 * scale_factor;

				draw_element(
					m,
					a,
					a_cx,
					a_cy,
					scale_factor,
					rgb,
					canvas,
					false);

				draw_element_side(
					m,
					a,
					m.common_side_num(a, element_id),
					a_cx,
					a_cy,
					scale_factor,
					e_rgb,
					canvas,
					false
					);
			}
		}

	} else {
		for(Integer s = 0; s < n_sides(e); ++s) {

			const auto &adj = m.dual_graph().adj(element_id);
			const auto a = adj[s];

			bool skip = a == INVALID_INDEX && skip_invalid_adj_sides;
			if(!skip) {
				auto rgb = RGB(0, 0, 0);

				if(a != INVALID_INDEX) {
					rgb.r = hsv[a*3];
					rgb.g = hsv[a*3 + 1];
					rgb.b = hsv[a*3 + 2];
				}

				draw_element_side(
					m,
					element_id,
					s,
					cx,
					cy,
					scale_factor,
					rgb,
					canvas,
					false
					);

				if(adj[s] != INVALID_INDEX) {
					e.side(s, side);

					Integer k = 0;
					for(; k < ManifoldDim + 1; ++k) {
						auto en = sorted_nodes[k];
						if(std::find(side.nodes.begin(), side.nodes.end(), en) == side.nodes.end()) {
							break;
						}
					}

					const auto a_cx = cx + cos(k * dAngle + M_PI) * 2.5 * scale_factor;
					const auto a_cy = cy + sin(k * dAngle + M_PI) * 2.5 * scale_factor;

					draw_element(
						m,
						a,
						a_cx,
						a_cy,
						scale_factor,
						rgb,
						canvas,
						false);

					draw_element_side(
						m,
						a,
						m.common_side_num(a, element_id),
						a_cx,
						a_cy,
						scale_factor,
						e_rgb,
						canvas,
						false
						);
				}
			}
		}
			// }
	}



	draw_element(
		m,
		element_id,
		cx,
		cy,
		scale_factor,
		e_rgb,
		canvas,
		false);

	return canvas.write(path);
}



	template<Integer Dim, Integer ManifoldDim>
bool write_element_with_subsurfaces(
	const std::string &path,
	const Mesh<Dim, ManifoldDim> &m,
	const Integer element_id,
	const Real scale_factor)
{
	// const auto &m = rgr.get_mesh();
	const auto &e = m.elem(element_id);
	auto sorted_nodes = e.nodes;
		// std::sort(sorted_nodes.begin(), sorted_nodes.end());

	MultilevelElementMap<ManifoldDim, 2> mlem;
	mlem.update(m);

	std::vector<Real> hsv;
	PlotOpts opts;
	opts.plot_fun = PLOT_ID;
	opts.active_only = false;
	mesh_color(m, opts, hsv);

	const Real margin = 1.3;
	const Real node_size = 0.5;
	const Real line_width_e = 0.04;

	const Real canvas_w = (ManifoldDim-1) * 3 * 2. * (margin * scale_factor);
	const Real canvas_h = (ManifoldDim-1) * 3 * 2. * (margin * scale_factor);

	moonolith::EPSCanvas canvas;

	canvas.update_box(-canvas_w/2, -canvas_h/2);
	canvas.update_box(canvas_w/2, canvas_h/2);

	std::array<Real, 2*(ManifoldDim+1)> xy;

	const Real dAngle = 2.*M_PI/n_nodes(e);
	const Real cx = 0.;
	const Real cy = 0.;

	auto e_rgb = RGB(hsv[element_id*3], hsv[element_id*3+1], hsv[element_id*3+2]);

	Integer n_circle = 1;
	for(Integer sub_manifold = ManifoldDim; sub_manifold >= 2; --sub_manifold, ++n_circle) {
		std::vector<Integer> neigs;
		
		mlem.adj(
			sub_manifold,
			e,
			neigs
		);

		if(neigs.empty()) continue;

		const auto dNeighAngle = 2.*M_PI/neigs.size();

		for(Integer s = 0; s < neigs.size(); ++s) {
			const auto a = neigs[s];
			if(a == element_id) continue;

			const auto &e_neigh = m.elem(a);

			auto rgb = RGB(hsv[a*3], hsv[a*3 + 1], hsv[a*3 + 2]);
				
			Integer k = s;
			// Integer k = 0;
			// for(; k < ManifoldDim + 1; ++k) {
			// 	auto en = sorted_nodes[k];
			// 	if(std::find(e_neigh.nodes.begin(), e_neigh.nodes.end(), en) == e_neigh.nodes.end()) {
			// 		break;
			// 	}
			// }

			const auto a_cx = cx + cos(k * dNeighAngle + M_PI) * 2.5 * scale_factor * n_circle;
			const auto a_cy = cy + sin(k * dNeighAngle + M_PI) * 2.5 * scale_factor * n_circle;

			draw_element(
				m,
				a,
				a_cx,
				a_cy,
				scale_factor,
				rgb,
				canvas,
				false);

		}
	}



	draw_element(
		m,
		element_id,
		cx,
		cy,
		scale_factor,
		e_rgb,
		canvas,
		false);

	return canvas.write(path);
}
}

#endif //MARS_VISUALIZATION_HPP
