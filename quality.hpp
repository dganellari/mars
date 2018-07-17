#ifndef MARS_QUALITY_HPP
#define MARS_QUALITY_HPP


namespace mars {

	//from https://www.sciencedirect.com/science/article/pii/0168874X94900337
	template<Integer Dim, Integer ManifoldDim>
	class Quality {
	public:
		class Utils {
		public:
			static Real surface_area(
				const Simplex<Dim, ManifoldDim> &s,
			 	const std::vector<Vector<Real, Dim>> &points)
			{
				Simplex<Dim, ManifoldDim-1> side; 

				Real ret = 0.;
				for(Integer i = 0; i < n_sides(s); ++i) {
					s.side(i, side);
					ret += unsigned_volume(side, points);
				}

				return ret;
			}


			static Real avg_edge_len(
				const Simplex<Dim, ManifoldDim> &s,
			 	const std::vector<Vector<Real, Dim>> &points)
			{
				Integer v1, v2;
				Real ret = 0.;
				
				for(Integer i = 0; i < n_edges(s); ++i) {
					s.edge(i, v1, v2);
					ret += (points[v1] - points[v2]).norm();
				}

				return ret/n_edges(s);
			}

			static Real rms_edge_len(
				const Simplex<Dim, ManifoldDim> &s,
			 	const std::vector<Vector<Real, Dim>> &points)
			{
				Integer v1, v2;
				Real ret = 0.;
				
				for(Integer i = 0; i < n_edges(s); ++i) {
					s.edge(i, v1, v2);
					auto len = (points[v1] - points[v2]).norm();
					ret += len*len;
				}

				return std::sqrt(ret/n_edges(s));
			}

			static Real max_edge_len(
				const Simplex<Dim, ManifoldDim> &s,
			 	const std::vector<Vector<Real, Dim>> &points)
			{
				Integer v1, v2;
				Real ret = -std::numeric_limits<Real>::max();
				
				for(Integer i = 0; i < n_edges(s); ++i) {
					s.edge(i, v1, v2);
					ret = std::max((points[v1] - points[v2]).norm(), ret);
				}

				return ret;
			}

			static Real min_edge_len(
				const Simplex<Dim, ManifoldDim> &s,
			 	const std::vector<Vector<Real, Dim>> &points)
			{
				Integer v1, v2;
				Real ret = std::numeric_limits<Real>::max();
				
				for(Integer i = 0; i < n_edges(s); ++i) {
					s.edge(i, v1, v2);
					ret = std::min((points[v1] - points[v2]).norm(), ret);
				}

				return ret;
			}

			static Real max_side_area(
				const Simplex<Dim, ManifoldDim> &s,
			 	const std::vector<Vector<Real, Dim>> &points)
			{
				Simplex<Dim, ManifoldDim-1> side; 

				Real ret = -std::numeric_limits<Real>::max();
				for(Integer i = 0; i < n_sides(s); ++i) {
					s.side(i, side);
					ret = std::max(unsigned_volume(side, points)), ret;
				}

				return ret;
			}

			static Real min_side_area(
				const Simplex<Dim, ManifoldDim> &s,
			 	const std::vector<Vector<Real, Dim>> &points)
			{
				Simplex<Dim, ManifoldDim-1> side; 

				Real ret = std::numeric_limits<Real>::max();
				for(Integer i = 0; i < n_sides(s); ++i) {
					s.side(i, side);
					ret = std::min(volume(side, points), ret);
				}

				return ret;
			}

			static Real avg_side_area(
				const Simplex<Dim, ManifoldDim> &s,
			 	const std::vector<Vector<Real, Dim>> &points)
			{
				return surface_area(s, points)/n_sides(s);
			}

		};

		class Metric {
		public:
			virtual ~Metric() {}
			
			virtual void init(const Mesh<Dim, ManifoldDim> &m)
			{
				q_.resize(m.n_elements());
			}

			virtual void compute_and_store(
				const Simplex<Dim, ManifoldDim> &s,
			 	const std::vector<Vector<Real, Dim>> &points) 
			{
				q_[s.id] = compute(s, points);
			}

			virtual Real compute(
				const Simplex<Dim, ManifoldDim> &s,
			 	const std::vector<Vector<Real, Dim>> &points) = 0;

			virtual std::string name() const = 0;
			// virtual Real optimum() const = 0;

			virtual void finalize(const Mesh<Dim, ManifoldDim> &m)
			{
				q_avg = 0;
				q_min = std::numeric_limits<Real>::max();
				q_max = -std::numeric_limits<Real>::max();

				for(Integer i = 0; i < m.n_elements(); ++i) {
					if(!m.is_active(i)) continue;

					auto q = q_[i];

					q_avg += q;
					q_min = std::min(q_min, q);
					q_max = std::max(q_max, q);
				}

				q_avg /= m.n_active_elements();
			}

			static void describe_header(std::ostream &os)
			{
				os << "metric\t";
				os << "avg\t"; 
				os << "min\t"; 
				os << "max\n"; 
			}

			void describe(std::ostream &os) const
			{
				os << name() << "\t";
				os << q_avg << "\t";
				os << q_min << "\t";
				os << q_max << "\n";
			}

			Real q_avg;
			Real q_min;
			Real q_max;

			std::vector<Real> q_;
		};


		class Alpha : public Metric {
		public:
			 Real compute(
			 	const Simplex<Dim, ManifoldDim> &s,
			 	const std::vector<Vector<Real, Dim>> &points) override
			 {
			 	auto S_avg = Utils::avg_edge_len(s, points);
			 	auto V     = unsigned_volume(s, points);

			 	return pow(S_avg, ManifoldDim)/V;
			 }

			 std::string name() const override
			 {
			 	return "Alpha";
			 }

			 // inline Real optimum() const override {
			 // 	return 8.48528;
			 // }

		};

		class Gamma : public Metric {
		public:
			 Real compute(
			 	const Simplex<Dim, ManifoldDim> &s,
			 	const std::vector<Vector<Real, Dim>> &points) override
			 {
			 	auto S_rms = Utils::rms_edge_len(s, points);
			 	auto V     = unsigned_volume(s, points);

			 	return pow(S_rms, ManifoldDim)/V;
			 }

			 std::string name() const override
			 {
			 	return "Gamma";
			 }

			 // inline Real optimum() const override {
			 // 	return 8.48528;
			 // }

		};

		class Tau : public Metric {
		public:
			 Real compute(
			 	const Simplex<Dim, ManifoldDim> &s,
			 	const std::vector<Vector<Real, Dim>> &points) override
			 {
			 	auto S_max = Utils::max_edge_len(s, points);
			 	auto S_min = Utils::min_edge_len(s, points);

			 	return S_max/S_min;
			 }

			 std::string name() const override
			 {
			 	return "Tau";
			 }

			 // inline Real optimum() const override {
			 // 	return 1.;
			 // }

		};


		class Report {
		public:
			// class Point {
			// public:
			// 	Integer n_elements;
			// 	Real q_measure;

			// 	Point(const Integer n_elements, const Real q_measure)
			// 	: n_elements(n_elements), q_measure(q_measure)
			// 	{}
			// };

			// typedef std::shared_ptr<std::vector<Point>> VecPtr;

				typedef std::shared_ptr<std::vector<Real>> VecPtr;

			std::map<std::string, VecPtr> x;
			std::map<std::string, VecPtr> y;

			void normalize_data_points()
			{
				for(auto &yi : y) {
					auto &v = *yi.second;
					auto first = v.front();

					for(auto &vi : v) {
						vi /= first;
					}
				}
			}	


			void add_data_point(
				const std::string &name,
				const Integer n_elements,
				const Real value)
			{
				auto &x_ptr = x[name];
				auto &y_ptr = y[name];

				if(!x_ptr) {
					x_ptr = std::make_shared< std::vector<Real> >();
				}

				if(!y_ptr) {
					y_ptr = std::make_shared< std::vector<Real> >();
				}

				x_ptr->push_back(n_elements);
				y_ptr->push_back(value);
			}

			bool save_csv(
				const std::string &name,
				const std::string &path,
				const bool print_header = true) const
			{
				std::ofstream os;
				os.open(path);
				if(!os.good()) {
					return false;
				}

				if(print_header) {
					os << "name,type,x,y\n";
				}

				Integer i = 0;
				for(auto xi : x) {
					auto yi_iter = y.find(xi.first);
					assert(yi_iter != y.end());

					for(Integer j = 0; j < xi.second->size(); ++j) {
						const auto x_val = xi.second->at(j);
						const auto y_val = yi_iter->second->at(j);

						os << name << "," << xi.first << "," << x_val << "," << y_val << "\n";
					}
				}

				os.close();
				return true;
			}

			bool save(const std::string &path) const
			{
#ifdef WITH_MOONOLITH
				using namespace moonolith;

				const Integer n_data = x.size();
				std::vector<plot_data> data(n_data);

				RGB colors[] = {
					RGB(1., 0., 0.),
					RGB(1., 0.5, 0.),
					RGB(1., 0.5, 0.5),
					RGB(1., 0.5, 1.),
					RGB(1., 1., 0.5),

					RGB(0., 1., 0.),
					RGB(0., 0.5, 0.),
					RGB(0.5, 1., 1.),
					
					RGB(0.5, 0.5, 1.),
					RGB(0., 0.5, 1.),
					RGB(0., 0., 1.)
				};

				Real min_y_val = std::numeric_limits<Real>::max();
				Real max_y_val = -std::numeric_limits<Real>::max();

				Real min_x_val = std::numeric_limits<Real>::max();
				Real max_x_val = -std::numeric_limits<Real>::max();

				Integer i = 0;
				for(auto xi : x) {
					auto yi_iter = y.find(xi.first);
					assert(yi_iter != y.end());

					auto yi = *yi_iter;

					data[i].x = &(*xi.second)[0];
					data[i].y = &(*yi.second)[0];
					data[i].name = xi.first;
					data[i].n_data = xi.second->size();

					data[i].line_color[0] = colors[i].r;
					data[i].line_color[1] = colors[i].g;
					data[i].line_color[2] = colors[i].b;

					data[i].show_marker = true;
					data[i].marker_size = 2;

					min_y_val = std::min(min_y_val, *std::min_element(yi.second->begin(), yi.second->end()));
					max_y_val = std::max(max_y_val, *std::max_element(yi.second->begin(), yi.second->end()));

					min_x_val = std::min(min_x_val, *std::min_element(xi.second->begin(), xi.second->end()));
					max_x_val = std::max(max_x_val, *std::max_element(xi.second->begin(), xi.second->end()));

					++i;
				}			

				axis_data axis;
				axis.axis[0]=min_x_val;
				axis.axis[1]=max_x_val;

				axis.axis[2]=min_y_val;
				axis.axis[3]=max_y_val;
				axis.use_integral_values_x = true;

				axis.init_ticks(5, 5);
				// axis.log_x = true;
				axis.log_y = true;

				legend legend;
				legend.show = true;
				legend.text_width = 40;
				legend.position = legend::NORTH;
				legend.font_size = 6;

				plot_options opts;

				opts.fit_to_size(axis, 150, 150);

				SVGCanvas canvas;
				plot(axis, legend, &data[0], n_data, canvas, opts);

				return canvas.write(path);
#else
				return false;
#endif //WITH_MOONOLITH
			}

		};

		Quality(const Mesh<Dim, ManifoldDim> &mesh)
		: mesh(mesh)
		{
			init();
		}

		void init()
		{
			measures.push_back(std::make_shared<Alpha>());
			// measures.push_back(std::make_shared<Gamma>());
			measures.push_back(std::make_shared<Tau>());
		}

		void compute_and_describe(std::ostream &os)
		{
			compute();

			Metric::describe_header(os);
			for(auto m : measures) {
				m->describe(os);
			}
		}

		void compute()
		{
			for(auto m : measures) {
				m->init(mesh);
			}

			for(Integer i = 0; i < mesh.n_elements(); ++i) {
				if(!mesh.is_active(i)) continue;

				for(auto m : measures) {
					m->compute_and_store(mesh.elem(i), mesh.points());
				}
			}

			for(auto m : measures) {
				m->finalize(mesh);
			}

			Integer na = mesh.n_active_elements();
			for(auto m : measures) {
				report.add_data_point(
					m->name() + "max",
					na,
					m->q_max);

				report.add_data_point(
					m->name() + "min",
					na,
					m->q_min);

				report.add_data_point(
					m->name() + "avg",
					na,
					m->q_avg);
			}
		}

		bool save_csv(const std::string &name, const std::string &path, const bool print_header = true)
		{
			return report.save_csv(name, path, print_header);
		}

		bool save_report(const std::string &path)
		{
			return report.save(path);
		}

		const Mesh<Dim, ManifoldDim> &mesh;
		std::vector< std::shared_ptr<Metric> > measures;
		Report report;


	};
}

#endif //MARS_QUALITY_HPP
