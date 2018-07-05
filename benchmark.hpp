#ifndef MARS_BENCHMARK_HPP
#define MARS_BENCHMARK_HPP

namespace mars {
	template<Integer Dim, Integer ManifoldDim>
	class Benchmark {
	public:
		typedef std::shared_ptr<EdgeSelect<Dim, ManifoldDim>> EdgeSelectPtr;

		void run(
			const Integer n_levels,
			const Mesh<Dim, ManifoldDim> &mesh_in,
			const std::string &output_path)
		{
			std::vector<EdgeSelectPtr> edge_selects;

			//recursive 
			edge_selects.push_back(std::make_shared<LongestEdgeSelect<Dim, ManifoldDim>>());
			edge_selects.push_back(std::make_shared<NewestVertexEdgeSelect<Dim, ManifoldDim>>());
			edge_selects.push_back(std::make_shared<NewestVertexAndLongestEdgeSelect<Dim, ManifoldDim>>());

			//non-recursive
			// edge_selects.push_back(std::make_shared<LongestEdgeSelect<Dim, ManifoldDim>>(false));
			// edge_selects.push_back(std::make_shared<NewestVertexEdgeSelect<Dim, ManifoldDim>>(false));
			// edge_selects.push_back(std::make_shared<NewestVertexAndLongestEdgeSelect<Dim, ManifoldDim>>(false));

			//refine once for creating nice intial set-up for newest vertex algorithm
			auto mesh = mesh_in;
			Bisection<Dim, ManifoldDim> b(mesh);
			b.uniform_refine(1);

			Integer exp_num = 0;
			for(auto es : edge_selects) {
				run_benchmark(n_levels, es, mesh, output_path, exp_num++);
			}
		}

		void run_benchmark(
			const Integer n_levels,
			const EdgeSelectPtr &edge_select,
			const Mesh<Dim, ManifoldDim> &mesh_in,
			const std::string &output_path,
			const Integer exp_num
			)
		{
			using namespace mars;
			std::cout << "======================================\n";
			
			//copy mesh
			auto mesh = mesh_in;

			Quality<Dim, ManifoldDim> q(mesh);
			q.compute();

			mark_boundary(mesh);

			std::cout << mesh.n_boundary_sides() << std::endl;
			std::cout << "volume: " << mesh.volume() << std::endl;
			std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;

			Bisection<Dim, ManifoldDim> b(mesh);
			b.set_edge_select(edge_select);
			b.uniform_refine(2);

			for(Integer i = 0; i < n_levels; ++i) {
				std::vector<mars::Integer> elements;
				
				Vector<Real, Dim> center;
				center.set(0.5);

				mark_hypersphere_for_refinement(
					mesh,
					center,
					0.25,
					elements
					);

				std::cout << "n_marked(" << i << "/" << n_levels << ") : " << elements.size() << std::endl;

				b.refine(elements);
				q.compute();
			}

			std::cout << "volume: " << mesh.volume() << std::endl;
			std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;

			mesh.update_dual_graph();
			// q.report.normalize_data_points();

			std::string suffix;
			if(edge_select->is_recursive()) {
				suffix = "_Recursive";
			}

			q.save_csv(edge_select->name() + suffix, output_path + "/Quality_" + std::to_string(exp_num) + "_" + edge_select->name() + suffix + ".csv", true); //exp_num == 0
			q.save_report(output_path + "/Quality_" + edge_select->name() + suffix + ".svg");

			std::cout << "======================================\n";
		}
	};
}

#endif //MARS_BENCHMARK_HPP
