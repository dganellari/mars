#ifndef MARS_BENCHMARK_HPP
#define MARS_BENCHMARK_HPP

#include "mars_newest_vertex.hpp"
#include "mars_longest_edge.hpp"
#include "mars_bisection.hpp"
#include "mars_quality.hpp"
#include <chrono>

namespace mars {
	template<class Mesh>
	class Benchmark {
	public:
		using EdgeSelectPtr = std::shared_ptr<EdgeSelect<Mesh>>;
		
		static const Integer Dim 		 = Mesh::Dim;
		static const Integer ManifoldDim = Mesh::ManifoldDim;

		void run(
			const Integer n_levels,
			const Mesh &mesh_in,
			const std::string &output_path)
		{
			std::vector<EdgeSelectPtr> edge_selects;

			//recursive 
			edge_selects.push_back(std::make_shared<LongestEdgeSelect<Mesh>>());
			// edge_selects.push_back(std::make_shared<NewestVertexEdgeSelect<Mesh>>());
			// edge_selects.push_back(std::make_shared<NewestVertexAndLongestEdgeSelect<Mesh>>());

			//non-recursive
			// edge_selects.push_back(std::make_shared<LongestEdgeSelect<Mesh>>(false));
			// edge_selects.push_back(std::make_shared<NewestVertexEdgeSelect<Mesh>>(false));
			// edge_selects.push_back(std::make_shared<NewestVertexAndLongestEdgeSelect<Mesh>>(false));

			//refine once for creating nice intial set-up for newest vertex algorithm
			auto mesh = mesh_in;
			Bisection<Mesh> b(mesh);
			b.uniform_refine(1);

			Integer exp_num = 0;
			for(auto es : edge_selects) {
				run_benchmark(n_levels, es, mesh, output_path, exp_num++);
			}
		}

		void run_benchmark(
			const Integer n_levels,
			const EdgeSelectPtr &edge_select,
			const Mesh &mesh_in,
			const std::string &output_path,
			const Integer exp_num
			)
		{
			using namespace mars;
			std::cout << "======================================\n";
			
			//copy mesh
			auto mesh = mesh_in;

			Quality<Mesh> q(mesh);
			q.compute();

			mark_boundary(mesh);

			std::cout << mesh.n_boundary_sides() << std::endl;
			std::cout << "volume: " << mesh.volume() << std::endl;
			std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;

			Bisection<Mesh> b(mesh);
			b.set_edge_select(edge_select);
			b.uniform_refine(2);

			std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

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

				std::cout << "n_marked(" << (i+1) << "/" << n_levels << ") : " << elements.size() << std::endl;

				b.refine(elements);
				q.compute();
			}

			std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

			auto duration = std::chrono::duration_cast < std::chrono::seconds > (t2 - t1).count();
			std::cout << "LE Refinement took: "<< duration<<" seconds."<<std::endl;

			std::cout << "volume: " << mesh.volume() << std::endl;
			std::cout << "n_active_elements: " << mesh.n_active_elements() << std::endl;

			if(n_levels <= 20){
				VTKMeshWriter<Mesh> w;
				w.write(edge_select->name() + std::to_string(n_levels) + ".vtu", mesh);
			}

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
