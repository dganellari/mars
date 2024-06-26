#ifndef MARS_LEPP_BENCHMARK_KOKKOS_HPP_
#define MARS_LEPP_BENCHMARK_KOKKOS_HPP_

#include "mars_base.hpp"

#include <chrono>

#include "mars_bisection_kokkos.hpp"
#include "mars_longest_edge_kokkos.hpp"
#include "mars_mark_kokkos.hpp"
#include "mars_sub_view.hpp"
#include "mars_utils_kokkos.hpp"

namespace mars {
    template <class Mesh>
    class ParallelLeppBenchmark {
    public:
        using EdgeSelectPtr = std::shared_ptr<EdgeSelect<Mesh>>;

        static constexpr Integer Dim = Mesh::Dim;
        static constexpr Integer ManifoldDim = Mesh::ManifoldDim;
        using Elem = mars::Simplex<Dim, ManifoldDim, KokkosImplementation>;

        void run(const Integer n_levels, const Mesh &mesh_in, const std::string &output_path) {
            auto mesh = mesh_in;

            ParallelBisection<Mesh> b(&mesh);
            /* b.set_verbose(true); */

            b.uniform_refine(1);

            Integer exp_num = 0;
            run_benchmark(n_levels, mesh, output_path, exp_num++);
        }

        void run_benchmark(const Integer n_levels, Mesh &mesh, const std::string &output_path, const Integer exp_num) {
            using namespace mars;
            using namespace Kokkos;

            std::cout << "======================================\n";

            /*	Quality<Mesh> q(mesh);
                    q.compute();

                    mark_boundary(mesh);

                    std::cout << mesh.n_boundary_sides() << std::endl;
                    std::cout << "volume: " << mesh.volume() << std::endl; */

            ParallelBisection<Mesh> b(&mesh);
            /* b.set_verbose(true); */

            b.uniform_refine(1);

            ViewVectorTypeC<Real, Dim> center("center");
            fill_view(center, 0.5);

            Kokkos::Timer timer_refine;

            for (Integer i = 0; i < n_levels; ++i) {
                std::cout << "======================================\n";

                ViewVectorType<Integer> elements =
                    mark_hypersphere_for_refinement<Mesh>(b.get_mesh(), center, 0.5, b.get_host_mesh()->n_elements());

                std::cout << "n_marked(" << (i + 1) << "/" << n_levels << ") : " << elements.extent(0) << std::endl;

                b.refine(elements);
                // q.compute();
            }

            double time = timer_refine.seconds();
            std::cout << "parallel LEPP Refinement took: " << time << " seconds." << std::endl;

            // typename Mesh::SerialMesh sMesh;
            // convert_parallel_mesh_to_serial(sMesh, mesh);

            std::cout << "======================================\n";

            // std::cout << "volume: " << sMesh.volume() << std::endl;
            // std::cout << "n_active_elements: " << sMesh.n_active_elements() << std::endl;

            /* const Integer meshsize = mesh.get_view_elements().extent(0);
            Kokkos::parallel_for("p", meshsize,
                    MARS_LAMBDA(const Integer i){
                    [>Elem e = mesh.elem_with_children(i);<]
                    [>SubView<Integer, 2> children = e.children;<]
                    [>SubView<Integer, 3> nodes = e.nodes;<]
                    Integer children[2];

                        [>printf("eid: %li, nodes: %li, %li, %li\n", i, nodes(0), nodes(1), nodes(2));<]
                    [>if(children.is_valid())<]
                    if(mesh.get_children(children, i))
                    {
                        [>printf("eid: %li, children: %li, %li\n", i, children(0), children(1));<]
                        printf("eid: %li, children: %li, %li\n", i, children[0], children[1]);
                        }
                        }); */
            /*
                                    if(n_levels <= 20 && mesh.ManifoldDim <4){
                                            VTKMeshWriter<typename Mesh::SerialMesh> w;
                                            w.write("Parallel_LEPP_" + std::to_string(n_levels) + "_" +
               std::to_string(mesh.ManifoldDim) +".vtu", sMesh);
                                    }
            */
            /*

                                    mesh.update_dual_graph();
                                    // q.report.normalize_data_points();

                                    std::string suffix;
                                    if(edge_select->is_recursive()) {
                                            suffix = "_Recursive";
                                    }

                                    q.save_csv(edge_select->name() + suffix, output_path + "/Quality_" +
               std::to_string(exp_num) + "_" + edge_select->name() + suffix + "_lepp.csv", true); //exp_num == 0
                                    q.save_report(output_path + "/Quality_" + edge_select->name() + suffix +
               "_lepp.svg");
            */

            std::cout << "======================================\n";
        }
    };
}  // namespace mars

#endif  // MARS_BENCHMARK_HPP
