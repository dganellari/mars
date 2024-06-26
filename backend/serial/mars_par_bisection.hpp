#ifndef MARS_PAR_BISECTION_HPP
#define MARS_PAR_BISECTION_HPP

#include "mars_communicator.hpp"
#include "mars_fwd.hpp"
#include "mars_longest_edge.hpp"
#include "mars_oldest_edge.hpp"
#include "mars_par_edge_split_pool.hpp"
#include "mars_par_mesh.hpp"

namespace mars {

    template <class Mesh, class EdgeSelect>
    class Bisection;

    template <class ParMesh>
    class ParBisection {
    public:
        template <typename T>
        using ptr = std::shared_ptr<T>;

        static const Integer Dim = ParMesh::Dim;
        static const Integer ManifoldDim = ParMesh::ManifoldDim;
        using Mesh = typename ParMesh::Mesh;

        using ES = EdgeSelect<Mesh>;
        using Bisection = mars::Bisection<Mesh, ES>;
        using DefaultES = OldestEdgeSelect<Mesh>;

        ParBisection(ParMesh &mesh) : mesh_(mesh), bisection_(mesh.get_serial_mesh()), edge_split_pool_(mesh.comm()) {
            bisection_.set_edge_select(std::make_shared<DefaultES>(mesh_.node_map()));
        }

        void set_edge_select(const std::shared_ptr<ES> &edge_select) { bisection_.set_edge_select(edge_select); }

        // 1) (parallel)
        void refine_elements(const std::vector<mars::Integer> &elements) {
            if (verbose) {
                std::cout << "------------------------------\n";
            }

            build_edge_interface();
            bisection_.refine(elements);
        }

        void build_edge_interface() {
            if (bisection_.edge_element_map().empty()) {
                bisection_.edge_element_map().update(mesh_.get_serial_mesh());
            }

            edge_split_pool_.build_edge_interface(mesh_, bisection_);
        }

        void exchange_edge_pools() { edge_split_pool_.synchronize(); }

        // 2) (synchronized)
        void update_edge_pools() {
            edge_split_pool_.update(
                mesh_, bisection_.edge_element_map(), bisection_.edge_node_map(), bisection_.bisected_edges());

            exchange_edge_pools();

            bisection_.clear_bisected_edges();
        }

        // 3) (synchronized)
        void update_global_ids() {
            edge_split_pool_.write_to_mesh(bisection_.edge_node_map(), mesh_);

            mesh_.synchronize();

            edge_split_pool_.read_from_mesh(bisection_.edge_node_map(), mesh_);

            exchange_edge_pools();

            edge_split_pool_.write_to_mesh(bisection_.edge_node_map(), mesh_);

            edge_split_pool_.update_midpoint_parts(bisection_.edge_node_map(), mesh_);

            build_edge_interface();
        }

        // 4) (parallel)
        bool conform_interfaces() {
            // refine edges and go to 2)
            std::vector<Edge> splits;
            edge_split_pool_.collect_splits_to_local_edges(mesh_, splits);

            // refine edges
            bisection_.refine_edges(splits);
            int has_more = !edge_split_pool_.empty();

            mesh_.comm().all_reduce(&has_more, 1, MPIMax());
            return !has_more;
        }

        void uniform_refine(const Integer n_levels) {
            for (Integer l = 0; l < n_levels; ++l) {
                if (verbose) {
                    std::cout << "------------------------------\n";
                }

                build_edge_interface();

                bisection_.uniform_refine(1);

                green_refinement();
            }
        }

        void green_refinement() {
            bool complete = false;

            Integer max_loops = 20;
            Integer loops = 0;

            while (!complete) {
                update_edge_pools();
                update_global_ids();
                complete = conform_interfaces();

                std::cout << "loop " << loops++ << " done. complete = " << complete << std::endl;

                if (loops >= max_loops && !complete) {
                    assert(false);
                }
            }
        }

        void refine(std::vector<std::vector<mars::Integer>> &elements) {
            refine_elements(elements);
            green_refinement();
        }

        ParMesh &mesh_;
        Bisection bisection_;
        ParEdgeSplitPool edge_split_pool_;
        bool verbose = false;
    };
}  // namespace mars

#endif  // MARS_PAR_BISECTION_HPP
