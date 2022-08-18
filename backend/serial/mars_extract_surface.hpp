#ifndef MARS_EXTRACT_SURFACE_HPP
#define MARS_EXTRACT_SURFACE_HPP

#include <numeric>
#include <vector>
#include "mars_mesh.hpp"

namespace mars {

    template <Integer Dim, Integer ManifoldDim>
    bool extract_surface(const mars::Mesh<Dim, ManifoldDim> &vol, mars::Mesh<Dim, ManifoldDim - 1> &surf) {
        using namespace mars;

        assert(!vol.dual_graph().empty());

        if (vol.dual_graph().empty()) {
            std::cerr << "[Error] extract_surface(,) requires the dual_graph to be built in advance" << std::endl;
            return false;
        }

        std::vector<std::pair<Integer, Integer>> boundary_elems;

        Simplex<Dim, ManifoldDim - 1> side;
        for (Integer i = 0; i < vol.n_elements(); ++i) {
            if (vol.is_active(i)) {
                const auto &e = vol.elem(i);
                const auto &adj = vol.dual_graph().adj(i);

                for (Integer side_num = 0; side_num < n_sides(e); ++side_num) {
                    auto a = adj[side_num];

                    if (a < 0) {
                        boundary_elems.emplace_back(i, side_num);
                    }
                }
            }
        }

        Integer surf_n_elems = boundary_elems.size();
        std::vector<Integer> node_index(vol.n_nodes(), 0);

        for (const auto &b : boundary_elems) {
            const auto &e = vol.elem(b.first);
            const auto side_num = b.second;

            e.side(side_num, side);

            // std::cout << e.side_tags[side_num] << std::endl;

            for (auto n : side.nodes) {
                node_index[n] = 1;
            }
        }

        Integer surf_n_nodes = std::accumulate(std::begin(node_index), std::end(node_index), 0);

        Integer cum_sum = 0;
        for (auto &n : node_index) {
            if (n == 1) {
                cum_sum += n;
                n = cum_sum;
            }
        }

        std::for_each(std::begin(node_index), std::end(node_index), [](Integer &v) { v -= 1; });

        surf.reserve(surf_n_elems, surf_n_nodes);

        const Integer vol_n_nodes = vol.n_nodes();
        for (Integer i = 0; i < vol_n_nodes; ++i) {
            if (node_index[i] >= 0) {
                surf.add_point(vol.point(i));
            }
        }

        for (const auto &b : boundary_elems) {
            const auto &vol_e = vol.elem(b.first);
            vol_e.side(b.second, side);

            for (auto &n : side.nodes) {
                assert(node_index[n] >= 0);
                n = node_index[n];
            }

            surf.add_elem(side);
            surf.tags().push_back(vol_e.side_tags[b.second]);
        }

        return true;
    }
}  // namespace mars

#endif  // MARS_EXTRACT_SURFACE_HPP
