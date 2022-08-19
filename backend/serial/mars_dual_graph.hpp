#ifndef MARS_DUAL_GRAPH_HPP
#define MARS_DUAL_GRAPH_HPP

#include <algorithm>
#include <array>
#include <vector>
#include "mars_fwd.hpp"

namespace mars {

    template <Integer ManifoldDim>
    class DualGraph {
    public:
        inline Integer n_adjacients(const Integer id) const {
            assert(id >= 0);
            assert(id < dual_graph_.size());

            Integer ret = 0;
            for (auto a : dual_graph_[id]) {
                ret += a != INVALID_INDEX;
            }

            return ret;
        }

        inline const std::array<Integer, ManifoldDim + 1> &adj(const Integer id) const {
            assert(id >= 0);
            assert(id < dual_graph_.size());
            return dual_graph_[id];
        }

        inline std::array<Integer, ManifoldDim + 1> &adj(const Integer id) {
            assert(id >= 0);
            assert(id < dual_graph_.size());
            return dual_graph_[id];
        }

        inline const std::array<Integer, ManifoldDim + 1> &safe_adj(const Integer id) const {
            if (dual_graph_.size() <= id) {
                dual_graph_.resize(id + 1);
            }

            return dual_graph_[id];
        }

        inline std::array<Integer, ManifoldDim + 1> &safe_adj(const Integer id) {
            if (dual_graph_.size() <= id) {
                dual_graph_.resize(id + 1);
            }

            return dual_graph_[id];
        }

        template <Integer Dim>
        void update(const SimplicialMesh<Dim, ManifoldDim> &mesh, const bool force = false) {
            const Integer n_nodes = mesh.n_nodes();
            const Integer n_elements = mesh.n_elements();
            Integer el_index_size = 0;

            if (dual_graph_.size() == n_elements && !force) {
                return;
            }

            std::vector<std::vector<Integer>> node_2_element(n_nodes);
            dual_graph_.resize(n_elements);

            for (Integer i = 0; i < n_elements; ++i) {
                if (!mesh.is_active(i)) continue;

                const auto &e = mesh.elem(i);
                const Integer nn = ManifoldDim + 1;

                std::fill(std::begin(dual_graph_[i]), std::end(dual_graph_[i]), INVALID_INDEX);

                for (Integer k = 0; k < nn; ++k) {
                    node_2_element[e.nodes[k]].push_back(i);
                }
            }

            std::vector<Integer> el_offset(n_elements, 0);

            for (Integer i = 0; i < n_nodes; ++i) {
                const auto &elements = node_2_element[i];

                for (std::size_t e_i = 0; e_i < elements.size(); ++e_i) {
                    const Integer e = elements[e_i];

                    for (std::size_t e_i_adj = 0; e_i_adj < elements.size(); ++e_i_adj) {
                        if (e_i == e_i_adj) continue;

                        const Integer e_adj = elements[e_i_adj];

                        bool must_add = true;
                        for (Integer k = 0; k < el_offset[e]; ++k) {
                            if (e_adj == dual_graph_[e][k]) {
                                must_add = false;
                            }
                        }

                        if (must_add && mesh.have_common_side(e, e_adj)) {
                            assert(el_offset[e] < ManifoldDim + 1);
                            dual_graph_[e][el_offset[e]] = e_adj;
                            ++el_offset[e];
                        }
                    }
                }
            }

            for (Integer i = 0; i < n_elements; ++i) {
                if (!mesh.is_active(i)) continue;

                const std::array<Integer, ManifoldDim + 1> neighs = dual_graph_[i];

                std::fill(std::begin(dual_graph_[i]), std::end(dual_graph_[i]), INVALID_INDEX);

                for (Integer j = 0; j < neighs.size(); ++j) {
                    if (neighs[j] == INVALID_INDEX) break;

                    const auto s = mesh.common_side_num(i, neighs[j]);
                    assert(s != INVALID_INDEX);

                    if (dual_graph_[i][s] != INVALID_INDEX) {
                        std::cerr << "bad side numbering or creation" << std::endl;
                        assert(dual_graph_[i][s] == INVALID_INDEX);
                    }

                    dual_graph_[i][s] = neighs[j];
                }
            }
        }

        inline bool empty() const { return dual_graph_.empty(); }

        inline Integer size() const { return dual_graph_.size(); }

        void describe(std::ostream &os) const {
            for (std::size_t i = 0; i < dual_graph_.size(); ++i) {
                describe_adj(i, os);
            }
        }

        void describe_adj(const Integer i, std::ostream &os) const {
            os << "[" << i << "]:";
            for (std::size_t j = 0; j < adj(i).size(); ++j) {
                os << " " << adj(i)[j];
            }
            os << "\n";
        }

        void clear() { dual_graph_.clear(); }

    private:
        std::vector<std::array<Integer, ManifoldDim + 1>> dual_graph_;
    };
}  // namespace mars

#endif  // MARS_DUAL_GRAPH_HPP
