#ifndef MARS_RED_GREEN_REFINEMENT_HPP
#define MARS_RED_GREEN_REFINEMENT_HPP

#include <algorithm>
#include <array>
#include <deque>
#include <set>
#include <vector>
#include "mars_fwd.hpp"

#include "mars_simplex.hpp"

#include "mars_edge_element_map.hpp"
#include "mars_edge_node_map.hpp"

namespace mars {

    enum RefinementFlag {
        NONE = 0,
        RED = 1,
        GREEN_1 = 2,
        GREEN_2 = 3,
        GREEN_3 = 4,
        CHILD_OF_GREEN = 5,
        PARENT_PROMOTED_TO_RED = 6,
        BISECTION = 7
    };

    template <class Mesh>
    class RedGreenRefinement {
    public:
        static const Integer Dim = Mesh::Dim;
        static const Integer ManifoldDim = Mesh::ManifoldDim;

        class Refinement {
        public:
            class Trigger {
            public:
                // red|green|...
                Integer flag;

                // element id
                Integer element;

                // self=Manifold|side=Manifold-1|sub-side=Manifold-2|....|edge=1
                Integer type;

                // face index|edge index
                Integer index;
            };

            std::vector<Trigger> triggers;
        };

        RedGreenRefinement(Mesh &mesh) : mesh(mesh) {}

        void update_children(const Integer id) {
            for (auto c : mesh.elem(id).children) {
                update_element(c);
            }
        }

        void update_element(const Integer id, const bool promoted = false) {
            const auto &e = mesh.elem(id);
            const auto p_id = e.parent_id;
            // const auto &parent_e = mesh.elem(p_id);

            // FIXME?
            if (p_id == INVALID_INDEX) return;

            mesh.elem(id).block = mesh.elem(p_id).block;

            // emplace if not there
            mesh.dual_graph().safe_adj(id);

            std::fill(mesh.dual_graph().adj(id).begin(), mesh.dual_graph().adj(id).end(), INVALID_INDEX);

            std::map<Integer, std::vector<Integer>> local_node_2_element;

            for (auto c : mesh.elem(p_id).children) {
                if (mesh.elem(c).children.empty() || !promoted) {
                    for (auto n : mesh.elem(c).nodes) {
                        local_node_2_element[n].push_back(c);
                    }
                } else {
                    for (auto cc : mesh.elem(c).children) {
                        for (auto n : mesh.elem(cc).nodes) {
                            local_node_2_element[n].push_back(cc);
                        }
                    }
                }
            }

            for (auto a : mesh.dual_graph().adj(p_id)) {
                if (a == INVALID_INDEX) continue;

                for (auto c : mesh.elem(a).children) {
                    for (auto n : mesh.elem(c).nodes) {
                        local_node_2_element[n].push_back(c);
                    }
                }
            }

            bool updated = false;

            for (Integer i = 0; i < ManifoldDim + 1; ++i) {
                auto it = local_node_2_element.find(e.nodes[i]);
                if (it == local_node_2_element.end()) continue;

                for (auto other : it->second) {
                    if (id == other) continue;

                    if (mesh.have_common_side(id, other)) {
                        updated = true;
                        auto &e_adj = mesh.dual_graph().safe_adj(id);
                        e_adj[mesh.common_side_num(id, other)] = other;

                        auto &other_adj = mesh.dual_graph().safe_adj(other);
                        other_adj[mesh.common_side_num(other, id)] = id;
                    }
                }
            }

            if (!updated) {
                std::cerr << "element " << id << " with parent " << p_id << " not updated " << std::endl;
                assert(updated);
            }

            if (id == 9) {
                mesh.dual_graph().describe_adj(9, std::cout);
            }
        }

        void promote_to_red(const Integer green_parent) {
            std::cout << "promoted to red " << green_parent << std::endl;
            assert(is_green(green_parent));
            set_refinement_flag(green_parent, RED);
            mesh.deactivate_children(green_parent);
            mesh.set_active(green_parent, true);
        }

        void breadth_first_promote_green(const Integer id, std::set<Integer> &promoted) {
            std::deque<Integer> to_visit;
            to_visit.push_back(id);
            std::set<Integer> visited;

            promoted.clear();

            while (!to_visit.empty()) {
                auto f = to_visit.front();
                to_visit.pop_front();
                auto parent = mesh.elem(f).parent_id;
                visited.insert(f);

                if (is_green(parent)) {
                    promoted.insert(parent);
                    mesh.set_active(f, false);
                }

                const auto &adj = mesh.dual_graph().adj(f);

                for (auto a : adj) {
                    if (a != INVALID_INDEX && mesh.is_active(a) && is_green(mesh.elem(a).parent_id)) {
                        if (visited.find(a) == visited.end()) {
                            to_visit.push_back(a);
                        }
                    }
                }
            }

            for (auto p : promoted) {
                promote_to_red(p);
            }
        }

        bool red_refine(const std::vector<Integer> &elements_to_refine) {
            update_dual_graph();
            assert(mesh.n_elements() == mesh.dual_graph().size());

            if (elements_to_refine.empty()) return false;

            if (refinement_flag_.empty()) {
                refinement_flag_.resize(mesh.n_elements(), NONE);
            }

            std::set<Integer> promoted_neighs;
            std::vector<Integer> promoted;
            std::vector<Integer> red_elements;

            bool has_refinement = false;
            for (auto e : elements_to_refine) {
                if (!mesh.is_valid(e)) {
                    std::cerr << e << " skipping invalid element" << e << std::endl;
                    continue;
                }

                if (!mesh.is_active(e)) {
                    std::cerr << e << " cannot refine inactive" << std::endl;
                    continue;
                }

                if (!is_green(mesh.elem(e).parent_id)) {
                    set_refinement_flag(e, RED);
                    red_elements.push_back(e);
                }

                breadth_first_promote_green(e, promoted_neighs);
                promoted.insert(promoted.end(), promoted_neighs.begin(), promoted_neighs.end());
                has_refinement = true;
            }

            for (auto e : promoted) {
                red_elements.push_back(e);
                breadth_first_promote_green(e, promoted_neighs);
                promoted.insert(promoted.end(), promoted_neighs.begin(), promoted_neighs.end());
            }

            if (!has_refinement) {
                std::cerr << "did not refine" << std::endl;
                return false;
            }

            for (auto e : red_elements) {
                std::cout << e << " ";
            }

            std::cout << std::endl;

            for (auto &e : red_elements) {
                red_refine_element(e);
            }

            for (auto &e : red_elements) {
                update_element(e);
            }

            // for(auto &e : promoted) {
            // 	update_element(e, true);
            // }

            for (auto e : red_elements) {
                update_children(e);
            }

            // update_dual_graph();
            assert(mesh.n_elements() == mesh.dual_graph().size());

            mesh.tags() = refinement_flag_;
            return true;
        }

        void green_refine() {
            std::vector<Integer> green_elements, promoted_elements;
            for (Integer i = 0; i < mesh.n_elements(); ++i) {
                if (mesh.is_active(i) && !is_red(i)) {
                    const auto &adj = mesh.dual_graph().adj(i);

                    Integer n_red_neighs = 0;
                    for (auto a : adj) {
                        if (a == INVALID_INDEX || mesh.is_active(a)) continue;
                        if (is_red(a)) {
                            n_red_neighs++;
                        }
                    }

                    if (n_red_neighs == 0) continue;

                    switch (n_red_neighs) {
                        case 1: {
                            set_refinement_flag(i, GREEN_1);
                            green_elements.push_back(i);
                            break;
                        }

                        case 2: {
                            set_refinement_flag(i, GREEN_2);
                            green_elements.push_back(i);
                            break;
                        }

                        default: {
                            set_refinement_flag(i, RED);
                            promoted_elements.push_back(i);
                            break;
                        }
                    }
                }
            }

            if (!promoted_elements.empty()) {
                for (auto e : green_elements) {
                    set_refinement_flag(e, NONE);
                }

                std::cout << "handling promoted elements" << std::endl;
                red_refine(promoted_elements);
                green_refine();
                return;
            }

            for (auto e : green_elements) {
                green_refine_element(e);
                update_element(e);
            }

            update_dual_graph();

            mesh.tags() = refinement_flag_;
        }

        bool refine(const std::vector<Integer> &elements_to_refine) {
            if (!red_refine(elements_to_refine)) return false;
            green_refine();
            return true;
        }

        inline void green_refine_element(const Integer element_id) {
            mesh.set_active(element_id, false);

            const auto &adj = mesh.dual_graph().adj(element_id);
            std::array<Integer, ManifoldDim> side_flags;
            std::vector<Integer> red_side_index;
            std::vector<Integer> green_side_index;

            Integer n_red_neighs = 0;
            Integer n_green_neighs = 0;

            Integer k = 0;
            for (auto a : adj) {
                if (a == INVALID_INDEX) {
                    side_flags[k++] = INVALID_INDEX;
                    continue;
                }

                if (refinement_flag(a) == RED) {
                    red_side_index.push_back(k);
                    ++n_red_neighs;
                } else if (refinement_flag(a) == GREEN_1) {
                    ++n_green_neighs;
                    green_side_index.push_back(k);
                }

                side_flags[k++] = refinement_flag(a);
            }

            Simplex<Dim, ManifoldDim - 1> side_1, side_2;
            Simplex<Dim, ManifoldDim> child;
            child.parent_id = element_id;
            mesh.elem(element_id).children.clear();

            switch (n_red_neighs) {
                case 1: {
                    mesh.elem(element_id).side(red_side_index[0], side_1);

                    Integer n0 = side_1.nodes[0];
                    Integer n1 = side_1.nodes[1];
                    const Integer midpoint = edge_node_map_.get(n0, n1);
                    const Integer opposite = mesh.elem(element_id).vertex_opposite_to_side(red_side_index[0]);

                    child.nodes[0] = n0;
                    child.nodes[1] = opposite;
                    child.nodes[2] = midpoint;

                    mesh.elem(element_id).children.push_back(add_elem(child));

                    child.nodes[0] = n1;
                    child.nodes[1] = midpoint;
                    child.nodes[2] = opposite;

                    mesh.elem(element_id).children.push_back(add_elem(child));
                    break;
                }

                case 2: {
                    mesh.elem(element_id).side(red_side_index[0], side_1);

                    Integer n1_0 = side_1.nodes[0];
                    Integer n1_1 = side_1.nodes[1];
                    const Integer midpoint_1 = edge_node_map_.get(n1_0, n1_1);
                    const Integer opposite_1 = mesh.elem(element_id).vertex_opposite_to_side(red_side_index[0]);

                    mesh.elem(element_id).side(red_side_index[1], side_2);

                    Integer n2_0 = side_2.nodes[0];
                    Integer n2_1 = side_2.nodes[1];
                    const Integer midpoint_2 = edge_node_map_.get(n2_0, n2_1);
                    const Integer opposite_2 = mesh.elem(element_id).vertex_opposite_to_side(red_side_index[1]);

                    child.nodes[0] = midpoint_1;
                    child.nodes[1] = opposite_2;
                    child.nodes[2] = opposite_1;

                    mesh.elem(element_id).children.push_back(add_elem(child));

                    child.nodes[0] = midpoint_1;
                    child.nodes[1] = midpoint_2;
                    child.nodes[2] = n2_0;

                    mesh.elem(element_id).children.push_back(add_elem(child));

                    child.nodes[0] = midpoint_1;
                    child.nodes[1] = n2_1;
                    child.nodes[2] = midpoint_2;

                    mesh.elem(element_id).children.push_back(add_elem(child));
                    break;
                }

                case 3: {
                    std::cerr << "[" << element_id << "] should have been promoted to RED" << std::endl;
                    break;
                }
                default: {
                    assert(false);
                }
            }
        }

        void update_dual_graph(const bool force = false) { mesh.update_dual_graph(force); }

        inline void red_refine_element(const Integer element_id) {
            static const Integer NSubs = NSubSimplices<ManifoldDim>::value;
            static_assert(NSubSimplices<ManifoldDim>::value > 0, "!");

            std::vector<Vector<Real, Dim>> parent_points;
            mesh.points(element_id, parent_points);

            std::array<Simplex<Dim, ManifoldDim>, NSubs> children;
            std::vector<Vector<Real, Dim>> children_points;
            auto interp = std::make_shared<SimplexInterpolator<ManifoldDim>>();

            Simplex<Dim, ManifoldDim> modified_e = mesh.elem(element_id);

            if (ManifoldDim == 4) {
                // 4D hack
                std::sort(modified_e.nodes.begin(), modified_e.nodes.end());
            }

            red_refinement<Dim, ManifoldDim, NSubs>(modified_e, parent_points, children, children_points, *interp);

            if (interp_.size() <= element_id) {
                interp_.resize(element_id + 1);
            }

            std::vector<Integer> point_ids(interp->rows(), INVALID_INDEX);

            for (Integer i = 0; i < ManifoldDim + 1; ++i) {
                point_ids[i] = modified_e.nodes[i];

                for (Integer j = i + 1; j < ManifoldDim + 1; ++j) {
                    Integer offset = midpoint_index<ManifoldDim>(i, j);
                    point_ids[offset] = edge_node_map_.get(modified_e.nodes[i], modified_e.nodes[j]);

                    if (point_ids[offset] == INVALID_INDEX) {
                        const auto new_id = mesh.add_point(children_points[offset]);
                        edge_node_map_.update(modified_e.nodes[i], modified_e.nodes[j], new_id);

                        point_ids[offset] = new_id;
                        assert(new_id < mesh.n_nodes());
                    }
                }
            }

            interp_[element_id] = interp;
            mesh.elem(element_id).children.clear();

            for (auto &c : children) {
                for (Integer i = 0; i < ManifoldDim + 1; ++i) {
                    c.nodes[i] = point_ids[c.nodes[i]];
                }

                // std::sort(c.nodes.begin(), c.nodes.end()); add_elem(c);

                // 4D hack
                const auto c_id = add_elem(c);
                mesh.repair_element(c_id, ManifoldDim != 4);
                mesh.elem(element_id).children.push_back(c_id);
            }

            assert(mesh.elem(element_id).children.size() == NSubs);
            mesh.set_active(element_id, false);
        }

        void uniformly_refine(const Integer n_levels = 1) {
            for (Integer l = 0; l < n_levels; ++l) {
                auto ne = mesh.n_elements();

                for (Integer i = 0; i < ne; ++i) {
                    if (mesh.is_active(i)) {
                        red_refine_element(i);
                    }
                }
            }
        }

        inline bool is_green(const Integer id) const {
            if (INVALID_INDEX == id) return false;

            switch (refinement_flag(id)) {
                case GREEN_1:
                    return true;
                case GREEN_2:
                    return true;
                case GREEN_3:
                    return true;
                default:
                    return false;
            }
        }

        inline bool is_red(const Integer id) const { return refinement_flag(id) == RED; }

        inline void set_refinement_flag(const Integer &element_id, const Integer flag) {
            refinement_flag_[element_id] = flag;
        }

        inline Integer refinement_flag(const Integer &element_id) const { return refinement_flag_[element_id]; }

        inline Mesh &get_mesh() { return mesh; }

        inline const Mesh &get_mesh() const { return mesh; }

        const EdgeNodeMap &edge_node_map() const { return edge_node_map_; }

    private:
        inline Integer add_elem(const Simplex<Dim, ManifoldDim> &elem) {
            refinement_flag_.push_back(NONE);
            return mesh.add_elem(elem);
        }

        Mesh &mesh;

        // refinement
        std::vector<Integer> refinement_flag_;
        std::vector<std::shared_ptr<SimplexInterpolator<ManifoldDim>>> interp_;
        EdgeNodeMap edge_node_map_;

        MultilevelElementMap<ManifoldDim, 2> mlem;
    };
}  // namespace mars

#endif  // MARS_RED_GREEN_REFINEMENT_HPP
