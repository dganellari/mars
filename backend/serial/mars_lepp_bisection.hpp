#ifndef MARS_LEPP_BISECTION_HPP_
#define MARS_LEPP_BISECTION_HPP_

#include "mars_bisection.hpp"

#include <stack>
#include "mars_err.hpp"

#define select_err 1

namespace mars {

    template <typename T>
    class TreeNode {
    public:
        std::vector<TreeNode>& getChildren() { return children; }

        const std::vector<TreeNode>& getChildren() const { return children; }

        const T get_id() const { return id; }

        void set_id(T id) { this->id = id; }

        void add_child(const T element_id) { children.emplace_back(element_id, std::vector<TreeNode>()); }

        TreeNode() = default;

        TreeNode(T elem_id, std::vector<TreeNode> elem_children)
            : id(std::move(elem_id)), children(std::move(elem_children)) {}

    private:
        T id;
        std::vector<TreeNode> children;
    };

    template <class Mesh_>
    class LeppBisection : public Bisection<Mesh_> {
        // int global =0;

    public:
        using Mesh = Mesh_;

        LeppBisection(Mesh& mesh) : Bisection<Mesh>(mesh) {}

        void refine_element(const Integer element_id) {
            if (!Bisection<Mesh>::edge_select()->can_refine(Bisection<Mesh>::get_mesh(), element_id)) {
                Bisection<Mesh>::get_incomplete_elements().push_back(element_id);
                assert(!Bisection<Mesh>::get_fail_if_not_refine());
                return;
            }

            TreeNode<Integer> lepp_tree;
            lepp_tree.set_id(element_id);

            while (Bisection<Mesh>::get_mesh().is_active(element_id)) {
                depth_first(lepp_tree);
            }
        }

    private:
        void set_edge_select(const std::shared_ptr<EdgeSelect<Mesh>>& edge_select) {
            if (edge_select->name() != "LongestEdge")
                errorx(select_err,
                       "%s %d %s Error: Calling set_edge_select on LeppBisection with %s is not supported!",
                       __FILE__,
                       __LINE__,
                       __func__,
                       edge_select->name().c_str());
        }

        void compute_lepp(TreeNode<Integer>& node) {
            Integer edge_num =
                Bisection<Mesh>::edge_select()->stable_select(Bisection<Mesh>::get_mesh(), node.get_id());
            // Integer edge_num = Bisection<Mesh>::edge_select()->select(Bisection<Mesh>::get_mesh(), element_id);
            // Integer edge_num = Bisection<Mesh>::edge_select()->select(Bisection<Mesh>::get_mesh(), element_id,
            // Bisection<Mesh>::edge_element_map());

            Edge edge;
            Bisection<Mesh>::get_mesh().elem(node.get_id()).edge(edge_num, edge.nodes[0], edge.nodes[1]);
            edge.fix_ordering();

            auto incidents = Bisection<Mesh>::edge_element_map().elements(edge);

            if (is_terminal(edge, node, incidents)) {
                for (const Integer element : incidents) {
                    if (Bisection<Mesh>::get_mesh().is_active(element) &&
                        Bisection<Mesh>::edge_select()->can_refine(Bisection<Mesh>::get_mesh(), element))
                        Bisection<Mesh>::bisect_element(element, edge);
                }
            }
        }

        bool is_terminal(const Edge& edge, TreeNode<Integer>& node, const std::vector<Integer>& incidents) {
            bool terminal = true;  // in case the elements share the longest edge or there is only one incident (itself)
                                   // - meaning the longest edge is on the boundary.

            for (auto i : incidents) {
                if (Bisection<Mesh>::get_mesh().is_active(i) &&
                    Bisection<Mesh>::edge_select()->can_refine(Bisection<Mesh>::get_mesh(), i)) {
                    Edge new_edge;
                    const Integer edge_num =
                        Bisection<Mesh>::edge_select()->stable_select(Bisection<Mesh>::get_mesh(), i);
                    // const Integer edge_num = Bisection<Mesh>::edge_select()->select(Bisection<Mesh>::get_mesh(),
                    // edge, i);
                    Bisection<Mesh>::get_mesh().elem(i).edge(edge_num, new_edge.nodes[0], new_edge.nodes[1]);
                    new_edge.fix_ordering();

                    if (edge != new_edge) {
                        terminal = false;
                        node.add_child(i);
                    }
                }
            }

            return terminal;
        }

        bool is_leaf(const TreeNode<Integer>& node) {
            for (const auto& i : node.getChildren()) {
                if (Bisection<Mesh>::get_mesh().is_active(i.get_id())) {
                    return false;
                }
            }

            return true;
        }

        void depth_first(TreeNode<Integer>& node) {
            if (is_leaf(node)) {
                compute_lepp(node);
            }

            for (auto& child : node.getChildren()) {
                if (Bisection<Mesh>::get_mesh().is_active(
                        child.get_id()))  // cases when the child is added as active but on the way it is bisected and
                                          // becomes inactive.
                    depth_first(child);
            }
        }
    };

    // optimized using a stack instead of the TreeNode.
    template <>
    class LeppBisection<Mesh2> : public Bisection<Mesh2, LongestEdgeSelect<Mesh2>> {
    public:
        using Mesh = Mesh2;
        template <class T>
        using Bisection = Bisection<T, LongestEdgeSelect<T>>;

        LeppBisection(Mesh& mesh) : Bisection<Mesh>(mesh) {}

        Integer longest_edge_neighbor(Integer element_id, Edge edge) {
            auto incidents = Bisection<Mesh>::edge_element_map().elements(edge);

            Integer index = INVALID_INDEX;  // in case a boundary longest edge

            for (auto i : incidents) {
                assert(has_edge(Bisection<Mesh>::get_mesh().elem(i), edge.nodes[0], edge.nodes[1]));

                if (!Bisection<Mesh>::get_mesh().is_active(i) || i == element_id) continue;

                if (Bisection<Mesh>::edge_select()->can_refine(Bisection<Mesh>::get_mesh(), i)) index = i;
            }

            return index;
        }

        void refine_element(const Integer element_id) {
            if (!Bisection<Mesh>::edge_select()->can_refine(Bisection<Mesh>::get_mesh(), element_id)) {
                Bisection<Mesh>::get_incomplete_elements().push_back(element_id);
                assert(!Bisection<Mesh>::get_fail_if_not_refine());
                return;
            }

            std::stack<Integer> s;
            s.push(element_id);

            while (!s.empty()) {
                Integer top = s.top();

                Integer edge_num = Bisection<Mesh>::edge_select()->stable_select(Bisection<Mesh>::get_mesh(), top);

                Edge edge;
                Bisection<Mesh>::get_mesh().elem(top).edge(edge_num, edge.nodes[0], edge.nodes[1]);
                edge.fix_ordering();

                Integer le_nhb = longest_edge_neighbor(top, edge);

                Edge new_edge;

                // If the longest edge of top is not boundary edge (if there exists a longest edge neighbor).
                if (le_nhb != INVALID_INDEX) {
                    // get the longest edge of that neighbor.
                    edge_num = Bisection<Mesh>::edge_select()->select(Bisection<Mesh>::get_mesh(), edge, le_nhb);
                    //				edge_num =
                    // Bisection<Mesh>::edge_select()->select(Bisection<Mesh>::get_mesh(), edge, le_nhb);

                    Bisection<Mesh>::get_mesh().elem(le_nhb).edge(edge_num, new_edge.nodes[0], new_edge.nodes[1]);
                    new_edge.fix_ordering();
                }

                // if top is terminal triangle (if the longest edge is shared or the longest edge is boundary).
                if (le_nhb == INVALID_INDEX || edge == new_edge) {
                    // bisect top;
                    Bisection<Mesh>::bisect_element(top, edge);

                    if (le_nhb != INVALID_INDEX) {
                        // bisect le_nhb
                        Bisection<Mesh>::bisect_element(le_nhb, new_edge);
                    }
                    s.pop();
                } else
                    s.push(le_nhb);
            }
        }

    private:
        void set_edge_select(const std::shared_ptr<LongestEdgeSelect<Mesh>>& edge_select) final {
            if (edge_select->name() != "LongestEdge")
                errorx(select_err,
                       "%s %d %s Error: Calling set_edge_select on LeppBisection with %s is not supported!",
                       __FILE__,
                       __LINE__,
                       __func__,
                       edge_select->name().c_str());
        }
    };

    // LEPP re-computation becomes the bottleneck of this algorithm.
    template <class Mesh_>
    class RecursiveLeppBisection : public Bisection<Mesh_> {
    public:
        using Mesh = Mesh_;

        RecursiveLeppBisection(Mesh& mesh) : Bisection<Mesh>(mesh) {}

        void refine_element(const Integer element_id) {
            if (!Bisection<Mesh>::edge_select()->can_refine(Bisection<Mesh>::get_mesh(), element_id)) {
                Bisection<Mesh>::get_incomplete_elements().push_back(element_id);
                assert(!Bisection<Mesh>::get_fail_if_not_refine());
                return;
            }

            while (Bisection<Mesh>::get_mesh().is_active(element_id)) {
                std::map<Edge, std::vector<Integer>> terminal_edges;
                compute_lepp(terminal_edges, element_id);

                for (auto const& star : terminal_edges) {
                    for (const Integer element : star.second) {
                        if (Bisection<Mesh>::get_mesh().is_active(element) &&
                            Bisection<Mesh>::edge_select()->can_refine(Bisection<Mesh>::get_mesh(), element))
                            Bisection<Mesh>::bisect_element(element, star.first);
                    }
                }
            }
        }

    private:
        void set_edge_select(const std::shared_ptr<EdgeSelect<Mesh_>>& edge_select) final {
            if (edge_select->name() != "LongestEdge")
                errorx(select_err,
                       "%s %d %s Error: Calling set_edge_select on LeppBisection with %s is not supported!",
                       __FILE__,
                       __LINE__,
                       __func__,
                       edge_select->name().c_str());
        }

        void compute_lepp(std::map<Edge, std::vector<Integer>>& lepp, const Integer element_id) {
            Integer edge_num = Bisection<Mesh>::edge_select()->stable_select(Bisection<Mesh>::get_mesh(), element_id);
            // Integer edge_num = Bisection<Mesh>::edge_select()->select(Bisection<Mesh>::get_mesh(), element_id);
            // Integer edge_num = Bisection<Mesh>::edge_select()->select(Bisection<Mesh>::get_mesh(), element_id,
            // Bisection<Mesh>::edge_element_map());

            Edge edge;
            Bisection<Mesh>::get_mesh().elem(element_id).edge(edge_num, edge.nodes[0], edge.nodes[1]);
            edge.fix_ordering();

            auto incidents = Bisection<Mesh>::edge_element_map().elements(edge);

            if (is_terminal(edge, lepp, incidents)) {
                lepp[edge] = incidents;
            }
        }

        bool is_terminal(const Edge& edge,
                         std::map<Edge, std::vector<Integer>>& lepp,
                         const std::vector<Integer>& incidents) {
            bool terminal = true;  // in case the elements share the longest edge or there is only one incident (itself)
                                   // - meaning the longest edge is on the boundary.

            for (auto i : incidents) {
                if (Bisection<Mesh>::get_mesh().is_active(i) &&
                    Bisection<Mesh>::edge_select()->can_refine(Bisection<Mesh>::get_mesh(), i)) {
                    Edge new_edge;
                    //				const Integer edge_num =
                    // Bisection<Mesh>::edge_select()->stable_select(Bisection<Mesh>::get_mesh(),  i);
                    const Integer edge_num =
                        Bisection<Mesh>::edge_select()->select(Bisection<Mesh>::get_mesh(), edge, i);
                    Bisection<Mesh>::get_mesh().elem(i).edge(edge_num, new_edge.nodes[0], new_edge.nodes[1]);
                    new_edge.fix_ordering();

                    if (edge != new_edge) {
                        terminal = false;
                        compute_lepp(lepp, i);
                    }
                }
            }

            return terminal;
        }
    };

}  // namespace mars

#endif /* MARS_LEPP_BISECTION_HPP_ */
