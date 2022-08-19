#ifndef MARS_EDGE_NODE_MAP_HPP
#define MARS_EDGE_NODE_MAP_HPP

#include "mars_edge.hpp"

#include <cassert>
#include <map>
#include <memory>
#include <vector>

namespace mars {
    class EdgeNodeMap {
    public:
        inline void update(const Integer an_edge_node, const Integer another_edge_node, const Integer midpoint_node) {
            mapping_[Edge(an_edge_node, another_edge_node)] = midpoint_node;
        }

        inline void update(const Edge &edge, const Integer midpoint_node) {
            assert(edge.is_valid());
            mapping_[edge] = midpoint_node;
        }

        inline Integer get(const Integer an_edge_node, const Integer another_edge_node) const {
            return get(Edge(an_edge_node, another_edge_node));
        }

        inline Integer get(const Edge &edge) const {
            auto it = mapping_.find(edge);
            if (it == mapping_.end()) {
                return INVALID_INDEX;
            }

            return it->second;
        }

        void describe(std::ostream &os) const {
            for (const auto &m : mapping_) {
                os << "(" << m.first.nodes[0] << "," << m.first.nodes[1] << ") -> " << m.second << "\n";
            }
        }

        inline std::map<Edge, Integer>::const_iterator begin() const { return mapping_.begin(); }

        inline std::map<Edge, Integer>::const_iterator end() const { return mapping_.end(); }

        void clear() { mapping_.clear(); }

        inline Integer size() const { return mapping_.size(); }

        inline bool empty() const { return mapping_.empty(); }

        std::map<Edge, Integer> &mapping() { return mapping_; }

        std::map<Edge, Integer> mapping_;
    };
}  // namespace mars

#endif  // MARS_EDGE_NODE_MAP_HPP
