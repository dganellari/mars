#ifndef MARS_EDGE_NODE_MAP_HPP
#define MARS_EDGE_NODE_MAP_HPP

#include "edge.hpp"

#include <map>
#include <vector>
#include <memory>
#include <cassert>

namespace mars {
	class EdgeNodeMap {
	public:
		inline void update(const Integer an_edge_node,
			        const Integer another_edge_node,
			        const Integer midpoint_node)
		{
			mapping_[Edge(an_edge_node, another_edge_node)] = midpoint_node;
		}

		inline Integer get(const Integer an_edge_node,
			        const Integer another_edge_node) const
		{
			auto it = mapping_.find(Edge(an_edge_node, another_edge_node));
			if(it == mapping_.end()) {
				return INVALID_INDEX;
			}

			return it->second;
		}

		std::map<Edge, Integer> mapping_;
	};
}

#endif //MARS_EDGE_NODE_MAP_HPP
