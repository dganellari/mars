#ifndef MARS_EDGE_HPP
#define MARS_EDGE_HPP

#include "base.hpp"
#include <array>

namespace mars {
	class Edge {
	public:
		std::array<Integer, 2> nodes;

		Edge(const Integer a_node, const Integer another_node)
		{
			if(a_node < another_node) {
				nodes[0] = a_node;
				nodes[1] = another_node;
			} else {
				nodes[0] = another_node;
				nodes[1] = a_node;
			}
		}

		inline bool operator<(const Edge &other) const
		{
			if(nodes[0] < other.nodes[0]) {
				return true;
			}

			if(nodes[0] > other.nodes[0]) {
				return false;
			}

			return nodes[1] < other.nodes[1];
		}
	};
}

#endif //MARS_EDGE_HPP
