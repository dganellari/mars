#ifndef MARS_EDGE_SPLIT_HPP
#define MARS_EDGE_SPLIT_HPP

#include "base.hpp"
#include "edge.hpp"

namespace mars {
	class EdgeSplit {
	public:
		EdgeSplit(const Edge &edge, const Integer midpoint = INVALID_INDEX)
		: edge(edge), midpoint(midpoint)
		{}

		Edge edge;
		Integer midpoint;
	};
}

#endif //MARS_EDGE_SPLIT_HPP
