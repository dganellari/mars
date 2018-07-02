#ifndef MARS_EDGE_SPLIT_HPP
#define MARS_EDGE_SPLIT_HPP

#include "base.hpp"
#include "edge.hpp"

#include <ostream>

namespace mars {
	class EdgeSplit {
	public:
		EdgeSplit(const Edge &edge = Edge(), const Integer midpoint = INVALID_INDEX, const Integer owner = INVALID_INDEX)
		: edge(edge), midpoint(midpoint), owner(owner)
		{}

		Edge edge;
		Integer midpoint;
		Integer owner;
		std::set<Integer> partitions;

		void describe(std::ostream &os) const
		{
			os << "split";
			edge.describe(os); 
			os << "\nmidpoint:   " << midpoint << "\n";
			os << "owner:      " << owner << std::endl;
			os << "partitions:";
			for(auto p : partitions) {
				os << " " << p;
			}

			os << "\n";

		}
	};
}

#endif //MARS_EDGE_SPLIT_HPP
