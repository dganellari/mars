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

		inline bool is_valid() const 
		{
			return edge.is_valid();
		}

		inline bool only_on_partition(const Integer partition_id) const
		{
			if(partitions.size() == 1) {
				return *partitions.begin() == partition_id;
			}

			return false;
		}

		// void determine_owner()
		// {
		// 	assert(!partitions.empty());
		// 	owner = *partitions.begin();
		// }

		inline void describe(std::ostream &os) const
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
