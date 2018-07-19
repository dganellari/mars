#ifndef MARS_EDGE_SPLIT_HPP
#define MARS_EDGE_SPLIT_HPP

#include "mars_base.hpp"
#include "mars_edge.hpp"

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

	template<class OutputStream>
	void write(
	    const EdgeSplit &edge_split,
	    OutputStream &os)
	{
	    write(edge_split.edge, os);
	    write(edge_split.midpoint, os);
	    write(edge_split.owner, os);
	    write(edge_split.partitions, os);
	}

	template<class InputStream>
	void read(
	    EdgeSplit &edge_split,
	    InputStream &is)
	{
	    read(edge_split.edge, is);
	    read(edge_split.midpoint, is);
	    read(edge_split.owner, is);
	    read(edge_split.partitions, is);
	}
}

#endif //MARS_EDGE_SPLIT_HPP
