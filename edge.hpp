#ifndef MARS_EDGE_HPP
#define MARS_EDGE_HPP

#include "base.hpp"
#include <array>
#include <vector>
#include <algorithm>
#include <initializer_list>

namespace mars {

	template<Integer N>
	class Side {
	public:
		static_assert(N > 0, "N cannot be zero");

		std::array<Integer, N> nodes;
		
		virtual ~Side() {}

		Side()
		{
			std::fill(nodes.begin(), nodes.end(), INVALID_INDEX);
		}

		Side(const std::array<Integer, N> &in)
		: nodes(in)
		{
			fix_ordering();
		}

		inline bool is_valid() const
		{
			for(auto n : nodes) {
				if(n == INVALID_INDEX) return false;
			}

			return true;
		}

		void fix_ordering()
		{
			std::sort(std::begin(nodes), std::end(nodes));
		}

		Side(const std::vector<Integer> &in)
		{
			assert(N == in.size());

			std::copy(std::begin(in), std::end(in), std::begin(nodes));
			std::sort(std::begin(nodes), std::end(nodes));
		}

		Side(std::initializer_list<Integer> in)
		{
			assert(N == in.size());

			std::copy(std::begin(in), std::end(in), std::begin(nodes));
			std::sort(std::begin(nodes), std::end(nodes));
		}

		inline bool operator==(const Side &other) const
		{
			for(Integer i = 0; i < N; ++i) {
				if(nodes[i] != other.nodes[i]) return false;
			}

			return true;
		}

		inline bool operator!=(const Side &other) const
		{
			return !((*this) == other);
		}

		inline bool operator<(const Side &other) const
		{
			for(Integer i = 0; i < N-1; ++i) {
				if(nodes[i] < other.nodes[i]) {
					return true;
				}

				if(nodes[i] > other.nodes[i]) {
					return false;
				}
			}

			return nodes[N-1] < other.nodes[N-1];
		}

		void describe(std::ostream &os) const
		{
			os << "(";
			for(Integer i = 0; i < N-1; ++i) {
				os << nodes[i] << ",";
			}

			os << nodes[N-1] << ")";
		}
	};

	class Edge : public Side<2> {
	public:
		Edge() : Side<2>() {}
		Edge(const Integer a_node, const Integer another_node)
		: Side({a_node, another_node})
		{}
	};
}

#endif //MARS_EDGE_HPP
