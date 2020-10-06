#ifndef GENERATION_NON_SIMPLEX_MARS_NON_SIMPLEX_HPP_
#define GENERATION_NON_SIMPLEX_MARS_NON_SIMPLEX_HPP_

#include "mars_fwd.hpp"
#include <vector>

namespace mars
{
	template<Integer Type, class Implementation_>
	class NonSimplex final : public IElem
	{
	public:
		std::array<Integer,Type> nodes;
		std::array<Integer, Type> side_tags;

		Integer id = INVALID_INDEX;
		Integer parent_id = INVALID_INDEX;
		Integer block = INVALID_INDEX;

	//	std::vector<Integer> children;

        static constexpr Integer ElemType = Type;

		inline void get_nodes(std::vector<Integer> &nodes_copy) const override
		{
			nodes_copy.resize(nodes.size());
			std::copy(std::begin(nodes), std::end(nodes), std::begin(nodes_copy));
		}

		inline Integer type() const override
		{
			return Type;
		}

		inline Integer n_nodes() const override
		{
			return nodes.size();
		}

		inline Integer node(const Integer idx) const override
		{
			assert(idx < nodes.size());
			return nodes[idx];
		}

		inline Integer get_block() const override
		{
			return block;
		}

		inline void set_block(const Integer block_id) override
		{
			block = block_id;
		}

	};


	template<Integer Type>
	inline constexpr static Integer n_nodes(const NonSimplex<Type> &e)
	{
		return e.nodes.size();
	}


	template<Integer Type>
	inline bool has_node(const NonSimplex<Type> &s, const Integer &node)
	{
		return std::find(s.nodes.begin(), s.nodes.end(), node) != s.nodes.end();
	}


	using Quad4Elem = NonSimplex<ElementType::Quad4>;
	using Hex8Elem  = NonSimplex<ElementType::Hex8>;
}

#endif /* GENERATION_NON_SIMPLEX_MARS_NON_SIMPLEX_HPP_ */
