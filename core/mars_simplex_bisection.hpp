#pragma once

#if defined(__CUDACC__) || defined(__HIPCC__)
#define MARS_BISECTION_HD __host__ __device__
#else
#define MARS_BISECTION_HD
#endif

namespace mars::bisection
{

template<class KeyType>
struct Edge
{
    KeyType first{}, second{};

    MARS_BISECTION_HD Edge() = default;
    MARS_BISECTION_HD Edge(KeyType a, KeyType b)
        : first(a < b ? a : b)
        , second(a < b ? b : a)
    {
    }

    MARS_BISECTION_HD bool operator==(const Edge& other) const
    {
        return first == other.first && second == other.second;
    }

    MARS_BISECTION_HD bool operator!=(const Edge& other) const { return !(*this == other); }

    MARS_BISECTION_HD bool operator<(const Edge& other) const
    {
        return first < other.first || (first == other.first && second < other.second);
    }
};

// Identities, rather than local connectivity indices, resolve equal lengths.
template<class KeyType, class RealType>
struct LongestEdge
{
    Edge<KeyType> edge;
    RealType squared_length = RealType(-1);
    int index               = -1;

    MARS_BISECTION_HD void consider(RealType length, KeyType a, KeyType b, int candidate_index)
    {
        Edge<KeyType> candidate(a, b);
        if (length > squared_length || (length == squared_length && candidate < edge))
        {
            squared_length = length;
            edge           = candidate;
            index          = candidate_index;
        }
    }
};

// Visit every incident: the legacy visitor also records the nonterminal dependencies.
template<class EdgeType, class Count, class Selection>
MARS_BISECTION_HD bool terminal_edge(const EdgeType& edge, Count count, Selection selection)
{
    bool terminal = true;
    for (Count i = 0; i < count; ++i)
        if (selection(i) != edge) terminal = false;
    return terminal;
}

// Endpoint replacement preserves the parent's signed orientation.
template<int Nodes, class Parent, class KeyType, class Child>
MARS_BISECTION_HD void child_connectivity(const Parent& parent, int replaced, KeyType midpoint, Child&& child)
{
    for (int i = 0; i < Nodes; ++i)
        child[i] = i == replaced ? midpoint : parent[i];
}

template<class RealType>
MARS_BISECTION_HD RealType midpoint_value(RealType a, RealType b)
{
    return RealType(0.5) * a + RealType(0.5) * b;
}

// A face containing the split edge belongs to both children. The other two
// faces each survive in the child that keeps their opposite endpoint.
MARS_BISECTION_HD inline unsigned face_children(int opposite, int edge_a, int edge_b)
{
    return opposite == edge_a ? 2u : opposite == edge_b ? 1u : 3u;
}

} // namespace mars::bisection

#undef MARS_BISECTION_HD
