#pragma once

namespace mars {

    class DefaultImplementation {};

    template <Integer Type, class Implementation_ = DefaultImplementation>
    class NonSimplex;

    template <Integer Dim_, Integer Manifold_ = Dim_, class Implementation_ = DefaultImplementation>
    class Simplex;

    template <Integer Dim_,
              Integer Manifold_ = Dim_,
              class Implementation_ = DefaultImplementation,
              class Simplex_ = Simplex<Dim_, Manifold_, Implementation_>>
    class Mesh;

    template <Integer N, Integer K, class Implementation_ = DefaultImplementation>
    class CombinationsAux;

    template <Integer N, Integer ChooseM, class Implementation_ = DefaultImplementation>
    class Combinations;

    template <class Mesh_, class EdgeSelect_>
    class Bisection;

    template <Integer N, class Implementation_ = DefaultImplementation>
    class SubManifoldElementMap;

    template <Integer N, class Implementation_ = DefaultImplementation>
    class Side;

    template <class Mesh, class Implementation_ = DefaultImplementation>
    class EdgeSelect;

    template <class Mesh, class Implementation_ = DefaultImplementation>
    class LongestEdgeSelect;
}  // namespace mars
