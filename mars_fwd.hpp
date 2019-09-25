#pragma once

namespace mars{

	class DefaultImplementation {};

	template<Integer Dim_, Integer Manifold_ = Dim_, class Implementation_ = DefaultImplementation>
	class Mesh;


	template<Integer Dim_, Integer Manifold_ = Dim_, class Implementation_ = DefaultImplementation>
	class Simplex;

	template<Integer N, Integer K, class Implementation_= DefaultImplementation>
	class CombinationsAux;

	template<Integer N, Integer ChooseM, class Implementation_= DefaultImplementation>
	class Combinations;
}



