#pragma once

namespace mars{

	class DefaultImplementation {};

	template<Integer Dim_, Integer Manifold_ = Dim_, class Implementation_ = DefaultImplementation>
	class Mesh;
}



