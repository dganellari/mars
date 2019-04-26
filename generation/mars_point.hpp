#pragma once

namespace mars {
namespace unit_generation {

#include "mars_vector.hpp"
template<typename T, Integer Dim>
class Point: public Vector<T,Dim> {
public:

	bool isActive() const {
		return active;
	}

	void setActive() {
		this->active = true;
	}

	Point() {
	}

	Point(std::initializer_list<T> values, bool active = false) :
			Vector<T,Dim>(values) {
		this->active = active;
	}

private:
	bool active = false;

};

}
}
