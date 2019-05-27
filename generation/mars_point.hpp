#pragma once

#include "mars_vector.hpp"

namespace mars {
namespace generation {

template<typename T, Integer Dim>
class Point: public Vector<T,Dim> {
public:

	bool is_active() const {
		return active_;
	}

	void set_active() {
		this->active_ = true;
	}

	Point() {
	}

	Point(std::initializer_list<T> values, bool active = false) :
			Vector<T,Dim>(values) {
		this->active_ = active;
	}

private:
	bool active_ = false;

};

}
}
