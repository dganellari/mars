#ifndef MARS_I_MESH_HPP
#define MARS_I_MESH_HPP

#include "mars_vector.hpp"

#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <memory>
#include <algorithm>  

namespace mars {

	class IElem {
	public:
		virtual ~IElem() {}
		virtual void get_nodes(std::vector<Integer> &nodes) const = 0;
		virtual Integer n_nodes() const = 0;
		virtual Integer node(const Integer idx) const = 0;
        virtual Integer get_block() const = 0;
        virtual void set_block(const Integer block_id) = 0;

		//for the moment it just returns the simplex type (e.g. 2 for triangles)
		virtual Integer type() const = 0;
	};

	inline Integer n_nodes(const IElem &elem)
	{
		return elem.n_nodes();
	}

	template<Integer Dim_>
	class IMesh {
	public:
		static const Integer Dim = Dim_;
		using Point = mars::Vector<Real, Dim>;
		using Elem  = mars::IElem;

		virtual ~IMesh() {}
		virtual void points(const Integer id, std::vector<Point> &pts) const = 0;
		virtual IElem &elem(const Integer id) = 0;
		virtual const IElem &elem(const Integer id) const = 0;
		virtual bool is_active(const Integer id) const = 0;
		virtual Integer n_nodes() const = 0;
		virtual Integer n_elements() const = 0;
		virtual Integer n_active_elements() const = 0;
		virtual Integer add_point(const Point &point) = 0;
		virtual Point &point(const Integer i) = 0;		
		virtual const Point &point(const Integer i) const = 0;
		virtual Integer add_elem(const IElem &elem) = 0;
		virtual Integer add_elem(const std::vector<Integer> &nodes) = 0;
		// virtual const std::vector<Point> &points() const = 0;

		virtual Integer type() const = 0;
		virtual void reserve(
			const std::size_t n_elements,
			const std::size_t n_points) = 0;
	};

}

#endif //MARS_I_MESH_HPP
