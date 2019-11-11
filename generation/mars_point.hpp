#pragma once

#include "mars_device_vector.hpp"
#include "mars_sub_view.hpp"

namespace mars
{

template<typename T, Integer Dim>
class Point
{

public:

	MARS_INLINE_FUNCTION Point()
	{

	}

	MARS_INLINE_FUNCTION Point(SubView<T, Dim> pv) :
			pointView(pv)
	{
		temp_ = false;
	}

	MARS_INLINE_FUNCTION Point(ViewMatrixType<T> v, Integer id) :
			Point(SubView<T, Dim>(v, id))
	{

	}

	MARS_INLINE_FUNCTION Point(ViewVectorTypeC<T, Dim> v)
	{
		for(Integer i=0; i<Dim; ++i){
			pointArray[i] = v(i);
		}
	}

	MARS_INLINE_FUNCTION
	T &operator[](const Integer i)
	{
		assert(i < Dim);

		if (temp_)
			return this->pointArray[i];

		return this->pointView[i];
	}

	MARS_INLINE_FUNCTION
	const T &operator[](const Integer i) const
	{
		assert(i < Dim);

		if (temp_)
			return this->pointArray[i];

		return this->pointView[i];
	}

	MARS_INLINE_FUNCTION
	Point operator-(const Point &right) const
	{
		Point ret;
		for (Integer i = 0; i < Dim; ++i)
		{
			ret[i] = this->operator[](i) - right[i];
		}

		return ret;
	}

	MARS_INLINE_FUNCTION
	Point operator+(const Point &right) const
	{
		Point ret;
		for (Integer i = 0; i < Dim; ++i)
		{
			ret[i] = this->operator[](i) + right[i];

		}

		return ret;
	}

	MARS_INLINE_FUNCTION
	Point operator*(const Point &right) const
	{
		Point ret;
		for (Integer i = 0; i < Dim; ++i)
		{
			ret[i] = this->operator[](i) * right[i];
		}

		return ret;
	}

	MARS_INLINE_FUNCTION
	const Point &operator/(const T divident)
	{
		for (Integer i = 0; i < Dim; ++i)
		{
			this->operator[](i) /= divident;
		}

		return *this;
	}

	//just for the purpose of the stable select. Stable on the coordinates comp.
	MARS_INLINE_FUNCTION
	bool operator<(const Point &right) const
	{
		for (Integer i = 0; i < Dim; ++i)
		{
			if(this->operator[](i) < right[i])
				return true;
			if(this->operator[](i) > right[i])
				return false;
		}

		return false;
	}

	/*
	 void describe(std::ostream &os) const
	 {
	 for (Integer i = 0; i < Dim; ++i)
	 {
	 os << this->view(this->index, i) << " ";
	 }

	 os << "\n";
	 }

	 friend std::ostream &operator<<(std::ostream &os, const SubView &v)
	 {
	 v.describe(os);
	 return os;
	 }*/

	MARS_INLINE_FUNCTION
	T squared_norm() const
	{
		T sqn = this->operator[](0) * this->operator[](0);
		for (Integer i = 1; i < Dim; ++i)
		{
			sqn += this->operator[](i) * this->operator[](i);
		}

		return sqn;
	}

	MARS_INLINE_FUNCTION
	T norm() const
	{
		return std::sqrt(squared_norm());
	}

	MARS_INLINE_FUNCTION
	Point &normalize()
	{
		T len = norm();

		for (Integer i = 0; i < Dim; ++i)
		{
			this->operator[](i) /= len;
		}

		return *this;
	}

	bool is_temp() const
	{
		return temp_;
	}

	void set_temp(bool temp = false)
	{
		temp_ = temp;
	}

private:
	SubView<T, Dim> pointView;
	TempArray<T, Dim> pointArray;
	bool temp_ = true;

};

}
