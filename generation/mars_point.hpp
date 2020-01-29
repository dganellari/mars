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
	}

	MARS_INLINE_FUNCTION Point(const ViewMatrixType<T>& v, Integer id) :
			Point(SubView<T, Dim>(&v, id))
	{
	}


	MARS_INLINE_FUNCTION
	T &operator[](const Integer i)
	{
		assert(i < Dim);
		return this->pointView[i];
	}

	MARS_INLINE_FUNCTION
	const T &operator[](const Integer i) const
	{
		assert(i < Dim);
		return this->pointView[i];
	}

	MARS_INLINE_FUNCTION
	TempArray<T,Dim> operator-(const Point &right) const
	{
		TempArray<T,Dim> ret;
		for (Integer i = 0; i < Dim; ++i)
		{
			ret[i] = this->operator[](i) - right[i];
		}

		return ret;
	}

	MARS_INLINE_FUNCTION
	TempArray<T,Dim> operator-(const TempArray<T,Dim>  &right) const
	{
		TempArray<T,Dim> ret;
		for (Integer i = 0; i < Dim; ++i)
		{
			ret[i] = this->operator[](i) - right[i];
		}

		return ret;
	}

	MARS_INLINE_FUNCTION
	TempArray<T,Dim> operator-(const ViewVectorTypeC<T, Dim>& right) const
	{
		TempArray<T,Dim> ret;
		for (Integer i = 0; i < Dim; ++i)
		{
			ret[i] = this->operator[](i) - right(i);
		}

		return ret;
	}

	MARS_INLINE_FUNCTION
	TempArray<T,Dim> operator+(const Point &right) const
	{
		TempArray<T,Dim> ret;
		for (Integer i = 0; i < Dim; ++i)
		{
			ret[i] = this->operator[](i) + right[i];
		}

		return ret;
	}


	MARS_INLINE_FUNCTION
	TempArray<T,Dim> operator*(const Point &right) const
	{
		TempArray<T,Dim> ret;
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

	//just for the sake of the stable select. Stable on the coordinates comp.
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

	/*MARS_INLINE_FUNCTION
	T squared_norm(const Point& right) const
	{
		T tmp[Dim];

		for (Integer i = 0; i < Dim; ++i)
		{
			tmp[i] = this->operator[](i) - right[i];
		}

		T sqn = tmp[0] * tmp[0];
		for (Integer i = 1; i < Dim; ++i)
		{
			sqn += tmp[i] * tmp[i];
		}

		return sqn;
	}*/

	MARS_INLINE_FUNCTION
	T squared_distance(const Point& right) const
	{
		T sqn = (this->operator[](0) - right[0])
				* (this->operator[](0) - right[0]);
		for (Integer i = 1; i < Dim; ++i)
		{
			sqn += (this->operator[](i) - right[i])
					* (this->operator[](i) - right[i]);
		}
		return sqn;
	}

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


private:
	SubView<T, Dim> pointView;
};

}
