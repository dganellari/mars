#ifndef MARS_VECTOR_KOKKOS_HPP
#define MARS_VECTOR_KOKKOS_HPP

#include <cmath>
#include <iostream>
#include "mars_utils_kokkos.hpp"

namespace mars {

	template<typename T, Integer Dim>
	class TempArray {
	public:
	    T  values[Dim];
	    TempArray() {}

	    MARS_INLINE_FUNCTION friend TempArray operator*(const T &factor, const TempArray &v)
	    {
	    	TempArray ret;
	    	for(Integer i = 0; i < Dim; ++i) {
	    	    ret(i) = factor * v(i);
	    	}
	    	
	    	return ret;
	    }

	    MARS_INLINE_FUNCTION friend TempArray operator/(const TempArray &v, const T &factor)
	    {
	    	TempArray ret;
	    	for(Integer i = 0; i < Dim; ++i) {
	    	    ret(i) = v(i)/factor;
	    	}
	    	
	    	return ret;
	    }
	    
	    MARS_INLINE_FUNCTION TempArray operator-(const TempArray &right) const
	    {
	        TempArray ret;
	        for(Integer i = 0; i < Dim; ++i) {
	            ret(i) = (*this)(i) - right(i);
	        }
	        
	        return ret;
	    }
	    
	    MARS_INLINE_FUNCTION TempArray operator+(const TempArray &right) const
	    {
	        TempArray ret;
	        for(Integer i = 0; i < Dim; ++i) {
	            ret(i) = (*this)(i) + right(i);
	        }
	        
	        return ret;
	    }

	    MARS_INLINE_FUNCTION TempArray operator+=(const TempArray &right)
	    {
	        for(Integer i = 0; i < Dim; ++i) {
	            (*this)(i) += right(i);
	        }
	        
	        return *this;
	    }

	    MARS_INLINE_FUNCTION TempArray operator*=(const T &right)
	    {
	        for(Integer i = 0; i < Dim; ++i) {
	            (*this)(i) *= right;
	        }
	        
	        return *this;
	    }

	    MARS_INLINE_FUNCTION TempArray operator/=(const T &right)
	    {
	        for(Integer i = 0; i < Dim; ++i) {
	            (*this)(i) /= right;
	        }
	        
	        return *this;
	    }

	    MARS_INLINE_FUNCTION TempArray operator*(const TempArray &right) const
	    {
	        TempArray ret;
	        for(Integer i = 0; i < Dim; ++i) {
	            ret(i) = (*this)(i) * right(i);
	        }
	        
	        return ret;
	    }
	    
	    MARS_INLINE_FUNCTION T &operator()(const Integer i)
	    {
	        assert(i < Dim);
	        return values[i];
	    }
	    
	    MARS_INLINE_FUNCTION const T &operator()(const Integer i) const
	    {
	        assert(i < Dim);
	        return values[i];
	    }
	    
	    MARS_INLINE_FUNCTION T &operator[](const Integer i)
	    {
	        assert(i < Dim);
	        return values[i];
	    }
	    
	    MARS_INLINE_FUNCTION const T &operator[](const Integer i) const
	    {
	        assert(i < Dim);
	        return values[i];
	    }


	    MARS_INLINE_FUNCTION void describe(std::ostream &os) const
	    {
	        for(Integer i = 0; i < Dim; ++i) {
	            os << (*this)(i) << " ";
	        }

	        os << "\n";
	    }

	    MARS_INLINE_FUNCTION friend std::ostream &operator<<(std::ostream &os, const TempArray &v)
	    {
	        v.describe(os);
	        return os;
	    }

	    MARS_INLINE_FUNCTION T squared_norm() const
	    {
	        T sqn = (*this)(0) * (*this)(0);
	        for(Integer i = 1; i < Dim; ++i)
	        {
	            sqn += (*this)(i) * (*this)(i);
	        }

	        return sqn;
	    }

	    MARS_INLINE_FUNCTION T norm() const
	    {
	        return std::sqrt(squared_norm());
	    }

	    MARS_INLINE_FUNCTION TempArray &normalize() {
	        T len = norm();

	        for(Integer i = 0; i < Dim; ++i) {
	            (*this)(i) /= len;
	        }

	        return *this;
	    }

	    MARS_INLINE_FUNCTION TempArray &zero()
	    {
	        for(Integer i = 0; i < Dim; ++i) {
	           (*this)(i) = 0;
	        }

	    	return *this;
	    }

	    MARS_INLINE_FUNCTION TempArray &set(const Real value)
		{
			for (Integer i = 0; i < Dim; ++i) {
				(*this)(i) = value;
			}

			return *this;
		}
	};

	template<typename T, Integer Dim>
	MARS_INLINE_FUNCTION T dot(const TempArray<T, Dim> &left, const TempArray<T, Dim> &right)
	{
		T ret = 0.;
		for(Integer d = 0; d < Dim; ++d) {
			ret += left(d) * right(d);
		}

		return ret;
	}
}

#endif //MARS_VECTOR_HPP
