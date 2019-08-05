#ifndef MARS_VECTOR_HPP
#define MARS_VECTOR_HPP
#include "mars_tensor_base.hpp"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Written by Gabriele Rovi (July 2019)                                                                    ////////
////// We define the Vector class:                                                                             ////////
////// Vector<typename T, Integer Dim>                                                                         ////////
////// We want vectors to be eventually constexpr and static. A proper constructor is then needed              ////////
////// To build it, we must define a VectorBase class from which Vector inherit                                ////////
////// In this way we can use Vector<T,Dim> as we would normally do, but we can also define:                   ////////
////// static constexpr Vector<T,Dim> static_mat{T1,T2....};                                                   ////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace mars {


	template<typename T, Integer Dim_>
	class Vector: public TensorBase<T, std::make_index_sequence<Dim_>> 
	{
	public:
		static constexpr Integer Dim=Dim_;
		using type=Vector<T,Dim>;
	    using subtype=T;
		using MB = TensorBase<T, std::make_index_sequence<Dim>>;
		using MB::MB;
		using MB::values;


	    friend constexpr Vector operator*(const Real &alpha, const Vector &v)
	    {
	    	Vector ret;
	    	for(Integer i = 0; i < Dim; ++i) {
	    	    ret(i) = alpha * v(i);
	    	}


	    	return ret;
	    }

	    friend constexpr Vector operator/(const Vector &v, const Real &alpha)
	    {
	    	Vector ret;
	    	for(Integer i = 0; i < Dim; ++i) {
	    	    ret(i) = v(i)/alpha;
	    	}
	    	
	    	return ret;
	    }

	    constexpr Vector operator-() const
	    {
	        Vector ret;
	        for(Integer i = 0; i < Dim; ++i) {
	            ret(i) = -(*this)(i) ;
	        }
	        
	        return ret;
	    }

	    constexpr Vector operator-(const Vector &right) const
	    {
	        Vector ret;
	        for(Integer i = 0; i < Dim; ++i) {
	            ret(i) = (*this)(i) - right(i);
	        }
	        
	        return ret;
	    }
	    
	    constexpr Vector operator+(const Vector &right) const
	    {
	        Vector ret;
	        for(Integer i = 0; i < Dim; ++i) {
	            ret(i) = (*this)(i) + right(i);
	        }
	        
	        return ret;
	    }


         

	     constexpr Vector& operator=(const Vector &right)
	    {
	        for(Integer i = 0; i < Dim; ++i) {
	            // const_cast<T&>(static_cast<const std::array<T,Dim>& >((*this)())[i] )=right(i);
	            (*this)()[i]=right(i);
	        }
	        
	        return *this;
	    }

	    
	    constexpr Vector &operator=(const T &value) const
	    {
	        for(Integer i = 0; i < Dim; ++i) {
	            (*this)(i) = value;
	        }
	        
	        return *this;
	    }

	    constexpr Vector &operator+=(const Vector &right)
	    {
	        for(Integer i = 0; i < Dim; ++i) {
	            (*this)(i) += right(i);
	        }
	        
	        return *this;
	    }

	    constexpr Vector &operator-=(const Vector &right)
	    {
	        for(Integer i = 0; i < Dim; ++i) {
	            (*this)(i) -= right(i);
	        }
	        
	        return *this;
	    }

	    constexpr Vector &operator*=(const T &right)
	    {
	        for(Integer i = 0; i < Dim; ++i) {
	            (*this)(i) *= right;
	        }
	        
	        return *this;
	    }

	    constexpr Vector &operator/=(const T &right)
	    {
	        for(Integer i = 0; i < Dim; ++i) {
	            (*this)(i) /= right;
	        }
	        
	        return *this;
	    }

	    constexpr Vector operator*(const Vector &right) const
	    {
	        Vector ret;
	        for(Integer i = 0; i < Dim; ++i) {
	            ret(i) = (*this)(i) * right(i);
	        }
	        
	        return ret;
	    }
	
		inline constexpr std::array<T,Dim> &operator()()
		{
			return values;
		}

		inline constexpr const std::array<T,Dim> &operator()()const
		{
			return values;
		}

	    inline constexpr T &operator()(const Integer i)
	    {
	        assert(i < Dim);
	        return values[i];
	    }
	    
	    inline constexpr const T &operator()(const Integer i) const
	    {
	        assert(i < Dim);
	        return values[i];
	    }
	    
	    inline constexpr T &operator[](const Integer i)
	    {
	        assert(i < Dim);
	        return values[i];
	    }
	    
	    inline constexpr const T &operator[](const Integer i) const
	    {
	        assert(i < Dim);
	        return values[i];
	    }



	    void describe(std::ostream &os) const
	    {
	        for(Integer i = 0; i < Dim; ++i) {
	            os << (*this)(i) << " ";
	        }

	        os << "\n";
	    }

	    friend std::ostream &operator<<(std::ostream &os, const Vector &v)
	    {
	        v.describe(os);
	        return os;
	    }

	    inline constexpr T squared_norm() const
	    {
	        T sqn = (*this)(0) * (*this)(0);
	        for(Integer i = 1; i < Dim; ++i)
	        {
	            sqn += (*this)(i) * (*this)(i);
	        }

	        return sqn;
	    }

	    inline constexpr T norm() const
	    {
	        return std::sqrt(squared_norm());
	    }

	    constexpr Vector &normalize() {
	        T len = norm();

	        for(Integer i = 0; i < Dim; ++i) {
	            (*this)(i) /= len;
	        }

	        return *this;
	    }

	    inline constexpr Vector &zero()
	    {
	        std::fill(begin(values), end(values), 0.);
	        return *this;
	    }

	    inline constexpr Vector &set(const Real value)
	    {
	    	std::fill(begin(values), end(values), value);
	    	return *this;
	    }



	    inline constexpr T sum() const
	    {
	        T s=(*this)(0);
	        for(Integer i = 1; i < Dim; ++i) {
	            s+=(*this)(i);
	        }
	        return s;
	    }
        
        // this works if T is Real, vector, matrix.
        // It does not for integers, in general (it should return a Real)
        // this is why we call it Tmean and not mean
	    inline constexpr const T Tmean() 
	    {
	        return sum()/Dim;
 	    }
        
        inline constexpr const Integer size()const {return Dim;};


	};

	template<typename T, Integer Dim>
	inline constexpr T dot(const Vector<T, Dim> &left, const Vector<T, Dim> &right)
	{
		T ret = 0.;
		for(Integer d = 0; d < Dim; ++d) {
			ret += left(d) * right(d);
		}

		return ret;
	}



}

#endif //MARS_VECTOR_HPP
