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

// <<<<<<< HEAD
// =======
	    // Vector(std::initializer_list<T> values)
	    // {
	    // 	assert(values.size() == Dim);
	    //     std::copy(std::begin(values), std::begin(values) + Dim, std::begin(this->values));
	    // }
// >>>>>>> remotes/origin/master

	    friend constexpr Vector operator*(const Real &alpha, const Vector &v)
	    {
	    	Vector ret;
	    	for(Integer i = 0; i < Dim; ++i) {
	    	    ret[i] = alpha * v[i];
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


        static constexpr Vector eye_vector(const T& t){
        	Vector tmp;
        	for(Integer i=0; i<Dim;i++)
        		tmp[i]=t;
        	return tmp;
        }

        constexpr bool is_zero(){
        	Vector tmp;
        	for(Integer i=0; i<Dim;i++)
        	  if((*this)[i]!=0)
        		return false;
        	return true;
        }

	};

	template<typename T, Integer Dim>
	inline constexpr T dot(const Vector<T, Dim> &left, const Vector<T, Dim> &right)
	{
		T ret = left(0) * right(0);
		for(Integer d = 1; d < Dim; ++d) {
			ret += left(d) * right(d);
		}

		return ret;
	}



 template<Integer Dim>
 bool are_parallel_vectors(const Vector<Real,Dim>& p_a,const Vector<Real,Dim>& p_b,const Vector<Real,Dim>& q_a,const Vector<Real,Dim>& q_b)
 {
 	Real k;
 	Real toll=1e-8;
 	Real a_dot_b=(p_a[0]-p_b[0])*(q_a[0]-q_b[0]);
 	Real a_norm=(p_a[0]-p_b[0])*(p_a[0]-p_b[0]);
 	Real b_norm=(q_a[0]-q_b[0])*(q_a[0]-q_b[0]);

 	for(Integer k=1;k<Dim;k++)
 	{
 		a_dot_b+=(p_a[k]-p_b[k])*(q_a[k]-q_b[k]);
 		a_norm+=(p_a[k]-p_b[k])*(p_a[k]-p_b[k]);
 		b_norm+=(q_a[k]-q_b[k])*(q_a[k]-q_b[k]);
 	}

 	// std::cout<<"a_dot_b="<<a_dot_b<<std::endl;
 	// std::cout<<"a_norm="<<a_norm<<std::endl;
 	// std::cout<<"b_norm="<<b_norm<<std::endl;



 	if(abs(a_dot_b*a_dot_b-a_norm*b_norm)<toll)
 		return true;
 	else
 		return false;

	// if(abs(p_a[0]-p_b[0])>toll)
	// {
	// 	k=(q_a[0]-q_b[0])/(p_a[0]-p_b[0]);
	// }
	// else
	// {
	// 	// if p_a[0]-p_b[0]=0 but q_a[0]-q_b[0]!=0, then they are not parallel
	// 	if(abs(q_a[0]-q_b[0])>toll)
	// 		return false;
	// 	else
	// 		k=0;
	// }

 // 	for(Integer i=1;i<Dim;i++)
 // 	{
 // 		if(k==0)
 // 		{
	// 		if(abs(p_a[i]-p_b[i])>toll)
	// 		{
	// 			k=(q_a[i]-q_b[i])/(p_a[i]-p_b[i]);
	// 		}
	// 		else
	// 		{
	// 			// if p_a[0]-p_b[0]=0 but q_a[0]-q_b[0]!=0, then they are not parallel
	// 			if(abs(q_a[0]-q_b[0])>toll)
	// 				return false;
	// 			else
	// 				k=0;
	// 		}

 // 		}
 // 		else
 // 		{
 // 			// k!= 0, but p_a[i]-p_b[i]==0
 // 		 	if(abs(p_a[i]-p_b[i])<toll)
	// 		{
	// 			return false;
	// 		}	
	// 		else 
	// 		{
	// 			if(abs(k-(q_a[i]-q_b[i])/(p_a[i]-p_b[i]))>toll )
	// 				return false;
	// 		}



 // 		}

 // 	}
 // 	// if k==0 they are both null vectors, so they are parallel
 // 	// if k!=0 and we have not exited yet, they are parallel
 //     return true;
 }

  template<Integer Dim>
  bool does_point_belong_to_segment(const Vector<Real,Dim>& q_a,const Vector<Real,Dim>& q_b,const Vector<Real,Dim>& w)
  {

  	Real toll=1e-8;
  	// for q_b-q_a to be a subsegment of p_b-p_a
  	// it must occur that in the line defined by
	// p= p_a + alpha (p_b-p_a)
  	// if we substitue q_a and q_b to p, the equation must be satisfied with alpha in [0,1]
    // since we know the two vectors are parallel (p_b-p_a)= k (q_b-q_a)
    Real alpha;
    Real alpha_tmp;
    bool found_non_zero=false;
    // std::cout<<"---------------"<<std::endl;

    // std::cout<<"q_a"<<std::endl;
    // std::cout<<q_a<<std::endl;
    // std::cout<<"q_b"<<std::endl;
    // std::cout<<q_b<<std::endl;
    // std::cout<<"w"<<std::endl;
    // std::cout<<w<<std::endl;

    for(Integer i=0;i<Dim;i++)
    {
    	// std::cout<<"i = "<< i<<std::endl;
        // std::cout<<"abs((q_b[i]-q_a[i]))="<< abs((q_b[i]-q_a[i]))<<std::endl;
    	if(!found_non_zero)
    	{
			if(abs((q_b[i]-q_a[i]))>toll)//abs(w[i]-q_a[i])>toll)
			{
				alpha=(w[i]-q_a[i])/(q_b[i]-q_a[i]);
				found_non_zero=true;
				// std::cout<<"if(abs(w[i]-q_a[i])>toll) alpha"<< alpha<<std::endl;
				if(alpha<-toll || alpha>1+toll)
					return false;
			}
			else
			{
				// std::cout<<"else"<<std::endl;
				if(abs((w[i]-q_a[i]))>toll)
					return false;
				// found_non_zero is still false
			}

    	}
    	else
    	{
    		if(abs(q_b[i]-q_a[i])>toll)
    		  {
			  	alpha_tmp=(w[i]-q_a[i])/(q_b[i]-q_a[i]);
			  	// std::cout<<"alpha_tmp"<< alpha_tmp<<std::endl;
			
				if(abs(alpha-alpha_tmp)>toll || alpha_tmp<-toll || alpha_tmp>1+toll)
					return false;
    		  }
		    else
    		 {
    		 	// std::cout<<"q_b[i]-q_a[i]"<< q_b[i]-q_a[i]<<std::endl;
    		 	if(abs(w[i]-q_a[i])>toll)
    		 		return false;
    		 	else
    		 	{
    		 		//  q_b[i]-q_a[i]=w[i]-q_a[i] =0
    		 	}

    		 }



    	}

    // std::cout<<" alpha = "<<alpha<<std::endl;

    }  
    // std::cout<< "return true"<<std::endl;	
    return true;
  }





  template<Integer Dim>
  bool is_subsegment(const Vector<Real,Dim>& p_a,const Vector<Real,Dim>& p_b,const Vector<Real,Dim>& q_a,const Vector<Real,Dim>& q_b)
  {
    // if the two vectors are not parallel, for sure p_b-p_a cannot be subsegment of q_b-q_a
  	// std::cout<<"are_parallel_vectors(p_a,p_b,q_a,q_b)"<<std::endl;
  	// std::cout<<are_parallel_vectors(p_a,p_b,q_a,q_b)<<std::endl;

  	if(!are_parallel_vectors(p_a,p_b,q_a,q_b))
  	   return false;
  	// std::cout<<"are_parallel_vectors(p_a,p_b,q_a,q_b)"<<std::endl;
  	std::cout<<are_parallel_vectors(p_a,p_b,q_a,q_b)<<std::endl;

  	// std::cout<<"does_point_belong_to_segment(p_a,p_b,q_a)"<<std::endl;
  	bool b1=does_point_belong_to_segment(p_a,p_b,q_a);
  	// std::cout<<b1<<std::endl;
  	// std::cout<<"does_point_belong_to_segment(p_a,p_b,q_b)"<<std::endl;
  	bool b2=does_point_belong_to_segment(p_a,p_b,q_b);
  	std::cout<<b2<<std::endl;

    return b1*b2;
    // return does_point_belong_to_segment(p_a,p_b,q_a)*does_point_belong_to_segment(p_a,p_b,q_b);
  }

  template<Integer Dim>
  Real segment_length(const Vector<Real,Dim>& p_a,const Vector<Real,Dim>& p_b)
  {
  	Real length;
  	length= (p_a[0]-p_b[0])*(p_a[0]-p_b[0]);
  	for(Integer i=1;i<Dim;i++)
  		length+=(p_a[i]-p_b[i])*(p_a[i]-p_b[i]);

  	return sqrt(length);

  }

  template<Integer Dim>
  Real segment_fraction_length(const Vector<Real,Dim>& p_a,const Vector<Real,Dim>& p_b,const Vector<Real,Dim>& q_a,const Vector<Real,Dim>& q_b)
  {

  	return segment_length(q_a,q_b)/segment_length(p_a,p_b);

  }


  template<Integer Dim>
  bool are_vectors_equal(const Vector<Real,Dim>& a,const Vector<Real,Dim>& b)
  {
  	Real toll=1e-8;

    for(Integer i=0;i<Dim;i++)
    {
    	if(abs(a[i]-b[i])>toll)
    		return false;

    }
  	return true;

  }


}

#endif //MARS_VECTOR_HPP
