#ifndef MARS_Array_HPP
#define MARS_Array_HPP
#include "mars_tensor_base.hpp"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Written by Gabriele Rovi (July 2019)                                                                    ////////
////// We define the Array class:                                                                              ////////
////// Array<typename T, Integer Dim>                                                                          ////////
////// We want Arrays to be eventually constexpr and static. A proper constructor is then needed               ////////
////// To build it, we must define a ArrayBase class from which Array inherit                                  ////////
////// In this way we can use Array<T,Dim> as we would normally do, but we can also define:                    ////////
////// static constexpr Array<T,Dim> static_mat{T1,T2....};                                                    ////////
////// With respect to Vector, Array does not provide +,-,*,/ operators                                       /////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace mars {


    template<typename T, Integer Dim_>
    class Array: public TensorBase<T, std::make_index_sequence<Dim_>> 
    {
    public:
        static constexpr Integer Dim=Dim_;
        using type=Array<T,Dim>;
        using value_type=T;
        using size_type= Integer;
        using difference_type=std::ptrdiff_t;
        using reference=value_type&;
        using const_reference=const value_type&;
        using pointer=value_type*;
        using const_pointer=const value_type*;


        using MB = TensorBase<T, std::make_index_sequence<Dim>>;
        using MB::MB;
        using MB::values;

        constexpr Array& operator=(const std::vector<T> &right)
        {
            for(Integer i = 0; i < Dim; ++i) {
                // const_cast<T&>(static_cast<const std::array<T,Dim>& >((*this)())[i] )=right(i);
                (*this)()[i]=right[i];
            }
            
            return *this;
        } 

        constexpr Array& operator=(const std::array<T,Dim> &right)
        {
            for(Integer i = 0; i < Dim; ++i) {
                // const_cast<T&>(static_cast<const std::array<T,Dim>& >((*this)())[i] )=right(i);
                (*this)()[i]=right[i];
            }
            
            return *this;
        } 

        constexpr Array& operator=(const Array &right)
        {
            for(Integer i = 0; i < Dim; ++i) {
                // const_cast<T&>(static_cast<const std::array<T,Dim>& >((*this)())[i] )=right(i);
                (*this)()[i]=right(i);
            }
            
            return *this;
        }

        
        constexpr Array &operator=(const T &value) const
        {
            for(Integer i = 0; i < Dim; ++i) {
                (*this)(i) = value;
            }
            
            return *this;
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

        inline constexpr const T & front()const 
        {
           assert(0 < Dim);
           return values[0];           
       }

       inline constexpr T & front() 
       {
           assert(0 < Dim);
           return values[0];           
       }

       inline constexpr const T & back()const 
       {
           assert(0 < Dim);
           return values[Dim-1];           
       }

       inline constexpr T & back() 
       {
           assert(0 < Dim);
           return values[Dim-1];           
       }

       inline constexpr T& at(const Integer i )
       {
           assert(i < Dim);
           return values[i];
       }

       inline constexpr const T& at(const Integer i ) const
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

    friend std::ostream &operator<<(std::ostream &os, const Array &v)
    {
        v.describe(os);
        return os;
    }


    inline constexpr Array &zero()
    {
        std::fill(begin(values), end(values), 0.);
        return *this;
    }

    inline constexpr Array &set(const Real value)
    {
        std::fill(begin(values), end(values), value);
        return *this;
    }

    inline constexpr void fill(const Real value)
    {
       for(Integer ii=0;ii<Dim;ii++)
           values[ii]=value;           
   } 

   inline constexpr const Integer size()const {return Dim;};

};


template<typename S, typename T,Integer Dim1,Integer Dim2>
constexpr auto subarray(const Array<S,Dim1>&v,const Array<T,Dim2>&w)
{
  Array<S,Dim2> u;
  assert((Dim1>=Dim2)&&" cannot extract a subarray with components not belonging to the original one");
  for(Integer ii=0;ii<Dim2;ii++)
    u[ii]=v[w[ii]];
  return u;

}

//     template<Integer N,Integer K>
// constexpr void combinations_generate_aux(
//     Array<Integer, K> &data,
//     const Integer index, 
//     const Integer i,
//     Array<Array<Integer, K>, binomial_coefficient(N,K)> &combs,
//     Integer &comb_index)
// {
//     if(index == K) {
//         for(Integer ii=0;ii<data.size();ii++)
//             combs[comb_index][ii]=data[ii];
//         comb_index++;
//         return;
//     }

//     if(i >= N) {
//         return;
//     }

//     data[index] = i;

//     combinations_generate_aux<N,K>(data, index+1, i+1, combs, comb_index);

//             // current is excluded, replace it with next (Note that
//             // i+1 is passed, but index is not changed)
//     combinations_generate_aux<N,K>(data, index, i+1, combs, comb_index);
// }

//         template<Integer N,Integer K >
// constexpr Array<Array<Integer, K>, binomial_coefficient(N,K)> combinations_generate()
// {
//     Array<Array<Integer, K>, binomial_coefficient(N,K)> combs;
//     Array<Integer, K> data;
//     Integer comb_index = 0;
//     combinations_generate_aux<N,K>(data, 0, 0, combs, comb_index);
//     return combs;
// }

//     template<Integer Dim, Integer ManifoldDim,Integer Npoints>
// inline constexpr auto jacobian_faces(Integer face)
// {
//     static_assert(Dim >= ManifoldDim, "Dim must be greater or equal ManifoldDim");
//     static_assert(Npoints == ManifoldDim+1, "Npoints must be equal to ManifoldDim+1");
//     Array<Array<Real, Dim>,Npoints> points;
//     const auto &n_faces=Npoints;
//     const auto combs=combinations_generate<ManifoldDim+1,ManifoldDim>(); 
//     Matrix<Real, Dim, ManifoldDim> J;
//     Array<Matrix<Real, Dim, ManifoldDim>,n_faces> Jmat;
//     // loop on all the faces
//     for(Integer ii=0;ii<n_faces;ii++)
//     {
//         // take the indices of the reference simplex related to the face ii
//         const auto &comb_ii=combs[ii];
//         // fill points with the corresponding reference simplex face 
//         for(Integer jj=0;jj<ManifoldDim;jj++)
//           for(Integer kk=0;kk<ManifoldDim;kk++)
//             points[jj][kk]=Simplex<Dim,ManifoldDim>::reference[comb_ii[jj]][kk];
        
//         // compute the jacobian of the given face
//         for(Integer ii=0;ii<Npoints;ii++)
//         {
//             Array<Real, Dim> v0 = points[0];
//             for(Integer i = 1; i < Npoints; ++i) {
//                 const auto &vi = points[i];
//                 J.col(i-1, vi - v0);
//             }
//         }
//         Jmat[ii]=J;

//     }
//     return Jmat;
// }

}

#endif //MARS_Array_HPP
