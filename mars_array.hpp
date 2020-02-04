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

template<typename T>
void subvector(std::vector<T>& u,const std::vector<T>&v,const std::vector<Integer>&w)
{
  Integer size=w.size();
  u.resize(size);
  for(Integer ii=0;ii<size;ii++)
    u[ii]=v[w[ii]];
}

template<typename S, typename T,Integer Dim1,Integer Dim2>
constexpr auto subarray(const Array<S,Dim1>&v,const Array<T,Dim2>&w)
{
  Array<S,Dim2> u;
  assert((Dim1>=Dim2)&&" cannot extract a subarray with components not belonging to the original one");
  for(Integer ii=0;ii<Dim2;ii++)
    u[ii]=v[w[ii]];
  return u;

}


template<typename S, typename T,Integer Dim1,Integer Dim2>
constexpr auto subarray(Array<S,Dim2>&u , const Array<S,Dim1>&v,const Array<T,Dim2>&w)
{

  assert((Dim1>=Dim2)&&" cannot extract a subarray with components not belonging to the original one");
  for(Integer ii=0;ii<Dim2;ii++)
    u[ii]=v[w[ii]];

}

template<typename S, typename T,Integer Dim>
constexpr auto subarray(Array<S,Dim>&u , const std::vector<S>&v,const Array<T,Dim>&w)
{

  // assert((v.size()<Dim)&&" cannot extract a subarray with components not belonging to the std::Vector");
  for(Integer ii=0;ii<Dim;ii++)
    u[ii]=v[w[ii]];

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


template<Integer N,Integer M>
constexpr auto cumulative_array
(const Array<Integer, N>& t1,const Array<Integer, M>& t2)
{
   Array<Integer, N+1> t3;

   for(Integer i=0;i<t2[0];i++)
      {
       t3[i]= 0;
      }
   Integer cont=t2[0];
   Integer dofs_count=t1[0];
   for(Integer i=1;i<t2.size();i++)
      {
       for(Integer j=0;j<t2[i];j++)
          {
           t3[cont]=dofs_count;
           cont++;}

        dofs_count+=t1[cont-1];
      }

  return t3;
}

template<Integer N,Integer M>
auto cumulative_array
(const Array<std::vector<Integer>, N>& t1,const Array<Integer, M>& t2)
{

// std::cout<<"N==="<<N<<std::endl;
//   std::cout<<"<<<<<<<t1>>>>>>>>>>>2"<<std::endl;
//    for(Integer i=0;i<t1.size();i++)
//    {
//     for(Integer j=0;j<t1[i].size();j++)
//     std::cout<<t1[i][j]<<" ";
// std::cout<<std::endl;
    
//    }

  // std::cout<<"<<<<<<<t2>>>>>>>>>>>2"<<std::endl;
  //  for(Integer i=0;i<t2.size();i++)
  //  {
  //   std::cout<<t2[i]<<" ";
  //  }
  // std::cout<<std::endl;


  // std::cout<<"<<<<<<<initialize t3 with zeros>>>>>>>>>>>1"<<std::endl;
  
   auto levels=t1[0].size();

   Array<std::vector<Integer>, N+1> t3;
   // std::cout<<"levels="<<levels<<std::endl;
   for(Integer i=0;i<N+1;i++)
      {
       t3[i].resize(levels, 0);
       // for(Integer j=0;j<levels;j++)
       // std::cout<<t3[i][j]<<" ";
       // std::cout<<std::endl;
      }


  // std::cout<<"<<<<<<<cumulative_array>>>>>>>>>>>2"<<std::endl;
   Integer cont;
   cont=t2[0];
   std::vector<Integer> dofs_count(levels);
   
   for(Integer i=0;i<levels;i++)
   {
    // std::cout<<i<<" "<<levels<<std::endl;
    dofs_count[i]=t1[0][i];
    // std::cout<<dofs_count[i]<<" ";
    
   }
   // std::cout<<std::endl;
   // std::cout<<"cont="<<cont<<std::endl;
  // std::cout<<"<<<<<<<cumulative_array>>>>>>>>>>>3 t2size="<<t2.size()<<std::endl;
   for(Integer i=1;i<M;i++)
      {
        // std::cout<<"t2[i]"<<t2[i]<<std::endl;

       for(Integer j=0;j<t2[i];j++)
          {
            // std::cout<<"j="<<j<<std::endl;
            // t3[cont].resize(levels,0);
            // std::cout<<"levels="<<levels<<std::endl;
             for(Integer s=0;s<levels;s++)
                { 
                  // std::cout<<cont<<" "<<s<<" ";
                  // std::cout<<"cont="<<cont<<std::endl;
                  // std::cout<<"s="<<s<<std::endl;
                  // std::cout<<"dofs_count[s]="<<dofs_count[s]<<std::endl;
                  t3[cont][s]=dofs_count[s];
                  // std::cout<<t3[cont][s]<<" "<<std::endl;
                  // std::cout<<"t3[cont][s]="<<t3[cont][s]<<std::endl;

                          // for(Integer j=0;j<t3[cont].size();j++)
                          //  std::cout<<t3[cont][j]<<" "<<std::endl;
                          //  std::cout<<std::endl;

                }
              cont++;

                  // std::cout<<"AFTER I T3"<<std::endl;
                  //        for(Integer h=0;h<t3.size();h++)
                  //        {
                  //         for(Integer j=0;j<t3[h].size();j++)
                  //          std::cout<<t3[h][j]<<" "<<std::endl;
                  //          std::cout<<std::endl;
                  //         }



        }

        // if(i<M-1)
        // {
        for(Integer s=0;s<levels;s++)
        dofs_count[s]+=t1[cont-1][s];
        
        // std::cout<<std::endl;
        // for(Integer s=0;s<levels;s++)
        //   std::cout<<t1[cont-1][s]<<" ";
        // std::cout<<std::endl;
        // }


      }


      for(Integer s=0;s<levels;s++)
          t3[cont][s]=dofs_count[s];   
  // std::cout<<"<<<<<<<cumulative_array>>>>>>>>>>>4"<<std::endl;


  // std::cout<<"<<<<<<<t3>>>>>>>>>>>"<<std::endl;
//    for(Integer i=0;i<t3.size();i++)
//    {
//     for(Integer j=0;j<t3[i].size();j++)
//     std::cout<<t3[i][j]<<" ";
// std::cout<<std::endl;
    
//    }
  return t3;
}


template<Integer N,Integer M>
constexpr auto cumulative_array_and_zero
(const Array<Integer, N>& t1,const Array<Integer, M>& t2)
{

  Array<Integer, N+M> t3;

  for(Integer ii=0;ii<N;ii++)
       t3[ii]=t1[ii];

  for(Integer ii=N-1;ii<N+M+1;ii++)
       t3[ii+1]=0;

  return t3;
}

template<Integer N,Integer M>
constexpr auto cumulative_array_and_zero
(const Array<std::vector<Integer>, N>& t1,const Array<std::vector<Integer>, M>& t2)
{

  Array<std::vector<Integer>, N+M> t3;

  // std::cout<<"<<<<<<<cumulative_array AND ZERO >>>>>>>>>>>5"<<std::endl;
  // std::cout<<"N+M"<<N+M<<std::endl;
  Integer levels=t1[0].size();
  // std::cout<<"t1"<<std::endl;
  for(Integer ii=0;ii<N;ii++)
    {
      t3[ii].resize(levels);
      // std::cout<<ii<<std::endl;
      for(Integer s=0;s<t1[ii].size();s++)
         t3[ii][s]=t1[ii][s];
      // std::cout<<std::endl;
   }
  // std::cout<<"end"<<std::endl;

  for(Integer ii=N;ii<N+M;ii++)
    {
      t3[ii].resize(levels);
      // std::cout<<ii<<std::endl;
      for(Integer s=0;s<levels;s++)
       t3[ii][s]=0;
     }
// std::cout<<"end"<<std::endl;

 //  std::cout<<"t2"<<std::endl;
 //  for(Integer ii=0;ii<t2.size();ii++)
 //    {
 //      std::cout<<std::endl;
 //      for(Integer s=0;s<t2[ii].size();s++)
 //      std::cout<<t2[ii][s]<<" ";
 //      std::cout<<std::endl;
 //   }



 //  std::cout<<"qui11"<<std::endl;
 //  auto levels=t2[0].size();



 //  std::cout<<"qui22"<<std::endl;
 //  for(Integer ii=0;ii<N+M;ii++)
 //    {
 //      std::cout<<"ii="<<ii<<std::endl;
 //      t3[ii].resize(levels);
 //   }


 //  std::cout<<"qui22"<<std::endl;
 //  for(Integer ii=0;ii<N;ii++)
 //    {
 //      std::cout<<"N="<<N<<std::endl;
 //      std::cout<<"ii="<<ii<<std::endl;
 //      std::cout<<"levels="<<levels<<"   "<<t3[ii].size()<<"  "<<t1[ii].size()<<std::endl;
 //      for(Integer s=0;s<levels;s++)
 //       t3[ii][s]=t1[ii][s];
 //   }
 //  std::cout<<"qui33"<<std::endl;
 //  for(Integer ii=N-1;ii<N+M+1;ii++)
 //    {
 //      for(Integer s=0;s<levels;s++)
 //       t3[ii+1][s]=0;
 //   }
 // std::cout<<"qui44"<<std::endl;
  return t3;
}

template<typename NewType, typename Arr>
class ArrayChangeTypeHelper;


template<typename NewType, typename T,Integer Dim>
class ArrayChangeTypeHelper<NewType, Array<T,Dim>>
{
public:
  using type=Array<NewType,Dim>;
};


template<typename NewType, typename Arr>
using ArrayChangeType=typename ArrayChangeTypeHelper<NewType, Arr>::type;




  template<Integer Ndofs>
  constexpr auto zero_array()
  {
   Array<Real,Ndofs> arr;
   for(Integer i=0;i<Ndofs;i++)
    arr[i]=0;
  return arr;
  }

  template<typename T,Integer Rows,Integer Cols, Integer NonZeroRows>
  class Matrix;

  template<typename T, Integer Dim1,Integer Dim2>
  constexpr auto ArrayOfArray2Matrix(const Array<Array<T,Dim2>,Dim1>& arr)
  {
    Matrix<T,Dim1,Dim2,-1> mat;
    for(Integer i=0;i<Dim1;i++)
      for(Integer j=0;j<Dim2;j++)
        mat(i,j)=arr[i][j];

    return mat;
  }




}

#endif //MARS_Array_HPP
