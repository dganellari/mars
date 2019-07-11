#ifndef MARS_NUMBER_HPP
#define MARS_NUMBER_HPP

#include "mars_base.hpp"

namespace mars{

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// Number<N>
//////// Class used to use Numbers as types and avoid constexpr arrays
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<Integer N>
class Number
{
public:
  static constexpr Integer value=N;
};
using Zero=Number<0>;
using One=Number<1>;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// TupleOfNumbers<Nmin,Nmax>
//////// Create a tuple of the Numbers between Nmin and Nmax
//////// Example:
//////// TupleOfNumbers<4,6>=tuple<Number<4>, Number<5>, Number<6>>
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<Integer Nmin,Integer Nmax,Integer N>
class TupleNumberCreateHelper;

template<Integer Nmin,Integer Nmax>
class TupleNumberCreateHelper<Nmin,Nmax,Nmax>
{
public:
  using type=std::tuple<Number<Nmax>>;
};
template<Integer Nmin,Integer Nmax,Integer N>
class TupleNumberCreateHelper
{
 public:
  static_assert(Nmin<=Nmax,"  In TupleNumberCreate Nmin<=Nmax");
  using single_type =std::tuple<Number<N>>;
  using type=decltype(std::tuple_cat(std::declval<single_type>(),
                      std::declval<typename TupleNumberCreateHelper<Nmin,Nmax,N+1>::type>()));
};
template<Integer Nmin,Integer Nmax>
using TupleOfNumbers=typename TupleNumberCreateHelper<Nmin,Nmax,Nmin>::type;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// TupleOfNumbers<Nmin,Nshigt,Nmax>
//////// Create a tuple of the Numbers between Nmin and Nmax, shifted by Nshift
//////// Example:
//////// TupleOfShiftedNumbers<2,3,10>=tuple<Number<2>, Number<5>, Number<8>>
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<Integer Nmin,Integer Nshift,Integer Nmax,Integer N>
class TupleShiftedNumberCreateHelper;

template<Integer Nmin,Integer Nshift,Integer Nmax>
class TupleShiftedNumberCreateHelper<Nmin,Nshift,Nmax,Nmax>
{
public:
  using type=std::tuple<Number<Nmax>>;
};
template<Integer Nmin,Integer Nshift,Integer Nmax,Integer N>
class TupleShiftedNumberCreateHelper
{
 public:
  static_assert(Nmin<=Nmax,"  In TupleNumberCreate Nmin<=Nmax");
  using single_type =std::tuple<Number<N>>;
  using type=decltype(std::tuple_cat(std::declval<single_type>(),
                      std::declval<typename TupleShiftedNumberCreateHelper<Nmin,Nshift,Nmax,(N+Nshift)>::type>()));
};

template<Integer Nmin,Integer Nshift,Integer Nmax>
using TupleOfShiftedNumbers=typename TupleShiftedNumberCreateHelper<Nmin,Nshift,((Nmax-Nmin)/Nshift)*Nshift+Nmin,Nmin>::type;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// FindNonZeroNumbers<TupleOfNumbers>
//////// Create a tuple of the positions of the TupleOfNumbers which are non zero
//////// Example:
//////// TupleOfNumbers=tuple<Number<0>, Number<4>, Number<6>,Number<0>>
//////// FindNonZeroNumbers=tuple<Number<1>, Number<2>>
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename TupleOfNumbers,Integer Nmax,Integer N>
class FindNonZeroNumbersHelper;

template <typename TupleOfNumbers,Integer Nmax>
class FindNonZeroNumbersHelper<TupleOfNumbers,Nmax,Nmax>
{
 public:
  using number=GetType<Nmax,TupleOfNumbers>;
  using type =typename std::conditional<std::is_same<number,Number<0>>::value,std::tuple<>,std::tuple<Number<Nmax>> >::type;
};
template <typename TupleOfNumbers,Integer Nmax,Integer N>
class FindNonZeroNumbersHelper
{
 public:
  using number=GetType<N,TupleOfNumbers>;
  using single_type =typename std::conditional<std::is_same<number,Number<0>>::value,std::tuple<>,std::tuple<Number<N>> >::type;
  using type=decltype(std::tuple_cat(std::declval<single_type>(),
                      std::declval<typename FindNonZeroNumbersHelper<TupleOfNumbers,Nmax,N+1>::type>()));
};

template <typename TupleOfNumbers>
using FindNonZeroNumbers=typename FindNonZeroNumbersHelper<TupleOfNumbers,TupleTypeSize<TupleOfNumbers>::value-1,0>::type;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// NumbersToArray<std::tuple<Number<Ns>...>>
//////// Convert tuple of numbers in the corresponding array of numbers
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T  >
class NumbersToArray;

template<Integer...Ns  >
class NumbersToArray<std::tuple<Number<Ns>...>>
{
public:
  static constexpr Integer value[]={Ns...};
};
template<Integer...Ns  >
constexpr Integer NumbersToArray<std::tuple<Number<Ns>...>>::value[];

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// TupleAllToUniqueMap<All,Unique>
//////// Unique is a subtuple of All, with no repetition
//////// TupleAllToUniqueMap creates a tuple of Numbers, with the same length as All.
//////// Each Number<N> defines the corresponding position of GetType<N,All> in Unique
//////// Example:
//////// All=std::tuple< A, B, C, D, A, B, E> 
//////// Unique=std::tuple<A,B,C,D,E>
//////// TupleAllToUniqueMap=std::tuple<Number<0>,Number<1>,Number<2>,Number<3>,Number<0>,Number<1>,Number<4>  > 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename All,typename Unique, Integer Nmax,Integer N>
class TupleAllToUniqueMapHelper;

template<typename All,typename Unique, Integer Nmax>
class TupleAllToUniqueMapHelper<All,Unique,Nmax,Nmax>
{
  public:
  using type=std::tuple<Number<ElementPosition<Nmax,All,Unique>::value>>;
};
template<typename All,typename Unique, Integer Nmax,Integer N=0>
class TupleAllToUniqueMapHelper
{
  public:
  using single_type=std::tuple<Number<ElementPosition<N,All,Unique>::value>>;
  using type=decltype(std::tuple_cat(std::declval<single_type>(),std::declval<typename TupleAllToUniqueMapHelper<All,Unique,Nmax,N+1>::type>()) );
};

template<typename All,typename Unique>
using TupleAllToUniqueMap=typename TupleAllToUniqueMapHelper<All,Unique,TupleTypeSize<All>::value-1,0>::type;
 



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// TupleOfTupleToBooleanTuple<Tuple>
//////// Transform a Tuple of tuples into a tuple of Zero (if the component is an empty tuple) or One (otherwise)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename Tuple, Integer Nmax, Integer N>
class TupleOfTupleToBooleanTupleHelper;

template <typename Tuple, Integer Nmax>
class TupleOfTupleToBooleanTupleHelper<Tuple,Nmax,Nmax>
{
  public:  
    using type=std::tuple<typename std::conditional< std::is_same< GetType<Nmax,Tuple>, std::tuple<> >::value, Number<0>,Number<1>  >::type>;
  };

template <typename Tuple, Integer Nmax, Integer N>
class TupleOfTupleToBooleanTupleHelper
{
  public:
    using single_type=std::tuple<typename std::conditional< std::is_same< GetType<N,Tuple>, std::tuple<> >::value, Number<0>,Number<1>  >::type>;
    using type=decltype(std::tuple_cat(std::declval<single_type>(),
                                       std::declval<typename TupleOfTupleToBooleanTupleHelper<Tuple,Nmax,N+1>::type>()));
};

template <typename Tuple>
using TupleOfTupleToBooleanTuple=typename TupleOfTupleToBooleanTupleHelper<Tuple, TupleTypeSize<Tuple>::value-1,0>::type;




















//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// StaticScalarProductForVector<NonZeroNumbers>(v1,v2)
//////// It takes the vectors v1, v2 and the NonZeroNumbers tuples, which contains the position where both components of v1i and v2i are non-zero
//////// Then it makes the sum of all this products.
//////// Example:
//////// A=[1 0 3 0 0 6], B=[1 1 1 1 1 0]; NonZeroNumbers=tuple<Number<0>,Number<2>>
//////// StaticScalarProductForVector=1 * 1 + 3 * 1;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename NonZeroNumbers, Integer Nmax,Integer N,typename T, Integer Dim>
constexpr typename std::enable_if< (std::is_same<NonZeroNumbers,std::tuple<>>::value),T>::type
StaticScalarProductForVectorHelper(const Vector<T,Dim> & v1,const  Vector<T,Dim>& v2 )
{
  return 0;
}
template<typename NonZeroNumbers, Integer Nmax,Integer N,typename T, Integer Dim>
constexpr typename std::enable_if< (N==Nmax) && !(std::is_same<NonZeroNumbers,std::tuple<>>::value),T>::type
StaticScalarProductForVectorHelper(const Vector<T,Dim> & v1,const  Vector<T,Dim>& v2 )
{
  return v1[GetType<N,NonZeroNumbers>::value]*v2[GetType<N,NonZeroNumbers>::value];
}
template<typename NonZeroNumbers, Integer Nmax,Integer N,typename T, Integer Dim>
constexpr typename std::enable_if< (N<Nmax)&& !(std::is_same<NonZeroNumbers,std::tuple<>>::value),T>::type
StaticScalarProductForVectorHelper(const Vector<T,Dim> & v1,const  Vector<T,Dim>& v2 )
{
  return v1[GetType<N,NonZeroNumbers>::value]*v2[GetType<N,NonZeroNumbers>::value]+ 
         StaticScalarProductForVectorHelper<NonZeroNumbers,Nmax,N+1>(v1,v2);
}
template<typename NonZeroNumbers,typename T, Integer Dim>
constexpr T StaticScalarProductForVector(const Vector<T,Dim> & v1,const  Vector<T,Dim>& v2 )
{
  return StaticScalarProductForVectorHelper<NonZeroNumbers,TupleTypeSize<NonZeroNumbers>::value-1,0,T,Dim>(v1,v2);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// StaticVectorContractionFindNonZeros<Tuple1,Tuple2>
//////// It takes Tuple1, Tuple2, i.e. the tuples of Zero,One related to two vectors A and B
//////// It checks which product of A*B=A0*B0+A1*B1 + A2*B2 +...  is zero (for example A1*B1=0, because A1=0 or B1=0 or both are zero ) 
//////// It builds a tuple whose elements are pairs of the surviving C1i, C2i 
//////// Example:
//////// A=[1 0 3 0 0 6], B=[1 1 1 1 1 0];
//////// then the result is: std::tuple<Number<0>, Number<2>  > 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename Tuple1, typename Tuple2, Integer Nmax, Integer N>
class StaticVectorContractionFindNonZerosHelper;

template <typename Tuple1, typename Tuple2, Integer Nmax>
class StaticVectorContractionFindNonZerosHelper<Tuple1,Tuple2,Nmax,Nmax>
{
  public:  
    static_assert(TupleTypeSize<Tuple1>::value==TupleTypeSize<Tuple2>::value,"In StaticVectorContractionFindNonZerosHelper tuple1 and tuple2 must have same size ");
    using type1=GetType<Nmax,Tuple1>;
    using type2=GetType<Nmax,Tuple2>;
    using type=typename std::conditional<Greater(type1::value*type2::value,0), std::tuple<Number<Nmax>>, std::tuple<> >::type;
};
template <typename Tuple1, typename Tuple2, Integer Nmax, Integer N>
class StaticVectorContractionFindNonZerosHelper
{
  public:
    static_assert(TupleTypeSize<Tuple1>::value==TupleTypeSize<Tuple2>::value,"In StaticVectorContractionFindNonZerosHelper tuple1 and tuple2 must have same size ");
    using type1=GetType<N,Tuple1>;
    using type2=GetType<N,Tuple2>;
    using single_type=typename std::conditional<Greater(type1::value*type2::value,0), std::tuple<Number<N>>, std::tuple<> >::type;
    using type=decltype(std::tuple_cat(std::declval<single_type>(),
                                       std::declval<typename StaticVectorContractionFindNonZerosHelper<Tuple1,Tuple2,Nmax,N+1>::type>()));
};
template <typename Tuple1, typename Tuple2>
using StaticVectorContractionFindNonZeros=typename StaticVectorContractionFindNonZerosHelper<Tuple1,Tuple2,TupleTypeSize<Tuple1>::value-1,0>::type;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// StaticMatrixContractionFindNonZeros<Tuple1,C1,Tuple2,C2>
//////// It takes Tuple1, Tuple2, i.e. the tuples of Zero,One related to two matrices A and B
//////// It also takes the indices of the row related to A and the column related to B, i.e. C1 and C2
//////// It checks which product of C1*C2=C10*C20+C11*C21 + C12*C22 +...  is zero (for example C11*C21=0, because C11=0 or C21=0 or both ) 
//////// It builds a tuple whose elements are pairs of the surviving C1i, C2i 
//////// Example:
//////// A=[1 2 3;4 5 6; 7 8 9], B=[1 0 0; 0 1 0; 0 0 1];
//////// Tuple1=[1,...,1],  Tuple2=[1 0 0; 0 1 0; 0 0 1];
//////// C1=[3 4 5] (second row), C2=[2 5 8] (last column)
//////// then the result is: std::tuple< std::tuple<Number<4>, Number<5> >  > 
//////// because in C1*C2=4 * 0 + 5 * 1 + 6 * 0 only the second product survives 
//////// and is related to the 4-th element in tuple1 and to the 5-th in tuple2
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename Tuple1, typename C1, typename Tuple2, typename C2, Integer Nmax, Integer N>
class StaticVectorContractionFindNonZerosHelper2;

template <typename Tuple1, typename C1, typename Tuple2, typename C2, Integer Nmax>
class StaticVectorContractionFindNonZerosHelper2<Tuple1,C1,Tuple2,C2,Nmax,Nmax>
{
  public:  
    static_assert(TupleTypeSize<C1>::value==TupleTypeSize<C2>::value,"In StaticVectorContractionFindNonZerosHelper2 row and col must have same size ");
    static constexpr Integer N1=GetType<Nmax,C1>::value;
    static constexpr Integer N2=GetType<Nmax,C2>::value;
    using type1=GetType< N1 ,Tuple1>;
    using type2=GetType< N2 ,Tuple2>;
    using type=typename std::conditional<Greater(type1::value*type2::value,0), std::tuple<std::tuple<Number<N1> , Number<N2>> >, std::tuple<> >::type;
};

template <typename Tuple1, typename C1, typename Tuple2, typename C2, Integer Nmax, Integer N>
class StaticVectorContractionFindNonZerosHelper2
{
  public:
    static_assert(TupleTypeSize<C1>::value==TupleTypeSize<C2>::value,"In StaticVectorContractionFindNonZerosHelper2 row and col must have same size ");
    static constexpr Integer N1=GetType<N,C1>::value;
    static constexpr Integer N2=GetType<N,C2>::value;
    using type1=GetType< N1 ,Tuple1>;
    using type2=GetType< N2 ,Tuple2>;
    using single_type=typename std::conditional<Greater(type1::value*type2::value,0), std::tuple<std::tuple<Number<N1> , Number<N2>> >, std::tuple<> >::type;
    using type=decltype(std::tuple_cat(std::declval<single_type>(),
                                       std::declval<typename StaticVectorContractionFindNonZerosHelper2<Tuple1,C1,Tuple2,C2,Nmax,N+1>::type>()));
};

template <typename Tuple1, typename C1, typename Tuple2, typename C2>
using StaticMatrixContractionFindNonZeros=typename StaticVectorContractionFindNonZerosHelper2<Tuple1,C1,Tuple2,C2,TupleTypeSize<C1>::value-1,0>::type;



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// StaticScalarProductForMatrixProduct<NonZeroNumbers2>(m1,m2)
//////// It returns the product of the components of the matrices m1 and m2 identified by NonZeroNumbers2
//////// Example:
//////// A=[1 2 3;4 5 6; 7 8 9], B=[1 0 1; 0 1 0; 0 0 1];
//////// NonZeroNumbers2=tuple< tuple<Number<3>,Number<2>> , tuple<Number<5>,Number<8>> > (second row times third column) 
//////// StaticScalarProductForMatrixProduct= 4 * 1 + 6 * 1
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename NonZeroNumbers2, Integer Nmax,Integer N,typename T, Integer Rows,Integer CommonDim,Integer Cols>
constexpr typename std::enable_if< (std::is_same<NonZeroNumbers2,std::tuple<>>::value),T>::type
StaticScalarProductForMatrixProductHelper(const Matrix<T,Rows,CommonDim> & m1,const  Matrix<T,CommonDim,Cols>& m2 )
{
  return 0;
}
template<typename NonZeroNumbers2, Integer Nmax,Integer N,typename T, Integer Rows,Integer CommonDim,Integer Cols>
constexpr typename std::enable_if< (N==Nmax) && !(std::is_same<NonZeroNumbers2,std::tuple<>>::value),T>::type
StaticScalarProductForMatrixProductHelper(const Matrix<T,Rows,CommonDim> & m1,const  Matrix<T,CommonDim,Cols>& m2 )
{
  return m1(GetType<0,GetType<N,NonZeroNumbers2>>::value)*m2(GetType<1,GetType<N,NonZeroNumbers2>>::value);
}
template<typename NonZeroNumbers2, Integer Nmax,Integer N,typename T, Integer Rows,Integer CommonDim,Integer Cols>
constexpr typename std::enable_if< (N<Nmax)&& !(std::is_same<NonZeroNumbers2,std::tuple<>>::value),T>::type
StaticScalarProductForMatrixProductHelper(const Matrix<T,Rows,CommonDim> & m1,const  Matrix<T,CommonDim,Cols>& m2 )
{
  return m1(GetType<0,GetType<N,NonZeroNumbers2>>::value)*m2(GetType<1,GetType<N,NonZeroNumbers2>>::value)+ 
         StaticScalarProductForMatrixProductHelper<NonZeroNumbers2,Nmax,N+1>(m1,m2);
}
template<typename NonZeroNumbers2,typename T, Integer Rows,Integer CommonDim,Integer Cols>
constexpr T StaticScalarProductForMatrixProduct(const Matrix<T,Rows,CommonDim> & m1,const  Matrix<T,CommonDim,Cols>& m2 )
{
  return StaticScalarProductForMatrixProductHelper<NonZeroNumbers2,TupleTypeSize<NonZeroNumbers2>::value-1,0,T,Rows,CommonDim,Cols>(m1,m2);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// StaticBooleanMatrixMatrixMultiplicationFindNonZeros<Rows1,Cols1,Tuple1,Rows2,Cols2,Tuple2>
//////// Tuple1/2 contains Zero and One for the corresponding matrix Mat1/2 of dimension Rows1/2 x Cols1/2
//////// For brevity, we call Mat the matrix result of the product of Mat1*Mat2
//////// StaticBooleanMatrixMatrixMultiplicationFindNonZeros returns a tuple of subtuples.
//////// A subtuple is empty if, in that position, Mat is known at compile-time to be zero
//////// Otherwise, a subtuple at position (i,j) contains the pairs of the components of Mat1 and Mat2 
////////                  used in the dot product for computing Mat(i,j) which are both non-zero
//////// Example:
//////// A=[1 2 3;4 5 6; 0 0 1], B=[1 0 1; 0 1 1; 0 0 1];
//////// NonZeroProducts=tuple< tuple< tuple< Number<0>,Number<0> > > ,
////////                        tuple< tuple< Number<1>,Number<4> > > 
////////                        tuple< tuple< Number<0>,Number<2> > , tuple< Number<1>,Number<5> >, tuple< Number<2>,Number<8> > > 
////////                        tuple< tuple< Number<3>,Number<0> > >, 
////////                        tuple< tuple< Number<4>,Number<4> > > ,
////////                        tuple< tuple< Number<3>,Number<2> >, tuple< Number<4>,Number<5> >, tuple< Number<5>,Number<8> > > ,
////////                        tuple<>,
////////                        tuple<>,
////////                        tuple< tuple<Number<8>,Number<8> >> >
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <Integer Rows1,Integer Cols1,typename Tuple1, Integer Rows2, Integer Cols2, typename Tuple2,Integer N1,Integer N2 >
class StaticBooleanMatrixMatrixMultiplicationFindNonZerosHelper;

template <Integer Rows1,Integer Cols1,typename Tuple1, Integer Rows2, Integer Cols2, typename Tuple2>
class StaticBooleanMatrixMatrixMultiplicationFindNonZerosHelper<Rows1,Cols1,Tuple1,Rows2,Cols2,Tuple2,Rows1-1,Cols2-1>
{
  public:  
    static_assert(TupleTypeSize<Tuple1>::value==Rows1*Cols1,"In StaticBooleanMatrixMatrixMultiplicationFindNonZerosHelper tuple1 must have length=Rows1*Cols1");
    static_assert(TupleTypeSize<Tuple2>::value==Rows2*Cols2,"In StaticBooleanMatrixMatrixMultiplicationFindNonZerosHelper tuple2 must have length=Rows2*Cols2");
    static_assert(Cols1==Rows2,"In StaticBooleanMatrixMatrixMultiplicationFindNonZerosHelper Cols1==Rows2");

    using R1=TupleOfNumbers<Cols1*(Rows1-1),Cols1*Rows1-1>;
    using C1=TupleOfShiftedNumbers<Cols2-1,Cols2, Rows2*Cols2-1>;
    using type=std::tuple< StaticMatrixContractionFindNonZeros<Tuple1,R1,Tuple2,C1> >;
};


template <Integer Rows1,Integer Cols1,typename Tuple1, Integer Rows2, Integer Cols2, typename Tuple2, Integer N1>
class StaticBooleanMatrixMatrixMultiplicationFindNonZerosHelper<Rows1,Cols1,Tuple1,Rows2,Cols2,Tuple2,N1,Cols2-1>
{
  public:
    static_assert(TupleTypeSize<Tuple1>::value==Rows1*Cols1,"In StaticBooleanMatrixMatrixMultiplicationFindNonZerosHelper tuple1 must have length=Rows1*Cols1");
    static_assert(TupleTypeSize<Tuple2>::value==Rows2*Cols2,"In StaticBooleanMatrixMatrixMultiplicationFindNonZerosHelper tuple2 must have length=Rows2*Cols2");
    static_assert(Cols1==Rows2,"In StaticBooleanMatrixMatrixMultiplicationFindNonZerosHelper Cols1==Rows2");

    using R1=TupleOfNumbers<Cols1*N1,Cols1*(N1+1)-1>;
    using C1=TupleOfShiftedNumbers<Cols2-1,Cols2, Rows2*Cols2-1>;
    using single_type=std::tuple< StaticMatrixContractionFindNonZeros<Tuple1,R1,Tuple2,C1> >;
    using type=decltype(std::tuple_cat(std::declval<single_type>(),
                                       std::declval<typename StaticBooleanMatrixMatrixMultiplicationFindNonZerosHelper<Rows1,Cols1,Tuple1,Rows2,Cols2,Tuple2,N1+1,0>::type>()));
};



template <Integer Rows1,Integer Cols1,typename Tuple1, Integer Rows2, Integer Cols2, typename Tuple2, Integer N1, Integer N2>
class StaticBooleanMatrixMatrixMultiplicationFindNonZerosHelper
{
  public:
    static_assert(TupleTypeSize<Tuple1>::value==Rows1*Cols1,"In StaticBooleanMatrixMatrixMultiplicationFindNonZerosHelper tuple1 must have length=Rows1*Cols1");
    static_assert(TupleTypeSize<Tuple2>::value==Rows2*Cols2,"In StaticBooleanMatrixMatrixMultiplicationFindNonZerosHelper tuple2 must have length=Rows2*Cols2");
    static_assert(Cols1==Rows2,"In StaticBooleanMatrixMatrixMultiplicationFindNonZerosHelper Cols1==Rows2");

    using R1=TupleOfNumbers<Cols1*N1,Cols1*(N1+1)-1>;
    using C1=TupleOfShiftedNumbers<N2,Cols2, Rows2*Cols2-1>;
    using single_type=std::tuple< StaticMatrixContractionFindNonZeros<Tuple1,R1,Tuple2,C1> >;
    using type=decltype(std::tuple_cat(std::declval<single_type>(),
                                       std::declval<typename StaticBooleanMatrixMatrixMultiplicationFindNonZerosHelper<Rows1,Cols1,Tuple1,Rows2,Cols2,Tuple2,N1,N2+1>::type>()));
};


template <Integer Rows1,Integer Cols1,typename Tuple1, Integer Rows2, Integer Cols2, typename Tuple2>
using StaticBooleanMatrixMatrixMultiplicationFindNonZeros=typename StaticBooleanMatrixMatrixMultiplicationFindNonZerosHelper<Rows1,Cols1,Tuple1,Rows2,Cols2,Tuple2,0,0>::type;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// StaticMatrixProduct<NonZeroNumbers, NonZeroProducts>(m1,m2,m3)
//////// It computes the product m3=m1*m2, omitting all the computations where zeros appear
//////// NonZeroNumbers: positions in m3 which are non-zero
//////// NonZeroProducts: for each non-zero component of m3, contains the tuple of the pairs of the components of m1 and m2 
////////                  used in the dot product for computing m3 which are both non-zero 
////////                  (see StaticBooleanMatrixMatrixMultiplicationFindNonZeros)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename NonZeroNumbers, typename NonZeroProducts,Integer N,typename T, Integer Rows,Integer Cols,Integer CommonDim>
inline typename std::enable_if<(N==TupleTypeSize<NonZeroNumbers>::value-1)  , void>::type
StaticMatrixProductAux(const Matrix<T,Rows,CommonDim> & m1,const  Matrix<T,CommonDim,Cols>& m2, Matrix<T,Rows,Cols> & m3)
{
 constexpr Integer Tot=GetType<N,NonZeroNumbers>::value;

 // constexpr Integer J=Modulo<Tot,Rows>::value;
 // constexpr Integer I=Tot/Rows;
 // m3(I,J)=0;
 //  for(Integer kk=0;kk<CommonDim;kk++)
    // m3(I,J)+=m1(I,kk)*m2(kk,J);

m3(Tot)=StaticScalarProductForMatrixProduct<GetType<Tot,NonZeroProducts>>(m1,m2);

}

template<typename NonZeroNumbers, typename NonZeroProducts,Integer N,typename T, Integer Rows,Integer Cols,Integer CommonDim>
inline typename std::enable_if<(N<TupleTypeSize<NonZeroNumbers>::value-1)  , void>::type
StaticMatrixProductAux(const Matrix<T,Rows,CommonDim> & m1,const  Matrix<T,CommonDim,Cols>& m2, Matrix<T,Rows,Cols> & m3)
{
   constexpr Integer Tot=GetType<N,NonZeroNumbers>::value;

 // constexpr Integer J=Modulo<Tot,Rows>::value;
 // constexpr Integer I=Tot/Rows;
  // m3(I,J)=0;
  // for(Integer kk=0;kk<CommonDim;kk++)
  //   m3(I,J)+=m1(I,kk)*m2(kk,J);

   m3(Tot)=StaticScalarProductForMatrixProduct<GetType<Tot,NonZeroProducts>>(m1,m2);


  StaticMatrixProductAux<NonZeroNumbers,NonZeroProducts,N+1>(m1,m2,m3);
}



template<typename NonZeroNumbers, typename NonZeroProducts,typename T, Integer Rows,Integer Cols,Integer CommonDim>
inline typename std::enable_if<!std::is_same<NonZeroNumbers,std::tuple<>>::value  , void>::type
StaticMatrixProduct(const Matrix<T,Rows,CommonDim> & m1,const  Matrix<T,CommonDim,Cols>& m2, Matrix<T,Rows,Cols> & m3)
{
 
 StaticMatrixProductAux<NonZeroNumbers,NonZeroProducts,0>(m1,m2,m3);
}


template<typename NonZeroNumbers, typename NonZeroProducts,typename T, Integer Rows,Integer Cols,Integer CommonDim>
inline typename std::enable_if<std::is_same<NonZeroNumbers,std::tuple<>>::value  , void>::type
 StaticMatrixProduct(const Matrix<T,Rows,CommonDim> & m1,const  Matrix<T,CommonDim,Cols>& m2, Matrix<T,Rows,Cols> & m3)
{}










// template<Integer M,Integer N>
// class Multiply2< Number<M>, Number<N> >
// {
// public:
//   static constexpr Integer value=M*N;
//   using type=Number<value>;
// };

// template<Integer M,Integer N>
// class Addition2< Number<M>, Number<N> >
// {
// public:
//   static constexpr Integer value=M+N;
//   using type=Number<value>;
// };





// template<typename Vector1, typename Vector2, Integer Nmax, Integer N>
// class StaticBooleanContractionHelper;

// template<typename Vector1, typename Vector2, Integer Nmax>
// class StaticBooleanContractionHelper<Vector1,Vector2,Nmax,Nmax>
// {
//  public:
//  static_assert(TupleTypeSize<Vector1>::value==TupleTypeSize<Vector2>::value," Static boolean contraction requires the two vectors to have the same length ");
//  using type = typename Multiply2<GetType<Nmax,Vector1>,GetType<Nmax,Vector2>>::type;
// };

// template<typename Vector1, typename Vector2, Integer Nmax,Integer N>
// class StaticBooleanContractionHelper
// {
//  public:
//  static_assert(TupleTypeSize<Vector1>::value==TupleTypeSize<Vector2>::value," Static boolean contraction requires the two vectors to have the same length ");
//  using single_type = typename Multiply2<GetType<N,Vector1>,GetType<N,Vector2>>::type;
//  using type= typename Addition2<single_type, typename StaticBooleanContractionHelper<Vector1,Vector2,Nmax,N+1>::type >::type;

// };


// template<typename Vector1, typename Vector2>
// using StaticBooleanContraction=typename StaticBooleanContractionHelper<Vector1,Vector2,TupleTypeSize<Vector1>::value-1,0>::type;




// template <Integer Rows,Integer Cols,typename Tuple1, Integer Dim, typename Tuple2, Integer Nmax, Integer N>
// class StaticBooleanMatrixVectorMultiplicationFindNonZerosHelper;

// template <Integer Rows,Integer Cols,typename Tuple1, Integer Dim, typename Tuple2, Integer Nmax>
// class StaticBooleanMatrixVectorMultiplicationFindNonZerosHelper<Rows,Cols,Tuple1,Dim,Tuple2,Nmax,Nmax>
// {
//   public:  
//     static_assert(TupleTypeSize<Tuple1>::value==Rows*Cols,"In StaticBooleanMatrixVectorMultiplicationFindNonZerosHelper tuple1 must have length=Rows*Cols");
//     static_assert(TupleTypeSize<Tuple2>::value==Dim,"In StaticBooleanMatrixVectorMultiplicationFindNonZerosHelper tuple2 must have length=Dim");
//     static_assert(Dim==Cols,"In StaticBooleanMatrixVectorMultiplicationFindNonZerosHelper Dim must be equal to Cols");

//     using type1=SubTupleType<Cols*Nmax,Cols*(Nmax+1)-1,Tuple1>;
//     using type2=Tuple2;
//     using type=std::tuple< StaticVectorContractionFindNonZeros<type1,type2> >;
// };


// template <Integer Rows,Integer Cols,typename Tuple1, Integer Dim, typename Tuple2, Integer Nmax, Integer N>
// class StaticBooleanMatrixVectorMultiplicationFindNonZerosHelper
// {
//   public:
//     static_assert(TupleTypeSize<Tuple1>::value==Rows*Cols,"In StaticBooleanMatrixVectorMultiplicationFindNonZerosHelper tuple1 must have length=Rows*Cols");
//     static_assert(TupleTypeSize<Tuple2>::value==Dim,"In StaticBooleanMatrixVectorMultiplicationFindNonZerosHelper tuple2 must have length=Dim");
//     static_assert(Dim==Cols,"In StaticBooleanMatrixVectorMultiplicationFindNonZerosHelper Dim must be equal to Cols");

//     using type1=SubTupleType<Cols*N,Cols*(N+1)-1,Tuple1>;
//     using type2=Tuple2;
//     using single_type=std::tuple< StaticVectorContractionFindNonZeros<type1,type2> >;

//     using type=decltype(std::tuple_cat(std::declval<single_type>(),
//                                        std::declval<typename StaticBooleanMatrixVectorMultiplicationFindNonZerosHelper<Rows,Cols,Tuple1,Dim,Tuple2,Nmax,N+1>::type>()));
// };



// template <Integer Rows,Integer Cols,typename Tuple1, Integer Dim, typename Tuple2>
// using StaticBooleanMatrixVectorMultiplicationFindNonZeros=typename StaticBooleanMatrixVectorMultiplicationFindNonZerosHelper<Rows,Cols,Tuple1,Dim,Tuple2,Rows-1,0>::type;







// template<Integer Rows_,Integer Cols_,typename...Args>
// class StaticBooleanMatrix
// {
//  public:
//   static constexpr Integer Rows=Rows_;
//   static constexpr Integer Cols=Cols_;
//   static_assert(sizeof...(Args)==Rows*Cols,"Wrong boolean input wrt the dimension of the matrix");
//   using type=std::tuple<Args...>;
// };

// template<Integer Dim_,typename...Args>//,typename...Args>
// class StaticBooleanVector
// {
//  public:
//   static constexpr Integer Dim=Dim_;
//   using type=std::tuple<Args...>;
//   static_assert(sizeof...(Args)==Dim,"Wrong boolean input wrt the dimension of the vector");

// };


// template<Integer Rows,Integer Cols,typename Tuple1, Integer Dim,typename Tuple2, Integer N>
// class StaticBooleanMatrixVectorMultiplicationHelper;


// template<Integer Rows,Integer Cols,typename Tuple1, Integer Dim,typename Tuple2>
// class StaticBooleanMatrixVectorMultiplicationHelper<Rows,Cols,Tuple1,Dim,Tuple2,Rows-1>
// {
// public:
// static constexpr Integer N=Rows-1;
// using Row=SubTupleType<N*Cols,(N+1)*Cols-1,Tuple1 >;
// using type=std::tuple< StaticBooleanContraction<Row,Tuple2>   >;
// };


// template<Integer Rows,Integer Cols,typename Tuple1, Integer Dim,typename Tuple2, Integer N=0>
// class StaticBooleanMatrixVectorMultiplicationHelper
// {
//  public:
//  using Row=SubTupleType<N*Cols,(N+1)*Cols-1,Tuple1 >;
//  using single_type=std::tuple< StaticBooleanContraction<Row,Tuple2>   >;
//  using type=decltype(std::tuple_cat(std::declval<single_type>(),
//                            std::declval<typename StaticBooleanMatrixVectorMultiplicationHelper<Rows,Cols,Tuple1,Dim,Tuple2,N+1>::type>()));
// };

// template<Integer Rows,Integer Cols,typename Tuple1, Integer Dim,typename Tuple2>
// using StaticBooleanMatrixVectorMultiplication=typename StaticBooleanMatrixVectorMultiplicationHelper<Rows,Cols,Tuple1,Dim,Tuple2,0>::type;


// template<Integer Rows,Integer Cols,typename...Args1, Integer Dim,typename...Args2>
// class Multiply2< StaticBooleanMatrix<Rows,Cols,Args1...> , StaticBooleanVector<Dim, Args2... > >
// {
// public:
//   static_assert(Cols==Dim,"static matrix vector multiplication requires Cols==Dim");
//   // using type=typename StaticBooleanVector<Rows,>::type ;

// };













}

#endif