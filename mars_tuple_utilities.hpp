


#ifndef MARS_TUPLE_UTILITIES_HPP
#define MARS_TUPLE_UTILITIES_HPP

#include "mars_base.hpp"
#include "mars_static_math.hpp"

namespace mars{




template<typename T,Integer Nmax,Integer N>
class TupleOfTypeTCreateHelper;

template<typename T,Integer Nmax>
class TupleOfTypeTCreateHelper<T,Nmax,Nmax>
{
public:
  using type=std::tuple<T>;
};

template<typename T,Integer Nmax,Integer N>
class TupleOfTypeTCreateHelper
{
 public:
  static_assert(N<=Nmax,"  In TupleNumberCreate Nmin<=Nmax");
  using single_type =std::tuple<T>;
  using type=decltype(std::tuple_cat(std::declval<single_type>(),
                      std::declval<typename TupleOfTypeTCreateHelper<T,Nmax,N+1>::type>()));
};
template<typename T,Integer Size>
using TupleOfTypeTCreate=typename TupleOfTypeTCreateHelper<T,Size-1,0>::type;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// Get the N-th type of the tuple std::tuple<T, Ts...>
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <Integer N, typename... Ts>
class GetHelper;

template <Integer N, typename T, typename... Ts>
class GetHelper<N, std::tuple<T, Ts...>>
{   public: 
    using type = typename GetHelper<N - 1, std::tuple<Ts...>>::type;
};

template <typename T, typename... Ts>
class GetHelper<0, std::tuple<T, Ts...>>
{   public: 
    using type = T;
};

template <>
class GetHelper<0, std::tuple<>>
{   public: 
    using type = std::tuple<>;
};

template <Integer N, typename... Ts>
using GetTypeAux=typename GetHelper<N,Ts...>::type;




template <Integer N, typename T, typename... Ts>
class GetHelper<N, const std::tuple<T, Ts...>>
{   public: 
    using type = typename GetHelper<N - 1, std::tuple<Ts...>>::type;
};

template <typename T, typename... Ts>
class GetHelper<0, const std::tuple<T, Ts...>>
{   public: 
    using type = T;
};

template <>
class GetHelper<0, const std::tuple<>>
{   public: 
    using type = std::tuple<>;
};

template <Integer N, typename... Ts>
using GetTypeAux=typename GetHelper<N,Ts...>::type;


template <typename Tuple, Integer N, Integer...Ns >
class GetTypeHelper;

template <typename Tuple,Integer N>
class GetTypeHelper< Tuple,N>
{
public:
  using type = GetTypeAux<N,Tuple>;
};


template <typename Tuple,Integer N, Integer...Ns>
class GetTypeHelper
{
public:
  using type =typename GetTypeHelper<GetTypeAux<N,Tuple>,Ns...>::type;
};

template <typename Tuple,Integer N, Integer...Ns>
using GetType=typename GetTypeHelper<Tuple,N,Ns...>::type;


// template <Integer N,typename Tuple>
// constexpr auto get(const Tuple& tuple){return std::get<N>(tuple);}


// template <Integer...N,Integer M,typename Tuple>
// constexpr auto get(const Tuple& tuple)
// {return std::get<M>(std::get<N>(tuple));}






template <Integer N,Integer...Ns,typename Tuple>
constexpr std::enable_if_t<(0==sizeof...(Ns)),const GetType<Tuple,N>&> 
tuple_get(const Tuple& tuple)
{return std::get<N>(tuple);}


template <Integer N,Integer...Ns,typename Tuple>
constexpr std::enable_if_t<(0<sizeof...(Ns)),const GetType<Tuple,N,Ns...>&> 
tuple_get(const Tuple& tuple)
{return tuple_get<Ns...>(std::get<N>(tuple));}




template <Integer N,Integer...Ns,typename Tuple>
constexpr std::enable_if_t<(0==sizeof...(Ns)), GetType<Tuple,N>&> 
tuple_get(Tuple& tuple)
{return std::get<N>(tuple);}


template <Integer N,Integer...Ns,typename Tuple>
constexpr std::enable_if_t<(0<sizeof...(Ns)), GetType<Tuple,N,Ns...>&> 
tuple_get(Tuple& tuple)
{return tuple_get<Ns...>(std::get<N>(tuple));}


template<typename...Args>
std::tuple<Args...> add_costant (const std::tuple<Args...>& t1,const Integer& value);


template<std::size_t Nelem_dofs>
std::vector<std::array<Integer, Nelem_dofs>> add_costant
(const std::vector<std::array<Integer, Nelem_dofs>>& t1,const Integer& value)
{
  std::vector<std::array<Integer, Nelem_dofs>> t2;
  const auto& t1size=t1.size();
  t2.resize(t1size);

  std::cout<<"value=="<<value<<std::endl;
  for(Integer ii=0;ii<t1size;ii++)
    for(Integer jj=0;jj<t1[ii].size();jj++)
       t2[ii][jj]=t1[ii][jj]+value;
  return t2;
}

template<Integer N,typename...Args>
typename std::enable_if<sizeof...(Args)==N+1, void>::type
 add_costant_tuple_loop(const std::tuple<Args...>& t1, std::tuple<Args...>& t2,const Integer& value)
{
std::get<N>(t2)=add_costant(std::get<N>(t1),value);
}

template<Integer N,typename...Args>
typename std::enable_if< N+1<sizeof...(Args), void>::type
 add_costant_tuple_loop(const std::tuple<Args...>& t1, std::tuple<Args...>& t2,const Integer& value)
{
  std::get<N>(t2)=add_costant(std::get<N>(t1),value);
  add_costant_tuple_loop<N+1,Args...>(t1,t2,value);
}


template<typename...Args>
std::tuple<Args...> add_costant (const std::tuple<Args...>& t1,const Integer& value)
{
  std::tuple<Args...> t2=t1;

  add_costant_tuple_loop<0,Args...>(t1,t2,value);

  return t2;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////// Remove all duplicates types in the tuple: 
///////// std::tuple<int,int,char,int,double> -> std::tuple<int,char,double>
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <class Haystack, class Needle>
struct ContainsHelper;

template <class Car, class... Cdr, class Needle>
struct ContainsHelper<std::tuple<Car, Cdr...>, Needle> : ContainsHelper<std::tuple<Cdr...>, Needle>
{};

template <class... Cdr, class Needle>
struct ContainsHelper<std::tuple<Needle, Cdr...>, Needle> : std::true_type
{};

template <class Needle>
struct ContainsHelper<std::tuple<>, Needle> : std::false_type
{};



template <class Out, class In>
struct RemoveTupleDuplicatesHelper;

template <class... Out, class InCar, class... InCdr>
struct RemoveTupleDuplicatesHelper<std::tuple<Out...>, std::tuple<InCar, InCdr...>>
{
  using type = typename std::conditional<
    ContainsHelper<std::tuple<Out...>, InCar>::value
    , typename RemoveTupleDuplicatesHelper<std::tuple<Out...>, std::tuple<InCdr...>>::type
    , typename RemoveTupleDuplicatesHelper<std::tuple<Out..., InCar>, std::tuple<InCdr...>>::type
  >::type;
};

template <class Out>
struct RemoveTupleDuplicatesHelper<Out, std::tuple<>>
{
  using type = Out;
};


template <class T>
using RemoveTupleDuplicates = typename RemoveTupleDuplicatesHelper<std::tuple<>, T>::type;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// Cat tuple types: TupleCatType< std::tuple<Ts1...>, std::tuple<Ts2...>,....  > =TupleCatType< std::tuple<Ts1...,Ts2...> > 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename ... input_t>
using TupleCatType=
decltype(std::tuple_cat(std::declval<input_t>()...));

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// Remove all the occurences of the type T from std::tuple<Ts...> 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename...Ts>
class TupleRemoveTypeHelper;

template<typename T, typename...Ts>
class TupleRemoveTypeHelper<T,std::tuple<Ts...>> 
{
public:
using type= TupleCatType<
    typename std::conditional<
        std::is_same<T, Ts>::value,
        std::tuple<>,
        std::tuple<Ts>
    >::type...
>;
};

template<typename T, typename...Ts>
using TupleRemoveType=typename TupleRemoveTypeHelper<T,Ts...>::type;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// Remove the first occurence of the type Tremove from std::tuple<Ts...> 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename...Ts>
class TupleRemovesingleTypeHelper;

template<typename Tremove, typename T>
class TupleRemovesingleTypeHelper<Tremove,std::tuple<T> > 
{
public:
static constexpr bool boolval=std::is_same<T, Tremove>::value;
using type= typename std::conditional<boolval, std::tuple<>, std::tuple<T> >::type;
};

template<typename Tremove, typename T,typename...Ts>
class TupleRemovesingleTypeHelper<Tremove,std::tuple<T,Ts...>> 
{
public:
static constexpr bool boolval=std::is_same<T, Tremove>::value;
using single_type= typename std::conditional<boolval, std::tuple<>, std::tuple<T> >::type;

using type=typename std::conditional<boolval, 
                                     // if we have found the type to remvoe, just add all the other ones
                                     TupleCatType< std::tuple<Ts...> >,
                                     // otherwise add T and check to the next one 
                                     TupleCatType<single_type, typename TupleRemovesingleTypeHelper< Tremove,std::tuple<Ts...> >::type >
                                    >::type;
};

template<typename Tremove, typename...Ts>
using TupleRemovesingleType=typename TupleRemovesingleTypeHelper<Tremove,Ts...>::type;



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// Cat tuple of tuple types: TupleCatType< std::tuple<std::tuple<T11>...>, std::tuple<std::tuple<T21>...>,....  > =
////////                                         std::tuple<std::tuple<T11,T12,...>,std::tuple<T21,T22...>....> > 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<Integer N,Integer Nmax,typename ... Ts>
class TupleOfTupleCatTypeHelper;


template<Integer Nmax, typename ... Ts>
class TupleOfTupleCatTypeHelper<Nmax,Nmax,Ts...>
{
 public:
 using type= std::tuple<TupleCatType<GetType<Ts,Nmax>...>>;
};

template<Integer N,Integer Nmax, typename ... Ts>
class TupleOfTupleCatTypeHelper
{
public:
using type= 
decltype(std::tuple_cat(std::declval<std::tuple<TupleCatType< GetType<Ts,N>...>>>(),
                           std::declval<typename TupleOfTupleCatTypeHelper<N+1,Nmax,Ts...>::type>()));
};

template<Integer Nmax, typename ... Ts>
using TupleOfTupleCatType=typename TupleOfTupleCatTypeHelper<0,Nmax,Ts ... >::type;




template<Integer N,Integer Nmax,typename ... Ts>
class RemoveTupleOfTupleDuplicatesHelper;


template<Integer Nmax, typename ... Ts>
class RemoveTupleOfTupleDuplicatesHelper<Nmax,Nmax,Ts...>
{
 public:
 using type= std::tuple<RemoveTupleDuplicates<TupleCatType<GetType<Ts,Nmax>...>>>;
};

template<Integer N,Integer Nmax, typename ... Ts>
class RemoveTupleOfTupleDuplicatesHelper
{
public:
    //TupleRemovesingleType<std::tuple<>,>
using type= 
decltype(std::tuple_cat(std::declval<std::tuple< RemoveTupleDuplicates<TupleCatType< GetType<Ts,N>...>>>>(),
                           std::declval<typename RemoveTupleOfTupleDuplicatesHelper<N+1,Nmax,Ts...>::type>()));
};


template<typename T,typename ... Ts>
using RemoveTupleOfTupleDuplicates=typename RemoveTupleOfTupleDuplicatesHelper<0,std::tuple_size<T>::value-1,T,Ts ... >::type;



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// Create tuple of dimension N of given type
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <size_t I,typename T> 
struct tuple_n{
    template< typename...Args> using type = typename tuple_n<I-1, T>::template type<T, Args...>;
};

template <typename T> 
struct tuple_n<0, T> {
    template<typename...Args> using type = std::tuple<Args...>;   
};
template <size_t I,typename T>  using TupleOfType = typename tuple_n<I,T>::template type<>;







template <class T, class Tuple>
struct TypeToTupleElementPosition;

template <class T, class... Types>
struct TypeToTupleElementPosition<T, std::tuple<T, Types...>> {
    static const std::size_t value = 0;
};

template <class T, class U, class... Types>
struct TypeToTupleElementPosition<T, std::tuple<U, Types...>> {
    static const std::size_t value = 1 + TypeToTupleElementPosition<T, std::tuple<Types...>>::value;
};







template <typename = void>
constexpr std::size_t TupleSizeHelper ()
 { return 0u; }

template <std::size_t I0, std::size_t ... Is>
constexpr std::size_t TupleSizeHelper ()
 { return I0 + TupleSizeHelper<Is...>(); }

template <typename ... Ts>
constexpr std::size_t TupleSize (std::tuple<Ts...> const &)
 { return TupleSizeHelper<sizeof(Ts)...>(); }


template <typename T,typename ... Types>
class TupleTypeSize;

template <>
class TupleTypeSize<std::tuple<>>
{
public:
 static constexpr std::size_t value = 1;
};

template <typename T>
class TupleTypeSize<std::tuple<T>>
{
public:
 static constexpr std::size_t value = 1;
};

template <typename T,typename ... Types>
class TupleTypeSize<std::tuple<T,Types...>>
{
public:
 static constexpr std::size_t value = 1 + TupleTypeSize<std::tuple<Types...>>::value;
};


template <typename T,typename ... Types>
class TupleTypeSize<std::tuple<T,Types...>&>
{
public:
 static constexpr std::size_t value = TupleTypeSize<std::tuple<T,Types...>>::value;
};

// for 
template<Integer N,typename All, typename Unique>
 class ElementPosition
{public:
 static_assert(N>=0, " negative integer");
 static_assert(N<TupleTypeSize<All>::value, " exceeding tuple dimension");
 using AllType=GetType<All,N>;
 static constexpr Integer value=TypeToTupleElementPosition<AllType,Unique>::value;
};



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// Find the subtuple from Nmin to Nmax of the given std::tuple<Ts...>.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<Integer N,Integer Nmin,Integer Nmax,typename...Ts>
class SubTupleHelper;

template<Integer N,Integer Nmin,Integer Nmax>
class SubTupleHelper<N,Nmin,Nmax,std::tuple<>>
{
public:
 using type=std::tuple<>;
};

template<Integer N,Integer Nmin,Integer Nmax,typename T,typename...Ts>
class SubTupleHelper<N,Nmin,Nmax,std::tuple<T,Ts...>>
{
public:
 using single_type=typename std::conditional< (N<Nmin || N>Nmax), std::tuple<>, std::tuple<T> >::type;
 using type=TupleCatType< single_type, typename SubTupleHelper<N+1,Nmin,Nmax,std::tuple<Ts...> >::type >;
};

template<Integer Nmin,Integer Nmax,typename...Ts>
using SubTupleType=typename SubTupleHelper<0,Nmin,Nmax,Ts... >::type;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// Find the subtuple starting from Nstart and containing all the elements, shifted of Nshift 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<Integer N,Integer Nstart,Integer Nshift,Integer Nmax,typename Tuple>
class SubTupleShiftHelper;

template<Integer Nstart,Integer Nshift,Integer Nmax,typename Tuple>
class SubTupleShiftHelper<Nmax,Nstart,Nshift,Nmax,Tuple>
{
public:
 static_assert(Nstart<TupleTypeSize<Tuple>::value," In SubTupleShiftHelper, Nstart exceeds tuple length");
 static constexpr Integer M=Nstart+Nmax*Nshift;
 using type=std::tuple< GetType<Tuple,M> >;
};

template<Integer N,Integer Nstart,Integer Nshift,Integer Nmax,typename Tuple>
class SubTupleShiftHelper
{
public:
 static constexpr Integer M=Nstart+N*Nshift;
 static_assert(Nstart<TupleTypeSize<Tuple>::value," In SubTupleShiftHelper, Nstart exceeds tuple length");
 static_assert(N<Nmax," In SubTupleShiftHelper, M exceeds tuple length");
 using single_type=std::tuple< GetType<Tuple,M> >;
 using type=TupleCatType< single_type, typename SubTupleShiftHelper<N+1,Nstart,Nshift,Nmax,Tuple >::type >;
};

template<Integer Nstart,Integer Nshift,typename Tuple>
using SubTupleShift=typename SubTupleShiftHelper<0,Nstart,Nshift, (TupleTypeSize<Tuple>::value-1-Nstart)/Nshift,Tuple>::type;




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// Add Tadd in std::tuple<Ts...> in position N.
///// For example, Tadd=float, std::tuple<char,int>:
/////                          N=0: std::tuple<float,char,int>:
/////                          N=1: std::tuple<char,float,int>:
/////                          N=2: std::tuple<char,int,float>:
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<Integer N,typename Tadd,typename...Ts>
class TupleAddingTypeHelper;


// template<typename Tadd,typename...Ts>
// class TupleAddingTypeHelper<0,Tadd,std::tuple<Ts...>> 
// {
// public:

//  using type=std::tuple<Tadd,Ts...>;
// };



template<Integer N,typename Tadd,typename...Ts>
class TupleAddingTypeHelper<N,Tadd,std::tuple<Ts...>> 
{
public:
 using typeLeft=typename std::conditional< (N==0), std::tuple<>, SubTupleType<0,N-1 , Ts... > >::type;
 using typeRight=SubTupleType<N,sizeof...(Ts)+1, Ts...>;//typename std::conditional< (N==sizeof...(Ts)+1), std::tuple<>, SubTupleType<N,sizeof...(Ts)+1, Ts...> >::type;
 using type=TupleCatType<typeLeft,std::tuple<Tadd>,typeRight>;
};

template<Integer N,typename Tadd, typename...Ts>
using TupleAddingType=typename TupleAddingTypeHelper< N,Tadd,std::tuple<Ts...> >::type;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// Change the N-th type of the tuple to type Tchange
///// For example, Tchange=float, std::tuple<char,int>:
/////                          N=0: std::tuple<float,int>:
/////                          N=1: std::tuple<char,float>:
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<Integer N,typename Tchange,typename...Ts>
class TupleChangeTypeHelper;

template<Integer N,typename Tchange,typename...Ts>
class TupleChangeTypeHelper<N,Tchange,std::tuple<Ts...>> 
{
public:
 using typeLeft=SubTupleType<0,N-1, Ts...>;//typename std::conditional< (N==0), std::tuple<>, SubTupleType<0,N-1 , Ts... > >::type;
 using typeRight=SubTupleType<N+1,sizeof...(Ts)+2, Ts...>;//typename std::conditional< (N==sizeof...(Ts)+1), std::tuple<>, SubTupleType<N,sizeof...(Ts)+1, Ts...> >::type;
 using type=TupleCatType<typeLeft,std::tuple<Tchange>,typeRight>;
};

template<Integer N,typename Tchange, typename...Ts>
using TupleChangeType=typename TupleChangeTypeHelper< N,Tchange,std::tuple<Ts...> >::type;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// Change to Tchange type all the N-th types of the tuples contained in the main tuple
///// For N=1, Tchange=char, std::tuple< std::tuple<  std::tuple<int,std::tuple<>>  >, 
/////                                        std::tuple<  std::tuple<int,char>, std::tuple<int,int>>, 
/////                                        std::tuple<  std::tuple<double,std::tuple<>>, std::tuple<std::tuple<>,std::tuple<>> > >:
///// The new tuple becomes:
/////                            std::tuple< std::tuple<  std::tuple<int,char>  >, 
/////                                        std::tuple<  std::tuple<int,char>, std::tuple<int,char>>, 
/////                                        std::tuple<  std::tuple<double,char>, std::tuple<std::tuple<>,char> > >:
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////










template<Integer N,Integer Nmax,typename Tuple,Integer Nchange,typename Tchange>
class provahelper;

template<Integer Nmax,typename Tuple,Integer Nchange,typename Tchange>
class provahelper<Nmax,Nmax,Tuple,Nchange,Tchange>
{
public:
using type=std::tuple<TupleChangeType<Nchange,Tchange,GetType<Tuple,Nmax>>>;
};

template<Integer N,Integer Nmax,typename Tuple,Integer Nchange,typename Tchange>
class provahelper
{
 public:
 using single_type=std::tuple<TupleChangeType<Nchange,Tchange,GetType<Tuple,N>>>;

 using type=decltype(std::tuple_cat(std::declval<single_type>(),
                           std::declval<typename provahelper<N+1,Nmax,Tuple,Nchange,Tchange>::type>()));


};

template<typename Tuple,Integer Nchange,typename Tchange>
using prova=typename provahelper<0,TupleTypeSize<Tuple>::value-1,Tuple,Nchange,Tchange>::type;



template<Integer N,Integer Nmax,Integer Nchange,typename Tchange, typename TupleOfTuple>
class TupleOfTupleChangeTypeHelper;


template<Integer Nmax,Integer Nchange,typename Tchange, typename TupleOfTuple>
class TupleOfTupleChangeTypeHelper<Nmax,Nmax,Nchange,Tchange,TupleOfTuple>
{
 public:
 using T=GetType<TupleOfTuple,Nmax>;
 using type=typename std::conditional< std::is_same< T,std::tuple<> >::value , std::tuple<T>, std::tuple<prova<T,Nchange,Tchange >>>::type;
};



template<Integer N,Integer Nmax,Integer Nchange,typename Tchange, typename TupleOfTuple>
class TupleOfTupleChangeTypeHelper
{
public:
using T=GetType<TupleOfTuple,N>;
using single_type=typename std::conditional< std::is_same< T,std::tuple<> >::value , std::tuple<T>, std::tuple<prova<T,Nchange,Tchange >>>::type;
using type=decltype(std::tuple_cat(std::declval<single_type>(),
                           std::declval<typename TupleOfTupleChangeTypeHelper<N+1,Nmax,Nchange,Tchange,TupleOfTuple>::type>()));
};


template<Integer Nchange,typename Tchange,typename TupleOfTuple>
using TupleOfTupleChangeType=typename TupleOfTupleChangeTypeHelper<0,TupleTypeSize<TupleOfTuple>::value-1,Nchange,Tchange,TupleOfTuple>::type;









template<Integer N,Integer Nmax, typename Tuple>
class Tuple2ToTuple1Helper;


template<Integer Nmax, typename Tuple>
class Tuple2ToTuple1Helper<Nmax,Nmax,Tuple>
{
 public:
 using type=std::tuple<
 std::tuple<
          GetType<GetType<Tuple,Nmax>,0>,
 typename GetType<GetType<Tuple,Nmax>,1>::Elem
 >
 > ;
};



template<Integer N,Integer Nmax, typename Tuple>
class Tuple2ToTuple1Helper
{
public:
using single_type=std::tuple<
 std::tuple<
          GetType<GetType<Tuple,N>,0>,
 typename GetType<GetType<Tuple,N>,1>::Elem >
 >;
using rest=typename Tuple2ToTuple1Helper<N+1,Nmax,Tuple>::type;
using type=decltype(std::tuple_cat(std::declval<single_type>(),
                                   std::declval<rest>()));
};

template<typename Tuple>
using Tuple2ToTuple1=typename Tuple2ToTuple1Helper<0,TupleTypeSize<Tuple>::value-1,Tuple>::type;





template<Integer N,Integer Nmax, typename TupleOfTuple>
class TupleOfTupleRemoveQuadratureHelper;


template<Integer Nmax, typename TupleOfTuple>
class TupleOfTupleRemoveQuadratureHelper<Nmax,Nmax,TupleOfTuple>
{
 public:
 using T=GetType<TupleOfTuple,Nmax>;
 using type=std::tuple<RemoveTupleDuplicates<Tuple2ToTuple1<T>>>;
};



template<Integer N,Integer Nmax, typename TupleOfTuple>
class TupleOfTupleRemoveQuadratureHelper
{
public:
using T=GetType<TupleOfTuple,N>;
using single_type=std::tuple<RemoveTupleDuplicates<Tuple2ToTuple1<T>>>;
using type=decltype(std::tuple_cat(std::declval<single_type>(),
                           std::declval<typename TupleOfTupleRemoveQuadratureHelper<N+1,Nmax,TupleOfTuple>::type>()));
};


template<typename TupleOfTuple>
using TupleOfTupleRemoveQuadrature=typename TupleOfTupleRemoveQuadratureHelper<0,TupleTypeSize<TupleOfTuple>::value-1,TupleOfTuple>::type;














template<typename TupleOfSpaces, Integer Nmax,Integer N>
class UniqueFEFamiliesHelper;

template<typename TupleOfSpaces, Integer Nmax>
class UniqueFEFamiliesHelper<TupleOfSpaces,Nmax,Nmax>
{
 public:

  using Space=GetType<GetType<TupleOfSpaces,Nmax>,1>;
  using type=std::tuple< Number<Space::FEFamily> >;
};

template<typename TupleOfSpaces, Integer Nmax,Integer N=0>
class UniqueFEFamiliesHelper
{
 public:
  using Space=GetType<GetType<TupleOfSpaces,N>,1>;
  using single_type=std::tuple< Number<Space::FEFamily> >;
  using type=decltype(std::tuple_cat(std::declval<single_type>(),
                             std::declval<typename UniqueFEFamiliesHelper<TupleOfSpaces,Nmax,N+1>::type>()));
};


template<typename TupleOfSpaces>
using UniqueFEFamilies=RemoveTupleDuplicates<typename UniqueFEFamiliesHelper<TupleOfSpaces,TupleTypeSize<TupleOfSpaces>::value-1,0>::type>;


template <typename TupleOfNumbers,Integer Nmax,Integer N>
class MaxNumberInTupleHelper;

template<typename TupleOfSpaces,Integer Nmax,Integer N>
class SpacesToUniqueFEFamiliesHelper;

template<typename TupleOfSpaces,Integer Nmax>
class SpacesToUniqueFEFamiliesHelper<TupleOfSpaces,Nmax,Nmax>
{
 public:
  using TupleFEFamilies=UniqueFEFamilies<TupleOfSpaces>;
  using FamilyNumber=Number<GetType<GetType<TupleOfSpaces,Nmax>,1>::FEFamily>;
  using PositionNumber=Number<TypeToTupleElementPosition<FamilyNumber,TupleFEFamilies>::value>;
  using type=std::tuple< PositionNumber >;
};

template<typename TupleOfSpaces,Integer Nmax,Integer N=0>
class SpacesToUniqueFEFamiliesHelper
{
 public:
  using TupleFEFamilies=UniqueFEFamilies<TupleOfSpaces>;
  using FamilyNumber=Number<GetType<GetType<TupleOfSpaces,N>,1>::FEFamily>;
  using PositionNumber=Number<TypeToTupleElementPosition<FamilyNumber,TupleFEFamilies>::value>;
  using single_type=std::tuple< PositionNumber  >;
  using type=decltype(std::tuple_cat(std::declval<single_type>(),
                             std::declval<typename SpacesToUniqueFEFamiliesHelper<TupleOfSpaces,Nmax,N+1>::type>()));
};



template<typename TupleOfSpaces>
using SpacesToUniqueFEFamilies=typename SpacesToUniqueFEFamiliesHelper<TupleOfSpaces,TupleTypeSize<TupleOfSpaces>::value-1,0>::type;








template<typename NumberTuple,typename Maps,Integer M, Integer Nmax,Integer N>
class UniqueMapSingleSpaceHelper;

template<typename NumberTuple,typename Maps,Integer M,Integer Nmax>
class UniqueMapSingleSpaceHelper<NumberTuple,Maps,M,Nmax,Nmax>
{
 public:
  using type= typename std::conditional< GetType<NumberTuple,Nmax>::value == M , GetType<Maps,Nmax> , std::tuple<>>::type;
};

template<typename NumberTuple,typename Maps,Integer M,Integer Nmax,Integer N=0>
class UniqueMapSingleSpaceHelper
{
 public:
  using single_type= typename std::conditional< GetType<NumberTuple,N>::value == M , GetType<Maps,N> , std::tuple<>>::type;
  using type=decltype(std::tuple_cat(std::declval<single_type>(),
                             std::declval<typename UniqueMapSingleSpaceHelper<NumberTuple,Maps,M,Nmax,N+1>::type>()));
};

template<typename NumberTuple,typename Maps,Integer M>
using UniqueMapSingleSpace=RemoveTupleDuplicates< typename UniqueMapSingleSpaceHelper<NumberTuple,Maps,M,TupleTypeSize<NumberTuple>::value-1,0>::type>;







template <typename Tuple,typename Type,Integer Nmax,Integer N>
class HowMayTypesDifferentFromTypeHelper;


template <typename Tuple,typename Type,Integer Nmax>
class HowMayTypesDifferentFromTypeHelper<Tuple,Type,Nmax,Nmax>
{
public:
  static constexpr Integer value=IsDifferent<Type,GetType<Tuple,Nmax>>::value;
};

template <typename Tuple,typename Type,Integer Nmax,Integer N>
class HowMayTypesDifferentFromTypeHelper
{
public:
  static constexpr Integer value=IsDifferent<Type,GetType<Tuple,N>>::value+HowMayTypesDifferentFromTypeHelper<Tuple,Type,Nmax,N+1>::value;
};

template <typename Tuple,typename Type>
static constexpr Integer HowMayTypesDifferentFromType=HowMayTypesDifferentFromTypeHelper<Tuple,Type,TupleTypeSize<Tuple>::value-1,0>::value;





template <typename TupleOfNumbers>
class MaxNumberInTuple;

template<typename NumberTuple,typename Maps,Integer Nmax,Integer N>
class UniqueMapHelper;


template<typename NumberTuple,typename Maps,Integer Nmax>
class UniqueMapHelper<NumberTuple,Maps,Nmax,Nmax>
{
 public:
  using type= std::tuple<UniqueMapSingleSpace<NumberTuple,Maps,Nmax>>;
};

template<typename NumberTuple,typename Maps,Integer Nmax,Integer N=0>
class UniqueMapHelper
{
 public:
  using single_type= std::tuple<UniqueMapSingleSpace<NumberTuple,Maps,N>>;
  using type=decltype(std::tuple_cat(std::declval<single_type>(),
                             std::declval<typename UniqueMapHelper<NumberTuple,Maps,Nmax,N+1>::type>()));
};



template<typename NumberTuple,typename Maps>
using UniqueMap=typename UniqueMapHelper<NumberTuple,Maps,MaxNumberInTuple<NumberTuple>::value,0>::type;






template<typename Operator, typename Elem,Integer FEFamily> //  BaseFunctioSpace>
class MapFromReference5;

template<typename Tuple, Integer FEFamily, Integer Nmax,Integer N>
class MapOperatorTupleHelper;

template<typename Tuple, Integer FEFamily, Integer Nmax>
class MapOperatorTupleHelper<Tuple,FEFamily,Nmax,Nmax>
{
 public:

  using Nthelem=GetType<Tuple,Nmax>;
  using Operator=GetType<Nthelem,0>;
  using Elem=GetType<Nthelem,1>;
  using type=std::tuple< MapFromReference5<Operator,Elem,FEFamily> >;
};

template<typename Tuple, Integer FEFamily,Integer Nmax,Integer N=0>
class MapOperatorTupleHelper
{
 public:
  using Nthelem=GetType<Tuple,N>;
  using Operator=GetType<Nthelem,0>;
  using Elem=GetType<Nthelem,1>;
  using single_type=std::tuple< MapFromReference5<Operator,Elem,FEFamily> >;
  using type=decltype(std::tuple_cat(std::declval<single_type>(),
                             std::declval<typename MapOperatorTupleHelper<Tuple,FEFamily,Nmax,N+1>::type>()));
};


template<typename Tuple, Integer FEFamily,Integer Nmax>
using MapOperatorTuple=typename MapOperatorTupleHelper<Tuple,FEFamily,Nmax,0>::type;


template<typename TupleOfTuple, typename TupleSpaces, Integer Nmax,Integer N>
class MapOperatorTupleOfTupleHelper;


template<typename TupleOfTuple,typename TupleSpaces, Integer Nspaces>
class MapOperatorTupleOfTupleHelper<TupleOfTuple,TupleSpaces, Nspaces, Nspaces>
{
public:
using Tuple=GetType<TupleOfTuple,Nspaces>;
static constexpr Integer FEFamily=GetType<GetType<TupleSpaces,Nspaces>,1>::FEFamily;
static constexpr Integer Nmax=TupleTypeSize<Tuple>::value-1;
using type=std::tuple< MapOperatorTuple< Tuple,FEFamily,Nmax> > ;
};

template<typename TupleOfTuple,typename TupleSpaces, Integer Nspaces, Integer N=0>
class MapOperatorTupleOfTupleHelper
{
public:
using Tuple=GetType<TupleOfTuple,N>;
static constexpr Integer FEFamily=GetType<GetType<TupleSpaces,N>,1>::FEFamily;
static constexpr Integer Nmax=TupleTypeSize<Tuple>::value-1;
using single_type=std::tuple< MapOperatorTuple< Tuple, FEFamily, Nmax> > ;
using type=decltype(std::tuple_cat(std::declval<single_type>(),
                           std::declval<typename MapOperatorTupleOfTupleHelper<TupleOfTuple,TupleSpaces,Nspaces,N+1>::type>()));
};

template<typename TupleOfTuple,typename TupleSpaces>
using MapOperatorTupleOfTuple=typename MapOperatorTupleOfTupleHelper<TupleOfTuple,TupleSpaces,TupleTypeSize<TupleSpaces>::value-1,0>::type;







//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// Create Tuple (for each function space) of tuples of shape functions related to the pair of pair<operator,quadrature_rule>
///// The input is the tuple of tuple of pair<operator,quadrature_rule> 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename Elem, typename BaseFunctioSpace, typename Operator, typename QuadratureRule>
class ShapeFunctionDependent;

template<typename Elem, typename BaseFunctioSpace, typename Tuple, Integer Nmax,Integer N>
class TupleShapeFunctionCreate;

template<typename Elem, typename BaseFunctioSpace, typename Tuple, Integer Nmax>
class TupleShapeFunctionCreate<Elem,BaseFunctioSpace,Tuple,Nmax,Nmax>
{
 public:

  using Nelem=GetType<Tuple,Nmax>;
  using Operator=GetType<Nelem,0>;
  using QuadratureRule=GetType<Nelem,1>;
  using type=std::tuple< ShapeFunctionDependent<Elem, BaseFunctioSpace,Operator,QuadratureRule > >;
};



template<typename Elem, typename BaseFunctioSpace, typename Tuple, Integer Nmax,Integer N=0>
class TupleShapeFunctionCreate
{
 public:
  using Nelem=GetType<Tuple,N>;
  using Operator=GetType<Nelem,0>;
  using QuadratureRule=GetType<Nelem,1>;
  using single_type=std::tuple< ShapeFunctionDependent<Elem, BaseFunctioSpace,Operator,QuadratureRule > >;
  using type=decltype(std::tuple_cat(std::declval<single_type>(),
                             std::declval<typename TupleShapeFunctionCreate<Elem,BaseFunctioSpace,Tuple,Nmax,N+1>::type>()));
};


template<typename TupleSpaces, typename TupleOperatorsAndQuadrature, Integer Nmax,Integer N>
class TupleOfTupleShapeFunctionCreate;



template<typename TupleSpaces,typename TupleOperatorsAndQuadrature, Integer Nspaces>
class TupleOfTupleShapeFunctionCreate<TupleSpaces,TupleOperatorsAndQuadrature, Nspaces, Nspaces>
{
public:
using OperatorAndQuadrature=GetType<TupleOperatorsAndQuadrature,Nspaces>;
static constexpr Integer Nmax=TupleTypeSize<OperatorAndQuadrature>::value-1;
using Space=GetType<TupleSpaces,Nspaces>;
using Elem=GetType<Space,0>;
using FunctionSpace=GetType<Space,1>;
using type=std::tuple<typename TupleShapeFunctionCreate< Elem,FunctionSpace, OperatorAndQuadrature,Nmax>::type>;
};

template<typename TupleSpaces,typename TupleOperatorsAndQuadrature, Integer Nspaces, Integer N=0>
class TupleOfTupleShapeFunctionCreate
{
public:
using OperatorAndQuadrature=GetType<TupleOperatorsAndQuadrature,N>;
static constexpr Integer Nmax=TupleTypeSize<OperatorAndQuadrature>::value-1;
using Space=GetType<TupleSpaces,N>;
using Elem=GetType<Space,0>;
using FunctionSpace=GetType<Space,1>;
using single_type=std::tuple<typename TupleShapeFunctionCreate< Elem,FunctionSpace, OperatorAndQuadrature,Nmax>::type>;
using type=decltype(std::tuple_cat(std::declval<single_type>(),
                           std::declval<typename TupleOfTupleShapeFunctionCreate<TupleSpaces,TupleOperatorsAndQuadrature,Nspaces,N+1>::type>()));
};



template<typename TupleSpaces,typename TupleOperatorsAndQuadrature>
using TupleOfTupleShapeFunctionType=typename TupleOfTupleShapeFunctionCreate<TupleSpaces,TupleOperatorsAndQuadrature,TupleTypeSize<TupleSpaces>::value-1,0>::type;



template<typename Operator,typename Tuple,Integer Nmax,Integer N>
class OperatorToMapHelper;

template<typename Operator,typename Tuple,Integer Nmax>
class OperatorToMapHelper<Operator,Tuple,Nmax,Nmax>
{
public:
  using type=Number<Nmax>;
};

template<typename Operator,typename Tuple,Integer Nmax,Integer N>
class OperatorToMapHelper
{
public:
  using OperatorMap=typename GetType<Tuple,N>::Operator;
  using type =typename std::conditional< std::is_same<OperatorMap,Operator>::value, 
                                         Number<N>, 
                                         typename OperatorToMapHelper<Operator,Tuple,Nmax,N+1>::type 
                                       >::type ;
};

template<typename Operator,typename Tuple>
using OperatorToMap=typename OperatorToMapHelper<Operator,Tuple,TupleTypeSize<Tuple>::value-1,0>::type;





template<typename TupleSpaces,typename Map, Integer Nmax, Integer N>
class MapTupleInit2Helper;

template<typename TupleSpaces,typename Map, Integer Nmax>
class MapTupleInit2Helper<TupleSpaces,Map, Nmax, Nmax>
{
public:
      using operator_space=typename GetType<TupleSpaces,Nmax>::Operator;
      using type=std::tuple<OperatorToMap<operator_space,Map>>;
};

template<typename TupleSpaces,typename Map, Integer Nmax, Integer N>
class MapTupleInit2Helper
{
public:
      using operator_space=typename GetType<TupleSpaces,N>::Operator;
      using single_type=std::tuple<OperatorToMap<operator_space,Map>>;
      using type=decltype(std::tuple_cat(std::declval<single_type>(),
                                 std::declval<typename MapTupleInit2Helper<TupleSpaces,Map,Nmax,N+1>::type>()));
};

template<typename TupleSpaces,typename Map>
using MapTupleInit2=typename MapTupleInit2Helper<TupleSpaces,Map,TupleTypeSize<TupleSpaces>::value-1,0>::type;




template<typename TupleOfTupleSpaces,typename SpaceToMap, typename TupleOfMaps, Integer Nmax, Integer N>
class MapTupleInitHelper;

template<typename TupleOfTupleSpaces,typename SpaceToMap,typename TupleOfMaps, Integer Nmax>
class MapTupleInitHelper<TupleOfTupleSpaces,SpaceToMap, TupleOfMaps,Nmax, Nmax>
{
public:
      using TupleSpaces=GetType<TupleOfTupleSpaces,Nmax>;
      using Maps=GetType< TupleOfMaps, GetType<SpaceToMap,Nmax>::value>;
      using type=std::tuple< MapTupleInit2<TupleSpaces,Maps> >;
};

template<typename TupleOfTupleSpaces,typename SpaceToMap,typename TupleOfMaps, Integer Nmax, Integer N>
class MapTupleInitHelper
{
public:
      using TupleSpaces=GetType<TupleOfTupleSpaces,N>;
      using Maps=GetType<TupleOfMaps,GetType<SpaceToMap,N>::value>;
      using single_type=std::tuple< MapTupleInit2<TupleSpaces,Maps> >;
      using type=decltype(std::tuple_cat(std::declval<single_type>(),
                                 std::declval<typename MapTupleInitHelper<TupleOfTupleSpaces,SpaceToMap,TupleOfMaps,Nmax,N+1>::type>()));
};


template<typename TupleOfTupleSpaces,typename SpaceToMap,typename TupleOfMaps>
using MapTupleInit=typename MapTupleInitHelper<TupleOfTupleSpaces,SpaceToMap,TupleOfMaps,TupleTypeSize<TupleOfTupleSpaces>::value-1,0>::type;









  
  template<typename FunctionSpaces>
  class ShapeFunctionCoefficient;



  template<typename Elem, Integer FEFamily,Integer Order, Integer Ndofs,typename Shape>
  void shape_function_init_aux_aux_aux(Shape& shape, const Vector<Real,Ndofs> &alpha);

  template<typename Ele, Integer FEFamily,Integer Order, Integer Ndofs,typename Shape>
  typename std::enable_if_t< FEFamily==LagrangeFE,void> 
  shape_function_init_aux_aux_aux(Shape& shape, const Vector<Real,Ndofs> &alpha)
  {
    shape.init();
  }


  template<typename Elem, Integer FEFamily,Integer Order, Integer Ndofs,typename Shape>
  typename std::enable_if_t<FEFamily==RaviartThomasFE,void> 
  shape_function_init_aux_aux_aux(Shape& shape, const Vector<Real,Ndofs> &alpha)
  {
    shape.init(alpha);
  }




  template<Integer Nmax,Integer N,Integer Ndofs,typename Tuple>
  typename std::enable_if_t<(N>Nmax),void> 
  shape_function_init_aux_aux(Tuple& tuple, const Vector<Real,Ndofs> &alpha)
  {}

  template<Integer Nmax,Integer N,Integer Ndofs,typename Tuple>
  typename std::enable_if_t< (N<=Nmax),void> 
  shape_function_init_aux_aux(Tuple& tuple, const Vector<Real,Ndofs> &alpha)
  {
    shape_function_init_aux_aux<Nmax,N+1>(tuple,alpha);
    tuple_get<N>(tuple).init(alpha);
  }

  template<Integer Nmax,Integer N, typename Tuple,typename...Args>
  typename std::enable_if_t< (N>Nmax),void> 
  shape_function_init_aux(Tuple& tuple, const ShapeFunctionCoefficient<Args...> &coefficients)
  {}

  template<Integer Nmax,Integer N, typename Tuple,typename...Args>
  typename std::enable_if_t<N<=Nmax,void> 
  shape_function_init_aux(Tuple& tuple, const ShapeFunctionCoefficient<Args...> &coefficients)
  {
    constexpr Integer Nmax_aux=TupleTypeSize<GetType<Tuple,N>>::value-1;
    // shape_function_init_aux_aux<Nmax_aux,0>(tuple_get<N>(tuple),coefficients);
    shape_function_init_aux<Nmax,N+1>(tuple,coefficients);
  }



  template< typename Tuple,typename...Args>
  constexpr void init(Tuple& tuple, const ShapeFunctionCoefficient<Args...> &coefficients)
  {
   constexpr Integer Nmax=TupleTypeSize<Tuple>::value-1;

   // CORREGGI, TOGLI NMAX-1 CAMBIA IN NMAX
   shape_function_init_aux<Nmax,0>(tuple,coefficients);
  }



















template<typename T,Integer Nmax,Integer N>
class UniqueElementFEFamilyHelper;

template<typename T,Integer Nmax>
class UniqueElementFEFamilyHelper<T,Nmax,Nmax>
{
public:
  using TN=GetType<T,Nmax>;
  using type=std::tuple<std::tuple<GetType<TN,0>,Number<GetType<TN,1>::FEFamily>>>;
};

template<typename T,Integer Nmax,Integer N>
class UniqueElementFEFamilyHelper
{
 public:
  using TN=GetType<T,N>;
  using single_type=std::tuple<std::tuple<GetType<TN,0>,Number<GetType<TN,1>::FEFamily>>>;
  using type=decltype(std::tuple_cat(std::declval<single_type>(),
                      std::declval<typename UniqueElementFEFamilyHelper<T,Nmax,N+1>::type>()));
};
template<typename T>
using UniqueElementFEFamily=RemoveTupleDuplicates<typename UniqueElementFEFamilyHelper<T,TupleTypeSize<T>::value-1,0>::type>;




// template<typename TupleSpaces,typename TupleOperatorsAndQuadrature,typename UniqueMapping>
// class MapTuple
// {
// public:
//       using TupleOfTupleNoQuadrature=TupleOfTupleRemoveQuadrature<TupleOperatorsAndQuadrature>;
//       using Unique=SpacesToUniqueFEFamilies<TupleSpaces>;
//       using Map=MapOperatorTupleOfTuple<TupleOfTupleNoQuadrature,TupleSpaces>;
//       using UniqueMap=UniqueMap<Unique,Map> ;
      
//       template<Integer Nmax,Integer N,Integer Mmax,Integer M,typename Jacobian>
//       typename std::enable_if<(N>Nmax), void >::type
//       init_aux(const Jacobian& J){}

//       template<Integer Nmax,Integer N,Integer Mmax,Integer M,typename Jacobian>
//       typename std::enable_if<(N<=Nmax && M==Mmax), void >::type
//       init_aux(const Jacobian& J)
//       {
//         const auto map_n_m_th=std::get<M>(std::get<N>(tuple_of_maps_));
//         map_n_m_th.init(J);
//         init_aux<Nmax,N+1,0,0> (J);
//       }

//       template<Integer Nmax,Integer N,Integer Mmax,Integer M,typename Jacobian>
//       typename std::enable_if<(N<=Nmax && M<Mmax), void >::type
//       init_aux(const Jacobian& J)
//       {
//         const auto map_n_m_th=std::get<M>(std::get<N>(tuple_of_maps_));
//         map_n_m_th.init(J);
//         init_aux<Nmax,N,TupleTypeSize<GetType<UniqueMapping,N>>::value-1,M+1> (J);
//       }


//       // TODO: IN CASE THE MAP HAS DIFFERENT ELEMENTS, THEN ALSO DIFFERENT JACOBIANS ARE NEEDED
//       template<Integer Dim, Integer ManifoldDim>
//       void init(const Matrix<Real, Dim, ManifoldDim>& J)
//       {init_aux<TupleTypeSize<UniqueMapping>::value-1,0,0,0>(J);}

//       UniqueMapping operator()(){return tuple_of_maps_;}
//     private:
//       UniqueMapping tuple_of_maps_;
// };



// template<typename TupleSpaces,typename TupleOperatorsAndQuadrature,typename UniqueMapping>
// class ShapeFunctionTuple
// {
//  using TupleOfTupleShapeFunction=TupleOfTupleShapeFunctionType<TupleSpaces,TupleOperatorsAndQuadrature>;
//  using SpacesToUniqueSpaces=SpacesToUniqueFEFamilies<TupleSpaces>;
//  using MapTupleNumbersW1=MapTupleInit<TupleOfTupleShapeFunction, SpacesToUniqueSpaces , UniqueMapping>;
// };



// template<typename TupleSpaces,typename TupleOperatorsAndQuadrature,typename Maps>
// class ShapeFunctionAssembly
// {
//  public:
//    using type=TupleOfTupleShapeFunctionType<TupleSpaces,TupleOperatorsAndQuadrature>;
//    static constexpr Integer Nmax=TupleTypeSize<TupleSpaces>::value-1;









//    ShapeFunctionAssembly(const Maps& maps)
//    {
//     //init_spaces(maps)
//    }




//    // /////////////// reference initialization of shape functions of all spaces, all operators and all quadrature rules used 

//    // template<Integer Nmax_aux,Integer N,typename Tuple>
//    // typename std::enable_if< (N>Nmax_aux) ,void>::type 
//    // init_reference_aux_aux(Tuple& t){}

//    // template<Integer Nmax_aux,Integer N,typename Tuple>
//    // typename std::enable_if< (N<=Nmax_aux) ,void>::type 
//    // init_reference_aux_aux(Tuple& t) 
//    // {auto& t_nth=std::get<N>(t); 
//    //  t_nth.init_reference();
//    //  init_reference_aux_aux<Nmax_aux,N+1,Tuple>(t);}

//    // template<Integer N>
//    // typename std::enable_if< (N>Nmax) ,void>::type 
//    // init_reference_aux(){}

//    // template<Integer N>
//    // typename std::enable_if< (N<=Nmax) ,void>::type 
//    // init_reference_aux() 
//    // {auto& tuple_nth=std::get<N>(tuple_);
//    //  using tuple_nth_type=decltype(tuple_nth);
//    //  init_reference_aux_aux<TupleTypeSize<tuple_nth_type>::value-1,0,tuple_nth_type>(tuple_nth);
//    //  init_reference_aux<N+1>(); }

//    // void init_reference(){init_reference_aux<0>();}
   




//    // template<Integer Nmax_aux,Integer N,typename Tuple,typename Jacobian>
//    // typename std::enable_if< (N>Nmax_aux) ,void>::type 
//    // init_aux_aux(const Jacobian&J,Tuple& t){}

//    // template<Integer Nmax_aux,Integer N,typename Tuple,typename Jacobian>
//    // typename std::enable_if< (N<=Nmax_aux) ,void>::type 
//    // init_aux_aux(const Jacobian&J,Tuple& t) 
//    // {auto& t_nth=std::get<N>(t); 
//    //  t_nth.init(J);
//    //  init_aux_aux<Nmax_aux,N+1,Tuple,Jacobian>(J,t);}

//    // template<Integer N,typename Jacobian>
//    // typename std::enable_if< (N>Nmax) ,void>::type 
//    // init_aux(const Jacobian&J){}

//    // template<Integer N,typename Jacobian>
//    // typename std::enable_if< (N<=Nmax) ,void>::type 
//    // init_aux(const Jacobian&J) 
//    // {auto& tuple_nth=std::get<N>(tuple_);
//    //  using tuple_nth_type=decltype(tuple_nth);
//    //  init_aux_aux<TupleTypeSize<tuple_nth_type>::value-1,0,tuple_nth_type>(J,tuple_nth);
//    //  init_aux<N+1>(J); }

//    // template<typename Jacobian>
//    // void init(const Jacobian&J){init_aux<0>(J);}




//    const type& get()const{return tuple_;};
//          type&  get()     {return tuple_;};

//    template<Integer N,Integer M>
//    const GetType<type,N,M> & get()const
//    {const auto& tuple_nth=std::get<N>(tuple_);
//     return std::get<M>(tuple_nth);};

//    template<Integer N,Integer M>
//          GetType<type,N,M>&  get()    
//    { auto& tuple_nth=std::get<N>(tuple_);
//     return std::get<M>(tuple_nth);};




//  private:
//    type tuple_;
// };




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// Compile-time Max, Min, Equal, Greater, Lesser
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
constexpr T Max (const T& a,const T& b) 
{
  return a > b ? a : b;
}

template<typename T>
constexpr T Min (const T& a,const T& b) 
{
  return a < b ? a : b;
}

template<typename T,typename S>
constexpr bool Equal (const T& a,const S& b) 
{
  return a == b ? 1 : 0;
}


template<typename T,typename S>
constexpr bool Greater (const T& a,const S& b) 
{
  return a > b ? 1 : 0;
}

template<typename T,typename S>
constexpr bool Lesser (const T& a,const S& b) 
{
  return a < b ? 1 : 0;
}





















}

#endif