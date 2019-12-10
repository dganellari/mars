


#ifndef MARS_TUPLE_UTILITIES_HPP
#define MARS_TUPLE_UTILITIES_HPP

#include "mars_base.hpp"
#include "mars_static_math.hpp"
#include "mars_array.hpp"
#include "mars_compiletime_sqrt.hpp"


namespace mars{


class EmptyClass;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// Compile-time Max, Min, Equal, Greater, Lesser
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
constexpr T Max(const T&t)
{
  return t;
}


template<typename T>
constexpr T Max (const T& a,const T& b) 
{
  return a > b ? a : b;
}

template<typename T,typename...Ts>
constexpr T Max(const T&t,const Ts&...ts)
{
  return Max(t,Max(ts...));
}



template<typename T>
constexpr T Min (const T& a) 
{
  return a ;
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

template<typename T>
constexpr T Sum(const T&t)
{
  return t;
}


template<typename T,typename...Ts>
constexpr T Sum(const T&t,const Ts&...ts)
{
  return t+Sum(ts...);
}


// returns the sign of val (1,0,-1)
template <typename T> 
constexpr Integer Sign(const T& val) {
    return (T(0) < val) - (val < T(0));
}

// returns the sign of val (1,0,-1)
template <typename T, Integer Dim> 
constexpr Array<T,Dim> Sign(const Array<T,Dim>& val) {
    Array<T,Dim> output;
    for(Integer ii=0;ii<Dim;ii++)
       output[ii]=(T(0) < val[ii]) - (val[ii] < T(0));
    return output;
}


template<typename T>
constexpr T Heaviside (const T& a) 
{
  return a > 0 ? a : 0;
}

template<typename T>
constexpr bool IsNegative (const T& a) 
{
  return a < 0 ? true : false;
}

template<typename T>
constexpr bool IsPositive (const T& a) 
{
  return a > 0 ? true : false;
}


template<typename T, Integer Dim >
constexpr bool IsPositive (const Array<T,Dim>& a) 
{
  bool output=true;
  for(Integer ii=0;ii<Dim;ii++)
      if(IsNegative(a[ii]))
       {output=false;
        break;}
  return output;
}


template <typename T>
T argsort(const T &v) {
    // initialize original index locations
    T idx;
    std::iota(idx.begin(), idx.end(), 0);
    // sort indexes based on comparing values in v
    std::sort(idx.begin(), idx.end(),
              [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
    
    return idx;
}

// template <typename T>
// void argsort(T& idx,const T &v) {
//     // initialize original index locations
//     std::iota(idx.begin(), idx.end(), 0);
//     // sort indexes based on comparing values in v
//     std::sort(idx.begin(), idx.end(),
//               [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
    
// }
// template <typename T>
// T argsort(const T &v) {
//     T idx;
//     argsort(idx,v);
//     return idx;
// }

// returns v, resorted based on index 
template <typename T, typename S>
T sort_by_index(const T &v,const S& index) {
    T sortv;
    assert(index.size()==v.size() && "sort_by_index: v and index must have the same length, otherwise last elements in v are not initialized");
    for(Integer ii=0;ii<index.size();ii++)
        sortv[ii]=v[index[ii]];
    return sortv;
}
// returns the sub array of v of indices = index
template <typename S,std::size_t SDim, typename T,std::size_t TDim>
std::array<S,TDim> sub_array(const std::array<S,SDim> &v,const std::array<T,TDim>& index) {
    
    static_assert(TDim<=SDim,"in sub_array the length of the index vector must be smaller than the one of the vector");
    std::array<S,TDim> subvec;
    for(Integer ii=0;ii<TDim;ii++)
        subvec[ii]=v[index[ii]];
    
    return subvec;
}

template<Integer N,Integer K>
 constexpr void combinations_generate_aux(
            Array<Integer, K> &data,
            const Integer index, 
            const Integer i,
            Array<Array<Integer, K>, binomial_coefficient(N,K)> &combs,
            Integer &comb_index)
        {
            if(index == K) {
                for(Integer ii=0;ii<data.size();ii++)
                    combs[comb_index][ii]=data[ii];
                comb_index++;
                return;
            }

            if(i >= N) {
                return;
            }

            data[index] = i;

            combinations_generate_aux<N,K>(data, index+1, i+1, combs, comb_index);
            
            // current is excluded, replace it with next (Note that
            // i+1 is passed, but index is not changed)
            combinations_generate_aux<N,K>(data, index, i+1, combs, comb_index);
        }

        template<Integer N,Integer K >
        constexpr Array<Array<Integer, K>, binomial_coefficient(N,K)> combinations_generate()
        {
            Array<Array<Integer, K>, binomial_coefficient(N,K)> combs;
            Array<Integer, K> data;
            Integer comb_index = 0;
            combinations_generate_aux<N,K>(data, 0, 0, combs, comb_index);
            return combs;
        }
        

template<typename Left,typename Right>
class ChooseHelper
{
 public: 
  using type=typename std::conditional<IsSame<Left,EmptyClass>::value, Right,Left>::type;
};

template<typename Left,typename Right>
using Choose=typename ChooseHelper<Left,Right>::type;


template<typename ...Args>
class TupleOfSubTypesHelper;


template<typename ...Args>
class TupleOfSubTypesHelper<std::tuple<Args...>>
{
 public:
  using type=std::tuple<typename Args::value_type...> ;
};

template<typename ...Args>
class TupleOfSubTypesHelper<const std::tuple<Args...>>
{
 public:
  using type=std::tuple<typename Args::value_type...> ;
};
template<typename ...Args>
using TupleOfSubTypes=typename TupleOfSubTypesHelper<std::remove_const<std::tuple<Args...>>>::type;




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

template <Integer N>
class GetHelper<N, std::tuple<>>
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

template <Integer N>
class GetHelper<N, const std::tuple<>>
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


template <typename Tuple>
class GetTypeHelper<Tuple,-1>
{
public:
  using type =std::tuple<>;
};

template <Integer N>
class GetTypeHelper<std::tuple<>,N>
{
public:
  using type =std::tuple<>;
};

template <typename Tuple,Integer N=0, Integer...Ns>
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


// template <std::size_t N, typename... Ts, typename T>
// auto change_tuple_element(const std::tuple<Ts...>& tuple,const T& t)
// {
//  return std::tuple_cat(sub_tuple<>);
// }



// template<typename...Args>
// std::tuple<Args...> add_costant (const std::tuple<Args...>& t1,const Integer& value);


template<Integer Nelem_dofs>
auto add_costant
(std::vector<Array<Integer, Nelem_dofs>>& t2,const std::vector<Array<Integer, Nelem_dofs>>& t1,const Integer& value)
{
  const auto& t1size=t1.size();
  t2.resize(t1size);

  for(Integer ii=0;ii<t1size;ii++)
    for(Integer jj=0;jj<t1[ii].size();jj++)
       t2[ii][jj]=t1[ii][jj]+value;
  return t2;
}

template<Integer N,typename...Args>
typename std::enable_if<sizeof...(Args)==N+1, void>::type
 add_costant_tuple_loop(const std::tuple<Args...>& t1, std::tuple<Args...>& t2,const Integer& value)
{
  add_costant(std::get<N>(t2),std::get<N>(t1),value);
}

template<Integer N,typename...Args>
typename std::enable_if< N+1<sizeof...(Args), void>::type
 add_costant_tuple_loop(const std::tuple<Args...>& t1, std::tuple<Args...>& t2,const Integer& value)
{
  add_costant(std::get<N>(t2),std::get<N>(t1),value);
  add_costant_tuple_loop<N+1>(t1,t2,value);

}


template<typename...Args>
std::tuple<Args...> add_costant (const std::tuple<Args...>& t1,const Integer& value)
{
  std::tuple<Args...> t2=t1;
  add_costant_tuple_loop<0>(t1,t2,value);

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
    // using fixme=RemoveTupleDuplicates<TupleCatType< GetType<Ts,Nmax>...>>;

  // using type= std::tuple<std::tuple<int>>;
};

template<Integer N,Integer Nmax, typename ... Ts>
class RemoveTupleOfTupleDuplicatesHelper
{
public:
//     //TupleRemovesingleType<std::tuple<>,>
//   using fixme=RemoveTupleDuplicates<TupleCatType< GetType<Ts,N>...>>;

// using type= 
// decltype(std::tuple_cat(std::declval<std::tuple<fixme>>(),
//                            std::declval<typename RemoveTupleOfTupleDuplicatesHelper<N+1,Nmax,Ts...>::type>()));
using type= 
decltype(std::tuple_cat(std::declval<std::tuple< RemoveTupleDuplicates<TupleCatType< GetType<Ts,N>...>>>>(),
                           std::declval<typename RemoveTupleOfTupleDuplicatesHelper<N+1,Nmax,Ts...>::type>()));

};


template<typename T,typename ... Ts>
using RemoveTupleOfTupleDuplicates=typename RemoveTupleOfTupleDuplicatesHelper<0,Max(std::tuple_size<T>::value,std::tuple_size<Ts>::value...)-1,T,Ts ... >::type;
// template<typename T,typename ... Ts>
// using RemoveTupleOfTupleDuplicates=typename RemoveTupleOfTupleDuplicatesHelper<0,std::tuple_size<T>::value-1,T,Ts ... >::type;



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





//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// TupleTypeSize<Tuple>::value== number of elements of the tuple
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

template <>
class TupleTypeSize<const std::tuple<>>
{
public:
 static constexpr std::size_t value = 1;
};

template <>
class TupleTypeSize<std::tuple<>&>
{
public:
 static constexpr std::size_t value = 1;
};

template <>
class TupleTypeSize<const std::tuple<>&>
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

template <typename T,typename ... Types>
class TupleTypeSize<const std::tuple<T,Types...>>
{
public:
 static constexpr std::size_t value = TupleTypeSize<std::tuple<T,Types...>>::value;
};

template <typename T,typename ... Types>
class TupleTypeSize<const std::tuple<T,Types...>&>
{
public:
 static constexpr std::size_t value = TupleTypeSize<std::tuple<T,Types...>>::value;
};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// TypeToTupleElementPosition<T,Tuple>::value==position of T in Tuple. If not present, return -1
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template <class T,Integer Nmax, class Tuple>
struct TypeToTupleElementPositionHelper;

template <class T,Integer Nmax>
struct TypeToTupleElementPositionHelper<T, Nmax,std::tuple<>> {
    static const Integer value = -Nmax-1;
};

template <class T,Integer Nmax, class... Types>
struct TypeToTupleElementPositionHelper<T, Nmax, std::tuple<T, Types...>> {
    static const Integer value = 0;
};

template <class T,Integer Nmax, class U, class... Types>
struct TypeToTupleElementPositionHelper<T, Nmax, std::tuple<U, Types...>> {
    static const Integer value = 1 + TypeToTupleElementPositionHelper<T, Nmax, std::tuple<Types...>>::value;
};

template <class T, class Tuple>
using TypeToTupleElementPosition=TypeToTupleElementPositionHelper<T,TupleTypeSize<Tuple>::value,Tuple>;






template<typename T, typename TupleOfTuple,Integer M,Integer Nmax,Integer N>
constexpr std::enable_if_t<(N==Nmax+1),Integer > NthTupleOfTupleTypePositionHelper()
{return -1;};

template<typename T, typename TupleOfTuple,Integer M,Integer Nmax,Integer N>
constexpr std::enable_if_t<(N<Nmax+1),Integer > NthTupleOfTupleTypePositionHelper()
{
  if(IsSame<T,GetType<TupleOfTuple,N,M> >::value)
    return N;
  else 
    return NthTupleOfTupleTypePositionHelper<T,TupleOfTuple,M,Nmax,N+1>();
};

template<typename T, typename TupleOfTuple,Integer M>
constexpr Integer NthTupleOfTupleTypePosition()
{
  return NthTupleOfTupleTypePositionHelper<T,TupleOfTuple,M,TupleTypeSize<TupleOfTuple>::value-1,0>();
}







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

// template<Integer N,Integer Nmin,Integer Nmax>
// class SubTupleHelper<N,Nmin,Nmax,std::tuple<>>
// {
// public:
//  using type=std::tuple<>;
// };

// template<Integer N,Integer Nmin,Integer Nmax,typename T,typename...Ts>
// class SubTupleHelper<N,Nmin,Nmax,std::tuple<T,Ts...>>
// {
// public:
//  using single_type=typename std::conditional< (N<Nmin || N>Nmax), std::tuple<>, std::tuple<T> >::type;
//  using type=TupleCatType< single_type, typename SubTupleHelper<N+1,Nmin,Nmax,std::tuple<Ts...> >::type >;
// };

// template<Integer Nmin,Integer Nmax,typename...Ts>
// using SubTupleType=typename SubTupleHelper<0,Nmin,Nmax,Ts... >::type;



template<Integer N,Integer Nmin,Integer Nmax>
class SubTupleHelper<N,Nmin,Nmax,std::tuple<>>
{
public:
 using type=std::tuple<>;
};


template<Integer Nmin,Integer Nmax,typename...Ts>
class SubTupleHelper<sizeof...(Ts)-1,Nmin,Nmax,std::tuple<Ts...>>
{
public:
 // using single_type=std::conditional_t< (N<Nmin || N>Nmax), std::tuple<>, std::tuple<GetType<std::tuple<Ts...>,N>>> ;
 static constexpr Integer N=sizeof...(Ts)-1;
 using type=std::conditional_t< (N<Nmin || N>Nmax),
                                 std::tuple<>,
                                 std::tuple<GetType<std::tuple<Ts...>,N>>>;
};

template<Integer N,Integer Nmin,Integer Nmax,typename...Ts>
class SubTupleHelper<N,Nmin,Nmax,std::tuple<Ts...>>
{
public:
 // using single_type=std::conditional_t< (N<Nmin || N>Nmax), std::tuple<>, std::tuple<GetType<std::tuple<Ts...>,N>>> ;
// if N<Nmin, we just move on, without adding types
// if N>Nmax, we stop, returning an empty tuple
// otherwise add the component
using type=
std::conditional_t<(N<Nmin),
                   typename SubTupleHelper<N+1,Nmin,Nmax,std::tuple<Ts...> >::type,
                   std::conditional_t< (N>Nmax),
                                       std::tuple<>,
                                       TupleCatType<std::tuple<GetType<std::tuple<Ts...>,N>>, 
                                                    typename SubTupleHelper<N+1,Nmin,Nmax,std::tuple<Ts...> >::type>
                                     >

                  >;


 // using type=std::conditional_t< (N<Nmin || N>Nmax),
 //                                 std::tuple<>,
 //                                 TupleCatType< std::tuple<GetType<std::tuple<Ts...>,N>>, 
 //                                               typename SubTupleHelper<N+1,Nmin,Nmax,std::tuple<Ts...> >::type >>;
};



template<Integer Nmin,Integer Nmax,typename...Ts>
using SubTupleType=typename SubTupleHelper<0,Nmin,Nmax,std::tuple<Ts...> >::type;


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
///// TupleRemoveNthPosition<N,Tuple> remove the n-th element of the tuple from Tuple
///// If N>Size(Tuple) or N<0, nothing is removed
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<Integer N, typename...Ts>
using TupleRemoveNthPosition=TupleCatType<SubTupleType<0,N-1,Ts...>,SubTupleType<N+1,TupleTypeSize<std::tuple<Ts...>>::value-1,Ts...>>;


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
 using typeRight=SubTupleType<N+1,TupleTypeSize<std::tuple<Ts...>>::value-1, Ts...>;//typename std::conditional< (N==sizeof...(Ts)+1), std::tuple<>, SubTupleType<N,sizeof...(Ts)+1, Ts...> >::type;
 // using type=Number<TupleTypeSize<std::tuple<Ts...>>::value-1>;//
 using type=TupleCatType<typeLeft,std::tuple<Tchange>,typeRight>;
};

template<Integer N,typename Tchange, typename...Ts>
using TupleChangeType=typename TupleChangeTypeHelper< N,Tchange,Ts... >::type;


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
class TupleChangeNthTypeHelper;

template<Integer Nmax,typename Tuple,Integer Nchange,typename Tchange>
class TupleChangeNthTypeHelper<Nmax,Nmax,Tuple,Nchange,Tchange>
{
public:
using type=std::tuple<TupleChangeType<Nchange,Tchange,GetType<Tuple,Nmax>>>;
};

template<Integer N,Integer Nmax,typename Tuple,Integer Nchange,typename Tchange>
class TupleChangeNthTypeHelper
{
 public:
 using single_type=std::tuple<TupleChangeType<Nchange,Tchange,GetType<Tuple,N>>>;

 using type=decltype(std::tuple_cat(std::declval<single_type>(),
                           std::declval<typename TupleChangeNthTypeHelper<N+1,Nmax,Tuple,Nchange,Tchange>::type>()));


};

template<typename Tuple,Integer Nchange,typename Tchange>
using TupleChangeNthType=typename TupleChangeNthTypeHelper<0,TupleTypeSize<Tuple>::value-1,Tuple,Nchange,Tchange>::type;



template<Integer N,Integer Nmax,Integer Nchange,typename Tchange, typename TupleOfTuple>
class TupleOfTupleChangeTypeHelper;


template<Integer Nmax,Integer Nchange,typename Tchange, typename TupleOfTuple>
class TupleOfTupleChangeTypeHelper<Nmax,Nmax,Nchange,Tchange,TupleOfTuple>
{
 public:
 using T=GetType<TupleOfTuple,Nmax>;
 using type=typename std::conditional< std::is_same< T,std::tuple<> >::value , std::tuple<T>, std::tuple<TupleChangeNthType<T,Nchange,Tchange >>>::type;
};



template<Integer N,Integer Nmax,Integer Nchange,typename Tchange, typename TupleOfTuple>
class TupleOfTupleChangeTypeHelper
{
public:
using T=GetType<TupleOfTuple,N>;
using single_type=typename std::conditional< std::is_same< T,std::tuple<> >::value , std::tuple<T>, std::tuple<TupleChangeNthType<T,Nchange,Tchange >>>::type;
using type=decltype(std::tuple_cat(std::declval<single_type>(),
                           std::declval<typename TupleOfTupleChangeTypeHelper<N+1,Nmax,Nchange,Tchange,TupleOfTuple>::type>()));
};


template<Integer Nchange,typename Tchange,typename TupleOfTuple>
using TupleOfTupleChangeType=typename TupleOfTupleChangeTypeHelper<0,TupleTypeSize<TupleOfTuple>::value-1,Nchange,Tchange,TupleOfTuple>::type;









template<Integer N,Integer Nmax, typename Tuple>
class Tuple2ToTuple1Helper;

template<Integer Nmax>
class Tuple2ToTuple1Helper<Nmax,Nmax,std::tuple<>>
{
 public:
 using type=std::tuple<> ;
};


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

template<Integer N,Integer Nmax>
class Tuple2ToTuple1Helper<N,Nmax,std::tuple<>>
{
 public:
 using type=std::tuple<> ;
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

















template<typename TupleOfSpaces, Integer Nmax,Integer N>
class UniqueFEFamiliesHelper2;

template<typename TupleOfSpaces, Integer Nmax>
class UniqueFEFamiliesHelper2<TupleOfSpaces,Nmax,Nmax>
{
 public:
  using Space=GetType<TupleOfSpaces,Nmax>;
  using type=std::tuple< Number<Space::FEFamily> >;
};

template<typename TupleOfSpaces, Integer Nmax,Integer N=0>
class UniqueFEFamiliesHelper2
{
 public:
  using Space=GetType<TupleOfSpaces,N>;
  using single_type=std::tuple< Number<Space::FEFamily> >;
  using type=decltype(std::tuple_cat(std::declval<single_type>(),
                             std::declval<typename UniqueFEFamiliesHelper2<TupleOfSpaces,Nmax,N+1>::type>()));
};


template<typename TupleOfSpaces>
using UniqueFEFamilies2=RemoveTupleDuplicates<typename UniqueFEFamiliesHelper2<TupleOfSpaces,TupleTypeSize<TupleOfSpaces>::value-1,0>::type>;








template<typename TupleOfSpaces,Integer Nmax,Integer N>
class SpacesToUniqueFEFamiliesHelper2;

template<typename TupleOfSpaces,Integer Nmax>
class SpacesToUniqueFEFamiliesHelper2<TupleOfSpaces,Nmax,Nmax>
{
 public:
  using TupleFEFamilies=UniqueFEFamilies2<TupleOfSpaces>;
  using FamilyNumber=Number<GetType<TupleOfSpaces,Nmax>::FEFamily>;
  using PositionNumber=Number<TypeToTupleElementPosition<FamilyNumber,TupleFEFamilies>::value>;
  using type=std::tuple< PositionNumber >;
};

template<typename TupleOfSpaces,Integer Nmax,Integer N=0>
class SpacesToUniqueFEFamiliesHelper2
{
 public:
  using TupleFEFamilies=UniqueFEFamilies2<TupleOfSpaces>;
  using FamilyNumber=Number<GetType<TupleOfSpaces,N>::FEFamily>;
  using PositionNumber=Number<TypeToTupleElementPosition<FamilyNumber,TupleFEFamilies>::value>;
  using single_type=std::tuple< PositionNumber  >;
  using type=decltype(std::tuple_cat(std::declval<single_type>(),
                             std::declval<typename SpacesToUniqueFEFamiliesHelper2<TupleOfSpaces,Nmax,N+1>::type>()));
};



template<typename TupleOfSpaces>
using SpacesToUniqueFEFamilies2=typename SpacesToUniqueFEFamiliesHelper2<TupleOfSpaces,TupleTypeSize<TupleOfSpaces>::value-1,0>::type;












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







// template <typename Tuple,typename Type,Integer Nmax,Integer N>
// class HowMayTypesDifferentFromTypeHelper;


// template <typename Tuple,typename Type,Integer Nmax>
// class HowMayTypesDifferentFromTypeHelper<Tuple,Type,Nmax,Nmax>
// {
// public:
//   static constexpr Integer value=IsDifferent<Type,GetType<Tuple,Nmax>>::value;
// };

// template <typename Tuple,typename Type,Integer Nmax,Integer N>
// class HowMayTypesDifferentFromTypeHelper
// {
// public:
//   static constexpr Integer value=IsDifferent<Type,GetType<Tuple,N>>::value+HowMayTypesDifferentFromTypeHelper<Tuple,Type,Nmax,N+1>::value;
// };

// template <typename Tuple,typename Type>
// static constexpr Integer HowMayTypesDifferentFromType=HowMayTypesDifferentFromTypeHelper<Tuple,Type,TupleTypeSize<Tuple>::value-1,0>::value;





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
class MapFromReference;

template<typename Tuple, Integer FEFamily, Integer Nmax,Integer N>
class MapOperatorTupleHelper;

template<typename Tuple, Integer FEFamily, Integer Nmax>
class MapOperatorTupleHelper<Tuple,FEFamily,Nmax,Nmax>
{
 public:

  using Nthelem=GetType<Tuple,Nmax>;
  using Operator=GetType<Nthelem,0>;
  using Elem=GetType<Nthelem,1>;
  using type=std::tuple< MapFromReference<Operator,Elem,FEFamily> >;
};

template<typename Tuple, Integer FEFamily,Integer Nmax,Integer N=0>
class MapOperatorTupleHelper
{
 public:
  using Nthelem=GetType<Tuple,N>;
  using Operator=GetType<Nthelem,0>;
  using Elem=GetType<Nthelem,1>;
  using single_type=std::tuple< MapFromReference<Operator,Elem,FEFamily> >;
  using type=decltype(std::tuple_cat(std::declval<single_type>(),
                             std::declval<typename MapOperatorTupleHelper<Tuple,FEFamily,Nmax,N+1>::type>()));
};


template<typename Tuple, Integer FEFamily,Integer Nmax>
using MapOperatorTuple=typename MapOperatorTupleHelper<Tuple,FEFamily,Nmax,0>::type;


template<typename TupleOfTuple, typename TupleSpaces, Integer Nmax,Integer N>
class MapOperatorTupleOfTupleHelper;


template<typename Tuple,typename TupleSpaces, Integer N>
class MapOperatorTupleOfTupleHelperAux;

template<typename TupleSpaces, Integer N>
class MapOperatorTupleOfTupleHelperAux<std::tuple<>,TupleSpaces,N>
{
public:
  using type=std::tuple<std::tuple<>>;
};

template<typename Tuple,typename TupleSpaces, Integer N>
class MapOperatorTupleOfTupleHelperAux
{
public:
static constexpr Integer FEFamily=GetType<TupleSpaces,N>::FEFamily;
static constexpr Integer Nmax=TupleTypeSize<Tuple>::value-1;
using type=std::tuple< MapOperatorTuple< Tuple,FEFamily,Nmax> > ;
};

template<typename TupleOfTuple,typename TupleSpaces, Integer Nspaces>
class MapOperatorTupleOfTupleHelper<TupleOfTuple,TupleSpaces, Nspaces, Nspaces>
{
public:
using Tuple=GetType<TupleOfTuple,Nspaces>;
using type=typename MapOperatorTupleOfTupleHelperAux<Tuple,TupleSpaces,Nspaces>::type;
// static constexpr Integer FEFamily=GetType<TupleSpaces,Nspaces>::FEFamily;
// static constexpr Integer Nmax=TupleTypeSize<Tuple>::value-1;
// using type=std::tuple< MapOperatorTuple< Tuple,FEFamily,Nmax> > ;
};

template<typename TupleOfTuple,typename TupleSpaces, Integer Nspaces, Integer N=0>
class MapOperatorTupleOfTupleHelper
{
public:
using Tuple=GetType<TupleOfTuple,N>;
using single_type=typename MapOperatorTupleOfTupleHelperAux<Tuple,TupleSpaces,N>::type;
// static constexpr Integer FEFamily=GetType<TupleSpaces,N>::FEFamily;
// static constexpr Integer Nmax=TupleTypeSize<Tuple>::value-1;
// using single_type=std::tuple< MapOperatorTuple< Tuple, FEFamily, Nmax> > ;
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
class ShapeFunction;

template<typename Elem, typename BaseFunctioSpace, typename Operator, typename QuadratureRule,Integer M>
class ShapeFunctionCombinations;

// template<typename Elem, typename BaseFunctioSpace, typename Tuple, Integer Nmax,Integer N>
// class TupleShapeFunctionCreate;

// template<typename Elem, typename BaseFunctioSpace, typename Tuple, Integer Nmax>
// class TupleShapeFunctionCreate<Elem,BaseFunctioSpace,Tuple,Nmax,Nmax>
// {
//  public:

//   using Nelem=GetType<Tuple,Nmax>;
//   using Operator=GetType<Nelem,0>;
//   using QuadratureRule=GetType<Nelem,1>;
//   using type=std::tuple< ShapeFunction<Elem, BaseFunctioSpace,Operator,QuadratureRule > >;
// };



// template<typename Elem, typename BaseFunctioSpace, typename Tuple, Integer Nmax,Integer N=0>
// class TupleShapeFunctionCreate
// {
//  public:
//   using Nelem=GetType<Tuple,N>;
//   using Operator=GetType<Nelem,0>;
//   using QuadratureRule=GetType<Nelem,1>;
//   using single_type=std::tuple< ShapeFunction<Elem, BaseFunctioSpace,Operator,QuadratureRule > >;
//   using type=decltype(std::tuple_cat(std::declval<single_type>(),
//                              std::declval<typename TupleShapeFunctionCreate<Elem,BaseFunctioSpace,Tuple,Nmax,N+1>::type>()));
// };


// template<typename TupleSpaces, typename TupleOperatorsAndQuadrature, Integer Nmax,Integer N>
// class TupleOfTupleShapeFunctionCreate;



// template<typename TupleSpaces,typename TupleOperatorsAndQuadrature, Integer Nspaces>
// class TupleOfTupleShapeFunctionCreate<TupleSpaces,TupleOperatorsAndQuadrature, Nspaces, Nspaces>
// {
// public:
// using OperatorAndQuadrature=GetType<TupleOperatorsAndQuadrature,Nspaces>;
// static constexpr Integer Nmax=TupleTypeSize<OperatorAndQuadrature>::value-1;
// using Space=GetType<TupleSpaces,Nspaces>;
// using Elem=GetType<Space,0>;
// using FunctionSpace=GetType<Space,1>;
// using type=std::tuple<typename TupleShapeFunctionCreate< Elem,FunctionSpace, OperatorAndQuadrature,Nmax>::type>;
// };

// template<typename TupleSpaces,typename TupleOperatorsAndQuadrature, Integer Nspaces, Integer N=0>
// class TupleOfTupleShapeFunctionCreate
// {
// public:
// using OperatorAndQuadrature=GetType<TupleOperatorsAndQuadrature,N>;
// static constexpr Integer Nmax=TupleTypeSize<OperatorAndQuadrature>::value-1;
// using Space=GetType<TupleSpaces,N>;
// using Elem=GetType<Space,0>;
// using FunctionSpace=GetType<Space,1>;
// using single_type=std::tuple<typename TupleShapeFunctionCreate< Elem,FunctionSpace, OperatorAndQuadrature,Nmax>::type>;
// using type=decltype(std::tuple_cat(std::declval<single_type>(),
//                            std::declval<typename TupleOfTupleShapeFunctionCreate<TupleSpaces,TupleOperatorsAndQuadrature,Nspaces,N+1>::type>()));
// };



// template<typename TupleSpaces,typename TupleOperatorsAndQuadrature>
// using TupleOfTupleShapeFunctionType=typename TupleOfTupleShapeFunctionCreate<TupleSpaces,TupleOperatorsAndQuadrature,TupleTypeSize<TupleSpaces>::value-1,0>::type;
























template<typename Elem, typename BaseFunctioSpace, typename Tuple, Integer Nmax,Integer N>
class TupleShapeFunctionCreate2;

template<typename Elem, typename BaseFunctioSpace, typename Tuple, Integer Nmax>
class TupleShapeFunctionCreate2<Elem,BaseFunctioSpace,Tuple,Nmax,Nmax>
{
 public:

  using Nelem=GetType<Tuple,Nmax>;
  using Operator=GetType<Nelem,0>;
  using QuadratureRule=GetType<Nelem,1>;
  using type=std::tuple< ShapeFunction<Elem, BaseFunctioSpace,Operator,QuadratureRule > >;
};



template<typename Elem, typename BaseFunctioSpace, typename Tuple, Integer Nmax,Integer N=0>
class TupleShapeFunctionCreate2
{
 public:
  using Nelem=GetType<Tuple,N>;
  using Operator=GetType<Nelem,0>;
  using QuadratureRule=GetType<Nelem,1>;
  using single_type=std::tuple< ShapeFunction<Elem, BaseFunctioSpace,Operator,QuadratureRule > >;
  using type=decltype(std::tuple_cat(std::declval<single_type>(),
                             std::declval<typename TupleShapeFunctionCreate2<Elem,BaseFunctioSpace,Tuple,Nmax,N+1>::type>()));
};

template<typename TupleSpaces,typename OperatorAndQuadrature, Integer N>
class TupleOfTupleShapeFunctionCreate2Aux;


template<typename TupleSpaces, Integer N>
class TupleOfTupleShapeFunctionCreate2Aux<TupleSpaces,std::tuple<>,N>
{
 public:
  using type=std::tuple<std::tuple<>>;
};

template<typename TupleSpaces,typename OperatorAndQuadrature, Integer N>
class TupleOfTupleShapeFunctionCreate2Aux
{
public:
  static constexpr Integer Nmax=TupleTypeSize<OperatorAndQuadrature>::value-1;
  using Space=GetType<TupleSpaces,N>;
  using Elem=typename Space::Elem;
  using FunctionSpace=Elem2FunctionSpace<Space>;
  using type=std::tuple<typename TupleShapeFunctionCreate2< Elem,FunctionSpace, OperatorAndQuadrature,Nmax>::type>;
};


template<typename TupleSpaces, typename TupleOperatorsAndQuadrature, Integer Nmax,Integer N>
class TupleOfTupleShapeFunctionCreate2;



template<typename TupleSpaces,typename TupleOperatorsAndQuadrature, Integer Nspaces>
class TupleOfTupleShapeFunctionCreate2<TupleSpaces,TupleOperatorsAndQuadrature, Nspaces, Nspaces>
{
public:
using OperatorAndQuadrature=GetType<TupleOperatorsAndQuadrature,Nspaces>;
using type=typename TupleOfTupleShapeFunctionCreate2Aux<TupleSpaces,OperatorAndQuadrature,Nspaces>::type;
// static constexpr Integer Nmax=TupleTypeSize<OperatorAndQuadrature>::value-1;
// using Space=GetType<TupleSpaces,Nspaces>;
// using Elem=typename Space::Elem;
// using FunctionSpace=Elem2FunctionSpace<Space>;
// using type=std::tuple<typename TupleShapeFunctionCreate2< Elem,FunctionSpace, OperatorAndQuadrature,Nmax>::type>;
};

template<typename TupleSpaces,typename TupleOperatorsAndQuadrature, Integer Nspaces, Integer N=0>
class TupleOfTupleShapeFunctionCreate2
{
public:
using OperatorAndQuadrature=GetType<TupleOperatorsAndQuadrature,N>;
using single_type=typename TupleOfTupleShapeFunctionCreate2Aux<TupleSpaces,OperatorAndQuadrature,N>::type;

// static constexpr Integer Nmax=TupleTypeSize<OperatorAndQuadrature>::value-1;
// using Space=GetType<TupleSpaces,N>;
// using Elem=typename Space::Elem;
// using FunctionSpace=Elem2FunctionSpace<Space>;
// using single_type=std::tuple<typename TupleShapeFunctionCreate2< Elem,FunctionSpace, OperatorAndQuadrature,Nmax>::type>;
using type=decltype(std::tuple_cat(std::declval<single_type>(),
                           std::declval<typename TupleOfTupleShapeFunctionCreate2<TupleSpaces,TupleOperatorsAndQuadrature,Nspaces,N+1>::type>()));
};



template<typename TupleSpaces,typename TupleOperatorsAndQuadrature>
using TupleOfTupleShapeFunctionType2=typename TupleOfTupleShapeFunctionCreate2<TupleSpaces,TupleOperatorsAndQuadrature,TupleTypeSize<TupleSpaces>::value-1,0>::type;





























template<typename Elem, typename BaseFunctioSpace, typename Tuple, Integer Nmax,Integer N,Integer M>
class TupleShapeFunctionCreate3;

template<typename Elem, typename BaseFunctioSpace, typename Tuple, Integer Nmax,Integer M>
class TupleShapeFunctionCreate3<Elem,BaseFunctioSpace,Tuple,Nmax,Nmax,M>
{
 public:

  using Nelem=GetType<Tuple,Nmax>;
  using Operator=GetType<Nelem,0>;
  using QuadratureRule=GetType<Nelem,1>;
  using type=std::tuple<ShapeFunctionCombinations<Elem, BaseFunctioSpace,Operator,QuadratureRule,M> >;
};



template<typename Elem, typename BaseFunctioSpace, typename Tuple, Integer Nmax,Integer N,Integer M>
class TupleShapeFunctionCreate3
{
 public:
  using Nelem=GetType<Tuple,N>;
  using Operator=GetType<Nelem,0>;
  using QuadratureRule=GetType<Nelem,1>;
  using single_type=std::tuple<ShapeFunctionCombinations<Elem, BaseFunctioSpace,Operator,QuadratureRule,M>>;
  using type=decltype(std::tuple_cat(std::declval<single_type>(),
                             std::declval<typename TupleShapeFunctionCreate3<Elem,BaseFunctioSpace,Tuple,Nmax,N+1,M>::type>()));
};

template<Integer N,Integer Nmax,typename...Ts>
class VoidTupleOrShapeFunctionCombination;


template<Integer N,Integer Nmax,typename Space>
class VoidTupleOrShapeFunctionCombination<N,Nmax,Space,std::tuple<>>
{
 public:
  using type=std::tuple<std::tuple<>>;
};

template<Integer N,Integer Nmax,typename Space,typename OperatorAndQuadrature>
class VoidTupleOrShapeFunctionCombination<N,Nmax,Space,OperatorAndQuadrature>
{
 public:
using Elem=typename Space::Elem;
using FunctionSpace=Elem2FunctionSpace<Space>;
using type=std::tuple<typename TupleShapeFunctionCreate3< Elem,FunctionSpace, OperatorAndQuadrature,Nmax,0,N>::type>;
};


template<typename TupleSpaces, typename TupleOperatorsAndQuadrature, Integer Nmax,Integer N>
class TupleOfTupleShapeFunctionCreate3;

template<typename TupleSpaces,typename TupleOperatorsAndQuadrature, Integer Nspaces>
class TupleOfTupleShapeFunctionCreate3<TupleSpaces,TupleOperatorsAndQuadrature, Nspaces, Nspaces>
{
public:
using OperatorAndQuadrature=GetType<TupleOperatorsAndQuadrature,Nspaces>;
static constexpr Integer Nmax=TupleTypeSize<OperatorAndQuadrature>::value-1;
using Space=GetType<TupleSpaces,Nspaces>;
using type=typename VoidTupleOrShapeFunctionCombination<Nspaces,Nmax,Space,OperatorAndQuadrature>::type; 
};

template<typename TupleSpaces,typename TupleOperatorsAndQuadrature, Integer Nspaces, Integer N=0>
class TupleOfTupleShapeFunctionCreate3
{
public:
using OperatorAndQuadrature=GetType<TupleOperatorsAndQuadrature,N>;
static constexpr Integer Nmax=TupleTypeSize<OperatorAndQuadrature>::value-1;
using Space=GetType<TupleSpaces,N>;
using single_type=typename VoidTupleOrShapeFunctionCombination<N,Nmax,Space,OperatorAndQuadrature>::type; 

using type=decltype(std::tuple_cat(std::declval<single_type>(),
                           std::declval<typename TupleOfTupleShapeFunctionCreate3<TupleSpaces,TupleOperatorsAndQuadrature,Nspaces,N+1>::type>()));
};



template<typename TupleSpaces,typename TupleOperatorsAndQuadrature>
using TupleOfTupleShapeFunctionType3=typename TupleOfTupleShapeFunctionCreate3<TupleSpaces,TupleOperatorsAndQuadrature,TupleTypeSize<TupleSpaces>::value-1,0>::type;

































template<typename...Ts>
class TupleRemoveNumber0Helper;

template<typename...Ts>
class TupleRemoveNumber0Helper<std::tuple<Ts...>> 
{
public:
using Tuple=std::tuple<Ts...>;
static constexpr Integer size=TupleTypeSize<Tuple>::value;
using type= GetType<typename std::conditional<size==1,
                                             Tuple, 
                                             TupleRemoveType<Number<0>,Tuple>
                                            >::type>;
};


template<typename...Ts>
using TupleRemoveNumber0=typename TupleRemoveNumber0Helper<std::tuple<Ts...>>::type;













































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

template<typename Map, Integer Nmax>
class MapTupleInit2Helper<std::tuple<>,Map, Nmax, Nmax>
{
public:
      using type=std::tuple<std::tuple<>>;
};

template<typename Map, Integer Nmax,Integer N>
class MapTupleInit2Helper<std::tuple<>,Map, Nmax, N>
{
public:
      using type=std::tuple<std::tuple<>>;
};

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









  
  template<typename GeneralForm,typename...GeneralForms>
  class ShapeFunctionCoefficientsCollection;



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
  shape_function_init_aux(Tuple& tuple, const ShapeFunctionCoefficientsCollection<Args...> &coefficients)
  {}

  template<Integer Nmax,Integer N, typename Tuple,typename...Args>
  typename std::enable_if_t<N<=Nmax,void> 
  shape_function_init_aux(Tuple& tuple, const ShapeFunctionCoefficientsCollection<Args...> &coefficients)
  {
    constexpr Integer Nmax_aux=TupleTypeSize<GetType<Tuple,N>>::value-1;
    // shape_function_init_aux_aux<Nmax_aux,0>(tuple_get<N>(tuple),coefficients);
    shape_function_init_aux<Nmax,N+1>(tuple,coefficients);
  }



  template< typename Tuple,typename...Args>
  constexpr void init(Tuple& tuple, const ShapeFunctionCoefficientsCollection<Args...> &coefficients)
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







template<typename T,Integer Nmax,Integer N>
class UniqueElementFEFamilyHelper2;

template<typename T,Integer Nmax>
class UniqueElementFEFamilyHelper2<T,Nmax,Nmax>
{
public:
  using TN=GetType<T,Nmax>;
  using type=std::tuple<std::tuple<typename TN::Elem,Number<TN::FEFamily>>>;
};

template<typename T,Integer Nmax,Integer N>
class UniqueElementFEFamilyHelper2
{
 public:
  using TN=GetType<T,N>;
  using single_type=std::tuple<std::tuple<typename TN::Elem,Number<TN::FEFamily>>>;
  using type=decltype(std::tuple_cat(std::declval<single_type>(),
                      std::declval<typename UniqueElementFEFamilyHelper2<T,Nmax,N+1>::type>()));
};
template<typename T>
using UniqueElementFEFamily2=RemoveTupleDuplicates<typename UniqueElementFEFamilyHelper2<T,TupleTypeSize<T>::value-1,0>::type>;
















template<typename ...Ts> 
class BubbleSortTupleOfPairsNumbersAuxHelper;

template<typename ...Ts> 
using BubbleSortTupleOfPairsNumbersAux=typename BubbleSortTupleOfPairsNumbersAuxHelper<Ts...>::type;

// linear form, main case (only one given number)
template<Integer N>
class BubbleSortTupleOfPairsNumbersAuxHelper<std::tuple<std::tuple<Number<N>>> >
{
    public:
    using type=std::tuple<std::tuple<Number<N>>>;
};

// bilinear form, main case (only a pairs of given numbers)
template<Integer M1, Integer M2>
class BubbleSortTupleOfPairsNumbersAuxHelper<std::tuple<Number<M1>,Number<M2>>>
{
    public:
    using type=typename std::tuple<std::tuple<Number<M1>,Number<M2>>>;
};

// linear form, main reordering case
template<Integer M, Integer N>
class BubbleSortTupleOfPairsNumbersAuxHelper<std::tuple<Number<M>>, 
                                             std::tuple<Number<N>> >
{
    public:
    using T1=std::tuple<Number<M>>;
    using T2=std::tuple<Number<N>>;
    using type=typename std::conditional<(M<N),
                                         std::tuple<T1,T2> ,
                                         std::tuple<T2,T1>>::type;
};

// template<Integer...Ns>
// class BubbleSortTupleOfPairsNumbersAuxHelper<std::tuple<>,
//                                              std::tuple< std::tuple<Number<Ns>>...>>
// {
//     public:
//     using type=BubbleSortTupleOfPairsNumbersAux<std::tuple< std::tuple<Number<Ns>>...>>;
// };

// template<Integer...Ns>
// class BubbleSortTupleOfPairsNumbersAuxHelper<std::tuple<>, 
//                                              std::tuple<std::tuple<Number<Ns>>...> >
// {
//     public:
//     using type=BubbleSortTupleOfPairsNumbersAuxHelper<std::tuple<Number<Ns>>...>;
// };

// template<Integer...Ns>
// class BubbleSortTupleOfPairsNumbersAuxHelper<std::tuple<std::tuple<Number<Ns>>...>,
//                                              std::tuple<> >
// {
//     public:
//     using type=BubbleSortTupleOfPairsNumbersAuxHelper<std::tuple<Number<Ns>>...>;
// };



template<Integer M1, Integer M2, Integer N1, Integer N2>
class BubbleSortTupleOfPairsNumbersAuxHelper<std::tuple<Number<M1>,Number<M2>>, 
                                          std::tuple<Number<N1>,Number<N2>> >
{
    public:
    using T1=std::tuple<Number<M1>,Number<M2>>;
    using T2=std::tuple<Number<N1>,Number<N2>>;
    using type=typename std::conditional<(M1<N1)||(M1==N1 && M2<N2),
                                         std::tuple<T1,T2> ,
                                         std::tuple<T2,T1>>::type;
};

template<Integer M1, Integer M2>
class BubbleSortTupleOfPairsNumbersAuxHelper<std::tuple< std::tuple<Number<M1>,Number<M2>>>>
{
    public:

    using type=std::tuple< std::tuple<Number<M1>,Number<M2>>>;
};

template<Integer M1, Integer M2, Integer N1, Integer N2>
class BubbleSortTupleOfPairsNumbersAuxHelper<std::tuple<std::tuple<Number<M1>,Number<M2>>, 
                                                        std::tuple<Number<N1>,Number<N2>>> >
{
    public:
    using T1=std::tuple<Number<M1>,Number<M2>>;
    using T2=std::tuple<Number<N1>,Number<N2>>;
    using type=BubbleSortTupleOfPairsNumbersAux<T1,T2>;
};

template<typename S, typename T,typename ...Ss> 
class BubbleSortTupleOfPairsNumbersAuxHelper<std::tuple<Ss...>, std::tuple<S, T>>
{
 
 public:
  using switch_type=BubbleSortTupleOfPairsNumbersAux<T,S>;
  using type=TupleCatType<std::tuple<Ss...>,switch_type>;

};

template<typename S, typename T,typename ...Ts,typename ...Ss> 
class BubbleSortTupleOfPairsNumbersAuxHelper<std::tuple<Ss...>, std::tuple<S, T, Ts ...>>
{
 
 public:
  using switch_type=BubbleSortTupleOfPairsNumbersAux<T,S>;
  using tmp_type=TupleCatType<switch_type,std::tuple<Ts...>>;
  using T1=TupleCatType<std::tuple<Ss...>,std::tuple<GetType<tmp_type,0>>>;
  using T2=TupleCatType<std::tuple<GetType<tmp_type,1>>, std::tuple<Ts...>>;
  using type=BubbleSortTupleOfPairsNumbersAux<T1,T2> ;
};


// template<typename S,typename T,typename ...Ts> 
// class BubbleSortTupleOfPairsNumbersAuxHelper<std::tuple<>, std::tuple<S, T, Ts ...>>
// {
 
//  public:
//   using switch_type=BubbleSortTupleOfPairsNumbersAux<S,T>;
//   using tmp_type=TupleCatType<switch_type,std::tuple<Ts...>>;
//   using T1=std::tuple<GetType<tmp_type,0>>;
//   using T2=TupleCatType<std::tuple<GetType<tmp_type,1>>,std::tuple<Ts...>>;
//   using type=BubbleSortTupleOfPairsNumbersAux<T1,T2> ;
// };

template<typename R,typename S> 
class BubbleSortTupleOfPairsNumbersAuxHelper<std::tuple<R, S>>
{
 
 public:
  using switch_type=BubbleSortTupleOfPairsNumbersAux<R,S>;
  using type=TupleCatType<switch_type>;
};


template<typename R,typename S,typename ...Ts> 
class BubbleSortTupleOfPairsNumbersAuxHelper<std::tuple<R, S, Ts ...>>
{
 
 public:
  using switch_type=BubbleSortTupleOfPairsNumbersAux<R,S>;
  using tmp_type=TupleCatType<switch_type,std::tuple<Ts...>>;
  using T1=std::tuple<GetType<tmp_type,0>>;
  using T2=TupleCatType<std::tuple<GetType<tmp_type,1>>,std::tuple<Ts...>>;
  using type=BubbleSortTupleOfPairsNumbersAux<T1,T2> ;
};

///////////////////////////////////////////////////////////////////////////////////////////////
///// This is the first function which is called
///// We insert a void tuple 
///////////////////////////////////////////////////////////////////////////////////////////////
// template<typename ...Ts> 
// class BubbleSortTupleOfPairsNumbersAuxHelper<std::tuple<Ts...>>
// {
 
//  public:
//   using type=BubbleSortTupleOfPairsNumbersAux< std::tuple<>, std::tuple<Ts...> >;
  
// };

template<typename Tuple, Integer Nmax, Integer N>
class BubbleSortTupleOfPairsNumbersHelper;




///////////////////////////////////////////////////////////////////////////////////////////////////////////
///// We have to reorder the X=tuple<x0,x1,x2,x3...>. We 
///// We Start by  checking each pair:
///// Y_1=(sort(x0,x1),x2,x3...)=(y0,y1,y2...)
///// Z_1=(y0,sort(y1,y2),y3...)=(z0,z1,z2...)
///// Until we get to the end. Then, we only know that if x_0 was the largest, now it is at the end.
///// But the other elements still have to be reordered. 
///// So we call BubbleSortTupleOfPairsNumbersHelper Nmax times, where Nmax is the size of the tuple.
///////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Tuple>
using BubbleSortTupleOfPairsNumbers=typename BubbleSortTupleOfPairsNumbersHelper<Tuple,TupleTypeSize<Tuple>::value-1,0>::type;


template<typename Tuple, Integer Nmax>
class BubbleSortTupleOfPairsNumbersHelper<Tuple,Nmax,Nmax>
{
 public:
  using type=BubbleSortTupleOfPairsNumbersAux<Tuple>;
};

template<typename Tuple, Integer Nmax, Integer N>
class BubbleSortTupleOfPairsNumbersHelper
{
 public:
  using type= typename BubbleSortTupleOfPairsNumbersHelper<BubbleSortTupleOfPairsNumbersAux<Tuple>,Nmax,N+1>::type;

};

template<Integer Nmax>
class BubbleSortTupleOfPairsNumbersHelper<std::tuple<std::tuple<>>,Nmax,Nmax>
{
 public:
  using type=std::tuple<>;
};



template <Integer N,std::size_t ... Is>
constexpr auto index_sequence_shift (std::index_sequence<Is...> const &)
   -> decltype( std::index_sequence<(Is+N)...>{} );

template <std::size_t N,std::size_t ShiftN>
using make_index_sequence_shift = decltype(index_sequence_shift<ShiftN>(std::make_index_sequence<N>{}));



template <typename... T, std::size_t... I>
auto sub_tuple_aux(const std::tuple<T...>& t, std::index_sequence<I...>)
-> decltype(std::make_tuple(std::get<I>(t)...))
{
  return std::make_tuple(std::get<I>(t)...);
}

template <std::size_t N_start,std::size_t N_end, typename... Ts>
// std::enable_if_t< (N_start<=N_end && 0<=N_end && 0<=N_start && N_end<sizeof...(Ts) && N_end>=0) , SubTupleType<N_start,N_end,std::tuple<Ts...>>> 
std::enable_if_t< (N_start<=N_end && 0<=N_end && 0<=N_start && N_end<sizeof...(Ts) && N_end>=0) , SubTupleType<N_start,N_end,Ts...>> 
sub_tuple(const std::tuple<Ts...>& t) 
{
  return sub_tuple_aux(t, make_index_sequence_shift<N_end-N_start+1,N_start>());
}

template <Integer N_start,Integer N_end, typename... Ts>
// std::enable_if_t< (N_start<=N_end && 0<=N_end && 0>N_start && N_end<sizeof...(Ts)&& N_end>=0), SubTupleType<0,N_end,std::tuple<Ts...>>> 
std::enable_if_t< (N_start<=N_end && 0<=N_end && 0>N_start && N_end<sizeof...(Ts)&& N_end>=0), SubTupleType<0,N_end,Ts...>> 
sub_tuple(const std::tuple<Ts...>& t) 
{
  return sub_tuple_aux(t, make_index_sequence_shift<N_end+1,0>());
}

template <Integer N_start,Integer N_end, typename... Ts>
// std::enable_if_t< (N_start<=N_end && 0<=N_end && 0<=N_start && N_end>=sizeof...(Ts)), SubTupleType<N_start,sizeof...(Ts)-1,std::tuple<Ts...>>> 
std::enable_if_t< (N_start<=N_end && 0<=N_end && 0<=N_start && N_end>=sizeof...(Ts)), SubTupleType<N_start,sizeof...(Ts)-1,Ts...>> 
sub_tuple(const std::tuple<Ts...>& t) 
{
  return sub_tuple_aux(t, make_index_sequence_shift<sizeof...(Ts)-N_start,N_start>());
}

template <Integer N_start,Integer N_end, typename... Ts>
// std::enable_if_t< (N_start<=N_end && 0<=N_end  && 0>N_start && N_end>=sizeof...(Ts)), SubTupleType<0,sizeof...(Ts)-1,std::tuple<Ts...>>> 
std::enable_if_t< (N_start<=N_end && 0<=N_end  && 0>N_start && N_end>=sizeof...(Ts)), SubTupleType<0,sizeof...(Ts)-1,Ts...>> 
sub_tuple(const std::tuple<Ts...>& t) 
{
  return sub_tuple_aux(t, make_index_sequence_shift<sizeof...(Ts),0>());
}


template <Integer N_start,Integer N_end, typename... Ts>
std::enable_if_t< ((N_start<=N_end &&  N_end<0 )||
                   // (N_start<=N_end && N_start>=sizeof...(Ts)-1)||
                   (N_start>N_end)), 
                  std::tuple<>> 
sub_tuple(const std::tuple<Ts...>& t) 
{
  return std::tuple<>();
}

template <Integer N,typename T, typename... Ts>
auto tuple_change_element(const std::tuple<Ts...>& tuple,const T& t) 
{
  static_assert( (0<=N && N<sizeof...(Ts))&&"In tuple_change_element, N must be between 0 and tuple_size-1");
  return std::tuple_cat(sub_tuple<0,N-1>(tuple),
                        std::tuple<T>(t),
                        sub_tuple<N+1,sizeof...(Ts)-1>(tuple) );
}


template<typename T, typename Tuple>
constexpr auto& get_tuple_element_of_type(const Tuple& tuple)
{
  // std::cout<< TypeToTupleElementPosition<T,Tuple>::value <<std::endl;
  return tuple_get< TypeToTupleElementPosition<T,Tuple>::value >(tuple);
}



template<Integer N, typename T>
std::enable_if_t<(N==-1),std::tuple<T>> tuple_cat_unique_aux(const T&t)
{return std::tuple<T>(t);}

template<Integer N,typename T>
std::enable_if_t<(N>=0),std::tuple<>> tuple_cat_unique_aux(const T&t)
{return std::tuple<>();}

template<typename T, typename Tuple>
constexpr auto tuple_cat_unique( Tuple& tuple, const T&t)
{
  return std::tuple_cat(tuple,tuple_cat_unique_aux<TypeToTupleElementPosition<T,Tuple>::value>(t));
}

template<typename T>
constexpr auto tuple_cat_unique( std::tuple<> tuple, const T&t)
{
  return std::tuple<T>(t);
}



template<typename Derived>
class Expression;

template<typename T>
class UnaryPlus;

template<typename T>
class UnaryMinus;

template<typename Left,typename Right>
class Addition;

template<typename Left,typename Right>
class Subtraction;

template<typename Left,typename Right>
class Multiplication;

template<typename Left,typename Right>
class Division;

template<typename...Parameters>
class CompositeOperator;


template<typename...Parameters>
class Evaluation;

// template<typename Left,typename Right,Integer QR>
// class L2DotProductIntegral; 

// template<typename Form, typename...Forms>
// class ShapeFunctionsCollection;


// template<typename Left,typename Right,Integer QR, typename...Forms>
// constexpr auto build_tuple_of_evals_aux_aux(const std::tuple<>& null, 
//                            const L2DotProductIntegral<Left,Right,QR>& l2prod,
//                                  ShapeFunctionsCollection<Forms...>&shape_functions)
// {
//   using L2dot=L2DotProductIntegral<Left,Right,QR>;

//   return Evaluation<Expression<L2dot>,ShapeFunctionsCollection<Forms...>>
//         // (Eval(L2dot(l2prod.left(),l2prod.right()),shape_functions));
//         (Eval(l2prod,shape_functions));
// }

// template<typename Left,typename Right,Integer QR, typename...Forms>
// constexpr auto build_tuple_of_evals_aux_aux(const L2DotProductIntegral<Left,Right,QR>& l2prod,
//                                             const std::tuple<>& null,
//                                                   ShapeFunctionsCollection<Forms...>&shape_functions)
// {
//   using L2dot=L2DotProductIntegral<Left,Right,QR>;

//   return Evaluation<Expression<L2dot>,ShapeFunctionsCollection<Forms...>>
// //         (Eval(L2dot(l2prod.left(),l2prod.right()),shape_functions));
//   (Eval(l2prod,shape_functions));

// }

// template<typename Left1,typename Right1,Integer QR1,
//          typename Left2,typename Right2,Integer QR2, typename...Forms>
// constexpr auto build_tuple_of_evals_aux_aux(const L2DotProductIntegral<Left1,Right1,QR1>& l2prod1, 
//                            const L2DotProductIntegral<Left2,Right2,QR2>& l2prod2,
//                                  ShapeFunctionsCollection<Forms...>&shape_functions)
// {
//   // using L2dot1=L2DotProductIntegral<Left1,Right1,QR1>;
//   // using L2dot2=L2DotProductIntegral<Left2,Right2,QR2>;


//   return Eval(l2prod1+l2prod2,shape_functions);
//         // Addition<Expression<Evaluation<Expression<L2dot1>,ShapeFunctionsCollection<Forms...> > >,
//         //           Expression<Evaluation<Expression<L2dot2>,ShapeFunctionsCollection<Forms...> > > >
//         // (Eval(L2dot1(l2prod1.left(),l2prod1.right()),shape_functions),
//         //  Eval(L2dot2(l2prod2.left(),l2prod2.right()),shape_functions));
// }

// template<typename Left, typename Left1,typename Right1,Integer QR1, typename...Forms>
// constexpr auto build_tuple_of_evals_aux_aux(const Left& left, 
//                            const L2DotProductIntegral<Left1,Right1,QR1>& l2prod,
//                                  ShapeFunctionsCollection<Forms...>&shape_functions)
// {
//   using L2dot1=L2DotProductIntegral<Left1,Right1,QR1>;
//   // return Addition<Expression<Left>,
//   //                 Expression< Evaluation<Expression<L2dot1>,ShapeFunctionsCollection<Forms...> >>>
//   //       (left,
//   //        Eval(L2dot1(l2prod.left(),l2prod.right()),shape_functions));
//   const auto& left_=left.expression();
//   return Eval(left_+l2prod,shape_functions);
// }

// // template<typename Left1,typename Right1,Integer QR1,typename Right, typename...Forms>
// // constexpr auto build_tuple_of_evals_aux_aux(const L2DotProductIntegral<Left1,Right1,QR1>& l2prod,
// //                            const Right& right,
// //                                  ShapeFunctionsCollection<Forms...>&shape_functions)
// // {
// //   using L2dot1=L2DotProductIntegral<Left1,Right1,QR1>;
// //   std::cout<<"---------"<<l2prod.left().eval()<<std::endl;
// //   return Addition<Expression<Evaluation<Expression<L2dot1>,ShapeFunctionsCollection<Forms...>>>,
// //                   Expression<Right>>
// //         (Eval(L2dot1(l2prod.left(),l2prod.right()),shape_functions),
// //          right());
// // }






// template<typename TupleOfPairsNumbers,typename Tuple, typename Left,typename Right,Integer QR, typename...Forms>
// constexpr auto build_tuple_of_evals_aux(const Tuple& tuple,
//                           const L2DotProductIntegral<Left,Right,QR>& l2prod,
//                                 ShapeFunctionsCollection<Forms...>&shape_functions)
// {
//   using L2=L2DotProductIntegral<Left,Right,QR>;
//   using TestTrialNumbers=typename L2::TestTrialNumbers;
//   auto tuple_nth=tuple_get< TypeToTupleElementPosition<TestTrialNumbers,TupleOfPairsNumbers>::value >(tuple);
//   auto new_elem=build_tuple_of_evals_aux_aux(tuple_nth,l2prod,shape_functions);
//   // todo fixme tuple add std::tuple_cat
//   return tuple_change_element<TypeToTupleElementPosition<TestTrialNumbers,TupleOfPairsNumbers>::value>
//    (tuple, decltype(new_elem)(new_elem));
// }


// template<typename TupleOfPairsNumbers,typename Tuple, typename Left,typename Right, typename...Forms>
// constexpr auto build_tuple_of_evals_aux(const Tuple& tuple,
//                           const Addition<Expression<Left>,Expression<Right>>& addition,
//                                 ShapeFunctionsCollection<Forms...>&shape_functions)
// {
//   auto tuple_new=build_tuple_of_evals_aux<TupleOfPairsNumbers>(tuple,addition.left(),shape_functions);
//   return build_tuple_of_evals_aux<TupleOfPairsNumbers>(tuple_new,addition.right(),shape_functions);
// }


// template<typename TupleOfPairsNumbers, typename Expr, typename...Forms>
// constexpr auto build_tuple_of_evals(const Expr& expr,ShapeFunctionsCollection<Forms...>&shape_functions)
// {
//   using emptytuple=TupleOfType<TupleTypeSize<TupleOfPairsNumbers>::value,std::tuple<> > ;
//   return build_tuple_of_evals_aux<TupleOfPairsNumbers,emptytuple>(emptytuple(),expr,shape_functions);
// }
















































// template<typename T>
// class Transposed;

template<typename MixedSpace, Integer N, typename OperatorType>
class Test;

template<typename MixedSpace, Integer N, typename OperatorType>
class Trial;

template<typename ConstType, typename...Inputs>
class ConstantTensor;

template<typename FullSpace, Integer M,typename Operator,typename FuncType>
class Function;









// template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Operator,typename Expr>
// constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,Operator>& t, 
//                                               const Expr& expr)
// {
//     using T=TestOrTrial<MixedSpace,N,Operator>;
//     using FunctionSpace=typename T::FunctionSpace;
//     using FromElementFunctionSpacesToFirstSpaceTupleType=typename FunctionSpace::FromElementFunctionSpacesToFirstSpaceTupleType;
//     constexpr Integer FirstSpace=GetType<FromElementFunctionSpacesToFirstSpaceTupleType,T::value>::value;
//     return TestOrTrial<MixedSpace,FirstSpace,Operator>(t.spaces_ptr());
// }

// template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename Operator>
// constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
//                                               const Operator& op)
// {
//     using T=TestOrTrial<MixedSpace,N,Operator>;
//     using FunctionSpace=typename T::FunctionSpace;
//     using FromElementFunctionSpacesToFirstSpaceTupleType=typename FunctionSpace::FromElementFunctionSpacesToFirstSpaceTupleType;
//     constexpr Integer FirstSpace=GetType<FromElementFunctionSpacesToFirstSpaceTupleType,T::value>::value;
//     return TestOrTrial<MixedSpace,FirstSpace,Operator>(t.spaces_ptr());
// }

// // template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename Operator>
// // constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
// //                                               const Transposed<Expression<Operator>>& expr)
// // {
// //   auto e=form_of_composite_operator_aux(t,expr.derived());
// //   // decltype(expr.derived()) eee(6);
// //   return Transpose(e);

// //   // using T=TestOrTrial<MixedSpace,N,Operator>;
// //   // using FunctionSpace=typename T::FunctionSpace;
// //   // using FromElementFunctionSpacesToFirstSpaceTupleType=typename FunctionSpace::FromElementFunctionSpacesToFirstSpaceTupleType;
// //   // constexpr Integer FirstSpace=GetType<FromElementFunctionSpacesToFirstSpaceTupleType,T::value>::value;  
// //   // return Transpose(TestOrTrial<MixedSpace,FirstSpace,Operator>(t.spaces_ptr()));
// // }

// template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename ConstType,typename...Inputs>
// constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
//                                               const ConstantTensor<ConstType,Inputs...>& constant)
// {return constant;}

// // template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename ConstType,typename...Inputs>
// // constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
// //                                               const Transposed<Expression<ConstantTensor<ConstType,Inputs...>>>& transposed_constant)
// // {
// //   return transposed_constant;
// // }

// // template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename ConstType,typename...Inputs>
// // constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
// //                                               const TraceOperator<Expression<ConstantTensor<ConstType,Inputs...>>>& trace_constant)
// // {
// //   return trace_constant;
// // }


// template<template<class,Integer,class > class TestOrTrial, template<class>class Unary,
// typename MixedSpace,Integer N,typename Expr,typename ConstType,typename...Inputs>
// constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
//                                               const Unary<Expression<ConstantTensor<ConstType,Inputs...>>>& unary_operator_applied_to_constant)
// {
//     return unary_operator_applied_to_constant;
// }

// template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename FullSpace, Integer M,typename Operator,typename FuncType>
// constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
//                                               const Function<FullSpace,M,Operator,FuncType>& func)
// {return func;}

// // template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename FullSpace, Integer M,typename Operator,typename FuncType>
// // constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
// //                                               const Transposed<Expression<Function<FullSpace,M,Operator,FuncType>>>& transposed_func)
// // {return transposed_func;}


// template<template<class,Integer,class > class TestOrTrial, template<class>class Unary,
// typename MixedSpace,Integer N,typename Expr,typename FullSpace, Integer M,typename Operator,typename FuncType>
// constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
//                                               const Unary<Expression<Function<FullSpace,M,Operator,FuncType>>>& unary_operator_applied_to_func)
// {return unary_operator_applied_to_func;}












// template<template<class,Integer,class > class TestOrTrial, template<class>class Unary,
// typename MixedSpace,Integer N,typename Expr,typename T>
// constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
//                                               const Unary<Expression<T>>& expr)
// {
//     auto e=form_of_composite_operator_aux(t,expr.derived());
//     // decltype(expr.derived()) eee(6);
//     return Unary<Expression<decltype(e)>>(e);
// }

// // template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename T>
// // constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
// //                                               const UnaryPlus<Expression<T>>& expr)
// // {
// //   auto e=form_of_composite_operator_aux(t,expr.derived());
// //   // decltype(expr.derived()) eee(6);
// //   return +e;
// // }

// // template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename T>
// // constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
// //                                               const UnaryMinus<Expression<T>>& expr)
// // {
// //   auto e=form_of_composite_operator_aux(t,expr.derived());
// //   return -e;
// // }



// template<template<class,Integer,class > class TestOrTrial, template<class,class>class Binary,
// typename MixedSpace,Integer N,typename Expr,typename Left,typename Right>
// constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
//                                               const Binary<Expression<Left>,Expression<Right>>& expr)
// {
//     auto left=form_of_composite_operator_aux(t,expr.left());
//     auto right=form_of_composite_operator_aux(t,expr.right());
    
//     return Binary<Expression<decltype(left)>,Expression<decltype(right)>>(left,right);
// }


// // template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename Left,typename Right>
// // constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
// //                                               const Addition<Expression<Left>,Expression<Right>>& expr)
// // {
// //   auto left=form_of_composite_operator_aux(t,expr.left());
// //   auto right=form_of_composite_operator_aux(t,expr.right());

// //   return left+right;
// // }

// // template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename Left,typename Right>
// // constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
// //                                               const Subtraction<Expression<Left>,Expression<Right>>& expr)
// // {
// //   auto left=form_of_composite_operator_aux(t,expr.left());
// //   auto right=form_of_composite_operator_aux(t,expr.right());

// //   return left-right;
// // }

// // template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename Left,typename Right>
// // constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
// //                                               const Multiplication<Expression<Left>,Expression<Right>>& expr)
// // {
// //   auto left=form_of_composite_operator_aux(t,expr.left());
// //   auto right=form_of_composite_operator_aux(t,expr.right());

// //   return left*right;
// // }

// // template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename Left,typename Right>
// // constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
// //                                               const Division<Expression<Left>,Expression<Right>>& expr)
// // {
// //   auto left=form_of_composite_operator_aux(t,expr.left());
// //   auto right=form_of_composite_operator_aux(t,expr.right());

// //   return left/right;
// // }


// template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr>
// constexpr auto form_of_composite_operator(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t)
// {
//     return form_of_composite_operator_aux(t,t.composite_operator().composite_operator());
// }

// template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Operator>
// constexpr auto form_of_composite_operator(const TestOrTrial<MixedSpace,N,Operator>& t)
// {
//     return t;
// }




// template<typename QuadratureRule, typename Tuple,typename Other>
// constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple, const Other& other)
// {
//  return tuple;
// }

// // template<typename QuadratureRule, typename Tuple,typename ConstType,typename...Inputs>
// // constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple, const ConstantTensor<ConstType,Inputs...>& other)
// // {
// //  return tuple;
// // }


// // template<typename QuadratureRule, typename Tuple,typename FullSpace, Integer M,typename Operator,typename FuncType>
// // constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple, const Function<FullSpace,M,Operator,FuncType>& other)
// // {
// //  return tuple;
// // }

// // template<typename QuadratureRule, typename Tuple,typename MixedSpace, Integer N, typename Operator>
// // constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple, const Trial<MixedSpace,N,Operator>& other)
// // {
// //  return tuple;
// // }
// // template<typename QuadratureRule, typename Tuple,typename MixedSpace, Integer N, typename Operator>
// // constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple, const Test<MixedSpace,N,Operator>& other)
// // {
// //  return tuple;
// // }

// template<typename QuadratureRule, typename Tuple,
//          template<class,Integer,class > class TestOrTrial_,typename MixedSpace, Integer N, typename Expr>
// constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple, 
//                            const TestOrTrial_<MixedSpace,N,CompositeOperator<Expression<Expr>>>& testortrial)
// {
//  // check if already exists test or trial of the same input with quadrature ruel
//   using TestOrTrial=TestOrTrial_<MixedSpace,N,CompositeOperator<Expression<Expr>>>;
//   auto composite_testortrial=
//        form_of_composite_operator(Test<MixedSpace,N,CompositeOperator<Expression<Expr>>>(testortrial.spaces_ptr(),testortrial.composite_operator()));
//   auto eval_composite=Evaluation<Expression<decltype(composite_testortrial)>,QuadratureRule>(composite_testortrial);
//   auto tuple_nth=tuple_get<TestOrTrial::value>(tuple);
  
//   // decltype(composite_testortrial) eee(5);
//   return tuple_change_element<TestOrTrial::value>(tuple,tuple_cat_unique(tuple_nth,decltype(eval_composite)(eval_composite)));
// }









// template<typename QuadratureRule, template<class>class Unary, typename Tuple,typename T>
// constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple,
//                            const Unary<Expression<T>>& unary)
// {
//   return build_tuple_of_combination_functions_aux_aux<QuadratureRule>(tuple,unary.derived());
// }

// // template<typename QuadratureRule, typename Tuple,typename T>
// // constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple,
// //                            const UnaryPlus<Expression<T>>& plusexpr)
// // {
// //   return build_tuple_of_combination_functions_aux_aux<QuadratureRule>(tuple,plusexpr.derived());
// // }

// // template<typename QuadratureRule, typename Tuple,typename T>
// // constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple,
// //                            const UnaryMinus<Expression<T>>& minusexpr)
// // {
// //   return build_tuple_of_combination_functions_aux_aux<QuadratureRule>(tuple,minusexpr.derived());
// // }


// template<typename QuadratureRule, template<class,class>class Binary, typename Tuple,typename Left,typename Right>
// constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple,
//                            const Binary<Expression<Left>,Expression<Right>>& binary)
// {
//   auto tuple_new=build_tuple_of_combination_functions_aux_aux<QuadratureRule>(tuple,binary.left());
//   return build_tuple_of_combination_functions_aux_aux<QuadratureRule>(tuple_new,binary.right());
// }


// // template<typename QuadratureRule, typename Tuple,typename Left,typename Right>
// // constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple,
// //                            const Addition<Expression<Left>,Expression<Right>>& addition)
// // {
// //   auto tuple_new=build_tuple_of_combination_functions_aux_aux<QuadratureRule>(tuple,addition.left());
// //   return build_tuple_of_combination_functions_aux_aux<QuadratureRule>(tuple_new,addition.right());
// // }

// // template<typename QuadratureRule, typename Tuple,typename Left,typename Right>
// // constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple,
// //                            const Subtraction<Expression<Left>,Expression<Right>>& subtraction)
// // {
// //   auto tuple_new=build_tuple_of_combination_functions_aux_aux<QuadratureRule>(tuple,subtraction.left());
// //   return build_tuple_of_combination_functions_aux_aux<QuadratureRule>(tuple_new,subtraction.right());
// // }

// // template<typename QuadratureRule, typename Tuple,typename Left,typename Right>
// // constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple,
// //                            const Multiplication<Expression<Left>,Expression<Right>>& multiplication)
// // {
// //   auto tuple_new=build_tuple_of_combination_functions_aux_aux<QuadratureRule>(tuple,multiplication.left());
// //   return build_tuple_of_combination_functions_aux_aux<QuadratureRule>(tuple_new,multiplication.right());
// // }

// // template<typename QuadratureRule, typename Tuple,typename Left,typename Right>
// // constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple,
// //                            const Division<Expression<Left>,Expression<Right>>& division)
// // {
// //   auto tuple_new=build_tuple_of_combination_functions_aux_aux<QuadratureRule>(tuple,division.left());
// //   return build_tuple_of_combination_functions_aux_aux<QuadratureRule>(tuple_new,division.right());
// // }


// // template<typename QuadratureRule, typename Tuple,typename T>
// // constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple,
// //                            const Transposed<Expression<T>>& transposed_expr)
// // {
// //   return build_tuple_of_combination_functions_aux_aux<QuadratureRule>(tuple,transposed_expr.derived());
// // }





// template<typename Tuple, typename Left,typename Right,Integer QR>
// constexpr auto build_tuple_of_combination_functions_aux(const Tuple& tuple,
//                                                         const L2DotProductIntegral<Left,Right,QR>& l2prod)
// {
//   using L2=L2DotProductIntegral<Left,Right,QR>;
//   using QRule=typename L2::QRule;
//   auto tuple_new=build_tuple_of_combination_functions_aux_aux<QRule>(tuple,l2prod.left());
//   return build_tuple_of_combination_functions_aux_aux<QRule>(tuple_new,l2prod.right());
// }


// template<typename Tuple, typename Left,typename Right>
// constexpr auto build_tuple_of_combination_functions_aux(const Tuple& tuple,
//                           const Addition<Expression<Left>,Expression<Right>>& addition)
// {
//   auto tuple_new=build_tuple_of_combination_functions_aux(tuple,addition.left());
//   return build_tuple_of_combination_functions_aux(tuple_new,addition.right());
// }


// // template<template<class,class>class AdditionOrSubtraction,typename Tuple, typename Left,typename Right>
// // constexpr auto build_tuple_of_combination_functions_aux(const Tuple& tuple,
// //                           const AdditionOrSubtraction<Expression<Left>,Expression<Right>>& addition)
// // {
// //   auto tuple_new=build_tuple_of_combination_functions_aux(tuple,addition.left());
// //   return build_tuple_of_combination_functions_aux(tuple_new,addition.right());
// // }




// template<typename Tuple,typename Form>
// constexpr auto build_tuple_of_combination_functions_form(const Tuple& tuple,const Form& form)
// {
//   return build_tuple_of_combination_functions_aux(tuple,form);
// }

// template<typename Tuple,typename Form,typename...Forms>
// constexpr auto build_tuple_of_combination_functions_form(const Tuple& tuple,const Form& form, const Forms&...forms)
// {
//   auto tuple_new= build_tuple_of_combination_functions_aux(tuple,form);
//   return build_tuple_of_combination_functions_form(tuple_new,forms...);
// }


// template<Integer Nmax,typename Form,typename...Forms>
// constexpr auto build_tuple_of_combination_functions(const Form& form, const Forms&...forms)
// {
//   using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
//   return build_tuple_of_combination_functions_form(emptytuple(),form,forms...);
// }
















template<typename ConstType,typename...Inputs>
class ConstantTensor;

template<typename FullSpace,Integer N,typename Operator_,typename FuncType>
class Function;

// template<typename...Ts>
// class FormOfCompositeOperatorAuxType;


// template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename Operator>
// class FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,Operator>,Expr>
// {
// public:
//   // using type=TestOrTrial<MixedSpace,N,Operator>;
//   using type=Test<MixedSpace,N,Operator>;
// };

// template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename Operator>
// class FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Operator>
// {
// public:
//   // using type=TestOrTrial<MixedSpace,N,Operator>;
//   using type=Test<MixedSpace,N,Operator>;
// };

// template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename ConstType,typename...Inputs>
// class FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,ConstantTensor<ConstType,Inputs...>>
// {
// public:
//   using type=ConstantTensor<ConstType,Inputs...>;
// };

// template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename FullSpace, Integer M,typename Operator,typename FuncType>
// class FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Function<FullSpace,M,Operator,FuncType>>
// {
// public:
//   using type=Function<FullSpace,M,Operator,FuncType>;
// };

// template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename T>
// class FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,UnaryPlus<Expression<T>>>
// {
// public:
//   using TT=typename FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,T>::type;
//   using type=UnaryPlus<Expression<TT>>;
// };

// template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename T>
// class FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,UnaryMinus<Expression<T>>>
// {
// public:
//   using TT=typename FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,T>::type;
//   using type=UnaryMinus<Expression<TT>>;
// };

// template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename Left,typename Right>
// class FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Addition<Expression<Left>,Expression<Right>>>
// {
// public:
//   using LeftT=typename FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Left>::type;
//   using RightT=typename FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Right>::type;
//   using type=Addition<Expression<LeftT>,Expression<RightT>>;
// };

// template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename Left,typename Right>
// class FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Subtraction<Expression<Left>,Expression<Right>>>
// {
// public:
//   using LeftT=typename FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Left>::type;
//   using RightT=typename FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Right>::type;
//   using type=Subtraction<Expression<LeftT>,Expression<RightT>>;
// };

// template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename Left,typename Right>
// class FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Multiplication<Expression<Left>,Expression<Right>>>
// {
// public:
//   using LeftT=typename FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Left>::type;
//   using RightT=typename FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Right>::type;
//   using type=Multiplication<Expression<LeftT>,Expression<RightT>>;
// };

// template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename Left,typename Right>
// class FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Division<Expression<Left>,Expression<Right>>>
// {
// public:
//   using LeftT=typename FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Left>::type;
//   using RightT=typename FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Right>::type;
//   using type=Division<Expression<LeftT>,Expression<RightT>>;
// };




// template<typename...Ts>
// class FormOfCompositeOperatorType;

// template<template<class,Integer,class > class TestOrTrial_, typename MixedSpace,Integer N,typename Expr>
// class FormOfCompositeOperatorType<TestOrTrial_< MixedSpace,N,CompositeOperator<Expression<Expr>> >>
// {
// public:
//   using type=typename FormOfCompositeOperatorAuxType<TestOrTrial_<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Expr>::type;
// };























// template<typename Number1,typename Number2,Integer X>
// class MaxNumberHelper;

// template<Integer M,Integer N,Integer X>
// class MaxNumberHelper<Number<M>,Number<N>,X>
// {
//  public:
//   using type=typename std::conditional<(M==N), Number<M>,Number<X>>::type;
// };

// template<typename Number1,typename Number2, Integer X>
// using MaxNumber=typename MaxNumberHelper<Number1,Number2,X>::type;


template<typename Tuple,typename NewTuple,Integer Nmax, Integer N>
class FromElementFunctionSpacesToFirstSpaceTupleTypeAux;


template<typename Tuple>
class FromElementFunctionSpacesToFirstSpaceTupleTypeAux<Tuple,std::tuple<>,0,0>
{
public:
  using type=std::tuple<Number<0>>;
};


template<typename Tuple,typename NewTuple,Integer Nmax>
class FromElementFunctionSpacesToFirstSpaceTupleTypeAux<Tuple,NewTuple,Nmax,Nmax>
{
public:
  static constexpr Integer N1=GetType<NewTuple,TupleTypeSize<NewTuple>::value-1>::value;
  static constexpr Integer N2=GetType<Tuple,Nmax>::value;
  using type=typename std::conditional<(N1<N2) , 
                                           TupleCatType< NewTuple, std::tuple<Number<Nmax>> >,
                                           NewTuple>::type;
};

template<typename Tuple,Integer Nmax>
class FromElementFunctionSpacesToFirstSpaceTupleTypeAux<Tuple,std::tuple<>,Nmax,0>
{
public:
  using NewTuple=std::tuple<Number<0>>;
  using type=typename FromElementFunctionSpacesToFirstSpaceTupleTypeAux<Tuple,NewTuple,Nmax,1>::type;
};


template<typename Tuple,typename NewTuple_,Integer Nmax, Integer N>
class FromElementFunctionSpacesToFirstSpaceTupleTypeAux
{
public:
  static constexpr Integer N1=GetType<NewTuple_,TupleTypeSize<NewTuple_>::value-1>::value;
  static constexpr Integer N2=GetType<Tuple,N>::value;
  using NewTuple=typename std::conditional<(N1<N2) , 
                                           TupleCatType< NewTuple_, std::tuple<Number<N>> >,
                                           NewTuple_>::type;
  using type=typename FromElementFunctionSpacesToFirstSpaceTupleTypeAux<Tuple,NewTuple,Nmax,N+1>::type;
};




template<typename Tuple>
class FromElementFunctionSpacesToFirstSpaceTupleTypeHelper
{
public:
  static constexpr Integer Nmax=TupleTypeSize<Tuple>::value-1;
  using type=typename FromElementFunctionSpacesToFirstSpaceTupleTypeAux<Tuple,std::tuple<>,Nmax,0>::type;
};

template<typename Tuple>
using FromElementFunctionSpacesToFirstSpaceTupleType=typename FromElementFunctionSpacesToFirstSpaceTupleTypeHelper<Tuple>::type;





//////////////////////////////////////////////////////////////////////////
////// Compile-time cat of Vector or Arrays
//////////////////////////////////////////////////////////////////////////


template<Integer... Is> struct seq{};
template<Integer N, Integer... Is>
struct gen_seq : gen_seq<N-1, N-1, Is...>{};
template<Integer... Is>
struct gen_seq<0, Is...> : seq<Is...>{};

template<typename T, Integer N>
constexpr Array<T, N> concat(const Array<T, N>& a){
  return a;
}

template<typename T, Integer N1, Integer... I1, Integer N2, Integer... I2>
// Expansion pack
constexpr Array<T, N1+N2> concat(const Array<T, N1>& a1, const Array<T, N2>& a2, seq<I1...>, seq<I2...>){
  return { a1[I1]..., a2[I2]... };
}

template<typename T, Integer N1, Integer N2>
// Initializer for the recursion
constexpr Array<T, N1+N2> concat(const Array<T, N1>& a1, const Array<T, N2>& a2){
  return concat(a1, a2, gen_seq<N1>{}, gen_seq<N2>{});
}


template<typename T, Integer N1, Integer N2, Integer...Ns>
// Initializer for the recursion
constexpr auto concat(const Array<T, N1>& a1, const Array<T, N2>& a2, const Array<T,Ns>&...vectors){
  return concat(concat(a1, a2, gen_seq<N1>{}, gen_seq<N2>{}),vectors... ) ;
}


template<typename T, Integer N>
constexpr Vector<T, N> concat(const Vector<T, N>& a){
  return a;
}

template<typename T, Integer N1, Integer... I1, Integer N2, Integer... I2>
// Expansion pack
constexpr Vector<T, N1+N2> concat(const Vector<T, N1>& a1, const Vector<T, N2>& a2, seq<I1...>, seq<I2...>){
  return { a1[I1]..., a2[I2]... };
}

template<typename T, Integer N1, Integer N2>
// Initializer for the recursion
constexpr Vector<T, N1+N2> concat(const Vector<T, N1>& a1, const Vector<T, N2>& a2){
  return concat(a1, a2, gen_seq<N1>{}, gen_seq<N2>{});
}


template<typename T, Integer N1, Integer N2, Integer...Ns>
// Initializer for the recursion
constexpr auto concat(const Vector<T, N1>& a1, const Vector<T, N2>& a2, const Vector<T,Ns>&...vectors){
  return concat(concat(a1, a2, gen_seq<N1>{}, gen_seq<N2>{}),vectors... ) ;
}












}

#endif