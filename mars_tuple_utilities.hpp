


#ifndef MARS_TUPLE_UTILITIES_HPP
#define MARS_TUPLE_UTILITIES_HPP

#include "mars_base.hpp"

namespace mars{


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

template <Integer N, typename... Ts>
using GetType=typename GetHelper<N,Ts...>::type;




template <Integer N, Integer M, typename... Ts>
class GetHelper2;

template <Integer M, typename T, typename... Ts>
class GetHelper2< 0, M, std::tuple<T,Ts...> >
{
public:
  using type =GetType<M, T >;
};


template <Integer N, Integer M, typename T, typename... Ts>
class GetHelper2< N, M, std::tuple<T,Ts...> >
{
public:
  using type =typename GetHelper2<N-1,M, std::tuple<Ts...> >::type;
};



template <Integer N,Integer M, typename... Ts>
using GetType2=typename GetHelper2<N,M,Ts...>::type;




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
class TupleRemoveType;

template<typename T, typename...Ts>
class TupleRemoveType<T,std::tuple<Ts...>> 
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
 using type= std::tuple<TupleCatType<GetType<Nmax,Ts>...>>;
};

template<Integer N,Integer Nmax, typename ... Ts>
class TupleOfTupleCatTypeHelper
{
public:
using type= 
decltype(std::tuple_cat(std::declval<std::tuple<TupleCatType< GetType<N,Ts>...>>>(),
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
 using type= std::tuple<RemoveTupleDuplicates<TupleCatType<GetType<Nmax,Ts>...>>>;
};

template<Integer N,Integer Nmax, typename ... Ts>
class RemoveTupleOfTupleDuplicatesHelper
{
public:
    //TupleRemovesingleType<std::tuple<>,>
using type= 
decltype(std::tuple_cat(std::declval<std::tuple< RemoveTupleDuplicates<TupleCatType< GetType<N,Ts>...>>>>(),
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
 class ElementPositiion
{public:
 static_assert(N>=0, " negative integer");
 static_assert(N<TupleTypeSize<All>::value, " exceeding tuple dimension");
 using AllType=GetType<N,All>;
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
 using typeLeft=typename std::conditional< (N==0), std::tuple<>, SubTupleType<0,N-1 , Ts... > >::type;
 using typeRight=SubTupleType<N+1,sizeof...(Ts)+1, Ts...>;//typename std::conditional< (N==sizeof...(Ts)+1), std::tuple<>, SubTupleType<N,sizeof...(Ts)+1, Ts...> >::type;
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
using type=std::tuple<TupleChangeType<Nchange,Tchange,GetType<Nmax,Tuple>>>;
};

template<Integer N,Integer Nmax,typename Tuple,Integer Nchange,typename Tchange>
class provahelper
{
 public:
 using single_type=std::tuple<TupleChangeType<Nchange,Tchange,GetType<N,Tuple>>>;

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
 using T=GetType<Nmax,TupleOfTuple>;
 using type= std::tuple<prova<T,Nchange,Tchange >>; 
};



template<Integer N,Integer Nmax,Integer Nchange,typename Tchange, typename TupleOfTuple>
class TupleOfTupleChangeTypeHelper
{
public:
using T=GetType<N,TupleOfTuple>;
using single_type= std::tuple<prova<T,Nchange,Tchange >>;
using type=decltype(std::tuple_cat(std::declval<single_type>(),
                           std::declval<typename TupleOfTupleChangeTypeHelper<N+1,Nmax,Nchange,Tchange,TupleOfTuple>::type>()));
};


template<Integer Nchange,typename Tchange,typename TupleOfTuple>
using TupleOfTupleChangeType=typename TupleOfTupleChangeTypeHelper<0,TupleTypeSize<TupleOfTuple>::value-1,Nchange,Tchange,TupleOfTuple>::type;







//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// Create Tuple (for each function space) of tuples of shape functions related to the pair of pair<operator,quadrature_rule>
///// The input is the tuple of tuple of pair<operator,quadrature_rule> 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename Elem, typename BaseFunctioSpace, typename Operator, typename QuadratureRule>
class ShapeFunctionOperatorDependent;

template<typename Elem, typename BaseFunctioSpace, typename Tuple, Integer Nmax,Integer N>
class TupleShapeFunctionCreate;

template<typename Elem, typename BaseFunctioSpace, typename Tuple, Integer Nmax>
class TupleShapeFunctionCreate<Elem,BaseFunctioSpace,Tuple,Nmax,Nmax>
{
 public:

  using Nelem=GetType<Nmax,Tuple>;
  using Operator=GetType<0,Nelem>;
  using QuadratureRule=GetType<1,Nelem>;
  using type=std::tuple< ShapeFunctionOperatorDependent<Elem, BaseFunctioSpace,Operator,QuadratureRule > >;
};



template<typename Elem, typename BaseFunctioSpace, typename Tuple, Integer Nmax,Integer N=0>
class TupleShapeFunctionCreate
{
 public:
  using Nelem=GetType<N,Tuple>;
  using Operator=GetType<0,Nelem>;
  using QuadratureRule=GetType<1,Nelem>;
  using single_type=std::tuple< ShapeFunctionOperatorDependent<Elem, BaseFunctioSpace,Operator,QuadratureRule > >;
  using type=decltype(std::tuple_cat(std::declval<single_type>(),
                             std::declval<typename TupleShapeFunctionCreate<Elem,BaseFunctioSpace,Tuple,Nmax,N+1>::type>()));
};


template<typename TupleSpaces, typename TupleOperatorsAndQuadrature, Integer Nmax,Integer N>
class TupleOfTupleShapeFunctionCreate;



template<typename TupleSpaces,typename TupleOperatorsAndQuadrature, Integer Nspaces>
class TupleOfTupleShapeFunctionCreate<TupleSpaces,TupleOperatorsAndQuadrature, Nspaces, Nspaces>
{
public:
using OperatorAndQuadrature=GetType<Nspaces,TupleOperatorsAndQuadrature>;
static constexpr Integer Nmax=TupleTypeSize<OperatorAndQuadrature>::value-1;
using Space=GetType<Nspaces,TupleSpaces>;
using Elem=GetType<0,Space>;
using FunctionSpace=GetType<1,Space>;
using type=std::tuple<typename TupleShapeFunctionCreate< Elem,FunctionSpace, OperatorAndQuadrature,Nmax>::type>;
};

template<typename TupleSpaces,typename TupleOperatorsAndQuadrature, Integer Nspaces, Integer N=0>
class TupleOfTupleShapeFunctionCreate
{
public:
using OperatorAndQuadrature=GetType<N,TupleOperatorsAndQuadrature>;
static constexpr Integer Nmax=TupleTypeSize<OperatorAndQuadrature>::value-1;
using Space=GetType<N,TupleSpaces>;
using Elem=GetType<0,Space>;
using FunctionSpace=GetType<1,Space>;
using single_type=std::tuple<typename TupleShapeFunctionCreate< Elem,FunctionSpace, OperatorAndQuadrature,Nmax>::type>;
using type=decltype(std::tuple_cat(std::declval<single_type>(),
                           std::declval<typename TupleOfTupleShapeFunctionCreate<TupleSpaces,TupleOperatorsAndQuadrature,Nspaces,N+1>::type>()));
};



template<typename TupleSpaces,typename TupleOperatorsAndQuadrature>
using TupleOfTupleShapeFunctionType=typename TupleOfTupleShapeFunctionCreate<TupleSpaces,TupleOperatorsAndQuadrature,TupleTypeSize<TupleSpaces>::value-1,0>::type;

template<typename TupleSpaces,typename TupleOperatorsAndQuadrature>
class ShapeFunctionAssembly
{
 public:
   using type=TupleOfTupleShapeFunctionType<TupleSpaces,TupleOperatorsAndQuadrature>;
   static constexpr Integer Nmax=TupleTypeSize<TupleSpaces>::value-1;

   /////////////// reference initialization of shape functions of all spaces, all operators and all quadrature rules used 

   template<Integer Nmax_aux,Integer N,typename Tuple>
   typename std::enable_if< (N>Nmax_aux) ,void>::type 
   init_reference_aux_aux(Tuple& t){}

   template<Integer Nmax_aux,Integer N,typename Tuple>
   typename std::enable_if< (N<=Nmax_aux) ,void>::type 
   init_reference_aux_aux(Tuple& t) 
   {auto& t_nth=std::get<N>(t); 
    t_nth.init_reference();
    init_reference_aux_aux<Nmax_aux,N+1,Tuple>(t);}

   template<Integer N>
   typename std::enable_if< (N>Nmax) ,void>::type 
   init_reference_aux(){}

   template<Integer N>
   typename std::enable_if< (N<=Nmax) ,void>::type 
   init_reference_aux() 
   {auto& tuple_nth=std::get<N>(tuple_);
    using tuple_nth_type=decltype(tuple_nth);
    init_reference_aux_aux<TupleTypeSize<tuple_nth_type>::value-1,0,tuple_nth_type>(tuple_nth);
    init_reference_aux<N+1>(); }

   void init_reference(){init_reference_aux<0>();}
   




   template<Integer Nmax_aux,Integer N,typename Tuple,typename Jacobian>
   typename std::enable_if< (N>Nmax_aux) ,void>::type 
   init_aux_aux(const Jacobian&J,Tuple& t){}

   template<Integer Nmax_aux,Integer N,typename Tuple,typename Jacobian>
   typename std::enable_if< (N<=Nmax_aux) ,void>::type 
   init_aux_aux(const Jacobian&J,Tuple& t) 
   {auto& t_nth=std::get<N>(t); 
    t_nth.init(J);
    init_aux_aux<Nmax_aux,N+1,Tuple,Jacobian>(J,t);}

   template<Integer N,typename Jacobian>
   typename std::enable_if< (N>Nmax) ,void>::type 
   init_aux(const Jacobian&J){}

   template<Integer N,typename Jacobian>
   typename std::enable_if< (N<=Nmax) ,void>::type 
   init_aux(const Jacobian&J) 
   {auto& tuple_nth=std::get<N>(tuple_);
    using tuple_nth_type=decltype(tuple_nth);
    init_aux_aux<TupleTypeSize<tuple_nth_type>::value-1,0,tuple_nth_type>(J,tuple_nth);
    init_aux<N+1>(J); }

   template<typename Jacobian>
   void init(const Jacobian&J){init_aux<0>(J);}




   const type& get()const{return tuple_;};
         type&  get()     {return tuple_;};

   template<Integer N,Integer M>
   const GetType2<N,M,type> & get()const
   {const auto& tuple_nth=std::get<N>(tuple_);
    return std::get<M>(tuple_nth);};
   template<Integer N,Integer M>
         GetType2<N,M,type>&  get()    
   { auto& tuple_nth=std::get<N>(tuple_);
    return std::get<M>(tuple_nth);};



 private:
   type tuple_;
};

// template<typename TupleSpaces,typename TupleOperatorsAndQuadrature, Integer Nspaces>
// init_reference()
// sfod_grad.init_reference();



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// Compile-time Max and Min.
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


}

#endif