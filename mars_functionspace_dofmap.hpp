#ifndef MARS_FUNCTIONSPACE_DOFMAP_HPP
#define MARS_FUNCTIONSPACE_DOFMAP_HPP


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Written by Gabriele Rovi (April 2019)                                                                                            ////////
////// We define the following class:                                                                                                   ////////
////// 1) EntitiesOfFunctionSpace<Integer Dim, Integer ManifoldDim, Integer FEFamily, Integer Order>:                                   ////////                                        
//////    given the function space FEFamily of order Order,                                                                             ////////
//////    it builds a tuple of the entities related to its Dim_dofs                                                                     ////////
////// 2) EntitiesOfFunctionSpaceType is the type of EntitiesOfFunctionSpace                                                            ////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "mars_base.hpp"
#include "mars_elementfunctionspace.hpp"
#include "mars_elem_to_sub_elem.hpp"
#include "mars_tuple_utilities.hpp"
#include "mars_array.hpp"

#include "mars_tracker.hpp"

namespace mars{
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// FunctionSpaceDofsPerElem:                                                                                     //////////////
/////// For the given entity of dimension entity[N] of FunctionSpace, returns the corresponding number of dofs        //////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename FunctionSpace,Integer N=FunctionSpace::entity.size()-1>
class 
FunctionSpaceDofsPerElem
{
     public: 
     using Elem = typename FunctionSpace::Elem; 
     static_assert(N>0," FunctionSpaceDofsPerElem N >0");
     static constexpr std::size_t ManifoldDim=FunctionSpace::ManifoldDim;
     static constexpr std::size_t n_components=FunctionSpace::NComponents;
     static constexpr std::size_t entity_dim=FunctionSpace::entity[N];
     static constexpr std::size_t dofs_per_entity=FunctionSpace::dofs_per_entity[N];
     static constexpr std::size_t dofs_per_elem=n_components * 
                                                dofs_per_entity   * 
                                                ElemEntityCombinations<Elem,entity_dim>::value;
     
     static constexpr std::size_t value=FunctionSpaceDofsPerElem<FunctionSpace,N-1>::value+dofs_per_elem;
};

template<typename FunctionSpace>
class
FunctionSpaceDofsPerElem<FunctionSpace,0>
{
     public:  
     using Elem = typename FunctionSpace::Elem; 
     static constexpr std::size_t N=0;
     static constexpr std::size_t ManifoldDim=FunctionSpace::ManifoldDim;
     static constexpr std::size_t n_components=FunctionSpace::NComponents;
     static constexpr std::size_t entity_dim=FunctionSpace::entity[N];
     static constexpr std::size_t dofs_per_entity=FunctionSpace::dofs_per_entity[N];
     static constexpr std::size_t value=n_components * 
                                        dofs_per_entity *
                                        ElemEntityCombinations<Elem,entity_dim>::value;
};





///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// FunctionSpaceDofsPerSubEntityElem:                                                                       //////////////
/////// For the given entity of dimension entity[N] of FunctionSpace, returns the corresponding number of dofs   //////////////
/////// for a given subentity of the element                                                                     //////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename FunctionSpace,Integer SubElemDim,Integer N=FunctionSpace::entity.size()-1>
class 
FunctionSpaceDofsPerSubEntityElem
{
     public: 
     using Elem = typename FunctionSpace::Elem; 
     static_assert(N>0," FunctionSpaceDofsPerElem N >0");
     static constexpr std::size_t ManifoldDim=FunctionSpace::ManifoldDim;
     static constexpr std::size_t n_components=FunctionSpace::NComponents;
     static constexpr std::size_t entity_dim=FunctionSpace::entity[N];
     static constexpr std::size_t dofs_per_entity=FunctionSpace::dofs_per_entity[N];
     static constexpr std::size_t dofs_per_elem=n_components * 
                                                dofs_per_entity   * 
                                                ElemEntityCombinations<typename ChangeElemDim<Elem,SubElemDim>::type,entity_dim>::value;
     
     static constexpr std::size_t value=FunctionSpaceDofsPerSubEntityElem<FunctionSpace,SubElemDim,N-1>::value+dofs_per_elem;
};

template<typename FunctionSpace,Integer SubElemDim>
class
FunctionSpaceDofsPerSubEntityElem<FunctionSpace,SubElemDim,0>
{
     public:  
     using Elem = typename FunctionSpace::Elem; 
     static constexpr std::size_t N=0;
     static constexpr std::size_t ManifoldDim=FunctionSpace::ManifoldDim;
     static constexpr std::size_t n_components=FunctionSpace::NComponents;
     static constexpr std::size_t entity_dim=FunctionSpace::entity[N];
     static constexpr std::size_t dofs_per_entity=FunctionSpace::dofs_per_entity[N];
     static constexpr std::size_t value=n_components * 
                                        dofs_per_entity *
                                        ElemEntityCombinations<typename ChangeElemDim<Elem,SubElemDim>::type,entity_dim>::value;
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// EntitiesOfFunctionSpaceType:                                                                                  //////////////
/////// For given FunctionSpaces on a given Elem, returns the entities corresponding to its dofs                      //////////////                                        
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// template<typename Elem, Integer FEFamily, Integer Order,Integer N>
// struct EntitiesOfFunctionSpaceType
// {
//       using rest = typename EntitiesOfFunctionSpaceType<Elem,FEFamily,Order,N-1>::type; 
//       using ens = ElemEntity<Elem,ElementFunctionSpace<Elem,FEFamily,Order>::entity[N]>;
//       using tuple_ens=std::tuple<ens>;
//       using type = decltype( std::tuple_cat( std::declval< rest >(), std::declval< tuple_ens >() ) );
// };


// template<typename Elem,Integer FEFamily, Integer Order>
// struct EntitiesOfFunctionSpaceType<Elem,FEFamily,Order,0>
// {
//      using ens = ElemEntity<Elem,ElementFunctionSpace<Elem,FEFamily,Order>::entity[0]>;
//      using type = typename std::tuple<ens>;
// };


// template<typename Elem,Integer FEFamily, 
//          Integer Order, 
//          Integer N=(ElementFunctionSpace<Elem,FEFamily,Order>::entity.size()-1),
//          typename MeshT>
// typename std::enable_if<N==0, 
//                         typename EntitiesOfFunctionSpaceType<Elem,FEFamily,Order,N>::type>::type 
// EntitiesOfFunctionSpace(const MeshT& mesh, 
//                         const std::vector< std::vector<Integer> > &node_2_element)
// {   
//       using ens=ElemEntity<Elem,ElementFunctionSpace<Elem,FEFamily,Order>::entity[N]>;
//       return std::tuple<ens>(ens(mesh,node_2_element));
// }


// template<typename Elem,Integer FEFamily, 
//          Integer Order, 
//          Integer N=(ElementFunctionSpace<Elem,FEFamily,Order>::entity.size()-1),
//          typename MeshT>
// typename std::enable_if< 0<N, 
//                          typename EntitiesOfFunctionSpaceType<Elem,FEFamily,Order,N>::type >::type  
// EntitiesOfFunctionSpace(const MeshT& mesh, 
//                         const std::vector< std::vector<Integer> >&node_2_element )
// {

//       using ens=ElemEntity<Elem,ElementFunctionSpace<Elem,FEFamily,Order>::entity[N]>;
//       return std::tuple_cat(EntitiesOfFunctionSpace<Elem,FEFamily,Order,N-1>(mesh,node_2_element),
//                             std::tuple<ens>(ens(mesh,node_2_element)));
// }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// tuple2array_dofs_offset:                                                                                  //////////////
/////// Transform the dofs_offset tuple to a more suitable array for the functionspace                            //////////////                                        
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// template<Integer N, typename T, typename A>
// typename std::enable_if< 0==N, void>::type
// tuple2array_dofs_offset(const T &tuple,A& array)
// {
//       auto tuple_N=std::get<N>(tuple);
//       auto tuple_N_size=tuple_N.size();
//       array[N].resize(tuple_N_size);
//       for(Integer nn=0;nn<tuple_N_size;nn++)
//           array[N][nn]=tuple_N[nn];
// }


// template<Integer N, typename T, typename A>
// typename std::enable_if< 0<N, void>::type
// tuple2array_dofs_offset(const T &tuple,A& array)
// {
//       auto tuple_N=std::get<N>(tuple);
//       auto tuple_N_size=tuple_N.size();
//       array[N].resize(tuple_N_size);
//       for(Integer nn=0;nn<tuple_N_size;nn++)
//           array[N][nn]=tuple_N[nn];
//       tuple2array_dofs_offset<N-1>(tuple,array);
// }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// initialize_flag_tuple_entities, FlagTuple:                                                                    //////////////
/////// We create a flag tuple (for each space) of tuples (for each entity of the dofs of the FunctionSpace)          ////////////// 
/////// So we can remember, for each space, the entities already visited and we can create the dofmap                 //////////////                                       
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// template<typename Elem, Integer FEFamily, Integer Order, Integer M  , typename ...Args>
// typename std::enable_if< M==sizeof ...(Args),void >::type
// initialize_flag_tuple_entities(std::tuple< Args...> const &tuple, 
//                                std::array<std::vector< std::array<Integer,2> >, 
//                                ElementFunctionSpace<Elem,FEFamily,Order>::entity.size() > & entity_found)
//     {static_assert(M>=0," the tuple must have non negative length ");};

// template<typename Elem, Integer FEFamily, Integer Order, Integer M = 0, typename...Args>
// typename std::enable_if< M<sizeof ...(Args),void >::type
// initialize_flag_tuple_entities(std::tuple< Args...> const &tuple,
//                                std::array<std::vector< std::array<Integer,2> >, 
//                                ElementFunctionSpace<Elem,FEFamily,Order>::entity.size() > & entity_found)
// {
//      static_assert(M>=0," the tuple must have non negative length ");
//      Integer entity_length=std::get<M>(tuple).size();
//      entity_found[M].resize(entity_length,{-1,-1});
//      initialize_flag_tuple_entities<Elem,FEFamily,Order,M+1,Args...>(tuple,entity_found);       
// };

// template<typename Elem,typename FunctionSpace,typename...FunctionSpaces>
// struct FlagTupleType
// {
//      using FS=ElemFunctionSpace<Elem,FunctionSpace>;
//      static constexpr Integer FEFamily=FS::FEFamily;
//      static constexpr Integer Order=FS::Order; 
//      static constexpr Integer entities_nums=FS::entity.size();
//      using rest = typename FlagTupleType<Elem,FunctionSpaces...>::type;
//      using ens  = typename std::array<std::vector<std::array<Integer,2>>, entities_nums>;
//      using type = decltype( std::tuple_cat( std::declval< rest >(), 
//                                             std::declval< std::tuple<ens> >() ) );
// };


// template<typename Elem,typename FunctionSpace>
// struct FlagTupleType<Elem,FunctionSpace>
// {
//  using FS=ElemFunctionSpace<Elem,FunctionSpace>;
//  static constexpr Integer FEFamily=FS::FEFamily;
//  static constexpr Integer Order=FS::Order; 
//  static constexpr Integer entities_nums=FS::entity.size();
//  using ens =  typename std::array<std::vector<std::array<Integer,2>>, entities_nums>;
//  using type = typename std::tuple<ens>;
// };


// template<typename Elem,typename FunctionSpace,typename...FunctionSpaces, typename...Args>
// typename std::enable_if< 0<sizeof...(FunctionSpaces), 
//                          typename FlagTupleType<Elem,FunctionSpace,FunctionSpaces...>::type>::type
// FlagTuple(std::tuple<Args...> tuple)
// {    
//       using FS=ElemFunctionSpace<Elem,FunctionSpace>;
//       constexpr Integer FEFamily=FS::FEFamily;
//       constexpr Integer Order=FS::Order;
//       static constexpr Integer M=sizeof...(FunctionSpaces);
//       static constexpr Integer entities_nums=FS::entity.size();
//       using type = typename std::array<std::vector<std::array<Integer,2>>, entities_nums>;
//       type ens;

//       initialize_flag_tuple_entities<Elem,FEFamily,Order>(std::get<M>(tuple),ens);
//       return std::tuple_cat(FlagTuple<Elem,FunctionSpaces...>(tuple),
//                             std::tuple<type>(ens)
//                            );
// }

// template<typename Elem,typename FunctionSpace, typename...FunctionSpaces, typename...Args>
//  typename std::enable_if< 0==sizeof...(FunctionSpaces), 
//                           typename FlagTupleType<Elem,FunctionSpace>::type>::type
// FlagTuple(std::tuple<Args...> tuple)
// {   
//       using FS=ElemFunctionSpace<Elem,FunctionSpace>;
//       constexpr Integer FEFamily=FS::FEFamily;
//       constexpr Integer Order=FS::Order;
//       static constexpr Integer entities_nums=FS::entity.size();
//       using type = typename std::array<std::vector<std::array<Integer,2>>, entities_nums>;
//       type ens;
//       initialize_flag_tuple_entities<Elem,FEFamily,Order>(std::get<0>(tuple),ens);
//       return std::tuple<type>(ens);
// }

// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// /////// EntitiesOfFunctionSpaceTuple:                                                                             //////////////
// /////// Tuple of EntitiesOfFunctionSpace for more FunctionSpaces on a given Elem                                  //////////////                                        
// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// template<typename Elem,typename FunctionSpace,typename...FunctionSpaces>
// struct EntitiesOfFunctionSpaceTupleType
// {
//      using FS=ElemFunctionSpace<Elem,FunctionSpace>;
//      static constexpr Integer FEFamily=FS::FEFamily;
//      static constexpr Integer Order=FS::Order; 
//      static constexpr Integer N=FS::entity.size()-1;
//      using rest = typename EntitiesOfFunctionSpaceTupleType<Elem,FunctionSpaces...>::type;
//      using ens  = typename EntitiesOfFunctionSpaceType<Elem,FEFamily,Order,N>::type;
//      using type = decltype( std::tuple_cat( std::declval< rest >(), 
//                                             std::declval< std::tuple<ens> >() ) );
// };


// template<typename Elem,typename FunctionSpace>
// struct EntitiesOfFunctionSpaceTupleType<Elem,FunctionSpace>
// {
//      using FS=ElemFunctionSpace<Elem,FunctionSpace>;
//      static constexpr Integer FEFamily=FS::FEFamily;
//      static constexpr Integer Order=FS::Order; 
//      static constexpr Integer N=FS::entity.size()-1;
//      using ens =  typename EntitiesOfFunctionSpaceType<Elem,FEFamily,Order,N>::type ;
//      using type = typename std::tuple<ens>;
// };


// template<typename Elem,typename FunctionSpace,typename...FunctionSpaces, typename MeshT>
// typename std::enable_if< 0<sizeof...(FunctionSpaces), 
//                          typename EntitiesOfFunctionSpaceTupleType<Elem,FunctionSpace,FunctionSpaces...>::type>::type
// EntitiesOfFunctionSpaceTuple(const MeshT& mesh, 
//                              const std::vector< std::vector<Integer> > &node_2_element)
// {    
//       using FS=ElemFunctionSpace<Elem,FunctionSpace>;
//       static constexpr Integer FEFamily=FS::FEFamily;
//       static constexpr Integer Order=FS::Order;
//       static constexpr Integer N=FS::entity.size()-1;
//       using type = typename EntitiesOfFunctionSpaceType<Elem,FEFamily,Order,N>::type;
//       return std::tuple_cat(EntitiesOfFunctionSpaceTuple<Elem,FunctionSpaces...>(mesh,node_2_element),
//                             std::tuple<type>(EntitiesOfFunctionSpace<Elem,FEFamily,Order>(mesh,node_2_element))
//                            );
// }

// template<typename Elem,typename FunctionSpace, typename...FunctionSpaces, typename MeshT>
//  typename std::enable_if< 0==sizeof...(FunctionSpaces), 
//                           typename EntitiesOfFunctionSpaceTupleType<Elem,FunctionSpace>::type>::type
// EntitiesOfFunctionSpaceTuple(const MeshT& mesh, 
//                              const std::vector< std::vector<Integer> > &node_2_element)
// {   
//       using FS=ElemFunctionSpace<Elem,FunctionSpace>;
//       static constexpr Integer FEFamily=FS::FEFamily;
//       static constexpr Integer Order=FS::Order;
//       static constexpr Integer N=FS::entity.size()-1;
//       using type = typename EntitiesOfFunctionSpaceType<Elem,FEFamily,Order,N>::type;
//       return std::tuple<type>(EntitiesOfFunctionSpace<Elem,FEFamily,Order>(mesh,node_2_element));
// }
// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// /////// ElementDofMap, ElementDofMap_LoopEntities:                                                                    //////////////
// /////// We create the element dofmap by looping on the functionspaces and the corresponding entities                  //////////////                                       
// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// ElementDofMap, ElementDofMap_LoopEntities:                                                                    //////////////
/////// We create the element dofmap by looping on the functionspaces and the corresponding entities                  //////////////                                       
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



template<Integer EntityDim,typename ...T>
class DofsOrdering;
template<Integer EntityDim,typename ...T>
class DofsOrdered;





////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// DofsPerElemNums:                                                                                              //////////////
/////// For given FunctionSpaces, returns the total number of dofs for the givene element                             //////////////                                        
/////// For P3=LagrangeSpace of Order=3 and P1=LagrangeSpace of Order=1 in 5D                                         //////////////
/////// arr=[0,5,25,35], because dofs_per_entity(P3)=[1,2,1] and dofs_per_entity(P1)=[1]                              //////////////
/////// We have DofsPerElemNums=5 nodes * 1, 10 edges * 2, 10 triangles * 1, 5 nodes * 1 = 40                         //////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Elem, typename FunctionSpace,typename...FunctionSpaces>
class DofsPerElemNums
{ 
 public: 
 static constexpr Integer value= (DofsPerElemNums<Elem,FunctionSpaces...>::value + FunctionSpaceDofsPerElem<ElemFunctionSpace<Elem,FunctionSpace>>::value);

};

template<typename Elem,typename FunctionSpace>
class DofsPerElemNums<Elem,FunctionSpace>
{
 public: 
 static constexpr Integer value= FunctionSpaceDofsPerElem< ElemFunctionSpace<Elem,FunctionSpace>>::value;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// FunctionSpaceOffSetDofs:                                                                                     //////////////
/////// For FunctionSpace, returns an array with offsets                                                             //////////////
/////// For LagrangeSpace of Order=3 in 5D, arr=[0,5,25]+ addConst, because dofs_per_entity=[1,2,1]                  //////////////
/////// We have 5 nodes * 1, 10 edges * 2, 10 triangles * 1                                                          //////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<Integer MaxSize,typename FunctionSpace,Integer N,Integer AddConst=0>
typename  std::enable_if< (0<N)&&(N==MaxSize-1), void>::type
FunctionSpaceOffSetDofs(std::array<Integer , MaxSize> &arr)
 {   
    arr[N] = AddConst + FunctionSpaceDofsPerElem<FunctionSpace,N-1>::value;
 }

template<Integer MaxSize,typename FunctionSpace,Integer N,Integer AddConst=0>
typename  std::enable_if< (0<N)&&(N<MaxSize-1), void>::type
FunctionSpaceOffSetDofs(std::array<Integer , MaxSize> &arr)
{   
    arr[N] = AddConst + FunctionSpaceDofsPerElem<FunctionSpace,N-1>::value;
    FunctionSpaceOffSetDofs<MaxSize,FunctionSpace,N+1,AddConst>(arr);    
 }



template<Integer MaxSize,typename FunctionSpace,Integer N=0,Integer AddConst=0>
typename  std::enable_if< 0==N && 1<MaxSize, void>::type
FunctionSpaceOffSetDofs(std::array<Integer , MaxSize> &arr)
{   
    arr[0] = AddConst + 0;
    FunctionSpaceOffSetDofs<MaxSize,FunctionSpace,N+1,AddConst>(arr);
}
template<Integer MaxSize,typename FunctionSpace,Integer N=0,Integer AddConst=0>
typename  std::enable_if< 0==N && 1==MaxSize, void>::type
FunctionSpaceOffSetDofs(std::array<Integer , MaxSize> &arr)
{   
    arr[0] = AddConst + 0;
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// OffSetDofs:                                                                                                   //////////////
/////// Generalize FunctionSpaceOffSetDofs to more FunctionSpaces , returns an array with offsets                     //////////////                                        
/////// For P3=LagrangeSpace of Order=3 and P1=LagrangeSpace of Order=1 in 5D                                         //////////////
/////// arr=[0,5,25,35], because dofs_per_entity(P3)=[1,2,1] and dofs_per_entity(P1)=[1]                              //////////////
/////// We have 5 nodes * 1, 10 edges * 2, 10 triangles * 1, 5 nodes * 1                                              //////////////
/////// PS: removed the add_const (not anymore necessary, each space is considered separately)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Elem,typename FunctionSpace,typename...FunctionSpaces>
struct OffsetDofsType
{
    using FS=ElemFunctionSpace<Elem,FunctionSpace>;
    static constexpr Integer entities_nums=FS::entity.size();
    using rest = typename OffsetDofsType<Elem, FunctionSpaces...>::type;
    using ens = std::array<Integer,entities_nums+1>;
    using tuple_ens=std::tuple<ens>;
    using type = decltype( std::tuple_cat(std::declval< tuple_ens >(), std::declval< rest >() ) );
};

template<typename Elem,typename FunctionSpace>
struct OffsetDofsType<Elem,FunctionSpace>
{
 using FS=ElemFunctionSpace<Elem,FunctionSpace>;
 static constexpr Integer entities_nums=FS::entity.size();
 using ens = std::array<Integer,entities_nums+1>;
 using type = typename std::tuple<ens>;
};


template<Integer M=0,typename Elem,typename FunctionSpace,typename...FunctionSpaces>
typename std::enable_if< 0==sizeof...(FunctionSpaces),
                         typename OffsetDofsType<Elem,FunctionSpace>::type >::type
OffsetDofs()
{   
    using FS=ElemFunctionSpace<Elem,FunctionSpace>;
    constexpr const auto entities_nums=FS::entity.size();
    std::array<Integer,entities_nums+1> arr;
    FunctionSpaceOffSetDofs<entities_nums+1,FS,0,M>(arr);
    return std::tuple<decltype(arr)>(arr);
}

template<Integer M=0,typename Elem, typename FunctionSpace,typename...FunctionSpaces>
typename std::enable_if< 0<sizeof...(FunctionSpaces),
                         typename OffsetDofsType<Elem,FunctionSpace,FunctionSpaces...>::type >::type 
OffsetDofs()
{
    using FS=ElemFunctionSpace<Elem,FunctionSpace>;
    constexpr const auto entities_nums=FS::entity.size();
    std::array<Integer,entities_nums+1> arr;

    FunctionSpaceOffSetDofs<entities_nums+1,FS,0,M>(arr);
    // return std::tuple_cat(std::tuple<decltype(arr)>(arr),
    //                       OffsetDofs<M+FunctionSpaceDofsPerElem<FS,entities_nums-1>::value,
    //                                   Elem,FunctionSpaces...>());
    return std::tuple_cat(std::tuple<decltype(arr)>(arr),
                          OffsetDofs<M,
                                      Elem,FunctionSpaces...>());
}




template<Integer N,Integer M, typename FunctionSpace, typename...FunctionSpaces>
typename std::enable_if< 0==sizeof...(FunctionSpaces) , void >::type
function_space_info(Array<Array<Integer,4>, M> &array)
{
 array[N][0]=FunctionSpace::FEFamily;
 array[N][1]=FunctionSpace::Order;
 array[N][2]=FunctionSpace::Continuity;
 array[N][3]=FunctionSpace::NComponents;
};


template<Integer N,Integer M, typename FunctionSpace, typename...FunctionSpaces>
typename std::enable_if< 0<sizeof...(FunctionSpaces) , void >::type
function_space_info(Array<Array<Integer,4>, M> &array)
{
 array[N][0]=FunctionSpace::FEFamily;
 array[N][1]=FunctionSpace::Order;
 array[N][2]=FunctionSpace::Continuity;
 array[N][3]=FunctionSpace::NComponents;
 function_space_info<N+1,M,FunctionSpaces...>(array);
};





template<Integer N=0, typename...Args>
std::enable_if_t<(N==sizeof...(Args)),void>
resize_tuple_of_vector(std::tuple<Args...>& tuple, const Integer size)
{}

template<Integer N=0, typename...Args>
std::enable_if_t<(N<sizeof...(Args)),void>
resize_tuple_of_vector(std::tuple<Args...>& tuple, const Integer size)
{
 tuple_get<N>(tuple).resize(size);
 resize_tuple_of_vector<N+1>(tuple,size);
}
















template<typename Elem>
class Node2ElemMap;








    template<typename FunctionSpace, typename Value, Integer Nmax,Integer N>
    struct TupleMapConstructorHelper;

    template<typename FunctionSpace, typename Value,Integer Nmax>
    struct TupleMapConstructorHelper<FunctionSpace,Value,Nmax,Nmax-1>
    {

      using Elem=typename FunctionSpace::Elem;
      static constexpr Integer entity_points=ElemEntityNPoints<Elem,FunctionSpace::entity[Nmax-1]>::value;
      using Key=std::array<Integer,entity_points>;
      using ens = std::map<Key, Value>;
      using type = typename std::tuple<ens>;
    };

    template<typename FunctionSpace, typename Value,Integer Nmax,Integer N>
    struct TupleMapConstructorHelper
    {
      using Elem=typename FunctionSpace::Elem;
      static constexpr Integer entity_points=ElemEntityNPoints<Elem,FunctionSpace::entity[N]>::value;
      using Key=std::array<Integer,entity_points>;
      using ens = std::map<Key, Value>;
      using tuple_ens=std::tuple<ens>;
      using rest = typename TupleMapConstructorHelper<FunctionSpace,Value,Nmax,N+1>::type; 
      using type = decltype( std::tuple_cat( std::declval< tuple_ens >(),std::declval< rest >() ) );
    };

    template<typename Value,typename FunctionSpace>
    using TupleOfMapConstructorType=typename TupleMapConstructorHelper<FunctionSpace,Value,FunctionSpace::entity.size(),0>::type;
    

    // using TupleParentMap=TupleOfMap<bool,n_entity>;
    // using TupleMap=TupleOfMap<std::shared_ptr<IntegerVector>,n_entity>;
    // using TupleMapDof=TupleOfMap<Integer,n_entity>;


 




    template<typename Value, typename Elem, typename BaseFunctionSpace,typename...BaseFunctionSpaces>
    struct TupleOfTupleMapConstructorHelper;

    template<typename Value, typename Elem, typename BaseFunctionSpace>
    struct TupleOfTupleMapConstructorHelper<Value,Elem,BaseFunctionSpace>
    {
      using ens = TupleOfMapConstructorType<Value,ElemFunctionSpace<Elem,BaseFunctionSpace>>;
      using type = typename std::tuple<ens>;
    };

    template<typename Value, typename Elem, typename BaseFunctionSpace,typename...BaseFunctionSpaces>
    struct TupleOfTupleMapConstructorHelper
    {
      using ens = TupleOfMapConstructorType<Value,ElemFunctionSpace<Elem,BaseFunctionSpace>>;
      using tuple_ens=std::tuple<ens>;
      using rest = typename TupleOfTupleMapConstructorHelper<Value,Elem,BaseFunctionSpaces...>::type; 
      using type = decltype( std::tuple_cat( std::declval< tuple_ens >(),std::declval< rest >() ) );
    };

    
    template<typename Value, typename Elem,typename...BaseFunctionSpaces>
    using TupleOfTupleMapConstructor=typename TupleOfTupleMapConstructorHelper<Value,Elem,BaseFunctionSpaces...>::type;









template<typename MeshT,typename FunctionSpace>
class DofMapSingleSpace
{
public:
   
   // using MeshT=typename FunctionSpace::MeshT;
   using BisectionT=Bisection<MeshT>;
   using Elem=typename MeshT::Elem;
   static constexpr Integer ManifoldDim=MeshT::ManifoldDim;
   static constexpr Integer Dim=MeshT::Dim;
   static constexpr const auto entity= FunctionSpace::entity;
   static constexpr auto n_entity=entity.size();
   static constexpr Integer NComponents=FunctionSpace::NComponents;
   static constexpr Integer n_components_entity=NComponents*n_entity;

   // static constexpr auto Npoints=ArrayEntityNPoints<Elem,FunctionSpace>::value();
   // static constexpr auto EntityCombinations=ArrayEntityCombinations<Elem,FunctionSpace>::value();
   using IntegerVector=std::vector<Integer>;
   // using Key=std::array<Integer,Npoints>;
   // using Value=std::shared_ptr<IntegerVector>;
   // using Map=std::map<Key, Value>;
   // using MapDof=std::map<Key,Integer>;





 // std::tuple<std::map<Key, bool>

 //   using Key=std::array<Integer,Npoints>;
 //   using Value=std::shared_ptr<IntegerVector>;
 //   using Map=std::map<Key, Value>;
 //   using MapDof=std::map<Key,Integer>;

 //        Map map_;
 //      MapDof map_dof_;
 //      std::map<Key, bool> parent_map_;
 //     std::vector<std::vector<Integer>>










   DofMapSingleSpace(MeshT& mesh,Node2ElemMap<MeshT>& node2elem,Bisection<MeshT>& bisection):
   mesh_(mesh),
   node2elem_(node2elem),
   bisection_(bisection)
   {}
   

   template<Integer N, typename Elem2Dofs, typename Map, typename ParentMap,typename MapDof>
   void init_elem_aux_aux(
             Elem2Dofs& dofmap_vec,
             Map& map,
             ParentMap& parent_map,
             MapDof& map_dof,
             Integer& loc_dof_count,
             IntegerVector& count_n_entity_vec,
             const Elem elem )
   {
   // // std::cout<<"___________LEVEL="<<0<<std::endl;

    // // std::cout<< " ElementDofMap_LoopEntities5 N = "<< N <<std::endl;
    using Elem = typename FunctionSpace::Elem;
    constexpr const auto ManifoldDim=FunctionSpace::ManifoldDim;
    constexpr const auto continuity=FunctionSpace::Continuity;
    constexpr const auto n_components=FunctionSpace::NComponents;
    constexpr const auto entity_dim=FunctionSpace::entity[N];
    constexpr const auto dofs_per_entity=FunctionSpace::dofs_per_entity[N];
    constexpr auto entity_points=entity_dim+1;
    constexpr auto manifold_points=ManifoldDim+1;
    constexpr auto combinations_nums=ElemEntityCombinations<Elem,entity_dim>::value;


    std::array<Integer,entity_points> entity_nodes;

    Integer cont;

    const auto& nodes=elem.nodes;


    Integer entity_[entity_points];
    Integer global_old_dof;


    auto& global_dof_count=count_n_entity_vec[0];
    auto& elem_id=elem.id;
    // loop on all the entities of a given entity_dim
    for(Integer entity_iter=0;entity_iter<combinations_nums;entity_iter++)
       { 

        Combinations<manifold_points, entity_points>::generate(entity_iter,entity_);
        for(std::size_t i=0;i<entity_points;i++)
          entity_nodes[i]=nodes[entity_[i]];

        // // std::cout<<"entity_nodes="<<std::endl;
        //   for(int s=0;s<entity_nodes.size();s++)
        //   // std::cout<<entity_nodes[s]<<std::endl;
        //   // std::cout<<std::endl;

        std::sort(entity_nodes.begin(),entity_nodes.end());

        // // std::cout<<"must_reorder="<<DofsOrdering<entity_dim,FunctionSpace>::must_reorder<<std::endl;

        auto ordered_entity_nodes=DofsOrdering<entity_dim,FunctionSpace>::value(nodes,entity_iter);
        
        auto& dof=map_dof[entity_nodes];

    //     // consider more dofs per entity (RT1)
        if(DofsOrdering<entity_dim,FunctionSpace>::must_reorder)
        {
          auto ordered_dofs=DofsOrdered<entity_dim,FunctionSpace>::local_dofs(ordered_entity_nodes,loc_dof_count);
          cont=0;
          // if entity not found yet
          if(!dof || continuity==Discontinuous )
          {
          dof=std::make_shared<Integer>(global_dof_count);
          for(Integer m=0;m<ordered_dofs.size();m++)
          {
            dofmap_vec[elem_id][ordered_dofs[cont]]=global_dof_count;
            // // std::cout<<"(*dofmap_vec)[elem_id][ordered_dofs[cont]]="<<dofmap_vec[elem_id][ordered_dofs[cont]]<<std::endl;
            cont++;
            loc_dof_count++;
            global_dof_count++;
          }
          }
          else
          {    
          global_old_dof=*dof;
          auto ordered_dofs=DofsOrdered<entity_dim,FunctionSpace>::local_dofs(ordered_entity_nodes,loc_dof_count);
          cont =0;
          for(Integer m=0;m<ordered_dofs.size();m++)
          {
            dofmap_vec[elem_id][ordered_dofs[cont]]=global_old_dof;
            // // std::cout<< "(*dofmap_vec)[elem_id][ordered_dofs[cont]]="<<dofmap_vec[elem_id][ordered_dofs[cont]]<<std::endl;
            cont++;
            loc_dof_count++;
            global_old_dof++;
          }
          }
        }
        // otherwise
        else
        {
          // // std::cout<<"do not reorder"<<std::endl;
    //     // if the entity has not been already visited, then create n_components new dofs
        if(!dof || continuity==Discontinuous)
        {
          // // std::cout<<"!dof || continuity==Discontinuous"<<std::endl;
          dof=std::make_shared<Integer>(global_dof_count);
          for(Integer entity_dofs_iter=0;entity_dofs_iter<dofs_per_entity;entity_dofs_iter++)
             for(Integer fs_dim=0;fs_dim<n_components;fs_dim++)
                {            
                 dofmap_vec[elem_id][loc_dof_count]=global_dof_count;

                 // to add ??
                 // space_dofs->push_back(global_dof_count);
                 // // std::cout<<global_dof_count <<std::endl;
                 loc_dof_count++;
                 global_dof_count++;

                 }
        }
        else
        {
          // // std::cout<<"else"<<std::endl;
         global_old_dof=*dof;
         for(Integer entity_dofs_iter=0;entity_dofs_iter<dofs_per_entity;entity_dofs_iter++)
            for(Integer fs_dim=0;fs_dim<n_components;fs_dim++)
               {           
                dofmap_vec[elem_id][loc_dof_count]=global_old_dof;
                loc_dof_count++;
                // // std::cout<<global_old_dof <<std::endl;
                global_old_dof++;
               }
        }
       }
      }


      //  // std::cout<<" DOFMAP VEC5 ==" <<std::endl;
      //  // std::cout<<" dofmap_vec.size() ==" <<dofmap_vec.size()<<std::endl;
     
      //        for(int i=0;i<dofmap_vec.size();i++)
      //         {
      //           for(int j=0;j<dofmap_vec[i].size();j++)
      //         // std::cout<<dofmap_vec[i][j]<<" ";
      //         // std::cout<<std::endl;
      //       }
      // // std::cout<<"after  dofmap_vec.size() ==" <<dofmap_vec.size()<<std::endl;
}
   
































   template<Integer N, typename Elem2Dofs, typename Map, typename ParentMap,typename MapDof>
   void init_fine_elem_aux_aux(
             const Integer& cont_tot,
             Integer& cont_new,
             std::vector<Integer>& entity_used,
             Elem2Dofs& dofmap_vec,
             Map& map,
             ParentMap& parent_map,
             MapDof& map_dof,
             Integer& loc_dof_count,
             IntegerVector& count_n_entity_vec,
             const Elem elem,
             const Integer level )
   {
   // // std::cout<<"___________LEVEL="<<level<<std::endl;
   // // std::cout<<"cont_new="<<cont_new<<std::endl;
   // // std::cout<<"cont_tot="<<cont_tot<<std::endl;

    // // std::cout<< " ElementDofMap_LoopEntities5 N = "<< N <<std::endl;
    using Elem = typename FunctionSpace::Elem;
    constexpr const auto ManifoldDim=FunctionSpace::ManifoldDim;
    constexpr const auto continuity=FunctionSpace::Continuity;
    constexpr const auto n_components=FunctionSpace::NComponents;
    constexpr const auto entity_dim=FunctionSpace::entity[N];
    constexpr const auto dofs_per_entity=FunctionSpace::dofs_per_entity[N];
    constexpr auto entity_points=entity_dim+1;
    constexpr auto manifold_points=ManifoldDim+1;
    constexpr auto combinations_nums=ElemEntityCombinations<Elem,entity_dim>::value;


        // std::cout  << "ELEM===="<<elem.id<<"     with entity==="<<entity_dim<<std::endl;

      // for(std::size_t i=0;i<count_n_entity_vec.size();i++)
      //   // std::cout<<count_n_entity_vec[i]<<" ";

      //   // std::cout  << "\n-----------------------------------\n";

      //   for(const auto &m : map_dof) {
      //     const auto& nodes=m.first;
      //     for(const auto n : nodes) {
      //       // std::cout << n << " ";
      //     }

      //     // std::cout << "-> ";

      //       // std::cout << *m.second << " ";
          

      //     // std::cout << "\n";
      //   }



    std::array<Integer,entity_points> entity_nodes;

    Integer cont;

    const auto& nodes=elem.nodes;


    Integer entity_[entity_points];
    Integer global_old_dof;

    // // std::cout<<"level="<<level<<std::endl;
    // // std::cout<<"count_n_entity_vec"<<std::endl;
    // for(Integer i=0; i<count_n_entity_vec.size();i++ )
    //   // std::cout<<count_n_entity_vec[i]<<" ";
    // // std::cout<<std::endl;

    auto& global_dof_count=count_n_entity_vec[level];
    auto& elem_id=elem.id;
    // loop on all the entities of a given entity_dim
    for(Integer entity_iter=0;entity_iter<combinations_nums;entity_iter++)
       { 

        Combinations<manifold_points, entity_points>::generate(entity_iter,entity_);
        for(std::size_t i=0;i<entity_points;i++)
          entity_nodes[i]=nodes[entity_[i]];


        // std::cout<<"entity_dim="<<entity_dim<<std::endl;

        // std::cout<<"entity_nodes="<<std::endl;
        //   for(int s=0;s<entity_nodes.size();s++)
        //   std::cout<<entity_nodes[s]<<std::endl;
        //   std::cout<<std::endl;

        std::sort(entity_nodes.begin(),entity_nodes.end());

        // // std::cout<<"must_reorder="<<DofsOrdering<entity_dim,FunctionSpace>::must_reorder<<std::endl;
        // std::cout<<"nodes="<<std::endl;
        //   for(int s=0;s<nodes.size();s++)
        //   std::cout<<nodes[s]<<std::endl;
        //   std::cout<<std::endl;


        auto ordered_entity_nodes=DofsOrdering<entity_dim,FunctionSpace>::value(nodes,entity_iter);
        
        auto& dof=map_dof[entity_nodes];

                      // if(cont_new<cont_tot)
                      // {
                      //  child_dof=entity_used[cont_new];
                      //  cont_new++; 
                      // }
                      // else
                      // {child_dof=count_n_entity_vec[level];
                      // global_dof_count++;}

   Integer global_dof_count_tmp;
          // std::cout<<"must reorder="<<DofsOrdering<entity_dim,FunctionSpace>::must_reorder<<std::endl;
    //     // consider more dofs per entity (RT1)
        if(DofsOrdering<entity_dim,FunctionSpace>::must_reorder)
        {
          cont=0;

          auto ordered_dofs=DofsOrdered<entity_dim,FunctionSpace>::local_dofs(ordered_entity_nodes,loc_dof_count);
        // std::cout<<"ordered_entity_nodes="<<std::endl;
        //   for(int s=0;s<ordered_entity_nodes.size();s++)
        //   std::cout<<ordered_entity_nodes[s]<<std::endl;
        //   std::cout<<std::endl;
        //  std::cout<<"ordered_dofs="<<std::endl;
        //   for(int s=0;s<ordered_dofs.size();s++)
        //   std::cout<<ordered_dofs[s]<<std::endl;
        //   std::cout<<std::endl;         

          auto order_of_entity_nodes=argsort(entity_nodes);
          // if entity not found yet
          if(!dof || continuity==Discontinuous )
          {
            // std::cout<<"!dof ="<<std::endl;

           if(cont_new<cont_tot)
           {


             // // std::cout<<"cont_new<cont_tot=>>"<<entity_used[cont_new]<<std::endl;
             // for(int s=0;s<entity_used.size();s++)
             //  // std::cout<<entity_used[s]<<" ";
             // // std::cout<<std::endl;
             dof=std::make_shared<Integer>(entity_used[cont_new]);
             cont_new++; 
             // // std::cout<<"cont_new="<<cont_new <<std::endl;
             global_dof_count_tmp=*dof;
            for(Integer m=0;m<ordered_dofs.size();m++)
            {
              dofmap_vec[elem_id][ordered_dofs[cont]]=global_dof_count_tmp;
              // std::cout<<"(dofmap_vec)[elem_id]["<<ordered_dofs[cont]<<"]="<<dofmap_vec[elem_id][ordered_dofs[cont]]<<std::endl;
              cont++;
              loc_dof_count++;
              global_dof_count_tmp++;
            }

           }
           else
           {

            // std::cout<<"==dof ="<<std::endl;
            dof=std::make_shared<Integer>(global_dof_count);
            // // std::cout<<"cont_new>=cont_tot=>>"<<global_dof_count<<std::endl;
            for(Integer m=0;m<ordered_dofs.size();m++)
            {
              dofmap_vec[elem_id][ordered_dofs[cont]]=global_dof_count;
              // std::cout<<"(dofmap_vec)[elem_id]["<<ordered_dofs[cont]<<"]="<<dofmap_vec[elem_id][ordered_dofs[cont]]<<std::endl;
              cont++;
              loc_dof_count++;
              global_dof_count++;
            }

           }


        


          }
          else
          {    
          global_old_dof=*dof;
          auto ordered_dofs=DofsOrdered<entity_dim,FunctionSpace>::local_dofs(ordered_entity_nodes,loc_dof_count);
          cont =0;

          // std::cout<<"global_old_dof= "<<global_old_dof<<std::endl;
          for(Integer m=0;m<ordered_dofs.size();m++)
          {
            dofmap_vec[elem_id][ordered_dofs[cont]]=global_old_dof;
            // std::cout<< "(*dofmap_vec)[elem_id]["<<ordered_dofs[cont]<<"]="<<dofmap_vec[elem_id][ordered_dofs[cont]]<<std::endl;
            cont++;
            loc_dof_count++;
            global_old_dof++;
          }
          }
        }
        // otherwise
        else
        {
    //     // if the entity has not been already visited, then create n_components new dofs
        if(!dof || continuity==Discontinuous)
        {


          // std::cout  << "ELEM===="<<elem.id<<"     with entity==="<<entity_dim<<std::endl;

          //  std::cout<<"!dof ="<<std::endl;
          //  std::cout<<"cont_new ="<<cont_new<<std::endl;
          //  std::cout<<"cont_tot ="<<cont_tot<<std::endl;
           if(cont_new<cont_tot)
           {
            // // std::cout<<"entity_used[cont_new]="<<entity_used[cont_new]<<std::endl;
             dof=std::make_shared<Integer>(entity_used[cont_new]);
             cont_new++; 
             global_dof_count_tmp=*dof;
              for(Integer entity_dofs_iter=0;entity_dofs_iter<dofs_per_entity;entity_dofs_iter++)
                 for(Integer fs_dim=0;fs_dim<n_components;fs_dim++)
                    {            
                     dofmap_vec[elem_id][loc_dof_count]=global_dof_count_tmp;
                     // std::cout<<"if loc_dof_count="<<loc_dof_count<<",  global_dof_count_tmp="<<global_dof_count_tmp<<"   dofmap_vec[elem_id][loc_dof_count]="<<dofmap_vec[elem_id][loc_dof_count]<<std::endl;
                     // to add ??
                     // space_dofs->push_back(global_dof_count);
                     loc_dof_count++;
                     global_dof_count_tmp++;
                     }

           }
           else
           {
            dof=std::make_shared<Integer>(global_dof_count);
            // // std::cout<<"global_dof_count="<<global_dof_count<<std::endl;
            for(Integer entity_dofs_iter=0;entity_dofs_iter<dofs_per_entity;entity_dofs_iter++)
               for(Integer fs_dim=0;fs_dim<n_components;fs_dim++)
                  {            
                   dofmap_vec[elem_id][loc_dof_count]=global_dof_count;
                   // std::cout<<"else loc_dof_count="<<loc_dof_count<<", global_dof_count="<<global_dof_count<<"   dofmap_vec[elem_id][loc_dof_count]="<<dofmap_vec[elem_id][loc_dof_count]<<std::endl;
                   // to add ??
                   // space_dofs->push_back(global_dof_count);
                   loc_dof_count++;
                   global_dof_count++;
                   }
           }

        }
        else
        {
          // // std::cout<<"===dof ="<<std::endl;
         global_old_dof=*dof;
         for(Integer entity_dofs_iter=0;entity_dofs_iter<dofs_per_entity;entity_dofs_iter++)
            for(Integer fs_dim=0;fs_dim<n_components;fs_dim++)
               {           
                dofmap_vec[elem_id][loc_dof_count]=global_old_dof;
                // // std::cout<<"dofmap_vec[elem_id][loc_dof_count]="<<dofmap_vec[elem_id][loc_dof_count]<<std::endl;
                loc_dof_count++;
                global_old_dof++;
               }
        }
       }
      }


      //  // std::cout<<" DOFMAP VEC5 ==" <<std::endl;
      //  // std::cout<<" dofmap_vec.size() ==" <<dofmap_vec.size()<<std::endl;
     
      //        for(int i=0;i<dofmap_vec.size();i++)
      //         {
      //           for(int j=0;j<dofmap_vec[i].size();j++)
      //         // std::cout<<dofmap_vec[i][j]<<" ";
      //         // std::cout<<std::endl;
      //       }
      // // std::cout<<"after  dofmap_vec.size() ==" <<dofmap_vec.size()<<std::endl;
}




   template<Integer N=0, typename Elem2Dofs, typename Map, typename ParentMap,typename MapDof>
   std::enable_if_t<(N==n_entity),void >
   init_fine_elem_aux(
             const Integer& cont_tot,
             Integer& cont_new,
             std::vector<Integer>& entity_used,
             Elem2Dofs& dofmap_vec,
             Map& tuple_map,
             ParentMap& tuple_parent_map,
             MapDof& tuple_map_dof,
             Integer& loc_dof_count,
             IntegerVector& count_n_entity_vec,
             const Elem elem,const Integer level )
   {}


   template<Integer N=0, typename Elem2Dofs, typename Map, typename ParentMap,typename MapDof>
   std::enable_if_t<(N<n_entity),void >
   init_fine_elem_aux(
             const Integer& cont_tot,
             Integer& cont_new,
             std::vector<Integer>& entity_used,
             Elem2Dofs& dofmap_vec,
             Map& tuple_map,
             ParentMap& tuple_parent_map,
             MapDof& tuple_map_dof,
             Integer& loc_dof_count,
             IntegerVector& count_n_entity_vec,
             const Elem elem,const Integer level )
   {

    // std::cout<<"init_fine_elem_aux N="<<N<<std::endl;
    // std::cout<<"count_n_entity_vec ="<<count_n_entity_vec[level]<<std::endl;
    auto& map=tuple_get<N>(tuple_map);
    auto& parent_map=tuple_get<N>(tuple_parent_map);
    auto& map_dof=tuple_get<N>(tuple_map_dof);    
    init_fine_elem_aux_aux<N>(cont_tot,cont_new,entity_used,dofmap_vec,
                          map,parent_map,map_dof,
                          loc_dof_count,count_n_entity_vec,elem,level);
    init_fine_elem_aux<N+1>(cont_tot,cont_new,entity_used,dofmap_vec,
                       tuple_map,tuple_parent_map,tuple_map_dof,
                       loc_dof_count,count_n_entity_vec,elem,level);
   }


   template<typename Elem2Dofs, typename Map, typename ParentMap,typename MapDof>
   void init_fine_elem(
             const Integer& cont_tot,  
             Integer& cont_new,         
             std::vector<Integer>& entity_used,
             Elem2Dofs& dofmap_vec,
             Map& tuple_map,
             ParentMap& tuple_parent_map,
             MapDof& tuple_map_dof,
             IntegerVector& count_n_entity_vec,
             const Elem elem,const Integer level )

   {
    // Integer cont_new=0;
    Integer loc_dof_count=0;
    // // std::cout<<"--------------count_n_entity_vec="<<std::endl;
      // for(std::size_t i=0;i<count_n_entity_vec.size();i++)
      //   // std::cout<<count_n_entity_vec[i]<<" ";
      // // std::cout<<"\n"<<std::endl;
    // loop on spaces
    // // std::cout<<"init_fine_elem level = "<<level<<std::endl;
    init_fine_elem_aux(cont_tot,cont_new,entity_used,dofmap_vec,tuple_map,tuple_parent_map,tuple_map_dof,loc_dof_count,count_n_entity_vec,elem,level);
   }



   template<Integer N=0, typename Map, typename ParentMap,typename MapDof>
   std::enable_if_t<(N==n_entity),void> count_new_dofs_aux(
             Integer& cont,
             std::vector<Integer>& entity_used,
             Map& tuple_map,
             ParentMap& tuple_parent_map,
             MapDof& tuple_map_dof,
             const Elem elem, const Integer level,const MeshT& mesh ,const Tracker&tracker)
   {}



   template<Integer N=0, typename Map, typename ParentMap,typename MapDof>
   std::enable_if_t<(N<n_entity),void> count_new_dofs_aux(
             Integer& cont,
             std::vector<Integer>& entity_used,
             Map& tuple_map,
             ParentMap& tuple_parent_map,
             MapDof& tuple_map_dof,
             const Elem elem, const Integer level ,const MeshT& mesh, const Tracker&tracker)
{
      constexpr const auto entity_dim=FunctionSpace::entity[N];
      constexpr const auto dofs_per_entity=FunctionSpace::dofs_per_entity[N];
      constexpr auto entity_points=entity_dim+1;
      constexpr auto manifold_points=ManifoldDim+1;
      constexpr auto combinations_nums=ElemEntityCombinations<Elem,entity_dim>::value;

      auto& map=tuple_get<N>(tuple_map);
      auto& parent_map=tuple_get<N>(tuple_parent_map);
      auto& map_dof=tuple_get<N>(tuple_map_dof);

      Integer parent_id=get_used_parent_elem_id(elem.id,mesh,tracker);
      // const auto& parent_id=elem.parent_id;

      const auto& parent_elem=mesh_.elem(parent_id);
      const auto& parent_nodes=parent_elem.nodes;
      std::array<Integer,entity_points> parent_entity_nodes;
      std::array<Integer,entity_points> child_entity_nodes;
      
          bool found_parent_entity;
      Integer entity_[entity_points];


        // // std::cout  << "ELEM===="<<elem.id<<"     with entity==="<<entity_dim<<std::endl;

        // // std::cout  << "\n------------------------parent_map-----------\n";

        // for(const auto &m : parent_map) {
        //   const auto& nodes=m.first;
        //   for(const auto n : nodes) {
        //     // std::cout << n << " ";
        //   }

        //   // std::cout << "-> ";

        //     // std::cout << m.second << " ";
          

        //   // std::cout << "\n";
        // }



     
       // find the parent element
       // loop on all its entities parent_entity of dimension EntityDim
       // given entity parent_entity:
       //     loop on all the child elements child_elem
       //          loop on all the entities child_entity
       //               if child_entity==parent_entity
           //                  do nothing, such entity already exists
           //      if no agreement has been found,  cont++
           // then we start creating new dofs for the entities on this level
           // entities which belong to both coarse and fine level are untouched
           // the new ones are created, counting them from n_coarse_entities-cont
       // we loop on all the entities of dimension EntityDim of the element
   

       // const auto& child=parent_elem.children;
       auto child=get_used_children_elem_id(parent_id,mesh,tracker);

     
           for(std::size_t parent_entity=0;parent_entity<combinations_nums;parent_entity++)
            {
         Combinations<manifold_points, entity_points>::generate(parent_entity,entity_);
               // // std::cout<<"\n count_new_dofs_aux parent_entity_nodes"<<std::endl;
             for(std::size_t i=0;i<entity_points;i++)
                    {
                      parent_entity_nodes[i]=parent_nodes[entity_[i]];
                      // // std::cout<<parent_entity_nodes[i]<<" ";
                    }
               // // std::cout<<std::endl;
               std::sort(parent_entity_nodes.begin(),parent_entity_nodes.end());

              // if(!FunctionSpace::Continuity)
              // {
              //   cont++;

              // }
            
              // else 
                if(!parent_map[parent_entity_nodes])
              {
                      // // std::cout<<" entro in !parent_map"<<std::endl;
                    found_parent_entity=false;
                    // loop on all the child elements child_elem          
              for(std::size_t i=0;i<child.size();i++)
              {
               const auto& child_elem=mesh_.elem(child[i]);
               const auto& child_id=child_elem.id;
               const auto& child_nodes=child_elem.nodes;
                     // loop on all the entities child_entity
                   for(std::size_t child_entity=0;child_entity<combinations_nums;child_entity++)
                    {
                 Combinations<manifold_points, entity_points>::generate(child_entity,entity_);

                     for(std::size_t i=0;i<entity_points;i++)
                            child_entity_nodes[i]=child_nodes[entity_[i]];
                       std::sort(child_entity_nodes.begin(),child_entity_nodes.end());
                         // parent_map[parent_entity_nodes]=true;
                       found_parent_entity=std::equal(child_entity_nodes.begin(),child_entity_nodes.end(),parent_entity_nodes.begin());
                   
                   if(found_parent_entity)             
                     {
                      // // std::cout<<"child_entity_nodes"<<std::endl;
                      // for(Integer i=0;i<child_entity_nodes.size();i++)
                        // // std::cout<<child_entity_nodes[i]<<" ";
                      // // std::cout<<std::endl;
                      goto label;
                    }
                    
                   }
                }
                label:
                if(found_parent_entity)
                {
                }
                else
                {

                    parent_map[parent_entity_nodes]=true;
                    entity_used.push_back(*map_dof[parent_entity_nodes]);
                    // // std::cout<<"\n entity_used  "<<*map_dof[parent_entity_nodes]<<std::endl;
                    cont++;
                }
                {}
          }
        }
   // // std::cout<<"count_new_dofs_aux==="<<cont<<std::endl;
    count_new_dofs_aux<N+1>(cont,entity_used,tuple_map,tuple_parent_map,tuple_map_dof,elem,level,mesh,tracker);

}





   template <typename Map, typename ParentMap,typename MapDof>
   void count_new_dofs(
             Integer& cont,
             std::vector<Integer>& entity_used,
             Map& tuple_map,
             ParentMap& tuple_parent_map,
             MapDof& tuple_map_dof,
             const Elem elem, const Integer level,const MeshT& mesh,const Tracker&tracker )
   {
     count_new_dofs_aux(cont,entity_used,tuple_map,tuple_parent_map,tuple_map_dof,elem,level,mesh,tracker);
   }





   // // init refined elem for entitydim
   template<typename Elem2Dofs, typename Map, typename ParentMap,typename MapDof>
   void init_fine_elem_aux(
             Elem2Dofs& elem2dofs,
             Map& tuple_map,
             ParentMap& tuple_parent_map,
             MapDof& tuple_map_dof,
             Integer& loc_dof_count,
             IntegerVector& count_n_entity_vec,
             const Elem elem, const Integer level,const MeshT& mesh, const Tracker&tracker )
   {



          
          // // std::cout<<"count_n_entity_vec.size()-1="<<count_n_entity_vec.size()-1<<std::endl;
      if(count_n_entity_vec.size()-1<level)
        count_n_entity_vec.push_back(count_n_entity_vec[count_n_entity_vec.size()-1]);
      // // std::cout<<"count_n_entity_vec[END]="<<std::endl;

      // for(std::size_t i=0;i<count_n_entity_vec.size();i++)
      //   // std::cout<<count_n_entity_vec[i]<<" ";

      // const auto& parent_id=elem.parent_id;
      Integer parent_id=get_used_parent_elem_id(elem.id,mesh,tracker);
      const auto& parent_elem=mesh_.elem(parent_id);
      Integer cont_tot;
      std::vector<Integer> entity_used;
      cont_tot=0;
      count_new_dofs(cont_tot,entity_used,tuple_map,tuple_parent_map,tuple_map_dof,elem,level,mesh,tracker);
      // // std::cout<<"CONT TOT ==="<<cont_tot<<std::endl;

      // loop on child nodes
      Integer cont_new=0;
      auto child=get_used_children_elem_id(parent_id,mesh,tracker);
      // const auto& child=parent_elem.children;
      // std::cout<<"init_fine_elem_aux ___________LEVEL="<<level<<std::endl;
      // std::cout<<"parent_id ="<<parent_id<<std::endl;




      // for(Integer i=0;i<count_n_entity_vec.size();i++)
      // std::cout<<"count_n_entity_vec[level] ="<<count_n_entity_vec[i]<<std::endl;

      // for(std::size_t i=0;i<elem2dofs.size();i++)
      //   std::cout<<"elem2dofs ="<<elem2dofs[i]<<std::endl;

      Integer tmp=0;
      for(std::size_t i=0;i<child.size();i++)
      {
        for(Integer j=0;j<elem2dofs[child[i]].size();j++)
          tmp+=elem2dofs[child[i]][j];

      }


      if(!FunctionSpace::Continuity )
      {

        if(tmp==0)
        {
        // std::cout<<"elem2dofs[parent_id]="<<elem2dofs[parent_id]<<std::endl;
      

       for(Integer i=0;i<elem2dofs[parent_id].size();i++)
        {
          elem2dofs[child[0]][i]=elem2dofs[parent_id][i];
           // count_n_entity_vec[level]++;
        }



        for(std::size_t i=1;i<child.size();i++)
        {
        for(Integer j=0;j<elem2dofs[parent_id].size();j++)
          {
            elem2dofs[child[i]][j]=count_n_entity_vec[level];
            // std::cout<<"count_n_entity_vec[level]====="<<count_n_entity_vec[level]<<std::endl;
          count_n_entity_vec[level]++;
        }
        }
        }


      }
      else
      {
        for(std::size_t i=0;i<child.size();i++)
        {
         const auto& child_elem=mesh_.elem(child[i]);
         init_fine_elem(cont_tot,cont_new,entity_used,elem2dofs,tuple_map,tuple_parent_map,tuple_map_dof,count_n_entity_vec,child_elem,level);
       }
      }







     
       // find the parent element
       // loop on all its entities parent_entity of dimension EntityDim
       // given entity parent_entity:
       //     loop on all the child elements child_elem
       //          loop on all the entities child_entity
       //               if child_entity==parent_entity
           //                  do nothing, such entity already exists
           //      if no agreement has been found,  cont++
       //     // then we start creating new dofs for the entities on this level
       //     // entities which belong to both coarse and fine level are untouched
       //     // the new ones are created, counting them from n_coarse_entities-cont
       // // we loop on all the entities of dimension EntityDim of the element
       //  // std::cout<<"qui"<<std::endl;

       //        count_n_entity_vec[level]=count_n_entity_vec[level];//-cont;
       //  // std::cout<<"qui2"<<std::endl;
       //        // loop on all the child elements child_elem          
       //  for(std::size_t i=0;i<child.size();i++)
       //  {
       //   const auto& child_elem=mesh_.elem(child[i]);
       //   // std::cout<<"qui3"<<std::endl;
       //   const auto& child_id=child_elem.id;
       //   // std::cout<<"qui4"<<std::endl;
       //   const auto& child_nodes=child_elem.nodes; 
       //   // std::cout<<"qui4.5  ->"<<elem2dofs.size()<<std::endl;         
       //   if(elem2dofs[child_id].size()==0)
       //   {






       //  // std::cout<<"qui5"<<std::endl;
       //         // loop on all the entities child_entity
       //       for(std::size_t child_entity=0;child_entity<combinations_nums;child_entity++)
       //        {
       //     Combinations<manifold_points, entity_points>::generate(child_entity,entity_);
       //         for(std::size_t i=0;i<entity_points;i++)
       //                child_entity_nodes[i]=child_nodes[entity_[i]];
       //           std::sort(child_entity_nodes.begin(),child_entity_nodes.end());
                   
       //             auto& child_map=map[child_entity_nodes];

       //             // if the entity has not been visited yet
       //             if(!child_map)
       //             {
       //              // create the vector
       //                child_map=std::make_shared<IntegerVector>();
       //                auto& child_dof=*map_dof[child_entity_nodes];
       //                if(cont2<cont)
       //                {
       //                 child_dof=entity_used[cont2];
       //                 // // std::cout<<" cont2<cont "<<child_dof<<std::endl;
       //                 cont2++; 
       //                }
       //                else
       //                {child_dof=count_n_entity_vec[level];
       //                 // // std::cout<<" cont2>=cont "<<child_dof<<std::endl;

       //       count_n_entity_vec[level]++;}
       //             }
       //             // // std::cout<<"count_n_entity_vec[level]="<<count_n_entity_vec[level]<<std::endl;
       //             child_map->push_back(child_id);
       //             // std::cout<<" child_id = "<<child_id<<" with "<<*map_dof[child_entity_nodes]<<std::endl;
                   



       //            // DECOMMENTA. RICORDA CHE ARRAY. RICORDA CHE NON BASTA PUSH BACK PERCHE' HO PIU' DOFS E NCOMPONENTS
       //             elem2dofs[child_id].push_back(*map_dof[child_entity_nodes]);
       //    }
       //   }
       //  }
       }














































   // // init refined elem for entitydim
   template<Integer N, typename Elem2Dofs, typename Map, typename ParentMap,typename MapDof>
   void init_elem_aux_aux(
             Elem2Dofs& elem2dofs,
             Map& map,
             ParentMap& parent_map,
             MapDof& map_dof,
             Integer& loc_dof_count,
             IntegerVector& count_n_entity_vec,
             const Elem elem, const Integer level, const MeshT& mesh, const Tracker& tracker )
   {

      constexpr const auto entity_dim=FunctionSpace::entity[N];
      constexpr const auto dofs_per_entity=FunctionSpace::dofs_per_entity[N];
      constexpr auto entity_points=entity_dim+1;
      constexpr auto manifold_points=ManifoldDim+1;
      constexpr auto combinations_nums=ElemEntityCombinations<Elem,entity_dim>::value;


          // // std::cout<<"___________LEVEL="<<level<<std::endl;
          // // std::cout<<"count_n_entity_vec.size()-1="<<count_n_entity_vec.size()-1<<std::endl;
      if(count_n_entity_vec.size()-1<level)
        count_n_entity_vec.push_back(count_n_entity_vec[count_n_entity_vec.size()-1]);
      // // std::cout<<"count_n_entity_vec[END]="<<std::endl;

      // for(std::size_t i=0;i<count_n_entity_vec.size();i++)
        // // std::cout<<count_n_entity_vec[i]<<" ";
      // // std::cout<<std::endl;
      Integer cont=0;
      Integer cont1=0;
      Integer cont2=0;

      // const auto& parent_id=elem.parent_id;
      Integer parent_id=get_used_parent_elem_id(elem.id,mesh,tracker);
      const auto& parent_elem=mesh_.elem(parent_id);
      const auto& parent_nodes=parent_elem.nodes;
      std::array<Integer,entity_points> parent_entity_nodes;
      std::array<Integer,entity_points> child_entity_nodes;
      std::array<Integer,entity_points> entity_used;
          bool found_parent_entity;
      Integer entity_[entity_points];






     
       // find the parent element
       // loop on all its entities parent_entity of dimension EntityDim
       // given entity parent_entity:
       //     loop on all the child elements child_elem
       //          loop on all the entities child_entity
       //               if child_entity==parent_entity
           //                  do nothing, such entity already exists
           //      if no agreement has been found,  cont++
           // then we start creating new dofs for the entities on this level
           // entities which belong to both coarse and fine level are untouched
           // the new ones are created, counting them from n_coarse_entities-cont
       // we loop on all the entities of dimension EntityDim of the element
       // const auto& child=parent_elem.children;
       auto child=get_used_children_elem_id(parent_id,mesh,tracker);
     
           for(std::size_t parent_entity=0;parent_entity<combinations_nums;parent_entity++)
            {
         Combinations<manifold_points, entity_points>::generate(parent_entity,entity_);
               // // std::cout<<std::endl;
             for(std::size_t i=0;i<entity_points;i++)
                    {
                      parent_entity_nodes[i]=parent_nodes[entity_[i]];
                      // // std::cout<<parent_entity_nodes[i]<<" ";
                    }
               // // std::cout<<std::endl;
               std::sort(parent_entity_nodes.begin(),parent_entity_nodes.end());
              

              if(!parent_map[parent_entity_nodes])
              {
                // // std::cout<<" entro in !parent_map"<<std::endl;
              found_parent_entity=false;
              // loop on all the child elements child_elem          
        for(std::size_t i=0;i<child.size();i++)
        {
         const auto& child_elem=mesh_.elem(child[i]);
         const auto& child_id=child_elem.id;
         const auto& child_nodes=child_elem.nodes;
               // loop on all the entities child_entity
             for(std::size_t child_entity=0;child_entity<combinations_nums;child_entity++)
              {
           Combinations<manifold_points, entity_points>::generate(child_entity,entity_);

               for(std::size_t i=0;i<entity_points;i++)
                      child_entity_nodes[i]=child_nodes[entity_[i]];
                 std::sort(child_entity_nodes.begin(),child_entity_nodes.end());
                   // parent_map[parent_entity_nodes]=true;
                 found_parent_entity=std::equal(child_entity_nodes.begin(),child_entity_nodes.end(),parent_entity_nodes.begin());
             
             if(found_parent_entity)             
               goto label;
              
             }
          }
          label:
          if(found_parent_entity)
          {
            // parent_map[parent_entity_nodes]=true;
            // // std::cout<<"found_parent_entity="<<found_parent_entity<<std::endl;
            // entity_used[cont1]=map_dof[parent_entity_nodes];
            // cont1++;
          }
          else
          {

              // // std::cout<<"entity_used["<<cont1<<"]="<<map_dof_[parent_entity_nodes]<<std::endl;
                      parent_map[parent_entity_nodes]=true;
                entity_used[cont]=*map_dof[parent_entity_nodes];
              cont=cont+n_components_entity;
          }
          {}
          }
        }
        // // std::cout<<"qui"<<std::endl;

              count_n_entity_vec[level]=count_n_entity_vec[level];//-cont;
        // // std::cout<<"qui2"<<std::endl;
              // loop on all the child elements child_elem          
        for(std::size_t i=0;i<child.size();i++)
        {
         const auto& child_elem=mesh_.elem(child[i]);
         // // std::cout<<"qui3"<<std::endl;
         const auto& child_id=child_elem.id;
         // // std::cout<<"qui4"<<std::endl;
         const auto& child_nodes=child_elem.nodes; 
         // // std::cout<<"qui4.5  ->"<<elem2dofs.size()<<std::endl;         
         if(elem2dofs[child_id].size()==0)
         {
        // // std::cout<<"qui5"<<std::endl;
               // loop on all the entities child_entity
             for(std::size_t child_entity=0;child_entity<combinations_nums;child_entity++)
              {
           Combinations<manifold_points, entity_points>::generate(child_entity,entity_);
               for(std::size_t i=0;i<entity_points;i++)
                      child_entity_nodes[i]=child_nodes[entity_[i]];
                 std::sort(child_entity_nodes.begin(),child_entity_nodes.end());
                   
                   auto& child_map=map[child_entity_nodes];

                   // if the entity has not been visited yet
                   if(!child_map)
                   {
                    // create the vector
                      child_map=std::make_shared<IntegerVector>();
                      auto& child_dof=*map_dof[child_entity_nodes];
                      if(cont2<cont)
                      {
                       child_dof=entity_used[cont2];
                       // // std::cout<<" cont2<cont "<<child_dof<<std::endl;
                       cont2++; 
                      }
                      else
                      {child_dof=count_n_entity_vec[level];
                       // // std::cout<<" cont2>=cont "<<child_dof<<std::endl;

             count_n_entity_vec[level]++;}
                   }
                   // // std::cout<<"count_n_entity_vec[level]="<<count_n_entity_vec[level]<<std::endl;
                   child_map->push_back(child_id);
                   // // std::cout<<" child_id = "<<child_id<<" with "<<*map_dof[child_entity_nodes]<<std::endl;
                   



                  // DECOMMENTA. RICORDA CHE ARRAY. RICORDA CHE NON BASTA PUSH BACK PERCHE' HO PIU' DOFS E NCOMPONENTS
                   // elem2dofs[child_id].push_back(*map_dof[child_entity_nodes]);
          }
         }
        }
       }
    



   template<Integer N=0, typename Elem2Dofs, typename Map, typename ParentMap,typename MapDof>
   std::enable_if_t<(N==n_entity),void> 
   init_elem_aux(
             Elem2Dofs& elem2dofs,
             Map& map,
             ParentMap& parent_map,
             MapDof& map_dof,
             Integer& local_dof_count,
             IntegerVector& count_n_entity_vec,
             const Elem elem,const Integer level, const MeshT& mesh, const Tracker& tracker )
   {}

   template<Integer N=0, typename Elem2Dofs, typename Map, typename ParentMap,typename MapDof>
   std::enable_if_t<(N<n_entity),void> 
   init_elem_aux(
             Elem2Dofs& elem2dofs,
             Map& map,
             ParentMap& parent_map,
             MapDof& map_dof,
             Integer& local_dof_count,
             IntegerVector& count_n_entity_vec,
             const Elem elem,const Integer level, const MeshT& mesh, const Tracker& tracker )
   {
    auto& m1=tuple_get<N>(map);
    auto& m2=tuple_get<N>(parent_map);
    auto& m3=tuple_get<N>(map_dof);
    // // std::cout<<"INIT ELEM aux ="<<std::endl;
    init_elem_aux_aux<N>(elem2dofs,m1,m2,m3,local_dof_count,count_n_entity_vec,elem,level,mesh,tracker);
    init_elem_aux<N+1>(elem2dofs,map,parent_map,map_dof,local_dof_count,count_n_entity_vec,elem,level,mesh,tracker);
   }


   template<Integer N=0, typename Elem2Dofs, typename Map, typename ParentMap,typename MapDof>
   std::enable_if_t<(N==n_entity),void> 
   init_coarse_elem_aux(
             Elem2Dofs& elem2dofs,
             Map& map,
             ParentMap& parent_map,
             MapDof& map_dof,
             Integer& local_dof_count,
             IntegerVector& count_n_entity_vec,
             const Elem elem)
   {}

   template<Integer N=0, typename Elem2Dofs, typename Map, typename ParentMap,typename MapDof>
   std::enable_if_t<(N<n_entity),void> 
   init_coarse_elem_aux(
             Elem2Dofs& elem2dofs,
             Map& map,
             ParentMap& parent_map,
             MapDof& map_dof,
             Integer& local_dof_count,
             IntegerVector& count_n_entity_vec,
             const Elem elem)

   {
    // // std::cout<<"INIT COARSE ELEM ="<<std::endl;
    auto& m1=tuple_get<N>(map);
    auto& m2=tuple_get<N>(parent_map);
    auto& m3=tuple_get<N>(map_dof);
    // init_elem_aux_aux<N>(elem2dofs,tuple_get<N>(map),tuple_get<N>(parent_map),tuple_get<N>(map_dof),count_n_entity_vec,elem);
    init_elem_aux_aux<N>(elem2dofs,m1,m2,m3,local_dof_count,count_n_entity_vec,elem);
    init_coarse_elem_aux<N+1>(elem2dofs,map,parent_map,map_dof,local_dof_count,count_n_entity_vec,elem);
   }



   template< typename Elem2Dofs, typename Map, typename ParentMap,typename MapDof>
   void init_elem(
             Elem2Dofs& elem2dofs,
             Map& map,
             ParentMap& parent_map,
             MapDof& map_dof,
             IntegerVector& count_n_entity_vec,
             const Elem elem,const Integer level,const MeshT& mesh, const Tracker&tracker)

   {
    Integer local_dof_count=0;
    if(level==0)
      init_coarse_elem_aux(elem2dofs,map,parent_map,map_dof,local_dof_count,count_n_entity_vec,elem);
    else
      {
        // // std::cout<<"INIT ELEM count_n_entity_vec="<<std::endl;

        // for(int i=0;i<count_n_entity_vec.size();i++)
        //   // std::cout<<count_n_entity_vec[i]<<" ";
        // // std::cout<<std::endl;
        // // std::cout<<"level="<<level<<std::endl;

        init_fine_elem_aux(elem2dofs,map,parent_map,map_dof,local_dof_count,count_n_entity_vec,elem,level,mesh,tracker);
        // init_elem_aux(elem2dofs,map,parent_map,map_dof,local_dof_count,count_n_entity_vec,elem,level);
      // Integer cont3;
      // std::vector<Integer> vec;
      // count_new_dofs(cont3,vec,map,parent_map,map_dof,elem,level);
      // // std::cout<<"CONT 3 ==="<<cont3<<std::endl;
      }
    //init_elem_aux1<N,FunctionSpace,FunctionSpaces...>(t1,t2,t3,t4,count_n_entity_vec,elem,level);
   }




   // template< typename Elem2Dofs, typename Map, typename ParentMap,typename MapDof>
   // void init_elem6(
   //           Integer& cont_tot,
   //           Integer& cont_new,
   //           std::vector<Integer>& entity_used,     
   //           Elem2Dofs& elem2dofs,
   //           Map& map,
   //           ParentMap& parent_map,
   //           MapDof& map_dof,
   //           Integer& loc_dof_count,
   //           IntegerVector& count_n_entity_vec,
   //           const Elem elem,const Integer level)

   // {

   //      // init_fine_elem(cont_tot,cont_new,elem2dofs,map,parent_map,map_dof,loc_dof_count,count_n_entity_vec,elem,level);

   //        init_fine_elem_aux(cont_tot,cont_new,entity_used,elem2dofs,map,parent_map,map_dof,loc_dof_count,count_n_entity_vec,elem,level);

   //  //init_elem_aux1<N,FunctionSpace,FunctionSpaces...>(t1,t2,t3,t4,count_n_entity_vec,elem,level);
   // }

private:
    MeshT& mesh_;
    Node2ElemMap<MeshT>& node2elem_;
    Bisection<MeshT>& bisection_;
    // Integer last_n_elements_;
    // std::vector<std::vector<Integer>> elem2dofs;

};








template<typename MeshT,typename...BaseFunctionSpaces>
class DofMapFESpace
{
public:
   using Elem=typename MeshT::Elem;
   using EntitiesTuple=std::tuple<DofMapSingleSpace<MeshT,ElemFunctionSpace<Elem,BaseFunctionSpaces>>...>;
   // static constexpr Integer ManifoldDim=MeshT::ManifoldDim;
   // static constexpr Integer Dim=MeshT::Dim;
   // static constexpr auto entity= FunctionSpace::entity;
   // static constexpr auto n_entity=entity.size();
   // static constexpr auto Npoints=array_entity_n_points<Elem>(entity);
   // static constexpr auto EntityCombinations=array_entity_combinations(entity);
   using IntegerVector=std::vector<Integer>;
   // using Key=std::array<Integer,Npoints>;
   // using Value=std::shared_ptr<IntegerVector>;
   // using Map=std::map<Key, Value>;
   // using MapDof=std::map<Key,Integer>;


    template<typename BFS,typename...BFSs>
    struct ConstructorTupleHelper;

    template<typename BFS>
    struct ConstructorTupleHelper<BFS>
    {
      using ens = DofMapSingleSpace<MeshT,ElemFunctionSpace<Elem,BFS>>;
      using type = typename std::tuple<ens>;
    };

    template<typename BFS,typename...BFSs>
    struct ConstructorTupleHelper
    {
      using rest = typename ConstructorTupleHelper<BFSs...>::type; 
      using ens = DofMapSingleSpace<MeshT,ElemFunctionSpace<Elem,BFS>>;
      using tuple_ens=std::tuple<ens>;
      using type = decltype( std::tuple_cat( std::declval< tuple_ens >(),std::declval< rest >() ) );
    };

    template<typename...BFSs>
    using ConstructorTuple=typename ConstructorTupleHelper<BFSs...>::type;




    template<typename BFS,typename...BFSs>
    std::enable_if_t<(0==sizeof...(BFSs)), ConstructorTuple<BFS> >
    construct_tuple(MeshT& mesh,Node2ElemMap<MeshT>& node2elem,Bisection<MeshT>& bisection)
    {   
      using type=DofMapSingleSpace<MeshT,ElemFunctionSpace<Elem,BFS>>;
      return std::tuple<type>(type(mesh,node2elem,bisection));
    }


    template<typename BFS,typename...BFSs>
    std::enable_if_t<(0<sizeof...(BFSs)), ConstructorTuple<BFS,BFSs...> >
    construct_tuple(MeshT& mesh,Node2ElemMap<MeshT>& node2elem,Bisection<MeshT>& bisection)
    {
      using type=DofMapSingleSpace<MeshT,ElemFunctionSpace<Elem,BFS>>;
      return std::tuple_cat(std::tuple<type>(type(mesh,node2elem,bisection)),
                            construct_tuple<BFSs...>(mesh,node2elem,bisection));
    }







   DofMapFESpace(MeshT& mesh,Node2ElemMap<MeshT>& node2elem,Bisection<MeshT>& bisection):
   mesh_(mesh),
   node2elem_(node2elem),
   bisection_(bisection),
   tuple_dofmapfespace_(construct_tuple<BaseFunctionSpaces...>(mesh,node2elem,bisection))
   {}



template<Integer N, typename BFS,typename...BFSs,typename Elem2Dofs, typename Map, typename ParentMap,typename MapDof>
 std::enable_if_t<(sizeof...(BFSs)==0),void> 
 init_elem_aux(
             Elem2Dofs& elem2dofs,
             Map& map,
             ParentMap& parent_map,
             MapDof& map_dof,
             IntegerVector& count_n_entity_vec,
             const Elem elem,const Integer level,const MeshT& mesh, const Tracker&tracker)

{ 
  auto& dofmapfespace=tuple_get<N>(tuple_dofmapfespace_);
  auto& m0=*tuple_get<N>(elem2dofs);
  auto& m1=tuple_get<N>(map);
  auto& m2=tuple_get<N>(parent_map);
  auto& m3=tuple_get<N>(map_dof);
  // std::cout<<"DofMapFESpace init_elem_aux N== "<<N<<std::endl;
  // // std::cout<<"DofMapFESpace aux level = " << level<<std::endl;
  dofmapfespace.init_elem(m0,m1,m2,m3,count_n_entity_vec,elem,level,mesh,tracker);
}


template<Integer N, typename BFS,typename...BFSs,typename Elem2Dofs, typename Map, typename ParentMap,typename MapDof>
 std::enable_if_t<(sizeof...(BFSs)>0),void> 
 init_elem_aux(             
             Elem2Dofs& elem2dofs,
             Map& map,
             ParentMap& parent_map,
             MapDof& map_dof,
             IntegerVector& count_n_entity_vec,
             const Elem elem,const Integer level,const MeshT& mesh, const Tracker&tracker)

{ 
  auto& dofmapfespace=tuple_get<N>(tuple_dofmapfespace_); 
  auto& m0=*tuple_get<N>(elem2dofs);
  auto& m1=tuple_get<N>(map);
  auto& m2=tuple_get<N>(parent_map);
  auto& m3=tuple_get<N>(map_dof);
  // // std::cout<<"DofMapFESpace aux level = " << level<<std::endl;
  // std::cout<<"DofMapFESpace init_elem_aux N== "<<N<<std::endl;
  dofmapfespace.init_elem(m0,m1,m2,m3,count_n_entity_vec,elem,level,mesh,tracker);
  init_elem_aux<N+1,BFSs...>(elem2dofs,map,parent_map,map_dof,count_n_entity_vec,elem,level,mesh,tracker);
}

 template<typename Elem2Dofs, typename Map, typename ParentMap,typename MapDof>
 void init_elem(
             Elem2Dofs& elem2dofs,
             Map& map,
             ParentMap& parent_map,
             MapDof& map_dof,
             IntegerVector& count_n_entity_vec,
             const Elem elem,const Integer level,const MeshT& mesh, const Tracker&tracker)
{

  // QUI DISTINGUO SE E' LEVEL = 0 OPPURE NO
  // // std::cout<<"dofmap init elem"<<std::endl;
  // // std::cout<<"DofMapFESpace level = " << level<<std::endl;
  init_elem_aux<0,BaseFunctionSpaces...>(elem2dofs,map,parent_map,map_dof,count_n_entity_vec,elem,level,mesh, tracker);

}













// template<Integer N, typename BFS,typename...BFSs,typename Elem2Dofs, typename Map, typename ParentMap,typename MapDof>
//  std::enable_if_t<(sizeof...(BFSs)==0),void> 
//  init_elem_aux6(
//              Integer& cont_tot,
//              Integer& cont_new,
//              std::vector<Integer>& entity_used,  
//              Elem2Dofs& elem2dofs,
//              Map& map,
//              ParentMap& parent_map,
//              MapDof& map_dof,
//              Integer& loc_dof_count,
//              IntegerVector& count_n_entity_vec,
//              const Elem elem,const Integer level)

// { 
//   auto& dofmapfespace=tuple_get<N>(tuple_dofmapfespace_);
//   auto& m0=*tuple_get<N>(elem2dofs);
//   auto& m1=tuple_get<N>(map);
//   auto& m2=tuple_get<N>(parent_map);
//   auto& m3=tuple_get<N>(map_dof);
//   dofmapfespace.init_elem6(cont_tot,cont_new,entity_used,m0,m1,m2,m3,loc_dof_count,count_n_entity_vec,elem,level);
// }


// template<Integer N, typename BFS,typename...BFSs,typename Elem2Dofs, typename Map, typename ParentMap,typename MapDof>
//  std::enable_if_t<(sizeof...(BFSs)>0),void> 
//  init_elem_aux6(    
//              Integer& cont_tot,
//              Integer& cont_new,
//              std::vector<Integer>& entity_used,          
//              Elem2Dofs& elem2dofs,
//              Map& map,
//              ParentMap& parent_map,
//              MapDof& map_dof,
//              Integer& loc_dof_count,
//              IntegerVector& count_n_entity_vec,
//              const Elem elem,const Integer level)

// { 
//   auto& dofmapfespace=tuple_get<N>(tuple_dofmapfespace_); 
//   auto& m0=*tuple_get<N>(elem2dofs);
//   auto& m1=tuple_get<N>(map);
//   auto& m2=tuple_get<N>(parent_map);
//   auto& m3=tuple_get<N>(map_dof);
//   dofmapfespace.init_elem6(cont_tot,cont_new,entity_used,m0,m1,m2,m3,loc_dof_count,count_n_entity_vec,elem,level);
//   init_elem_aux6<N+1,BFSs...>(cont_tot,cont_new,entity_used,elem2dofs,map,parent_map,map_dof,loc_dof_count,count_n_entity_vec,elem,level);
// }

//  template<typename Elem2Dofs, typename Map, typename ParentMap,typename MapDof>
//  void init_elem6(
//              Integer& cont_tot,
//              Integer& cont_new,
//              std::vector<Integer>& entity_used,  
//              Elem2Dofs& elem2dofs,
//              Map& map,
//              ParentMap& parent_map,
//              MapDof& map_dof,
//              Integer& loc_dof_count,
//              IntegerVector& count_n_entity_vec,
//              const Elem elem,const Integer level)
// {

//   // QUI DISTINGUO SE E' LEVEL = 0 OPPURE NO
  
//   init_elem_aux6<0,BaseFunctionSpaces...>(cont_tot,cont_new,entity_used,elem2dofs,map,parent_map,map_dof,loc_dof_count,count_n_entity_vec,elem,level);

// }






//    template<Integer N=0, typename Map, typename ParentMap,typename MapDof>
//    std::enable_if_t<(N==sizeof...(BaseFunctionSpaces)),void> count_new_dofs_aux(
//              Integer& cont,
//              std::vector<Integer>& entity_used,
//              Map& tuple_map,
//              ParentMap& tuple_parent_map,
//              MapDof& tuple_map_dof,
//              const Elem elem, const Integer level )
//   {}

//    template<Integer N=0, typename Map, typename ParentMap,typename MapDof>
//    std::enable_if_t<(N<sizeof...(BaseFunctionSpaces)),void> count_new_dofs_aux(
//              Integer& cont_tot,
//              std::vector<Integer>& entity_used,
//              Map& tuple_map,
//              ParentMap& tuple_parent_map,
//              MapDof& tuple_map_dof,
//              const Elem elem, const Integer level )
//   {

//    auto& dofmapfespace=tuple_get<N>(tuple_dofmapfespace_); 
//    auto& m1=tuple_get<N>(tuple_map);
//    auto& m2=tuple_get<N>(tuple_parent_map);
//    auto& m3=tuple_get<N>(tuple_map_dof);
//    dofmapfespace.count_new_dofs(cont_tot,entity_used,m1,m2,m3,elem,level);
//    // // std::cout<<"cont_tot=="<<cont_tot<<std::endl;
//    count_new_dofs_aux<N+1>(cont_tot,entity_used,tuple_map,tuple_parent_map,tuple_map_dof,elem,level);
//   }










// template<typename Elem2Dofs, typename Map, typename ParentMap,typename MapDof>
//  void init_elem5(
//              Elem2Dofs& elem2dofs,
//              Map& tuple_map,
//              ParentMap& tuple_parent_map,
//              MapDof& tuple_map_dof,
//              IntegerVector& count_n_entity_vec,
//              const Elem elem,const Integer level)
// {

// if(level==0)
// init_elem_aux<0,BaseFunctionSpaces...>(elem2dofs,tuple_map,tuple_parent_map,tuple_map_dof,count_n_entity_vec,elem,level);

// else{
//       if(count_n_entity_vec.size()-1<level)
//         count_n_entity_vec.push_back(count_n_entity_vec[count_n_entity_vec.size()-1]);
//       // // std::cout<<"count_n_entity_vec[END]="<<std::endl;

//       // for(std::size_t i=0;i<count_n_entity_vec.size();i++)
//       //   // std::cout<<count_n_entity_vec[i]<<" ";

//       const auto& parent_id=elem.parent_id;
//       const auto& parent_elem=mesh_.elem(parent_id);
//       Integer cont_tot=0;
//       std::vector<Integer> entity_used;
//       count_new_dofs_aux(cont_tot,entity_used,tuple_map,tuple_parent_map,tuple_map_dof,elem,level);
//       // // std::cout<<"CONT TOT ==="<<cont_tot<<std::endl;

//       // loop on child nodes
//       Integer cont_new=0;
//       Integer loc_dof_count;
//       const auto& child=parent_elem.children;
//         for(std::size_t i=0;i<child.size();i++)
//         {
//           loc_dof_count=0;
//          const auto& child_elem=mesh_.elem(child[i]);
//          init_elem6(cont_tot,cont_new,entity_used,elem2dofs,tuple_map,tuple_parent_map,tuple_map_dof,loc_dof_count,count_n_entity_vec,child_elem,level);
//        }
//   }
// }





private:
    MeshT& mesh_;
    Node2ElemMap<MeshT>& node2elem_;
    Bisection<MeshT>& bisection_;
    EntitiesTuple tuple_dofmapfespace_;

};






// template<typename MeshT,Integer N>
// class ConnectivitySimpliacialMap;

// template<typename MeshT>
// class ConnectivitySimpliacialMapCollection;



// template<Integer N,typename FunctionSpace,typename FlagArray,typename Element,typename MeshT, Integer P, long M,typename OS, typename SpaceDofs>
// // template<Integer N,typename FunctionSpace,typename FlagArray,typename...Args1, long M,typename OS>
// typename std::enable_if< -1<N,void>::type
// ElementDofMap_LoopEntities5(const ConnectivitySimpliacialMap<MeshT,P>& entitiestuple,
//                            FlagArray& flagtuples,
//                            Element& elem,
//                            const Integer& elem_id, 
//                            std::vector<Array<Integer, M>>&  dofmap_vec,
//                            Integer& global_dof_count,
//                            Integer& loc_dof_count,
//                            const OS &dofs_offset,
//                            SpaceDofs& space_dofs,
//                            Integer& dofs_count)
// {

// // // std::cout<< " ElementDofMap_LoopEntities5 N = "<< N <<std::endl;
// using Elem = typename FunctionSpace::Elem;
// constexpr auto ManifoldDim=FunctionSpace::ManifoldDim;
// constexpr auto continuity=FunctionSpace::Continuity;
// constexpr auto n_components=FunctionSpace::NComponents;
// constexpr auto entity_dim=FunctionSpace::entity[N];
// constexpr auto dofs_per_entity=FunctionSpace::dofs_per_entity[N];
// constexpr auto entity_points=entity_dim+1;
// constexpr auto manifold_points=ManifoldDim+1;
// constexpr auto combinations_nums=ElemEntityCombinations<Elem,entity_dim>::value;
// // const     auto& entity=std::get<entity_dim>(entitiestuple);
// // const     auto& elem2entity=entity.elem_2_entity(elem_id);
//           auto& flag=tuple_get<N>(flagtuples);


// std::array<Integer,entity_points> entity_nodes;

// Integer cont;

// // move to the entity of dimension entity[N-1]
// ElementDofMap_LoopEntities5<N-1,FunctionSpace>
//              (entitiestuple,flagtuples,elem,elem_id,dofmap_vec,global_dof_count,loc_dof_count,dofs_offset,space_dofs,dofs_count);

// // // std::cout<<"pre nodes ="<<std::endl;
// const auto& nodes=elem.nodes;

// // // std::cout<<"after nodes ="<<std::endl;

// Integer entity_[entity_points];

// // // std::cout<<"after entity_ ="<<std::endl;


// // loop on all the entities of a given entity_dim
// for(Integer entity_iter=0;entity_iter<combinations_nums;entity_iter++)
//    { 

//     Combinations<manifold_points, entity_points>::generate(entity_iter,entity_);
//     for(std::size_t i=0;i<entity_points;i++)
//       entity_nodes[i]=nodes[entity_[i]];

//     // // std::cout<<"entity_nodes="<<std::endl;
//     //   for(int s=0;s<entity_nodes.size();s++)
//     //   // std::cout<<entity_nodes[s]<<std::endl;
//     //   // std::cout<<std::endl;

//     std::sort(entity_nodes.begin(),entity_nodes.end());

//     // // std::cout<<"must_reorder="<<DofsOrdering<entity_dim,FunctionSpace>::must_reorder<<std::endl;

//     auto ordered_entity_nodes=DofsOrdering<entity_dim,FunctionSpace>::value(nodes,entity_iter);

//     // consider more dofs per entity (RT1)
//     if(DofsOrdering<entity_dim,FunctionSpace>::must_reorder)
//     {

//       // // std::cout<<"ordered_entity_nodes="<<std::endl;
//       // for(int s=0;s<ordered_entity_nodes.size();s++)
//       // // std::cout<<ordered_entity_nodes[s]<<" "<<std::endl;
//       auto ordered_dofs=DofsOrdered<entity_dim,FunctionSpace>::local_dofs(ordered_entity_nodes,loc_dof_count);
//       // // std::cout<<std::endl;
//       // // std::cout<<"ordered_dofs="<<std::endl;
//       // for(int s=0;s<ordered_dofs.size();s++)
//       // // std::cout<<ordered_dofs[s]<<std::endl;
//       // // std::cout<<std::endl;
//       cont=0;
//       // if entity not found yet
//       if(!flag[entity_nodes] || continuity==Discontinuous )
//       {

//       // // std::cout<<"!flag[ordered_entity_nodes]"<<std::endl;
//       flag[entity_nodes]=std::make_shared<std::array<Integer,2>>(std::array<Integer,2>{elem_id,entity_iter});
//       // // std::cout<<"after !flag[ordered_entity_nodes]"<<std::endl;
//       // flag[ordered_entity_nodes][0]=elem_id;
//       // flag[ordered_entity_nodes][1]=entity_iter;
//       for(Integer m=0;m<ordered_dofs.size();m++)
//       {
//         // // std::cout<<"m="<<m<<std::endl;
//         // // std::cout<<"dofmap_vec[elem_id]"<<dofmap_vec.size()<<std::endl;
//         // // std::cout<<"dofmap_vec[elem_id]"<<dofmap_vec[elem_id].size()<<std::endl;
        
//         dofmap_vec[elem_id][ordered_dofs[cont]]=global_dof_count;
//         // // std::cout<<"(*dofmap_vec)[elem_id][ordered_dofs[cont]]="<<dofmap_vec[elem_id][ordered_dofs[cont]]<<std::endl;
//         cont++;
//         loc_dof_count++;
//         // // std::cout<<global_dof_count<<std::endl;
//         global_dof_count++;
//         dofs_count++;
//       }

//       }
//       else
//       {

//        // // std::cout<<"flag[entity_nodes]"<<std::endl;
//       const auto& elem_id_tmp=(*flag[entity_nodes])[0];
//       const auto& iter_tmp=(*flag[entity_nodes])[1];
//             auto& elem_dm=dofmap_vec[elem_id_tmp];
//             auto ex_loc_dof_count=(dofs_offset[N]+dofs_per_entity*iter_tmp*FunctionSpace::NComponents);
//             auto global_old_dof=elem_dm[ex_loc_dof_count];
     
//       // // std::cout<<"elem_id_tmp="<<elem_id_tmp<<std::endl;
//       // // std::cout<<"global_old_dof="<<global_old_dof<<std::endl;
//       for (int i = 0; i < n_components*dofs_per_entity; ++i)
//       {
//         const auto m= i +ex_loc_dof_count;
//         if(elem_dm[m]<global_old_dof)
//         {
//           global_old_dof=elem_dm[m];
//         }
//       }

//       auto ordered_dofs=DofsOrdered<entity_dim,FunctionSpace>::local_dofs(ordered_entity_nodes,loc_dof_count);

//       cont =0;
//       // // std::cout<<"ordered_dofs="<<ordered_dofs<<std::endl;
//       for(Integer m=0;m<ordered_dofs.size();m++)
//       {
//         // // std::cout<<"dofmap_vec[elem_id].size()="<<dofmap_vec[elem_id].size()<<std::endl;
//         // // std::cout<<"ordered_dofs[cont]="<<ordered_dofs[cont]<<std::endl;
//         dofmap_vec[elem_id][ordered_dofs[cont]]=global_old_dof;
//         // // std::cout<< "(*dofmap_vec)[elem_id][ordered_dofs[cont]]="<<dofmap_vec[elem_id][ordered_dofs[cont]]<<std::endl;
//         cont++;
//         loc_dof_count++;
//         global_old_dof++;
//         dofs_count++;
//       }



//       }
     


//     }
//     // otherwise
//     else
//     {

//     // if the entity has not been already visited, then create n_components new dofs
//     if(!flag[entity_nodes] || continuity==Discontinuous)
//     {
//       // // std::cout<<"!flag[ordered_entity_nodes]"<<std::endl;
//       flag[entity_nodes]=std::make_shared<std::array<Integer,2>>(std::array<Integer,2>{elem_id,entity_iter});

//       // flag[ordered_entity_nodes][0]=elem_id;
//       // flag[ordered_entity_nodes][1]=entity_iter;


//       // // std::cout<<"elem_id = "<<elem_id<<std::endl;
//       // // std::cout<<"entity_iter = "<<entity_iter<<std::endl;
//       for(Integer entity_dofs_iter=0;entity_dofs_iter<dofs_per_entity;entity_dofs_iter++)
//          for(Integer fs_dim=0;fs_dim<n_components;fs_dim++)
//             {            
//              dofmap_vec[elem_id][loc_dof_count]=global_dof_count;
//              space_dofs->push_back(global_dof_count);
//              // // std::cout<< "(dofmap_vec)[elem_id][loc_dof_count]="<<dofmap_vec[elem_id][loc_dof_count]<<std::endl;
//              // // std::cout<<global_dof_count<<std::endl;
//              // // std::cout<<"entity_dofs_iter = "<<entity_dofs_iter<<std::endl;
//              // // std::cout<<"global_dof_count = "<<global_dof_count<<std::endl;
//              // // std::cout<<"(dofsmap)="<<(*dofmap_vec)[elem_id][loc_dof_count]<<std::endl;
//              // // std::cout<<" global_dof_count (new)= "<<global_dof_count<<std::endl;
//              loc_dof_count++;
//              global_dof_count++;
//              dofs_count++;
//              }
//       // // std::cout<<"global_dof_count = "<<global_dof_count<<std::endl;
//     }
//     else
//     {
//      // // std::cout<<"flag[ordered_entity_nodes]"<<std::endl;
//      // if the entity has been already visited, find the element (elem_id_tmp) in which was found 
//      // and the corresponding entity_iter (iter_tmp)
//      // we do not need any offset, because dofmap_vec contains  of the entity N, we move to the iter_tmp entity and find in dofmap_vec the already numbered dof
//      const auto& elem_id_tmp=(*flag[entity_nodes])[0];
//      const auto& iter_tmp=(*flag[entity_nodes])[1];
//            auto old_dof=dofmap_vec[elem_id_tmp][dofs_offset[N]+dofs_per_entity*iter_tmp*FunctionSpace::NComponents];
//       // // std::cout<<" elem_id = "<<elem_id<<std::endl;
//       // // std::cout<<" entity_iter = "<<flag[entity_id][1]<<std::endl  ;  
     
//      for(Integer entity_dofs_iter=0;entity_dofs_iter<dofs_per_entity;entity_dofs_iter++)
//         for(Integer fs_dim=0;fs_dim<n_components;fs_dim++)
//            {           
//             dofmap_vec[elem_id][loc_dof_count]=old_dof;
//             // // std::cout<<old_dof<<std::endl;
//              // // std::cout<<"entity_dofs_iter = "<<entity_dofs_iter<<std::endl;
//              // // std::cout<<"old global_dof_count = "<<global_dof_count<<std::endl;
//             // // std::cout<< "(dofmap_vec)[elem_id][loc_dof_count]="<<dofmap_vec[elem_id][loc_dof_count]<<std::endl;
//              loc_dof_count++;
//              old_dof++;
//            }
//       // // std::cout<<"global_dof_count = "<<global_dof_count<<std::endl;
//     }
 
//     }

//     // // std::cout<<std::endl;
//     // for(int s=0;s<dofmap_vec[elem_id].size();s++)
//       // // std::cout<<dofmap_vec[elem_id][s]<<" ";
//     // // std::cout<<std::endl;
//     // // std::cout<<"loc_dof_count="<<loc_dof_count<<std::endl;
//     // // std::cout<<"elem_id="<<elem_id<<std::endl;

//   }


//    // // std::cout<<" DOFMAP VEC5 ==" <<std::endl;
//    // // std::cout<<" dofmap_vec.size() ==" <<dofmap_vec.size()<<std::endl;
 
//         //  for(int i=0;i<dofmap_vec.size();i++)
//         //   {
//         //     for(int j=0;j<dofmap_vec[i].size();j++)
//         //   // std::cout<<dofmap_vec[i][j]<<" ";
//         //   // std::cout<<std::endl;
//         // }
//   // // std::cout<<"after  dofmap_vec.size() ==" <<dofmap_vec.size()<<std::endl;
// }

// template<Integer N,typename FunctionSpace, typename FlagArray,typename Element,typename MeshT, Integer P, long M,typename OS, typename SpaceDofs>
// typename std::enable_if< -1==N,void>::type
// ElementDofMap_LoopEntities5(const ConnectivitySimpliacialMap<MeshT,P>& entitiestuple,
//               FlagArray& flagtuples,
//               Element& elem,
//               const Integer& elem_id,
//               std::vector<Array<Integer, M>>&  dofmap_vec,
//               Integer& global_dof_count,
//               Integer &loc_dof_count,
//               const OS &dofs_offset,
//               SpaceDofs& space_dofs,
//               Integer& dofs_count  )
// {
//   // // std::cout<<"ElementDofMap_LoopEntities5 K ==-1" <<std::endl;
// }






// template<Integer K=0,typename Elem,typename FunctionSpace,typename...FunctionSpaces, typename Element, typename...Args1,typename...Args2,typename...Args3,typename OS, typename SpaceDofs,typename ArrayNdofs>
// typename std::enable_if< 0<sizeof...(FunctionSpaces),void>::type
// ElementDofMap5(const std::tuple<Args1...>& entitiestuple,
//                    std::tuple<Args2...>& flagtuples,
//                    Element& elem,
//                    const Integer& elem_id,
//                    std::tuple<Args3...>& dofmap_vec,
//                    Integer& global_dof_count,
//                    Integer& loc_dof_count,
//                    OS &dofs_offset,
//                    SpaceDofs& space_dofs,
//                    ArrayNdofs& array_ndofs,
//                    const Integer level   )
// {





//     // // std::cout<<"__________qui0111111=2=" << K <<std::endl;

//  static constexpr Integer M=sizeof...(FunctionSpaces);
//  // const auto& entity_tuple=std::get<M>(entitiestuple);

//  const auto& entity_tuple=std::get<K>(entitiestuple);
//  // auto &flagtuple=std::get<M>(flagtuples);

//   auto &flagtuple=std::get<K>(flagtuples);
//  static constexpr Integer NN=std::tuple_size<remove_all_t<decltype(flagtuple)>>::value-1;
//  loc_dof_count=0;

//  ElementDofMap_LoopEntities5<NN,ElemFunctionSpace<Elem,FunctionSpace>>
//                             (entity_tuple,flagtuple,elem,elem_id,(*tuple_get<K>(dofmap_vec)),global_dof_count,loc_dof_count,std::get<K>(dofs_offset),space_dofs[K],array_ndofs[K][level]);

//  // // std::cout<<"__________qui==" << K <<std::endl;
//  auto & dm=*tuple_get<K>(dofmap_vec);
//  // // std::cout<<" elem dofmap 5 K==" << K <<std::endl;
//  // // std::cout<<" dm.size()==" << dm.size() <<std::endl;

// const auto& m3=entity_tuple.elem2dofs(elem_id);
//    // // std::cout<<"elem2dofs=="<<std::endl;
//     // for(int i=0;i<m3.size();i++)
//       // // std::cout<<m3[i]<<" "<<std::endl;
//    // // std::cout<<std::endl;

// //  for(int i=0;i<dm.size();i++)
// //   {
// //     for(int j=0;j<dm[i].size();j++)
// //   // std::cout<<dm[i][j]<<" ";
// //   // std::cout<<std::endl;

// // }

//   ElementDofMap5<K+1,Elem,FunctionSpaces...>(entitiestuple,flagtuples,elem,elem_id,dofmap_vec,global_dof_count,loc_dof_count,dofs_offset,space_dofs,array_ndofs,level);

// };




// template<Integer K=0,typename Elem,typename FunctionSpace,typename...FunctionSpaces, typename Element, typename...Args1,typename...Args2,typename...Args3,typename OS, typename SpaceDofs,typename ArrayNdofs>
// typename std::enable_if< 0==sizeof...(FunctionSpaces),void>::type
// ElementDofMap5(const std::tuple<Args1...>& entitiestuple,
//                    std::tuple<Args2...>& flagtuples,
//                    Element& elem,
//                    const Integer& elem_id,
//                    std::tuple<Args3...>& dofmap_vec,
//                    Integer& global_dof_count, 
//                    Integer& loc_dof_count,
//                    OS &dofs_offset,
//                    SpaceDofs& space_dofs,
//                    ArrayNdofs& array_ndofs,
//                    const Integer level )
// {


//   // // std::cout<<"__________qui00000=2=" << K <<std::endl;
//  // const auto& entity_tuple=std::get<0>(entitiestuple);
//   const auto& entity_tuple=std::get<K>(entitiestuple);
//  auto& flagtuple=std::get<K>(flagtuples);
//  static constexpr Integer NN=std::tuple_size<remove_all_t<decltype(flagtuple)>>::value-1;
//  loc_dof_count=0;


//  ElementDofMap_LoopEntities5<NN,ElemFunctionSpace<Elem,FunctionSpace>>
//                            (entity_tuple,flagtuple,elem,elem_id,(*tuple_get<K>(dofmap_vec)),global_dof_count,loc_dof_count,std::get<K>(dofs_offset),space_dofs[K],array_ndofs[K][level]);

//  // // std::cout<<"__________qui2=2=" << K <<std::endl;

//  auto & dm=*tuple_get<K>(dofmap_vec);
//  // // std::cout<<" elem dofmap 4 K==" << K <<std::endl;
//  // // std::cout<<" dm.size()==" << dm.size() <<std::endl;
//  const auto& m3=entity_tuple.elem2dofs(elem_id);

//    // // std::cout<<"elem2dofs=="<<std::endl;
//    //  for(int i=0;i<m3.size();i++)
//    //    // std::cout<<m3[i]<<" "<<std::endl;
//    // // std::cout<<std::endl;

// //  for(int i=0;i<dm.size();i++)
// //   {
// //     for(int j=0;j<dm[i].size();j++)
// //   // std::cout<<dm[i][j]<<" ";
// //   // std::cout<<std::endl;
// // }


// };








// template<typename FunctionSpace,  Integer Mmax, Integer M,  typename ...Args>
// typename std::enable_if< (M>Mmax),void >::type
// initialize_flag_tuple_entities_aux5(std::tuple< Args...> const &tuple, 
//                                std::array<std::vector< std::array<Integer,2> >, 
//                                FunctionSpace::entity.size() > & entity_found)
//     {static_assert(M>=0," the tuple must have non negative length ");};

// template<typename FunctionSpace,  Integer Mmax,Integer M, typename...Args>
// typename std::enable_if< (M<=Mmax),void >::type
// initialize_flag_tuple_entities_aux5(std::tuple< Args...> const &tuple,
//                                std::array<std::vector< std::array<Integer,2> >, FunctionSpace::entity.size() > & entity_found)
// {
//      static_assert(M>=0," the tuple must have non negative length ");
//      static constexpr Integer entity_dim=FunctionSpace::entity[M];
//      // maybe 
//      Integer entity_length=tuple_get<entity_dim>(tuple).size();

//      entity_found[M].resize(entity_length,{-1,-1});

//      // // std::cout<<"entity_found[M] size==" <<entity_found[M].size()<<std::endl;
//      initialize_flag_tuple_entities_aux5<FunctionSpace,Mmax,M+1>(tuple,entity_found);       
// };


// template<typename FunctionSpace, typename...Args>
//  void initialize_flag_tuple_entities5(
//   std::tuple< Args...> const &tuple,
//   std::array<std::vector< std::array<Integer,2> >, FunctionSpace::entity.size() > & entity_found)
// {
//      static constexpr Integer entities_nums=FunctionSpace::entity.size();
//      initialize_flag_tuple_entities_aux5<FunctionSpace,entities_nums-1,0>(tuple,entity_found);       
// };




// template<typename Elem,typename FunctionSpace,typename...FunctionSpaces>
// struct FlagTupleType5
// {
//      using FS=ElemFunctionSpace<Elem,FunctionSpace>;
//      static constexpr Integer FEFamily=FS::FEFamily;
//      static constexpr Integer Order=FS::Order; 
//      static constexpr Integer entities_nums=FS::entity.size();
//      using rest = typename FlagTupleType5<Elem,FunctionSpaces...>::type;
//      using ens  = typename std::array<std::vector<std::array<Integer,2>>, entities_nums>;
//      using type = decltype( std::tuple_cat( std::declval< rest >(), 
//                                             std::declval< std::tuple<ens> >() ) );
// };


// template<typename Elem,typename FunctionSpace>
// struct FlagTupleType5<Elem,FunctionSpace>
// {
//  using FS=ElemFunctionSpace<Elem,FunctionSpace>;
//  static constexpr Integer FEFamily=FS::FEFamily;
//  static constexpr Integer Order=FS::Order; 
//  static constexpr Integer entities_nums=FS::entity.size();
//  using ens =  typename std::array<std::vector<std::array<Integer,2>>, entities_nums>;
//  using type = typename std::tuple<ens>;
// };


// template<typename Elem,typename BaseFunctionSpace,typename...BaseFunctionSpaces, typename MeshT>
// typename std::enable_if< 0<sizeof...(BaseFunctionSpaces), 
//                          typename FlagTupleType4<Elem,BaseFunctionSpace,BaseFunctionSpaces...>::type>::type
// FlagTuple5(ConnectivitySimpliacialMapCollection<MeshT>& entities)
// {    
//       using FunctionSpace=ElemFunctionSpace<Elem,BaseFunctionSpace>;
//       static constexpr Integer entities_nums=FunctionSpace::entity.size();
//       using type = typename std::array<std::vector<std::array<Integer,2>>, entities_nums>;
//       type ens;

//       initialize_flag_tuple_entities5<FunctionSpace>(entities.tuple_entities(),ens);
//       return std::tuple_cat(FlagTuple5<Elem,BaseFunctionSpaces...>(entities),
//                             std::tuple<type>(ens));
// }

// template<typename Elem,typename BaseFunctionSpace, typename...BaseFunctionSpaces, typename MeshT>
//  typename std::enable_if< 0==sizeof...(BaseFunctionSpaces), 
//                           typename FlagTupleType4<Elem,BaseFunctionSpace>::type>::type
// FlagTuple5(ConnectivitySimpliacialMapCollection<MeshT>& entities)
// {   
//       using FunctionSpace=ElemFunctionSpace<Elem,BaseFunctionSpace>;
//       static constexpr Integer entities_nums=FunctionSpace::entity.size();
//       using type = typename std::array<std::vector<std::array<Integer,2>>, entities_nums>;
//       type ens;
//       initialize_flag_tuple_entities5<FunctionSpace>(entities.tuple_entities(),ens);
//       return std::tuple<type>(ens);
// }










  // template <Integer N>
  // struct FlagEntitiesTuple;

  // template<Integer N>
  // struct FlagSimplexEntitiesTuple
  //  {
  //   using type=std::map<std::array<Integer,N+1>,std::shared_ptr<std::array<Integer,2>>>;
  //  };


  // template<typename FunctionSpace,Integer Nmax,Integer N>
  // struct FlagSimplexEntitiesFunctionSpaceTupleHelper;


  // template<typename FunctionSpace,Integer Nmax>
  // struct FlagSimplexEntitiesFunctionSpaceTupleHelper<FunctionSpace,Nmax,Nmax>
  //  {
  //   using single_type=typename FlagSimplexEntitiesTuple<FunctionSpace::entity[Nmax]>::type;
  //   using type=std::tuple<single_type>;
  //  };



  // template<typename FunctionSpace,Integer Nmax,Integer N>
  // struct FlagSimplexEntitiesFunctionSpaceTupleHelper
  //  {
  //    using single_type=typename FlagSimplexEntitiesTuple<FunctionSpace::entity[N]>::type;
  //    using rest=typename FlagSimplexEntitiesFunctionSpaceTupleHelper<FunctionSpace,Nmax,N+1>::type;
  //    using type = decltype( std::tuple_cat( std::declval< std::tuple<single_type> >() ,
  //                                           std::declval< rest >() 
  //                                           ) );

  //  };


  // template<typename FunctionSpace>
  // struct FlagSimplexEntitiesFunctionSpaceTuple
  //  {
  //   using type=typename FlagSimplexEntitiesFunctionSpaceTupleHelper<FunctionSpace,FunctionSpace::entity.size()-1,0>::type;
  //  };
 //  template <Integer N>
 //  struct FlagEntitiesTuple
 //  {
 //   using rest = typename FlagEntitiesTuple<N-1>::type; 
 //   using single_type=std::map<std::array<Integer,N+1>,Integer>;
 //   using tuple_single_type=std::tuple<single_type>;
 //   using type = decltype( std::tuple_cat(std::declval< rest >() , std::declval< tuple_single_type >()) );
 //  };

 //  template <Integer N>
 // using FlagEntities=typename FlagEntitiesTuple<N>::type;



// template<typename Elem,typename BaseFunctionSpace,typename...BaseFunctionSpaces>
// struct FlagTupleType6
// {
//      using FunctionSpace=ElemFunctionSpace<Elem,BaseFunctionSpace>;
//      using rest = typename FlagTupleType6<Elem,BaseFunctionSpaces...>::type;
//      using ens  = typename FlagSimplexEntitiesFunctionSpaceTuple<FunctionSpace>::type;
//      using type = decltype( std::tuple_cat( std::declval< std::tuple<ens> >(),std::declval< rest >()
//                                              ) );
// };


// template<typename Elem,typename BaseFunctionSpace>
// struct FlagTupleType6<Elem,BaseFunctionSpace>
// {
//  using FunctionSpace=ElemFunctionSpace<Elem,BaseFunctionSpace>;
//  using single_type = typename FlagSimplexEntitiesFunctionSpaceTuple<FunctionSpace>::type;
//  using type = typename std::tuple<single_type>;
// };


// template<typename Elem,typename...BaseFunctionSpaces>
// using FlagTuple6=typename FlagTupleType6<Elem,BaseFunctionSpaces...>::type;






















template<typename BaseFunctionSpace, typename...BaseFunctionSpaces, typename MeshT, typename Dofmap, typename OffSetT,typename T1,typename T2,typename T3>//, typename DofMapT, typename OffSetT>
void dofmap_fespace5( 
             MeshT& mesh,
             Bisection<MeshT>& bisection,
             Node2ElemMap<MeshT>& node2elem,
             Dofmap& dofsdm_,
             OffSetT& dofs_offset_arr,
             Integer& n_elements,//Array<std::vector<Integer>,1 + sizeof...(BaseFunctionSpaces)>& level_array_ndofs
             T1& tuple_map,
             T2& tuple_parent_map,
             T3& tuple_map_dofs,
             bool update=false
             )
{

    using     Elem = typename MeshT::Elem; 
    constexpr auto Dim=MeshT::Dim;
    constexpr auto ManifoldDim=MeshT::ManifoldDim;

    constexpr auto dofs_per_elem=DofsPerElemNums<Elem,BaseFunctionSpace,BaseFunctionSpaces...>::value;


    // this shuld be passed from outside
    // FlagTuple6<Elem,BaseFunctionSpace,BaseFunctionSpaces...> flag_entities;
    // auto n_levels=entities.n_levels();
    // auto mesh=entities.mesh();
    // auto& space_dofs=dofsdm_.space_dofs();
    auto& dofmap_vec=dofsdm_.dofmap();


    // // std::cout<<"dofs dm reserve begin"<<std::endl;
    dofsdm_.init_dofmap(mesh.n_elements());
    // // std::cout<<"dofs dm reserve end"<<std::endl;
    // compute the connection node to elem (auxiliary tool, often used)
    clock_t begin = clock();
    
    // auto& entities_collection=entities.tuple_entities();
    // init<  BaseFunctionSpace,BaseFunctionSpaces...>(entities);
//     NodeToElem<MeshT> node_2_elem(mesh);
//     const auto& node2elem=node_2_elem.val();
//     const auto& n_elements=mesh.n_elements();
    const auto dofs_offset=OffsetDofs<0,Elem,BaseFunctionSpace,BaseFunctionSpaces...>();
//     const auto n_spaces=1+sizeof...(BaseFunctionSpaces); 
//           auto entitiestuple=EntitiesOfFunctionSpaceTuple<Elem,BaseFunctionSpace,BaseFunctionSpaces...>(mesh,node2elem);
          // auto flagtuples= FlagTuple5<Elem,BaseFunctionSpace,BaseFunctionSpaces...>(entities);
    
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    // // std::cout<<"dofmap entities tuple elapsed_secs="<<elapsed_secs<<std::endl;
    // resize_tuple_of_ptr_vector(dofmap_vec,n_elements); 
    


    auto& tracker=bisection.tracker();
    auto n_levels= tracker.current_iterate();


    std::vector<Integer> global_dof_count;
    if(update==false)  
    {
    // // std::cout<<"false"<<std::endl;
    global_dof_count.push_back(0);
    } 
    else
    {
      // // std::cout<<"true"<<std::endl;
    // for(Integer i=0;i<dofsdm_.level_n_dofs_array().size();i++)
    // {
    //   for(Integer j=0;j<dofsdm_.level_n_dofs_array()[i].size();j++)
    //     // std::cout<<dofsdm_.level_n_dofs_array()[i][j]<<" ";
    //   // std::cout<<std::endl;
    //     // space_dofs[space_id]->resize(space_components[space_id][3]);
    // } 

    // // std::cout<<"global_dof_count"<<std::endl;
    for(Integer j=0;j<dofsdm_.level_n_dofs_array()[0].size();j++)
    {
        global_dof_count.push_back(dofsdm_.level_n_dofs_array()[0][j]);
        // // std::cout<<global_dof_count[j]<<" ";
        // space_dofs[space_id]->resize(space_components[space_id][3]);
    } 

    // std::cout<<std::endl;
    }
    // loop on all the elements
     for(Integer j=0;j<global_dof_count.size();j++)
    {
        // // std::cout<<global_dof_count[j]<<" ";
        // space_dofs[space_id]->resize(space_components[space_id][3]);
    }  


    // for(std::size_t i=0;i<level_array_ndofs.size();i++)
    //   {
    //     level_array_ndofs[i].resize(n_levels,0);
    //   //   for(std::size_t j=0;j<n_levels;j++)
    //   // level_array_ndofs[i][j]=0;
    //   }
    // // std::cout<<"dofmap_fespace5 2="<<elapsed_secs<<std::endl;
    // auto& bisection=entities.bisection();

    

  // for(Integer ll=0;ll<n_levels;ll++)
  //   {
  //       if(global_dof_count.size()<ll)
  //       global_dof_count[ll]=dofsdm_.level_n_dofs_array()[0][j];
  //       // space_dofs[space_id]->resize(space_components[space_id][3]);
  //   }



    // auto& node2elem=entities.node2elem();
    // auto level=entities.level();
    begin = clock();

   DofMapFESpace<MeshT,BaseFunctionSpace,BaseFunctionSpaces...> dofmap(mesh,node2elem,bisection);
   

   std::vector<std::vector<Integer>> elem2dofs;
   // TupleOfTupleMapConstructor<std::shared_ptr<std::vector<Integer>>,Elem,BaseFunctionSpace,BaseFunctionSpaces...> tuple_map;
   // TupleOfTupleMapConstructor<bool,Elem,BaseFunctionSpace,BaseFunctionSpaces...> tuple_parent_map;
   // TupleOfTupleMapConstructor<std::shared_ptr<Integer>,Elem,BaseFunctionSpace,BaseFunctionSpaces...> tuple_map_dofs;

    // // std::cout<<"levels="<<n_levels<<std::endl;
    // // std::cout<<"PRE LEVEL N DOFS ARRAY ="<<std::endl;
    // loop on the spaces
    // for(std::size_t i=0; i<dofsdm_.n_dofs().size() ;i++)
    //   { dofsdm_.level_n_dofs_array()[i].resize(n_levels);
    //     // loop on the levels
    //     for(std::size_t j=0;j<n_levels;j++)
    //        {dofsdm_.level_n_dofs_array()[i][j]=global_dof_count[j];
    //         // std::cout<<dofsdm_.level_n_dofs_array()[i][j]<<" ";}
    //      // std::cout<<std::endl;
    //    }
    // // std::cout<<"dofmap_fespace5 3="<<elapsed_secs<<std::endl;
    // Integer n_levels=bisection.tracker().current_iterate()+1;
    // decltype(flag_entities) ffe(5);
    Integer level=0;
    bool has_tracked=tracker.has_tracked();
    for(Integer elem_iter=n_elements;elem_iter<mesh.n_elements();elem_iter++)
    {


      // if(!elem_belongs_to_level(mesh,elem_iter,level,bisection)) continue;
       
       // change it for smarter algorithms   
      // // std::cout<<" dofmap elem == "<< elem_iter <<" /"<< mesh.n_elements()<< std::endl;
      // std::cout<<"dofmap5 elem_iter="<<elem_iter<<std::endl;
      // for(Integer i=0;i<global_dof_count.size();i++)
      //   std::cout<<"global_dof_count="<<global_dof_count[i]<<std::endl;

       auto& elem=mesh.elem(elem_iter);
       auto &elem_id=elem.id;
       // // std::cout<<"has_tracked="<<has_tracked<<std::endl;
       if(has_tracked)
          level=tracker.get_iterate(elem_id);
       // // std::cout<<"level="<<level<<std::endl;
       // // std::cout<<"0<level && global_dof_count.size()-1<level="<<(0<level && global_dof_count.size()-1<level)<<std::endl;
       
      if(0<level && global_dof_count.size()-1<level)
         global_dof_count.push_back(global_dof_count[global_dof_count.size()-1]);

 
       Integer loc_dof_count=0;
       // // std::cout<<"DOFMAPSPACE5 ELEM ITER="<<elem_iter<<" level=="<<level<<std::endl;
       if(level>-1)
       dofmap.init_elem(dofmap_vec,tuple_map,tuple_parent_map,tuple_map_dofs,
                        global_dof_count,elem,level,mesh,tracker);

       // ElementDofMap5<0,Elem,BaseFunctionSpace,BaseFunctionSpaces...>
       //                   (entities.tuple_entities(),flag_entities,elem,elem_id,dofmap_vec,global_dof_count[level],loc_dof_count,dofs_offset,space_dofs,level_array_ndofs,level);

 
      // }

    }
    n_elements=mesh.n_elements();
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    // // std::cout<<"dofmap tuple elapsed_secs="<<elapsed_secs<<std::endl;
    
  //   for(std::size_t i=0; i<dofsdm_.n_dofs().size() ;i++)
  //     { dofsdm_.level_n_dofs_array()[i].resize(n_levels);
  //       for(std::size_t j=0;j<n_levels;j++)
  //     dofsdm_.level_n_dofs_array()[i][j]=global_dof_count[j];
  // }


  // auto dm=(*tuple_get<0>(dofmap_vec));
  // for(std::size_t i=0;i<dm.size();i++)
  //   {for(std::size_t j=0;j<dm[i].size();j++)
  //   // std::cout<< dm[i][j]<<" ";
  //    // std::cout<<std::endl;}
// // std::cout<<" level_array_ndofs[i][j]"<<std::endl;

//   for(std::size_t i=0;i<level_array_ndofs.size();i++)
//     {for(std::size_t j=0;j<level_array_ndofs[i].size();j++)
//     // std::cout<< level_array_ndofs[i][j]<<" ";
//      // std::cout<<std::endl;}

// // std::cout<<" global_dof_count"<<std::endl;

//   for(std::size_t i=0;i<global_dof_count.size();i++)
//     {
//     // std::cout<< global_dof_count[i]<<" ";
//      // std::cout<<std::endl;}
// 
    // // std::cout<<"AFTER LEVEL N DOFS ARRAY ="<<std::endl;
    // loop on the spaces
    for(std::size_t i=0; i<dofsdm_.n_dofs().size() ;i++)
      { dofsdm_.level_n_dofs_array()[i].resize(n_levels);
        // loop on the levels
        for(std::size_t j=0;j<n_levels;j++)
           {
            dofsdm_.level_n_dofs_array()[i][j]=global_dof_count[j];
            // // std::cout<<dofsdm_.level_n_dofs_array()[i][j]<<" ";
          }
         // // std::cout<<std::endl;
       }

    // // std::cout<<"AFTER LEVEL N DOFS ="<<std::endl;

        for(Integer i=0;i<dofsdm_.n_dofs().size();i++)
        {
          dofsdm_.n_dofs()[i]=dofsdm_.level_n_dofs_array()[i][n_levels-1];
          // // std::cout<<dofsdm_.n_dofs()[i]<<std::endl;
        }


 // // std::cout<<" after dofmap tuple elapsed_secs="<<elapsed_secs<<std::endl;
// tuple2array_dofs_offset<n_spaces-1>(dofs_offset,dofs_offset_arr);
}





  template<typename Elem>
  class FiniteElem;

  template<typename FullFunctionSpace,Integer EntityDim>
  class Entity2Dofs
  {
  public:
     using FunctionSpace=FullFunctionSpace;
     using TupleOfSpaces=typename FunctionSpace::FunctionSpace::TupleOfSpaces;
     using ElemDofMap = typename FunctionSpace::DofsDM::ElemDofMap;
     using Elem = typename FullFunctionSpace::Elem;     
     static constexpr Integer ManifoldDim=Elem::ManifoldDim;
     static constexpr Integer manifold_points=ElemEntityNPoints<Elem,ManifoldDim>::value;
     static constexpr Integer entity_points=ElemEntityNPoints<Elem,EntityDim>::value;
     static constexpr Integer combinations_nums=ElemEntityCombinations<Elem,EntityDim>::value;
     static constexpr auto NLocalDofs=FunctionSpace::DofsDM::NLocalDofs;
     using ArrEntityPoints=std::array<Integer,entity_points>;
     using Map=std::map<ArrEntityPoints, std::vector<std::vector<Integer>> >;
     using MapBoolPtr=std::map<std::array<Integer,entity_points+1>,bool>;
     using Vec=std::vector<Integer>;
     using VecVec=std::vector<Vec>;
     using VecVecVec=std::vector<VecVec>;
     using SortedVec=std::vector<std::vector<Integer>>;
     using SortedLevelsVec=std::vector<SortedVec>;



    Entity2Dofs(const std::shared_ptr<FullFunctionSpace> W_ptr):
    spaces_ptr_(W_ptr)
    {

    }



    template<bool Dofs2Elements,typename Space,Integer Nmax,Integer N=0,typename Tracker, typename DofsDM >
    inline std::enable_if_t<(N>=Nmax),void> 
    entity_loop(Integer& cont, Tracker& tracker, DofsDM& dm, const std::array<Integer,entity_points>& entity_nodes1, const Elem& elem, const FiniteElem<Elem>&FE)
    {}


    template<bool Dofs2Elements,typename Space,Integer Nmax,Integer N=0,typename Tracker, typename DofsDM >
    inline std::enable_if_t<(N<Nmax),void> 
    entity_loop(Integer& cont,Tracker& tracker, DofsDM& dm, const std::array<Integer,entity_points>& entity_nodes1, const Elem& elem, const FiniteElem<Elem>&FE)
    {

      // // std::cout<<" entity dim=="<<Space::entity[N]<<std::endl;

      constexpr Integer EntityDim2=Space::entity[N];

      constexpr Integer entity_points2=ElemEntityNPoints<Elem,EntityDim2>::value;

      constexpr Integer combinations_nums2=ElemEntityCombinations<Elem,EntityDim2>::value;



      std::array<Integer,entity_points2> entity_nodes2;

      Integer entity2[entity_points2];

      auto& level=FE.level();
      auto n_levels=tracker.current_iterate();


      const auto& nodes=elem.nodes;
      Integer cont_nodes;
      bool found;

      // // std::cout<<"nodes="<<std::endl;

      //        for(std::size_t i=0;i<nodes.size();i++)
      //       {
      //         // std::cout<<nodes[i]<<" ";
      //       }   
      // // std::cout<<" "<<std::endl;

      for(Integer entity_iter=0;entity_iter<combinations_nums2;entity_iter++)
         { 
          found=false;
          // // std::cout<<"entity_iter="<<entity_iter<<std::endl;
          Combinations<manifold_points, entity_points2>::generate(entity_iter,entity2);
          // for(std::size_t i=0;i<entity_points2;i++)
          //   {
          //     // std::cout<<entity2[i]<<" ";
          //   }

          // // std::cout<<" "<<std::endl;


          for(std::size_t i=0;i<entity_points2;i++)
            {
              entity_nodes2[i]=nodes[entity2[i]];
              // // std::cout<<entity_nodes2[i]<<" ";
            }
          // // std::cout<<" "<<std::endl;

          std::sort(entity_nodes2.begin(),entity_nodes2.end());


          // check if all the nodes entity_nodes1 belong to entity_nodes2
          cont_nodes=0;
          for(Integer i=0;i<entity_points;i++)
          {
            for(Integer j=0;j<entity_points2;j++)
            {
             if(entity_nodes1[i]==entity_nodes2[j])
             {
              cont_nodes++;
              break;
             }
            }
            if(cont_nodes==entity_points || Dofs2Elements)
              {
               // add dofmap
                // // std::cout<<"found on level "<<level<<std::endl;
                found=true;
                auto& dofs=map_[entity_nodes1];

                if(dofs.size()==0)
                {
                  dofs.resize(n_levels);
                }
                else if(dofs.size()<n_levels)
                {
                  for(Integer i=dofs.size();i<n_levels;i++)
                      dofs.push_back(std::vector<Integer>{});
                }

                
                
                // // std::cout<<" entity_nodes1 "<<std::endl;

                // for(Integer m=0;m<entity_nodes1.size();m++)
                // {
                //   // std::cout<<entity_nodes1[m]<<" ";
                // }
                // // std::cout<<" "<<std::endl;
                // // std::cout<<"  map_ "<<std::endl;

                for(Integer m=0;m<Space::dofs_per_entity[N];m++)
                  for(Integer k=0;k<Space::NComponents;k++)
                     {

                      dofs[level].push_back(dm[cont]);
                      // // std::cout<<dm[cont]<<" ";
                      cont++;
                      }
                // // std::cout<<" "<<std::endl;

                // std::sort(dofs.begin(),dofs.end());
                // dofs.erase( std::unique( dofs.begin(), dofs.end() ), dofs.end() );

                // // std::cout<<"dofs "<<std::endl;

                // for(Integer m=0;m< dofs[level].size();m++)
                // {
                //   // std::cout<<dofs[level][m]<<" "<<std::endl;
                // }
                // // std::cout<<" "<<std::endl;

              }
          }


          if(found==false)
          {
            cont=cont+Space::NComponents*Space::dofs_per_entity[N];
          }


       }


      entity_loop<Dofs2Elements,Space,Nmax,N+1>(cont,tracker,dm,entity_nodes1,elem,FE);

    }




    template<bool Dofs2Elements,Integer N,Integer Nmax=TupleTypeSize<TupleOfSpaces>::value,typename Tracker, typename DofsDM >
    inline std::enable_if_t<(N>=Nmax),void> 
    loop_all_spaces(Tracker& tracker, DofsDM& dm,const std::array<Integer,entity_points>& entity_nodes, const Elem& elem, const FiniteElem<Elem>&FE)
    {}

    template<bool Dofs2Elements, Integer N,Integer Nmax=TupleTypeSize<TupleOfSpaces>::value,typename Tracker, typename DofsDM >
    inline std::enable_if_t<(N<Nmax),void> 
    loop_all_spaces(Tracker& tracker, DofsDM& dm, const std::array<Integer,entity_points>& entity_nodes, const Elem& elem, const FiniteElem<Elem>&FE)
    {
      // // std::cout<<" space=="<<N<<std::endl;
      using Space=GetType<TupleOfSpaces,N>;
      constexpr auto entity=Space::entity;
      auto& single_dm=tuple_get<N>(elemdofmap_);
      dm.template dofmap_get<N>(single_dm,FE.elem_id(),FE.level());
      // // std::cout<<" single_dm=="<<single_dm<<std::endl;

      Integer cont=0;
      entity_loop<Dofs2Elements,Space,entity.size()>(cont,tracker,single_dm,entity_nodes,elem,FE);
      loop_all_spaces<Dofs2Elements,N+1>(tracker,dm,entity_nodes,elem,FE);
    }



    template<bool Dofs2Elements,Integer N,Integer...Ns,typename Tracker, typename DofsDM >
    inline std::enable_if_t<(sizeof...(Ns)==0),void> 
    loop_some_spaces(Tracker& tracker, DofsDM& dm,const std::array<Integer,entity_points>& entity_nodes, const Elem& elem, const FiniteElem<Elem>&FE)
    {
      using Space=GetType<TupleOfSpaces,N>;
      constexpr auto entity=Space::entity;
      auto& single_dm=tuple_get<N>(elemdofmap_);
      dm.template dofmap_get<N>(single_dm,FE.elem_id(),FE.level());
      Integer cont=0;
      entity_loop<Dofs2Elements,Space,entity.size()>(cont,tracker,single_dm,entity_nodes,elem,FE);

    }

    template<bool Dofs2Elements,Integer N,Integer...Ns,typename Tracker, typename DofsDM >
    inline std::enable_if_t<(sizeof...(Ns)>0),void> 
    loop_some_spaces(Tracker& tracker, DofsDM& dm, const std::array<Integer,entity_points>& entity_nodes, const Elem& elem, const FiniteElem<Elem>&FE)
    {
      using Space=GetType<TupleOfSpaces,N>;
      constexpr auto entity=Space::entity;
      auto& single_dm=tuple_get<N>(elemdofmap_);
      dm.template dofmap_get<N>(single_dm,FE.elem_id(),FE.level());

      Integer cont=0;
      entity_loop<Dofs2Elements,Space,entity.size()>(cont,tracker,single_dm,entity_nodes,elem,FE);
      loop_some_spaces<Dofs2Elements,Ns...>(tracker,dm,entity_nodes,elem,FE);
    }


    template<bool Dofs2Elements,Integer M,Integer... Ms,Integer Nmax=TupleTypeSize<TupleOfSpaces>::value,typename Tracker, typename DofsDM >
    std::enable_if_t<(M==-1),void>
    loop_spaces_aux(Tracker& tracker, DofsDM& dofsdofmap, const std::array<Integer,entity_points>& entity_nodes, const Elem& elem, const FiniteElem<Elem>&FE)
    {
      loop_all_spaces<Dofs2Elements,0>(tracker,dofsdofmap,entity_nodes,elem,FE);
    }
      template<bool Dofs2Elements,Integer M,Integer... Ms,Integer Nmax=TupleTypeSize<TupleOfSpaces>::value,typename Tracker, typename DofsDM >
    std::enable_if_t<(M!=-1),void>
    loop_spaces_aux(Tracker& tracker, DofsDM& dofsdofmap, const std::array<Integer,entity_points>& entity_nodes, const Elem& elem, const FiniteElem<Elem>&FE)
    {
      loop_some_spaces<Dofs2Elements,M,Ms...>(tracker,dofsdofmap,entity_nodes,elem,FE);
    }


    template<bool Dofs2Elements, Integer M=-1, Integer...Ms>
    void build()
    {
     Integer entity_[entity_points];
     std::array<Integer,entity_points> entity_nodes;




     auto &dofsdofmap=spaces_ptr_->dofsdofmap();
     auto &mesh=spaces_ptr_->mesh();
     auto& bisection=spaces_ptr_->bisection();
     auto& tracker=bisection.tracker();

     FiniteElem<Elem> FE(mesh);

     // // std::cout<<" Entity2Dofs"<<std::endl;
     for(Integer el=0;el<mesh.n_elements();el++)
     {

          // if(tracker.get_iterate(el)<0)continue;
          if(tracker.get_level(el)<0)continue;
          auto& elem=mesh.elem(el);
          const auto& nodes=elem.nodes;


          // // std::cout<<" el=="<<el<<std::endl;
          // FE.init(el,tracker.get_iterate(el));
          FE.init(el,tracker.get_level(el));

          for(Integer entity_iter=0;entity_iter<combinations_nums;entity_iter++)
             { 

              Combinations<manifold_points, entity_points>::generate(entity_iter,entity_);

              for(std::size_t i=0;i<entity_points;i++)
                entity_nodes[i]=nodes[entity_[i]];

              std::sort(entity_nodes.begin(),entity_nodes.end());


              loop_spaces_aux<Dofs2Elements,M,Ms...>(tracker,dofsdofmap,entity_nodes,elem,FE);
             }
         }


          for (auto it=map_.begin(); it!=map_.end(); ++it)
            {
              auto& second=it->second;
              for(Integer i=0;i<second.size();i++)
              {
                auto& dofs=second[i];
                std::sort(dofs.begin(),dofs.end());
                dofs.erase( std::unique( dofs.begin(), dofs.end() ), dofs.end() );
              }
            }

    //  // // std::cout<<"-----dofmap of space="<<0<<std::endl;
    // GetType<ElemDofMap,0> elemdm0;
    // Integer cont1=0;

    // for(Integer i=0;i<dofsdofmap.template dofmap_size<0>();i++)
    // {

    //     dofsdofmap.template dofmap_get<0>(elemdm0,i);
    //          // std::cout<<"dofmap == "<< i<<std::endl;
    //          auto& elem=mesh.elem(i);
    //          auto& nodes=elem.nodes;
    //        dofsdofmap.template dofmap_get<0>(elemdm0,i);
    //          // std::cout<<"dofmap == "<< i<<std::endl;
    //          // // std::cout<<elemdm<<std::endl;

    //   }
    }


    void print()
    {
          for (auto it=map_.begin(); it!=map_.end(); ++it)
            {
              auto& first=it->first;
              auto& second=it->second;

              // std::cout<<"NODES "<<std::endl;

              for(Integer i=0;i<first.size();i++)
              {
                  // std::cout<<first[i]<<" ";
              }
              // std::cout<<std::endl;
              // std::cout<<"DOFS "<<std::endl;
              for(Integer i=0;i<second.size();i++)
              {
                // std::cout<<"level =="<<i<<std::endl;
                for(Integer j=0;j<second[i].size();j++)
                {
                  // std::cout<<second[i][j]<<" ";
                }
                // std::cout<<std::endl;
              }
              // std::cout<<std::endl;
            }

    }


    auto& get()
    {

      auto& mesh=spaces_ptr_->mesh();
      auto& bisection=spaces_ptr_->bisection();
      auto& tracker=bisection.tracker();
      auto& n2e=spaces_ptr_->node2elem();

      Integer n_levels=tracker.current_iterate();
      if(map_vec_.size()==0)
        map_vec_.resize(n_levels);

      if(map_vec_.size()<n_levels)
      {
        for(Integer i=map_vec_.size();i<n_levels;i++)
        map_vec_.push_back(VecVec{});
      }


      for(Integer i=0;i<n_levels;i++)
      {
        map_vec_[i]=get(i);
      }
      return map_vec_;

    }

    auto& get(const Integer level)
    {

      auto& mesh=spaces_ptr_->mesh();
      auto& bisection=spaces_ptr_->bisection();
      auto& tracker=bisection.tracker();
      auto& n2e=spaces_ptr_->node2elem();

      if(map_vec_.size()==0)
        map_vec_.resize(tracker.current_iterate());

      auto& vec=map_vec_[level];
    
      Integer max_nodes=n2e.max_n_nodes();
      // vec.reserve(mesh.n_elements()*NLocalDofs*max_nodes);

          for (auto it=map_.begin(); it!=map_.end(); ++it)
              {
              // for(Integer i=0;i<it->first.size();i++)
              // {
              //     // std::cout<<it->first[i]<<" ";
              // }              
              //   // std::cout<< std::endl;
                vec.push_back(it->second[level]);
              }
      return vec;

    }

    auto& get(const std::vector<Integer>& levels)
    {
      levels_vec_.clear();
      Integer size=levels.size();
      levels_vec_.resize(size);
      for(Integer i=0;i<size;i++)
        levels_vec_[i]=get(levels[i]);


      return levels_vec_;

    }


      inline void ordered_entities_aux_aux(const Integer el,Integer& cont_elem,std::vector<bool>& elems_iterate)
      {
          auto& mesh=spaces_ptr_->mesh();
          auto& bisection=spaces_ptr_->bisection();
          auto& tracker=bisection.tracker();
          auto& n2e=spaces_ptr_->node2elem();
          auto& elem=mesh.elem(el);
          FiniteElem<Elem> FE(mesh);
          Integer entity_[entity_points];
          std::array<Integer,entity_points> entity_nodes;
          std::array<Integer,entity_points+1> entity_nodes_and_level;
          // // std::pair<std::array<Integer,entity_points>,Integer> pair;

          const auto& nodes=elem.nodes;
          auto level=tracker.get_iterate(el);
          FE.init(el,level);

          // // std::cout<< "elem="<<el <<std::endl;

          //// iterate on all the entities of the element
          for(Integer entity_iter=0;entity_iter<combinations_nums;entity_iter++)
             { 
              // generate the entity_iter entity in entity_;
              Combinations<manifold_points, entity_points>::generate(entity_iter,entity_);
              for(std::size_t i=0;i<entity_points;i++)
                {
                  entity_nodes[i]=nodes[entity_[i]];
                  // // std::cout<< entity_nodes[i]<<"";
                }


                // // std::cout<< std::endl;

              // sort the nodes of found entity
              std::sort(entity_nodes.begin(),entity_nodes.end());

              for(std::size_t i=0;i<entity_points;i++)
                {
                  entity_nodes_and_level[i]=entity_nodes[i];
                  // // std::cout<< entity_nodes[i]<<"";
                }
              entity_nodes_and_level[entity_points]=level;



          //     // here do something for the sorted dofmap
          //     // std::cout<<"pre pair" <<std::endl;
          //      // pair=std::make_pair(entity_nodes,level);
          //      // std::cout<<"after pre entity_nodes_and_level" <<std::endl;
          //      for(Integer i=0;i<entity_nodes_and_level.size();i++)
          //       // std::cout<<entity_nodes_and_level[i]<<" ";
          //     // std::cout<<" "<<std::endl;
          //     // std::cout<<"level="<<level<<std::endl;
          //     // if ( map_bool_.find(entity_nodes_and_level) == map_bool_.end() ) {
          //     //      // std::cout<<"not found="<<std::endl;
          //     //     } else {
          //     //      // std::cout<<" found="<<std::endl;
          //     //     }


          //        // std::cout<<std::endl;
          //      for(auto it=map_bool_.begin(); it!=map_bool_.end(); ++it)
          //      {
          //        for(Integer i=0;i<entity_nodes_and_level.size();i++)
          //            // std::cout<<it->first[i]<<" ";
          //          // std::cout<<std::endl;
              
          //      }
          //      // std::cout<<"end pre pair2" <<std::endl;
          //      auto m=map_bool_[entity_nodes_and_level];

          //      // std::cout<<"after pre pair2" <<std::endl;

              if(map_bool_[entity_nodes_and_level])
              {
                map_bool_[entity_nodes_and_level]=true;
                entities_ordered_[level].push_back(map_[entity_nodes][level]);

              }
            }
      elems_iterate[el]=true;
      cont_elem++;

    }









    inline bool find_next_elem(Integer& next_elem,const Integer&el,std::vector<bool>& elems_iterate)
    {
        auto& mesh=spaces_ptr_->mesh();
        auto& bisection=spaces_ptr_->bisection();
        auto& tracker=bisection.tracker();
        auto& n2e=spaces_ptr_->node2elem();
        auto& elem=mesh.elem(el);
        auto& nodes=elem.nodes;  
        auto level=tracker.get_iterate(el);
        bool found=false;
        for(Integer n=0;n<nodes.size();n++)
        {
          auto node2elems=n2e.get(nodes[n],level);
          for(Integer l=0;l<node2elems.size();l++)
          {
            if(elems_iterate[node2elems[l]]==false)
             {
              next_elem=node2elems[l];
              found=true;
              goto finish;

             } 
          }
        }
        finish:
        {}
        return found;
    }

    inline void ordered_entities_aux(Integer& next_elem,const Integer el,Integer& cont_elem,std::vector<bool>& elems_iterate)
    {


      // cont_elem++;// REMOVE <-----------------------------










        auto& mesh=spaces_ptr_->mesh();
        auto& bisection=spaces_ptr_->bisection();
        auto& tracker=bisection.tracker();
        auto& n2e=spaces_ptr_->node2elem();
        auto& elem=mesh.elem(el);
        auto& nodes=elem.nodes;  
        auto level=tracker.get_iterate(el);

        next_elem=-1;
        ordered_entities_aux_aux(el,cont_elem,elems_iterate);

              for(Integer n=0;n<nodes.size();n++)
              {
                auto node2elems=n2e.get(nodes[n],level);
                for(Integer l=0;l<node2elems.size();l++)
                {
                  if(elems_iterate[node2elems[l]]==false)
                   {
                    ordered_entities_aux_aux(node2elems[l],cont_elem,elems_iterate);

                   } 
                }
              }
              // for(Integer n=0;n<nodes.size();n++)
              // {
              //   auto node2elems=n2e.get(nodes[n],level);
              //   for(Integer l=0;l<node2elems.size();l++)
              //   {
              //     if(find_next_elem(next_elem,node2elems[l],elems_iterate))
              //       goto finish;
              //   }

                
              //  }
      // finish:
      // {}


 
    }
       




    inline void ordered_entities_aux(std::vector<bool>& elems_iterate ,const Integer level, const Integer level_n_elements)
    {
      auto& mesh=spaces_ptr_->mesh();
      auto& bisection=spaces_ptr_->bisection();
      auto& tracker=bisection.tracker();
      auto& n2e=spaces_ptr_->node2elem();
      auto n_levels=tracker.current_iterate();
      Integer cont_elem=0;
      Integer next_elem=0;
      Integer el=0;
      Integer n_elements=mesh.n_elements();

      // std::vector<bool> elems_iterate(level_n_elements,false);







      
      // loop on all the elements 

       while(cont_elem<level_n_elements)
       {

        if(!elem_belongs_to_level(mesh,el,level,tracker))
          { el++;
            continue;
          }

        ordered_entities_aux(next_elem,el,cont_elem,elems_iterate);
        if(next_elem<0)
        {
         for(Integer i=0;i<n_elements;i++)
          if(elems_iterate[i]==false)
          {
            el=i;
            break;
          }
        }
        else
        {
          el=next_elem;
        }

        // // std::cout<<"level="<<level<<", el="<<el<<", next_elem="<<next_elem<<", cont_elem="<<cont_elem<<"/"<<level_n_elements<<std::endl;


       }
      }
  

    inline auto& ordered_entities()
    {
      auto& mesh=spaces_ptr_->mesh();
      auto& bisection=spaces_ptr_->bisection();
      auto& tracker=bisection.tracker();
      auto& n2e=spaces_ptr_->node2elem();
      auto n_levels=tracker.current_iterate();
      Integer cont_elem=0;
      Integer next_elem=0;
      Integer el=0;
      Integer n_elements=mesh.n_elements();

      std::vector<Integer> level_n_elements(n_levels);
      std::vector<bool> elems_iterate(n_elements,false);


      if(entities_ordered_.size()<n_levels)
      {
        for(Integer i=entities_ordered_.size();i<n_levels;i++)
          entities_ordered_.push_back(SortedVec{});
      }

      for(Integer el=0;el<n_elements;el++)
      {
        level_n_elements[tracker.get_iterate(el)]++;
      }
      
      for(Integer i=0;i<n_levels;i++)
      {
        // // std::cout<<"level=="<<i<<std::endl;
        ordered_entities_aux(elems_iterate,i,level_n_elements[i]);
      }

      return entities_ordered_;
      

    }



  private:
    std::shared_ptr<FullFunctionSpace> spaces_ptr_;
    ElemDofMap elemdofmap_;
    Map map_;
    VecVecVec map_vec_;
    SortedLevelsVec entities_ordered_;
    MapBoolPtr map_bool_;
    VecVecVec levels_vec_;


  };




  template<Integer M=-1,Integer...Ms,typename FullFunctionSpace, Integer EntityDim>
  void build(Entity2Dofs<FullFunctionSpace,EntityDim>& entity2dofs)
  {

    entity2dofs.template build<false,M,Ms...>();

  }



  template<Integer M=-1,Integer...Ms,typename FullFunctionSpace, Integer EntityDim>
  void build_all_element_dofs(Entity2Dofs<FullFunctionSpace,EntityDim>& entity2dofs)
  {

    entity2dofs.template build<true,M,Ms...>();

  }




}


#endif