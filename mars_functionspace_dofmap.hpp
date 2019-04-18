#ifndef MARS_FUNCTIONSPACE_DOFMAP_HPP
#define MARS_FUNCTIONSPACE_DOFMAP_HPP


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Written by Gabriele Rovi (April 2019)                                                                                            ////////
////// We define the following class:                                                                                                   ////////
////// 1) EntitiesOfFunctionSpace<Integer Dim, Integer ManifoldDim, Integer SpaceID, Integer Order>:                              ////////                                        ////////
//////    given the function space SpaceID of order Order,                                                                        ////////
//////    it builds a tuple of the entities related to its Dim_dofs                                                                     ////////
////// 2) EntitiesOfFunctionSpaceType is the type of EntitiesOfFunctionSpace                                                            ////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "mars_base.hpp"
#include "mars_functionspace.hpp"
namespace mars{

template <int N, typename T, typename SeqWithArgs>      
struct Append;                                                                  

template <int N, typename T, template <typename...> class Seq, typename... Args >      
struct Append<N, T, Seq<Args...> >                                               
{                                                                               
    using type = typename Append<N-1, T, Seq<T,Args...> >::type;                         
};                                                                              

template <typename T, template<typename...> class Seq, typename... Args>              
struct Append<0, T, Seq<Args...> >                                               
{                                                                               
    using type = Seq<Args...>;                                                  
};                                                                              

// create a big tuple with N number of "int" types in it                           


template<Integer N,typename T, typename TwithArgs>
struct TupleType;

template<Integer N,typename T,template<typename...> class C, typename...Args>
struct TupleType<N, T,C<T,Args...> >
{
    using type= typename TupleType< N-1,T,C<T,Args...>>::type;
    // using rest = typename TupleType< N-1,T,C<T,Args...>>::type;
    // using ens =  C<T,Args...>;
    // using tuple_ens=std::tuple<ens>;
    // using type = decltype( std::tuple_cat( std::declval< rest >(), std::declval< tuple_ens >() ) );
};


template<typename T,template<typename...> class C, typename...Args>
struct TupleType<0,T,C<Args...> >
{
 using ens =  C<Args...>;
 using type = typename std::tuple<ens>;
};



 









template<Integer Dim, Integer ManifoldDim, Integer SpaceID, Integer Order,Integer N>
struct EntitiesOfFunctionSpaceType
{
    using rest = typename EntitiesOfFunctionSpaceType<Dim,ManifoldDim,SpaceID,Order,N-1>::type;
    using ens = Entity<Dim,ManifoldDim,FunctionSpace<Dim,ManifoldDim,SpaceID,Order>::entity[N]>;
    using tuple_ens=std::tuple<ens>;
    using type = decltype( std::tuple_cat( std::declval< rest >(), std::declval< tuple_ens >() ) );
};


template<Integer Dim, Integer ManifoldDim,Integer SpaceID, Integer Order>
struct EntitiesOfFunctionSpaceType<Dim,ManifoldDim,SpaceID,Order,0>
{
 using ens = Entity<Dim,ManifoldDim,FunctionSpace<Dim,ManifoldDim,SpaceID,Order>::entity[0]>;
 using type = typename std::tuple<ens>;
};















// specialization for N=0
template<Integer Dim, Integer ManifoldDim,Integer SpaceID, Integer Order, Integer N=(FunctionSpace<Dim,ManifoldDim,SpaceID,Order>::entities_nums-1)>
typename std::enable_if<N==0, 
                        typename EntitiesOfFunctionSpaceType<Dim,ManifoldDim,SpaceID,Order,N>::type>::type 
EntitiesOfFunctionSpace(const Mesh<Dim, ManifoldDim>&mesh, const std::vector< std::vector<Integer> > &node_2_element)
{    
    using ens=Entity<Dim,ManifoldDim,FunctionSpace<Dim,ManifoldDim,SpaceID,Order>::entity[N]>;
    return std::tuple<ens>(ens(mesh,node_2_element));
}


template<Integer Dim, Integer ManifoldDim,Integer SpaceID, Integer Order, Integer N=(FunctionSpace<Dim,ManifoldDim,SpaceID,Order>::entities_nums-1)>
typename std::enable_if< 0<N, 
                         typename EntitiesOfFunctionSpaceType<Dim,ManifoldDim,SpaceID,Order,N>::type >::type  
EntitiesOfFunctionSpace(const Mesh<Dim,ManifoldDim>&mesh, const std::vector< std::vector<Integer> >&node_2_element )
{
    using ens=Entity<Dim,ManifoldDim,FunctionSpace<Dim,ManifoldDim,SpaceID,Order>::entity[N]>;
    return std::tuple_cat(EntitiesOfFunctionSpace<Dim,ManifoldDim,SpaceID,Order,N-1>(mesh,node_2_element),
                          std::tuple<ens>(ens(mesh,node_2_element)));
}
















template<Integer Dim,Integer ManifoldDim,typename FunctionSpace,typename...FunctionSpaces>
struct EntitiesOfFunctionSpaceTupleType
{
 static constexpr Integer SpaceID=FunctionSpace::id;
 static constexpr Integer Order=FunctionSpace::order; 
 static constexpr Integer N=FunctionSpace::entities_nums-1;
 using rest = typename EntitiesOfFunctionSpaceTupleType<Dim,ManifoldDim,FunctionSpaces...>::type;
 using ens  = typename EntitiesOfFunctionSpaceType<Dim,ManifoldDim,SpaceID,Order,N>::type;
 using type = decltype( std::tuple_cat( std::declval< rest >(), 
                                        std::declval< std::tuple<ens> >() ) );
};


template<Integer Dim,Integer ManifoldDim,typename FunctionSpace>
struct EntitiesOfFunctionSpaceTupleType<Dim,ManifoldDim,FunctionSpace>
{
 static constexpr Integer SpaceID=FunctionSpace::id;
 static constexpr Integer Order=FunctionSpace::order; 
 static constexpr Integer N=FunctionSpace::entities_nums-1;
 using ens =  typename EntitiesOfFunctionSpaceType<Dim,ManifoldDim,SpaceID,Order,N>::type ;
 using type = typename std::tuple<ens>;
};




/////////////////////////////////////////////////////////////////////////////



// general case 
template<Integer Dim,Integer ManifoldDim,typename FunctionSpace,typename...FunctionSpaces>
typename std::enable_if< 0<sizeof...(FunctionSpaces), 
                         typename EntitiesOfFunctionSpaceTupleType<Dim,ManifoldDim,FunctionSpace,FunctionSpaces...>::type>::type
EntitiesOfFunctionSpaceTuple
(const Mesh<Dim, ManifoldDim>&mesh, const std::vector< std::vector<Integer> > &node_2_element)
{    
    constexpr Integer SpaceID=FunctionSpace::id;
    constexpr Integer Order=FunctionSpace::order;
    static constexpr Integer N=FunctionSpace::entities_nums-1;
    using type = typename EntitiesOfFunctionSpaceType<Dim,ManifoldDim,SpaceID,Order,N>::type;
    std::cout<<"mmm1111"<<std::endl;
    //return std::tuple<type>(EntitiesOfFunctionSpace<Dim,ManifoldDim,SpaceID,Order>(mesh,node_2_element));
    return std::tuple_cat(EntitiesOfFunctionSpaceTuple<Dim,ManifoldDim,FunctionSpaces...>(mesh,node_2_element),
                          std::tuple<type>(EntitiesOfFunctionSpace<Dim,ManifoldDim,SpaceID,Order>(mesh,node_2_element))
                         );
}

// base case
template<Integer Dim,Integer ManifoldDim,typename FunctionSpace, typename...FunctionSpaces>
 typename std::enable_if< 0==sizeof...(FunctionSpaces), 
                          typename EntitiesOfFunctionSpaceTupleType<Dim,ManifoldDim,FunctionSpace>::type>::type
//typename EntitiesOfFunctionSpaceTupleType<Dim,ManifoldDim,FunctionSpace>::type
EntitiesOfFunctionSpaceTuple
(const Mesh<Dim, ManifoldDim>&mesh, 
 const std::vector< std::vector<Integer> > &node_2_element)
{    
    constexpr Integer SpaceID=FunctionSpace::id;
    constexpr Integer Order=FunctionSpace::order;
    static constexpr Integer N=FunctionSpace::entities_nums-1;
    using type = typename EntitiesOfFunctionSpaceType<Dim,ManifoldDim,SpaceID,Order,N>::type;
    std::cout<<"mmm0000"<<std::endl;
    return std::tuple<type>(EntitiesOfFunctionSpace<Dim,ManifoldDim,SpaceID,Order>(mesh,node_2_element));
}





/////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////   INITIALIZE BOOL VECTOR ENTITIES   ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////


// base case (M==tuple size)
template<Integer Dim, Integer ManifoldDim, Integer SpaceID, Integer Order, Integer M  , typename ...Args>
typename std::enable_if< M==sizeof ...(Args),void >::type
    initialize_vector_entities(std::tuple< Args...> const &tuple, 
                               std::array<std::vector<bool>, FunctionSpace<Dim,ManifoldDim,SpaceID,Order>::entities_nums > & entity_found)
    {
     static_assert(M>=0," the tuple must have non negative length ");
     static_assert(ManifoldDim > 0," the ManifoldDim must be positive ");
        
    };
// general case (M<tuple size)
template<Integer Dim,Integer ManifoldDim, Integer SpaceID, Integer Order, Integer M = 0, typename...Args>
typename std::enable_if< M<sizeof ...(Args),void >::type
    initialize_vector_entities(std::tuple< Args...> const &tuple,
                               std::array<std::vector<bool>, FunctionSpace<Dim,ManifoldDim,SpaceID,Order>::entities_nums > & entity_found)
    {
     static_assert(M>=0," the tuple must have non negative length ");
     static_assert(ManifoldDim > 0," the ManifoldDim must be positive ");

     Integer entity_length=std::get<M>(tuple).entity_nums();
     entity_found[M].resize(entity_length,false);
     std::cout<<" (M,entity_length )=("<<M<<", "<<entity_length<<")"<<std::endl;
     // iterate to the next 
     initialize_vector_entities<Dim,ManifoldDim,SpaceID,Order,M+1,Args...>(tuple,entity_found);
        
    };








template<Integer Dim,Integer ManifoldDim,typename FunctionSpace,typename...FunctionSpaces>
struct FlagTupleType
{
 static constexpr Integer SpaceID=FunctionSpace::id;
 static constexpr Integer Order=FunctionSpace::order; 
 static constexpr Integer entities_nums=FunctionSpace::entities_nums;

 using rest = typename FlagTupleType<Dim,ManifoldDim,FunctionSpaces...>::type;
 using ens  = typename std::array<std::vector<bool>, entities_nums>;
 using type = decltype( std::tuple_cat( std::declval< rest >(), 
                                        std::declval< std::tuple<ens> >() ) );
};


template<Integer Dim,Integer ManifoldDim,typename FunctionSpace>
struct FlagTupleType<Dim,ManifoldDim,FunctionSpace>
{
 static constexpr Integer SpaceID=FunctionSpace::id;
 static constexpr Integer Order=FunctionSpace::order; 
 static constexpr Integer entities_nums=FunctionSpace::entities_nums;
 using ens =  typename std::array<std::vector<bool>, entities_nums>;
 using type = typename std::tuple<ens>;
};




template<Integer Dim,Integer ManifoldDim,typename FunctionSpace,typename...FunctionSpaces, typename...Args>
typename std::enable_if< 0<sizeof...(FunctionSpaces), 
                         typename FlagTupleType<Dim,ManifoldDim,FunctionSpace,FunctionSpaces...>::type>::type
FlagTuple
(std::tuple<Args...> tuple)
{    
    constexpr Integer SpaceID=FunctionSpace::id;
    constexpr Integer Order=FunctionSpace::order;
    static constexpr Integer M=sizeof...(FunctionSpaces);
    static constexpr Integer entities_nums=FunctionSpace::entities_nums;
    using type = typename std::array<std::vector<bool>, entities_nums>;
    type ens;

    initialize_vector_entities<Dim,ManifoldDim,SpaceID,Order>(std::get<M>(tuple),ens);
    std::cout<<"aaa1111"<<std::endl;
    return std::tuple_cat(FlagTuple<Dim,ManifoldDim,FunctionSpaces...>(tuple),
                          std::tuple<type>(ens)
                         );
}

// base case
template<Integer Dim,Integer ManifoldDim,typename FunctionSpace, typename...FunctionSpaces, typename...Args>
 typename std::enable_if< 0==sizeof...(FunctionSpaces), 
                          typename FlagTupleType<Dim,ManifoldDim,FunctionSpace>::type>::type
//typename EntitiesOfFunctionSpaceTupleType<Dim,ManifoldDim,FunctionSpace>::type
FlagTuple(std::tuple<Args...> tuple)
{    
    constexpr Integer SpaceID=FunctionSpace::id;
    constexpr Integer Order=FunctionSpace::order;
    static constexpr Integer entities_nums=FunctionSpace::entities_nums;

    using type = typename std::array<std::vector<bool>, entities_nums>;
    type ens;
    initialize_vector_entities<Dim,ManifoldDim,SpaceID,Order>(std::get<0>(tuple),ens);
    std::cout<<"aaa0000"<<std::endl;
    return std::tuple<type>(ens);
}


































template<Integer ManifoldDim, Integer N,typename FunctionSpace,typename Array,typename...Args1>
typename std::enable_if< 0<N,void>::type
Loop_Entities(std::tuple<Args1...> entitiestuple,
              Array& flagtuples)
{

Loop_Entities<ManifoldDim,N-1,FunctionSpace>(entitiestuple,flagtuples);

auto&entity=std::get<N-1>(entitiestuple);
auto&flag=flagtuples[N-1];
// I have to loop on all entities of the given type

constexpr Integer entity_dim=FunctionSpace::entity[N];
constexpr Integer entity_points=entity_dim+1;
//Combinations<ManifoldDim+1,entity_points>::value

std::cout<<"N=="<<N<<", "<<FunctionSpace::entity[N]<<std::endl;
// do stuff
 
}
template<Integer ManifoldDim,Integer N,typename FunctionSpace, typename Array,typename...Args1>
typename std::enable_if< 0==N,void>::type
Loop_Entities(std::tuple<Args1...> entitiestuple,
              Array& flagtuples)
{
std::cout<<"N=="<<N<<", "<<FunctionSpace::entity[N]<<std::endl;
// do stuff
}



template<Integer ManifoldDim,typename FunctionSpace,typename...FunctionSpaces, typename...Args1,typename...Args2>
typename std::enable_if< 0<sizeof...(FunctionSpaces),void>::type
Loop_functionspace(std::tuple<Args1...> entitiestuple,
                   std::tuple<Args2...> flagtuples)
{
 static constexpr Integer M=sizeof...(FunctionSpaces);
 const auto m1=std::get<M>(entitiestuple);
 const auto m2=std::get<M>(flagtuples);
 static constexpr Integer NN=std::tuple_size<decltype(m1)>::value-1;

 std::cout<<"Loop_functionspace B "<<std::tuple_size<decltype(m1)>::value<<", "<<m2.size()<< std::endl;
 Loop_Entities<ManifoldDim,NN,FunctionSpace>(m1,m2);
 Loop_functionspace<ManifoldDim,FunctionSpaces...>(entitiestuple,flagtuples);
};

template<Integer ManifoldDim,typename FunctionSpace,typename...FunctionSpaces, typename...Args1,typename...Args2>
typename std::enable_if< 0==sizeof...(FunctionSpaces),void>::type
Loop_functionspace(std::tuple<Args1...> entitiestuple,
                        std::tuple<Args2...> flagtuples)
{
 const auto m1=std::get<0>(entitiestuple);
 const auto m2=std::get<0>(flagtuples);
 static constexpr Integer NN=std::tuple_size<decltype(m1)>::value-1;
 std::cout<<"Loop_functionspace A "<<std::tuple_size<decltype(m1)>::value<<", "<<m2.size()<< std::endl;
 Loop_Entities<ManifoldDim,NN,FunctionSpace>(m1,m2);
};












template<Integer Dim,Integer ManifoldDim, typename...Functions>  //Integer Dim, Integer ManifoldDim,Integer SpaceID,Integer Order>//, typename...FunctionSpace>
typename std::enable_if<0==sizeof...(Functions), void>::type
dofmap(const Mesh<Dim,ManifoldDim>&mesh){};


template<Integer Dim,Integer ManifoldDim,typename FunctionSpace, typename...FunctionSpaces>
void dofmap(const Mesh<Dim,ManifoldDim>&mesh)
{
    // compute the connection node to elem (auxiliary tool, often used)
    NodeToElem<Dim,ManifoldDim> node_2_elem(mesh);
    const std::vector< std::vector<Integer> > & node2elem=node_2_elem.val();
    // compute 
    constexpr Integer SpaceID=FunctionSpace::id;
    constexpr Integer Order=FunctionSpace::order;
    //auto ehi3=EntitiesOfFunctionSpace<Dim,ManifoldDim,SpaceID,Order>(mesh,node2elem);
    constexpr Integer NN=FunctionSpace::entities_nums-1;
    // const auto entities_tuple=std::tuple <typename EntitiesOfFunctionSpaceType<Dim,ManifoldDim,SpaceID,Order,NN>::type,
    //                                       typename EntitiesOfFunctionSpaceType<Dim,ManifoldDim,SpaceID,Order,NN>::type  >
    //                                       (ehi3,ehi3);


    //std::tuple<decltype(entities_tuple)> menghia=std::tuple<decltype(entities_tuple)>(entities_tuple);
    //auto ok=std::tuple<decltype(EntitiesOfFunctionSpace<Dim,ManifoldDim,SpaceID,Order>)>(EntitiesOfFunctionSpace<Dim,ManifoldDim,SpaceID,Order>(mesh,node2elem));
    auto entitiestuple=EntitiesOfFunctionSpaceTuple<Dim,ManifoldDim,FunctionSpace,FunctionSpaces...>(mesh,node2elem);
    using eccola=typename FlagTupleType<Dim,ManifoldDim,FunctionSpace,FunctionSpaces...>::type;
//    auto flagtuples= FlagTuple<Dim,ManifoldDim,FunctionSpace>(entitiestuple);
    auto flagtuples= FlagTuple<Dim,ManifoldDim,FunctionSpace,FunctionSpaces...>(entitiestuple);

    auto first_flag=std::get<0>(flagtuples);
    auto second_flag=std::get<1>(flagtuples);

    std::cout<<"first_flag="<<first_flag.size()<<std::endl;
    std::cout<<first_flag[0].size()<<", "<<first_flag[1].size()<<", "<<first_flag[2].size()<<std::endl;
    std::cout<<"second_flag="<<second_flag.size()<<std::endl;
    std::cout<<second_flag[0].size()<<", "<<second_flag[1].size()<<", "<<second_flag[2].size()<<std::endl;

    auto first_space=std::get<0>(entitiestuple);
    auto first_space1=std::get<0>(first_space);
    auto first_space2=std::get<1>(first_space);
    auto first_space3=std::get<2>(first_space);

    auto second_entity=std::get<1>(entitiestuple);
     auto second_entity1=std::get<0>(second_entity);
    auto second_entity2=std::get<1>(second_entity);
    auto second_entity3=std::get<2>(second_entity);
    

    // auto third_entity=std::get<2>(entitiestuple);
    //   auto third_entity1=std::get<0>(third_entity);
    // auto third_entity2=std::get<1>(third_entity);
    // auto third_entity3=std::get<2>(third_entity);
    

    std::cout<<" tuple size "<<std::tuple_size<decltype(entitiestuple)>::value<<std::endl;
    //auto second_entity=std::get<1>(entitiestuple);
    // loop on all the elements
    for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
    {
     // change it for smarter algorithms   
     auto &elem_id=elem_iter;
     // loop on all the function spaces
     Loop_functionspace<ManifoldDim,FunctionSpace,FunctionSpaces...>(entitiestuple,flagtuples);
    }

};





























 
// base case (M==tuple size)
template<Integer ManifoldDim,Integer FunctionSpaceDim,Integer EntityDim, typename FunctionSpace,typename ...FunctionSpaces,typename ...TupleComponentTypes>
typename std::enable_if< 0==sizeof...(FunctionSpaces),void >::type
    Loop_functionspace_iter(Integer elem_id, 
                            std::tuple< TupleComponentTypes...> const &entity_tuple,
                            std::array<std::vector<bool>, ManifoldDim+1> const&entity_found,
                            FunctionSpace functionspace,
                            FunctionSpaces... otherfunctionspaces) 
    {
     static_assert(sizeof...(FunctionSpaces)>=0," the component must be greater than zero ");
    };


template<Integer ManifoldDim,Integer FunctionSpaceDim,Integer EntityDim, typename FunctionSpace, typename ...FunctionSpaces,typename ...TupleComponentTypes>
typename std::enable_if< 0<sizeof...(FunctionSpaces),void >::type
    Loop_functionspace_iter(Integer elem_id, 
                     std::tuple< TupleComponentTypes...> const &entity_tuple,
                     std::array<std::vector<bool>, ManifoldDim+1> const&entity_found,
                     FunctionSpace functionspace,
                     FunctionSpaces... otherfunctionspaces) 
     {  // if the current space has dofs on this entity
        if(functionspace.value[EntityDim]>0)
        {
        std::cout<<"Loop_functionspace_iter=="<<functionspace.space_dim<<std::endl;
        for(Integer fs_dim=0;fs_dim<functionspace.functionspace_dim;fs_dim++)
           {
            std::cout<<" CAZZO==="<<functionspace.value[EntityDim]<<std::endl;
           }
        }
        Loop_functionspace_iter<ManifoldDim,FunctionSpaceDim,EntityDim>(elem_id,entity_tuple,entity_found,otherfunctionspaces...);
     };






// base case (M==tuple size)
template<Integer ManifoldDim,Integer FunctionSpaceDim,Integer EntityDim, Integer M,  
         typename FunctionSpace,typename ...FunctionSpaces,typename ...TupleComponentTypes>
typename std::enable_if< M==Combinations<ManifoldDim+1,EntityDim+1>::value,void >::type
    Loop_entity_iter(Integer elem_id, 
                     std::tuple< TupleComponentTypes...> const &entity_tuple,
                     std::array<std::vector<bool>, ManifoldDim+1> const&entity_found,
                     FunctionSpace functionspace,
                     FunctionSpaces... otherfunctionspaces) 
    {
     static_assert(M>=0," the component must be greater than zero ");
    };



// general case (M<tuple size)
template<Integer ManifoldDim,Integer FunctionSpaceDim,Integer EntityDim, Integer  M = 0,  
         typename FunctionSpace,typename ...FunctionSpaces,typename ...TupleComponentTypes>
typename std::enable_if< M<Combinations<ManifoldDim+1,EntityDim+1>::value,void >::type
    Loop_entity_iter(Integer elem_id, 
                     std::tuple< TupleComponentTypes...> const &entity_tuple,
                     std::array<std::vector<bool>, ManifoldDim+1> const&entity_found,
                     FunctionSpace functionspace,
                     FunctionSpaces... otherfunctionspaces
                    ) 
    {
     static_assert(M>=0," the component must be greater than zero ");
     const auto& entity_id=std::get<EntityDim>(entity_tuple).elem_2_entity(elem_id)[M];

     // the entity has not been found yet
     if(entity_found[EntityDim][entity_id]==false)
     {
        std::cout<<"2 sizeof...(FunctionSpaces)="<<sizeof...(FunctionSpaces)<<std::endl;
        Loop_functionspace_iter<ManifoldDim,FunctionSpaceDim,EntityDim>(elem_id,entity_tuple,entity_found,functionspace,otherfunctionspaces...);



                  //     {// loop on all the finite element spaces
            //     for(Integer function_space = 0; function_space < N_function_space; function_space++)
                    // loop on all the dimension of the field of the given function space
                    for(Integer function_space_dim=0;function_space_dim<FunctionSpaceDim; function_space_dim++ )
                    {
                     std::cout<<"function_space_dim=="<<function_space_dim<<std::endl;   
                    }
                    //     // loop on the number of dofs for the given entity of the given function_space
                    //     for(Integer dof_dim = 0 ;dof_dim < Dim_dofs(function_space,entity); dof_dim++)  
                    //         {

                                
            //                 }
            //     }

     }
     else
     {}
     std::cout<<"entity_id======="<<entity_id<<std::endl;
     // iterate to the next 
     Loop_entity_iter<ManifoldDim,FunctionSpaceDim,EntityDim,M+1>(elem_id,entity_tuple,entity_found,functionspace,otherfunctionspaces...);

    };



}


#endif