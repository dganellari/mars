#ifndef MARS_FUNCTIONSPACE_DOFMAP_HPP
#define MARS_FUNCTIONSPACE_DOFMAP_HPP


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Written by Gabriele Rovi (April 2019)                                                                                            ////////
////// We define the following class:                                                                                                   ////////
////// 1) EntitiesOfFunctionSpace<Integer Dim, Integer ManifoldDim, Integer SpecificSpace, Integer Order>:                              ////////                                        ////////
//////    given the function space SpecificSpace of order Order,                                                                        ////////
//////    it builds a tuple of the entities related to its Dim_dofs                                                                     ////////
////// 2) EntitiesOfFunctionSpaceType is the type of EntitiesOfFunctionSpace                                                            ////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "mars_base.hpp"
#include "mars_functionspace.hpp"
namespace mars{


template<Integer Dim, Integer ManifoldDim, Integer SpecificSpace, Integer Order,Integer N>
struct EntitiesOfFunctionSpaceType
{
    using rest = typename EntitiesOfFunctionSpaceType<Dim,ManifoldDim,SpecificSpace,Order,N-1>::type;
    using ens = Entity<Dim,ManifoldDim,FunctionSpace<Dim,ManifoldDim,SpecificSpace,Order>::entity[N]>;
    using tuple_ens=std::tuple<ens>;
    using type = decltype( std::tuple_cat( std::declval< rest >(), std::declval< tuple_ens >() ) );
};


template<Integer Dim, Integer ManifoldDim,Integer SpecificSpace, Integer Order>
struct EntitiesOfFunctionSpaceType<Dim,ManifoldDim,SpecificSpace,Order,0>
{
 using ens = Entity<Dim,ManifoldDim,FunctionSpace<Dim,ManifoldDim,SpecificSpace,Order>::entity[0]>;
 using type = typename std::tuple<ens>;
};





// specialization for N=0
template<Integer Dim, Integer ManifoldDim,Integer SpecificSpace, Integer Order, Integer N>
typename std::enable_if<N==0, 
                        typename EntitiesOfFunctionSpaceType<Dim,ManifoldDim,SpecificSpace,Order,N> ::type>::type 
EntitiesOfFunctionSpace(const Mesh<Dim, ManifoldDim> mesh, const std::vector< std::vector<Integer> >node_2_element)
{    
    using ens=Entity<Dim,ManifoldDim,FunctionSpace<Dim,ManifoldDim,SpecificSpace,Order>::entity[N]>;
    return std::tuple<ens>(ens(mesh,node_2_element));
}


template<Integer Dim, Integer ManifoldDim,Integer SpecificSpace, Integer Order, Integer N=(FunctionSpace<Dim,ManifoldDim,SpecificSpace,Order>::entities_nums-1)>
typename std::enable_if< 0<N, 
                         typename EntitiesOfFunctionSpaceType<Dim,ManifoldDim,SpecificSpace,Order,N>::type >::type  
EntitiesOfFunctionSpace(const Mesh<Dim,ManifoldDim> mesh, const std::vector< std::vector<Integer> >node_2_element )
{
    using ens=Entity<Dim,ManifoldDim,FunctionSpace<Dim,ManifoldDim,SpecificSpace,Order>::entity[N]>;
    return std::tuple_cat(EntitiesOfFunctionSpace<Dim,ManifoldDim,SpecificSpace,Order,N-1>(mesh,node_2_element),
                          std::tuple<ens>(ens(mesh,node_2_element)));
}















/////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// COMPUTE_LOCAL_DOFMAP NON-CONST CASE ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////


// base case (M==tuple size)
template<std::size_t M = 0 , typename ...Args>
typename std::enable_if< M==sizeof ...(Args),void >::type
    compute_local_dofmap(std::tuple< Args...> const &tuple )
    {

        
    };
// general case (M<tuple size)
template<std::size_t M = 0, typename...Args>
typename std::enable_if< M<sizeof ...(Args),void >::type
    compute_local_dofmap(std::tuple< Args...> const &tuple )
    {
     // iterate to the next 
     compute_local_dofmap<M+1, Args...>(tuple);
        
    };



/////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////   INITIALIZE BOOL VECTOR ENTITIES   ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////


// base case (M==tuple size)
template<Integer ManifoldDim, Integer M  , typename ...Args>
typename std::enable_if< M==sizeof ...(Args),void >::type
    initialize_vector_entities(std::tuple< Args...> const &tuple, 
                               std::array<std::vector<bool>, ManifoldDim+1> & entity_found)
    {
     static_assert(M>=0," the tuple must have non negative length ");
     static_assert(ManifoldDim > 0," the ManifoldDim must be positive ");
        
    };
// general case (M<tuple size)
template<Integer ManifoldDim, Integer M = 0, typename...Args>
typename std::enable_if< M<sizeof ...(Args),void >::type
    initialize_vector_entities(std::tuple< Args...> const &tuple,
                               std::array<std::vector<bool>, ManifoldDim+1> & entity_found
                                )
    {
     static_assert(M>=0," the tuple must have non negative length ");
     static_assert(ManifoldDim > 0," the ManifoldDim must be positive ");

     std::size_t entity_length=std::get<M>(tuple).entity_nums();
     entity_found[M].resize(entity_length,false);
     std::cout<<" (M,entity_length )=("<<M<<", "<<entity_length<<")"<<std::endl;
     // iterate to the next 
     initialize_vector_entities<ManifoldDim,M+1,Args...>(tuple,entity_found);
        
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


//template<Integer Dim, Integer ManifoldDim, Integer N>
// struct AllEntitiesType
// {
//     using rest = typename AllEntitiesType<Dim,ManifoldDim, N-1>::type;
//     using type = decltype(std::tuple_cat(std::declval<rest>() ,
//                                          std::declval< std::tuple< Entity<Dim,ManifoldDim,N>>>())
//                           );
// };


// template<Integer Dim, Integer ManifoldDim>
// struct AllEntitiesType<Dim,ManifoldDim,0>
// {
//  using type = typename std::tuple< Entity<Dim,ManifoldDim,0> >;
// };



// // specialization for N=0
// template<Integer Dim, Integer ManifoldDim, Integer N>
// typename std::enable_if<(N == 0u), typename AllEntitiesType<Dim,ManifoldDim,0>::type>::type 
// AllEntities(const Mesh<Dim, ManifoldDim> mesh, const std::vector< std::vector<Integer> >node_2_element)
// {
//     return std::tuple<Entity<Dim, ManifoldDim,0>>(Entity<Dim, ManifoldDim,0>(mesh,node_2_element));
// }


// template<Integer Dim, Integer ManifoldDim, Integer N>
// typename std::enable_if<(N > 0u), typename AllEntitiesType<Dim,ManifoldDim,N>::type>::type  
// AllEntities(const Mesh<Dim,ManifoldDim> mesh, const std::vector< std::vector<Integer> >node_2_element )
// {
//     return std::tuple_cat(AllEntities<Dim,ManifoldDim,N-1>(mesh,node_2_element),
//                           std::tuple<Entity<Dim,ManifoldDim,N>>(Entity<Dim,ManifoldDim,N>(mesh,node_2_element)) 
//                              )   ;
// }
#endif