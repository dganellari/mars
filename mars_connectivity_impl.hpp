#ifndef MARS_CONNECTIVITY_IMPL_HPP
#define MARS_CONNECTIVITY_IMPL_HPP

#include "mars_connectivity.hpp"


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Written by Gabriele Rovi (April 2019)                                                                                            ////////
////// We define the following class:                                                                                                   ////////
////// Connectivity<Integer Dim, Integer ManifoldDim, Integer EntityDimFrom, Integer SubEntityDimFrom, Integer EntityDimTo>             ////////
////// Given a mesh, with elements whose dimension is ManifoldDim, defined in a Dim-dimensional space                                   ////////
////// We want to define the connectivity between different entities                                                                    ////////
////// Given an entity e1 of dimension EntityDimFrom, we want to find all the other entities e2 of dimension EntityDimTo                ////////
////// Such that e2 share at least SubEntityDimFrom+1 points with EntityDimFrom                                                         ////////
////// Example 1) Dim=ManifoldDim=3D, e1=triangle, SubEntityDimFrom=0, e2=triangle:                                                     ////////
//////            all the triangles e2 which share a node with the triangle e1                                                          ////////
////// Example 2) Dim=ManifoldDim=3D, e1=triangle, SubEntityDimFrom=1, e2=triangle:                                                     ////////
//////            all the triangles e2 which share an edge with the triangle e1                                                         ////////
////// Example 3) Dim=ManifoldDim=3D, e1=triangle, SubEntityDimFrom=2, e2=triangle:                                                     ////////
//////            all the triangles e2 which share a triangle with the triangle e1 -> e1=e2, known a priori                             ////////
////// Rules: 1) if SubEntityDimFrom > EntityDimFrom: the set is empty                                                                  ////////
//////        2) if SubEntityDimFrom = EntityDimFrom == EntityDimTo: the set is the entity e1 itself                                    ////////
//////        3) The class with EntityDimFrom=SubEntityDimFrom=0 and EntityDimTo=ManifoldDim is defined separately                      ////////
//////           Indeed it will store a vector of vectors, whose component (the node id) returns the neighborhood elements              ////////
//////        4) For all the other classes, given the index entity, the connectivity is computed on the fly                             ////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace mars{

//template <Integer Dim, Integer ManifoldDim, Integer EntityDimFrom, Integer SubEntityDimFrom, Integer EntityDimTo>
//class Connectivity;

template <Integer Dim, Integer ManifoldDim, Integer EntityDimFrom, Integer SubEntityDimFrom, Integer EntityDimTo>
inline std::vector<Integer> Connectivity<Dim,ManifoldDim,EntityDimFrom,SubEntityDimFrom,EntityDimTo>::compute
                                   (const Entity<Dim,ManifoldDim,EntityDimFrom> &entity_from,
                                    const Integer index_from,
                                    const Entity<Dim,ManifoldDim,EntityDimTo>   &entity_to)
{


    std::vector<Integer> fromto;
    
    static_assert(Dim>=0 , " the space dimension must be non negative ");
    static_assert(ManifoldDim>=0 , " the manifold dimension must be non negative ");
    static_assert(Dim>=ManifoldDim , " the space dimension must be greater than or equal to than the manifold dimension ");   
    static_assert(EntityDimFrom>=0 , " the entityfrom dimension must be non negative ");
    static_assert(ManifoldDim>=EntityDimFrom , " the manifold dimension must be greater than or equal to the entityfrom dimension ");
    static_assert(SubEntityDimFrom>=0 , " the subentityfrom dimension must be non negative ");
    static_assert(ManifoldDim>=SubEntityDimFrom , " the manifold dimension must be greater than or equal to the subentityfrom dimension ");   
    static_assert(EntityDimTo>=0 , " the entityto dimension must be non negative ");
    static_assert(ManifoldDim>=EntityDimTo , " the manifold dimension must be greater than or equal to the entityto dimension ");
    
    // rule 1)
    if(SubEntityDimFrom>EntityDimFrom || SubEntityDimFrom > EntityDimTo)
         return fromto;
         
    // rule 2)
    else if(SubEntityDimFrom==EntityDimFrom && EntityDimFrom==EntityDimTo)
        {
         fromto.resize(1);
         fromto[0]=index_from;
         return fromto;
         }
         
    else
       {
        Integer entity_nodes_from[EntityDimFrom+1];
        Integer entity_nodes_to[EntityDimTo+1];
        Integer common_nodes;
        std::vector<Integer> el_neighs;
        const auto& elem_index_from=entity_from.entity_2_elem(index_from)[0];
        const auto& iter_entity_from=entity_from.entity_2_elem(index_from)[1];
        const auto& nodes_from=mesh_.elem(elem_index_from).nodes;
    
        // entity_nodes_from contains the local indexes of the nodes of the entity
        Combinations<ManifoldDim + 1, EntityDimFrom+1>::generate(iter_entity_from,entity_nodes_from);

        // now entity_nodes_from contains the global indexes of the nodes of the entity
        for(Integer nn=0;nn<EntityDimFrom+1;nn++)
            entity_nodes_from[nn]=nodes_from[entity_nodes_from[nn]];
        
        // loop on all the nodes of the entity with index=index_from
        for(Integer nn=0;nn<EntityDimFrom + 1;nn++)
            {
             const Integer& nn_from=entity_nodes_from[nn];
             // loop on all the elements surrounding
             for(Integer mm=0;mm<node2elem_[nn_from].size();mm++)
                 el_neighs.push_back(node2elem_[nn_from][mm]);
            }
    
        // define all the neighborhood elements in el_neighs
        std::sort(el_neighs.begin(),el_neighs.end());
        el_neighs.erase( std::unique( el_neighs.begin(), el_neighs.end() ), el_neighs.end() );          
    
        // loop on all the surrounding elements
        for(Integer elem_iter_to=0;elem_iter_to< el_neighs.size();elem_iter_to++)
           {
            const auto &elem_index_to=el_neighs[elem_iter_to];
            const auto& nodes_to=mesh_.elem(elem_index_to).nodes;
            // loop on all the entities of dimension EntityDimTo of the element with index=elem_index_to
            for(Integer iter_entity_to=0;iter_entity_to< entity_to.entity_combinations();iter_entity_to++)
                {
                // consider the iter_entity nodes
                const auto& index_entity_to= entity_to.elem_2_entity(elem_index_to)[iter_entity_to];    
                // entity_nodes_to contains the local indexes of the nodes of the simplex   
                Combinations<ManifoldDim + 1, EntityDimTo+1>::generate(iter_entity_to,entity_nodes_to);
                // now entity_nodes_to contains the global indexes of the nodes of the entity
                for(Integer nn=0;nn<EntityDimTo+1;nn++)
                    entity_nodes_to[nn]=nodes_to[entity_nodes_to[nn]];          
                // check if the actual index_entity_to has at least (SubEntityDimFrom + 1) nodes shared with index_entity_from
                common_nodes=0;
                for(Integer nn=0;nn<EntityDimFrom + 1;nn++)
                   for(Integer mm=0;mm<EntityDimTo + 1;mm++)
                      if(entity_nodes_from[nn]==entity_nodes_to[mm]) 
                         common_nodes++;
                 
                if(common_nodes >= SubEntityDimFrom + 1)
                   fromto.push_back(index_entity_to);
                }
            }
        
        std::sort(fromto.begin(),fromto.end());
        fromto.erase( std::unique( fromto.begin(), fromto.end() ), fromto.end() );      

        return fromto;
        }
    };

template <Integer Dim, Integer ManifoldDim>
void Connectivity<Dim, ManifoldDim,0,0,ManifoldDim>::init(const Mesh<Dim,ManifoldDim> mesh)
        {
            const auto & n_nodes    = mesh.n_nodes();
            const auto & n_elements = mesh.n_elements();            
            val_.resize(n_nodes);
                // loop on all the nodes of all the active elements 
                // add to node_2_element the corresponding active element
                for(Integer i = 0; i < n_elements; ++i) 
                   {
                        if(!mesh.is_active(i)) continue;

                        const auto &e    = mesh.elem(i);
                        // this can be put outside the loop, right?
                        const auto & nn = ManifoldDim + 1;
                        for(Integer k = 0; k < nn; ++k) 
                            {
                                val_[e.nodes[k]].push_back(i);
                            }
                    }       
        };











template<Integer Dim, Integer ManifoldDim, Integer N>
struct EntitiesTupleType
{
    using rest = typename EntitiesTupleType<Dim,ManifoldDim, N-1>::type;
    using type = decltype(std::tuple_cat(std::declval<rest>() ,
                                         std::declval< std::tuple< Entity<Dim,ManifoldDim,N>>>())
                          );
};


template<Integer Dim, Integer ManifoldDim>
struct EntitiesTupleType<Dim,ManifoldDim,0>
{
 using type = typename std::tuple< Entity<Dim,ManifoldDim,0> >;
};



// specialization for N=0
template<Integer Dim, Integer ManifoldDim, Integer N>
typename std::enable_if<(N == 0u), typename EntitiesTupleType<Dim,ManifoldDim,0>::type>::type 
Entities(const Mesh<Dim, ManifoldDim> mesh, const std::vector< std::vector<Integer> >node_2_element)
{
    return std::tuple<Entity<Dim, ManifoldDim,0>>(Entity<Dim, ManifoldDim,0>(mesh,node_2_element));
}


template<Integer Dim, Integer ManifoldDim, Integer N>
typename std::enable_if<(N > 0u), typename EntitiesTupleType<Dim,ManifoldDim,N>::type>::type  
Entities(const Mesh<Dim,ManifoldDim> mesh, const std::vector< std::vector<Integer> >node_2_element )
{
    return std::tuple_cat(Entities<Dim,ManifoldDim,N-1>(mesh,node_2_element),
                          std::tuple<Entity<Dim,ManifoldDim,N>>(Entity<Dim,ManifoldDim,N>(mesh,node_2_element)) 
                             )   ;
}




/////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////        COMPUTE_LOCAL_DOFMAP         ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
// We iterate over the tuple from position 0 to the end
// If M < sizeof(Args), i.e. M < sizeof(tuple)
// Continue with M+1, otherwise end



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
template<std::size_t ManifoldDim, std::size_t M = 0 , typename ...Args>
typename std::enable_if< M==sizeof ...(Args),void >::type
    initialize_vector_entities(std::tuple< Args...> const &tuple, 
                               std::array<std::vector<bool>, ManifoldDim+1> & entity_found)
    {
     static_assert(M>=0," the tuple must have non negative length ");
     static_assert(ManifoldDim > 0," the ManifoldDim must be positive ");
        
    };
// general case (M<tuple size)
template<std::size_t ManifoldDim, std::size_t M = 0, typename...Args>
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


/////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////   INITIALIZE BOOL VECTOR ENTITIES   ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////


// base case 
template<std::size_t ManifoldDim, std::size_t M>
typename std::enable_if< M==1,void >::type
initialize_combinations_value(std::array<Integer, ManifoldDim+1> & combinations_nums)
    {

     static_assert(M>=1," the tuple must have non negative length ");
     static_assert(ManifoldDim > 0," the ManifoldDim must be positive ");
     combinations_nums[0]=Combinations<ManifoldDim+1,M>::value;
    };
// general case 
template<std::size_t ManifoldDim, std::size_t M = ManifoldDim+1>
typename std::enable_if< 1<M,void >::type
initialize_combinations_value(std::array<Integer, ManifoldDim+1> & combinations_nums)
    {
     static_assert(M>=1," the tuple must have non negative length ");
     static_assert(ManifoldDim > 0," the ManifoldDim must be positive ");
     combinations_nums[M-1]=Combinations<ManifoldDim+1,M>::value;
     // iterate to the next 
     initialize_combinations_value<ManifoldDim,M-1>(combinations_nums);       
    };



/////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////          M-th TUPLE COMPONENT       ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
// base case (M==tuple size)
template<std::size_t K,std::size_t M , typename ...Args>
typename std::enable_if< M==K,void >::type
    tuple_component(std::tuple< Args...> const &tuple) 
    {
     static_assert(M>=0," the component must be greater than zero ");
     static_assert(K<sizeof ...(Args)," the component must be inside the tuple ");
        
    };
// general case (M<tuple size)
template<std::size_t K,std::size_t M = 0, typename...Args>
typename std::enable_if< M<sizeof ...(Args),void >::type
    tuple_component(std::tuple< Args...> const &tuple) 
    {
     static_assert(M>=0," the component must be greater than zero ");
     static_assert(K<sizeof ...(Args)," the component must be inside the tuple ");
     // iterate to the next 
     tuple_component<K,M+1,Args...>(tuple);
        
    };
























using std::cout;
using std::endl;
void connectivity_example()
{




  // std::array<std::shared_ptr<BaseEntity>,3> entity_vec1;
  //   entity_vec1[0]=std::static_pointer_cast<BaseEntity> (std::make_shared< Entity<ManifoldDim,ManifoldDim,0>>(mesh,node2elem));
  //   entity_vec1[1]=std::static_pointer_cast<BaseEntity> (std::make_shared< Entity<ManifoldDim,ManifoldDim,1>>(mesh,node2elem));
  //   entity_vec1[2]=std::static_pointer_cast<BaseEntity> (std::make_shared< Entity<ManifoldDim,ManifoldDim,2>>(mesh,node2elem));
  
    constexpr Integer ManifoldDim=4;
    constexpr Integer Dim=4;
    mars::Mesh<ManifoldDim, ManifoldDim> mesh;
    read_mesh("../data/pentatope_1.MFEM", mesh);
    //read_mesh("../data/cube_6.MFEM", mesh);
    //read_mesh("../data/square_2_def.MFEM", mesh);
    
    
    auto nelements=mesh.n_elements();
//     for(auto i=0;i<nelements;i++)
//     {
//      cout<<" i="<<i<<endl;
//      auto elem0=mesh.elem(i);
//      auto vertices=elem0.nodes;
//      auto nnodes=n_nodes(elem0);
//      auto nsides=n_sides(elem0);
//      auto ndims=n_dims(elem0);
//      auto dual_graph=mesh.dual_graph();
//  
//  
//      mesh.build_dual_graph();
//      cout<<" dual_graph.size()="<<dual_graph.size()<<endl;
//      for(int ii=0;ii<dual_graph.size();++ii)
//      {
//      auto adj=dual_graph.adj(ii);
//      cout<<"dual_graph.n_adjacients="<<dual_graph.n_adjacients(ii)<<endl;
//      for(int jj=0;jj<ManifoldDim+1;++jj)
//          {
//          cout<<" dual_graph="<<adj[jj]<<endl;
//          }
//      }
//  
//  
//      for (auto row = vertices.begin(); row != vertices.end(); row++) 
//               {
//               //  for (auto col = row->begin(); col != row->end(); col++) {
//               // do stuff ...
//               auto p=(mesh.point(*row)).values;
//               for( auto point_it=p.begin(); point_it!=p.end();point_it++)
//                   cout<<" point("<<*row<<")="<<*point_it<<endl;
//                   //}
//                }
// 
//              //code here to print record wise.
// 
//      cout<<"n_nodes="<<nnodes<<endl;
//      cout<<"n_sides="<<nsides<<endl;
//      cout<<"n_dims="<<ndims<<endl;
//     }
    




    
    NodeToElem4 node2elem3(mesh);
    auto node2elem=node2elem3.val();
    





     const auto const_entities_tuple=Entities<Dim,ManifoldDim,ManifoldDim>(mesh,node2elem);
   
    const auto prova0=std::get<0>(const_entities_tuple);
    const auto edge=std::get<1>(const_entities_tuple);
    const auto triangle=std::get<2>(const_entities_tuple);
    const auto prova3=std::get<3>(const_entities_tuple);
    std::cout<<" tuple size =" << std::tuple_size<decltype(const_entities_tuple)>::value << std::endl ;
     compute_local_dofmap(const_entities_tuple); 



    std::array<std::vector<bool>, ManifoldDim+1> entity_found;
    std::array<Integer, ManifoldDim+1> combinations_nums;
    initialize_vector_entities<ManifoldDim>(const_entities_tuple,entity_found);
    initialize_combinations_value<ManifoldDim>(combinations_nums);


    for(Integer entity_dim = 0 ;entity_dim <= ManifoldDim; entity_dim++)
    {
        cout<<" entity_dim="<<entity_dim<<endl;

        for(Integer entity = 0 ;entity < combinations_nums[entity_dim]; entity++)
            {
            std::cout<<combinations_nums[entity_dim]<<" ";
            //entity_dof= qualcosa dell elemento[entity_dim]; 
            // // if the entity has not been find yet, then define dofs on this entity
            // if(alread_mapped_entity[entity_dim][entity_dof]==false)
            //     {// loop on all the finite element spaces
            //     for(Integer function_space = 0; function_space < N_function_space; function_space++)
            //         // loop on all the dimension of the field of the given function space
            //         for(Integer function_space_dim=0;function_space_dim<Dim_function_space; function_space_dim++ )
            //             // loop on the number of dofs for the given entity of the given function_space
            //             for(Integer dof_dim = 0 ;dof_dim < Dim_dofs(function_space,entity); dof_dim++)  
            //                 {

                                
            //                 }
            //     }
            // // otherwise, search for the already numbered dofs and assign them to the ones of this element
            // else
            //     {

            //     }
            }
            cout<<endl;


    }
    // for(Integer entity_iter=0; entity_iter < ManifoldDim+1 ; entity_iter++)
    // {
    //   entity_found[entity_iter].resize(std::get<entity_iter>(const_entities_tuple).entity_nums(),false);
    // }







    //const std::tuple<Entity<Dim,ManifoldDim,2>> entitytuple(Entity<Dim,ManifoldDim,2>(mesh,node2elem)); 
    
    //std::tuple<Entity<Dim,ManifoldDim,2>> entitytuple;
    //std::get<0>(entitytuple)=Entity<Dim,ManifoldDim,2>(mesh,node2elem); 
    //const auto edge=std::get<0>(entitytuple); 

    //TriangleMap4 edge(mesh,node2elem);
    //EdgeMap3 edge(mesh,node2elem);    
    //TriangleMap4 triangle(mesh,node2elem);  
    //EdgeMap3 triangle(mesh,node2elem); 

    constexpr Integer entity_from_index=7;
    constexpr Integer entitydim_from=2;
    constexpr Integer subentitydim_from=1;
    constexpr Integer entitydim_to=1;

    Integer entity_e[3];
    Integer entity_t[3];
    Connectivity<ManifoldDim,ManifoldDim,entitydim_from,subentitydim_from,entitydim_to> conn_e2t(mesh,node2elem);
        
        
        
        
        
        
    const auto &entity_2_elem_e=edge.entity_2_elem();
    const auto &elem2entity_e=edge.elem_2_entity();
        
    const auto &entity_2_elem_t=triangle.entity_2_elem();
    const auto &elem2entity_t=triangle.elem_2_entity();
    
    cout<<"ENTITY FROM #elems="<<triangle.entity_nums()<<endl;
    cout<<"ENTITY FROM dimension="<<entitydim_from<<endl;

//     for(int ii=entity_from_index;ii<entity_2_elem_t.size();ii++)
//     {
//  
//      cout<<endl;
//      for(int jj=0;jj<1;jj++)
//         cout<<entity_2_elem_t[ii][jj]<<" ";
//      cout<<"ENTITY FROM id="<<ii<<"   ";
        const auto & elemii=mesh.elem(entity_2_elem_t[entity_from_index][0]);
        Combinations<ManifoldDim + 1, triangle.num_of_points()>::generate(entity_2_elem_t[entity_from_index][1],entity_t);
        for(int jj=0;jj<triangle.num_of_points();jj++)
           cout<<elemii.nodes[entity_t[jj]]<<" ";    
        cout<<endl;
//  }  


    cout<<"SUBENTITY FROM dimension="<<subentitydim_from<<endl;
    cout<<"ENTITY TO dimension="<<entitydim_to<<endl;
    cout<<"ENTITY FROM of interest id="<<entity_from_index<<endl;
    
 

 
    cout<<"ENTITY TO #elems="<<edge.entity_nums()<<endl;
    for(int ii=0;ii<edge.entity_nums();ii++)
    {
    
        cout<<endl;
        for(int jj=0;jj<1;jj++)
           cout<<entity_2_elem_e[ii][jj]<<" ";
        cout<<"--------ENTITY TO elem id="<<ii<<"   ";
        const auto & elemii=mesh.elem(entity_2_elem_e[ii][0]);
        Combinations<ManifoldDim + 1, edge.num_of_points()>::generate(entity_2_elem_e[ii][1],entity_e);
        for(int jj=0;jj<edge.num_of_points();jj++)
           cout<<elemii.nodes[entity_e[jj]]<<" ";    
    }   

    const auto & connection=conn_e2t.compute(triangle,entity_from_index,edge);

    cout<<endl;
    for(Integer ii=0;ii<connection.size();ii++)
       cout<<" CONNECTIONS: ="<<connection[ii]<<endl;

   //int r=3;


// Integer N=3;

// for_<N>([&] (auto i) {      
//   std::get<i.value>(t); // do stuff
// });
    // std::shared_ptr<BaseEntity> ok= std::static_pointer_cast<BaseEntity> (std::make_shared< Entity<ManifoldDim,ManifoldDim,0>>(mesh,node2elem));


 

    //auto valore=entity_vec1[0]->val();

    // //for(Integer ii=0;ii<3;ii++)
    // //entity_vec1[ii]=std::static_pointer_cast<BaseEntity> (std::make_shared< Entity<ManifoldDim,ManifoldDim,ii>>(mesh,node2elem));
    // typedef std::make_shared< Entity<ManifoldDim,ManifoldDim,*it>> myType;
    //  typedef std::static_pointer_cast<BaseEntity> myStaticPointer;
    // std::array<Integer,3> myarray={0,1,2};
    //    for ( const_iterator it = myarray.begin(); it != myarray.end(); ++it )
    //    {
    //         entity_vec1[*it]= myStaticPointer(myType(mesh,node2elem));

    //    }

       // std::array<BaseEntity,2> entity_vec2;
    // entity_vec[0]=new Entity<ManifoldDim,ManifoldDim,0>(mesh,node2elem);
    // entity_vec[1]=new Entity<ManifoldDim,ManifoldDim,1>(mesh,node2elem);
    // std::array<BaseEntity,3> entity_vec3;
    // entity_vec[0]=new Entity<ManifoldDim,ManifoldDim,0>(mesh,node2elem);
    // entity_vec[1]=new Entity<ManifoldDim,ManifoldDim,1>(mesh,node2elem);
    // entity_vec[2]=new Entity<ManifoldDim,ManifoldDim,2>(mesh,node2elem);

     }
}

#endif