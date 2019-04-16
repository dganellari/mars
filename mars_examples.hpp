#ifndef MARS_EXAMPLES_HPP
#define MARS_EXAMPLES_HPP

#include "mars_connectivity.hpp"
#include "mars_functionspace_dofmap.hpp"



namespace mars{


using std::cout;
using std::endl;
void connectivity_example()
{

    constexpr Integer ManifoldDim=4;
    constexpr Integer FunctionSpaceDim=2;
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
    




    const Integer vec[4]={0,2,3,4};

    NodeToElem4 node2elem3(mesh);
    auto node2elem=node2elem3.val();
 





     static constexpr std::array<Integer,4> entity_used={0,2,3,4};


    //try_template<4,decltype(RT1), RT1>();//mesh,node2elem);
    // LagrangianSpace<Dim,ManifoldDim,1> P111;
    // static constexpr RaviartSpace<Dim,ManifoldDim,0> RT1;
      const auto const_entities_tuple=EntitiesOfFunctionSpace<Dim,ManifoldDim,GeneralSpace,0>(mesh,node2elem);

   
    const auto prova0=std::get<0>(const_entities_tuple);
    const auto edge=std::get<1>(const_entities_tuple);
    const auto triangle=std::get<2>(const_entities_tuple);
    const auto prova3=std::get<3>(const_entities_tuple);
    std::cout<<" tuple size =" << std::tuple_size<decltype(const_entities_tuple)>::value << std::endl ;
     compute_local_dofmap(const_entities_tuple); 


    std::array<std::vector<bool>, ManifoldDim+1> entity_found;
    initialize_vector_entities<ManifoldDim>(const_entities_tuple,entity_found);








    for(Integer entity_dim = 0 ;entity_dim <= ManifoldDim; entity_dim++)
    {
        cout<<" entity_dim="<<entity_dim<<endl;

        //for(Integer entity = 0 ;entity < combinations_nums[entity_dim]; entity++)
            // {
            // std::cout<<combinations_nums[entity_dim]<<" ";
            // //const auto& const_entities_tuple.elem_2_entity[elem_iter];
            // //entity_dof= qualcosa dell elemento[entity_dim]; 
            // // // if the entity has not been find yet, then define dofs on this entity
            // // if(alread_mapped_entity[entity_dim][entity_dof]==false)
            // //     {// loop on all the finite element spaces
            // //     for(Integer function_space = 0; function_space < N_function_space; function_space++)
            // //         // loop on all the dimension of the field of the given function space
            // //         for(Integer function_space_dim=0;function_space_dim<Dim_function_space; function_space_dim++ )
            // //             // loop on the number of dofs for the given entity of the given function_space
            // //             for(Integer dof_dim = 0 ;dof_dim < Dim_dofs(function_space,entity); dof_dim++)  
            // //                 {

                                
            // //                 }
            // //     }
            // // // otherwise, search for the already numbered dofs and assign them to the ones of this element
            // // else
            // //     {

            // //     }
            // }
            // cout<<endl;


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

 

     }


}



#endif