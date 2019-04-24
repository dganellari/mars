#ifndef MARS_EXAMPLES_HPP
#define MARS_EXAMPLES_HPP

#include "mars_connectivity.hpp"
#include "mars_functionspace_dofmap.hpp"
#include "mars_simplex.hpp"



namespace mars{


// template<Integer N>
// vofamily dimmi()
// {
//     std::cout<<"dimmi_"<<N<<std::endl;
// } 

// template<typename Mesh>
// vofamily qualcosa()
// {

//     constexpr Integer dim=Mesh::Dim;
//     std::cout<<"qualcosa_"<<dim<<std::endl;
//     dimmi<dim>();
// } 
//        qualcosa<mars::Mesh<ManifoldDim, ManifoldDim>>();

using std::cout;
using std::endl;
void connectivity_example()
{

    constexpr Integer ManifoldDim=4;
    constexpr Integer Dim=4;
    constexpr Integer FEFamily=LagrangeFE;
    constexpr Integer Order=3;
    mars::Mesh<ManifoldDim, ManifoldDim> mesh;
    read_mesh("../data/pentatope_2.MFEM", mesh);
    //read_mesh("../data/cube_6.MFEM", mesh);
    //read_mesh("../data/square_2_def.MFEM", mesh);

     std::cout<<"LagrangeFE="<<LagrangeFE<<std::endl;
     std::cout<<ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,1,1,1>::space_dim<<" "<<std::endl;
    std::cout<<ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,1,0,1>::manifold_dim<<" "<<std::endl;
     std::cout<<ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,1,0,1>::family<<" "<<std::endl;
    std::cout<<ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,1,0,1>::order<<" "<<std::endl;
    std::cout<<ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,1,0,1>::n_components<<" "<<std::endl;
    
    
    NodeToElem4 node2elem3(mesh);
    auto node2elem=node2elem3.val();
    dofmap< Dim,ManifoldDim, Lagrange1< Simplex<4,4>>   , Lagrange3< Simplex<4,4> ,1,1> >(mesh);//,Lagrange2_4D,Lagrange3_4D);


    
    const auto const_entities_tuple=EntitiesOfFunctionSpace<Dim,ManifoldDim,GeneralSpace,0>(mesh,node2elem);
    const auto node=std::get<0>(const_entities_tuple);
    const auto edge=std::get<1>(const_entities_tuple);
    const auto triangle=std::get<2>(const_entities_tuple);
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
    
    cout<<"ENTITY FROM #elems="<<triangle.size()<<endl;
    cout<<"ENTITY FROM dimension="<<entitydim_from<<endl;

        const auto & elemii=mesh.elem(entity_2_elem_t[entity_from_index][0]);
        Combinations<ManifoldDim + 1, triangle.num_of_points()>::generate(entity_2_elem_t[entity_from_index][1],entity_t);
        for(int jj=0;jj<triangle.num_of_points();jj++)
           cout<<elemii.nodes[entity_t[jj]]<<" ";    
        cout<<endl;


    cout<<"SUBENTITY FROM dimension="<<subentitydim_from<<endl;
    cout<<"ENTITY TO dimension="<<entitydim_to<<endl;
    cout<<"ENTITY FROM of interest family="<<entity_from_index<<endl;
    
 

 
    cout<<"ENTITY TO #elems="<<edge.size()<<endl;
    for(int ii=0;ii<edge.size();ii++)
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