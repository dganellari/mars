#ifndef MARS_EXAMPLES_HPP
#define MARS_EXAMPLES_HPP

#include "mars_simplex.hpp"
#include "mars_connectivity.hpp"
#include "mars_functionspace_dofmap.hpp"
#include "mars_functionspace.hpp"




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
    using MeshT=mars::Mesh<Dim, ManifoldDim>;
    MeshT mesh;
    using Elem = typename MeshT::Elem; 
    read_mesh("../data/pentatope_2.MFEM", mesh);
    //read_mesh("../data/cube_6.MFEM", mesh);
    //read_mesh("../data/square_2_def.MFEM", mesh);

     std::cout<<"LagrangeFE="<<LagrangeFE<<std::endl;
     std::cout<<"ElementFunctionSpace<Dim,ManifoldDim,LagrangeFE,1>::space_dim="<<ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,1>::space_dim<<std::endl;
     std::cout<<"ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,1>::continuity="<<ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,1>::continuity<<std::endl;
     std::cout<<"ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,1>::n_components="<<ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,1>::n_components<<std::endl;

     std::cout<<ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,1,1,1>::space_dim<<" "<<std::endl;
    std::cout<<ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,1,0,1>::manifold_dim<<" "<<std::endl;
     std::cout<<ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,1,0,1>::family<<" "<<std::endl;
    std::cout<<ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,1,0,1>::order<<" "<<std::endl;
    std::cout<<ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,1,0,1>::n_components<<" "<<std::endl;
    
    

    ElemNodeToElem<Simplex<4,4>> node2elem3(mesh);
    auto node2elem=node2elem3.val();
    //auto vec=dofmap<ElemLagrange1<Elem>   , ElemLagrange3<Elem,1,1> >(mesh);





    // constexpr auto n_spaces=2;
    // constexpr auto dofs_per_elem=DofsPerElemNums1<Elem,RT0<1>,Lagrange3<1>>::value;
    // std::vector<std::array<Integer,dofs_per_elem>> dofmap_vec;
    // std::array<std::vector<Integer>,n_spaces> offset;
    // dofmap1<RT0<1>,Lagrange3<1>>(mesh,dofmap_vec,offset);
    

    FunctionSpace< MeshT, Lagrange1<1>, Lagrange3<1> > FEspace(mesh);

   std::cout<<"n_dofs="<<FEspace.n_dofs()<< std::endl;
    std::cout<<std::endl;
    for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
    {
     auto &elem_id=elem_iter;
     std::cout<<"elem_id="<<elem_id<<", "<<FEspace.dofmap(elem_id).size()<<std::endl;
     for(Integer nn=0;nn<FEspace.dofmap(elem_id).size();nn++)
     {
        std::cout<<FEspace.dofmap(elem_id)[nn]<<" ";
     }
     std::cout<<std::endl;
    } 

for(Integer ss=0;ss<FEspace.n_subspaces();ss++)
{
   std::cout<<" dofs of space ="<<ss<<std::endl;
    auto& space0=FEspace.space_dofs(ss);
    for(Integer mm=0;mm<FEspace.space_dofs(ss).size();mm++)
       std::cout<<FEspace.space_dofs(ss)[mm]<<" ";
   std::cout<<std::endl;

}
 

    for(Integer mm=0;mm<FEspace.offset().size();mm++)
     {
      std::cout<<" offset space="<<mm<<std::endl;
      for(Integer nn=0;nn<FEspace.offset()[mm].size();nn++)
         {
            std::cout<< FEspace.offset()[mm][nn]<<" ";
         }

     }
     std::cout<<std::endl;


   std::cout<<"--size="<<FEspace.dofmap(1,0).size()<<std::endl;
   for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
    {std::cout<<std::endl;
      const auto size=FEspace.dofmap(1,elem_iter).size();
      std::cout<<"size="<<size<<std::endl;
      for(Integer nn=0;nn<size;nn++)
          std::cout<<FEspace.dofmap(1,elem_iter)[nn]<<" ";
      std::cout<<std::endl;
     }


    const auto const_entities_tuple=EntitiesOfFunctionSpace<Elem,GeneralSpace,0>(mesh,node2elem);
    const auto node=std::get<0>(const_entities_tuple);
    const auto edge=std::get<1>(const_entities_tuple);
    const auto triangle=std::get<2>(const_entities_tuple);
    constexpr Integer entity_from_index=7;
    constexpr Integer entitydim_from=2;
    constexpr Integer subentitydim_from=1;
    constexpr Integer entitydim_to=1;

    Integer entity_e[3];
    Integer entity_t[3];
     ElemConnectivity<Elem,entitydim_from,subentitydim_from,entitydim_to> conn_e2t(mesh,node2elem);
     
        
        
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