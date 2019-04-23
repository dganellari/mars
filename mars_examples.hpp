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
    constexpr Integer Dim=4;
    constexpr Integer SpaceID=Lagrange;
    constexpr Integer Order=3;
    mars::Mesh<ManifoldDim, ManifoldDim> mesh;
    read_mesh("../data/pentatope_2.MFEM", mesh);
    //read_mesh("../data/cube_6.MFEM", mesh);
    //read_mesh("../data/square_2_def.MFEM", mesh);
    
  
    NodeToElem4 node2elem3(mesh);
    auto node2elem=node2elem3.val();
    // std::cout<<" Lagr11 ord1 =" <<Lagrange1_1D::entities_nums<<", "<<Lagrange1_1D::id<< std::endl ;
    // std::cout<<" 11 ord1 =" <<FunctionSpace<1,1,Lagrange,1>::entities_nums<< std::endl ;
    // std::cout<<" 22 ord1 =" <<FunctionSpace<2,2,Lagrange,1>::entities_nums<< std::endl ;
    // std::cout<<" 33 ord1 =" <<FunctionSpace<3,3,Lagrange,1>::entities_nums<< std::endl ;
    // std::cout<<" 44 ord1 =" <<FunctionSpace<4,4,Lagrange,1>::entities_nums<< std::endl ;

    // std::cout<<" 22 ord2 =" <<FunctionSpace<2,2,Lagrange,2>::entities_nums<< ", "<<FunctionSpace<2,2,Lagrange,2>::entity[1]<< std::endl ;
    // std::cout<<" 33 ord2 =" <<FunctionSpace<3,3,Lagrange,2>::entities_nums<< ", "<<FunctionSpace<3,3,Lagrange,2>::entity[1]<< std::endl ;
    // std::cout<<" 44 ord2 =" <<FunctionSpace<4,4,Lagrange,2>::entities_nums<< ", "<<FunctionSpace<4,4,Lagrange,2>::entity[1]<< std::endl ;


    // std::cout<<" 22 ord3 =" <<FunctionSpace<2,2,Lagrange,3>::entities_nums<< ", "<<FunctionSpace<2,2,Lagrange,3>::entity[2]<< std::endl ;
    // std::cout<<" 33 ord3 =" <<FunctionSpace<3,3,Lagrange,3>::entities_nums<< ", "<<FunctionSpace<3,3,Lagrange,3>::entity[2]<< std::endl ;
    // std::cout<<" 44 ord3 =" <<FunctionSpace<4,4,Lagrange,3>::entities_nums<< ", "<<FunctionSpace<4,4,Lagrange,3>::entity[2]<< std::endl ;



    //  std::cout<<" Lagrange1_3D order =" <<FunctionSpace<3,3,Lagrange,1>::order<<", "<<FunctionSpace<3,3,Lagrange,1>::space_dim<< std::endl ;
    //  std::cout<<" Lagrange1_2D order =" <<FunctionSpace<2,2,Lagrange,1>::order<<", "<<FunctionSpace<2,2,Lagrange,1>::space_dim<< std::endl ;

    //  std::cout<<" Lagrange1_1D order =" <<FunctionSpace<1,1,Lagrange,1>::order<<", "<<FunctionSpace<1,1,Lagrange,1>::space_dim<< std::endl ;
    // std::cout<<" Lagrange1_4D order =" <<BaseFunctionSpace<Dim,ManifoldDim,Lagrange,Order>::order<< std::endl ;
    //     std::cout<<" Lagrange1_4D order =" <<FunctionSpace<Dim,ManifoldDim,Lagrange,2>::order<< std::endl ;

   MixedFunctionSpace<Lagrange1_4D,Lagrange2_4D,Lagrange3_4D> mixed_space;
    
   // constexpr std::array<std::array<Integer,2>,2> iu=ArrayHolder<Lagrange2_4D::dofs_per_entity, Lagrange2_4D::dofs_per_entity>::data;
     std::array<std::array<Integer, 5> , 2 > tmp2;
    ArrOfArr<0,2,5,Lagrange1_4D,Lagrange2_4D>(tmp2);
    std::cout<<tmp2[0][0]<<" "<<tmp2[0][1]<<" "<<tmp2[0][2]<<" "<<tmp2[0][3]<<" "<<" "<<tmp2[0][4]<<std::endl;
    std::cout<<tmp2[1][0]<<" "<<tmp2[1][1]<<" "<<tmp2[1][2]<<" "<<tmp2[1][3]<<" "<<" "<<tmp2[1][4]<<std::endl;
     //std::cout<<"xxx"<<Lagrange2_4D::dofs_per_entity2[1]<<std::endl;

    //constexpr std::array<Integer,2> eic=Lagrange2_4D::dofs_per_entity2;
    constexpr std::array<Integer, 4> tmp={0,4,23,2};
   std::cout << mixed_space.entities_nums(0) << std::endl;
   std::cout << mixed_space.entities_nums(1) << std::endl;
   std::cout << mixed_space.entities_nums(2) << std::endl;

   std::cout << mixed_space.dofs_per_entity<0>()[0] << std::endl;
   std::cout << mixed_space.dofs_per_entity<1>()[1] << std::endl;
   std::cout << mixed_space.dofs_per_entity<2>()[1] << std::endl;


//    MixedFunctionSpace<Lagrange1_4D, MixedFunctionSpace<Lagrange1_4D> > mixed_space2;
   //  std::cout << std::get<0>(kkk)[0] << std::endl;
   // std::cout << std::get<1>(kkk)[1] << std::endl;
   // std::cout << std::get<2>(kkk)[1] << std::endl;
   // // std::cout << mgs.kkk<1>() << std::endl;
   // std::cout << kkk<2>() << std::endl;
    // auto ooooo=TupleCreate<Lagrange3_4D,Lagrange3_4D>();//mesh);
    // auto ooooo1=std::get<1>(cane0);
     // std::cout<<" --------11---"<< ooooo1::entities_nums<<std::endl;
   std::array<Integer,Lagrange3_4D::entities_nums> arr3;
   std::array<Integer,Lagrange2_4D::entities_nums> arr2;
   FunctionSpaceOffSetDofs<Lagrange3_4D::entities_nums,Lagrange3_4D>(arr3);
   FunctionSpaceOffSetDofs<Lagrange2_4D::entities_nums,Lagrange2_4D>(arr2);

std::cout<<"Lagrange2_4D"<<std::endl;
std::cout<<FunctionSpaceDofsPerElem<Lagrange2_4D,0>::value<<std::endl;
std::cout<<FunctionSpaceDofsPerElem<Lagrange2_4D,1>::value<<std::endl;
std::cout<<"Lagrange3_4D"<<std::endl;
std::cout<<FunctionSpaceDofsPerElem<Lagrange3_4D,0>::value<<std::endl;
std::cout<<FunctionSpaceDofsPerElem<Lagrange3_4D,1>::value<<std::endl;
std::cout<<FunctionSpaceDofsPerElem<Lagrange3_4D,2>::value<<std::endl;

   // const auto& eee=TupleTemplate<0,Lagrange1_4D,Lagrange2_4D>();


   //  const auto eee0=std::get<0>(eee);
   //  const auto eee1=std::get<1>(eee);
   // // const auto eee2=std::get<2>(eee);
   // for(Integer nn=0;nn<eee0.size();nn++)
   //     std::cout<<"0__________=="<<eee0[nn]  <<std::endl;
   // for(Integer nn=0;nn<eee1.size();nn++)
   //      std::cout<<"1__________=="<<eee1[nn]  <<std::endl;
   // for(Integer nn=0;nn<eee0.size();nn++)
   //     std::cout<<"2__________=="<<eee2[nn]  <<std::endl;

   
   for(Integer nn=0;nn<arr2.size();nn++)
      std::cout<<"hei2=="<<arr2[nn]  <<std::endl;
   for(Integer nn=0;nn<arr3.size();nn++)
      std::cout<<"hei3=="<<arr3[nn]  <<std::endl;
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
dofmap<Dim,ManifoldDim, Lagrange1_4D,Lagrange3_4D>(mesh);//,Lagrange2_4D,Lagrange3_4D);
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


    const auto const_entities_tuple=EntitiesOfFunctionSpace<Dim,ManifoldDim,SpaceID,Order>(mesh,node2elem);

   
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
    cout<<"ENTITY FROM of interest id="<<entity_from_index<<endl;
    
 

 
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