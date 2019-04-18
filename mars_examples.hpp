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
    read_mesh("../data/pentatope_1.MFEM", mesh);
    //read_mesh("../data/cube_6.MFEM", mesh);
    //read_mesh("../data/square_2_def.MFEM", mesh);
    
  
    NodeToElem4 node2elem3(mesh);
    auto node2elem=node2elem3.val();
    std::cout<<" Lagr11 ord1 =" <<Lagrange1_1D::entities_nums<<", "<<Lagrange1_1D::id<< std::endl ;
    std::cout<<" 11 ord1 =" <<FunctionSpace<1,1,Lagrange,1>::entities_nums<< std::endl ;
    std::cout<<" 22 ord1 =" <<FunctionSpace<2,2,Lagrange,1>::entities_nums<< std::endl ;
    std::cout<<" 33 ord1 =" <<FunctionSpace<3,3,Lagrange,1>::entities_nums<< std::endl ;
    std::cout<<" 44 ord1 =" <<FunctionSpace<4,4,Lagrange,1>::entities_nums<< std::endl ;

    std::cout<<" 22 ord2 =" <<FunctionSpace<2,2,Lagrange,2>::entities_nums<< ", "<<FunctionSpace<2,2,Lagrange,2>::entity[1]<< std::endl ;
    std::cout<<" 33 ord2 =" <<FunctionSpace<3,3,Lagrange,2>::entities_nums<< ", "<<FunctionSpace<3,3,Lagrange,2>::entity[1]<< std::endl ;
    std::cout<<" 44 ord2 =" <<FunctionSpace<4,4,Lagrange,2>::entities_nums<< ", "<<FunctionSpace<4,4,Lagrange,2>::entity[1]<< std::endl ;


    std::cout<<" 22 ord3 =" <<FunctionSpace<2,2,Lagrange,3>::entities_nums<< ", "<<FunctionSpace<2,2,Lagrange,3>::entity[2]<< std::endl ;
    std::cout<<" 33 ord3 =" <<FunctionSpace<3,3,Lagrange,3>::entities_nums<< ", "<<FunctionSpace<3,3,Lagrange,3>::entity[2]<< std::endl ;
    std::cout<<" 44 ord3 =" <<FunctionSpace<4,4,Lagrange,3>::entities_nums<< ", "<<FunctionSpace<4,4,Lagrange,3>::entity[2]<< std::endl ;



     std::cout<<" Lagrange1_3D order =" <<FunctionSpace<3,3,Lagrange,1>::order<<", "<<FunctionSpace<3,3,Lagrange,1>::space_dim<< std::endl ;
     std::cout<<" Lagrange1_2D order =" <<FunctionSpace<2,2,Lagrange,1>::order<<", "<<FunctionSpace<2,2,Lagrange,1>::space_dim<< std::endl ;

     std::cout<<" Lagrange1_1D order =" <<FunctionSpace<1,1,Lagrange,1>::order<<", "<<FunctionSpace<1,1,Lagrange,1>::space_dim<< std::endl ;
    std::cout<<" Lagrange1_4D order =" <<BaseFunctionSpace<Dim,ManifoldDim,Lagrange,Order>::order<< std::endl ;
        std::cout<<" Lagrange1_4D order =" <<FunctionSpace<Dim,ManifoldDim,Lagrange,2>::order<< std::endl ;


using T = Append<2, Entity<Dim,ManifoldDim,1>, std::tuple<int> >::type; 
// il problema e' che 
    Lagrange3_4D mmm;
    dofmap<Dim,ManifoldDim, Lagrange3_4D,Lagrange3_4D>(mesh);//,Lagrange2_4D,Lagrange3_4D);
    const auto const_entities_tuple=EntitiesOfFunctionSpace<Dim,ManifoldDim,SpaceID,Order>(mesh,node2elem);

   
    const auto node=std::get<0>(const_entities_tuple);
    const auto edge=std::get<1>(const_entities_tuple);
    const auto triangle=std::get<2>(const_entities_tuple);
    std::cout<<" tuple size =" << std::tuple_size<decltype(const_entities_tuple)>::value << std::endl ;


    constexpr const auto & entities_nums=FunctionSpace<Dim,ManifoldDim,SpaceID,Order>::entities_nums;
    std::array<std::vector<bool>, entities_nums> entity_found;
    initialize_vector_entities<Dim,ManifoldDim,SpaceID,Order>(const_entities_tuple,entity_found);
   
   for(int ii=0;ii<entity_found.size();ii++)
    {std::cout<<std::endl<<"ii=="<<ii<<std::endl;
     for(int jj=0;jj<entity_found[ii].size();jj++)
     {
        std::cout<<entity_found[ii][jj]<<" ";
     }
    }



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

        const auto & elemii=mesh.elem(entity_2_elem_t[entity_from_index][0]);
        Combinations<ManifoldDim + 1, triangle.num_of_points()>::generate(entity_2_elem_t[entity_from_index][1],entity_t);
        for(int jj=0;jj<triangle.num_of_points();jj++)
           cout<<elemii.nodes[entity_t[jj]]<<" ";    
        cout<<endl;


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