#ifndef MARS_EXAMPLES_HPP
#define MARS_EXAMPLES_HPP

#include "mars_simplex.hpp"
#include "mars_connectivity.hpp"
#include "mars_functionspace_dofmap.hpp"
#include "mars_functionspace.hpp"
#include "mars_shape_function.hpp"
#include "mars_matrix.hpp"



namespace mars{


using std::cout;
using std::endl;














































void cents_example()
{
    using Point=Vector<Real,3>;
    Point point{0,1,2};

    constexpr Integer ManifoldDim=2;
    constexpr Integer Dim=2;
    using MeshT=mars::Mesh<Dim, ManifoldDim>;
    MeshT mesh;
    using Elem = typename MeshT::Elem; 
    read_mesh("../data/beam-tri.MFEM", mesh);
    constexpr Integer QPOrder=4;
    constexpr Integer NQPoints=GaussPoints<Elem,QPOrder>::NQPoints; 
    using QP=typename GaussPoints<Elem,QPOrder>::qp_points_type;
    GaussPoints<Elem,QPOrder> gauss;
    const auto& qp_points=gauss.qp_points();

    IdentityOperator identity;

    using FSspace= FunctionSpace< MeshT, Lagrange2<1>,RT0<1>>;
    FSspace FEspace(mesh);
    using FSspace2= FunctionSpace< MeshT, Lagrange2<1>,RT0<1>>;
    FSspace2 FEspace2(mesh);


 using fespace1= ElemFunctionSpace<Simplex<2,2>, Lagrange1<1>>;
 using fespace2= ElemFunctionSpace<Simplex<2,2>, Lagrange2<2>>;
 using fespace3= ElemFunctionSpace<Simplex<2,2>, RT0<1>>;
 std::cout<<"fespace"<<fespace2::FEFamily<<std::endl;
 std::cout<<"fespace"<<fespace2::Order<<std::endl;
 // BilinearFormIntegrator<4,FSspace,FSspace,GradientOperator,IdentityOperator>(FEspace,FEspace);

// std::shared_ptr<CollectorFunctionSpace> collfes=std::make_shared<CollectorFunctionSpace>(FEspace,FEspace);
// auto e1=collfes->get(0);
// const auto dm1=e1->dofmap();
const auto ee=FEspace.dofmap();



auto W=MixedFunctionSpace(FEspace,FEspace2);
auto W1=MixedFunctionSpace(W,FEspace,FEspace2);

std::cout<<FEspace.n_dofs()<<std::endl;
std::cout<<FEspace2.n_dofs()<<std::endl;
std::cout<<W.n_dofs()<<std::endl;
std::cout<<W1.n_dofs()<<std::endl;

auto& dmW =W.dofmap();
auto& dmW0 =std::get<0>(dmW);
auto& dmW01=std::get<1>(dmW);

auto& fe2map=FEspace2.dofmap();

auto& dm =W1.dofmap();
auto& dm0=std::get<0>(dm);
auto& dm01=std::get<0>(dm0);
auto& dm02=std::get<1>(dm0);
auto& dm1=std::get<1>(dm);
auto& dm2=std::get<2>(dm);

std::cout<<"dmW0"<<std::endl;
for(Integer ii=0;ii<dmW0.size();ii++)
{
   for(Integer jj=0;jj<dmW0[ii].size();jj++)
      std::cout<<dmW0[ii][jj]<<" ";
    std::cout<<std::endl;
}

std::cout<<"fe2map"<<std::endl;
for(Integer ii=0;ii<fe2map.size();ii++)
{
   for(Integer jj=0;jj<fe2map[ii].size();jj++)
      std::cout<<fe2map[ii][jj]<<" ";
    std::cout<<std::endl;
}

std::cout<<"dmW01"<<std::endl;
for(Integer ii=0;ii<dmW01.size();ii++)
{
   for(Integer jj=0;jj<dmW01[ii].size();jj++)
      std::cout<<dmW01[ii][jj]<<" ";
    std::cout<<std::endl;
}


std::cout<<"dm01"<<std::endl;
for(Integer ii=0;ii<dm01.size();ii++)
{
   for(Integer jj=0;jj<dm01[ii].size();jj++)
      std::cout<<dm01[ii][jj]<<" ";
    std::cout<<std::endl;
}
std::cout<<"dm02"<<std::endl;

for(Integer ii=0;ii<dm02.size();ii++)
{
   for(Integer jj=0;jj<dm02[ii].size();jj++)
      std::cout<<dm02[ii][jj]<<" ";
    std::cout<<std::endl;
}
std::cout<<"dm1"<<std::endl;

for(Integer ii=0;ii<dm1.size();ii++)
{
   for(Integer jj=0;jj<dm1[ii].size();jj++)
      std::cout<<dm1[ii][jj]<<" ";
    std::cout<<std::endl;
}
std::cout<<"dm2"<<std::endl;
for(Integer ii=0;ii<dm2.size();ii++)
{
   for(Integer jj=0;jj<dm2[ii].size();jj++)
      std::cout<<dm2[ii][jj]<<" ";
    std::cout<<std::endl;
}
auto TestW= Test(W);

};















void functionspaces_example2D()
{

    constexpr Integer ManifoldDim=2;
    constexpr Integer Dim=2;
    using MeshT=mars::Mesh<Dim, ManifoldDim>;
    MeshT mesh;
    using Elem = typename MeshT::Elem; 
    read_mesh("../data/beam-tri.MFEM", mesh);
 


     


    NodeToElem<Elem> node2elem3(mesh);
    auto node2elem=node2elem3.val();

    for(Integer nn=0; nn<node2elem.size();nn++)
        {std::cout<<"node=="<<nn<<"     ";
         for(Integer mm=0; mm<node2elem[nn].size();mm++)
          std::cout<<node2elem[nn][mm]<<" ";
         std::cout<<std::endl;
      }
    ElemEntity<Elem,0> nodes(mesh,node2elem);

    std::cout<<"node 2 elem size=="<<nodes.entity_2_elem().size()<<std::endl;
    std::cout<<"node 2 elem=="<<std::endl;
    for(Integer nn=0; nn<nodes.entity_2_elem().size();nn++)
         std::cout<<nodes.entity_2_elem(nn)[0]<<" ";
    std::cout<<std::endl;
    for(Integer nn=0; nn<nodes.elem_2_entity().size();nn++)
        {
         std::cout<<"elem="<<nn<< " made of nodes "<< std::endl;
         for(Integer mm=0; mm<nodes.elem_2_entity(nn).size();mm++)
            std::cout<<nodes.elem_2_entity(nn)[mm]<<" ";
         std::cout<< std::endl;
        } 


    using FSspace= FunctionSpace< MeshT, Lagrange2<1>,RT0<1>>;
    FSspace FEspace(mesh);


   const auto& P2_ens0=ElemEntity<Elem,ElementFunctionSpace<Elem,LagrangeFE,2>::entity[0]>(mesh,node2elem);
   const auto& P2_ens1=ElemEntity<Elem,ElementFunctionSpace<Elem,LagrangeFE,2>::entity[1]>(mesh,node2elem);

   auto ens2elem20=P2_ens0.entity_2_elem();
   auto ens2elem21=P2_ens1.entity_2_elem();

   auto elem2ens20=P2_ens0.elem_2_entity();
   auto elem2ens21=P2_ens1.elem_2_entity();

   std::cout<<"ens2elem 2 0="<< std::endl;
   for(Integer nn=0;nn<ens2elem20.size();nn++)
    {
        for(Integer mm=0;mm<ens2elem20[nn].size();mm++)
        std::cout<<ens2elem20[nn][mm]<<" ";
    std::cout<<std::endl;
    } 
   std::cout<<"ens2elem 2 1="<< std::endl;
   for(Integer nn=0;nn<ens2elem21.size();nn++)
    {
        for(Integer mm=0;mm<ens2elem21[nn].size();mm++)
        std::cout<<ens2elem21[nn][mm]<<" ";
    std::cout<<std::endl;
    } 
   std::cout<<"elem2ens20 2 0="<< std::endl;
   for(Integer nn=0;nn<elem2ens20.size();nn++)
    {
        for(Integer mm=0;mm<elem2ens20[nn].size();mm++)
        std::cout<<elem2ens20[nn][mm]<<" ";
    std::cout<<std::endl;
    } 
   std::cout<<"elem2ens21 2 1="<< std::endl;
   for(Integer nn=0;nn<elem2ens21.size();nn++)
    {
        for(Integer mm=0;mm<elem2ens21[nn].size();mm++)
        std::cout<<elem2ens21[nn][mm]<<" ";
    std::cout<<std::endl;
    } 
   
   std::cout<<"n_dofs="<<FEspace.n_dofs()<< std::endl;
    std::cout<<std::endl;
    for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
    {
     auto &elem_id=elem_iter;
     std::cout<<std::endl<<" elem =="<<elem_iter<<std::endl;
     for(Integer nn=0;nn<FEspace.dofmap(elem_id).size();nn++)
     {
        std::cout<<FEspace.dofmap(elem_id)[nn]<< "  ";
     }
    } 


   FEspace.set_new_start(4);
    std::cout<<"dofmap_new_start n_dofs="<<FEspace.n_dofs()<< std::endl;
    std::cout<<std::endl;
    for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
    {
     auto &elem_id=elem_iter;
     std::cout<<std::endl<<" elem =="<<elem_iter<<std::endl;
     for(Integer nn=0;nn<FEspace.dofmap_new_start(elem_id).size();nn++)
     {
        std::cout<<FEspace.dofmap_new_start(elem_id)[nn]<< "  ";
     }
    } 


std::cout<<std::endl;
for(Integer ss=0;ss<FEspace.n_subspaces();ss++)
       std::cout<<"components of space["<<ss<<"]=="<<FEspace.components(ss)<<std::endl;
std::cout<<std::endl;

for(Integer ss=0;ss<FEspace.n_subspaces();ss++)
{
   std::cout<<"dofs of space ="<<ss<<std::endl;
   for(Integer cc=0;cc<FEspace.components(ss);cc++)
    {std::cout<<"component ="<<cc<<"   "<<std::endl;
      auto& vec=FEspace.space_dofs(ss,cc);
      for(Integer mm=0;mm<FEspace.n_dofs(ss,cc);mm++)
          std::cout<<vec[mm]<<" ";
   }
   std::cout<<std::endl;

   std::cout<<"dofs of space ="<<ss<<std::endl;
   for(Integer cc=0;cc<FEspace.components(ss);cc++)
    {std::cout<<"component ="<<cc<<"   "<<std::endl;
      auto& vec=FEspace.space_dofs_new_start(ss,cc);
      for(Integer mm=0;mm<FEspace.n_dofs(ss,cc);mm++)
          std::cout<<vec[mm]<<" ";
   }
   std::cout<<std::endl;

}
 

    for(Integer mm=0;mm<FEspace.offset().size();mm++)
     {
      std::cout<<"OFFSET space ="<<mm<<std::endl;
      for(Integer nn=0;nn<FEspace.offset()[mm].size();nn++)
         {
            std::cout<< FEspace.offset()[mm][nn]<<" ";
         }
       std::cout<<std::endl;
     }
     std::cout<<std::endl;


 auto spaceinf=FEspace.space_info();
 for(Integer mm=0;mm<spaceinf.size();mm++)
     {
      std::cout<<"Space=="<<mm<<std::endl;
      for(Integer nn=0;nn<spaceinf[mm].size();nn++)
        std::cout<<spaceinf[mm][nn]<<" ";
      std::cout<<std::endl;
     }

 }


void normals_example3D()
{
    constexpr Integer ManifoldDim=3;
    constexpr Integer Dim=3;   
    using MeshT=mars::Mesh<Dim, ManifoldDim>;
    MeshT mesh;
    using Elem = typename MeshT::Elem;   
    read_mesh("../data/beam-tet.MFEM", mesh); 
    auto sn=SignedNormal<Simplex<Dim,ManifoldDim>>(mesh);
    auto n=sn();
    auto alpha=sn.sign();
    sn.print(mesh);
}

void normals_example4D()
{
    constexpr Integer ManifoldDim=4;
    constexpr Integer Dim=4;   
    using MeshT=mars::Mesh<Dim, ManifoldDim>;
    MeshT mesh;
    using Elem = typename MeshT::Elem;   
    read_mesh("../data/pentatope_2.MFEM", mesh);
    auto sn=SignedNormal<Simplex<Dim,ManifoldDim>>(mesh);
    auto n=sn();
    auto alpha=sn.sign();
    for(Integer e=0;e<mesh.n_elements();e++)
        {
           std::cout<<"elem=="<<e<<std::endl; 
           std::cout<<"normals:"<<std::endl; 
           for(Integer f=0;f<n[e].size();f++)
               n[e][f].describe(std::cout);
           std::cout<<"signs:"<<std::endl;
           std::cout<<alpha[e]<<std::endl; 
        }    
    sn.print(mesh);

}



void functionspaces_example4D()
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
   

    NodeToElem<Elem> node2elem3(mesh);
    auto node2elem=node2elem3.val();
    //auto vec=dofmap<ElemLagrange1<Elem>   , ElemLagrange3<Elem,1,1> >(mesh);

    // Primitive minusone(-1);
    // Primitive zero(0);
    // Composite first(1); 
    // Composite second(2); 
    // Composite third(3);
    // Composite fifth(5);  

    // first.add(zero);//std::make_shared<decltype(zero)>(zero));
    // third.add(minusone); 
    // third.add(zero); 
    // first.add(third);//std::make_shared<decltype(third)>(third));
    // second.add(first);
    // second.traverse(); 
    // fifth.add(minusone,zero,second,third);
    // std::cout<<std::endl<<"-----------------------"<<std::endl;
    // fifth.traverse();
    // std::cout<<std::endl<<"-----------------------"<<std::endl;
    // fifth[0]->traverse();
    // std::cout<<std::endl<<"-----------------------"<<std::endl;
    // fifth[1]->traverse();
    // std::cout<<std::endl<<"-----------------------"<<std::endl;
    // fifth[2]->traverse();
    // std::cout<<std::endl<<"-----------------------"<<std::endl;

    // auto fourth=second[0];
    // fourth.traverse(); 
    // cout<<endl<<" nnn = "<< *nnn<<endl;
     Leaf a0(3);
     Leaf a1(1);
     Leaf a2(2);
     Tree t1(-10);
     Tree t2(-20);
     a0.print();
     a1.print();
     a2.print();
     a2[0].print();
     t2.add(a0);
     t2.add(a1);

     t1.add(t2);
     t1.add(a0);
     t1.add(a1);
     t1.add(a2);
     
     auto ecco=t1[0][0];
     ecco.print();
    // ecco.print();
    // tree t1;
    
    // t1.add(a1);
    // auto ecc=t1[0];
    // std::cout<<"-----------------------"<<a1.val()<<"   "<<ecc->val()<<std::endl;

    // a1.print();
    // auto c1=t1[0];
    // std::cout<<"MA CI ARRIVO?"<<std::endl;
    // std::cout<<a1.val()<<" "<<a2.val()<<" "<<t3.val_vec().size()<<" "<<t3.val_vec()[0]->val()<<std::endl;
    // // constexpr auto n_spaces=2;
    // constexpr auto dofs_per_elem=DofsPerElemNums1<Elem,RT0<1>,Lagrange3<1>>::value;
    // std::vector<std::array<Integer,dofs_per_elem>> dofmap_vec;
    // std::array<std::vector<Integer>,n_spaces> offset;
    // dofmap1<RT0<1>,Lagrange3<1>>(mesh,dofmap_vec,offset);
    

    FunctionSpace< MeshT, Lagrange1<2>, Lagrange2<1> > FEspace(mesh);

    // auto eeh=FunctionSpaceSystem(FEspace);


   std::cout<<"n_dofs="<<FEspace.n_dofs()<< std::endl;
    std::cout<<std::endl;
    for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
    {
     auto &elem_id=elem_iter;
     std::cout<<"elem_id="<<elem_id<<", number of dofs=s"<<FEspace.dofmap(elem_id).size()<<std::endl;
     for(Integer nn=0;nn<FEspace.dofmap(elem_id).size();nn++)
     {
        std::cout<<FEspace.dofmap(elem_id)[nn]<<" ";
     }
     std::cout<<std::endl;
    } 
std::cout<<std::endl;
for(Integer ss=0;ss<FEspace.n_subspaces();ss++)
       std::cout<<"components of space["<<ss<<"]=="<<FEspace.components(ss)<<std::endl;
std::cout<<std::endl;

for(Integer ss=0;ss<FEspace.n_subspaces();ss++)
{
   std::cout<<"dofs of space ="<<ss<<std::endl;
   for(Integer cc=0;cc<FEspace.components(ss);cc++)
    {std::cout<<"component ="<<cc<<"   "<<std::endl;
      auto& vec=FEspace.space_dofs(ss,cc);
      for(Integer mm=0;mm<FEspace.n_dofs(ss,cc);mm++)
          std::cout<<vec[mm]<<" ";
      std::cout<<std::endl;
   }
   std::cout<<std::endl;

}
 

    for(Integer mm=0;mm<FEspace.offset().size();mm++)
     {
      std::cout<<"offset space="<<mm<<std::endl;
      for(Integer nn=0;nn<FEspace.offset()[mm].size();nn++)
         {
            std::cout<< FEspace.offset()[mm][nn]<<" ";
         }
     }
     std::cout<<std::endl;


   for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
    {std::cout<<std::endl;
      const auto size=FEspace.dofmap(1,elem_iter).size();
      std::cout<<"elem_iter="<<elem_iter<<std::endl;
      for(Integer nn=0;nn<size;nn++)
          std::cout<<FEspace.dofmap(1,elem_iter)[nn]<<" ";
      std::cout<<std::endl;
     }

 std::cout<<std::endl;
 }





void connectivity_example5D()
{

    constexpr Integer ManifoldDim=4;
    constexpr Integer Dim=4;
    using MeshT=mars::Mesh<Dim, ManifoldDim>;
    MeshT mesh;
    using Elem = typename MeshT::Elem; 
    read_mesh("../data/pentatope_2.MFEM", mesh);
    //read_mesh("../data/cube_6.MFEM", mesh);
    //read_mesh("../data/square_2_def.MFEM", mesh)    
    

    NodeToElem<Elem> node2elem3(mesh);
    const auto node2elem=node2elem3.val();
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




// Base Function: f(x) = x
// class Function
// {
//     protected:
//     struct Implementation
//     {
//         virtual ~Implementation() {}
//         virtual double evaluate(double x) const { return x; }
//     };

//     public:
//     Function()
//     :   self_(std::make_shared<Implementation>())
//     {}

//     double operator () (double x) const { return self_->evaluate(x); }

//     protected:
//     Function(std::shared_ptr<Implementation> self)
//     :   self_(self)
//     {}

//     private:
//     std::shared_ptr<Implementation> self_;
// };



// class Function2: public Function
// {
//      protected:
//     struct Implementation: Function::Implementation
//     {
//         Function f;
//         Implementation(Function f):   
//         f(f)
//         {};
//         virtual double evaluate(double x) const override{return 3.333*x; }
//     };  
// public: 
//     // double operator () (double x) const  { return 2*x; };
//     // Function2(Function f)
//     // :   Function(std::make_shared<Implementation>(f))
//     // {}
//     Function2()  
//     :   Function(std::make_shared<Implementation>(Function()))
//     {};
// };

// // Unary Function: u(-f(x))
// class UnaryMinus : public Function
// {
//     protected:
//     struct Implementation : Function::Implementation
//     {
//         Function f;


//         Implementation(Function f):   
//         f(f)
//         {};

//         virtual double evaluate(double x) const override { return -f(x); }
//     };

//     public:
//     UnaryMinus(Function f)
//     :   Function(std::make_shared<Implementation>(f))
//     {}
// };

// // Binary Function: u(f(x) + g(x))
// class BinaryAdd : public Function
// {
//     protected:
//     struct Implementation : Function::Implementation
//     {
//         Function f;
//         Function g;
//         Implementation(Function f, Function g)
//         :   f(f), g(g)
//         {};

//         virtual double evaluate(double x) const override { return f(x) + g(x); }
//     };

//     public:
//     BinaryAdd(Function f, Function g)
//     :   Function(std::make_shared<Implementation>(f, g))
//     {}
// };

// // Binary Function: u(f(x) * g(x))
// class BinaryMultiply : public Function
// {
//     protected:
//     struct Implementation : Function::Implementation
//     {
//         Function f;
//         Function g;
//         Implementation(Function f, Function g)
//         :   f(f), g(g)
//         {};

//         virtual double evaluate(double x) const override { return f(x) * g(x); }
//     };

//     public:
//     BinaryMultiply(Function f, Function g)
//     :   Function(std::make_shared<Implementation>(f, g))
//     {}
// };

// // Binary Function: u(f(x) * g(x))
// class scalarMultiply : public Function
// {
//     protected:
//     struct Implementation : Function::Implementation
//     {
//         Function f;
//         double alpha_;
//         Implementation(Function f, double alpha):   
//         f(f), 
//         alpha_(alpha)
//         {};


//         virtual double evaluate(double x) const override { return f(x) * alpha_; }
//     };

//     public:
//     scalarMultiply(Function f, double alpha)
//     :   Function(std::make_shared<Implementation>(f, alpha))
//     {}

// };

// inline scalarMultiply operator * (Function f,double alpha) { return scalarMultiply(f,alpha); }
// inline scalarMultiply operator * (double alpha,Function f) { return scalarMultiply(f,alpha); }
// inline UnaryMinus operator - (Function f) { return UnaryMinus(f); }
// inline BinaryAdd operator + (Function f, Function g) { return BinaryAdd(f, g); }
// inline BinaryMultiply operator * (Function f, Function g) { return BinaryMultiply(f, g); }
