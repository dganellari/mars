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

Elem elem;
Matrix<Real, Dim, ManifoldDim> J;

constexpr Integer Npoints=Elem::Npoints;
std::vector<Vector<Real,Dim>> points(Npoints);

for(Integer ii=0;ii<Npoints;ii++)
   elem.nodes[ii]=ii;

points[0][0]=0;
points[0][1]=0;
points[1][0]=2;
points[1][1]=3;
points[2][0]=1;
points[2][1]=4;
  jacobian(elem,points,J);

    using FSspace1= FunctionSpace< MeshT, Lagrange2<1>,RT0<1>>;
    FSspace1 FEspace1(mesh);
    using FSspace2= FunctionSpace< MeshT, Lagrange2<1>,RT0<1>>;
    FSspace2 FEspace2(mesh);
    using FSspace3= FunctionSpace< MeshT, Lagrange1<1>, RT0<1>,Lagrange1<1>>;
    FSspace3 FEspace3(mesh);
    using FSspace4= FunctionSpace< MeshT, Lagrange2<1>, Lagrange3<1>>;
    FSspace4 FEspace4(mesh);
    using FSspace5= FunctionSpace< MeshT, Lagrange2<1>, Lagrange3<2>>;

    using FSspace10= FunctionSpace< MeshT, Lagrange1<1>>;
    FSspace10 FEspace10(mesh);
      using FSspace11= FunctionSpace< MeshT, RT0<1>>;
    FSspace11 FEspace11(mesh);
// auto ecchime=FEspace3.set_quadrature_rule<0,GaussPoints<Elem,2>,GaussPoints<Elem,2>>();
 using fespace1= ElemFunctionSpace<Simplex<2,2>, Lagrange1<1>>;
 using fespace2= ElemFunctionSpace<Simplex<2,2>, Lagrange2<2>>;
 using fespace3= ElemFunctionSpace<Simplex<2,2>, RT0<1>>;
 
 BilinearFormIntegrator<4,FSspace1,FSspace1,GradientOperator,IdentityOperator>(FEspace1,FEspace1);


//  using FESpace4=FESpace<MeshT,Lagrange1<2>,GaussQP<4>,Lagrange1<1>,GaussQP<4>,Lagrange1<3>,GaussQP<4>>;
//  FESpace4 FE4(mesh);
//   std::cout<<"------------------------------------------------------------"<<std::endl;
//   std::cout<<"Nelem_dofs==="<<FESpace4::Nelem_dofs<<std::endl;
//   std::cout<<"------------------------------------------------------------"<<std::endl;

// // std::shared_ptr<CollectorFunctionSpace> collfes=std::make_shared<CollectorFunctionSpace>(FEspace,FEspace);
// // auto e1=collfes->get(0);
// // const auto dm1=e1->dofmap();
// const auto ee=FEspace.dofmap();


// // auto Wp=MixedFunctionSpace(FE3,FE3);

auto W=MixedFunctionSpace(FEspace3,FEspace4);
auto W1=MixedFunctionSpace(W,FEspace4);


constexpr Integer maximum=Max(0,2) ;

typename decltype(W)::type_unique_base_function_space unique0;
GetType<1,decltype(W)::type_unique_base_function_space> unique1;
decltype(W1)::type_unique_base_function_space unique2;



std::cout<<ElementPosition<0,std::tuple<int,int,double,int>,std::tuple<int,double>>::value<<std::endl;
std::cout<<ElementPosition<1,std::tuple<int,int,double,int>,std::tuple<int,double>>::value<<std::endl;
std::cout<<ElementPosition<2,std::tuple<int,int,double,int>,std::tuple<int,double>>::value<<std::endl;
std::cout<<ElementPosition<maximum,std::tuple<int,int,double,int>,std::tuple<int,double>>::value<<std::endl;
std::cout<<"TypeToTupleElementPosition=="<<TypeToTupleElementPosition<double,std::tuple<int,char,double,int,double>>::value<<std::endl;
std::cout<<maximum<<std::endl;

constexpr Integer NComponents=2;

using QR=GaussPoints<Elem,QPOrder>;
ShapeFunctionOperatorDependent<Simplex<2,2>, Lagrange1<NComponents>,IdentityOperator,QR >  sfod;
sfod.init_reference();
sfod.init(J);
sfod.function();
std::cout<<"sfod.function()=="<<sfod.function()<<std::endl;








auto W3=MixedFunctionSpace(FEspace10,FEspace11);


ShapeFunctionExpression<0,2,Lagrange1<1>> sfe1; 
GradientShapeFunctionExpression<0,2,Lagrange1<1>> sfe2;
ShapeFunctionExpression<1,2,RT0<1>> sfe3;


auto sfen=sfe1+sfe2+sfe3;

using TupleOperatorsAndQuadrature=TupleOfTupleChangeType<1,QR,OperatorTupleType<decltype(sfen)>::type>; 
using TupleSpaces=typename decltype(W3)::type_unique_base_function_space;

// OperatorTupleType<decltype(sfen)>::type e0(3);
// TupleOperatorsAndQuadrature e(4);
// typename TupleShapeFunctionCreate<Elem0,BaseFunctioSpace0,OperatorsAndQuadrature0,TupleTypeSize<OperatorsAndQuadrature0>::value-1>::type ef2(3);
// typename TupleShapeFunctionCreate<Elem1,BaseFunctioSpace1,OperatorsAndQuadrature1,TupleTypeSize<OperatorsAndQuadrature1>::value-1>::type ef3(3);

TupleOfTupleShapeFunctionType<TupleSpaces,TupleOperatorsAndQuadrature> ee;

ShapeFunctionAssembly<TupleSpaces,TupleOperatorsAndQuadrature> ee2;
ee2.init_reference();
ee2.init(J);
const auto& ee3=ee2.get<0,1>();
std::cout<<"examples reference"<<ee3.reference()<<std::endl;
std::cout<<"examples function"<<ee3.function()<<std::endl;
using All=std::tuple<int,char,double,int,char>;
using Unique=std::tuple<int,char,double>;
std::cout<<"examples value="<<ElementPosition<0,All,Unique>::value<<std::endl;
std::cout<<"examples value="<<ElementPosition<1,All,Unique>::value<<std::endl;
std::cout<<"examples value="<<ElementPosition<2,All,Unique>::value<<std::endl;
std::cout<<"examples value="<<ElementPosition<3,All,Unique>::value<<std::endl;
std::cout<<"examples value="<<ElementPosition<4,All,Unique>::value<<std::endl;

TupleAllToUniqueMap<All,Unique> tup;
const auto tup0=std::get<0>(tup);
std::cout<<"examples value="<<tup0.value<<std::endl;
auto& mm=W1.space_avatar();

// TupleTrialOrTest<decltype(W)::type_tuple_spaces,0> eeee(3);
// TupleTrialOrTest<decltype(W)::type_tuple_spaces,1> eeee2(3);


auto u =     MakeTrial<0>(W1);
auto sigma = MakeTrial<1>(W1);
auto p =     MakeTrial<2>(W1);
auto r =     MakeTrial<3>(W1);
auto rr1 =     MakeTrial<4>(W1);
auto rr2 =     MakeTrial<5>(W1);
auto rr3 =     MakeTrial<6>(W1);

auto v =   MakeTest<0>(W1);
auto tau = MakeTest<1>(W1);
auto q =   MakeTest<2>(W1);
auto s =   MakeTest<3>(W1);
auto ss1 =   MakeTest<4>(W1);
auto ss2 =   MakeTest<5>(W1);
auto ss3 =   MakeTest<6>(W1);

constexpr bool bool1[]={true,false,true};
constexpr bool bool2[]={true,false,true};
constexpr bool bool3[]={Max(bool1[0],bool2[0])};


auto l20=L2Inner(mesh,Div(sigma),Div(tau))+L2Inner(mesh,Div(sigma),Div(tau));
auto l21=L2Inner(mesh,Grad(u),Grad(v));
auto l22=L2Inner(mesh,Div(sigma),Div(tau))+
         L2Inner(mesh,Grad(u),Grad(v))+
         L2Inner(mesh,u,Div(tau))+
         L2Inner(mesh,rr1,q)+
         L2Inner(mesh,r,s)+
         L2Inner(mesh,Grad(rr1)+sigma,tau);

auto l23=L2Inner(mesh,Grad(u)+u,Grad(v));
auto l24=L2Inner(mesh,Div(sigma)+u,Div(tau));
auto l25=L2Inner(mesh,Grad(u)+sigma,Grad(v));
auto l26=L2Inner(mesh,Grad(u)+sigma,Grad(v)+v);
auto l27=L2Inner(mesh,Grad(u)+u+sigma+v,Grad(v)+Div(sigma)+v+u+tau);


// std::cout<<"tensorvector size="<<SupportOverlap<TensorVector1,1,TensorVector2,1,0,0>()<<std::endl;
auto s3 =   MakeTest<bool1[0]+bool2[0]>(W1);
std::cout<<"tensorvector size="<<Overlap<TensorVector1,TensorVector2>()<<std::endl;
// std::cout<<"bbilinear order="<<l20.Order<<std::endl;
std::cout<<"bbilinear order="<<l21.Order<<std::endl;
// std::cout<<"bbilinear order="<<l22.Order<<std::endl;
std::cout<<"bbilinear order="<<l23.Order<<std::endl;
std::cout<<"bbilinear order="<<l24.Order<<std::endl;
std::cout<<"bbilinear order="<<l25.Order<<std::endl;
std::cout<<"bbilinear order="<<l26.Order<<std::endl;
std::cout<<"bbilinear order="<<l27.Order<<std::endl;
std::cout<<"decltype(Grad(u))::Nmax="<<decltype(Grad(u))::Nmax<<std::endl;
std::cout<<"decltype(Grad(u))::N="<<decltype(Grad(u))::value<<std::endl;
std::cout<<"decltype(Grad(sigma))::N="<<decltype(Grad(sigma))::value<<std::endl;
std::cout<<"decltype(Grad(p))::N="<<decltype(Grad(p))::value<<std::endl;
std::cout<<"decltype(Grad(r))::N="<<decltype(Grad(r))::value<<std::endl;
std::cout<<"decltype(Grad(rr1))::N="<<decltype(Grad(rr1))::value<<std::endl;
std::cout<<"decltype(Grad(rr2))::N="<<decltype(Grad(rr2))::value<<std::endl;
std::cout<<"decltype(Grad(rr3))::N="<<decltype(Grad(rr3))::value<<std::endl;
std::cout<<" TupleTypeSize< decltype(W1)::type_unique_base_function_space>::value"<< TupleTypeSize< decltype(W1)::type_unique_base_function_space>::value<<std::endl;
std::cout<<"decltype(Grad(u))::Nmax="<<decltype(Grad(u))::Nmax<<std::endl;
std::cout<<"decltype(Grad(u))::N="<<decltype(Grad(u))::value<<std::endl;
std::cout<<"decltype(Grad(u))::N="<<decltype(Grad(sigma))::value<<std::endl;
std::cout<<"decltype(Grad(u))::N="<<decltype(Grad(p))::Nmax<<std::endl;
std::cout<<"decltype(Grad(u))::N="<<decltype(Grad(r))::Nmax<<std::endl;
std::cout<<"decltype(Grad(u))::N="<<decltype(Grad(rr1))::Nmax<<std::endl;
std::cout<<"decltype(Grad(u))::N="<<decltype(Grad(rr2))::Nmax<<std::endl;
std::cout<<"decltype(Grad(u))::N="<<decltype(Grad(rr3))::Nmax<<std::endl;


 using TupleSpacesW1=typename decltype(W1)::type_unique_base_function_space;
 using TupleOperatorsAndQuadratureW1= typename OperatorTupleType<decltype(l22)>::type;
 using TupleOfTupleNoQuadrature=TupleOfTupleRemoveQuadrature<TupleOperatorsAndQuadratureW1>;
 TupleOfTupleNoQuadrature feef;


TupleOfTupleNoQuadrature deel;
TupleSpacesW1 r3r;
MapOperatorTupleOfTuple<TupleOfTupleNoQuadrature,TupleSpacesW1> eeded;
TupleOfTupleNoQuadrature frfr;
UniqueFEFamilies<TupleSpacesW1> ufe;
SpacesToUniqueFEFamilies<TupleSpacesW1,UniqueFEFamilies<TupleSpacesW1>> mce;
// typename ShapeFunctionAssembly<TupleSpacesW1,TupleOperatorsAndQuadratureW1>::type fghi(4);




std::cout<<"MINCHIA="<<minchia()<<std::endl;



using tp5=std::tuple< Number<1>,Number<0>,Number<0>,Number<0>,Number<1>,Number<0>,Number<0>,Number<0>,Number<0>>;
using tp6=std::tuple< Number<1>,Number<0>,Number<0>,Number<0>,Number<1>,Number<0>,Number<0>,Number<0>,Number<0>>;
using nonzeroproduct=StaticBooleanMatrixMatrixMultiplicationFindNonZeros<3,3,tp5,3,3,tp6>;
using booleanmatrix=TupleOfTupleToBooleanTuple<nonzeroproduct>;
using findnonzeronumbers=FindNonZeroNumbers<booleanmatrix>;
Matrix<Real,3,3> mat1{1,2,3,4,5,6,7,8,9};
Matrix<Real,3,3> mat2{2,3,4,5,6,7,8,9,10};
Matrix<Real,3,3> mat3;//{0,0,0,0,0,0,0,0,0};
Matrix<Real,3,3> mat4;//{0,0,0,0,0,0,0,0,0};


constexpr Integer Rows=3;
constexpr Integer Cols=3;

using boolean1=std::tuple<One,Zero,Zero,One,One,Zero,Zero,Zero,One>;
StaticMatrix<Real,3,3,boolean1> ede;
StaticMatrix<Real,3,3> hfe;

Real nonzero;
Integer Constant=pow(10,7);




constexpr std::array<Integer,4> arr{1,2,3,4};
constexpr Integer b=arr[0];

Matrix<Real,2,2> mat5{2,2,2,2};
Matrix<Real,2,2> mat7;

constexpr std::array<Real, 4> a6{1,1,1,1};
static_assert(a6[0]==1,"erfrr");


std::cout<<"mat5="<<mat5<<std::endl;


ShapeFunctionDependent<Simplex<2,2>,Lagrange1<2>,IdentityOperator,GaussPoints<Elem,3>> sfd;

std::cout<<sfd.reference()<<std::endl;

constexpr const QPValues<int,2> qpc{1,2};
constexpr const FQPValues<int,1,2> fqpc{{{1},{2}}};
constexpr const Vector<Vector<Real,1>,2> fqpc2{{1.0},{2.0}};
// ShapeFunctionDependent<Simplex<2,2>,Lagrange1<2>,IdentityOperator,GaussPoints<Elem,3>> sfd4;

// constexpr auto romolo=ShapeFunctionDependent<Simplex<2,2>,Lagrange1<2>,IdentityOperator,GaussPoints<Elem,3>>::init_reference(GaussPoints<Elem,3>::qp_points);

// Matrix<double,2,2> dynamicmat3{0.1,0.2,0.3,0.4};

// static_assert(staticmat(0,0)==0.1," is not 2");
// static_assert(staticmat2(0,0)==0.1," is not 2");
// assert(dynamicmat3(0,0)==0.1 &&" is not 2");


// Matrix<Real, 2, 2> jac3{0.1,2,3,5};
// static constexpr Matrix<Real,1,1> f{0.1};
// static constexpr Matrix<Matrix<Real,1,1>, 2, 2> jk4{0.1,0.1,0.1,0.1};//,f,f,f};

// Matrix< Vector<Real,5>, 2, 2> jk5{0.1,0.1,0.1,0.1};
// // Matrix< Matrix<Real,1,1>, 2, 2> jk7{f,f,f,{0.1}};

//  MAtrixProva3< MAtrixProva3<Real,1,1>, 2, 2> jk6{ {0.1},{0.1},{0.1},{0.1} };


//  constexpr MAtrixProva3< Real, 2, 2> jac4{1,1,2,3};
// jacobian(elem,points,jac3);


   static constexpr Matrix<double,2,2> staticmat{0.1,0.2,0.3,0.4};
   Matrix<double,2,2> jac{0.1,0.2,0.3,0.4};

   constexpr int zero=0;
   constexpr auto staticrow=staticmat.get_row(zero);
   
   constexpr auto static00=staticmat(0,0);
   Matrix<Matrix<double,1,1>, 2, 2> a4{{0.1}, {0.1}, {0.1}, {0.1}};
   Matrix<Matrix<double,1,1>, 2, 2> b4{0.1, 0.1, 0.1, 0.1};

   static constexpr Vector<int,2> staticvec{1,2};
   jacobian(elem,points,jac);

   constexpr double somma=staticvec.sum();
   static_assert(somma==3,"wrong static vector sum");
   static_assert(QPSimplex1_qp_points(0,1)==0.33333333333333,"wrong static vector sum");

   static_assert(GaussPoints< Simplex<Dim,2> , 1>::qp_points(0,0)==0.33333333333333,"wroeeeeeeeeng static vector sum");
   static_assert(GaussPoints< Simplex<Dim,2> , 1>::qp_points(0,1)==0.33333333333333,"wroeeeeeeeeng static vector sum");
std::cout<<"points="<<points[0]<<std::endl;
std::cout<<"points="<<points[1]<<std::endl;
std::cout<<"GaussPoints="<<QPSimplex1_qp_points(0,0)<<std::endl;
std::cout<<"GaussPoints="<<QPSimplex1_qp_points(0,1)<<std::endl;
std::cout<<"GaussPoints="<<GaussPoints< Simplex<Dim,2> , 1>::qp_points(0,0)<<std::endl;
std::cout<<"GaussPoints="<<GaussPoints< Simplex<Dim,2> , 1>::qp_points(0,1)<<std::endl;
// std::cout<<"GaussPoints="<<GaussPoints< Simplex<Dim,2> , 1>::qp_points(0,1)<<std::endl;

constexpr auto mmm=staticmat;
 clock_t begin = clock();
  for(Integer ii=0;ii<Constant;ii++)
  {
    //     for(Integer i = 0; i < 3; ++i) {
    //       for(Integer j = 0; j < 3; ++j) {
    //           mat1(i, j) = mat1(i, j) /pow(1.0,1.0/ii);
    //         }
    //       }
        

    // StaticMatrixProduct<findnonzeronumbers,nonzeroproduct>(mat1,mat2,mat3);
    // StaticMatrixProduct<findnonzeronumbers,nonzeroproduct>(mat1,mat3,mat4);
    // mat3=mat1*mat2;
    

        // mat3=0;
        // for(Integer i = 0; i < 3; ++i) {
        //   for(Integer j = 0; j < 3; ++j) {
        //     for(Integer k = 0; k < 3; ++k) {
        //       mat3(i, k) += mat1(i, j) * mat2(j, k);
        //     }
        //   }
        // }

//         for(Integer i = 0; i < 3; ++i) {
//           for(Integer j = 0; j < 3; ++j) {
//             mat3(i, j) = mat1(i, 0) * mat2(0,j);
//             for(Integer k = 1; k < 1; ++k) {
//               mat3(i, j) += mat1(i, k) * mat2(k,j);
//             }
//           }
//         }


// nonzero+=mat3(0,0)+mat3(0,1)+mat3(0,2)+mat3(1,0)+mat3(1,1)+mat3(1,2)+mat3(2,0)+mat3(2,1)+mat3(2,2);

   // nonzero=StaticScalarProduct<Nonzero>(v1,v2);
    // StaticMatrixProduct<findnonzeronumbers,nonzeroproduct>(mat1,mat2,mat3,vec1,vec2);
    //     // nonzero=0;


              // mat3= mat1 * mat2;

     // nonzero=v1[3]*v2[2]+v1[4]*v2[5]+v1[5]*v2[8];
    // nonzero=0;
    // for(Integer ii=0;ii<3;ii++)
    //   nonzero+=mat1(1,ii)*mat2(ii,2);
  }
  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

std::cout<<"mat==="<<mat3<<std::endl;

std::cout<<" elapsed_secs------------>"<<elapsed_secs<<std::endl;

std::cout<<" nonzero------------>"<<nonzero<<std::endl;

auto& dm =W1.dofmap();
auto& dm0=std::get<0>(dm);
auto& dm01=std::get<0>(dm0);
auto& dm02=std::get<1>(dm0);
auto& dm1=std::get<1>(dm);
// auto& dm2=std::get<2>(dm);




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
// std::cout<<"dm2"<<std::endl;
// for(Integer ii=0;ii<dm2.size();ii++)
// {
//    for(Integer jj=0;jj<dm2[ii].size();jj++)
//       std::cout<<dm2[ii][jj]<<" ";
//     std::cout<<std::endl;
// }





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


   // FEspace.set_new_start(4);
   //  std::cout<<"dofmap_new_start n_dofs="<<FEspace.n_dofs()<< std::endl;
   //  std::cout<<std::endl;
   //  for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
   //  {
   //   auto &elem_id=elem_iter;
   //   std::cout<<std::endl<<" elem =="<<elem_iter<<std::endl;
   //   for(Integer nn=0;nn<FEspace.dofmap_new_start(elem_id).size();nn++)
   //   {
   //      std::cout<<FEspace.dofmap_new_start(elem_id)[nn]<< "  ";
   //   }
   //  } 


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
     // Leaf a0(3);
     // Leaf a1(1);
     // Leaf a2(2);
     // Tree t1(-10);
     // Tree t2(-20);
     // a0.print();
     // a1.print();
     // a2.print();
     // a2[0].print();
     // t2.add(a0);
     // t2.add(a1);

     // t1.add(t2);
     // t1.add(a0);
     // t1.add(a1);
     // t1.add(a2);
     
     // auto ecco=t1[0][0];
     // ecco.print();
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
