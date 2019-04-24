#ifndef MARS_FUNCTIONSPACE_HPP
#define MARS_FUNCTIONSPACE_HPP

#include "mars_base.hpp"
#include "mars_base_functionspace.hpp"
#include "mars_simplex.hpp"



namespace mars{

template <Integer Dim,Integer ManifoldDim, Integer FEFamily,Integer Order,Integer Continuity=Continuous,Integer NComponents=1>
class FunctionSpace : public BaseFunctionSpace<Dim,ManifoldDim,LagrangeFE, Order, Continuity,1 >
{};


template <typename Elem, Integer FEFamily, Integer Order,Integer Continuity, Integer NComponents>
class ElementFunctionSpace: public BaseElementFunctionSpace<Elem,FEFamily,Order,Continuity,NComponents>
{};





template <Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
class ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,1,Continuity,NComponents>
      : public BaseElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,1,Continuity,NComponents>
{
  public: 
     static constexpr Integer entity[]={0};
     static constexpr Integer dofs_per_entity[]={1};
     static constexpr Integer entities_nums=1;  
};

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,1,Continuity,NComponents>::entity[];

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,1,Continuity,NComponents>::dofs_per_entity[];

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,1,Continuity,NComponents>::entities_nums;



template<typename Elem,Integer Continuity=Continuous, Integer NComponents=1>
using Lagrange1=ElementFunctionSpace<Elem,LagrangeFE,1,Continuity,NComponents>;













//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
class ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,2,Continuity,NComponents>
      : public BaseElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,2,Continuity,NComponents>
{
  public: 
     static constexpr Integer entity[]={0,1};
     static constexpr Integer dofs_per_entity[]={1,1};
     static constexpr Integer entities_nums=2;
};

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,2,Continuity,NComponents>::entity[];

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,2,Continuity,NComponents>::dofs_per_entity[];

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,2,Continuity,NComponents>::entities_nums;



template<typename Elem,Integer Continuity=Continuous, Integer NComponents=1>
using Lagrange2=ElementFunctionSpace<Elem,LagrangeFE,2,Continuity,NComponents>;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
class ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,3,Continuity,NComponents>
      : public BaseElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,3,Continuity,NComponents>
{
  public: 
     static constexpr Integer entity[]={0,1,2};
     static constexpr Integer dofs_per_entity[]={1,2,1};     
     static constexpr Integer entities_nums=3;
};

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,3,Continuity,NComponents>::entity[];

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,3,Continuity,NComponents>::dofs_per_entity[];

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,3,Continuity,NComponents>::entities_nums;



template<typename Elem,Integer Continuity=Continuous, Integer NComponents=1>
using Lagrange3=ElementFunctionSpace<Elem,LagrangeFE,3,Continuity,NComponents>;















// template <typename Elem,Integer Order,Integer Continuity, Integer NComponents>
// class Lagrange: public BaseElementFunctionSpace<Elem,LagrangeFE,Order,Continuity,NComponents>
// {};

// template <Integer Dim, Integer ManifoldDim,Integer Order,Integer Continuity, Integer NComponents>
// class Lagrange<Simplex<Dim, ManifoldDim>,Order,Continuity, NComponents> : public BaseElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,Order,Continuity,NComponents>
// {
// public:
//     // static const Integer space_dim=Dim; 
//     // static const Integer manifold_dim=ManifoldDim;

//      static constexpr Integer entity[]={0};
//      static constexpr Integer dofs_per_entity[]={1};
//      static constexpr Integer entities_nums=1;
// };

// template <Integer Dim,Integer ManifoldDim,Integer Order,Integer Continuity, Integer NComponents>
// constexpr Integer Lagrange<Simplex<Dim, ManifoldDim>,Order,Continuity, NComponents>::entities_nums;

// template<Integer Dim, Integer ManifoldDim,Integer Order,Integer Continuity, Integer NComponents>
// constexpr Integer Lagrange<Simplex<Dim, ManifoldDim>,Order,Continuity, NComponents>::dofs_per_entity[];

// template<Integer Dim, Integer ManifoldDim,Integer Order,Integer Continuity, Integer NComponents>
// constexpr Integer Lagrange<Simplex<Dim, ManifoldDim>,Order,Continuity, NComponents>::entity[];



// template <typename Elem,Integer Continuity, Integer NComponents>
// class Lagrange1p: public Lagrange<Elem,1,Continuity,NComponents>
// {};

// template <Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
// class Lagrange1p<Simplex<Dim, ManifoldDim>,Continuity,NComponents>: public Lagrange<Simplex<Dim, ManifoldDim>,1,Continuity,NComponents>
// {
// public:
//      static constexpr Integer entity[]={0};
//      static constexpr Integer dofs_per_entity[]={1};
//      static constexpr Integer entities_nums=1;
// };

// template <Integer Dim,Integer ManifoldDim,Integer Continuity, Integer NComponents>
// constexpr Integer Lagrange1p<Simplex<Dim, ManifoldDim>,Continuity, NComponents>::entities_nums;

// template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
// constexpr Integer Lagrange1p<Simplex<Dim, ManifoldDim>,Continuity, NComponents>::dofs_per_entity[];

// template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
// constexpr Integer Lagrange1p<Simplex<Dim, ManifoldDim>,Continuity, NComponents>::entity[];



// template <typename Elem,Integer Continuity, Integer NComponents>
// class Lagrange2p: public Lagrange<Elem,1,Continuity,NComponents>
// {};

// template <Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
// class Lagrange2p<Simplex<Dim, ManifoldDim>,Continuity,NComponents>: public Lagrange<Simplex<Dim, ManifoldDim>,2,Continuity,NComponents>
// {
// public:
//      static constexpr Integer entity[]={0,1};
//      static constexpr Integer dofs_per_entity[]={1,1};
//      static constexpr Integer entities_nums=2;
// };

// template <Integer Dim,Integer ManifoldDim,Integer Continuity, Integer NComponents>
// constexpr Integer Lagrange2p<Simplex<Dim, ManifoldDim>,Continuity, NComponents>::entities_nums;

// template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
// constexpr Integer Lagrange2p<Simplex<Dim, ManifoldDim>,Continuity, NComponents>::dofs_per_entity[];

// template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
// constexpr Integer Lagrange2p<Simplex<Dim, ManifoldDim>,Continuity, NComponents>::entity[];





// template <typename Elem,Integer Continuity, Integer NComponents>
// class Lagrange3p: public Lagrange<Elem,1,Continuity,NComponents>
// {};

// template <Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
// class Lagrange3p<Simplex<Dim, ManifoldDim>,Continuity,NComponents>: public Lagrange<Simplex<Dim, ManifoldDim>,3,Continuity,NComponents>
// {
// public:
//      static constexpr Integer entity[]={0,1,2};
//      static constexpr Integer dofs_per_entity[]={1,2,1};     
//      static constexpr Integer entities_nums=3;
// };

// template <Integer Dim,Integer ManifoldDim,Integer Continuity, Integer NComponents>
// constexpr Integer Lagrange3p<Simplex<Dim, ManifoldDim>,Continuity, NComponents>::entities_nums;

// template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
// constexpr Integer Lagrange3p<Simplex<Dim, ManifoldDim>,Continuity, NComponents>::dofs_per_entity[];

// template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
// constexpr Integer Lagrange3p<Simplex<Dim, ManifoldDim>,Continuity, NComponents>::entity[];




















///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////------------------------------------ ///////////////////////////////
///////////////////////////////---------- LAGRANGE SPACE ---------- ///////////////////////////////
///////////////////////////////------------------------------------ /////////////////////////////// 
///////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////--------- LAGRANGE SPACE 1 --------- ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

// template<Integer Dim,Integer Continuity,Integer NComponents=1>
// using Lagrange1=FunctionSpace<Dim,Dim,LagrangeFE,1,Continuity,NComponents>;
//////////////////////////////////////////////////
///////---- CONTINUOUS LAGRANGE 1 ----////////////
//////////////////////////////////////////////////
// template<Integer NComponents>
// using Lagrange1_1D=Lagrange1<1,Continuous,NComponents>;
// template<Integer NComponents>
// using Lagrange1_2D=Lagrange1<2,Continuous,NComponents>;
// template<Integer NComponents>
// using Lagrange1_3D=Lagrange1<3,Continuous,NComponents>;
// template<Integer NComponents>
// using Lagrange1_4D=Lagrange1<4,Continuous,NComponents>;


// //////////////////////////////////////////////////
// //////---- DISCONTINUOUS LAGRANGE 1 ----//////////
// //////////////////////////////////////////////////
// template<Integer NComponents>
// using Lagrange1_1D_DG=Lagrange1<1,Discontinuous,NComponents>;
// template<Integer NComponents>
// using Lagrange1_2D_DG=Lagrange1<2,Discontinuous,NComponents>;
// template<Integer NComponents>
// using Lagrange1_3D_DG=Lagrange1<3,Discontinuous,NComponents>;
// template<Integer NComponents>
// using Lagrange1_4D_DG=Lagrange1<4,Discontinuous,NComponents>;


template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
class FunctionSpace<Dim,ManifoldDim,LagrangeFE,1,Continuity,NComponents>: public BaseFunctionSpace<Dim,ManifoldDim,LagrangeFE,1, Continuity,NComponents >
{
public:
     static constexpr Integer entity[]={0};
     static constexpr Integer dofs_per_entity[]={1};
     static constexpr Integer entities_nums=1;
     
 };

template <Integer Dim,Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer FunctionSpace<Dim,ManifoldDim,LagrangeFE,1,Continuity,NComponents>::entities_nums;
template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer FunctionSpace<Dim,ManifoldDim,LagrangeFE,1,Continuity,NComponents>::dofs_per_entity[];
template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer FunctionSpace<Dim,ManifoldDim,LagrangeFE,1,Continuity,NComponents>::entity[];


///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////--------- LAGRANGE SPACE 2 --------- ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
// template<Integer Dim,Integer Continuity, Integer NComponents>
// using Lagrange2=FunctionSpace<Dim,Dim,LagrangeFE,2,Continuity, NComponents>;
// //////////////////////////////////////////////////
// ///////---- CONTINUOUS LAGRANGE 2 ----////////////
// //////////////////////////////////////////////////
// template<Integer NComponents>
// using Lagrange2_1D=Lagrange2<1,Continuous,NComponents>;
// template<Integer NComponents>
// using Lagrange2_2D=Lagrange2<2,Continuous,NComponents>;
// template<Integer NComponents>
// using Lagrange2_3D=Lagrange2<3,Continuous,NComponents>;
// template<Integer NComponents>
// using Lagrange2_4D=Lagrange2<4,Continuous,NComponents>;
// //////////////////////////////////////////////////
// //////---- DISCONTINUOUS LAGRANGE 2 ----//////////
// //////////////////////////////////////////////////
// template<Integer NComponents>
// using Lagrange2_1D_DG=Lagrange2<1,Discontinuous,NComponents>;
// template<Integer NComponents>
// using Lagrange2_2D_DG=Lagrange2<2,Discontinuous,NComponents>;
// template<Integer NComponents>
// using Lagrange2_3D_DG=Lagrange2<3,Discontinuous,NComponents>;
// template<Integer NComponents>
// using Lagrange2_4D_DG=Lagrange2<4,Discontinuous,NComponents>;



template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
class FunctionSpace<Dim,ManifoldDim,LagrangeFE,2,Continuity,NComponents>: public BaseFunctionSpace<Dim,ManifoldDim,LagrangeFE,2, Continuity,NComponents >
{
public:
     static constexpr Integer entity[]={0,1};
     static constexpr Integer dofs_per_entity[]={1,1};
     static constexpr Integer entities_nums=2;
     
 };

template <Integer Dim,Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer FunctionSpace<Dim,ManifoldDim,LagrangeFE,2,Continuity,NComponents>::entities_nums;
template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer FunctionSpace<Dim,ManifoldDim,LagrangeFE,2,Continuity,NComponents>::dofs_per_entity[];
template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer FunctionSpace<Dim,ManifoldDim,LagrangeFE,2,Continuity,NComponents>::entity[];


///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////--------- LAGRANGE SPACE 3 --------- ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

// template<Integer Dim,Integer Continuity>
// using Lagrange3=FunctionSpace<Dim,Dim,LagrangeFE,3,Continuity>;
// //////////////////////////////////////////////////
// ///////---- CONTINUOUS LAGRANGE 3 ----////////////
// //////////////////////////////////////////////////
// using Lagrange3_1D=Lagrange3<1,Continuous>;
// using Lagrange3_2D=Lagrange3<2,Continuous>;
// using Lagrange3_3D=Lagrange3<3,Continuous>;
// using Lagrange3_4D=Lagrange3<4,Continuous>;
// //////////////////////////////////////////////////
// //////---- DISCONTINUOUS LAGRANGE 3 ----//////////
// //////////////////////////////////////////////////
// using Lagrange3_1D_DG=Lagrange3<1,Discontinuous>;
// using Lagrange3_2D_DG=Lagrange3<2,Discontinuous>;
// using Lagrange3_3D_DG=Lagrange3<3,Discontinuous>;
// using Lagrange3_4D_DG=Lagrange3<4,Discontinuous>;


template<Integer Dim, Integer ManifoldDim,Integer Continuity>
class FunctionSpace<Dim,ManifoldDim,LagrangeFE,3,Continuity>: public BaseFunctionSpace<Dim,ManifoldDim,LagrangeFE,3, Continuity,1 >
{
public:
     static constexpr Integer entity[]={0,1,2};
     static constexpr Integer dofs_per_entity[]={1,2,1};     
     static constexpr Integer entities_nums=3;
     
 };

template <Integer Dim,Integer ManifoldDim,Integer Continuity>
constexpr Integer FunctionSpace<Dim,ManifoldDim,LagrangeFE,3,Continuity>::entities_nums;
template<Integer Dim, Integer ManifoldDim,Integer Continuity>
constexpr Integer FunctionSpace<Dim,ManifoldDim,LagrangeFE,3,Continuity>::dofs_per_entity[];
template<Integer Dim, Integer ManifoldDim,Integer Continuity>
constexpr Integer FunctionSpace<Dim,ManifoldDim,LagrangeFE,3,Continuity>::entity[];





// ///////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////------------------------------------ ///////////////////////////////
// ///////////////////////////////------- RAVIART THOMAS SPACE ------- ///////////////////////////////
// ///////////////////////////////------------------------------------ /////////////////////////////// 
// ///////////////////////////////////////////////////////////////////////////////////////////////////


// ///////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////--------- RAVIART THOMAS 0 --------- ///////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////////////////////
// template<Integer Dim, Integer Continuity>
// using RT0=FunctionSpace<Dim,Dim,RaviartThomasFE,0,Continuity>;
// //////////////////////////////////////////////////
// ///////---- CONTINUOUS RAVIAR-THOMAS 0 ----///////
// //////////////////////////////////////////////////
// using RT0_2D=RT0<2,Continuous>;
// using RT0_3D=RT0<3,Continuous>;
// using RT0_4D=RT0<4,Continuous>;
// //////////////////////////////////////////////////
// /////---- DISCONTINUOUS RAVIAR-THOMAS 0 ----//////
// //////////////////////////////////////////////////
// using RT0_2D_DG=RT0<2,Discontinuous>;
// using RT0_3D_DG=RT0<3,Discontinuous>;
// using RT0_4D_DG=RT0<4,Discontinuous>;


template<Integer Dim, Integer ManifoldDim, Integer Continuity>
class FunctionSpace<Dim,ManifoldDim,RaviartThomasFE,0,Continuity>: public BaseFunctionSpace<Dim,ManifoldDim,RaviartThomasFE,0 , Continuity,1 >
{
public:
     static constexpr Integer entity[]={ManifoldDim-1};
     static constexpr Integer dofs_per_entity[]={1};     
     static constexpr Integer entities_nums=1;
     
 };

template <Integer Dim,Integer ManifoldDim, Integer Continuity>
constexpr Integer FunctionSpace<Dim,ManifoldDim,RaviartThomasFE,0, Continuity>::entities_nums;
template<Integer Dim, Integer ManifoldDim,  Integer Continuity>
constexpr Integer FunctionSpace<Dim,ManifoldDim,RaviartThomasFE,0, Continuity>::dofs_per_entity[];
template<Integer Dim, Integer ManifoldDim,  Integer Continuity>
constexpr Integer FunctionSpace<Dim,ManifoldDim,RaviartThomasFE,0, Continuity>::entity[];




///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////------------------------------------ ///////////////////////////////
///////////////////////////////----------- GENERAL SPACE ---------- ///////////////////////////////
///////////////////////////////------------------------------------ /////////////////////////////// 
///////////////////////////////////////////////////////////////////////////////////////////////////

template<>
class FunctionSpace<4,4,GeneralSpace,0>: public BaseFunctionSpace<4,4,GeneralSpace,0,Continuous >
{
public:
     static constexpr Integer entity[]={0,1,2,3,4};
     static constexpr Integer dofs_per_entity[]={1,1,1,1,1};
     static constexpr Integer entities_nums=5;
};

// constexpr Integer FunctionSpace<4,4,GeneralSpace,0>::entities_nums;
// constexpr Integer FunctionSpace<4,4,GeneralSpace,0>::entity[FunctionSpace<4,4,GeneralSpace,0>::entities_nums];


constexpr Integer FunctionSpace<4,4,GeneralSpace,0>::entities_nums;
constexpr Integer FunctionSpace<4,4,GeneralSpace,0>::dofs_per_entity[];
constexpr Integer FunctionSpace<4,4,GeneralSpace,0>::entity[];



}




#endif
