#ifndef MARS_FUNCTIONSPACE_HPP
#define MARS_FUNCTIONSPACE_HPP

#include "mars_base.hpp"
#include "mars_base_functionspace.hpp"
#include "mars_simplex.hpp"



namespace mars{

template <Integer Dim,Integer ManifoldDim, Integer FEFamily,Integer Order,Integer Continuity=Continuous,Integer NComponents=1>
class FunctionSpace : public BaseFunctionSpace<Dim,ManifoldDim,LagrangeFE, Order, Continuity,1 >
{};


template <typename Elem, Integer FEFamily, Integer Order,Integer Continuity=Continuous, Integer NComponents=1>
class ElementFunctionSpace: public BaseElementFunctionSpace<Elem,FEFamily,Order,Continuity,NComponents>
{};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////----------------------------------------------------------------------------/////////////////
///////////////////////////////---------- LAGRANGE SPACE ----- LAGRANGE SPACE ----- LAGRANGE SPACE --------/////////////////
///////////////////////////////----------------------------------------------------------------------------/////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////--------- LAGRANGE 1 -------- LAGRANGE 1 -------- LAGRANGE 1 --------/////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////--------- LAGRANGE 2 -------- LAGRANGE 2 -------- LAGRANGE 2 --------/////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////--------- LAGRANGE 3 -------- LAGRANGE 3 -------- LAGRANGE 3 --------/////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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







template<>
class ElementFunctionSpace<Simplex<4,4>,GeneralSpace,0,1,1>: 
      public BaseElementFunctionSpace<Simplex<4,4>,GeneralSpace,0,1,1>
{
public:
     static constexpr Integer entity[]={0,1,2,3,4};
     static constexpr Integer dofs_per_entity[]={1,1,1,1,1};
     static constexpr Integer entities_nums=5;
};


constexpr Integer ElementFunctionSpace<Simplex<4,4>,GeneralSpace,0,1,1>::entities_nums;
constexpr Integer ElementFunctionSpace<Simplex<4,4>,GeneralSpace,0,1,1>::dofs_per_entity[];
constexpr Integer ElementFunctionSpace<Simplex<4,4>,GeneralSpace,0,1,1>::entity[];








// template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
// class FunctionSpace<Dim,ManifoldDim,LagrangeFE,1,Continuity,NComponents>: public BaseFunctionSpace<Dim,ManifoldDim,LagrangeFE,1, Continuity,NComponents >
// {
// public:
//      static constexpr Integer entity[]={0};
//      static constexpr Integer dofs_per_entity[]={1};
//      static constexpr Integer entities_nums=1;
     
//  };

// template <Integer Dim,Integer ManifoldDim,Integer Continuity, Integer NComponents>
// constexpr Integer FunctionSpace<Dim,ManifoldDim,LagrangeFE,1,Continuity,NComponents>::entities_nums;
// template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
// constexpr Integer FunctionSpace<Dim,ManifoldDim,LagrangeFE,1,Continuity,NComponents>::dofs_per_entity[];
// template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
// constexpr Integer FunctionSpace<Dim,ManifoldDim,LagrangeFE,1,Continuity,NComponents>::entity[];





// template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
// class FunctionSpace<Dim,ManifoldDim,LagrangeFE,2,Continuity,NComponents>: public BaseFunctionSpace<Dim,ManifoldDim,LagrangeFE,2, Continuity,NComponents >
// {
// public:
//      static constexpr Integer entity[]={0,1};
//      static constexpr Integer dofs_per_entity[]={1,1};
//      static constexpr Integer entities_nums=2;
     
//  };

// template <Integer Dim,Integer ManifoldDim,Integer Continuity, Integer NComponents>
// constexpr Integer FunctionSpace<Dim,ManifoldDim,LagrangeFE,2,Continuity,NComponents>::entities_nums;
// template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
// constexpr Integer FunctionSpace<Dim,ManifoldDim,LagrangeFE,2,Continuity,NComponents>::dofs_per_entity[];
// template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
// constexpr Integer FunctionSpace<Dim,ManifoldDim,LagrangeFE,2,Continuity,NComponents>::entity[];



// template<Integer Dim, Integer ManifoldDim,Integer Continuity>
// class FunctionSpace<Dim,ManifoldDim,LagrangeFE,3,Continuity>: public BaseFunctionSpace<Dim,ManifoldDim,LagrangeFE,3, Continuity,1 >
// {
// public:
//      static constexpr Integer entity[]={0,1,2};
//      static constexpr Integer dofs_per_entity[]={1,2,1};     
//      static constexpr Integer entities_nums=3;
     
//  };

// template <Integer Dim,Integer ManifoldDim,Integer Continuity>
// constexpr Integer FunctionSpace<Dim,ManifoldDim,LagrangeFE,3,Continuity>::entities_nums;
// template<Integer Dim, Integer ManifoldDim,Integer Continuity>
// constexpr Integer FunctionSpace<Dim,ManifoldDim,LagrangeFE,3,Continuity>::dofs_per_entity[];
// template<Integer Dim, Integer ManifoldDim,Integer Continuity>
// constexpr Integer FunctionSpace<Dim,ManifoldDim,LagrangeFE,3,Continuity>::entity[];





// ///////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////------------------------------------ ///////////////////////////////
// ///////////////////////////////------- RAVIART THOMAS SPACE ------- ///////////////////////////////
// ///////////////////////////////------------------------------------ /////////////////////////////// 
// ///////////////////////////////////////////////////////////////////////////////////////////////////


// ///////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////--------- RAVIART THOMAS 0 --------- ///////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////////////////////



// template<Integer Dim, Integer ManifoldDim, Integer Continuity>
// class FunctionSpace<Dim,ManifoldDim,RaviartThomasFE,0,Continuity>: public BaseFunctionSpace<Dim,ManifoldDim,RaviartThomasFE,0 , Continuity,1 >
// {
// public:
//      static constexpr Integer entity[]={ManifoldDim-1};
//      static constexpr Integer dofs_per_entity[]={1};     
//      static constexpr Integer entities_nums=1;
     
//  };

// template <Integer Dim,Integer ManifoldDim, Integer Continuity>
// constexpr Integer FunctionSpace<Dim,ManifoldDim,RaviartThomasFE,0, Continuity>::entities_nums;
// template<Integer Dim, Integer ManifoldDim,  Integer Continuity>
// constexpr Integer FunctionSpace<Dim,ManifoldDim,RaviartThomasFE,0, Continuity>::dofs_per_entity[];
// template<Integer Dim, Integer ManifoldDim,  Integer Continuity>
// constexpr Integer FunctionSpace<Dim,ManifoldDim,RaviartThomasFE,0, Continuity>::entity[];




///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////------------------------------------ ///////////////////////////////
///////////////////////////////----------- GENERAL SPACE ---------- ///////////////////////////////
///////////////////////////////------------------------------------ /////////////////////////////// 
///////////////////////////////////////////////////////////////////////////////////////////////////



}




#endif
