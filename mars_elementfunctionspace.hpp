/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Written by Gabriele Rovi (April 2019)                                                                             ////////
////// We define the ElementFunctionSpace specialization for simplicial elements                                         ////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef MARS_FUNCTIONSPACE_HPP
#define MARS_FUNCTIONSPACE_HPP

#include "mars_base.hpp"
#include "mars_base_elementfunctionspace.hpp"
#include "mars_simplex.hpp"



namespace mars{
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

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////------------------------------------ ///////////////////////////////
///////////////////////////////------- RAVIART THOMAS SPACE ------- ///////////////////////////////
///////////////////////////////------------------------------------ /////////////////////////////// 
///////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////--------- RAVIART THOMAS 0 --------- ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
template <Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
class ElementFunctionSpace<Simplex<Dim, ManifoldDim>,RaviartThomasFE,0,Continuity,NComponents>
      : public BaseElementFunctionSpace<Simplex<Dim, ManifoldDim>,RaviartThomasFE,0,Continuity,NComponents>
{
  public: 
     static constexpr Integer entity[]={ManifoldDim-1};
     static constexpr Integer dofs_per_entity[]={1};     
     static constexpr Integer entities_nums=1; 
};

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer ElementFunctionSpace<Simplex<Dim, ManifoldDim>,RaviartThomasFE,0,Continuity,NComponents>::entity[];

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer ElementFunctionSpace<Simplex<Dim, ManifoldDim>,RaviartThomasFE,0,Continuity,NComponents>::dofs_per_entity[];

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer ElementFunctionSpace<Simplex<Dim, ManifoldDim>,RaviartThomasFE,0,Continuity,NComponents>::entities_nums;

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////------------------------------------ ///////////////////////////////
///////////////////////////////--------   GENERAL SPACE ----------  ///////////////////////////////
///////////////////////////////------------------------------------ /////////////////////////////// 
///////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////--------   GENERAL SPACE 0 ---------  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
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
}


#endif
