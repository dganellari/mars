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
///////////////////////////////--------- LAGRANGE 0 -------- LAGRANGE 0 -------- LAGRANGE 0 --------/////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
class ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,0,Continuity,NComponents>
      : public BaseElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,0,Continuity,NComponents>
{
  public: 
     static constexpr const std::array<Integer,1> entity={ManifoldDim};
     static constexpr const std::array<Integer,1> dofs_per_entity={1};
     static constexpr const Integer ShapeFunctionDim1=1;
     static constexpr const Integer ShapeFunctionDim2=1;
};

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr std::array<Integer,1> ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,0,Continuity,NComponents>::entity;

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr std::array<Integer,1> ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,0,Continuity,NComponents>::dofs_per_entity;

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,0,Continuity,NComponents>::ShapeFunctionDim1;

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,0,Continuity,NComponents>::ShapeFunctionDim2;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////--------- LAGRANGE 1 -------- LAGRANGE 1 -------- LAGRANGE 1 --------/////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
class ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,1,Continuity,NComponents>
      : public BaseElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,1,Continuity,NComponents>
{
  public: 
     static constexpr const std::array<Integer,1> entity={0};
     static constexpr const std::array<Integer,1> dofs_per_entity={1};
     static constexpr const Integer ShapeFunctionDim1=1;
     static constexpr const Integer ShapeFunctionDim2=1;
};

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr std::array<Integer,1> ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,1,Continuity,NComponents>::entity;

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr std::array<Integer,1> ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,1,Continuity,NComponents>::dofs_per_entity;

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,1,Continuity,NComponents>::ShapeFunctionDim1;

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,1,Continuity,NComponents>::ShapeFunctionDim2;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////--------- LAGRANGE 2 -------- LAGRANGE 2 -------- LAGRANGE 2 --------/////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
class ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,2,Continuity,NComponents>
      : public BaseElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,2,Continuity,NComponents>
{
  public: 
     static constexpr const std::array<Integer,2> entity{0,1};
     static constexpr const std::array<Integer,2> dofs_per_entity{1,1};
     // static constexpr Integer entities_nums=entity.size();
     static constexpr const Integer ShapeFunctionDim1=1;
     static constexpr const Integer ShapeFunctionDim2=1;
};

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr std::array<Integer,2> ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,2,Continuity,NComponents>::entity;

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr std::array<Integer,2> ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,2,Continuity,NComponents>::dofs_per_entity;

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,2,Continuity,NComponents>::ShapeFunctionDim1;

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,2,Continuity,NComponents>::ShapeFunctionDim2;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////--------- LAGRANGE 3 -------- LAGRANGE 3 -------- LAGRANGE 3 --------/////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
class ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,3,Continuity,NComponents>
      : public BaseElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,3,Continuity,NComponents>
{
  public: 
     static constexpr const std::array<Integer,3> entity{0,1,2};
     static constexpr const std::array<Integer,3> dofs_per_entity{1,2,1};     
     static constexpr const Integer ShapeFunctionDim1=1;
     static constexpr const Integer ShapeFunctionDim2=1;
};

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr std::array<Integer,3> ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,3,Continuity,NComponents>::entity;

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr std::array<Integer,3> ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,3,Continuity,NComponents>::dofs_per_entity;

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,3,Continuity,NComponents>::ShapeFunctionDim1;

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer ElementFunctionSpace<Simplex<Dim, ManifoldDim>,LagrangeFE,3,Continuity,NComponents>::ShapeFunctionDim2;
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////------------------------------------ ///////////////////////////////
///////////////////////////////------- RAVIART THOMAS SPACE ------- ///////////////////////////////
///////////////////////////////------------------------------------ /////////////////////////////// 
///////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////
//////// For dimension of the Raviart-Thomas space, check:
//////// https://www-dimat.unipv.it/boffi/teaching/download/fulltext_duran.pdf
//////// S=Simplex of dimension n, k=order of the space 
//////// dim( RT_k(S) )=n * binomial_coeff(k+n,k)+binomial_coeff(k+n-1,k)
//////// dim( RT_0(S) ) = n + 1 
//////// dim( RT_1(S) ) = n * (1+n) + n =n (2+n)
///////////////////////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////--------- RAVIART THOMAS 0 --------- ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
template <Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
class ElementFunctionSpace<Simplex<Dim, ManifoldDim>,RaviartThomasFE,0,Continuity,NComponents>
      : public BaseElementFunctionSpace<Simplex<Dim, ManifoldDim>,RaviartThomasFE,0,Continuity,NComponents>
{
  public: 
     static constexpr const std::array<Integer,1> entity={ManifoldDim-1};
     static constexpr const std::array<Integer,1> dofs_per_entity={1};     
     static constexpr const Integer ShapeFunctionDim1=Dim;
     static constexpr const Integer ShapeFunctionDim2=1;
};

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr std::array<Integer,1> ElementFunctionSpace<Simplex<Dim, ManifoldDim>,RaviartThomasFE,0,Continuity,NComponents>::entity;

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr std::array<Integer,1> ElementFunctionSpace<Simplex<Dim, ManifoldDim>,RaviartThomasFE,0,Continuity,NComponents>::dofs_per_entity;

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer ElementFunctionSpace<Simplex<Dim, ManifoldDim>,RaviartThomasFE,0,Continuity,NComponents>::ShapeFunctionDim1;

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer ElementFunctionSpace<Simplex<Dim, ManifoldDim>,RaviartThomasFE,0,Continuity,NComponents>::ShapeFunctionDim2;

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////--------- RAVIART THOMAS 1 --------- ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
template <Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
class ElementFunctionSpace<Simplex<Dim, ManifoldDim>,RaviartThomasFE,1,Continuity,NComponents>
      : public BaseElementFunctionSpace<Simplex<Dim, ManifoldDim>,RaviartThomasFE,1,Continuity,NComponents>
{
  public: 
     static constexpr const std::array<Integer,2> entity={ManifoldDim-1,ManifoldDim};
     static constexpr const std::array<Integer,2> dofs_per_entity={ManifoldDim,ManifoldDim};     
     static constexpr const Integer ShapeFunctionDim1=Dim;
     static constexpr const Integer ShapeFunctionDim2=1;
};

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr std::array<Integer,2> ElementFunctionSpace<Simplex<Dim, ManifoldDim>,RaviartThomasFE,1,Continuity,NComponents>::entity;

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr std::array<Integer,2> ElementFunctionSpace<Simplex<Dim, ManifoldDim>,RaviartThomasFE,1,Continuity,NComponents>::dofs_per_entity;

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer ElementFunctionSpace<Simplex<Dim, ManifoldDim>,RaviartThomasFE,1,Continuity,NComponents>::ShapeFunctionDim1;

template<Integer Dim, Integer ManifoldDim,Integer Continuity, Integer NComponents>
constexpr Integer ElementFunctionSpace<Simplex<Dim, ManifoldDim>,RaviartThomasFE,1,Continuity,NComponents>::ShapeFunctionDim2;


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
     static constexpr const std::array<Integer,5> entity={0,1,2,3,4};
     static constexpr const std::array<Integer,5> dofs_per_entity={1,1,1,1,1};
     // static constexpr Integer entities_nums=entity.size();
};


// constexpr Integer ElementFunctionSpace<Simplex<4,4>,GeneralSpace,0,1,1>::entities_nums;
constexpr std::array<Integer,5> ElementFunctionSpace<Simplex<4,4>,GeneralSpace,0,1,1>::dofs_per_entity;
constexpr std::array<Integer,5> ElementFunctionSpace<Simplex<4,4>,GeneralSpace,0,1,1>::entity;
}


#endif
