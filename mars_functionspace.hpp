#ifndef MARS_FUNCTIONSPACE_HPP
#define MARS_FUNCTIONSPACE_HPP

#include "mars_base.hpp"
#include "mars_base_functionspace.hpp"

namespace mars{

template <Integer Dim,Integer ManifoldDim, Integer SpecificSpace,Integer Order>
class FunctionSpace : public BaseFunctionSpace<Dim,ManifoldDim,Lagrange, Order, 1 >
 {};



///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////------------------------------------ ///////////////////////////////
///////////////////////////////---------- LAGRANGE SPACE ---------- ///////////////////////////////
///////////////////////////////------------------------------------ /////////////////////////////// 
///////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////--------- LAGRANGE SPACE 1 --------- ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

template<Integer Dim>
using Lagrange1=FunctionSpace<Dim,Dim,Lagrange,1>;
using Lagrange1_1D=Lagrange1<1>;
using Lagrange1_2D=Lagrange1<2>;
using Lagrange1_3D=Lagrange1<3>;
using Lagrange1_4D=Lagrange1<4>;

///////////////////////////////------------- DIM == 1 --------------///////////////////////////////

template<>
class FunctionSpace<1,1,Lagrange,1>{
public:
     static constexpr Integer entity[]={0};
     static constexpr Integer entities_nums=1;
 };

constexpr Integer FunctionSpace<1,1,Lagrange,1>::entities_nums;
constexpr Integer FunctionSpace<1,1,Lagrange,1>::entity[FunctionSpace<1,1,Lagrange,1>::entities_nums];


///////////////////////////////------------- DIM == 2 --------------///////////////////////////////

template<>
class FunctionSpace<2,2,Lagrange,1>{
public:
     static constexpr Integer entity[]={0};
     static constexpr Integer entities_nums=1;
 };

constexpr Integer FunctionSpace<2,2,Lagrange,1>::entities_nums;
constexpr Integer FunctionSpace<2,2,Lagrange,1>::entity[FunctionSpace<2,2,Lagrange,1>::entities_nums];


///////////////////////////////------------- DIM == 3 --------------///////////////////////////////

template<>
class FunctionSpace<3,3,Lagrange,1>{
public:
     static constexpr Integer entity[]={0};
     static constexpr Integer entities_nums=1;
 };

constexpr Integer FunctionSpace<3,3,Lagrange,1>::entities_nums;
constexpr Integer FunctionSpace<3,3,Lagrange,1>::entity[FunctionSpace<3,3,Lagrange,1>::entities_nums];


///////////////////////////////------------- DIM == 4 --------------///////////////////////////////

template<>
class FunctionSpace<4,4,Lagrange,1>{
public:
     static constexpr Integer entity[]={0};
     static constexpr Integer entities_nums=1;
 };

constexpr Integer FunctionSpace<4,4,Lagrange,1>::entities_nums;
constexpr Integer FunctionSpace<4,4,Lagrange,1>::entity[FunctionSpace<4,4,Lagrange,1>::entities_nums];




///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////--------- LAGRANGE SPACE 2 --------- ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
template<Integer Dim>
using Lagrange2=FunctionSpace<Dim,Dim,Lagrange,2>;
using Lagrange2_1D=Lagrange2<1>;
using Lagrange2_2D=Lagrange2<2>;
using Lagrange2_3D=Lagrange2<3>;
using Lagrange2_4D=Lagrange2<4>;

///////////////////////////////------------- DIM == 2 --------------///////////////////////////////
template<>
class FunctionSpace<2,2,Lagrange,2>{
public:
     static constexpr Integer entity[]={0,1};
     static constexpr Integer entities_nums=2;
 };

constexpr Integer FunctionSpace<2,2,Lagrange,2>::entities_nums;
constexpr Integer FunctionSpace<2,2,Lagrange,2>::entity[FunctionSpace<2,2,Lagrange,2>::entities_nums];

///////////////////////////////------------- DIM == 3 --------------///////////////////////////////
template<>
class FunctionSpace<3,3,Lagrange,2>{
public:
     static constexpr Integer entity[]={0,1};
     static constexpr Integer entities_nums=2;
 };

constexpr Integer FunctionSpace<3,3,Lagrange,2>::entities_nums;
constexpr Integer FunctionSpace<3,3,Lagrange,2>::entity[FunctionSpace<3,3,Lagrange,2>::entities_nums];

///////////////////////////////------------- DIM == 4 --------------///////////////////////////////
template<>
class FunctionSpace<4,4,Lagrange,2>{
public:
     static constexpr Integer entity[]={0,1};
     static constexpr Integer entities_nums=2;
 };

constexpr Integer FunctionSpace<4,4,Lagrange,2>::entities_nums;
constexpr Integer FunctionSpace<4,4,Lagrange,2>::entity[FunctionSpace<4,4,Lagrange,2>::entities_nums];


///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////--------- LAGRANGE SPACE 3 --------- ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

template<Integer Dim>
using Lagrange3=FunctionSpace<Dim,Dim,Lagrange,3>;
using Lagrange3_1D=Lagrange3<1>;
using Lagrange3_2D=Lagrange3<2>;
using Lagrange3_3D=Lagrange3<3>;
using Lagrange3_4D=Lagrange3<4>;

///////////////////////////////------------- DIM == 2 --------------///////////////////////////////
template<>
class FunctionSpace<2,2,Lagrange,3>{
public:
     static constexpr Integer entity[]={0,1,2};
     static constexpr Integer entities_nums=3;
 };

constexpr Integer FunctionSpace<2,2,Lagrange,3>::entities_nums;
constexpr Integer FunctionSpace<2,2,Lagrange,3>::entity[FunctionSpace<2,2,Lagrange,3>::entities_nums];

///////////////////////////////------------- DIM == 3 --------------///////////////////////////////
template<>
class FunctionSpace<3,3,Lagrange,3>{
public:
     static constexpr Integer entity[]={0,1,2};
     static constexpr Integer entities_nums=3;
 };

constexpr Integer FunctionSpace<3,3,Lagrange,3>::entities_nums;
constexpr Integer FunctionSpace<3,3,Lagrange,3>::entity[FunctionSpace<3,3,Lagrange,3>::entities_nums];




///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////------------------------------------ ///////////////////////////////
///////////////////////////////------- RAVIART THOMAS SPACE ------- ///////////////////////////////
///////////////////////////////------------------------------------ /////////////////////////////// 
///////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////--------- RAVIART THOMAS 0 --------- ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
template<Integer Dim>
using RT0=FunctionSpace<Dim,Dim,RaviartThomas,0>;
using RT0_2D=RT0<2>;
using RT0_3D=RT0<3>;
using RT0_4D=RT0<4>;


///////////////////////////////------------- DIM == 2 --------------///////////////////////////////

template<>
class FunctionSpace<2,2,RaviartThomas,0>{
public:
     static constexpr Integer entity[]={1};
     static constexpr Integer entities_nums=1;};

constexpr Integer FunctionSpace<2,2,RaviartThomas,0>::entities_nums;
constexpr Integer FunctionSpace<2,2,RaviartThomas,0>::entity[FunctionSpace<2,2,RaviartThomas,0>::entities_nums];


///////////////////////////////------------- DIM == 3 --------------///////////////////////////////

template<>
class FunctionSpace<3,3,RaviartThomas,0>{
public:
     static constexpr Integer entity[]={2};
     static constexpr Integer entities_nums=1;};

constexpr Integer FunctionSpace<3,3,RaviartThomas,0>::entities_nums;
constexpr Integer FunctionSpace<3,3,RaviartThomas,0>::entity[FunctionSpace<3,3,RaviartThomas,0>::entities_nums];


///////////////////////////////------------- DIM == 4 --------------///////////////////////////////

template<>
class FunctionSpace<4,4,RaviartThomas,0>{
public:
     static constexpr Integer entity[]={3};
     static constexpr Integer entities_nums=1;};

constexpr Integer FunctionSpace<4,4,RaviartThomas,0>::entities_nums;
constexpr Integer FunctionSpace<4,4,RaviartThomas,0>::entity[FunctionSpace<4,4,RaviartThomas,0>::entities_nums];



///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////------------------------------------ ///////////////////////////////
///////////////////////////////----------- GENERAL SPACE ---------- ///////////////////////////////
///////////////////////////////------------------------------------ /////////////////////////////// 
///////////////////////////////////////////////////////////////////////////////////////////////////

template<>
class FunctionSpace<4,4,GeneralSpace,0>{
public:
     static constexpr Integer entity[]={0,1,2,3,4};
     static constexpr Integer entities_nums=5;};

constexpr Integer FunctionSpace<4,4,GeneralSpace,0>::entities_nums;
constexpr Integer FunctionSpace<4,4,GeneralSpace,0>::entity[FunctionSpace<4,4,GeneralSpace,0>::entities_nums];



}


#endif
