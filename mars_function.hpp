#ifndef MARS_FUNCTION_HPP
#define MARS_FUNCTION_HPP

#include "mars_base.hpp"


namespace mars{


class Function1
{
    public: 
    using Point=Vector<Real,2>;
    using Output=Real;

    static Output eval(const Point& p)
    {
     return p[0]+p[1]; 
    }
};



// template<typename FuncType, typename FunctionSpace,Integer N_=0, typename TupleSpaces=FunctionSpace>
// class Function: public Expression<Function<FuncType,FunctionSpace,N_,TupleSpaces>>
// {
  
//   public:

//     using ElementFunctionSpace=GetType<typename FunctionSpace::TupleOfSpaces,0>;
//     using BaseFunctionSpace=Elem2FunctionSpace<ElementFunctionSpace>;
//     static constexpr Integer FEFamily=ElementFunctionSpace::FEFamily;
//     static constexpr Integer Order=ElementFunctionSpace::Order;
//     static constexpr Integer Continuity=ElementFunctionSpace::Continuity;
//     static constexpr Integer NComponents=ElementFunctionSpace::NComponents;
//     using Elem=typename ElementFunctionSpace::Elem;
//     static constexpr Integer Dim=Elem::Dim;
//     using sub_type=typename FuncType::Output;
//     using Operator=IdentityOperator;
//     static constexpr Integer value=N_;// function counter
    

//   void init(){}

//   template<typename...Inputs>
//   static auto eval(const Inputs&...inputs){FuncType::eval(inputs...);}


//   Function(){} 


// private:
//     FuncType func_;
// };

// template<typename FuncType, typename FunctionSpace>
// constexpr auto FunctionInit(){return Function<FuncType,FunctionSpace>();} 


// template<typename FuncType, typename FunctionSpace,Integer N, typename TupleSpaces>
// class QuadratureOrder<Function<FuncType,FunctionSpace,N,TupleSpaces>>
// { public:
//   using type=Function<FuncType,FunctionSpace,N,TupleSpaces>;
//   // using UniqueElementFunctionSpacesTupleType=typename Test::UniqueElementFunctionSpacesTupleType;
//   // using BaseFunctionSpaceAndElement=GetType<UniqueElementFunctionSpacesTupleType,N>;
//   using Operator=typename type::Operator;
//   using ElementFunctionSpace=typename type::ElementFunctionSpace;
//   using BaseFunctionSpace=Elem2FunctionSpace<ElementFunctionSpace>;
//   static constexpr Integer value=QuadratureOrder<Operator,BaseFunctionSpace>::value;

// };







template<typename FullSpace,Integer N,typename Operator_=IdentityOperator,typename FuncType=EmptyClass>
class Function2: public Expression<Function2<FullSpace,N,Operator_,FuncType>>
{
  
  public:

    using FunctionSpace=FullSpace;
    using FunctionType=FuncType;
    using Operator=Operator_;
    static constexpr Integer number=N;
    static constexpr Integer value=GetType<typename FunctionSpace::FromElementFunctionSpacesToUniqueNumbersTupleType,N>::value;
    using Elem=typename FunctionSpace::Elem;
    static constexpr Integer Dim=Elem::Dim;
    using UniqueElementFunctionSpacesTupleType=typename FunctionSpace::UniqueElementFunctionSpacesTupleType;
    // using ElementFunctionSpace=typename FunctionSpace::ElementFunctionSpace;
    static constexpr Integer Nmax=TupleTypeSize<UniqueElementFunctionSpacesTupleType>::value;
    static constexpr Integer Nelem_dofs=FunctionSpace::Nelem_dofs_array[N];
    using LocalDofs=Array<Real,Nelem_dofs>;
 
  // template<typename...Inputs>
  // static auto eval(const Inputs&...inputs){FuncType::eval(inputs...);}

  Function2(){}
  
  Function2(const FullSpace& AuxW):
  space_ptr_(std::make_shared<FullSpace>(AuxW))
  {} 


  const LocalDofs& local_dofs()const{return local_dofs_;}


private:
    std::shared_ptr<FullSpace> space_ptr_;
    LocalDofs local_dofs_;
};


template<Integer N,typename OperatorType,typename FullSpace,typename FuncType>
class QuadratureOrder<Function2<FullSpace,N,OperatorType,FuncType>>
{ public:
  using type=Function2<FullSpace,N,OperatorType,FuncType>;
  // using UniqueElementFunctionSpacesTupleType=typename Test::UniqueElementFunctionSpacesTupleType;
  // using BaseFunctionSpaceAndElement=GetType<UniqueElementFunctionSpacesTupleType,N>;
  using Operator=typename type::Operator;
  using UniqueElementFunctionSpacesTupleType=typename type::UniqueElementFunctionSpacesTupleType;
  using BaseFunctionSpaceAndElement=GetType<UniqueElementFunctionSpacesTupleType,type::value>;
  using BaseFunctionSpace=Elem2FunctionSpace<BaseFunctionSpaceAndElement>;
  static constexpr Integer value=QuadratureOrder<Operator,BaseFunctionSpace>::value;

};



template<Integer N,typename FuncType=EmptyClass,typename OperatorType=IdentityOperator,typename FullSpace>
auto MakeFunction(const FullSpace& AuxW){return Function2<FullSpace,N+FullSpace::TrialSpaceSize,OperatorType,FuncType>(AuxW);}


}




#endif
