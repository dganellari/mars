#ifndef MARS_EVALUATION_FUNCTION_HPP
#define MARS_EVALUATION_FUNCTION_HPP
#include "mars_base.hpp"
#include "mars_function.hpp"


namespace mars {



template<typename FuncType,Integer N,typename Operator_,typename FullSpace,typename...OtherTemplateArguments>
class Evaluation<Expression<Function<FullSpace,N,Operator_,FuncType>>,OtherTemplateArguments...>
{
 public:
 using type=Function<FullSpace,N,Operator_,FuncType>;
 using FunctionSpace=typename type::FunctionSpace;
 using FunctionSpaces=GetType<typename type::UniqueElementFunctionSpacesTupleType,type::value>;
 using Elem=typename FunctionSpaces::Elem;
 using BaseFunctionSpace=Elem2FunctionSpace<FunctionSpaces>;
 using Operator=typename type::Operator;
 using value_type=OperatorType<type,OtherTemplateArguments...>;
 Evaluation(){};

 Evaluation(const type& expr):
 eval_(expr)
 {};

 template<typename Shapes>
 void compute(value_type& value, type&eval,const Jacobian<Elem>&J, const Shapes& shapes )
 {
  using type1=typename Shapes::type; 
  using type2=typename Shapes::subtype;
  using type3=typename Shapes::subtype::subtype; 
   eval.local_dofs_update(J);
   const auto& local_dofs=eval.local_dofs();
    // loop on qp points
    for(Integer ii=0;ii< type2::Dim;ii++)
    {
      for(Integer mm=0;mm<type3::Rows;mm++)
        for(Integer nn=0;nn<type3::Cols;nn++)
           value[ii](mm,nn)=local_dofs[0]*shapes[0][ii](mm,nn);
    // loop on dofs
     for(Integer jj=1;jj< type1::Dim;jj++)
      // loop on components of subtype (is a matrix)
      for(Integer mm=0;mm<type3::Rows;mm++)
        for(Integer nn=0;nn<type3::Cols;nn++)
           value[ii](mm,nn)+=local_dofs[jj]*shapes[jj][ii](mm,nn);

    }


 }
 
 // template<typename...Forms, typename Jacobian,typename...Inputs>
 // constexpr void apply(value_type& value,const Jacobian& J, const typename ShapeFunctions2<Forms...>::TupleOfTupleShapeFunction& shape_functions,const Inputs&...inputs)
 template<typename...Args, typename Jacobian,typename...Inputs>
 constexpr void apply(value_type& value,const Jacobian& J, const std::tuple<Args...>& shape_functions,const Inputs&...inputs)

 {
  using TupleOfTupleShapeFunction=std::tuple<Args...>;
  using tuple_type=GetType<TupleOfTupleShapeFunction,type::value>;
  // const auto& tuple=tuple_get<type::value>(shape_functions());
  const auto& tuple=tuple_get<type::value>(shape_functions);

  // using tuple_type=GetType<std::tuple<Args...>,type::value>;
  // const auto& tuple=tuple_get<type::value>(tuple_of_tuple);
  constexpr Integer M=TypeToTupleElementPosition<ShapeFunction<Elem,BaseFunctionSpace,Operator,OtherTemplateArguments...>,tuple_type>::value;
  // FIXME


  const auto& shapes=tuple_get<M>(tuple).eval();
  compute(value,eval_,J,shapes);
  std::cout<<"Evaluation<Expression<Function "<<value<<std::endl;


 }


private:
 type eval_;
};





}
#endif