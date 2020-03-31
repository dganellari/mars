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
 using Space=typename type::Space; 
 using Elem=typename FunctionSpaces::Elem;
 using BaseFunctionSpace=Elem2FunctionSpace<FunctionSpaces>;
 using Operator=typename type::Operator;
 using value_type=OperatorType<type,OtherTemplateArguments...>;
 static constexpr auto trace=TraceDofs<Space>::dofs();//trace_dofs<FunctionSpaces>(); //TraceDofs<FunctionSpace>::dofs();//
 // using trace_type=decltype(trace);//typename decltype(trace)::value_type; 
 using trace_type_tmp=typename decltype(trace)::value_type;
 using trace_type=ArrayChangeType<Real,trace_type_tmp>; 

 // Evaluation(){};

 Evaluation(const type& expr):
 eval_ptr_(std::make_shared<type>(expr)) 
 {
  // std::cout<<"-------constructor global_dofs_ptr="<<std::endl;
 };

 template<typename OperatorType, typename Shapes>
 std::enable_if_t<IsSame<OperatorType,TraceOperator>::value,void> compute
 (value_type& value, type&eval, FiniteElem<Elem>&FE, const Shapes& shapes )
 {
  using type1=typename Shapes::type; 
  using type2=typename Shapes::subtype;
  using type3=typename Shapes::subtype::subtype; 
  const Integer face=FE.side_id();

     // std::cout<<"Evaluation<Expression<Function TraceOperator="<<std::endl;

   eval.local_dofs_update(FE);
   const auto& all_local_dofs=eval.local_dofs();
   // eval_ptr_->local_dofs_update(FE);
   // const auto& all_local_dofs=eval_ptr_->local_dofs();
   // auto local_dofs=subarray(all_local_dofs,trace[face]);
   // decltype(trace[face]) m1(1); //const mars::Array<long, 4>
   // decltype(all_local_dofs) m2(1);//'const mars::Array<double, 6>
   // decltype(local_dofs_) m3(1); //const mars::Array<mars::Array<long, 4>, 3>
   

   subarray(local_dofs_,all_local_dofs,trace[face]);
   // std::cout<<"local_dofs TraceOperator="<<std::endl;
   // for(Integer i=0;i<all_local_dofs.size();i++)
   //  std::cout<<all_local_dofs[i]<<std::endl;
   // std::cout<<"type2::Dim"<<type2::Dim<<std::endl;
   // std::cout<<"type3::Rows"<<type3::Rows<<std::endl;
   // std::cout<<"type3::Cols"<<type3::Cols<<std::endl;
   // std::cout<<"shapes"<<std::endl;
   // std::cout<<shapes<<std::endl;
   // std::cout<<"local_dofs_"<<std::endl;
   // std::cout<<local_dofs_<<std::endl;
   // std::cout<<"___"<<std::endl;
    // loop on qp points
    for(Integer ii=0;ii< type2::Dim;ii++)
    {
      for(Integer mm=0;mm<type3::Rows;mm++)
        for(Integer nn=0;nn<type3::Cols;nn++)
           {value[ii](mm,nn)=local_dofs_[0]*shapes[0][ii](mm,nn);

            // std::cout<<value[ii](mm,nn)<<std::endl;
            // std::cout<<local_dofs_[0]<<" "<<shapes[0][ii](mm,nn)<<std::endl;
           }
           // std::cout<<"here ::Cols"<<type3::Cols<<std::endl;
    // loop on dofs
     for(Integer jj=1;jj< type1::Dim;jj++)
      // loop on components of subtype (is a matrix)
      for(Integer mm=0;mm<type3::Rows;mm++)
        for(Integer nn=0;nn<type3::Cols;nn++)
           {
            value[ii](mm,nn)+=local_dofs_[jj]*shapes[jj][ii](mm,nn);
            // std::cout<<"<IsSame<OperatorType,TraceOperator> value[ii](mm,nn)="<<value[ii](mm,nn)<<std::endl;
            // std::cout<<"<IsSame<OperatorType,TraceOperator> local_dofs_[jj]="<<local_dofs_[jj]<<std::endl;
          }

    }
    // std::cout<<"___"<<std::endl;


 }


 template<typename OperatorType, typename Shapes>
 std::enable_if_t<IsDifferent<OperatorType,TraceOperator>::value,void> compute
 (value_type& value, type&eval, FiniteElem<Elem>&FE, const Shapes& shapes )
 {
  using type1=typename Shapes::type; 
  using type2=typename Shapes::subtype;
  using type3=typename Shapes::subtype::subtype; 
  // std::cout<< " shapes"<<shapes<<std::endl;
  // std::cout<< " value"<<value<<std::endl;
  // std::cout<< " eval function compute 1"<<std::endl;
   eval.local_dofs_update(FE);
  // std::cout<< " eval function compute 2"<<std::endl;
   const auto& local_dofs=eval.local_dofs();
  // std::cout<< " eval function compute 3"<<std::endl;
  // std::cout<< " shapes"<<shapes<<std::endl;
  // std::cout<< " local_dofs"<<local_dofs<<std::endl;
  // std::cout<< " value"<<value<<std::endl;

   // std::cout<<"local_dofs not TraceOperator="<<local_dofs<<std::endl;
  //  std::cout<<shapes<<std::endl;
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
           {
            value[ii](mm,nn)+=local_dofs[jj]*shapes[jj][ii](mm,nn);
            // std::cout<<"<IsDifferent<OperatorType,TraceOperator> value[ii](mm,nn)="<<value[ii](mm,nn)<<std::endl;
            // std::cout<<"<IsDifferent<OperatorType,TraceOperator> local_dofs_[jj]="<<local_dofs[jj]<<std::endl;

          }

    }
  // std::cout<< " value post"<<value<<std::endl;

 }
 
 // template<typename...Forms, typename FiniteElem,typename...Inputs>
 // constexpr void apply(value_type& value,const FiniteElem& J, const typename ShapeFunctions2<Forms...>::TupleOfTupleShapeFunction& shape_functions,const Inputs&...inputs)
 template<typename...Args, typename FiniteElem,typename...Inputs>
 constexpr void apply(value_type& value, FiniteElem& FE, const std::tuple<Args...>& shape_functions,const Inputs&...inputs)

 {
  using TupleOfTupleShapeFunction=std::tuple<Args...>;
  using tuple_type=GetType<TupleOfTupleShapeFunction,type::value>;
  // const auto& tuple=tuple_get<type::value>(shape_functions());
  // std::cout<< " eval function 1"<<std::endl;
  const auto& tuple=tuple_get<type::value>(shape_functions);
  // std::cout<< " eval function 2"<<std::endl;
  // using tuple_type=GetType<std::tuple<Args...>,type::value>;
  // const auto& tuple=tuple_get<type::value>(tuple_of_tuple);
  constexpr Integer M=TypeToTupleElementPosition<ShapeFunction<Elem,BaseFunctionSpace,Operator,OtherTemplateArguments...>,tuple_type>::value;
  // FIXME

  // std::cout<< " eval function 3"<<std::endl;
  const auto& shapes=tuple_get<M>(tuple).eval();
  // std::cout<< " eval function 4"<<std::endl;
  compute<Operator>(value,*eval_ptr_,FE,shapes);
  // std::cout<< " eval function 5"<<std::endl;
  // std::cout<<"Evaluation<Expression<Function "<<value<<std::endl;


 }


private:
 std::shared_ptr<type> eval_ptr_;
 trace_type local_dofs_;
};





}
#endif