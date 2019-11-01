#ifndef MARS_EVALUATION_GENERAL_FORM_UTILS_HPP
#define MARS_EVALUATION_GENERAL_FORM_UTILS_HPP
#include "mars_base.hpp"
#include "mars_shape_functions_collection.hpp"

namespace mars {

    // BILINEAR FORMS
    // DOT OF DOTS BETWEEN FQPVALUES AND FQPVALUES (same type)
    // template<Integer Dim>
    // inline constexpr Real dotofdots(const Vector<Real, Dim> &left, const Vector<Real,Dim> &right)
    // {
    //      // std::cout<<"dotofdots const Vector<Real, Dim> &left, const Vector<Real,Dim> &right"<<std::endl;

    //     return dot(left,right);
    // }
    

    // here we get also Tranposed<Matrix<T,Rows,Cols>> because Transposed inherits from Matrix<T,Cols,Rows>
    template<Integer Rows,Integer Cols>
    inline constexpr Real dotofdots(const Matrix<Real, Rows,Cols> &left, const Matrix<Real,Rows,Cols> &right)
    {

      // std::cout<<"dotofdots mat mat"<<std::endl;
        Real result = left(0,0)*right(0,0);
        
        for(Integer col = 1; col < Cols; ++col) 
            result += left(0,col)*right(0,col);
        for(Integer col = 0; col < Cols; ++col)
           for(Integer row = 1; row < Rows; ++row)  
            result += left(row,col)*right(row,col);  


         // std::cout<<"dotofdots const Matrix<Real, Rows,Cols> &left, const Matrix<Real,Rows,Cols> &right"<<std::endl;
         // std::cout<<left<<" "<<right<<std::endl;
         // std::cout<<left(0,0)<<" "<<right(0,0)<<std::endl;
        return result;
    }

    template<Integer Rows,Integer Cols>
    inline constexpr Real dotofdots(const Transposed<Matrix<Real,Cols,Rows>> &left, const Matrix<Real,Rows,Cols> &right)
    {
            // std::cout<<"dotofdots trans mat mat"<<std::endl;
            // std::cout<<left.rows()<<" "<<left.cols()<<std::endl;
            // std::cout<<right.rows()<<" "<<right.cols()<<std::endl;

        Real result = left(0,0)*right(0,0);
         std::cout<<0<<" "<<0<<std::endl;
        for(Integer col = 1; col < Cols; ++col) 
            {
               // std::cout<<"here1 "<<std::endl;
               // std::cout<<col<<", "<<0<<std::endl;
              // result += left(col,0)*right(0,col);
              result += left(0,col)*right(0,col);
            }
        for(Integer col = 0; col < Cols; ++col)
           for(Integer row = 1; row < Rows; ++row)  
            {
                             // std::cout<<"here2 "<<std::endl;

               // std::cout<<col<<", "<<row<<std::endl;
              // result += left(col,row)*right(row,col);  
              result += left(row,col)*right(row,col);  
            }

             
         // std::cout<<"dotofdots const Matrix<Real, Rows,Cols> &left, const Matrix<Real,Rows,Cols> &right"<<std::endl;
         // std::cout<<left<<" "<<right<<std::endl;
         // std::cout<<left(0,0)<<" "<<right(0,0)<<std::endl;
        return result;
    }

    template<Integer Rows,Integer Cols>
    inline constexpr Real dotofdots(const Matrix<Real,Rows,Cols> &left, const Transposed<Matrix<Real,Cols,Rows>> &right)
    {
          // std::cout<<"dotofdots mat trans_mat"<<std::endl;

        Real result = left(0,0)*right(0,0);
        
        for(Integer col = 1; col < Cols; ++col) 
            // result += left(0,col)*right(col,0);
          result += left(0,col)*right(0,col);
        for(Integer col = 0; col < Cols; ++col)
           for(Integer row = 1; row < Rows; ++row)  
            // result += left(row,col)*right(col,row);  
            result += left(row,col)*right(row,col);  

         // std::cout<<"dotofdots const Matrix<Real, Rows,Cols> &left, const Matrix<Real,Rows,Cols> &right"<<std::endl;
         // std::cout<<left<<" "<<right<<std::endl;
         // std::cout<<left(0,0)<<" "<<right(0,0)<<std::endl;
        return result;
    }

    // here we get also Tranposed<Matrix<T,Rows,Cols>> because Transposed inherits from Matrix<T,Cols,Rows>
    template<Integer Rows,Integer Cols>
    inline constexpr Real dotofdots(const Transposed<Matrix<Real,Cols, Rows>> &left, const Transposed<Matrix<Real,Cols,Rows>> &right)
    {
              // std::cout<<"dotofdots trans_mat trans_mat"<<std::endl;

        Real result = left(0,0)*right(0,0);
        
        for(Integer col = 1; col < Cols; ++col) 
            result += left(0,col)*right(0,col);
        for(Integer col = 0; col < Cols; ++col)
           for(Integer row = 1; row < Rows; ++row)  
            result += left(row,col)*right(row,col);  


         // std::cout<<"dotofdots const Matrix<Real, Rows,Cols> &left, const Matrix<Real,Rows,Cols> &right"<<std::endl;
         // std::cout<<left<<" "<<right<<std::endl;
         // std::cout<<left(0,0)<<" "<<right(0,0)<<std::endl;
        return result;
    }
    // template<Integer Rows,Integer Cols>
    // inline constexpr Real dotofdots(const Transposed<Matrix<Real, Cols,Rows>> &left, const Matrix<Real,Rows,Cols> &right)
    // {
    //     Real result = left(0,0)*right(0,0);
        
    //     for(Integer col = 1; col < Cols; ++col) 
    //         result += left(col,0)*right(0,col);
    //     for(Integer col = 0; col < Cols; ++col)
    //        for(Integer row = 1; row < Rows; ++row)  
    //         result += left(col,row)*right(row,col);  

    //     return result;
    // }

    template<Integer Dim, typename S,typename T>
    inline constexpr Real dotofdots(const Vector<S, Dim> &left, const Vector<T, Dim> &right)
    {
        Real ret = dotofdots(left(0),right(0));
        for(Integer d = 1; d < Dim; ++d) {
            ret += dotofdots(left(d),right(d));
        }

        return ret;
    }

    template<Integer Dim, typename S, typename T>
    inline constexpr Real dotofdots(const QPValues<S, Dim> &left, const Vector<T, Dim> &right)
    {


        Real ret = dotofdots(left[0],right[0]);

        for(Integer d = 1; d < Dim; ++d) {
            ret += dotofdots(left[d],right[d]);
        }
        // using ok1=decltype(left[0]);
        // using ok2=decltype(right[0]);
        // QPValues<S, Dim> ok1(1,2);
        // S ok2(2,2);
        // Vector<T, Dim> ok3(3,3);
        // T ok4(4,3);

        // decltype(right[0]) ok2(333.44,5);
        std::cout<<dotofdots(left[0],right[0])<<std::endl;
        return ret;
    }

    template<typename S, typename T, Integer Rows,Integer Cols>
    inline constexpr Real dotofdots(const Matrix<S, Rows,Cols> &left, const Matrix<T,Rows,Cols> &right)
    {
        Real result = dotofdots(left(0,0),right(0,0));
        
        for(Integer col = 1; col < Cols; ++col) 
            result += dotofdots(left(0,col),right(0,col));
        for(Integer col = 0; col < Cols; ++col)
           for(Integer row = 1; row < Rows; ++row)  
            result += dotofdots(left(row,col),right(row,col));  

        return result;
    }

    // DOT OF DOTS BETWEEN QPVALUES AND FQPVALUES (T and Tensor<T>)
    template<typename S, Integer Dim, typename T>
    inline constexpr Real dotofdots(const S &left, const Vector<T, Dim> &right)
    {
      std::cout<<" eee"<<std::endl;
        Real ret = dotofdots(left,right(0));
        for(Integer d = 1; d < Dim; ++d) {
            ret += dotofdots(left,right(d));
        }

        return ret;
    }

    template<typename S, Integer Dim, typename T>
    inline constexpr Real dotofdots(const Vector<S, Dim> &left, const T &right)
    {
        Real ret = dotofdots(left(0),right);
        for(Integer d = 1; d < Dim; ++d) {
            ret += dotofdots(left(d),right);
        }

        return ret;
    }

    template<typename S, typename T, Integer Rows,Integer Cols>
    inline constexpr Real dotofdots(const S &left, const Matrix<T,Rows,Cols> &right)
    {
        Real result = left*right(0,0);
        
        for(Integer col = 1; col < Cols; ++col) 
            result += left*right(0,col);
        for(Integer col = 0; col < Cols; ++col)
           for(Integer row = 1; row < Rows; ++row)  
            result += left*right(row,col);  

        return result;
    }

    template<typename S, typename T, Integer Rows,Integer Cols>
    inline constexpr Real dotofdots(const Matrix<S, Rows,Cols> &left, const T &right)
    {
        Real result = left(0,0)*right;
        
        for(Integer col = 1; col < Cols; ++col) 
            result += left(0,col)*right;
        for(Integer col = 0; col < Cols; ++col)
           for(Integer row = 1; row < Rows; ++row)  
            result += left(row,col)*right;  

        return result;
    }
    
template<bool PositiveWeights,typename TestTrialSpaces,typename L2,typename Form>
class LocalTensor;



template<typename TestTrialSpaces, typename Left,typename Right,bool VolumeIntegral, Integer QR,typename Form>
class LocalTensor<false,TestTrialSpaces,L2DotProductIntegral<Left,Right,VolumeIntegral,QR>,Form>
// template<typename TestTrialSpaces, typename MeshT, typename Left,typename Right,Integer QR,typename Form>
// class LocalTensor<false,TestTrialSpaces,L2DotProductIntegral<MeshT,Left,Right,QR>,Form>
{
    static_assert("negative weights not permitted");
};


// template<typename TestTrialSpaces, typename MeshT, typename Left,typename Right,Integer QR,typename Form>
// class LocalTensor<true,TestTrialSpaces,L2DotProductIntegral<MeshT,Left,Right,QR>,Form>
template<typename TestTrialSpaces, typename Left,typename Right,bool VolumeIntegral,Integer QR>
class LocalTensor<true,TestTrialSpaces,L2DotProductIntegral<Left,Right,VolumeIntegral,QR>,Number<1>>
{
 public:
 using type= L2DotProductIntegral<Left,Right,VolumeIntegral,QR>;
 using QRule=typename type::QRule;
 using EvalLeft=Evaluation<Expression<Left>,QRule>;
 using EvalRight=Evaluation<Expression<Right>,QRule>;
 using subtype= OperatorType<type,Number<1>>;
 
 LocalTensor(){}

 LocalTensor(const type& expr)
 :
 eval_left_(expr.left()),
 eval_right_(expr.right())
 {
  // Left ooo(5,6,7,5);
 }


 template<typename Elem, typename T,typename S, Integer NQPoints, Integer NComponents>
 void apply_aux(subtype& vec, const FiniteElem<Elem>& J, 
               const FQPValues<T,NQPoints,NComponents> & fqp_values,
               const  QPValues<S,NQPoints> & qp_values)
 {
  const auto& detJ=J.get_det();
  // for(Integer ii=0;ii<vec.size();ii++)
    // {   
    //     vec[ii]=detJ*dotofdots(qp_values,fqp_values[ii]);
    // }
  for(Integer ii=0;ii<vec.rows();ii++)
    {   
        vec(ii,0)=detJ*dotofdots(qp_values,fqp_values[ii]);
    }
 }

 template<typename Elem, typename T,typename S, Integer NQPoints, Integer NComponents>
 void apply_aux(subtype& vec, const FiniteElem<Elem>& J, 
               const  QPValues<T,NQPoints> & qp_values,
               const FQPValues<S,NQPoints,NComponents> & fqp_values)
 {
  const auto& detJ=J.get_det();
  // for(Integer ii=0;ii<vec.size();ii++)
    // {   
    //     vec[ii]=detJ*dotofdots(qp_values,fqp_values[ii]);
    // }
  for(Integer ii=0;ii<vec.rows();ii++)
    {   
        vec(ii,0)=detJ*dotofdots(qp_values,fqp_values[ii]);
    }
 }

 template<typename Elem,typename...Forms>
  void apply(subtype& vec, const FiniteElem<Elem>& J, ShapeFunctionsCollection<Forms...>& shape_functions)

 {
  eval_left_.apply(left_value_,J,shape_functions.tuple<VolumeIntegral>(),shape_functions.template composite_tensor<VolumeIntegral>(),shape_functions.template composite_shapes<VolumeIntegral>());
  eval_right_.apply(right_value_,J,shape_functions.tuple<VolumeIntegral>(),shape_functions.template composite_tensor<VolumeIntegral>(),shape_functions.template composite_shapes<VolumeIntegral>());
  apply_aux(vec,J,left_value_,right_value_);
  std::cout<<"detJ=="<<J.get_det()<<std::endl;
  std::cout<<"vec=="<<vec<<std::endl;
 }

  
 private:
 EvalLeft eval_left_;
 EvalRight eval_right_;
 OperatorType<Left,QRule> left_value_;
 OperatorType<Right,QRule> right_value_;

};







// template<typename TestTrialSpaces, typename MeshT, typename Left,typename Right,Integer QR,typename Form>
// class LocalTensor<true,TestTrialSpaces,L2DotProductIntegral<MeshT,Left,Right,QR>,Form>
template<typename TestTrialSpaces, typename Left,typename Right,bool VolumeIntegral,Integer QR>
class LocalTensor<true,TestTrialSpaces,L2DotProductIntegral<Left,Right,VolumeIntegral,QR>,Number<2>>
{
 public:
 using type= L2DotProductIntegral<Left,Right,VolumeIntegral,QR>;
 using QRule=typename type::QRule;
 using EvalLeft=Evaluation<Expression<Left>,QRule>;
 using EvalRight=Evaluation<Expression<Right>,QRule>;
 using subtype= OperatorType<type,Number<2>>;
 
 LocalTensor(){}

 LocalTensor(const type& expr)
 :
 // eval_left_<QRule>(Eval(expr.left())),
 // eval_right_<QRule>(Eval(expr.right()))
 eval_left_(expr.left()),
 eval_right_(expr.right())
 {}

 template<typename Elem,typename ShapeFunctions>
 void apply(subtype& mat, const FiniteElem<Elem>& J, const ShapeFunctions& shape_functions)
 {

  const auto& detJ=J.get_det();
  // todo fixme
  // eval_left_.apply(left_value_,J,shape_functions);
  // eval_right_.apply(right_value_,J,shape_functions);

  eval_left_.apply(left_value_,J,shape_functions.tuple<VolumeIntegral>(),shape_functions.template composite_tensor<VolumeIntegral>(),shape_functions.template composite_shapes<VolumeIntegral>());
  eval_right_.apply(right_value_,J,shape_functions.tuple<VolumeIntegral>(),shape_functions.template composite_tensor<VolumeIntegral>(),shape_functions.template composite_shapes<VolumeIntegral>());
// shape_functions_.composite_tensor(),shape_functions_.composite_shapes());
  std::cout<<"left_value_"<<std::endl;
  // std::cout<<left_value_<<std::endl;
  std::cout<<"right_value_"<<std::endl;
  // std::cout<<right_value_<<std::endl;
  // std::cout<<"quiiiiiiiii 111"<<std::endl;
 // decltype(left_value_[0]) ee(6);
 // decltype(right_value_[0]) e4e(6);
  for(Integer ii=0;ii<mat.rows();ii++)
    for(Integer jj=0;jj<mat.cols();jj++)
    {   

        // std::cout<<"(ii,jj)=("<<ii<<","<<jj<<")"<<std::endl;
        // left value is trial, right value is test

        mat(ii,jj)=detJ*dotofdots(right_value_[ii],left_value_[jj]);
    }
    // std::cout<<detJ<<std::endl;
    std::cout<<"mat="<<mat<<std::endl;
   // std::cout<<"quiiiiiiiii 222"<<std::endl;
 }

  
 private:
 EvalLeft eval_left_;
 EvalRight eval_right_;
 OperatorType<Left,QRule> left_value_;
 OperatorType<Right,QRule> right_value_;

};

template<typename Tuple,Integer Nmax,Integer N, typename ShapeFunctions>
class EvalOfL2InnersAuxHelper;


template<typename Tuple,Integer Nmax, typename ShapeFunctions>
class EvalOfL2InnersAuxHelper<Tuple,Nmax,Nmax,ShapeFunctions>
{
 public:
    using single_type=std::conditional_t<IsSame<GetType<Tuple,Nmax>,EmptyExpression>::value,
                                         std::tuple<>,
                                         Evaluation<Expression<GetType<Tuple,Nmax>>,ShapeFunctions>>;
    using type=std::tuple<single_type>;
};

template<typename Tuple,Integer Nmax,Integer N, typename ShapeFunctions>
class EvalOfL2InnersAuxHelper
{
 public:
    using single_type=std::conditional_t<IsSame<GetType<Tuple,N>,EmptyExpression>::value,
                                         std::tuple<>,
                                         Evaluation<Expression<GetType<Tuple,N>>,ShapeFunctions>>;
    using type=TupleCatType<std::tuple<single_type>,typename EvalOfL2InnersAuxHelper<Tuple,Nmax,N+1,ShapeFunctions>::type >;
};



template<typename Tuple,Integer N, typename ShapeFunctions>
using EvalOfL2InnersAux=typename EvalOfL2InnersAuxHelper<Tuple,TupleTypeSize<Tuple>::value-1,N,ShapeFunctions>::type;

// template<typename Tuple, typename ShapeFunctions>
// using EvalOfL2InnersType=EvalOfL2InnersAux<Tuple,0,ShapeFunctions>;


// template<typename Tuple,Integer N,typename FullSpace,typename ShapeFunctions>
// constexpr std::enable_if_t<(N==TupleTypeSize<Tuple>::value-1), EvalOfL2InnersAux<Tuple,N,ShapeFunctions> > 
// EvalOfL2InnersHelper(const FullSpace&W,ShapeFunctions& shape_functions)
// {   
//       using ens0=GetType<Tuple,N>;
//       using ens=Evaluation<Expression<ens0>,ShapeFunctions>;
//       return std::tuple<ens>(Eval(Constructor<ens0>::apply(W),shape_functions));
// }


// template<typename Tuple,Integer N,typename FullSpace,typename ShapeFunctions>
// constexpr std::enable_if_t< (N<TupleTypeSize<Tuple>::value-1), EvalOfL2InnersAux<Tuple,N,ShapeFunctions>>  
// EvalOfL2InnersHelper(const FullSpace&W,ShapeFunctions& shape_functions)
// {
//       using ens0=GetType<Tuple,N>;
//       using ens=Evaluation<Expression<ens0>,ShapeFunctions>;
//       return std::tuple_cat(std::tuple<ens>(Eval(Constructor<ens0>::apply(W),shape_functions)),
//                             EvalOfL2InnersHelper<Tuple,N+1>(W,shape_functions));
// }


// template<typename Tuple,typename FullSpace,typename ShapeFunctions>
// constexpr auto EvalOfL2Inners(const FullSpace&W,ShapeFunctions& shape_functions)
// {
//       return EvalOfL2InnersHelper<Tuple,0>(W,shape_functions);
// }




template<typename GeneralForm, Integer N >
class LocalTupleOfTensors;

template<typename GeneralForm >
class LocalTupleOfTensors<GeneralForm,1>
{
public:
    using TupleOfNumbers=typename GeneralForm::TupleOfPairsNumbers;
    static constexpr auto Nelem_dofs_array=GeneralForm::FunctionSpace::Nelem_dofs_array;

    template<Integer Nmax, Integer N>
    class AuxHelper;

    template<Integer Nmax>
    class AuxHelper<Nmax,Nmax>
    {
        public:
        using Numbers=GetType<TupleOfNumbers,Nmax>;
        static constexpr Integer Dim=Nelem_dofs_array[GetType<Numbers>::value];
        // using type= std::tuple< Vector<Real,Dim> >;
        using type= std::tuple< Matrix<Real,Dim,1> >;
    };

    template<Integer Nmax, Integer N>
    class AuxHelper
    {
        public:
        using Numbers=GetType<TupleOfNumbers,N>;
        static constexpr Integer Dim=Nelem_dofs_array[GetType<Numbers>::value];
        // using type=TupleCatType< std::tuple< Vector<Real,Dim>>, typename AuxHelper<Nmax,N+1>::type >;
        using type=TupleCatType< std::tuple< Matrix<Real,Dim,1>>, typename AuxHelper<Nmax,N+1>::type >;
    };

    using type=typename AuxHelper<TupleTypeSize<TupleOfNumbers>::value-1,0>::type;
    
};



template<typename GeneralForm >
class LocalTupleOfTensors<GeneralForm,2>
{
public:
    using TupleOfPairsNumbers=typename GeneralForm::TupleOfPairsNumbers;
    static constexpr auto Nelem_dofs_array=GeneralForm::FunctionSpace::Nelem_dofs_array;

    template<Integer Nmax, Integer N>
    class AuxHelper;

    template<Integer Nmax>
    class AuxHelper<Nmax,Nmax>
    {
        public:
        using Numbers=GetType<TupleOfPairsNumbers,Nmax>;
        static constexpr Integer Rows=Nelem_dofs_array[GetType<Numbers,0>::value];
        static constexpr Integer Cols=Nelem_dofs_array[GetType<Numbers,1>::value];
        using type= std::tuple< Matrix<Real,Rows,Cols>>;
    };

    template<Integer Nmax, Integer N>
    class AuxHelper
    {
        public:
        using Numbers=GetType<TupleOfPairsNumbers,N>;
        static constexpr Integer Rows=Nelem_dofs_array[GetType<Numbers,0>::value];
        static constexpr Integer Cols=Nelem_dofs_array[GetType<Numbers,1>::value];
        using type=TupleCatType< std::tuple< Matrix<Real,Rows,Cols>>, typename AuxHelper<Nmax,N+1>::type >;
    };

    using type=typename AuxHelper<TupleTypeSize<TupleOfPairsNumbers>::value-1,0>::type;
    
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//// build_tuple_of_evals<TupleOfPairsNumbers>(Expr expr , ShapeFunctionsCollection<...Forms>shapes) /////
//// -TupleOfPairsNumbers={(0,0),(0,1),(0,2),(1,0),(2,0)}                                            /////
////  where (i,j) represents the coupling of the i test with the j trial                             /////
//// -Expr is a form                                                                                 /////
//// We want to build a tuple which associate to each element of TupleOfPairsNumbers                 /////
//// the corresponding Evaluation of the form                                                        /////
////                                                                                                 /////                                                            /////
//// Example:                                                                                        /////    
//// X=[P0,P1,P2] -> u,v in P0; r,s in P1; p,q in P2                                                 /////
//// Expr= int u v + int gradu gradv + int u s + int u q  - int v r - int v p                        /////
//// TupleOfPairsNumbers={(0,0),(0,1),(0,2),(1,0),(2,0)}                                             /////
//// Then the output of the function is:                                                             /////  
////                                                                                                 /////  
//// { Eval(int u v + int gradu gradv),                                                              /////                                                               /////  
////   Eval( int u s),                                                                               /////                                                               /////  
////   Eval(int u q ),                                                                               /////  
////   Eval( int (-v r)),                                                                            /////  
////   Eval(int (-v p))  }                                                                           /////  
////                                                                                                 /////  
//// To create this, we move along TupleOfPairsNumbers. Se we fix for example T=(0,0).               /////  
//// Then we move along the whole expressions, searching for each term related to T.                 /////  
//// If it is there, we add its Evaluation. So for example Eval(int u s)                             /////                      /////  
//// However, if other terms are there, we do not simply sum them up, like:                          /////   
//// Eval(int u s)+Eval(int gradu gradv)                                                             /////  
//// Instead we create Eval(int u s+int gradu gradv)                                                 ///// 
//// Since, by inspectioning the form, we do not know how many other term related to T are there     /////
//// we must destroy and reconstruct the evaluation evrytime we find a new one                       ///// 
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//// todo fixme
//// we can specialise this to the case where we explicitly ask for a volume or surface integral





template<typename Left,typename Right,bool VolumeIntegral,Integer QR, typename...Forms>
constexpr auto build_tuple_of_evals_aux_aux(const L2DotProductIntegral<Left,Right,VolumeIntegral,QR>& l2prod,
                                            ShapeFunctionsCollection<Forms...>&shape_functions)
{
  using L2dot=L2DotProductIntegral<Left,Right,VolumeIntegral,QR>;

  return Evaluation<Expression<L2dot>,ShapeFunctionsCollection<Forms...>>
        (Eval(l2prod,shape_functions));
}


template<typename Left,typename Right,bool VolumeIntegral,Integer QR, typename...Forms>
constexpr auto build_tuple_of_evals_aux_aux(const std::tuple<>& null, 
                           const L2DotProductIntegral<Left,Right,VolumeIntegral,QR>& l2prod,
                                 ShapeFunctionsCollection<Forms...>&shape_functions)
{
  return build_tuple_of_evals_aux_aux(l2prod,shape_functions);
  // using L2dot=L2DotProductIntegral<Left,Right,VolumeIntegral,QR>;

  // return Evaluation<Expression<L2dot>,ShapeFunctionsCollection<Forms...>>
  //       (Eval(l2prod,shape_functions));
}

template<typename Left,typename Right,bool VolumeIntegral,Integer QR, typename...Forms>
constexpr auto build_tuple_of_evals_aux_aux(const L2DotProductIntegral<Left,Right,VolumeIntegral,QR>& l2prod,
                                            const std::tuple<>& null,
                                                  ShapeFunctionsCollection<Forms...>&shape_functions)
{
    return build_tuple_of_evals_aux_aux(l2prod,shape_functions);

  // using L2dot=L2DotProductIntegral<Left,Right,VolumeIntegral,QR>;

  // return Evaluation<Expression<L2dot>,ShapeFunctionsCollection<Forms...>>
  // (Eval(l2prod,shape_functions));

}

template<typename Left1,typename Right1,bool VolumeIntegral1, Integer QR1,
         typename Left2,typename Right2,bool VolumeIntegral2, Integer QR2, typename...Forms>
constexpr auto build_tuple_of_evals_aux_aux(const L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>& l2prod1, 
                           const L2DotProductIntegral<Left2,Right2,VolumeIntegral2,QR2>& l2prod2,
                                 ShapeFunctionsCollection<Forms...>&shape_functions)
{
   return Eval(l2prod1+l2prod2,shape_functions);
}

template<typename Left, typename Left1,typename Right1,bool VolumeIntegral1,Integer QR1, typename...Forms>
constexpr auto build_tuple_of_evals_aux_aux(const Left& left, 
                           const L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>& l2prod,
                                 ShapeFunctionsCollection<Forms...>&shape_functions)
{
  using L2dot1=L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>;
  const auto& left_=left.expression();
  return Eval(left_+l2prod,shape_functions);
}








template<typename TupleOfPairsNumbers,typename Tuple, typename Left,typename Right,bool VolumeIntegral,Integer QR, typename...Forms>
constexpr auto build_tuple_of_evals_aux(const Tuple& tuple,
                          const L2DotProductIntegral<Left,Right,VolumeIntegral,QR>& l2prod,
                                ShapeFunctionsCollection<Forms...>&shape_functions)
{
  using L2=L2DotProductIntegral<Left,Right,VolumeIntegral,QR>;
  using TestTrialNumbers=typename L2::TestTrialNumbers;
  auto tuple_nth=tuple_get< TypeToTupleElementPosition<TestTrialNumbers,TupleOfPairsNumbers>::value >(tuple);
  auto new_elem=build_tuple_of_evals_aux_aux(tuple_nth,l2prod,shape_functions);
  // todo fixme tuple add std::tuple_cat
  return tuple_change_element<TypeToTupleElementPosition<TestTrialNumbers,TupleOfPairsNumbers>::value>
   (tuple, decltype(new_elem)(new_elem));
}


template<typename TupleOfPairsNumbers,typename Tuple, typename Left,typename Right, typename...Forms>
constexpr auto build_tuple_of_evals_aux(const Tuple& tuple,
                          const Addition<Expression<Left>,Expression<Right>>& addition,
                                ShapeFunctionsCollection<Forms...>&shape_functions)
{
  auto tuple_new=build_tuple_of_evals_aux<TupleOfPairsNumbers>(tuple,addition.left(),shape_functions);
  return build_tuple_of_evals_aux<TupleOfPairsNumbers>(tuple_new,addition.right(),shape_functions);
}


template<typename TupleOfPairsNumbers, typename Expr, typename...Forms>
constexpr auto build_tuple_of_evals(const Expr& expr,ShapeFunctionsCollection<Forms...>&shape_functions)
{
  using emptytuple=TupleOfType<TupleTypeSize<TupleOfPairsNumbers>::value,std::tuple<> > ;
  return build_tuple_of_evals_aux<TupleOfPairsNumbers,emptytuple>(emptytuple(),expr,shape_functions);
}





















// template<Integer H,typename TupleOfPairsNumbers,typename Tuple, typename Left,typename Right,bool VolumeIntegral,Integer QR, typename...Forms>
// constexpr std::enable_if_t< ! ((H==-1)||
//                                (H==0&&VolumeIntegral==true)||
//                                (H==1&&VolumeIntegral==false) ),
// Tuple> 
// build_tuple_of_evals_aux2(const Tuple& tuple,
//                          const L2DotProductIntegral<Left,Right,VolumeIntegral,QR>& l2prod,
//                                ShapeFunctionsCollection<Forms...>&shape_functions)
// {return tuple;}







// template<Integer H,typename TupleOfPairsNumbers,typename Tuple, typename Left,typename Right,bool VolumeIntegral,Integer QR, typename...Forms>
// constexpr std::enable_if_t<(H==-1)||
//                            (H==0&&VolumeIntegral==true)||
//                            (H==1&&VolumeIntegral==false),
// Evaluation<Expression<Addition<
// // we take the n-th type (we remove Evaluation) of the tuple
// Expression<typename GetType<Tuple,TypeToTupleElementPosition<typename L2DotProductIntegral<Left,Right,VolumeIntegral,QR>::TestTrialNumbers,TupleOfPairsNumbers>::value>::type>,
// // and we add it to integral
// Expression<L2DotProductIntegral<Left,Right,VolumeIntegral,QR>>
//                               >>>
// > 
// build_tuple_of_evals_aux2(const Tuple& tuple,
//                          const L2DotProductIntegral<Left,Right,VolumeIntegral,QR>& l2prod,
//                                ShapeFunctionsCollection<Forms...>&shape_functions)
// {
//   using L2=L2DotProductIntegral<Left,Right,VolumeIntegral,QR>;
//   using TestTrialNumbers=typename L2::TestTrialNumbers;
//   auto tuple_nth=tuple_get< TypeToTupleElementPosition<TestTrialNumbers,TupleOfPairsNumbers>::value >(tuple);
//   auto new_elem=build_tuple_of_evals_aux_aux(tuple_nth,l2prod,shape_functions);
//   return tuple_change_element<TypeToTupleElementPosition<TestTrialNumbers,TupleOfPairsNumbers>::value>
//    (tuple, decltype(new_elem)(new_elem));
// }


template<Integer H,typename TupleOfPairsNumbers,typename Tuple, typename Left,typename Right,bool VolumeIntegral,Integer QR, typename...Forms>
constexpr std::enable_if_t<
                           !((H==-1)||
                            (H==0&&VolumeIntegral==true)||
                            (H==1&&VolumeIntegral==false))
                           ,
                           Tuple       
                          > 

build_tuple_of_evals_aux2(const Tuple& tuple,
                          const L2DotProductIntegral<Left,Right,VolumeIntegral,QR>& l2prod,
                                ShapeFunctionsCollection<Forms...>&shape_functions)
{
  return tuple;
}




template<Integer H,typename TupleOfPairsNumbers,typename Tuple, typename Left,typename Right,bool VolumeIntegral,Integer QR, typename...Forms>
constexpr std::enable_if_t<
                           ((H==-1)||
                            (H==0&&VolumeIntegral==true)||
                            (H==1&&VolumeIntegral==false))
                           && 
                           IsSame<GetType<Tuple,TypeToTupleElementPosition<typename L2DotProductIntegral<Left,Right,VolumeIntegral,QR>::TestTrialNumbers,TupleOfPairsNumbers>::value>,std::tuple<>>::value,


                           // (H==-1)||
                           // (H==0&&VolumeIntegral==true)||
                           // (H==1&&VolumeIntegral==false),
// we change the n-th component of tuple 
TupleChangeType<
                // given the integral, we search the position of TestTrialNumbers
                TypeToTupleElementPosition<typename L2DotProductIntegral<Left,Right,VolumeIntegral,QR>::TestTrialNumbers,TupleOfPairsNumbers>::value,
                Evaluation<                                                // and we add it to integral
                                                Expression<L2DotProductIntegral<Left,Right,VolumeIntegral,QR>>,
                           ShapeFunctionsCollection<Forms...>                         
                          >,
                Tuple   
                >
> 

build_tuple_of_evals_aux2(const Tuple& tuple,
                          const L2DotProductIntegral<Left,Right,VolumeIntegral,QR>& l2prod,
                                ShapeFunctionsCollection<Forms...>&shape_functions)
{
  using L2=L2DotProductIntegral<Left,Right,VolumeIntegral,QR>;
  using TestTrialNumbers=typename L2::TestTrialNumbers;
  auto tuple_nth=tuple_get< TypeToTupleElementPosition<TestTrialNumbers,TupleOfPairsNumbers>::value >(tuple);
  auto new_elem=build_tuple_of_evals_aux_aux(tuple_nth,l2prod,shape_functions);
  // todo fixme tuple add std::tuple_cat
  return tuple_change_element<TypeToTupleElementPosition<TestTrialNumbers,TupleOfPairsNumbers>::value>
   (tuple, decltype(new_elem)(new_elem));
}








template<Integer H,typename TupleOfPairsNumbers,typename Tuple, typename Left,typename Right,bool VolumeIntegral,Integer QR, typename...Forms>
constexpr std::enable_if_t<
                           ((H==-1)||
                            (H==0&&VolumeIntegral==true)||
                            (H==1&&VolumeIntegral==false))
                           && 
                           IsDifferent<GetType<Tuple,TypeToTupleElementPosition<typename L2DotProductIntegral<Left,Right,VolumeIntegral,QR>::TestTrialNumbers,TupleOfPairsNumbers>::value>,std::tuple<>>::value,



                           // (H==-1)||
                           // (H==0&&VolumeIntegral==true)||
                           // (H==1&&VolumeIntegral==false),
// we change the n-th component of tuple 
TupleChangeType<
                // given the integral, we search the position of TestTrialNumbers
                TypeToTupleElementPosition<typename L2DotProductIntegral<Left,Right,VolumeIntegral,QR>::TestTrialNumbers,TupleOfPairsNumbers>::value,
                Evaluation<Expression<Addition<
                                                // we take the n-th type (we remove Evaluation) of the tuple
                                                
                                                Expression<typename GetType<Tuple,TypeToTupleElementPosition<typename L2DotProductIntegral<Left,Right,VolumeIntegral,QR>::TestTrialNumbers,TupleOfPairsNumbers>::value>::type>,
                                                // and we add it to integral
                                                Expression<L2DotProductIntegral<Left,Right,VolumeIntegral,QR>>
                                               >>,
                           ShapeFunctionsCollection<Forms...>
                          >,
                Tuple   
                >
> 

build_tuple_of_evals_aux2(const Tuple& tuple,
                          const L2DotProductIntegral<Left,Right,VolumeIntegral,QR>& l2prod,
                                ShapeFunctionsCollection<Forms...>&shape_functions)
{
  using L2=L2DotProductIntegral<Left,Right,VolumeIntegral,QR>;
  using TestTrialNumbers=typename L2::TestTrialNumbers;
  auto tuple_nth=tuple_get< TypeToTupleElementPosition<TestTrialNumbers,TupleOfPairsNumbers>::value >(tuple);
  auto new_elem=build_tuple_of_evals_aux_aux(tuple_nth,l2prod,shape_functions);
  // todo fixme tuple add std::tuple_cat
  return tuple_change_element<TypeToTupleElementPosition<TestTrialNumbers,TupleOfPairsNumbers>::value>
   (tuple, decltype(new_elem)(new_elem));
}





// template<Integer H,typename TupleOfPairsNumbers,typename Tuple, typename Left,typename Right,bool VolumeIntegral,Integer QR, typename...Forms>
// constexpr auto build_tuple_of_evals_aux2(const Tuple& tuple,
//                           const L2DotProductIntegral<Left,Right,VolumeIntegral,QR>& l2prod,
//                                 ShapeFunctionsCollection<Forms...>&shape_functions)
// {
//   using L2=L2DotProductIntegral<Left,Right,VolumeIntegral,QR>;
//   using TestTrialNumbers=typename L2::TestTrialNumbers;
//   auto tuple_nth=tuple_get< TypeToTupleElementPosition<TestTrialNumbers,TupleOfPairsNumbers>::value >(tuple);
//   auto new_elem=build_tuple_of_evals_aux_aux(tuple_nth,l2prod,shape_functions);
//   // todo fixme tuple add std::tuple_cat
//   return tuple_change_element<TypeToTupleElementPosition<TestTrialNumbers,TupleOfPairsNumbers>::value>
//    (tuple, decltype(new_elem)(new_elem));
// }



template<Integer H,typename TupleOfPairsNumbers,typename Tuple, typename Left,typename Right, typename...Forms>
constexpr auto build_tuple_of_evals_aux2(const Tuple& tuple,
                                        const Addition<Expression<Left>,Expression<Right>>& addition,
                                              ShapeFunctionsCollection<Forms...>&shape_functions)
{
  auto tuple_new=build_tuple_of_evals_aux2<H,TupleOfPairsNumbers>(tuple,addition.left(),shape_functions);
  return build_tuple_of_evals_aux2<H,TupleOfPairsNumbers>(tuple_new,addition.right(),shape_functions);
}


template<Integer H,typename TupleOfPairsNumbers, typename Expr, typename...Forms>
constexpr auto build_tuple_of_evals2(const Expr& expr,ShapeFunctionsCollection<Forms...>&shape_functions)
{
  using emptytuple=TupleOfType<TupleTypeSize<TupleOfPairsNumbers>::value,std::tuple<> > ;
  return build_tuple_of_evals_aux2<H,TupleOfPairsNumbers>(emptytuple(),expr,shape_functions);
}

























// template<typename...Ts>
// class EvaluationOfL2Inners;

// template<typename Form,typename ShapeFunctions_>
// class EvaluationOfL2Inners<Evaluation<Expression<GeneralForm<Form>>,ShapeFunctions_>>
// {
// public: 
//     using EvaluationGeneralForm=Evaluation<Expression<GeneralForm<Form>>,ShapeFunctions_>;
//     using GeneralForm=typename EvaluationGeneralForm::type;
//     using ShapeFunctions=typename EvaluationGeneralForm::ShapeFunctions;
//     // using TupleOfTupleShapeFunction=typename ShapeFunctions::TupleOfTupleShapeFunctionV;
//     using LocalTupleOfTensors=typename LocalTupleOfTensors<GeneralForm,GetType<typename GeneralForm::form>::value>::type;    
//     using L2Products=typename EvaluationGeneralForm::L2Products;
//     using EvalOfL2InnersType=EvalOfL2InnersAux<L2Products,0,ShapeFunctions>;
//     using TupleOfPairsNumbers=typename GeneralForm::TupleOfPairsNumbers;
    

//     EvaluationOfL2Inners(const GeneralForm& form,ShapeFunctions& shapefunctions)
//     :
//     // eval_inners_(build_tuple_of_evals<TupleOfPairsNumbers>(form(),shapefunctions))
//     eval_inners_(build_tuple_of_evals2<0,TupleOfPairsNumbers>(form(),shapefunctions))
//     {}

//  template<Integer N,typename Output,typename FiniteElem>
//     constexpr void apply_aux_aux(Output& mat,FiniteElem&J)
//     {
//     std::cout<<"pre jacobian evalinners2"<<std::endl;
//      auto& local_mat=std::get<N>(tensor_tuple_);
//      local_mat.zero();
//      auto & eval_N=std::get<N>(eval_inners_);
//      eval_N.apply(local_mat,J);

//   std::cout<<"after jacobian evalinners2"<<std::endl;

//     }   

//  template<Integer Nmax,Integer N,typename Output,typename FiniteElem>
//     constexpr std::enable_if_t<(N==Nmax),void> apply_aux(Output& mat,FiniteElem&J)
//     {
//      apply_aux_aux<N>(mat,J);
//     }

//  template<Integer Nmax,Integer N,typename Output,typename FiniteElem>
//     constexpr std::enable_if_t<(N<Nmax),void> apply_aux(Output& mat,FiniteElem&J)
//     {
//      {apply_aux_aux<N>(mat,J);}
//       apply_aux<Nmax,N+1>(mat,J);
//      }

//  template<typename Output,typename FiniteElem>
//     void apply(Output& mat,FiniteElem&J)
//     {
//       apply_aux<TupleTypeSize<L2Products>::value-1,0>(mat,J);
//     }

// private:
//     LocalTupleOfTensors tensor_tuple_;
//     EvalOfL2InnersType eval_inners_;    
// };


template<Integer H,typename...Ts>
class EvaluationOfL2Inners;






template<Integer H, typename Form,typename ShapeFunctions_>
class EvaluationOfL2Inners<H, Evaluation<Expression<GeneralForm<Form>>,ShapeFunctions_>>
{
public: 
    using EvaluationGeneralForm=Evaluation<Expression<GeneralForm<Form>>,ShapeFunctions_>;
    using GeneralForm=typename EvaluationGeneralForm::type;
    using ShapeFunctions=typename EvaluationGeneralForm::ShapeFunctions;
    // using TupleOfTupleShapeFunction=typename ShapeFunctions::TupleOfTupleShapeFunctionV;
    using LocalTupleOfTensors=typename LocalTupleOfTensors<GeneralForm,GetType<typename GeneralForm::form>::value>::type;    
    using L2Products=typename EvaluationGeneralForm::template L2Products<H>;
    using EvalOfL2InnersType=EvalOfL2InnersAux<L2Products,0,ShapeFunctions>;
    using TupleOfPairsNumbers=typename GeneralForm::TupleOfPairsNumbers;
    

    EvaluationOfL2Inners(const GeneralForm& form,ShapeFunctions& shapefunctions)
    :
    eval_inners_(build_tuple_of_evals2<H,TupleOfPairsNumbers>(form(),shapefunctions))
    {}


 template<typename Eval,typename Output,typename FiniteElem>
    constexpr std::enable_if_t<IsDifferent<Eval,std::tuple<>>::value,void>
    apply_aux_aux_aux(Eval& eval,Output& local_mat,const FiniteElem&J)
    {
     eval.apply(local_mat,J);
    }  

 template<typename Eval,typename Output,typename FiniteElem>
    constexpr std::enable_if_t<IsSame<Eval,std::tuple<>>::value,void>
    apply_aux_aux_aux(Eval& eval,Output& local_mat,const FiniteElem&J)
    {}  

 template<Integer N,typename Output,typename FiniteElem>
    constexpr void apply_aux_aux(Output& mat,FiniteElem&J)
    {
    std::cout<<"pre jacobian evalinners2"<<std::endl;
     auto& local_mat=std::get<N>(tensor_tuple_);
     local_mat.zero();
     auto & eval_N=std::get<N>(eval_inners_);
     apply_aux_aux_aux(eval_N,local_mat,J);
     // eval_N.apply(local_mat,J);

  std::cout<<"after jacobian evalinners2"<<std::endl;

    }   

 template<Integer Nmax,Integer N,typename Output,typename FiniteElem>
    constexpr std::enable_if_t<(N==Nmax),void> apply_aux(Output& mat,FiniteElem&J)
    {
     apply_aux_aux<N>(mat,J);
    }

 template<Integer Nmax,Integer N,typename Output,typename FiniteElem>
    constexpr std::enable_if_t<(N<Nmax),void> apply_aux(Output& mat,FiniteElem&J)
    {
     {apply_aux_aux<N>(mat,J);}
      apply_aux<Nmax,N+1>(mat,J);
     }

 template<typename Output,typename FiniteElem>
    void apply(Output& mat,FiniteElem&J)
    {
      apply_aux<TupleTypeSize<L2Products>::value-1,0>(mat,J);
    }


private:
    LocalTupleOfTensors tensor_tuple_;
    EvalOfL2InnersType eval_inners_;    
};





}
#endif