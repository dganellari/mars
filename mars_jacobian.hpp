#ifndef MARS_JACOBIAN_HPP
#define MARS_JACOBIAN_HPP
#include "mars_simplex.hpp"



namespace mars{




template<typename Elem>
class Jacobian;

template<Integer Dim, Integer ManifoldDim>
class Jacobian<Simplex<Dim, ManifoldDim>>
{
  public:


  Jacobian(const Mesh<Dim,ManifoldDim>&mesh):
  already_set_(false),
  mesh_ptr_(std::make_shared<Mesh<Dim,ManifoldDim>>(mesh))
  {}
  
  constexpr void init_det(){detJ_=det(J_);}

  constexpr void init(const Integer id)
   {
   if(already_set_==false)
    {
      points_.resize(ManifoldDim+1);
      already_set_=true;
    }
   const auto& points=mesh_ptr_->points();
   elem_id_=id;
   elem_=mesh_ptr_->elem(id);
   auto n = n_nodes(elem_);
   points_[0] = points[elem_.nodes[0]];
   for(Integer i = 1; i < n; ++i) 
      points_[i] = points[elem_.nodes[i]];
   jacobian(elem_,points,J_);
   init_det();
   }

  constexpr auto & operator()()const {return J_;}
  constexpr auto & get_det()   const {return detJ_;}
  constexpr auto & elem_id()   const {return elem_id_;}
  constexpr auto & points()    const {return points_;}

  private:  
  bool already_set_;
  Integer elem_id_;
  Simplex<Dim, ManifoldDim> elem_;
  std::vector<Vector<Real,Dim>> points_;
  Matrix<Real, Dim, ManifoldDim> J_;
  Real detJ_;
  std::shared_ptr<Mesh<Dim,ManifoldDim>> mesh_ptr_;
};



    // DOT OF DOTS BETWEEN FQPVALUES AND FQPVALUES (same type)
    template<Integer Dim>
    inline constexpr Real dotofdots(const Vector<Real, Dim> &left, const Vector<Real,Dim> &right)
    {
        return dot(left,right);
    }

    template<Integer Rows,Integer Cols>
    inline constexpr Real dotofdots(const Matrix<Real, Rows,Cols> &left, const Matrix<Real,Rows,Cols> &right)
    {
        Real result = left(0,0)*right(0,0);
        
        for(Integer col = 1; col < Cols; ++col) 
            result += left(0,col)*right(0,col);
        for(Integer col = 0; col < Cols; ++col)
           for(Integer row = 1; row < Rows; ++row)  
            result += left(row,col)*right(row,col);  

        return result;
    }


    template<Integer Dim, typename T>
    inline constexpr Real dotofdots(const Vector<T, Dim> &left, const Vector<T, Dim> &right)
    {
        Real ret = dotofdots(left(0),right(0));
        for(Integer d = 1; d < Dim; ++d) {
            ret += dotofdots(left(d),right(d));
        }

        return ret;
    }

    template<Integer Dim, typename T>
    inline constexpr Real dotofdots(const QPValues<T, Dim> &left, const Vector<T, Dim> &right)
    {
        Real ret = dotofdots(left[0],right[0]);
        for(Integer d = 1; d < Dim; ++d) {
            ret += dotofdots(left[d],right[d]);
        }

        return ret;
    }

    template<typename T, Integer Rows,Integer Cols>
    inline constexpr Real dotofdots(const Matrix<T, Rows,Cols> &left, const Matrix<T,Rows,Cols> &right)
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
    template<Integer Dim, typename T>
    inline constexpr Real dotofdots(const T &left, const Vector<T, Dim> &right)
    {
        Real ret = dotofdots(left,right(0));
        for(Integer d = 1; d < Dim; ++d) {
            ret += dotofdots(left,right(d));
        }

        return ret;
    }

    template<Integer Dim, typename T>
    inline constexpr Real dotofdots(const Vector<T, Dim> &left, const T &right)
    {
        Real ret = dotofdots(left(0),right);
        for(Integer d = 1; d < Dim; ++d) {
            ret += dotofdots(left(d),right);
        }

        return ret;
    }

    template<typename T, Integer Rows,Integer Cols>
    inline constexpr Real dotofdots(const T &left, const Matrix<T,Rows,Cols> &right)
    {
        Real result = left*right(0,0);
        
        for(Integer col = 1; col < Cols; ++col) 
            result += left*right(0,col);
        for(Integer col = 0; col < Cols; ++col)
           for(Integer row = 1; row < Rows; ++row)  
            result += left*right(row,col);  

        return result;
    }

    template<typename T, Integer Rows,Integer Cols>
    inline constexpr Real dotofdots(const Matrix<T, Rows,Cols> &left, const T &right)
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



template<typename TestTrialSpaces, typename Left,typename Right,Integer QR,typename Form>
class LocalTensor<false,TestTrialSpaces,L2DotProductIntegral<Left,Right,QR>,Form>
// template<typename TestTrialSpaces, typename MeshT, typename Left,typename Right,Integer QR,typename Form>
// class LocalTensor<false,TestTrialSpaces,L2DotProductIntegral<MeshT,Left,Right,QR>,Form>
{
    static_assert("negative weights not permitted");
};


// template<typename TestTrialSpaces, typename MeshT, typename Left,typename Right,Integer QR,typename Form>
// class LocalTensor<true,TestTrialSpaces,L2DotProductIntegral<MeshT,Left,Right,QR>,Form>
template<typename TestTrialSpaces, typename Left,typename Right,Integer QR>
class LocalTensor<true,TestTrialSpaces,L2DotProductIntegral<Left,Right,QR>,Number<1>>
{
 public:
 using type= L2DotProductIntegral<Left,Right,QR>;
 using QRule=typename type::QRule;
 using EvalLeft=Evaluation<Expression<Left>,QRule>;
 using EvalRight=Evaluation<Expression<Right>,QRule>;
 using subtype= OperatorType<type,Number<1>>;
 
 LocalTensor(){}

 LocalTensor(const type& expr)
 :
 // eval_left_<QRule>(Eval(expr.left())),
 // eval_right_<QRule>(Eval(expr.right()))
 eval_left_(expr.left()),
 eval_right_(expr.right())
 {}


 template<typename Elem, typename T,Integer NQPoints, Integer NComponents>
 void apply_aux(subtype& vec, const Jacobian<Elem>& J, 
               const FQPValues<T,NQPoints,NComponents> & fqp_values,
               const  QPValues<T,NQPoints> & qp_values)
 {
  const auto& detJ=J.get_det();
  for(Integer ii=0;ii<vec.size();ii++)
    {   

        // std::cout<<"(ii,jj)=("<<ii<<","<<jj<<")"<<std::endl;
        // left value is trial, right value is test
      // std::cout<<ii<<std::endl;
      // std::cout<<"left="<<left_value_.size()<<std::endl;
      // std::cout<<"right="<<right_value_.size()<<std::endl;
      // std::cout<<"right[0]="<<right_value_[0].size()<<std::endl;

        vec[ii]=detJ*dotofdots(qp_values,fqp_values[ii]);
    }
 }

 template<typename Elem, typename T,Integer NQPoints, Integer NComponents>
 void apply_aux(subtype& vec, const Jacobian<Elem>& J, 
               const  QPValues<T,NQPoints> & qp_values,
               const FQPValues<T,NQPoints,NComponents> & fqp_values)
 {
  const auto& detJ=J.get_det();
  for(Integer ii=0;ii<vec.size();ii++)
    {   

        // std::cout<<"(ii,jj)=("<<ii<<","<<jj<<")"<<std::endl;
      //   // left value is trial, right value is test
      // std::cout<<ii<<std::endl;
      // std::cout<<"left="<<left_value_.size()<<std::endl;
      // std::cout<<"right="<<right_value_.size()<<std::endl;
      // std::cout<<"right[0]="<<right_value_[0].size()<<std::endl;

        vec[ii]=detJ*dotofdots(qp_values,fqp_values[ii]);
    }
 }

 template<typename Elem,typename ShapeFunctions>
 void apply(subtype& vec, const Jacobian<Elem>& J, const ShapeFunctions& shape_functions)
 {

  std::cout<<"jacobian tensor left1 "<<std::endl;
  eval_left_.apply(left_value_,J,shape_functions);
  std::cout<<"jacobian tensor right1 "<<std::endl;
  eval_right_.apply(right_value_,J,shape_functions);
  std::cout<<"left_value_"<<std::endl;
  std::cout<<left_value_<<std::endl;
  std::cout<<"right_value_"<<std::endl;
  std::cout<<right_value_<<std::endl;

   apply_aux(vec,J,left_value_,right_value_);

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
template<typename TestTrialSpaces, typename Left,typename Right,Integer QR>
class LocalTensor<true,TestTrialSpaces,L2DotProductIntegral<Left,Right,QR>,Number<2>>
{
 public:
 using type= L2DotProductIntegral<Left,Right,QR>;
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
 void apply(subtype& mat, const Jacobian<Elem>& J, const ShapeFunctions& shape_functions)
 {

  const auto& detJ=J.get_det();
  // decltype(eval_left_) eee(6);
  eval_left_.apply(left_value_,shape_functions);
  eval_right_.apply(right_value_,shape_functions);

  std::cout<<"left_value_"<<std::endl;
  std::cout<<left_value_<<std::endl;
  std::cout<<"right_value_"<<std::endl;
  std::cout<<right_value_<<std::endl;
  std::cout<<"quiiiiiiiii 111"<<std::endl;
  for(Integer ii=0;ii<mat.rows();ii++)
    for(Integer jj=0;jj<mat.cols();jj++)
    {   

        // std::cout<<"(ii,jj)=("<<ii<<","<<jj<<")"<<std::endl;
        // left value is trial, right value is test
        mat(ii,jj)=detJ*dotofdots(right_value_[ii],left_value_[jj]);
    }
    std::cout<<detJ<<std::endl;
    std::cout<<"mat="<<mat<<std::endl;
   std::cout<<"quiiiiiiiii 222"<<std::endl;
 }

  
 private:
 EvalLeft eval_left_;
 EvalRight eval_right_;
 OperatorType<Left,QRule> left_value_;
 OperatorType<Right,QRule> right_value_;

};





template<typename ...Ts>
class Constructor;

template<typename T>
class Constructor< UnaryPlus<Expression<T>> >  ;

template<typename T>
class Constructor< UnaryMinus<Expression<T>> >  ;

template<typename Left,typename Right>
class Constructor< Addition<Expression<Left>,Expression<Right>> >  ;

template<typename Left,typename Right>
class Constructor< Subtraction<Expression<Left>,Expression<Right>> >  ;

template<typename Left,typename Right>
class Constructor< Multiplication<Expression<Left>,Expression<Right>> >  ;




template<typename T>
class Constructor<T>
                                   
{
public:
   template<typename...Inputs>
   static auto apply(const Inputs&...inputs)
   {
    return T(inputs...);
    }

};




template<typename Left1,typename Right1,Integer QR1>
class Constructor<L2DotProductIntegral<Left1,Right1,QR1>>
                                   
{
public:
   using T=L2DotProductIntegral<Left1,Right1,QR1>;

   template<typename...Inputs>
   static auto apply(const Inputs&...inputs)
   {
    return L2Inner(Constructor<typename T::Left>::apply(inputs...), Constructor<typename T::Right>::apply(inputs...));
    }

};

// template<typename MeshT,typename Left1,typename Right1,Integer QR1,
//                         typename Left2,typename Right2,Integer QR2>
// class Constructor<Addition<Expression<L2DotProductIntegral<MeshT,Left1,Right1,QR1>>,
//                                 Expression<L2DotProductIntegral<MeshT,Left2,Right2,QR2>> > >
                                   
// {
// public:
//    using Left=L2DotProductIntegral<MeshT,Left1,Right1,QR1>;
//    using Right=L2DotProductIntegral<MeshT,Left2,Right2,QR2>;

//    static auto apply(const MeshT& mesh)
//    {
//     auto e1=L2Inner(mesh, typename Left::Left(), typename Left::Right());
//     auto e2=L2Inner(mesh, typename Right::Left(), typename Right::Right());
//    return e1+e2;
//     }

// };
// template<typename Left1,typename Right1,Integer QR1,
//          typename Left2,typename Right2,Integer QR2>
// class Constructor<Addition<Expression<L2DotProductIntegral<Left1,Right1,QR1>>,
//                                 Expression<L2DotProductIntegral<Left2,Right2,QR2>> > >
                                   
// {
// public:
//    using Left=L2DotProductIntegral<Left1,Right1,QR1>;
//    using Right=L2DotProductIntegral<Left2,Right2,QR2>;

//    template<typename...Inputs>
//    static auto apply(const Inputs&...inputs)
//    {
//     auto e1=L2Inner(Constructor<typename Left::Left>(inputs...), Constructor<typename Left::Right>(inputs...));
//     auto e2=L2Inner(Constructor<typename Right::Left>(inputs...), Constructor<typename Right::Right>(inputs...));
//    return e1+e2;
//     }

// };
                         


// template<typename Left,typename Right>
// class Constructor< Addition<Expression<Left>,Expression<Right>> >                             
// {
// public:
//    template<typename MeshT>
//    static auto apply(const MeshT& mesh)
//    {
//     auto e1=Constructor<Left>::apply(mesh);
//     auto e2=Constructor<Right>::apply(mesh);
//    return e1+e2;
//     }

// };

template<typename Type>
class Constructor< UnaryPlus<Expression<Type>> >                             
{
public:
   template<typename...Inputs>
   static auto apply(const Inputs&...inputs)
   {
    return Constructor<Type>::apply(inputs...);
   }

};

template<typename Type>
class Constructor< UnaryMinus<Expression<Type>> >                             
{
public:
   template<typename...Inputs>
   static auto apply(const Inputs&...inputs)
   {
    return Constructor<Type>::apply(inputs...);
   }

};

template<typename Left,typename Right>
class Constructor< Addition<Expression<Left>,Expression<Right>> >                             
{
public:
   template<typename...Inputs>
   static auto apply(const Inputs&...inputs)
   {
    auto e1=Constructor<Left>::apply(inputs...);
    auto e2=Constructor<Right>::apply(inputs...);
   return e1+e2;
    }

};

template<typename Left,typename Right>
class Constructor< Subtraction<Expression<Left>,Expression<Right>> >                             
{
public:
   template<typename...Inputs>
   static auto apply(const Inputs&...inputs)
   {
    auto e1=Constructor<Left>::apply(inputs...);
    auto e2=Constructor<Right>::apply(inputs...);
   return e1-e2;
    }

};

template<typename Left,typename Right>
class Constructor< Multiplication<Expression<Left>,Expression<Right>> >                             
{
public:
   template<typename...Inputs>
   static auto apply(const Inputs&...inputs)
   {
    auto e1=Constructor<Left>::apply(inputs...);
    auto e2=Constructor<Right>::apply(inputs...);
   return e1*e2;
    }

};


// template<typename T>
// constexpr auto Construct(const UnaryPlus<Expression<T>>& x )                             
// {return Construct(x.derived());};

// template<typename T>
// constexpr auto Construct(const UnaryPlus<Expression<T>>& x )                             
// {return Construct(x.derived());};

// template<typename T>
// constexpr auto Construct(const UnaryMinus<Expression<T>>& x )                             
// {return -Construct(x.derived());};

// template<typename Left,typename Right>
// constexpr auto Construct(const Subtraction<Expression<Left>,Expression<Right>>& x )                             
// {return Construct(x.left())-Construct(x.right());};

// template<typename Left,typename Right>
// constexpr auto Construct(const Addition<Expression<Left>,Expression<Right>>& x )                             
// {return Construct(x.left())+Construct(x.right());};

// template<typename Left,typename Right>
// constexpr auto Construct(const Multiplication<Expression<Left>,Expression<Right>>& x )                             
// {return Construct(x.left())*Construct(x.right());};







template<typename Tuple,Integer Nmax,Integer N, typename ShapeFunctions>
class EvalOfL2InnersAuxHelper;

template<typename Tuple,Integer Nmax, typename ShapeFunctions>
class EvalOfL2InnersAuxHelper<Tuple,Nmax,Nmax,ShapeFunctions>
{
 public:
    using single_type=Evaluation<Expression<GetType<Tuple,Nmax>>,ShapeFunctions>;
    using type=std::tuple<single_type>;
};

template<typename Tuple,Integer Nmax,Integer N, typename ShapeFunctions>
class EvalOfL2InnersAuxHelper
{
 public:
    using single_type=Evaluation<Expression<GetType<Tuple,N>>,ShapeFunctions>;
    using type=TupleCatType<std::tuple<single_type>,typename EvalOfL2InnersAuxHelper<Tuple,Nmax,N+1,ShapeFunctions>::type >;
};

template<typename Tuple,Integer N, typename ShapeFunctions>
using EvalOfL2InnersAux=typename EvalOfL2InnersAuxHelper<Tuple,TupleTypeSize<Tuple>::value-1,N,ShapeFunctions>::type;

template<typename Tuple, typename ShapeFunctions>
using EvalOfL2InnersType=EvalOfL2InnersAux<Tuple,0,ShapeFunctions>;

// template<typename Tuple,Integer N,typename MeshT, typename ShapeFunctions>
// constexpr std::enable_if_t<(N==TupleTypeSize<Tuple>::value-1), EvalOfL2InnersAux<Tuple,N,ShapeFunctions> > 
// EvalOfL2InnersHelper(const MeshT& mesh, ShapeFunctions& shape_functions)
// {   
//       using ens0=GetType<Tuple,N>;
//       using ens=Evaluation<Expression<ens0>,ShapeFunctions>;
//       return std::tuple<ens>(Eval(Constructor<ens0>::apply(mesh),shape_functions));
// }
template<typename Tuple,Integer N,typename FullSpace,typename ShapeFunctions>
constexpr std::enable_if_t<(N==TupleTypeSize<Tuple>::value-1), EvalOfL2InnersAux<Tuple,N,ShapeFunctions> > 
EvalOfL2InnersHelper(const FullSpace&W,ShapeFunctions& shape_functions)
{   
      using ens0=GetType<Tuple,N>;
      using ens=Evaluation<Expression<ens0>,ShapeFunctions>;
      return std::tuple<ens>(Eval(Constructor<ens0>::apply(W),shape_functions));
}
// template<typename Tuple,Integer N,typename MeshT, typename ShapeFunctions>
// constexpr std::enable_if_t< (N<TupleTypeSize<Tuple>::value-1), EvalOfL2InnersAux<Tuple,N,ShapeFunctions>>  
// EvalOfL2InnersHelper(const MeshT& mesh, ShapeFunctions& shape_functions)
// {
//       using ens0=GetType<Tuple,N>;
//       using ens=Evaluation<Expression<ens0>,ShapeFunctions>;
//       return std::tuple_cat(std::tuple<ens>(Eval(Constructor<ens0>::apply(mesh),shape_functions)),
//                             EvalOfL2InnersHelper<Tuple,N+1>(mesh,shape_functions));
// }

template<typename Tuple,Integer N,typename FullSpace,typename ShapeFunctions>
constexpr std::enable_if_t< (N<TupleTypeSize<Tuple>::value-1), EvalOfL2InnersAux<Tuple,N,ShapeFunctions>>  
EvalOfL2InnersHelper(const FullSpace&W,ShapeFunctions& shape_functions)
{
      using ens0=GetType<Tuple,N>;
      using ens=Evaluation<Expression<ens0>,ShapeFunctions>;
      return std::tuple_cat(std::tuple<ens>(Eval(Constructor<ens0>::apply(W),shape_functions)),
                            EvalOfL2InnersHelper<Tuple,N+1>(W,shape_functions));
}

// template<typename Tuple,typename MeshT, typename ShapeFunctions>
// constexpr auto EvalOfL2Inners(const MeshT& mesh, ShapeFunctions& shape_functions)
// {
//       return EvalOfL2InnersHelper<Tuple,0>(mesh,shape_functions);
// }
template<typename Tuple,typename FullSpace,typename ShapeFunctions>
constexpr auto EvalOfL2Inners(const FullSpace&W,ShapeFunctions& shape_functions)
{
      return EvalOfL2InnersHelper<Tuple,0>(W,shape_functions);
}




template<typename GeneralForm, Integer N >
class LocalTupleOfTensors;

template<typename GeneralForm >
class LocalTupleOfTensors<GeneralForm,1>
{
public:
    // using GeneralForm=typename EvaluationGeneralForm::type;
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
        using type= std::tuple< Vector<Real,Dim> >;
    };

    template<Integer Nmax, Integer N>
    class AuxHelper
    {
        public:
        using Numbers=GetType<TupleOfNumbers,N>;
        static constexpr Integer Dim=Nelem_dofs_array[GetType<Numbers>::value];
        using type=TupleCatType< std::tuple< Vector<Real,Dim>>, typename AuxHelper<Nmax,N+1>::type >;
    };

    using type=typename AuxHelper<TupleTypeSize<TupleOfNumbers>::value-1,0>::type;
    
};



template<typename GeneralForm >
class LocalTupleOfTensors<GeneralForm,2>
{
public:
    // using GeneralForm=typename EvaluationGeneralForm::type;
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


// Evaluation<Expression<GeneralForm<Form>>,ShapeFunctions>
template<typename Form,typename FullSpace,typename ShapeFunctions_>//,typename ShapeFunctions>
class EvaluationOfL2Inners<Evaluation<Expression<GeneralForm<Form,FullSpace>>,ShapeFunctions_>>
{
public: 
    using EvaluationGeneralForm=Evaluation<Expression<GeneralForm<Form,FullSpace>>,ShapeFunctions_>;
    using GeneralForm=typename EvaluationGeneralForm::type;
    using ShapeFunctions=typename EvaluationGeneralForm::ShapeFunctions;
    using LocalTupleOfTensors=typename LocalTupleOfTensors<GeneralForm,GetType<typename GeneralForm::form>::value>::type;    
    using L2Products=typename EvaluationGeneralForm::L2Products;
    using EvalOfL2InnersType=EvalOfL2InnersAux<L2Products,0,ShapeFunctions>;
    

    EvaluationOfL2Inners(const GeneralForm& form,ShapeFunctions& shapefunctions)
    :
    eval_inners_(EvalOfL2Inners<L2Products>(form.space_ptr(),shapefunctions))
    {}

 template<Integer N,typename Output,typename Jacobian>
    constexpr void apply_aux_aux(Output& mat,Jacobian&J)
    {
  std::cout<<"jacobian evalinners1"<<std::endl;
     auto& local_mat=std::get<N>(mat_tuple_);
     local_mat.zero();
     auto & eval_N=std::get<N>(eval_inners_);
     eval_N.apply(local_mat,J);
  std::cout<<"jacobian evalinners2"<<std::endl;

    }   

 template<Integer Nmax,Integer N,typename Output,typename Jacobian>
    constexpr std::enable_if_t<(N==Nmax),void> apply_aux(Output& mat,Jacobian&J)
    {
     apply_aux_aux<N>(mat,J);
    }

 template<Integer Nmax,Integer N,typename Output,typename Jacobian>
    constexpr std::enable_if_t<(N<Nmax),void> apply_aux(Output& mat,Jacobian&J)
    {
     {apply_aux_aux<N>(mat,J);}
      apply_aux<Nmax,N+1>(mat,J);
     }

 template<typename Output,typename Jacobian>
    void apply(Output& mat,Jacobian&J)
    {apply_aux<TupleTypeSize<L2Products>::value-1,0>(mat,J);}

private:
    LocalTupleOfTensors mat_tuple_;
    EvalOfL2InnersType eval_inners_;    
};



}
#endif