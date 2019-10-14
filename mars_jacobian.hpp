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



    // BILINEAR FORMS
    // DOT OF DOTS BETWEEN FQPVALUES AND FQPVALUES (same type)
    template<Integer Dim>
    inline constexpr Real dotofdots(const Vector<Real, Dim> &left, const Vector<Real,Dim> &right)
    {
         // std::cout<<"dotofdots const Vector<Real, Dim> &left, const Vector<Real,Dim> &right"<<std::endl;

        return dot(left,right);
    }
    

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

    template<Integer Dim, typename T>
    inline constexpr Real dotofdots(const Vector<T, Dim> &left, const Vector<T, Dim> &right)
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
      std::cout<<" eee"<<std::endl;
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
 {
  // decltype(eval_left_) qui(5,4,4,5,5,5,5,5);
  // decltype(eval_right_) qui2(5,4,4,5,5,5,5,5);
  // std::cout<<"qui"<<std::endl;
 }


 template<typename Elem, typename T,typename S, Integer NQPoints, Integer NComponents>
 void apply_aux(subtype& vec, const Jacobian<Elem>& J, 
               const FQPValues<T,NQPoints,NComponents> & fqp_values,
               const  QPValues<S,NQPoints> & qp_values)
 {
  const auto& detJ=J.get_det();
  for(Integer ii=0;ii<vec.size();ii++)
    {   

        // std::cout<<"(ii,jj)=("<<ii<<","<<jj<<")"<<std::endl;
        // left value is trial, right value is test
      // std::cout<<ii<<std::endl;
      // std::cout<<"apply_auxleft  qp="<<qp_values<<std::endl;
      // std::cout<<"apply_aux right fqp ="<<fqp_values<<std::endl;
      // // std::cout<<"right[0]="<<right_value_[0].size()<<std::endl;
      //    std::cout<<"vec.size()="<<vec.size()<<std::endl;
      //     std::cout<<"before dotofdots="<<qp_values<<std::endl;
      //    std::cout<<"before dotofdots="<<fqp_values[ii]<<std::endl;

        vec[ii]=detJ*dotofdots(qp_values,fqp_values[ii]);
    }
 }

 template<typename Elem, typename T,typename S, Integer NQPoints, Integer NComponents>
 void apply_aux(subtype& vec, const Jacobian<Elem>& J, 
               const  QPValues<T,NQPoints> & qp_values,
               const FQPValues<S,NQPoints,NComponents> & fqp_values)
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
       // QPValues<T,NQPoints> a91(2);
       // FQPValues<S,NQPoints,NComponents> ok;

       // auto es=ok[0];
       // decltype(es) mm(6);
       std::cout<<"apply_auxleft qp ="<<qp_values<<std::endl;
      std::cout<<"apply_aux right fqp="<<fqp_values<<std::endl; 
        vec[ii]=detJ*dotofdots(qp_values,fqp_values[ii]);
    }
 }

 // template<typename Elem,typename...Args1, typename...Args2, typename...Args3>
 // void apply(subtype& vec, const Jacobian<Elem>& J, const std::tuple<Args1...>& tuple_shape_functions, const std::tuple<Args2...>&tuple_tensor, const std::tuple<Args3...>&tuple_evals)
 template<typename Elem,typename...Forms>
  void apply(subtype& vec, const Jacobian<Elem>& J, ShapeFunctions2<Forms...>& shape_functions)

 {
  // Left k(6);
  // OperatorType<Left,QRule> t(6);
  // std::cout<<"----"<<t<<std::endl;

  // todo fixme
  // eval_left_.apply(left_value_,J,shape_functions);
  // eval_right_.apply(right_value_,J,shape_functions);

  // eval_left_.apply(left_value_,J,tuple_shape_functions,tuple_tensor,tuple_evals);
  // eval_right_.apply(right_value_,J,tuple_shape_functions,tuple_tensor,tuple_evals);
  std::cout<<"before eval_left_"<<std::endl;
  // EvalLeft ok(6);
  eval_left_.apply(left_value_,J,shape_functions(),shape_functions.composite_tensor(),shape_functions.composite_shapes());
  std::cout<<"before eval_right_"<<std::endl;
  eval_right_.apply(right_value_,J,shape_functions(),shape_functions.composite_tensor(),shape_functions.composite_shapes());


  std::cout<<"left_value_"<<std::endl;
  std::cout<<left_value_<<std::endl;
  std::cout<<"right_value_"<<std::endl;
  std::cout<<right_value_<<std::endl;
    std::cout<<"before vec=="<<vec<<std::endl;
    std::cout<<left_value_<<std::endl;
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
  // todo fixme
  // eval_left_.apply(left_value_,J,shape_functions);
  // eval_right_.apply(right_value_,J,shape_functions);

  eval_left_.apply(left_value_,J,shape_functions(),shape_functions.composite_tensor(),shape_functions.composite_shapes());
  eval_right_.apply(right_value_,J,shape_functions(),shape_functions.composite_tensor(),shape_functions.composite_shapes());
// shape_functions_.composite_tensor(),shape_functions_.composite_shapes());
  std::cout<<"left_value_"<<std::endl;
  // std::cout<<left_value_<<std::endl;
  std::cout<<"right_value_"<<std::endl;
  // std::cout<<right_value_<<std::endl;
  // std::cout<<"quiiiiiiiii 111"<<std::endl;
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


template<typename Tuple,Integer N,typename FullSpace,typename ShapeFunctions>
constexpr std::enable_if_t<(N==TupleTypeSize<Tuple>::value-1), EvalOfL2InnersAux<Tuple,N,ShapeFunctions> > 
EvalOfL2InnersHelper(const FullSpace&W,ShapeFunctions& shape_functions)
{   
      using ens0=GetType<Tuple,N>;
      using ens=Evaluation<Expression<ens0>,ShapeFunctions>;
      return std::tuple<ens>(Eval(Constructor<ens0>::apply(W),shape_functions));
}


template<typename Tuple,Integer N,typename FullSpace,typename ShapeFunctions>
constexpr std::enable_if_t< (N<TupleTypeSize<Tuple>::value-1), EvalOfL2InnersAux<Tuple,N,ShapeFunctions>>  
EvalOfL2InnersHelper(const FullSpace&W,ShapeFunctions& shape_functions)
{
      using ens0=GetType<Tuple,N>;
      using ens=Evaluation<Expression<ens0>,ShapeFunctions>;
      return std::tuple_cat(std::tuple<ens>(Eval(Constructor<ens0>::apply(W),shape_functions)),
                            EvalOfL2InnersHelper<Tuple,N+1>(W,shape_functions));
}


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


template<typename Form,typename ShapeFunctions_>
class EvaluationOfL2Inners<Evaluation<Expression<GeneralForm<Form>>,ShapeFunctions_>>
{
public: 
    using EvaluationGeneralForm=Evaluation<Expression<GeneralForm<Form>>,ShapeFunctions_>;
    using GeneralForm=typename EvaluationGeneralForm::type;
    using ShapeFunctions=typename EvaluationGeneralForm::ShapeFunctions;
    using TupleOfTupleShapeFunction=typename ShapeFunctions::TupleOfTupleShapeFunction;
    using LocalTupleOfTensors=typename LocalTupleOfTensors<GeneralForm,GetType<typename GeneralForm::form>::value>::type;    
    using L2Products=typename EvaluationGeneralForm::L2Products;
    using EvalOfL2InnersType=EvalOfL2InnersAux<L2Products,0,ShapeFunctions>;
    using TupleOfPairsNumbers=typename GeneralForm::TupleOfPairsNumbers;
    

    EvaluationOfL2Inners(const GeneralForm& form,ShapeFunctions& shapefunctions)
    :
    eval_inners_(build_tuple_of_evals<TupleOfPairsNumbers>(form(),shapefunctions))
    {}

 template<Integer N,typename Output,typename Jacobian>
    constexpr void apply_aux_aux(Output& mat,Jacobian&J)
    {
    std::cout<<"pre jacobian evalinners2"<<std::endl;
     auto& local_mat=std::get<N>(tensor_tuple_);
     local_mat.zero();
     auto & eval_N=std::get<N>(eval_inners_);
     eval_N.apply(local_mat,J);
  std::cout<<"after jacobian evalinners2"<<std::endl;

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
    {
      apply_aux<TupleTypeSize<L2Products>::value-1,0>(mat,J);
    }

private:
    LocalTupleOfTensors tensor_tuple_;
    EvalOfL2InnersType eval_inners_;    
};



}
#endif