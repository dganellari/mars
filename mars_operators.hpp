#ifndef MARS_OPERATORS_HPP
#define MARS_OPERATORS_HPP
#include "mars_base.hpp"
#include "mars_tuple_utilities.hpp"
#include "mars_quadrature_rules.hpp"
#include "mars_expression.hpp"


namespace mars {

////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////---------- OPERATOR IDs ---------- ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

class IdentityOperator: public Expression<IdentityOperator>  
{public: IdentityOperator(){};}; 
class DivergenceOperator: public Expression<DivergenceOperator>
{public: DivergenceOperator(){};};

class GradientOperator : public Expression<GradientOperator> 
{public: GradientOperator(){}};
class CurlOperator: public Expression<CurlOperator>  
{public: CurlOperator(){}};
template<Integer N_>
class ComponentOperator: public Expression<ComponentOperator<N_>>  
 {public: static constexpr Integer N=N_;
  ComponentOperator(){}};

// class Operator{public:
//                static  IdentityOperator identity_op;
//                static  DivergenceOperator divergence_op;
//                static  GradientOperator gradient_op;
//                static  CurlOperator curl_op; 
//                // static  ComponentOperator component_op; 
//                static  IdentityOperator id(){return identity_op;};
//                static  DivergenceOperator div(){return divergence_op;};
//                static  GradientOperator grad(){return gradient_op;};
//                static  CurlOperator curl(){return curl_op;};
//                // static  ComponentOperator component(){return component_op;};
//            };



///////////////////////////////////////////////////////////////
/////////////////////    MAIN OPERATORS    ////////////////////
///////////////////////////////////////////////////////////////

template<typename...Parameters>
class Evaluation;

template<typename...Parameters>
class UnaryPlus;

template<typename...Parameters>
class UnaryMinus;

template<typename...Parameters>
class Multiplication;

template<typename...Parameters>
class Contraction2;

template<typename...Parameters>
class Division;

template<typename...Parameters>
class Subtraction;

template<typename...Parameters>
class Addition;


template<typename T,Integer NQPoints,Integer NComponents>
class FQPValues;

template<typename T,Integer NQPoints>
class QPValues;

template<typename...Operators>
class OperatorTypeHelper;

// template<typename MeshT, typename Left,typename Right,Integer QR>
// class L2DotProductIntegral;
template<typename Left,typename Right,Integer QR>
class L2DotProductIntegral;

template<typename Left, typename Right, typename...Ts>
class IsAddable
{
 public:
 static constexpr bool value=false;
};

template<typename T, Integer Dim,typename...Ts>
class IsAddable<Vector<T,Dim>,Vector<T,Dim>,Ts...>
{
 public:
 static constexpr bool value=true;
};

template<typename T, Integer Rows,Integer Cols,typename...Ts>
class IsAddable<Matrix<T,Rows,Cols>,Matrix<T,Rows,Cols>,Ts...>
{
 public:
 static constexpr bool value=true;
}; 


// template<typename MeshT, typename Left1,typename Right1,Integer QR1,typename Left2,typename Right2,Integer QR2,typename...Ts>
// class IsAddable<L2DotProductIntegral<MeshT,Left1,Right1,QR1>,L2DotProductIntegral<MeshT,Left2,Right2,QR2>,Ts... >
// {
//  public:
//  static constexpr bool value=true;
// };
template<typename Left1,typename Right1,Integer QR1,typename Left2,typename Right2,Integer QR2,typename...Ts>
class IsAddable<L2DotProductIntegral<Left1,Right1,QR1>,L2DotProductIntegral<Left2,Right2,QR2>,Ts... >
{
 public:
 static constexpr bool value=true;
};




// type(T) = type(+T)
template<typename T, typename...Ts>
class OperatorTypeHelper<T,Ts...>
{ public:
  using type=T;
};

template<typename T,Integer Rows, Integer Cols, typename...Ts>
class OperatorTypeHelper< Multiplication< Matrix<T,Rows,Cols>, 
                                          Vector<T,Cols> >, Ts...>
{
  public:
  using type=Vector<T,Rows>;
};

template<typename T,Integer Rows,Integer CommonDim, Integer Cols, typename...Ts>
class OperatorTypeHelper< Multiplication< Matrix<T,Rows,CommonDim>, 
                                          Matrix<T,CommonDim,Cols> >, Ts...>
{
  public:
  using type=Matrix<T,Rows,Cols>;
};

template<typename T, typename S, Integer NQPoints, typename...Ts>
class OperatorTypeHelper< Multiplication< QPValues<T,NQPoints>, 
                                          QPValues<S,NQPoints> >, Ts...>
{
  public:
  using type=QPValues<typename OperatorTypeHelper<Multiplication<T,S>,Ts...>::type ,NQPoints>;
};

template<typename T, typename S, Integer NQPoints,Integer NComponents, typename...Ts>
class OperatorTypeHelper< Multiplication< QPValues<T,NQPoints>, 
                                          FQPValues<S,NQPoints,NComponents> >, Ts...>
{
  public:
  using type=FQPValues<typename OperatorTypeHelper<Multiplication<T,S>,Ts...>::type ,NQPoints,NComponents>;
};

// type(T) = type(+T)
template<typename T, typename...Ts>
class OperatorTypeHelper< UnaryPlus< Expression< T > >, Ts...>
{ public:
  using type=typename OperatorTypeHelper<T,Ts...>::type;
};



// type(T) = type(-T)
template<typename T, typename...Ts>
class OperatorTypeHelper< UnaryMinus< Expression< T > >, Ts...>
{ public:
  using type=typename OperatorTypeHelper<T,Ts...>::type;
};


// type(A)=type(B) = type(A+B)
template<typename Left, typename Right, typename...Ts>
class OperatorTypeHelper< Addition< Expression<Left>, Expression<Right > >, Ts...>
{ public:
  using type1=typename OperatorTypeHelper<Left,Ts...>::type;
  using type2=typename OperatorTypeHelper<Right,Ts...>::type;
  static_assert(IsAddable<typename OperatorTypeHelper<Left>::type,
                          typename OperatorTypeHelper<Right>::type
                         >::value, " In Addition, Left and Right types must be equal");
  // static_assert(IsSame<Left,Right>::value, " In Addition, Left and Right types must be equal");
  using type=typename OperatorTypeHelper<Left,Ts...>::type;
};


// type(A)=type(B) = type(A-B)
template<typename Left, typename Right, typename...Ts>
class OperatorTypeHelper< Subtraction< Expression<Left>, Expression<Right > >, Ts...>
{ public:
  static_assert(IsSame<typename OperatorTypeHelper<Left>::type,
                       typename OperatorTypeHelper<Right>::type
                      >::value, " In Subtraction, Left and Right types must be equal");

  using type=typename OperatorTypeHelper<Left,Ts...>::type;
};

template<typename Left, typename Right, typename...Ts>
class OperatorTypeHelper< Multiplication< Expression<Left>, Expression<Right > >, Ts...>
{ public:
  // static_assert(IsSame<typename OperatorTypeHelper<Left>::type,
  //                      typename OperatorTypeHelper<Right>::type
  //                     >::value, " In Subtraction, Left and Right types must be equal");

  using LeftType=typename OperatorTypeHelper<Left,Ts...>::type;
  using RightType=typename OperatorTypeHelper<Right,Ts...>::type;
  using type=typename OperatorTypeHelper<Multiplication<LeftType,RightType>,Ts...>::type;
};


// type(T) = type(T/REAL)
template<typename T, typename...Ts>
class OperatorTypeHelper< Division< Expression<T>, Real >, Ts...>
{ public:
  using type=typename OperatorTypeHelper<T,Ts...>::type;
};

// type(T) = type(T * REAL)
template<typename T, typename...Ts>
class OperatorTypeHelper< Multiplication< Expression<T>, Real  >, Ts...>
{ public:
  using type=typename OperatorTypeHelper<T,Ts...>::type;
};

// type(T) = type(REAL * T)
template<typename T, typename...Ts>
class OperatorTypeHelper< Multiplication< Real, Expression<T> >, Ts...>
{ public:
  using type=typename OperatorTypeHelper<T,Ts...>::type;
};

template<typename...Inputs>
using OperatorType=typename OperatorTypeHelper<Inputs...>::type;





template<typename...Ts>
class OperatorTypeHelper< Contraction2< Expression<Real>, Expression<Real > >, Ts...>
{ public:
  using type=Real;
};

template<typename T, Integer Dim, typename...Ts>
class OperatorTypeHelper< Contraction2< Expression<Vector<T,Dim>>, Expression<Vector<T,Dim>> >, Ts...>
{ public:
  using type=T;
};

template<typename T, Integer Rows, Integer Cols, typename...Ts>
class OperatorTypeHelper< Contraction2< Expression<Matrix<T,Rows,Cols>>, Expression<Matrix<T,Rows,Cols>> >, Ts...>
{ public:
  using type=T;
};



template<typename Left, typename Right, typename...Ts>
class OperatorTypeHelper< Contraction2< Expression<Left>, Expression<Right > >, Ts...>
{ public:
  using typeleft=typename OperatorTypeHelper<Left,Ts...>::type;
  using typeright=typename OperatorTypeHelper<Right,Ts...>::type;
  using type=typename OperatorTypeHelper<Contraction2<typeleft,typeright>,Ts...>::type;
};

// template<typename Left, typename Right, typename...Ts>
// class OperatorTypeHelper< Contraction2< Expression<Addition<Expression<Left1>,Expression<Right1>> >, 
//                                         Expression<Addition<Expression<Left2>,Expression<Right2>> > >, Ts...>
// { public:
//   using typeleft=typename OperatorTypeHelper<Left,Ts...>::type;
//   using typeright=typename OperatorTypeHelper<Right,Ts...>::type;
//   using type=typename OperatorTypeHelper<Contraction2<typeleft,typeright>,Ts...>::type;
// };

// // type(VECTOR) = type(MATRIX * VECTOR)
// template<typename T, Integer RowsLeft,Integer Common>
// class OperatorType< Multiply< Matrix<T, RowsLeft, Common>, Vector<T, Common> > >
// { public:
//   using type=Vector<T, RowsLeft>;
// };

// // type(MATRIX) = type(MATRIX * MATRIX)

// template<typename T, Integer RowsLeft,Integer Common,Integer ColsRight>
// class OperatorType< Multiply< Matrix<T, RowsLeft, Common>, Matrix<T, Common, ColsRight> > >
// { public:
//   using type=Matrix<T, RowsLeft, ColsRight>;
// };

// template<typename T, Integer RowsLeft,Integer Common,Integer ColsRight, Integer NQPoints, Integer NComponents>
// class OperatorType< Multiply< QPValues < Matrix<T, RowsLeft, Common>, NQPoints> ,
//                                    FQPValues< Matrix<T, Common, ColsRight>, NQPoints, NComponents> > >
// { public:
//    using type=FQPValues< Matrix<T, RowsLeft, ColsRight>, NQPoints, NComponents>;
// };

// template<typename T, Integer RowsLeft,Integer Common,Integer ColsRight, Integer NQPoints>
// class OperatorType< Multiply< QPValues< Matrix<T, RowsLeft, Common>, NQPoints> ,
//                                    QPValues< Matrix<T, Common, ColsRight>, NQPoints> > >
// { public:
//    using type=QPValues< Matrix<T, RowsLeft, ColsRight>, NQPoints>;
// };


// template< typename DerivedLeft,typename DerivedRight, typename  TLeft,typename TRight, Integer NQPoints,Integer NComponents,Integer Dim>
// class OperatorType< Multiply< QPprova < DerivedLeft, TLeft, NQPoints,Dim> ,
//                                    FQPprova< DerivedRight,TRight, NQPoints, NComponents,Dim> > >
// { public:
//    using type= FQPValues<typename OperatorType<Multiply<TLeft,TRight>>::type , NQPoints, NComponents>;
// };

// template< typename DerivedLeft,typename DerivedRight, typename  TLeft,typename TRight, Integer NQPoints,Integer Dim>
// class OperatorType< Multiply< QPprova< DerivedLeft, TLeft,  NQPoints,Dim> ,
//                                    QPprova< DerivedRight,TRight, NQPoints,Dim> > >
// { public:
//    using type= QPValues<typename OperatorType<Multiply<TLeft,TRight>>::type , NQPoints>;
// };


class Contraction
{
public:
 template<typename T,Integer Rows,Integer Cols>
 T operator()(const Matrix<T, Rows,Cols> &A,const Matrix<T, Rows,Cols> &B)
 {
       T result=0;
       for(Integer i = 0; i < Rows; ++i) 
         for(Integer j = 0; j < Cols; ++j) 
             result += A(i, j) * B(i,j);
       return result;
 }

 template<typename T,Integer Rows,Integer Cols>
 Matrix<T,1,Rows> operator()(const Matrix<T, Rows,Cols> &A,const Matrix<T, 1,Cols> &B)
 {
       Matrix<T,1,Rows> result;
       result.zero();
       for(Integer i = 0; i < Rows; ++i) 
         for(Integer j = 0; j < Cols; ++j) 
             result(0,i) += A(i, j) * B(0,j);
       return result;
 }
 template<typename T,Integer Rows,Integer Cols>
 Matrix<T,Rows,1> operator()(const Matrix<T, Rows,Cols> &A,const Matrix<T, Cols,1> &B)
 {
       Matrix<T,Rows,1> result=0;
       for(Integer i = 0; i < Rows; ++i) 
         for(Integer j = 0; j < Cols; ++j) 
             result(i,0) += A(i, j) * B(j,0);
       return result;
 }

 template<typename T>
 T operator()(const Matrix<T, 1,1> &A,const T &B)
 {
       return A(0,0)*B;
 }

 template<typename T>
 T operator()(const T &A,const Matrix<T, 1,1> &B)
 {
       return A*B(0,0);
 }

 template<typename T,Integer Dim>
 T operator()(const Vector<T, Dim> &A,const Vector<T, Dim> &B)
 {
       T result=0;
       for(Integer i = 0; i < Dim; ++i) 
             result += A[i] * B[i];
       return result;
 }

 Real operator()(const Real &A,const Real &B)
 {
       return A*B;
 }

 template<typename T,Integer Rows>
 T operator()(const Matrix<T, Rows,1> &A,const Vector<T, Rows> &B)
 {
       T result=0;
       for(Integer i = 0; i < Rows; ++i) 
             result += A(i, 0) * B[i];
       return result;
 }

 template<typename T,Integer Cols>
 T operator()(const Matrix<T, 1,Cols> &A,const Vector<T, Cols> &B)
 {
       T result=0;
       for(Integer i = 0; i < Cols; ++i) 
             result += A(0,i) * B[i];
       return result;
 }

 template<typename T,Integer Rows>
 T operator()(const Vector<T, Rows> &B,const Matrix<T, Rows,1> &A)
 {
       T result=0;
       for(Integer i = 0; i < Rows; ++i) 
             result += A(i, 0) * B[i];
       return result;
 }

 template<typename T,Integer Cols>
 T operator()(const Vector<T, Cols> &B,const Matrix<T, 1,Cols> &A)
 {
       T result=0;
       for(Integer i = 0; i < Cols; ++i) 
             result += A(0,i) * B[i];
       return result;
 }

};




































template<typename T,Integer NQPoints>
class QPValues;

template<typename T,Integer NQPoints,Integer NComponents>
class FQPValues;

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////			  MATRIX EXPRESSION ALGEBRAIC OPERATIONS		    ///////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename T,Integer Rows,Integer Cols>
class UnaryPlus<Matrix<T,Rows,Cols>> 
{    
 public:       
 inline static void apply(Matrix<T,Rows,Cols>& A)
	{
	 for(Integer ii=0;ii<Rows;ii++)
	 	for(Integer jj=0;jj<Cols;jj++)
	 	   	   A(ii,jj)=A(ii,jj);
	};
};


template<typename T,Integer Rows,Integer Cols>
class UnaryMinus<Matrix<T,Rows,Cols>> 
{    
 public:       
 inline static void apply(Matrix<T,Rows,Cols>& A)
  {
   for(Integer ii=0;ii<Rows;ii++)
    for(Integer jj=0;jj<Cols;jj++)
           A(ii,jj)=-A(ii,jj);
  };
};


template<typename T,Integer Rows,Integer Cols>
class Addition<Matrix<T,Rows,Cols>> 
{       
 public: 

 inline static void apply(      Matrix<T,Rows,Cols>& A,
                          const Matrix<T,Rows,Cols>& B)
  {
   for(Integer ii=0;ii<Rows;ii++)
    for(Integer jj=0;jj<Cols;jj++)
    {
     A(ii,jj)=A(ii,jj)+B(ii,jj);
    }
  };

 inline static void apply(      Matrix<T,Rows,Cols>& A,
                          const Matrix<T,Rows,Cols>& B,
                          const Matrix<T,Rows,Cols>& C)
  {
   for(Integer ii=0;ii<Rows;ii++)
    for(Integer jj=0;jj<Cols;jj++)
    {
     A(ii,jj)=B(ii,jj)+C(ii,jj);
    }
  };
};

template<typename T,Integer Rows,Integer Cols>
class Subtraction<Matrix<T,Rows,Cols>> 
{       
 public:     

 inline static void apply(      Matrix<T,Rows,Cols>& A,
                          const Matrix<T,Rows,Cols>& B)
  {
   for(Integer ii=0;ii<Rows;ii++)
    for(Integer jj=0;jj<Cols;jj++)
    {
     A(ii,jj)=A(ii,jj)-B(ii,jj);
    }
  };

 inline static void apply(      Matrix<T,Rows,Cols>& A,
                          const Matrix<T,Rows,Cols>& B,
                          const Matrix<T,Rows,Cols>& C)
  {
   for(Integer ii=0;ii<Rows;ii++)
    for(Integer jj=0;jj<Cols;jj++)
    {
     A(ii,jj)=B(ii,jj)-C(ii,jj);
    }
  };
};

template<Integer Rows,Integer CommonDim,Integer Cols>
class Multiplication<Matrix<Real,Rows,CommonDim>,
                     Matrix<Real,CommonDim,Cols>> 
{    
 public:            
 inline static void apply(      Matrix<Real,Rows,Cols>& A,
                          const Matrix<Real,Rows,CommonDim>& B,
                          const Matrix<Real,CommonDim,Cols>& C)
  {
   for(Integer ii=0;ii<Rows;ii++)
    for(Integer jj=0;jj<Cols;jj++)
       {
        A(ii,jj)=B(ii,0)*C(0,jj);
        for(Integer cc=1;cc<CommonDim;cc++)
           A(ii,jj)+=B(ii,cc)*C(cc,jj);
       }
  };
};



// template<typename T,Integer Rows,Integer CommonDim,Integer Cols>
// class Multiplication<Matrix<T,Rows,CommonDim>,
//                      Matrix<T,CommonDim,Cols>,
//                      Matrix<T,Rows,Cols> > 
// {    
//  public:            
//  inline static void apply(      Matrix<T,Rows,Cols>& A,
//                           const Matrix<T,Rows,CommonDim>& B,
//  	                        const Matrix<T,CommonDim,Cols>& C)
// 	{
// 	 for(Integer ii=0;ii<Rows;ii++)
// 	 	for(Integer jj=0;jj<Cols;jj++)
// 	 	   {
// 	 	   	A(ii,jj)=B(ii,0)*C(0,jj);
// 	 	   	for(Integer cc=1;cc<CommonDim;cc++)
// 	 	   	   A(ii,jj)+=B(ii,cc)*C(cc,jj);
// 	 	   }
// 	};
// };

template<typename S, typename T,Integer NQPoints>
class Multiplication<QPValues<S,NQPoints>,
                     QPValues<T,NQPoints>> 
{    
 public:            
 inline static void apply(      QPValues<OperatorType<Multiplication<S,T>>,NQPoints>& A,
                          const QPValues<S,NQPoints>& B,
                          const QPValues<T,NQPoints>& C)
  {
    std::cout<<"apply QPValues-QPValues"<<std::endl;
    for(Integer qp=0;qp<NQPoints;qp++)
       {
        //FIXME TODO
        Multiply(A[qp],B[qp],C[qp]);
        std::cout<<"B[qp]"<<B[qp]<<std::endl;
        std::cout<<"C[qp]"<<C[qp]<<std::endl;
        std::cout<<"A[qp]"<<A[qp]<<std::endl;
       }
  };
};

template<typename S, typename T,Integer NQPoints,Integer Ndofs>
class Multiplication<QPValues<T,NQPoints>,
                     FQPValues<S,NQPoints,Ndofs>
                     >
                      
{    
 public:            
 inline static void apply(      FQPValues<OperatorType<Multiplication<S,T>>,NQPoints,Ndofs>& A,
                          const QPValues<S,NQPoints>& B,
                          const FQPValues<T,NQPoints,Ndofs>& C)
  {
   for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
    for(Integer qp=0;qp<NQPoints;qp++)
       {
           //FIXME TODO
        //FIXME TODO
        Multiply(A[n_dof][qp],B[qp],C[n_dof][qp]);
        std::cout<<"B[qp]"<<B[qp]<<std::endl;
        std::cout<<"C[qp]"<<C[n_dof][qp]<<std::endl;
        std::cout<<"A[qp]"<<A[n_dof][qp]<<std::endl;
       }
  };
};

template<typename S, typename T,Integer NQPoints,Integer Ndofs>
class Multiplication<FQPValues<S,NQPoints,Ndofs>,
                     QPValues<T,NQPoints>> 
{    
 public:          
 inline static void apply(      FQPValues<OperatorType<Multiplication<S,T>>,NQPoints,Ndofs>& A,
                          const FQPValues<S,NQPoints,Ndofs>& B,
                          const QPValues<T,NQPoints>& C)
  {
   for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
    for(Integer qp=0;qp<NQPoints;qp++)
       {
        //FIXME TODO
        //FIXME TODO
        Multiply(A[n_dof][qp],B[qp][n_dof],C[qp]);
        std::cout<<"B[qp]"<<B[qp][n_dof]<<std::endl;
        std::cout<<"C[qp]"<<C[qp]<<std::endl;
        std::cout<<"A[qp]"<<A[n_dof][qp]<<std::endl;
       }
  };
};


template<typename T,Integer Rows,Integer Cols>
class Multiplication<Real,Matrix<T,Rows,Cols>> 
{       
 public:               

 // inline static void apply(Matrix<T,Rows,Cols>& A,const Real& alpha)
 //  {
 //   for(Integer ii=0;ii<Rows;ii++)
 //    for(Integer jj=0;jj<Cols;jj++)
 //      A(ii,jj)=A(ii,jj)*alpha;
 //  };

 inline static void apply(      Matrix<T,Rows,Cols>& A,
                          const Matrix<T,Rows,Cols>& B,
                          const Real& alpha)
  {
   for(Integer ii=0;ii<Rows;ii++)
    for(Integer jj=0;jj<Cols;jj++)
      A(ii,jj)=B(ii,jj)*alpha;
  };

};

template<typename T,Integer Rows,Integer Cols>
class Multiplication<Matrix<T,Rows,Cols>,Real> 
{       
 public:               

 // inline static void apply(Matrix<T,Rows,Cols>& A,const Real& alpha)
 //  {
 //   Multiplication<Real,Matrix<T,Rows,Cols>>(A,alpha); 
 //  };

 inline static void apply(      Matrix<T,Rows,Cols>& A,
                          const Matrix<T,Rows,Cols>& B,
                          const Real& alpha)
  {
   for(Integer ii=0;ii<Rows;ii++)
    for(Integer jj=0;jj<Cols;jj++)
      A(ii,jj)=B(ii,jj)*alpha;
  };

};


template<typename T,Integer Rows,Integer Cols>
class Division<Matrix<T,Rows,Cols>,Real> 
{    
 public:            
 inline static void apply(      Matrix<T,Rows,Cols>& A,
                          const Matrix<T,Rows,Cols>& B,
                          const Real& alpha)
  {
   for(Integer ii=0;ii<Rows;ii++)
    for(Integer jj=0;jj<Cols;jj++)
      A(ii,jj)=B(ii,jj)/alpha;
  };
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////			    QP   EXPRESSION ALGEBRAIC OPERATIONS		    ///////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


// template<typename U,Integer NQPoints>
// class UnaryPlus<QPValues<U,NQPoints>,
//                 QPValues<U,NQPoints>> 
// {       
//  public:               
//  inline static void apply(const QPValues<U,NQPoints>& A,
//                                 QPValues<U,NQPoints>& C)
// 	{
// 	 for(Integer qp=0;qp<NQPoints;qp++)
//         UnaryPlus<U,U>::apply(A[qp],C[qp]);
// 	};

// };

// template<typename U,Integer NQPoints>
// class UnaryMinus<QPValues<U,NQPoints>,
//                  QPValues<U,NQPoints>> 
// {       
//  public:               
//  inline static void apply(const QPValues<U,NQPoints>& A,
//                                 QPValues<U,NQPoints>& C)
// 	{
// 	 for(Integer qp=0;qp<NQPoints;qp++)
//         UnaryMinus<U,U>::apply(A[qp],C[qp]);
// 	};

// };

// template<typename U,Integer NQPoints>
// class Addition<QPValues<U,NQPoints>,
//                QPValues<U,NQPoints>,
//                QPValues<U,NQPoints> > 
// {       
//  public:               
//  inline static void apply(const QPValues<U,NQPoints>& A,
//  	                      const QPValues<U,NQPoints>& B, 
//                                 QPValues<U,NQPoints>& C)
// 	{
// 	 for(Integer qp=0;qp<NQPoints;qp++)
//         Addition<U,U,U>::apply(A[qp],B[qp],C[qp]);
// 	};

// };

// template<typename U,Integer NQPoints>
// class Subtraction<QPValues<U,NQPoints>,
//                   QPValues<U,NQPoints>,
//                   QPValues<U,NQPoints> > 
// {       
//  public:               
//  inline static void apply(const QPValues<U,NQPoints>& A,
//  	                      const QPValues<U,NQPoints>& B, 
//                                 QPValues<U,NQPoints>& C)
// 	{
// 	 for(Integer qp=0;qp<NQPoints;qp++)
//         Subtraction<U,U,U>::apply(A[qp],B[qp],C[qp]);
// 	};

// };

// template<typename U, Integer NQPoints>
// class Multiply<QPValues<U,NQPoints>,
//                Real,
//                QPValues<U,NQPoints> > 
// {       
//  public:               

//  inline static void apply(const QPValues<U,NQPoints>& A,
//  	                      const Real& B, 
//                                 QPValues<U,NQPoints>& C)
// 	{
// 	 for(Integer qp=0;qp<NQPoints;qp++)
//         Multiply<U,Real,U>::apply(A[qp],B,C[qp]);
// 	};

// };

// template<typename U,Integer NQPoints>
// class Multiply<Real,
//                QPValues<U,NQPoints>,
//                QPValues<U,NQPoints> > 
// {       
//  public:               

//  inline static void apply(const Real& A,
//  	                      const QPValues<U,NQPoints>& B, 
//                                 QPValues<U,NQPoints>& C)
// 	{
// 	 for(Integer qp=0;qp<NQPoints;qp++)
//         Multiply<Real,U,U>::apply(A,B[qp],C[qp]);
// 	};

// };

// template<typename U,typename V, Integer NQPoints>
// class Multiply<QPValues<U,NQPoints>,
//                QPValues<V,NQPoints>,
//                QPValues<typename OperatorType<Multiply<U,V>>::type,NQPoints> > 
// {       
//  public:               
//  using W=typename OperatorType<Multiply<U,V>>::type;

//  inline static void apply(const QPValues<U,NQPoints>& A,
//  	                      const QPValues<V,NQPoints>& B, 
//                                 QPValues<W,NQPoints>& C)
// 	{
// 	 for(Integer qp=0;qp<NQPoints;qp++)
//         Multiply<U,V,W>::apply(A[qp],B[qp],C[qp]);
// 	};

// };

// template<typename U, Integer NQPoints>
// class Divide<QPValues<U,NQPoints>,
//                Real,
//                QPValues<U,NQPoints> > 
// {       
//  public:               

//  inline static void apply(const QPValues<U,NQPoints>& A,
//  	                      const Real& B, 
//                                 QPValues<U,NQPoints>& C)
// 	{
// 	 for(Integer qp=0;qp<NQPoints;qp++)
//         Divide<U,Real,U>::apply(A[qp],B,C[qp]);
// 	};

// };


///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////			    FQP   EXPRESSION ALGEBRAIC OPERATIONS		    ///////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename U, Integer NQPoints,Integer Ndofs>
class UnaryPlus< FQPValues<U,NQPoints,Ndofs>> 
{       
 public:               
 inline static void apply(FQPValues<U,NQPoints,Ndofs>& A)
  {
for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
   for(Integer qp=0;qp<NQPoints;qp++)
    {
      UnaryPlus<U>::apply(A[n_dof][qp]);
    }
  };
};

template<typename U, Integer NQPoints,Integer Ndofs>
class UnaryMinus< FQPValues<U,NQPoints,Ndofs>> 
{       
 public:               
 inline static void apply(FQPValues<U,NQPoints,Ndofs>& A)
	{
for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
	 for(Integer qp=0;qp<NQPoints;qp++)
	 	{
	 		UnaryMinus<U>::apply(A[n_dof][qp]);
	 	}
	};
};

template<typename U, Integer NQPoints,Integer Ndofs>
class Addition< FQPValues<U,NQPoints,Ndofs>> 
{       
 public:               
 inline static void apply(      FQPValues<U,NQPoints,Ndofs>& A,
                          const FQPValues<U,NQPoints,Ndofs>& B,
 	                        const FQPValues<U,NQPoints,Ndofs>& C)
	{
  for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
	 for(Integer qp=0;qp<NQPoints;qp++)
	 	{
	 		Addition<U>::apply(A[n_dof][qp],B[n_dof][qp],C[n_dof][qp]);
	 	}
	};


};

template<typename U, Integer NQPoints,Integer Ndofs>
class Subtraction< FQPValues<U,NQPoints,Ndofs>> 
{       
 public:               
 inline static void apply(      FQPValues<U,NQPoints,Ndofs>& A,
                          const FQPValues<U,NQPoints,Ndofs>& B,
                          const FQPValues<U,NQPoints,Ndofs>& C)
  {
  for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
   for(Integer qp=0;qp<NQPoints;qp++)
    {
      Subtraction<U>::apply(A[n_dof][qp],B[n_dof][qp],C[n_dof][qp]);
    }
  };
};

template<typename U, typename V, typename W, Integer NQPoints,Integer Ndofs>
class Multiplication<FQPValues<U,NQPoints,Ndofs>,
                FQPValues<V,NQPoints,Ndofs>,
                FQPValues<W,NQPoints,Ndofs> > 
{       
 public:               

 inline static void apply(      FQPValues<U,NQPoints,Ndofs>& A,
                          const FQPValues<V,NQPoints,Ndofs>& B,
                          const FQPValues<W,NQPoints,Ndofs>& C)
  {
  for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
   for(Integer qp=0;qp<NQPoints;qp++)
      Multiplication<U,V,W>::apply(A[n_dof][qp],B[n_dof][qp],C[n_dof][qp]);
  };

};

template<typename U,Integer NQPoints,Integer Ndofs>
class Multiplication< Real,FQPValues<U,NQPoints,Ndofs>> 
{       
 public:               
 inline static void apply(FQPValues<U,NQPoints,Ndofs>& A,const Real& alpha)
	{
  for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
   for(Integer qp=0;qp<NQPoints;qp++)
	 		Multiplication<Real,U>::apply(A[n_dof][qp],alpha);
	};
};


template<typename U,Integer NQPoints,Integer Ndofs>
class Multiplication<FQPValues<U,NQPoints,Ndofs>,Real> 
{       
 public:               
 inline static void apply(FQPValues<U,NQPoints,Ndofs>& A,const Real& alpha)
  {
    Multiplication< Real,FQPValues<U,NQPoints,Ndofs>>::apply(A,alpha);
  };
};

template<typename U,Integer NQPoints,Integer Ndofs>
class Division<FQPValues<U,NQPoints,Ndofs>,Real> 
{       
 public:               
 inline static void apply(FQPValues<U,NQPoints,Ndofs>& A,const Real& alpha)
  {
  for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
   for(Integer qp=0;qp<NQPoints;qp++)
    Division<U,Real>::apply(A[n_dof][qp],alpha);
    
  };
};



///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////				    QP-FQP    ALGEBRAIC OPERATIONS			    ///////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////



// template<typename U,typename V, Integer NQPoints,Integer NComponents>
// class Multiply< QPValues<U,NQPoints>,
//                 FQPValues<V,NQPoints,NComponents>,
//                 FQPValues<typename OperatorType<Multiply<U,V>>::type,NQPoints,NComponents> > 
// {       
//  public:               
//  using W=typename OperatorType<Multiply<U,V>>::type;

//  inline static void apply(const QPValues <U,NQPoints>& A,
//  	                      const FQPValues<V,NQPoints,NComponents>& B, 
//                                 FQPValues<W,NQPoints,NComponents>& C)
// 	{
// 	 for(Integer qp=0;qp<NQPoints;qp++)
// 	 	for(Integer n_comp=0;n_comp<NComponents;n_comp++)
// 	 	{
// 	 		Multiply<U,V,W>::apply(A[qp],B[n_comp][qp],C[n_comp][qp]);
// 	 	}
// 	 std::cout<<"qp * fqp C=="<<C<<std::endl;
// 	};

// };

// template<typename U,typename V, Integer NQPoints,Integer NComponents>
// class Multiply< FQPValues<U,NQPoints,NComponents>,
//                 QPValues<V,NQPoints>,
//                 FQPValues<typename OperatorType<Multiply<U,V>>::type,NQPoints,NComponents> > 
// {       
//  public:               
//  using W=typename OperatorType<Multiply<U,V>>::type;

//  inline static void apply(const FQPValues <U,NQPoints,NComponents>& A,
//  	                      const QPValues<V,NQPoints>& B, 
//                                 FQPValues<W,NQPoints,NComponents>& C)
// 	{
// 	 for(Integer qp=0;qp<NQPoints;qp++)
// 	 	for(Integer n_comp=0;n_comp<NComponents;n_comp++)
// 	 	{
// 	 		Multiply<U,V,W>::apply(A[qp],B[n_comp][qp],C[n_comp][qp]);
// 	 	}
// 	};

// };


template<typename T>
constexpr void Add2(T& t){UnaryPlus<T>::apply(t);}

template<typename T>
constexpr void Minus2(T& t){UnaryMinus<T>::apply(t);}

template<typename T>
constexpr void Add(T& a, const T&b, const T&c){Addition<T>::apply(a,b,c);}

template<typename T>
constexpr void Subtract(T& a, const T&b, const T&c){Subtraction<T>::apply(a,b,c);}

template<typename U,typename V,typename W>
constexpr void Multiply(U& u, const V& v, const W& w){Multiplication<V,W>::apply(u,v,w);}


template<typename T>
constexpr void Multiply(T& a, const Real& alpha){Multiplication<Real,T>::apply(a,alpha);}

template<typename T>
constexpr void Multiply(const Real& alpha,T& a){Multiplication<Real,T>::apply(a,alpha);}

template<typename T>
constexpr void Divide(T& a, const Real& alpha){Division<T,Real>::apply(a,alpha);}







///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////// 				ALGEBRAIC EXPRESSION OPERATIONS 			///////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

template< typename Derived>
class UnaryPlus< Expression <Derived> > 
: public Expression< UnaryPlus< Expression <Derived> > >
{
  public:
    UnaryPlus(const Expression<Derived>& expr): value_(expr.derived()){};
    Derived& operator()(){return value_;};
    const Derived& operator()()const{return value_;};
  private:
  Derived value_;
};

template< typename Derived>
class UnaryMinus< Expression <Derived> > 
: public Expression< UnaryMinus< Expression <Derived> > >
{
  public:
    UnaryMinus(const Expression<Derived>& expr): value_(expr.derived()){};
    Derived& operator()(){return value_;};
    const Derived& operator()()const{return value_;};
  private:
  Derived value_;
};


template< typename DerivedLeft,typename DerivedRight>
class Addition< Expression <DerivedLeft>, Expression <DerivedRight> > 
: public Expression< Addition< Expression <DerivedLeft>, Expression <DerivedRight> > >
{
  public:
    using Left=DerivedLeft;
    using Right=DerivedRight;
    Addition(const Expression<Left>& left, const Expression<Right>&right)
    : 
    left_(left.derived()),
    right_(right.derived())
    {};
    const Left& left()const{return left_;};
    const Right& right()const{return right_;};
          auto operator()(){return Addition<Expression<Left>,Expression<Right>>(left_,right_);};
          auto operator()()const{return Addition<Expression<Left>,Expression<Right>>(left_,right_);};
  private:
  Left left_;
  Right right_;
};


template< typename DerivedLeft,typename DerivedRight>
class Subtraction< Expression <DerivedLeft>, Expression <DerivedRight> > 
: public Expression< Subtraction< Expression <DerivedLeft>, Expression <DerivedRight> > >
{
  public:
    using Left=DerivedLeft;
    using Right=DerivedRight;
    Subtraction(const Expression<DerivedLeft>& left, const Expression<DerivedRight>&right)
    : 
    left_(left.derived()),
    right_(right.derived())
    {};
    const DerivedLeft& left()const{return left_;};
    const DerivedRight& right()const{return right_;};
  private:
  DerivedLeft left_;
  DerivedRight right_;
};

template< typename DerivedLeft,typename DerivedRight>
class Division< Expression <DerivedLeft>, Expression <DerivedRight> > 
: public Expression< Division< Expression <DerivedLeft>, Expression <DerivedRight> > >
{
  public:
    using Left=DerivedLeft;
    using Right=DerivedRight;
    Division(const Expression<DerivedLeft>& left, const Expression<DerivedRight>&right)
    : 
    left_(left.derived()),
    right_(right.derived())
    {};
    const DerivedLeft& left()const{return left_;};
    const DerivedRight& right()const{return right_;};
  private:
  DerivedLeft left_;
  DerivedRight right_;
};

template< typename DerivedLeft>
class Division< Expression <DerivedLeft>, Real > 
: public Expression< Division< Expression <DerivedLeft>, Real > >
{
  public:
    using Left=DerivedLeft;
    using Right=Real;
    Division(const Expression<DerivedLeft>& left, const Real&right)
    : 
    left_(left.derived()),
    right_(right)
    {};
    const DerivedLeft& left()const{return left_;};
    const Real& right()const{return right_;};
  private:
  DerivedLeft left_;
  Real right_;
};


template< typename DerivedLeft,typename DerivedRight>
class Contraction2< Expression <DerivedLeft>, Expression <DerivedRight> > 
: public Expression< Contraction2< Expression <DerivedLeft>, Expression <DerivedRight> > >
{
  public:
    using Left=DerivedLeft;
    using Right=DerivedRight;

    Contraction2(const Expression<DerivedLeft>& left, const Expression<DerivedRight>&right)
    : 
    left_(left.derived()),
    right_(right.derived())
    {
      // static_assert(std::is_same<DerivedLeft,DerivedRight>::value,"in Contraction2, type(Left)==type(Right)");
    };
    const DerivedLeft& left()const{return left_;};
    const DerivedRight& right()const{return right_;};
  private:
  DerivedLeft left_;
  DerivedRight right_;
};

template< typename DerivedLeft,typename DerivedRight>
class Multiplication< Expression <DerivedLeft>, Expression <DerivedRight> > 
: public Expression< Multiplication< Expression <DerivedLeft>, Expression <DerivedRight> > >
{
  public:
    using Left=DerivedLeft;
    using Right=DerivedRight;


    Multiplication(const Expression<DerivedLeft>& left, const Expression<DerivedRight>&right)
    : 
    left_(left.derived()),
    right_(right.derived())
    {};
    const DerivedLeft& left()const{return left_;};
    const DerivedRight& right()const{return right_;};
  private:
  DerivedLeft left_;
  DerivedRight right_;
};


template< typename DerivedLeft>
class Multiplication< Expression <DerivedLeft>, Real > 
: public Expression< Multiplication< Expression <DerivedLeft>, Real > >
{
  public:
    using Left=DerivedLeft;
    using Right=Real;

    Multiplication(const Expression<DerivedLeft>& left, const Real&right)
    : 
    left_(left.derived()),
    right_(right)
    {};
    const DerivedLeft& left()const{return left_;};
    const Real& right()const{return right_;};

  private:
  DerivedLeft left_;
  Real right_;
};

template< typename DerivedRight>
class Multiplication< Real,Expression <DerivedRight> > 
: public Expression< Multiplication<Real, Expression <DerivedRight> > >
{
  public:
    using Left=Real;
    using Right=DerivedRight;
    Multiplication(const Real& left, const Expression <DerivedRight>&right)
    : 
    left_(left),
    right_(right.derived())
    {};
    const Real& left()const{return left_;};
    const DerivedRight& right()const{return right_;};
  private:
  DerivedRight right_;
  Real left_;
};



template< typename DerivedLeft,typename DerivedRight>
class Contraction2< Expression <DerivedLeft>, Expression <DerivedRight> >
Inner(const Expression<DerivedLeft>&left, const Expression<DerivedRight>&right)
{return Contraction2< Expression <DerivedLeft>, Expression <DerivedRight> >(left,right);}




template< typename DerivedLeft>
class Division< Expression <DerivedLeft>, Real > 
operator/(const Expression<DerivedLeft>&left, const Real&right)
{return Division< Expression <DerivedLeft>, Real > (left,right);}


template< typename DerivedLeft,typename DerivedRight>
class Multiplication< Expression <DerivedLeft>, Expression <DerivedRight> >
operator*(const Expression<DerivedLeft>&left, const Expression<DerivedRight>&right)
{return Multiplication< Expression <DerivedLeft>, Expression <DerivedRight> >(left,right);}

template< typename DerivedLeft>
class Multiplication< Expression <DerivedLeft>, Real > 
operator*(const Expression<DerivedLeft>&left, const Real&right)
{return Multiplication< Expression <DerivedLeft>, Real > (left,right);}

template< typename DerivedRight>
class Multiplication< Real, Expression <DerivedRight> >
operator*(const Real&left, const Expression<DerivedRight>&right)
{return Multiplication< Real, Expression <DerivedRight> >(left,right);}




template< typename Derived>
class UnaryPlus< Expression <Derived> >
operator+(const Expression<Derived>& expr)
{return UnaryPlus< Expression <Derived> >(expr);}

template< typename Derived>
class UnaryMinus< Expression <Derived> >
operator-(const Expression<Derived>& expr)
{return UnaryMinus< Expression <Derived> >(expr);}

template< typename DerivedLeft,typename DerivedRight>
class Addition< Expression <DerivedLeft>, Expression <DerivedRight> >
operator+(const Expression<DerivedLeft>&left, const Expression<DerivedRight>&right)
{return Addition< Expression <DerivedLeft>, Expression <DerivedRight> >(left,right);}

template< typename DerivedLeft,typename DerivedRight>
class Subtraction< Expression <DerivedLeft>, Expression <DerivedRight> >
operator-(const Expression<DerivedLeft>&left, const Expression<DerivedRight>&right)
{return Subtraction< Expression <DerivedLeft>, Expression <DerivedRight> >(left,right);}






template<typename T, typename...Ts>
class MultipleAdditionHelper;

template<typename T, typename...Ts>
using MultipleAddition=typename MultipleAdditionHelper<T,Ts...>::type;

template<typename T>
class MultipleAdditionHelper<T>
{
  public:
    using type=T;
};


template<typename T, typename...Ts>
class MultipleAdditionHelper
{
  public:
    using type=Addition<Expression<T>,Expression<typename MultipleAdditionHelper<Ts...>::type>>;
};



///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////// 				ALGEBRAIC EXPRESSION EVALUATIONS 			///////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename Derived,typename...OtherTemplateArguments>
class Evaluation<Expression<UnaryPlus< Expression<Derived> > >,OtherTemplateArguments... >
{
 public:
 using type=UnaryPlus<Expression<Derived>>;
 using Eval=Evaluation<Expression<Derived>,OtherTemplateArguments...>;
 using subtype= OperatorType<type,OtherTemplateArguments...>;
 Evaluation(){};
 

 Evaluation(const Expression<type>& expr):
 expr_(expr.derived()),
 eval_(Eval(expr_()))
 {};
 
 template<typename...Inputs>
 void apply(subtype& output,const Inputs&...inputs)
 {
  // compute output
  // eval_.template apply<OtherTemplateArguments2...>(output,inputs...);

  eval_.apply(output,inputs...);
  // apply plus to output (we assume it does not change, so we do not compute)
  // Derived(output);
 }
private:

 type expr_;
 Eval eval_;
};


template<typename Derived,typename...OtherTemplateArguments>
class Evaluation<Expression<UnaryMinus< Expression<Derived> > >,OtherTemplateArguments...>
{
 public:
 using type=UnaryMinus<Expression<Derived>>;
 using Eval=Evaluation<Expression<Derived>,OtherTemplateArguments...>;
 using subtype= OperatorType<type,OtherTemplateArguments...>;
 Evaluation(){};
 
 Evaluation(const Expression<type>& expr):
 expr_(expr.derived()),
 eval_(Eval(expr_()))
 {};
 
 template<typename...Inputs>
 void apply(subtype& output,const Inputs&...inputs)
 {
  
  // compute output
    // eval_.template apply<OtherTemplateArguments2...>(output,inputs...);
  eval_.apply(output,inputs...);

  Minus2(output);
  // apply minus to output
  // Minus(output,tmp);
 }

private:

 type expr_;
 Eval eval_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Specialization for -(-X)=X. The evaluation is simply X.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename Derived,typename...OtherTemplateArguments>
class Evaluation<Expression<UnaryMinus< Expression<UnaryMinus< Expression<Derived> >> > >,OtherTemplateArguments... >
{
 public:
 using type=UnaryPlus<Expression<Derived>>;
 using Eval=Evaluation<Expression<Derived>,OtherTemplateArguments...>;
 using subtype= OperatorType<type,OtherTemplateArguments...>;
 Evaluation(){};
 

 Evaluation(const Expression<type>& expr):
 expr_(expr.derived()),
 eval_(Eval(expr_()))
 {};
 
 template<typename...Inputs>
 void apply(subtype& output,const Inputs&...inputs)
 {
  // compute output
  // eval_.template apply<OtherTemplateArguments2...>(output,inputs...);

  eval_.apply(output,inputs...);
  // apply plus to output (we assume it does not change, so we do not compute)
  // Derived(output);
 }
private:

 type expr_;
 Eval eval_;
};




template<typename DerivedRight,typename DerivedLeft,typename...OtherTemplateArguments>
class Evaluation<Expression<Addition< Expression<DerivedLeft>  ,  
                                        Expression<DerivedRight>>>,
                 OtherTemplateArguments...>
{
 public:
 using type=Addition<Expression<DerivedLeft>,Expression<DerivedRight>>;
 using EvalLeft=Evaluation<Expression<DerivedLeft>,OtherTemplateArguments...>;
 using EvalRight=Evaluation<Expression<DerivedRight>,OtherTemplateArguments...>;
 using subtype=OperatorType<type,OtherTemplateArguments...>;

 Evaluation(){};
 

 // Evaluation(const Expression<type>& expr):
 Evaluation(const Expression<type>& expr,OtherTemplateArguments&...args):
 expr_(expr.derived()),
 // eval_left_(EvalLeft(expr_.left())),
 // eval_right_(EvalRight(expr_.right()))
 eval_left_(EvalLeft(expr_.left(),args...)),
 eval_right_(EvalRight(expr_.right(),args...))
{};
 
 template<typename...OtherTemplateArguments2,typename...Inputs>
 void apply(subtype& output,const Inputs&...inputs)
 {
  
  eval_left_.apply(left_value_,inputs...);
  eval_right_.apply(right_value_,inputs...);
  Add(output,left_value_,right_value_);
 }
private:

 type expr_;
 subtype left_value_;
 subtype right_value_;
 EvalLeft eval_left_;
 EvalRight eval_right_;
};

template<typename DerivedRight,typename DerivedLeft,typename...OtherTemplateArguments>
class Evaluation<Expression<Subtraction< Expression<DerivedLeft>  ,  
                                           Expression<DerivedRight>>>,
                 OtherTemplateArguments...>
{
 public:
 using type=Subtraction<Expression<DerivedLeft>,Expression<DerivedRight>>;
 using EvalLeft=Evaluation<Expression<DerivedLeft>,OtherTemplateArguments...>;
 using EvalRight=Evaluation<Expression<DerivedRight>,OtherTemplateArguments...>;
 using subtype=OperatorType<type,OtherTemplateArguments...>;
 Evaluation(){};
 

 Evaluation(const Expression<type>& expr):
 expr_(expr.derived()),
 eval_left_(EvalLeft(expr_.left())),
 eval_right_(EvalRight(expr_.right()))
{};
 
 template<typename...OtherTemplateArguments2,typename...Inputs>
 void apply(subtype& output,const Inputs&...inputs)
 {
  eval_left_.apply(left_value_,inputs...);
  eval_right_.apply(right_value_,inputs...);
  Subtract(output,left_value_,right_value_);
 }
private:

 type expr_;
 subtype left_value_;
 subtype right_value_;
 EvalLeft eval_left_;
 EvalRight eval_right_;
};


template<typename DerivedRight,typename DerivedLeft,typename...OtherTemplateArguments>
class Evaluation<Expression<Multiplication< Expression<DerivedLeft>  ,  
                                              Expression<DerivedRight>>>,
                 OtherTemplateArguments...>
{
 public:
 using type=Multiplication<Expression<DerivedLeft>,Expression<DerivedRight>>;
 using EvalLeft=Evaluation<Expression<DerivedLeft>,OtherTemplateArguments...>;
 using EvalRight=Evaluation<Expression<DerivedRight>,OtherTemplateArguments...>;
 using subtype=OperatorType<type,OtherTemplateArguments...>;
 using subtypeleft=OperatorType<DerivedLeft,OtherTemplateArguments...>;
 using subtyperight=OperatorType<DerivedRight,OtherTemplateArguments...>;
 // Evaluation(){};
 

 Evaluation(const Expression<type>& expr):
 expr_(expr.derived()),
 eval_left_(EvalLeft(expr_.left())),
 eval_right_(EvalRight(expr_.right()))
{};

 template<typename...OtherTemplateArguments2,typename...Inputs>
 void apply(subtype& output,const Inputs&...inputs)
 {
  eval_left_.apply(left_value_,inputs...);
  eval_right_.apply(right_value_,inputs...);
  Multiply(output,left_value_,right_value_);
  // std::cout<<"evaluation Multiplication"<<std::endl;
  // std::cout<<"left_value_"<<left_value_<<std::endl;
  // std::cout<<"right_value_"<<right_value_<<std::endl;

 }
private:

 type expr_;
 subtype output;
 subtypeleft left_value_;
 subtyperight right_value_;
 EvalLeft eval_left_;
 EvalRight eval_right_;
};

template<typename DerivedRight,typename...OtherTemplateArguments>
class Evaluation<Expression<Multiplication< Real, Expression<DerivedRight> >>,OtherTemplateArguments...>                                   
{
 public:
 using type=Multiplication<Real,Expression<DerivedRight>>;
 using EvalRight=Evaluation<Expression<DerivedRight>,OtherTemplateArguments...>;
 using subtype=OperatorType<type,OtherTemplateArguments...>;
 Evaluation(){};
 

 Evaluation(const Expression<type>& expr):
 expr_(expr.derived()),
 eval_left_(expr_.left()),
 eval_right_(EvalRight(expr_.right()))
{};
 
 template<typename...OtherTemplateArguments2,typename...Inputs>
 void apply(subtype& output,const Inputs&...inputs)
 {
  // std::cout<<"evaluation Multiplication(Real,)"<<std::endl;
  eval_right_.apply(output,inputs...);
  Multiply(output,eval_left_);
 }
private:

 type expr_;
 Real eval_left_;
 EvalRight eval_right_;
};




template< typename DerivedLeft,typename...OtherTemplateArguments>
class Evaluation<Expression<Multiplication< Expression<DerivedLeft>,Real>>,
                 OtherTemplateArguments...>                                       
{
 public:
 using type=Multiplication<Expression<DerivedLeft>,Real>;
 using EvalLeft=Evaluation<Expression<DerivedLeft>,OtherTemplateArguments...>;
 using subtype=OperatorType<type,OtherTemplateArguments...>;
 
 Evaluation(){};
 

 Evaluation(const Expression<type>& expr):
 expr_(expr.derived()),
 eval_left_(EvalLeft(expr_.left())),
 eval_right_(expr_.right())
{};
 
 template<typename...OtherTemplateArguments2,typename...Inputs>
 void apply(subtype& output,const Inputs&...inputs)
 {
  // std::cout<<"evaluation Multiplication(,Real)"<<std::endl;
  eval_left_.apply(output,inputs...);
  Multiply(output,eval_right_);

 }
private:
 type expr_;
 EvalLeft eval_left_;
 Real eval_right_;
};

template< typename DerivedLeft,typename...OtherTemplateArguments>
class Evaluation<Expression<Division< Expression<DerivedLeft>, Real> >,
                 OtherTemplateArguments...>                                       
{
 public:
 using type=Division<Expression<DerivedLeft>,Real>;
 using EvalLeft=Evaluation<Expression<DerivedLeft>,OtherTemplateArguments...>;
 using subtype=OperatorType<type,OtherTemplateArguments...>;
 
 Evaluation(){};
 

 Evaluation(const Expression<type>& expr):
 expr_(expr.derived()),
 eval_left_(EvalLeft(expr_.left())),
 eval_right_(expr_.right())
{};
 
 template<typename...OtherTemplateArguments2,typename...Inputs>
 void apply(subtype& output,const Inputs&...inputs)
 {
  // std::cout<<"evaluation Division(,Real)"<<std::endl;
  eval_left_.apply(output,inputs...);
  Divide(output,eval_right_);
 }

private:

 type expr_;
 EvalLeft eval_left_;
 Real eval_right_;
};

















// template<typename T,Integer NQPoints,Integer Ndofs1, Integer Ndofs2>
// class Evaluation<Expression<Contraction2<FQPValues<T,NQPoints,Ndofs1>,FQPValues<T,NQPoints,Ndofs2>>>>
// {
//  public:
//  using type=Contraction2<FQPValues<T,NQPoints,Ndofs1>,FQPValues<T,NQPoints,Ndofs2>>;
//  using subtype=OperatorType<type,OtherTemplateArguments...>;

//  Evaluation(){};
 
//  Evaluation(const Expression<type>& expr,OtherTemplateArguments&...args):
//  expr_(expr.derived()),
//  eval_left_(EvalLeft(expr_.left(),args...)),
//  eval_right_(EvalRight(expr_.right(),args...))
//  {};
 
//  template<typename...OtherTemplateArguments2,typename...Inputs>
//  void apply(subtype& output,const Inputs&...inputs)
//  {
//   std::cout<<"evaluation Addition"<<std::endl;

//   eval_left_.apply(left_value_,inputs...);
//   eval_right_.apply(right_value_,inputs...);
//   std::cout<<"left_value_="<<left_value_<<std::endl;
//   std::cout<<"right_value_="<<right_value_<<std::endl;
//   // Add(output,left_value_,right_value_);
//  }
// private:

//  type expr_;
//  subtype left_value_;
//  subtype right_value_;
//  EvalLeft eval_left_;
//  EvalRight eval_right_;
// };

template<typename DerivedRight,typename DerivedLeft,typename...OtherTemplateArguments>
class Evaluation<Expression<Contraction2<Expression<DerivedLeft>,Expression<DerivedRight>>>,
                 OtherTemplateArguments...>
{
 public:
 using type=Contraction2<Expression<DerivedLeft>,Expression<DerivedRight>>;
 using EvalLeft=Evaluation<Expression<DerivedLeft>,OtherTemplateArguments...>;
 using EvalRight=Evaluation<Expression<DerivedRight>,OtherTemplateArguments...>;
 using subtype=OperatorType<type,OtherTemplateArguments...>;

 Evaluation(){};
 
 Evaluation(const Expression<type>& expr,OtherTemplateArguments&...args):
 expr_(expr.derived()),
 eval_left_(EvalLeft(expr_.left(),args...)),
 eval_right_(EvalRight(expr_.right(),args...))
{};
 
 template<typename...OtherTemplateArguments2,typename...Inputs>
 void apply(subtype& output,const Inputs&...inputs)
 {
  // std::cout<<"evaluation Addition"<<std::endl;

  eval_left_.apply(left_value_,inputs...);
  eval_right_.apply(right_value_,inputs...);
  // std::cout<<"left_value_="<<left_value_<<std::endl;
  // std::cout<<"right_value_="<<right_value_<<std::endl;
  // Add(output,left_value_,right_value_);
 }
private:

 type expr_;
 subtype left_value_;
 subtype right_value_;
 EvalLeft eval_left_;
 EvalRight eval_right_;
};



















// // 
// template<typename...Forms>
// class ShapeFunctions;

// // template<typename MeshT, typename Left,typename Right,Integer QR, typename Form>
// // class Evaluation<Expression<L2DotProductIntegral<MeshT,Left,Right,QR>>, ShapeFunctions<Form>>;
// template<typename Left,typename Right,Integer QR, typename Form>
// class Evaluation<Expression<L2DotProductIntegral<Left,Right,QR>>, ShapeFunctions<Form>>;


// template<typename MeshT, typename Left,typename Right,Integer QR, typename Form,typename...OtherTemplateArguments>
// constexpr auto Eval(const L2DotProductIntegral<MeshT,Left,Right,QR>& t,const OtherTemplateArguments&...ts)
// {return Evaluation< Expression<remove_all_t<decltype(t)>>,OtherTemplateArguments...>(t,ts...);}
template<typename Left,typename Right,Integer QR, typename Form,typename...OtherTemplateArguments>
constexpr auto Eval(const L2DotProductIntegral<Left,Right,QR>& t,const OtherTemplateArguments&...ts)
{return Evaluation< Expression<remove_all_t<decltype(t)>>,OtherTemplateArguments...>(t,ts...);}


template<typename...OtherTemplateArguments, typename T>
constexpr auto Eval(const T& t){return Evaluation< Expression<remove_all_t<decltype(t)>>,OtherTemplateArguments...>(t);}

// first parameter is an expression, while the others input are utils 
template<typename...OtherTemplateArguments,typename T,typename ...Ts>
constexpr auto Eval(const T& t,const Ts&...ts){return Evaluation< Expression<remove_all_t<decltype(t)>>,
                                                                              remove_all_t<decltype(ts)>...,OtherTemplateArguments... >(t,ts...);}

template<typename...OtherTemplateArguments,typename T,typename ...Ts>
constexpr auto Eval(const T& t, Ts&...ts){return Evaluation< Expression<remove_all_t<decltype(t)>>,
                                                                         remove_all_t<decltype(ts)>...,OtherTemplateArguments... >(t,ts...);}


template<typename...Args>
class GeneralForm;

template<typename Form,typename...Forms>
class ShapeFunctions2;

template<typename Elem>
class Jacobian;

template<typename...Ts>
class TupleOfTestTrialPairsNumbers;

template<typename TupleOfPairsNumbers, typename Form>
class TupleOfL2Products;

template<typename ...Ts>
class TupleOfEvaluationExpressionOfTypesHelper;

template<typename L2Products,typename Form>
class TupleOfEvals;

template<typename L2Products,typename Form>
class TupleOfEvals<L2Products,ShapeFunctions2<GeneralForm<Form>>>
{
 public:
  using type=typename TupleOfEvaluationExpressionOfTypesHelper<L2Products,ShapeFunctions2<GeneralForm<Form>>>::type;
  


  template<Integer Nmax,Integer N>
  constexpr std::enable_if_t<(Nmax==N),void>
  construct_aux(ShapeFunctions2<GeneralForm<Form>>& shapes)
  {
  }

  template<Integer Nmax,Integer N>
  constexpr std::enable_if_t<(Nmax>N),void>
  construct_aux(ShapeFunctions2<GeneralForm<Form>>& shapes)
  {
  }


  TupleOfEvals(ShapeFunctions2<GeneralForm<Form>>& shapes):
  shapes_(shapes)
  // ,
  // tuple_
  {}

 private:
  // type tuple;
  ShapeFunctions2<GeneralForm<Form>>& shapes_;
  // type tuple_;
};


// template<typename EvaluationGeneralForm,typename ShapeFunctions>
// class EvaluationOfL2Inners;
template<typename...Ts>
class EvaluationOfL2Inners;

template<typename Form,typename GeneralForm_, typename...GeneralForms>
class Evaluation<Expression<GeneralForm<Form>>,ShapeFunctions2<GeneralForm_,GeneralForms...>>
{
 public:
 using type= GeneralForm<Form>;
 using FunctionSpace= typename type::FunctionSpace;
 using TupleFunctionSpace=typename type::TupleFunctionSpace;
 using ShapeFunctions=ShapeFunctions2<GeneralForm_,GeneralForms...>;
 using Shapes=typename ShapeFunctions::TupleOfTupleShapeFunction;
 using TupleOfPairsNumbers=typename GeneralForm<Form>::TupleOfPairsNumbers;
 using L2Products=typename TupleOfL2Products< TupleOfPairsNumbers, Form >::type;
 using TupleOfEvals=TupleOfEvals<L2Products,ShapeFunctions>;
 using EvaluationOfL2Inners=EvaluationOfL2Inners<Evaluation<Expression<GeneralForm<Form>>,ShapeFunctions>>;
 const Form& operator()()const{return general_form_();};

 Evaluation(const GeneralForm<Form>& general_form,ShapeFunctions& shapesform):
 general_form_(general_form),
 shapesform_(shapesform),
 // tuple_eval(faccio<TupleOfPairsNumbers>(general_form_())),
 eval_inners_(general_form_,shapesform)
 {}
 
 template<typename Elem>
 constexpr void apply(const Jacobian<Elem>& J)
 {
  int token_mat=1;
  eval_inners_.apply(token_mat,J);
 }

 private:
 const GeneralForm<Form>& general_form_;
 ShapeFunctions& shapesform_;
 EvaluationOfL2Inners eval_inners_;
 // EvaluationOfL2Inners tuple_eval;
};


template<typename Form,typename FullSpace, typename GeneralForm_, typename...GeneralForms>
constexpr auto Eval(const GeneralForm<Form>& form, ShapeFunctions2<GeneralForm_,GeneralForms...>& shapes)
{return Evaluation< Expression<GeneralForm<Form>>,ShapeFunctions2<GeneralForm_,GeneralForms...> >(form,shapes);}
















// template< typename Derived>
// class Evaluation< UnaryPlus< Expression<Derived> > >
// {   
//  public:
//  using type= typename Evaluation< Derived >::type;

//  Evaluation(const UnaryPlus< Expression<Derived> >& expr):
//  eval_(expr())
//  {};

//  type& operator()(){return self_;};
//  const type& operator()()const{return self_;};

//  template<typename...Parameters>
//  inline type& apply(const Parameters&...parameters)
//  {
//   self_ =eval_.apply(parameters...);
//   UnaryPlus<type,type>::apply(self_,value_);
//   return value_;
//  };
// private:
//  type value_;
//  type self_;
//  Evaluation< Derived > eval_;
// };


// template< typename Derived>
// class Evaluation< UnaryMinus< Expression<Derived> > >
// {   
//  public:
//  using type= typename Evaluation< Derived >::type;

//  Evaluation(const UnaryMinus< Expression<Derived> >& expr):
//  eval_(expr())
//  {};

//  type& operator()(){return self_;};
//  const type& operator()()const{return self_;};

//  template<typename...Parameters>
//  inline type& apply(const Parameters&...parameters)
//  {
//   self_ =eval_.apply(parameters...);
//   UnaryMinus<type,type>::apply(self_,value_);
//   return value_;
//  };
// private:
//  type value_;
//  type self_;
//  Evaluation< Derived > eval_;
// };


// template< typename DerivedLeft,typename DerivedRight>
// class Evaluation< Addition< Expression<DerivedLeft> , Expression<DerivedRight> > >
// {   
//  public:
//  using Left= typename Evaluation< DerivedLeft >::type;
//  using Right=typename Evaluation< DerivedRight >::type;
//  using type= typename OperatorType< Addition<  Left, Right > >::type; // SBBAGLIATO, il prodotto

//  // typename OperatorType< Addition< typename Left::type, typename Right::type > >::type; // SBBAGLIATO, il prodotto

//  Evaluation(const Addition< Expression<DerivedLeft> , Expression<DerivedRight> >& expr):
//  eval_left_(expr.left()),
//  eval_right_(expr.right())
//  {};

//  Left& left(){return left_;};
//  Right& right(){return right_;};
//  const Left& left()const{return left_;};
//  const Right& right()const{return right_;};

//  template<typename...Parameters>
//  inline type& apply(const Parameters&...parameters)
//  {
//   left_ =eval_left_.apply(parameters...);
//   right_=eval_right_.apply(parameters...);
//   Addition<Left,Right,type>::apply(left_,right_,value_);
//   return value_;
//  };
// private:
//  Left left_;
//  Right right_;
//  type value_;
//  Evaluation< DerivedLeft > eval_left_;
//  Evaluation< DerivedRight> eval_right_; 
// };



// template< typename DerivedLeft,typename DerivedRight>
// class Evaluation< Subtraction< Expression<DerivedLeft> , Expression<DerivedRight> > >
// {   
//  public:
//  using Left= typename Evaluation< DerivedLeft >::type;
//  using Right=typename Evaluation< DerivedRight >::type;
//  using type= typename OperatorType< Subtraction<  Left, Right > >::type;
//  // typename OperatorType< Subtraction< typename Left::type, typename Right::type > >::type; // SBBAGLIATO, il prodotto

//  Evaluation(const Subtraction< Expression<DerivedLeft> , Expression<DerivedRight> >& expr):
//  eval_left_(expr.left()),
//  eval_right_(expr.right())
//  {};

//  Left& left(){return left_;};
//  Right& right(){return right_;};
//  const Left& left()const{return left_;};
//  const Right& right()const{return right_;};

//  template<typename...Parameters>
//  inline type& apply(const Parameters&...parameters)
//  {
//   left_ =eval_left_.apply(parameters...);
//   right_=eval_right_.apply(parameters...);
//   Subtraction<Left,Right,type>::apply(left_,right_,value_);
//   return value_;
//  };
// private:
//  Left left_;
//  Right right_;
//  type value_;
//  Evaluation< DerivedLeft > eval_left_;
//  Evaluation< DerivedRight> eval_right_; 
// };


// template< typename DerivedLeft,typename DerivedRight>
// class Evaluation< Multiplication< Expression<DerivedLeft> , Expression<DerivedRight> > >
// {   
//  public:
//  using Left= typename Evaluation< DerivedLeft >::type;
//  using Right=typename Evaluation< DerivedRight >::type;
//  using type= typename OperatorType< Multiply<  Left, Right > >::type;
//  // typename OperatorType< Multiply< typename Left::type, typename Right::type > >::type; // SBBAGLIATO, il prodotto

//  Evaluation(const Multiplication< Expression<DerivedLeft> , Expression<DerivedRight> >& expr):
//  eval_left_(expr.left()),
//  eval_right_(expr.right())
//  {};

//  Left& left(){return left_;};
//  Right& right(){return right_;};
//  const Left& left()const{return left_;};
//  const Right& right()const{return right_;};

//  template<typename...Parameters>
//  inline type& apply(const Parameters&...parameters)
//  {
//   left_ =eval_left_.apply(parameters...);
//   right_=eval_right_.apply(parameters...);
//   Multiply<Left,Right,type>::apply(left_,right_,value_);
//   return value_;
//  };
// private:
//  Left left_;
//  Right right_;
//  type value_;
//  Evaluation< DerivedLeft > eval_left_;
//  Evaluation< DerivedRight> eval_right_; 
// };




// template< typename DerivedLeft>
// class Evaluation< Multiplication< Expression<DerivedLeft> , Real > >
// {   
//  public:
//  using Left= typename Evaluation< DerivedLeft >::type;
//  using Right=Real;
//  using type= typename OperatorType< Multiply<  Left, Right > >::type;
//  // typename OperatorType< Multiply< typename Left::type, Real> >::type; // SBBAGLIATO, il prodotto

//  Evaluation(const Multiplication< Expression<DerivedLeft> , Real >& expr):
//  eval_left_(expr.left()),
//  right_(expr.right())
//  {};

//  Left& left(){return left_;};
//  Real& right(){return right_;};
//  const Left& left()const{return left_;};
//  const Real& right()const{return right_;};

//  template<typename...Parameters>
//  inline type& apply(const Parameters&...parameters)
//  {
//   left_ =eval_left_.apply(parameters...);
//   Multiply<Left,Right,type>::apply(left_,right_,value_);
//   return value_;
//  };
// private:
//  Left left_;
//  Real right_;
//  type value_;
//  Evaluation< DerivedLeft > eval_left_;
// };


// template< typename DerivedRight>
// class Evaluation< Multiplication< Real, Expression<DerivedRight> > >
// {   
//  public:
//  using Left= Real;
//  using Right=typename Evaluation< DerivedRight >::type;
//  using type= typename OperatorType< Multiply<  Left, Right > >::type;
//  // typename OperatorType< Multiply< Real, typename Right::type> >::type; 

//  Evaluation(const Multiplication< Real, Expression<DerivedRight> >& expr):
//  left_(expr.left()),
//  eval_right_(expr.right())
//  {};

//  Left& left(){return left_;};
//  Real& right(){return right_;};
//  const Left& left()const{return left_;};
//  const Real& right()const{return right_;};

//  template<typename...Parameters>
//  inline type& apply(const Parameters&...parameters)
//  {
//   right_ =eval_right_.apply(parameters...);
//   Multiply<Left,Right,type>::apply(left_,right_,value_);
//   return value_;
//  };
// private:
//  Real left_;
//  Right right_;
//  type value_;
//  Evaluation< DerivedRight > eval_right_;
// };

// template< typename DerivedLeft>
// class Evaluation< Division< Expression<DerivedLeft> , Real > >
// {   
//  public:
//  using Left= typename Evaluation< DerivedLeft >::type;
//  using Right=Real;
//  using type= typename OperatorType< Divide<  Left, Right > >::type;
//  // typename OperatorType< Divide< typename Left::type, Real> >::type; // SBBAGLIATO, il prodotto

//  Evaluation(const Division< Expression<DerivedLeft> , Real >& expr):
//  eval_left_(expr.left()),
//  right_(expr.right())
//  {};

//  Left& left(){return left_;};
//  Real& right(){return right_;};
//  const Left& left()const{return left_;};
//  const Real& right()const{return right_;};

//  template<typename...Parameters>
//  inline type& apply(const Parameters&...parameters)
//  {
//   left_ =eval_left_.apply(parameters...);
//   Divide<Left,Right,type>::apply(left_,right_,value_);
//   return value_;
//  };
// private:
//  Left left_;
//  Real right_;
//  type value_;
//  Evaluation< DerivedLeft > eval_left_;
// };







template<typename T,Integer NQPoints>
class ExpressionQPValues: 
public Expression< ExpressionQPValues<T,NQPoints > >
{
public:
        using type=QPValues<T,NQPoints>;
        ExpressionQPValues(){};

        ExpressionQPValues(const type& value):value_(value){};

        type& operator()(){return value_;};
        const type& operator()()const{return value_;};
protected:
    type value_;

};

template<typename T,Integer NQPoints>
class Evaluation<ExpressionQPValues<T,NQPoints> >
{
 public:
 using type=QPValues<T,NQPoints>;

 Evaluation(const ExpressionQPValues<T,NQPoints>&expr):
 value_(expr.derived()()){};
 template<typename...Parameters>
 inline type& apply(const Parameters&...parameters)
 {
   return value_;
 };
 private:
 	type value_;
};




template<typename T,Integer NQPoints, Integer NComponents>
class ExpressionFQPValues: 
public Expression< ExpressionFQPValues<T,NQPoints,NComponents > >
{
public:
        using type=FQPValues<T,NQPoints,NComponents>;
        ExpressionFQPValues(){};

        ExpressionFQPValues(const type& value):value_(value){};

        type& operator()(){return value_;};
        const type& operator()()const{return value_;};
protected:
    type value_;

};

template<typename T,Integer NQPoints, Integer NComponents>
class Evaluation<ExpressionFQPValues<T,NQPoints,NComponents> >
{
 public:
 using type=FQPValues<T,NQPoints,NComponents>;

 Evaluation(const ExpressionFQPValues<T,NQPoints,NComponents>&expr):
 value_(expr.derived()()){};
 template<typename...Parameters>
 inline type& apply(const Parameters&...parameters)
 {
   return value_;
 };
 private:
 	type value_;
};




template<typename T,Integer Rows,Integer Cols>
class ExpressionMatrixVar: 
public Expression< ExpressionMatrixVar<T,Rows,Cols> >
{
public:
        using type=Matrix<T,Rows,Cols>;
        ExpressionMatrixVar(){};

        type& operator()(){return value_;};
        const type& operator()()const{return value_;};
protected:
    type value_;

};

template<typename T,Integer Rows,Integer Cols>
class Evaluation<ExpressionMatrixVar<T,Rows,Cols> >
{
 public:
 using type=Matrix<T,Rows,Cols>;

 Evaluation(const ExpressionMatrixVar<T,Rows,Cols>&expr):
 value_(expr.derived()()){};
 template<typename QP>
 inline type& apply(const QP&qp_points)
 {

   for(Integer ii=0;ii<Rows;ii++)	
   	  for(Integer jj=0;jj<Cols;jj++)
   	      value_(ii,jj)	= qp_points(ii,jj);
   return value_;
 };
 private:
 	type value_;
};





template<typename T,Integer Rows,Integer Cols>
class ExpressionMatrix: 
public Expression< ExpressionMatrix<T,Rows,Cols> >
{
public:
        using type=Matrix<T,Rows,Cols>;
        ExpressionMatrix(){};

        ExpressionMatrix(const type& value):value_(value){};

        type& operator()(){return value_;};
        const type& operator()()const{return value_;};
protected:
    type value_;

};

template<typename T,Integer Rows,Integer Cols>
class Evaluation<ExpressionMatrix<T,Rows,Cols> >
{
 public:
 using type=Matrix<T,Rows,Cols>;

 Evaluation(const ExpressionMatrix<T,Rows,Cols>&expr):
 value_(expr.derived()()){};
 template<typename...Parameters>
 inline type& apply(const Parameters&...parameters)
 {
   return value_;
 };
 private:
 	type value_;
};








































































template<typename...T>
class QuadratureOrder;



template<typename Elem, Integer QuadratureRuleType,Integer ActualOrder>
class CheckMaxQuadratureOrder
{
public:
  static constexpr Integer value=Min(MaxOrder<Elem,QuadratureRuleType>::value,ActualOrder);
};


// order(T) = order(+T)
template<typename T>
class QuadratureOrder< UnaryPlus< Expression<T> > >
{ public:
  static constexpr Integer value=QuadratureOrder<T>::value;
};


// order(T) = order(-T)
template<typename T>
class QuadratureOrder< UnaryMinus< Expression<T> > >
{ public:
  static constexpr Integer value=QuadratureOrder<T>::value;
};




// order(A+B) =max(order(A),order(B))
template<typename Left, typename Right>
class QuadratureOrder< Addition< Expression<Left>,Expression<Right> > >
{ public:
  static constexpr Integer value=Max(QuadratureOrder<Left>::value,QuadratureOrder<Right>::value);
};


// order(A+B) =max(order(A),order(B))
template<typename Left, typename Right>
class QuadratureOrder< Subtraction< Expression<Left>,Expression<Right> > >
{ public:
  static constexpr Integer value=Max(QuadratureOrder<Left>::value,QuadratureOrder<Right>::value);
};


// order(T) = order(T/REAL)
template<typename T>
class QuadratureOrder< Division<Expression<T>, Real > >
{ public:
  static constexpr Integer value=QuadratureOrder<T>::value;
};

// order(T) = order(T * REAL)
template<typename T>
class QuadratureOrder< Multiplication< Expression<T>, Real > >
{ public:
  static constexpr Integer value=QuadratureOrder<T>::value;
};

// order(T) = order(REAL * T)
template<typename T>
class QuadratureOrder< Multiplication< Real, Expression<T> > >
{ public:
  static constexpr Integer value=QuadratureOrder<T>::value;
};


// order(Left*Right) = order(Left) * order(Right)
template<typename Left, typename Right>
class QuadratureOrder< Multiplication< Expression<Left>, Expression<Right> > >
{ public:

  static constexpr Integer value=QuadratureOrder<Left>::value + QuadratureOrder<Right>::value;
};


// order(Left*Right) = order(Left) * order(Right)
template<typename Left, typename Right>
class QuadratureOrder< Contraction2< Expression<Left>, Expression<Right> > >
{ public:

  static constexpr Integer value=QuadratureOrder<Left>::value + QuadratureOrder<Right>::value;
};









template<Integer Order,Integer Continuity, Integer NComponents>
class QuadratureOrder<IdentityOperator,BaseFunctionSpace<LagrangeFE,Order,Continuity,NComponents> >
{ public:
  static constexpr Integer value=Order;
};

template<Integer Order,Integer Continuity, Integer NComponents>
class QuadratureOrder<GradientOperator, BaseFunctionSpace<LagrangeFE,Order,Continuity,NComponents> >
{ public:
  static constexpr Integer value=Order-1;
};


template<Integer Order,Integer Continuity, Integer NComponents>
class QuadratureOrder<IdentityOperator, BaseFunctionSpace<RaviartThomasFE,Order,Continuity,NComponents> >
{ public:
  static constexpr Integer value=Order+1;
};



template<Integer Order,Integer Continuity, Integer NComponents>
class QuadratureOrder<DivergenceOperator, BaseFunctionSpace<RaviartThomasFE,Order,Continuity,NComponents> >
{ public:
  static constexpr Integer value=Order;
};







template<Integer N,Integer Nmax,typename...Ts>
class ShapeFunctionExpression: public Expression<ShapeFunctionExpression<N,Nmax,Ts...>>
{
public:
	static constexpr Integer value=N;

};

template<Integer N,Integer Nmax,typename...Ts>
class GradientShapeFunctionExpression: public Expression<GradientShapeFunctionExpression<N,Nmax,Ts...>>
{
public:
	static constexpr Integer value=N;

};

template<typename...T>
class OperatorTupleType;


// order(T) = order(+T)
template<typename T>
class OperatorTupleType< UnaryPlus< Expression<T> > >
{ public:
  using type=typename OperatorTupleType<T>::type;
};


// order(T) = order(-T)
template<typename T>
class OperatorTupleType< UnaryMinus< Expression<T> > >
{ public:
  using type=typename OperatorTupleType<T>::type;
};




// order(A+B) =max(order(A),order(B))
template<typename Left, typename Right>
class OperatorTupleType< Addition< Expression<Left>,Expression<Right> > >
{ public:
  using	LeftT=typename OperatorTupleType<Left>::type;
  using	RightT=typename OperatorTupleType<Right>::type;
  using type=RemoveTupleOfTupleDuplicates< LeftT, RightT >;
};


// order(A+B) =max(order(A),order(B))
template<typename Left, typename Right>
class OperatorTupleType< Subtraction< Expression<Left>,Expression<Right> > >
{ public:
  using	LeftT=typename OperatorTupleType<Left>::type;
  using	RightT=typename OperatorTupleType<Right>::type;
  using type=RemoveTupleOfTupleDuplicates<  LeftT, RightT >;
};


// order(T) = order(T/REAL)
template<typename T>
class OperatorTupleType< Division<Expression<T>, Real > >
{ public:
  using type=typename OperatorTupleType<T>::type;
};

// order(T) = order(T * REAL)
template<typename T>
class OperatorTupleType< Multiplication< Expression<T>, Real > >
{ public:
  using type=typename OperatorTupleType<T>::type;
};

// order(T) = order(REAL * T)
template<typename T>
class OperatorTupleType< Multiplication< Real, Expression<T> > >
{ public:
  using type=typename OperatorTupleType<T>::type;
};


template<typename Left, typename Right>
class OperatorTupleType< Multiplication< Expression<Left>, Expression<Right> > >
{ public:
  using	LeftT=typename OperatorTupleType<Left>::type;
  using	RightT=typename OperatorTupleType<Right>::type;
  using type=RemoveTupleOfTupleDuplicates<  LeftT, RightT >;
};


// order(Left*Right) = order(Left) * order(Right)
template<typename Left, typename Right>
class OperatorTupleType< Contraction2< Expression<Left>, Expression<Right> > >
{ public:
  using LeftT=typename OperatorTupleType<Left>::type;
  using RightT=typename OperatorTupleType<Right>::type;
  using type=RemoveTupleOfTupleDuplicates<  LeftT, RightT >;
};



// order(Left*Right) = order(Left) * order(Right)









/////////////////////// TO REMOVE /////////////////////////////

template<Integer N,Integer Nmax,typename...T>
class OperatorTupleType<ShapeFunctionExpression<N,Nmax,T...> >
{ public:
  using single_type=std::tuple<std::tuple< IdentityOperator,std::tuple<> >>;
  static constexpr Integer value=ShapeFunctionExpression<N,Nmax,T...>::value;
  using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
  using type=TupleChangeType<value,single_type,emptytuple>;

};

template<Integer N,Integer Nmax,typename...T>
class OperatorTupleType<GradientShapeFunctionExpression<N,Nmax,T...> >
{ public:
  using single_type=std::tuple<std::tuple< GradientOperator,std::tuple<> >>;
  static constexpr Integer value=ShapeFunctionExpression<N,Nmax,T...>::value;
  using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
  using type=TupleChangeType<value,single_type,emptytuple>;
};

template<typename...T>
class OperatorTupleType<ExpressionQPValues<T...> >
{ public:
  using type=std::tuple<>;
  static constexpr Integer value=-1;
};

}

#endif //MARS_OPERATORS_HPP
