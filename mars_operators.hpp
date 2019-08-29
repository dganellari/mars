#ifndef MARS_OPERATORS_HPP
#define MARS_OPERATORS_HPP
#include "mars_base.hpp"
#include "mars_tuple_utilities.hpp"
#include "mars_quadrature_rules.hpp"

namespace mars {


///////////////////////////////////////////////////////////////////////////////////////
/////////////////////    EXPRESSIONS FOR TEMPLATE EXPRESSIONS      ////////////////////
///////////////////////////////////////////////////////////////////////////////////////
template<typename Derived>
class Expression2
{
public:
        inline Derived &derived() { return static_cast<Derived &>(*this); }
        inline constexpr const Derived &derived() const { return static_cast<const Derived &>(*this); }
        inline operator Derived &() {return derived();}
        inline constexpr  operator const Derived &() const {return derived();}

         ///@return a string with the name of the class
        virtual std::string getClass() const {
            return "Expression of ";
        }
};


///////////////////////////////////////////////////////////////
/////////////////////    MAIN OPERATORS    ////////////////////
///////////////////////////////////////////////////////////////

template<typename...Parameters>
class Evaluation;

template<typename...Parameters>
class UnaryPlus2;

template<typename...Parameters>
class UnaryMinus2;

template<typename...Parameters>
class Multiplication2;

template<typename...Parameters>
class Contraction2;

template<typename...Parameters>
class Division2;

template<typename...Parameters>
class Subtraction2;

template<typename...Parameters>
class Addition2;


template<typename T,Integer NQPoints,Integer NComponents>
class FQPValues;

template<typename T,Integer NQPoints>
class QPValues;

template<typename...Operators>
class OperatorTypeHelper;



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

template<typename MeshT, typename Left,typename Right,Integer QR>
class L2DotProductIntegral;

template<typename MeshT, typename Left1,typename Right1,Integer QR1,typename Left2,typename Right2,Integer QR2,typename...Ts>
class IsAddable<L2DotProductIntegral<MeshT,Left1,Right1,QR1>,L2DotProductIntegral<MeshT,Left2,Right2,QR2>,Ts... >
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

// type(T) = type(+T)
template<typename T, typename...Ts>
class OperatorTypeHelper< UnaryPlus2< Expression2< T > >, Ts...>
{ public:
  using type=typename OperatorTypeHelper<T,Ts...>::type;
};



// type(T) = type(-T)
template<typename T, typename...Ts>
class OperatorTypeHelper< UnaryMinus2< Expression2< T > >, Ts...>
{ public:
  using type=typename OperatorTypeHelper<T,Ts...>::type;
};


// type(A)=type(B) = type(A+B)
template<typename Left, typename Right, typename...Ts>
class OperatorTypeHelper< Addition2< Expression2<Left>, Expression2<Right > >, Ts...>
{ public:
  using type1=typename OperatorTypeHelper<Left,Ts...>::type;
  using type2=typename OperatorTypeHelper<Right,Ts...>::type;
  // static_assert(IsAddable<typename OperatorTypeHelper<Left>::type,
  //                         typename OperatorTypeHelper<Right>::type
  //                        >::value, " In Addition, Left and Right types must be equal");

  using type=typename OperatorTypeHelper<Left,Ts...>::type;
};


// type(A)=type(B) = type(A-B)
template<typename Left, typename Right, typename...Ts>
class OperatorTypeHelper< Subtraction2< Expression2<Left>, Expression2<Right > >, Ts...>
{ public:
  static_assert(IsSame<typename OperatorTypeHelper<Left>::type,
                       typename OperatorTypeHelper<Right>::type
                      >::value, " In Subtraction, Left and Right types must be equal");

  using type=typename OperatorTypeHelper<Left,Ts...>::type;
};

// type(T) = type(T/REAL)
template<typename T, typename...Ts>
class OperatorTypeHelper< Division2< Expression2<T>, Real >, Ts...>
{ public:
  using type=typename OperatorTypeHelper<T,Ts...>::type;
};

// type(T) = type(T * REAL)
template<typename T, typename...Ts>
class OperatorTypeHelper< Multiplication2< Expression2<T>, Real  >, Ts...>
{ public:
  using type=typename OperatorTypeHelper<T,Ts...>::type;
};

// type(T) = type(REAL * T)
template<typename T, typename...Ts>
class OperatorTypeHelper< Multiplication2< Real, Expression2<T> >, Ts...>
{ public:
  using type=typename OperatorTypeHelper<T,Ts...>::type;
};

template<typename...Inputs>
using OperatorType=typename OperatorTypeHelper<Inputs...>::type;





template<typename...Ts>
class OperatorTypeHelper< Contraction2< Expression2<Real>, Expression2<Real > >, Ts...>
{ public:
  using type=Real;
};

template<typename T, Integer Dim, typename...Ts>
class OperatorTypeHelper< Contraction2< Expression2<Vector<T,Dim>>, Expression2<Vector<T,Dim>> >, Ts...>
{ public:
  using type=T;
};

template<typename T, Integer Rows, Integer Cols, typename...Ts>
class OperatorTypeHelper< Contraction2< Expression2<Matrix<T,Rows,Cols>>, Expression2<Matrix<T,Rows,Cols>> >, Ts...>
{ public:
  using type=T;
};



template<typename Left, typename Right, typename...Ts>
class OperatorTypeHelper< Contraction2< Expression2<Left>, Expression2<Right > >, Ts...>
{ public:
  using typeleft=typename OperatorTypeHelper<Left,Ts...>::type;
  using typeright=typename OperatorTypeHelper<Right,Ts...>::type;
  using type=typename OperatorTypeHelper<Contraction2<typeleft,typeright>,Ts...>::type;
};

// template<typename Left, typename Right, typename...Ts>
// class OperatorTypeHelper< Contraction2< Expression2<Addition2<Expression2<Left1>,Expression2<Right1>> >, 
//                                         Expression2<Addition2<Expression2<Left2>,Expression2<Right2>> > >, Ts...>
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
class UnaryPlus2<Matrix<T,Rows,Cols>> 
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
class UnaryMinus2<Matrix<T,Rows,Cols>> 
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
class Addition2<Matrix<T,Rows,Cols>> 
{       
 public:               
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
class Subtraction2<Matrix<T,Rows,Cols>> 
{       
 public:               
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


template<typename T,Integer Rows,Integer CommonDim,Integer Cols>
class Multiplication2<Matrix<T,Rows,CommonDim>,
                Matrix<T,CommonDim,Cols>,
                Matrix<T,Rows,Cols> > 
{    
 public:            
 inline static void apply(      Matrix<T,Rows,Cols>& A,
                          const Matrix<T,Rows,CommonDim>& B,
 	                        const Matrix<T,CommonDim,Cols>& C)
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


template<typename T,Integer Rows,Integer Cols>
class Multiplication2<Real,Matrix<T,Rows,Cols>> 
{       
 public:               

 inline static void apply(Matrix<T,Rows,Cols>& A,const Real& alpha)
  {
   for(Integer ii=0;ii<Rows;ii++)
    for(Integer jj=0;jj<Cols;jj++)
      A(ii,jj)=A(ii,jj)*alpha;
  };
};

template<typename T,Integer Rows,Integer Cols>
class Multiplication2<Matrix<T,Rows,Cols>,Real> 
{       
 public:               

 inline static void apply(Matrix<T,Rows,Cols>& A,const Real& alpha)
  {
   Multiplication2<Real,Matrix<T,Rows,Cols>>(A,alpha); 
  };
};


template<typename T,Integer Rows,Integer Cols>
class Division2<Matrix<T,Rows,Cols>,Real> 
{    
 public:            
 inline static void apply(Matrix<T,Rows,Cols>& A,const Real& alpha)
	{
	 for(Integer ii=0;ii<Rows;ii++)
	 	for(Integer jj=0;jj<Cols;jj++)
	 	   	   A(ii,jj)=A(ii,jj)/alpha;   
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
class UnaryPlus2< FQPValues<U,NQPoints,Ndofs>> 
{       
 public:               
 inline static void apply(FQPValues<U,NQPoints,Ndofs>& A)
  {
for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
   for(Integer qp=0;qp<NQPoints;qp++)
    {
      UnaryPlus2<U>::apply(A[n_dof][qp]);
    }
  };
};

template<typename U, Integer NQPoints,Integer Ndofs>
class UnaryMinus2< FQPValues<U,NQPoints,Ndofs>> 
{       
 public:               
 inline static void apply(FQPValues<U,NQPoints,Ndofs>& A)
	{
for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
	 for(Integer qp=0;qp<NQPoints;qp++)
	 	{
	 		UnaryMinus2<U>::apply(A[n_dof][qp]);
	 	}
	};
};

template<typename U, Integer NQPoints,Integer Ndofs>
class Addition2< FQPValues<U,NQPoints,Ndofs>> 
{       
 public:               
 inline static void apply(      FQPValues<U,NQPoints,Ndofs>& A,
                          const FQPValues<U,NQPoints,Ndofs>& B,
 	                        const FQPValues<U,NQPoints,Ndofs>& C)
	{
  for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
	 for(Integer qp=0;qp<NQPoints;qp++)
	 	{
	 		Addition2<U>::apply(A[n_dof][qp],B[n_dof][qp],C[n_dof][qp]);
	 	}
	};
};

template<typename U, Integer NQPoints,Integer Ndofs>
class Subtraction2< FQPValues<U,NQPoints,Ndofs>> 
{       
 public:               
 inline static void apply(      FQPValues<U,NQPoints,Ndofs>& A,
                          const FQPValues<U,NQPoints,Ndofs>& B,
                          const FQPValues<U,NQPoints,Ndofs>& C)
  {
  for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
   for(Integer qp=0;qp<NQPoints;qp++)
    {
      Subtraction2<U>::apply(A[n_dof][qp],B[n_dof][qp],C[n_dof][qp]);
    }
  };
};

template<typename U, typename V, typename W, Integer NQPoints,Integer Ndofs>
class Multiplication2<FQPValues<U,NQPoints,Ndofs>,
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
      Multiplication2<U,V,W>::apply(A[n_dof][qp],B[n_dof][qp],C[n_dof][qp]);
  };

};

template<typename U,Integer NQPoints,Integer Ndofs>
class Multiplication2< Real,FQPValues<U,NQPoints,Ndofs>> 
{       
 public:               
 inline static void apply(FQPValues<U,NQPoints,Ndofs>& A,const Real& alpha)
	{
  for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
   for(Integer qp=0;qp<NQPoints;qp++)
	 		Multiplication2<Real,U>::apply(A[n_dof][qp],alpha);
	};
};


template<typename U,Integer NQPoints,Integer Ndofs>
class Multiplication2<FQPValues<U,NQPoints,Ndofs>,Real> 
{       
 public:               
 inline static void apply(FQPValues<U,NQPoints,Ndofs>& A,const Real& alpha)
  {
    Multiplication2< Real,FQPValues<U,NQPoints,Ndofs>>::apply(A,alpha);
  };
};

template<typename U,Integer NQPoints,Integer Ndofs>
class Division2<FQPValues<U,NQPoints,Ndofs>,Real> 
{       
 public:               
 inline static void apply(FQPValues<U,NQPoints,Ndofs>& A,const Real& alpha)
  {
  for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
   for(Integer qp=0;qp<NQPoints;qp++)
    Division2<U,Real>::apply(A[n_dof][qp],alpha);
    
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
constexpr void Add2(T& t){UnaryPlus2<T>::apply(t);}

template<typename T>
constexpr void Minus2(T& t){UnaryMinus2<T>::apply(t);}

template<typename T>
constexpr void Add(T& a, const T&b, const T&c){Addition2<T>::apply(a,b,c);}

template<typename T>
constexpr void Subtract(T& a, const T&b, const T&c){Subtraction2<T>::apply(a,b,c);}

template<typename U,typename V,typename W>
constexpr void Multiply(U& u, const V& v, const W& w){Multiplication2<U,V,W>::apply(u,v,w);}


template<typename T>
constexpr void Multiply(T& a, const Real& alpha){Multiplication2<Real,T>::apply(a,alpha);}

template<typename T>
constexpr void Multiply(const Real& alpha,T& a){Multiplication2<Real,T>::apply(a,alpha);}

template<typename T>
constexpr void Divide(T& a, const Real& alpha){Division2<T,Real>::apply(a,alpha);}







///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////// 				ALGEBRAIC EXPRESSION OPERATIONS 			///////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

template< typename Derived>
class UnaryPlus2< Expression2 <Derived> > 
: public Expression2< UnaryPlus2< Expression2 <Derived> > >
{
  public:
    UnaryPlus2(const Expression2<Derived>& expr): value_(expr.derived()){};
    Derived& operator()(){return value_;};
    const Derived& operator()()const{return value_;};
  private:
  Derived value_;
};

template< typename Derived>
class UnaryMinus2< Expression2 <Derived> > 
: public Expression2< UnaryMinus2< Expression2 <Derived> > >
{
  public:
    UnaryMinus2(const Expression2<Derived>& expr): value_(expr.derived()){};
    Derived& operator()(){return value_;};
    const Derived& operator()()const{return value_;};
  private:
  Derived value_;
};


template< typename DerivedLeft,typename DerivedRight>
class Addition2< Expression2 <DerivedLeft>, Expression2 <DerivedRight> > 
: public Expression2< Addition2< Expression2 <DerivedLeft>, Expression2 <DerivedRight> > >
{
  public:
    using Left=DerivedLeft;
    using Right=DerivedRight;
    Addition2(const Expression2<DerivedLeft>& left, const Expression2<DerivedRight>&right)
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
class Subtraction2< Expression2 <DerivedLeft>, Expression2 <DerivedRight> > 
: public Expression2< Subtraction2< Expression2 <DerivedLeft>, Expression2 <DerivedRight> > >
{
  public:
    Subtraction2(const Expression2<DerivedLeft>& left, const Expression2<DerivedRight>&right)
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
class Division2< Expression2 <DerivedLeft>, Expression2 <DerivedRight> > 
: public Expression2< Division2< Expression2 <DerivedLeft>, Expression2 <DerivedRight> > >
{
  public:
    Division2(const Expression2<DerivedLeft>& left, const Expression2<DerivedRight>&right)
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
class Division2< Expression2 <DerivedLeft>, Real > 
: public Expression2< Division2< Expression2 <DerivedLeft>, Real > >
{
  public:
    Division2(const Expression2<DerivedLeft>& left, const Real&right)
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


template< typename DerivedLeft_,typename DerivedRight_>
class Contraction2< Expression2 <DerivedLeft_>, Expression2 <DerivedRight_> > 
: public Expression2< Contraction2< Expression2 <DerivedLeft_>, Expression2 <DerivedRight_> > >
{
  public:
    using DerivedLeft=DerivedLeft_;
    using DerivedRight=DerivedRight_;

    Contraction2(const Expression2<DerivedLeft>& left, const Expression2<DerivedRight>&right)
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

template< typename DerivedLeft_,typename DerivedRight_>
class Multiplication2< Expression2 <DerivedLeft_>, Expression2 <DerivedRight_> > 
: public Expression2< Multiplication2< Expression2 <DerivedLeft_>, Expression2 <DerivedRight_> > >
{
  public:
  	using DerivedLeft=DerivedLeft_;
  	using DerivedRight=DerivedRight_;

    Multiplication2(const Expression2<DerivedLeft>& left, const Expression2<DerivedRight>&right)
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


template< typename DerivedLeft_>
class Multiplication2< Expression2 <DerivedLeft_>, Real > 
: public Expression2< Multiplication2< Expression2 <DerivedLeft_>, Real > >
{
  public:
  	using DerivedLeft=DerivedLeft_;
  	using DerivedRight=Real;
    Multiplication2(const Expression2<DerivedLeft>& left, const Real&right)
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

template< typename DerivedRight_>
class Multiplication2< Real,Expression2 <DerivedRight_> > 
: public Expression2< Multiplication2<Real, Expression2 <DerivedRight_> > >
{
  public:
  	using DerivedLeft=Real;
  	using DerivedRight=DerivedRight_;
    Multiplication2(const Real& left, const Expression2 <DerivedRight>&right)
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
class Contraction2< Expression2 <DerivedLeft>, Expression2 <DerivedRight> >
Inner(const Expression2<DerivedLeft>&left, const Expression2<DerivedRight>&right)
{return Contraction2< Expression2 <DerivedLeft>, Expression2 <DerivedRight> >(left,right);}




template< typename DerivedLeft>
class Division2< Expression2 <DerivedLeft>, Real > 
operator/(const Expression2<DerivedLeft>&left, const Real&right)
{return Division2< Expression2 <DerivedLeft>, Real > (left,right);}


template< typename DerivedLeft,typename DerivedRight>
class Multiplication2< Expression2 <DerivedLeft>, Expression2 <DerivedRight> >
operator*(const Expression2<DerivedLeft>&left, const Expression2<DerivedRight>&right)
{return Multiplication2< Expression2 <DerivedLeft>, Expression2 <DerivedRight> >(left,right);}

template< typename DerivedLeft>
class Multiplication2< Expression2 <DerivedLeft>, Real > 
operator*(const Expression2<DerivedLeft>&left, const Real&right)
{return Multiplication2< Expression2 <DerivedLeft>, Real > (left,right);}

template< typename DerivedRight>
class Multiplication2< Real, Expression2 <DerivedRight> >
operator*(const Real&left, const Expression2<DerivedRight>&right)
{return Multiplication2< Real, Expression2 <DerivedRight> >(left,right);}




template< typename Derived>
class UnaryPlus2< Expression2 <Derived> >
operator+(const Expression2<Derived>& expr)
{return UnaryPlus2< Expression2 <Derived> >(expr);}

template< typename Derived>
class UnaryMinus2< Expression2 <Derived> >
operator-(const Expression2<Derived>& expr)
{return UnaryMinus2< Expression2 <Derived> >(expr);}

template< typename DerivedLeft,typename DerivedRight>
class Addition2< Expression2 <DerivedLeft>, Expression2 <DerivedRight> >
operator+(const Expression2<DerivedLeft>&left, const Expression2<DerivedRight>&right)
{return Addition2< Expression2 <DerivedLeft>, Expression2 <DerivedRight> >(left,right);}

template< typename DerivedLeft,typename DerivedRight>
class Subtraction2< Expression2 <DerivedLeft>, Expression2 <DerivedRight> >
operator-(const Expression2<DerivedLeft>&left, const Expression2<DerivedRight>&right)
{return Subtraction2< Expression2 <DerivedLeft>, Expression2 <DerivedRight> >(left,right);}












///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////// 				ALGEBRAIC EXPRESSION EVALUATIONS 			///////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename Derived,typename...OtherTemplateArguments>
class Evaluation<Expression2<UnaryPlus2< Expression2<Derived> > >,OtherTemplateArguments... >
{
 public:
 using type=UnaryPlus2<Expression2<Derived>>;
 using Eval=Evaluation<Expression2<Derived>,OtherTemplateArguments...>;
 using subtype= OperatorType<type,OtherTemplateArguments...>;
 Evaluation(){};
 

 Evaluation(const Expression2<type>& expr):
 expr_(expr.derived()),
 eval_(Eval(expr_()))
 {};
 
 template<typename...Inputs>
 void apply(subtype& output,const Inputs&...inputs)
 {
  std::cout<<"evaluation unaryplus"<<std::endl;
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
class Evaluation<Expression2<UnaryMinus2< Expression2<Derived> > >,OtherTemplateArguments...>
{
 public:
 using type=UnaryMinus2<Expression2<Derived>>;
 using Eval=Evaluation<Expression2<Derived>,OtherTemplateArguments...>;
 using subtype= OperatorType<type,OtherTemplateArguments...>;
 Evaluation(){};
 
 Evaluation(const Expression2<type>& expr):
 expr_(expr.derived()),
 eval_(Eval(expr_()))
 {};
 
 template<typename...Inputs>
 void apply(subtype& output,const Inputs&...inputs)
 {
  std::cout<<"evaluation unaryminus"<<std::endl;
  
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




template<typename DerivedRight,typename DerivedLeft,typename...OtherTemplateArguments>
class Evaluation<Expression2<Addition2< Expression2<DerivedLeft>  ,  
                                        Expression2<DerivedRight>>>,
                 OtherTemplateArguments...>
{
 public:
 using type=Addition2<Expression2<DerivedLeft>,Expression2<DerivedRight>>;
 using EvalLeft=Evaluation<Expression2<DerivedLeft>,OtherTemplateArguments...>;
 using EvalRight=Evaluation<Expression2<DerivedRight>,OtherTemplateArguments...>;
 using subtype=OperatorType<type,OtherTemplateArguments...>;

 Evaluation(){};
 

 Evaluation(const Expression2<type>& expr):
 expr_(expr.derived()),
 eval_left_(EvalLeft(expr_.left())),
 eval_right_(EvalRight(expr_.right()))
{};
 
 template<typename...OtherTemplateArguments2,typename...Inputs>
 void apply(subtype& output,const Inputs&...inputs)
 {
  std::cout<<"evaluation Addition2"<<std::endl;

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
class Evaluation<Expression2<Subtraction2< Expression2<DerivedLeft>  ,  
                                           Expression2<DerivedRight>>>,
                 OtherTemplateArguments...>
{
 public:
 using type=Subtraction2<Expression2<DerivedLeft>,Expression2<DerivedRight>>;
 using EvalLeft=Evaluation<Expression2<DerivedLeft>,OtherTemplateArguments...>;
 using EvalRight=Evaluation<Expression2<DerivedRight>,OtherTemplateArguments...>;
 using subtype=OperatorType<type,OtherTemplateArguments...>;
 Evaluation(){};
 

 Evaluation(const Expression2<type>& expr):
 expr_(expr.derived()),
 eval_left_(EvalLeft(expr_.left())),
 eval_right_(EvalRight(expr_.right()))
{};
 
 template<typename...OtherTemplateArguments2,typename...Inputs>
 void apply(subtype& output,const Inputs&...inputs)
 {
  std::cout<<"evaluation Subtraction2"<<std::endl;
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
class Evaluation<Expression2<Multiplication2< Expression2<DerivedLeft>  ,  
                                              Expression2<DerivedRight>>>,
                 OtherTemplateArguments...>
{
 public:
 using type=Multiplication2<Expression2<DerivedLeft>,Expression2<DerivedRight>>;
 using EvalLeft=Evaluation<Expression2<DerivedLeft>,OtherTemplateArguments...>;
 using EvalRight=Evaluation<Expression2<DerivedRight>,OtherTemplateArguments...>;
 using subtype=OperatorType<type,OtherTemplateArguments...>;
 using subtypeleft=OperatorType<DerivedLeft,OtherTemplateArguments...>;
 using subtyperight=OperatorType<DerivedRight,OtherTemplateArguments...>;
 Evaluation(){};
 

 Evaluation(const Expression2<type>& expr):
 expr_(expr.derived()),
 eval_left_(EvalLeft(expr_.left())),
 eval_right_(EvalRight(expr_.right()))
{};

 template<typename...OtherTemplateArguments2,typename...Inputs>
 void apply(subtype& output,const Inputs&...inputs)
 {
  std::cout<<"evaluation Multiplication2"<<std::endl;
  eval_left_.apply(left_value_,inputs...);
  eval_right_.apply(right_value_,inputs...);
  Multiply(output,left_value_,right_value_);
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
class Evaluation<Expression2<Multiplication2< Real, Expression2<DerivedRight> >>,OtherTemplateArguments...>                                   
{
 public:
 using type=Multiplication2<Real,Expression2<DerivedRight>>;
 using EvalRight=Evaluation<Expression2<DerivedRight>,OtherTemplateArguments...>;
 using subtype=OperatorType<type,OtherTemplateArguments...>;
 Evaluation(){};
 

 Evaluation(const Expression2<type>& expr):
 expr_(expr.derived()),
 eval_left_(expr_.left()),
 eval_right_(EvalRight(expr_.right()))
{};
 
 template<typename...OtherTemplateArguments2,typename...Inputs>
 void apply(subtype& output,const Inputs&...inputs)
 {
  std::cout<<"evaluation Multiplication2(Real,)"<<std::endl;
  eval_right_.apply(output,inputs...);
  Multiply(output,eval_left_);
 }
private:

 type expr_;
 Real eval_left_;
 EvalRight eval_right_;
};




template< typename DerivedLeft,typename...OtherTemplateArguments>
class Evaluation<Expression2<Multiplication2< Expression2<DerivedLeft>,Real>>,
                 OtherTemplateArguments...>                                       
{
 public:
 using type=Multiplication2<Expression2<DerivedLeft>,Real>;
 using EvalLeft=Evaluation<Expression2<DerivedLeft>,OtherTemplateArguments...>;
 using subtype=OperatorType<type,OtherTemplateArguments...>;
 
 Evaluation(){};
 

 Evaluation(const Expression2<type>& expr):
 expr_(expr.derived()),
 eval_left_(EvalLeft(expr_.left())),
 eval_right_(expr_.right())
{};
 
 template<typename...OtherTemplateArguments2,typename...Inputs>
 void apply(subtype& output,const Inputs&...inputs)
 {
  std::cout<<"evaluation Multiplication2(,Real)"<<std::endl;
  eval_left_.apply(output,inputs...);
  Multiply(output,eval_right_);

 }
private:
 type expr_;
 EvalLeft eval_left_;
 Real eval_right_;
};

template< typename DerivedLeft,typename...OtherTemplateArguments>
class Evaluation<Expression2<Division2< Expression2<DerivedLeft>, Real> >,
                 OtherTemplateArguments...>                                       
{
 public:
 using type=Division2<Expression2<DerivedLeft>,Real>;
 using EvalLeft=Evaluation<Expression2<DerivedLeft>,OtherTemplateArguments...>;
 using subtype=OperatorType<type,OtherTemplateArguments...>;
 
 Evaluation(){};
 

 Evaluation(const Expression2<type>& expr):
 expr_(expr.derived()),
 eval_left_(EvalLeft(expr_.left())),
 eval_right_(expr_.right())
{};
 
 template<typename...OtherTemplateArguments2,typename...Inputs>
 void apply(subtype& output,const Inputs&...inputs)
 {
  std::cout<<"evaluation Division2(,Real)"<<std::endl;
  eval_left_.apply(output,inputs...);
  Divide(output,eval_right_);
 }

private:

 type expr_;
 EvalLeft eval_left_;
 Real eval_right_;
};


// 
template<typename Form>
class ShapeFunctions;

template<typename MeshT, typename Left,typename Right,Integer QR, typename Form>
class Evaluation<Expression2<L2DotProductIntegral<MeshT,Left,Right,QR>>, ShapeFunctions<Form>>;


// template<typename...OtherTemplateArguments, typename T>
// constexpr auto Eval(const T& t){return Evaluation< Expression2<remove_all_t<decltype(t)>>,OtherTemplateArguments...>(t);}

// // first parameter is an expression, while the others input are utils 
// template<typename...OtherTemplateArguments,typename T,typename ...Ts>
// constexpr auto Eval(const T& t,const Ts&...ts){return Evaluation< Expression2<remove_all_t<decltype(t)>>,
//                                                                               remove_all_t<decltype(ts)>...,OtherTemplateArguments... >(t,ts...);}

// template<typename...OtherTemplateArguments,typename T,typename ...Ts>
// constexpr auto Eval(const T& t, Ts&...ts){return Evaluation< Expression2<remove_all_t<decltype(t)>>,
//                                                                          remove_all_t<decltype(ts)>...,OtherTemplateArguments... >(t,ts...);}


template<typename Form>
class GeneralForm;

template<typename Form>
class ShapeFunctions2;




template<typename Form>
class Evaluation<Expression2<GeneralForm<Form>>>
{
 public:

 Evaluation(const GeneralForm<Form>& generalform,ShapeFunctions2<GeneralForm<Form>>& shapes):
 generalform_(generalform),
 shapes_(shapes)
 {}

 private:
 const GeneralForm<Form>& generalform_;
 ShapeFunctions2<GeneralForm<Form>>& shapes_;
};

template<typename Form>
constexpr auto Eval(const GeneralForm<Form>& form,ShapeFunctions2<GeneralForm<Form>>& shapes)
{return Evaluation< Expression2<GeneralForm<Form>>>(form,shapes);}



// template< typename Derived>
// class Evaluation< UnaryPlus2< Expression2<Derived> > >
// {   
//  public:
//  using type= typename Evaluation< Derived >::type;

//  Evaluation(const UnaryPlus2< Expression2<Derived> >& expr):
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
// class Evaluation< UnaryMinus2< Expression2<Derived> > >
// {   
//  public:
//  using type= typename Evaluation< Derived >::type;

//  Evaluation(const UnaryMinus2< Expression2<Derived> >& expr):
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
// class Evaluation< Addition2< Expression2<DerivedLeft> , Expression2<DerivedRight> > >
// {   
//  public:
//  using Left= typename Evaluation< DerivedLeft >::type;
//  using Right=typename Evaluation< DerivedRight >::type;
//  using type= typename OperatorType< Addition<  Left, Right > >::type; // SBBAGLIATO, il prodotto

//  // typename OperatorType< Addition< typename Left::type, typename Right::type > >::type; // SBBAGLIATO, il prodotto

//  Evaluation(const Addition2< Expression2<DerivedLeft> , Expression2<DerivedRight> >& expr):
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
// class Evaluation< Subtraction2< Expression2<DerivedLeft> , Expression2<DerivedRight> > >
// {   
//  public:
//  using Left= typename Evaluation< DerivedLeft >::type;
//  using Right=typename Evaluation< DerivedRight >::type;
//  using type= typename OperatorType< Subtraction<  Left, Right > >::type;
//  // typename OperatorType< Subtraction< typename Left::type, typename Right::type > >::type; // SBBAGLIATO, il prodotto

//  Evaluation(const Subtraction2< Expression2<DerivedLeft> , Expression2<DerivedRight> >& expr):
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
// class Evaluation< Multiplication2< Expression2<DerivedLeft> , Expression2<DerivedRight> > >
// {   
//  public:
//  using Left= typename Evaluation< DerivedLeft >::type;
//  using Right=typename Evaluation< DerivedRight >::type;
//  using type= typename OperatorType< Multiply<  Left, Right > >::type;
//  // typename OperatorType< Multiply< typename Left::type, typename Right::type > >::type; // SBBAGLIATO, il prodotto

//  Evaluation(const Multiplication2< Expression2<DerivedLeft> , Expression2<DerivedRight> >& expr):
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
// class Evaluation< Multiplication2< Expression2<DerivedLeft> , Real > >
// {   
//  public:
//  using Left= typename Evaluation< DerivedLeft >::type;
//  using Right=Real;
//  using type= typename OperatorType< Multiply<  Left, Right > >::type;
//  // typename OperatorType< Multiply< typename Left::type, Real> >::type; // SBBAGLIATO, il prodotto

//  Evaluation(const Multiplication2< Expression2<DerivedLeft> , Real >& expr):
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
// class Evaluation< Multiplication2< Real, Expression2<DerivedRight> > >
// {   
//  public:
//  using Left= Real;
//  using Right=typename Evaluation< DerivedRight >::type;
//  using type= typename OperatorType< Multiply<  Left, Right > >::type;
//  // typename OperatorType< Multiply< Real, typename Right::type> >::type; 

//  Evaluation(const Multiplication2< Real, Expression2<DerivedRight> >& expr):
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
// class Evaluation< Division2< Expression2<DerivedLeft> , Real > >
// {   
//  public:
//  using Left= typename Evaluation< DerivedLeft >::type;
//  using Right=Real;
//  using type= typename OperatorType< Divide<  Left, Right > >::type;
//  // typename OperatorType< Divide< typename Left::type, Real> >::type; // SBBAGLIATO, il prodotto

//  Evaluation(const Division2< Expression2<DerivedLeft> , Real >& expr):
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
class Expression2QPValues: 
public Expression2< Expression2QPValues<T,NQPoints > >
{
public:
        using type=QPValues<T,NQPoints>;
        Expression2QPValues(){};

        Expression2QPValues(const type& value):value_(value){};

        type& operator()(){return value_;};
        const type& operator()()const{return value_;};
protected:
    type value_;

};

template<typename T,Integer NQPoints>
class Evaluation<Expression2QPValues<T,NQPoints> >
{
 public:
 using type=QPValues<T,NQPoints>;

 Evaluation(const Expression2QPValues<T,NQPoints>&expr):
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
class Expression2FQPValues: 
public Expression2< Expression2FQPValues<T,NQPoints,NComponents > >
{
public:
        using type=FQPValues<T,NQPoints,NComponents>;
        Expression2FQPValues(){};

        Expression2FQPValues(const type& value):value_(value){};

        type& operator()(){return value_;};
        const type& operator()()const{return value_;};
protected:
    type value_;

};

template<typename T,Integer NQPoints, Integer NComponents>
class Evaluation<Expression2FQPValues<T,NQPoints,NComponents> >
{
 public:
 using type=FQPValues<T,NQPoints,NComponents>;

 Evaluation(const Expression2FQPValues<T,NQPoints,NComponents>&expr):
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
class Expression2MatrixVar: 
public Expression2< Expression2MatrixVar<T,Rows,Cols> >
{
public:
        using type=Matrix<T,Rows,Cols>;
        Expression2MatrixVar(){};

        type& operator()(){return value_;};
        const type& operator()()const{return value_;};
protected:
    type value_;

};

template<typename T,Integer Rows,Integer Cols>
class Evaluation<Expression2MatrixVar<T,Rows,Cols> >
{
 public:
 using type=Matrix<T,Rows,Cols>;

 Evaluation(const Expression2MatrixVar<T,Rows,Cols>&expr):
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
class Expression2Matrix: 
public Expression2< Expression2Matrix<T,Rows,Cols> >
{
public:
        using type=Matrix<T,Rows,Cols>;
        Expression2Matrix(){};

        Expression2Matrix(const type& value):value_(value){};

        type& operator()(){return value_;};
        const type& operator()()const{return value_;};
protected:
    type value_;

};

template<typename T,Integer Rows,Integer Cols>
class Evaluation<Expression2Matrix<T,Rows,Cols> >
{
 public:
 using type=Matrix<T,Rows,Cols>;

 Evaluation(const Expression2Matrix<T,Rows,Cols>&expr):
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
class QuadratureOrder< UnaryPlus2< Expression2<T> > >
{ public:
  static constexpr Integer value=QuadratureOrder<T>::value;
};


// order(T) = order(-T)
template<typename T>
class QuadratureOrder< UnaryMinus2< Expression2<T> > >
{ public:
  static constexpr Integer value=QuadratureOrder<T>::value;
};




// order(A+B) =max(order(A),order(B))
template<typename Left, typename Right>
class QuadratureOrder< Addition2< Expression2<Left>,Expression2<Right> > >
{ public:
  static constexpr Integer value=Max(QuadratureOrder<Left>::value,QuadratureOrder<Right>::value);
};


// order(A+B) =max(order(A),order(B))
template<typename Left, typename Right>
class QuadratureOrder< Subtraction2< Expression2<Left>,Expression2<Right> > >
{ public:
  static constexpr Integer value=Max(QuadratureOrder<Left>::value,QuadratureOrder<Right>::value);
};


// order(T) = order(T/REAL)
template<typename T>
class QuadratureOrder< Division2<Expression2<T>, Real > >
{ public:
  static constexpr Integer value=QuadratureOrder<T>::value;
};

// order(T) = order(T * REAL)
template<typename T>
class QuadratureOrder< Multiplication2< Expression2<T>, Real > >
{ public:
  static constexpr Integer value=QuadratureOrder<T>::value;
};

// order(T) = order(REAL * T)
template<typename T>
class QuadratureOrder< Multiplication2< Real, Expression2<T> > >
{ public:
  static constexpr Integer value=QuadratureOrder<T>::value;
};


// order(Left*Right) = order(Left) * order(Right)
template<typename Left, typename Right>
class QuadratureOrder< Multiplication2< Expression2<Left>, Expression2<Right> > >
{ public:

  static constexpr Integer value=QuadratureOrder<Left>::value + QuadratureOrder<Right>::value;
};


// order(Left*Right) = order(Left) * order(Right)
template<typename Left, typename Right>
class QuadratureOrder< Contraction2< Expression2<Left>, Expression2<Right> > >
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
class ShapeFunctionExpression: public Expression2<ShapeFunctionExpression<N,Nmax,Ts...>>
{
public:
	static constexpr Integer value=N;

};

template<Integer N,Integer Nmax,typename...Ts>
class GradientShapeFunctionExpression: public Expression2<GradientShapeFunctionExpression<N,Nmax,Ts...>>
{
public:
	static constexpr Integer value=N;

};

template<typename...T>
class OperatorTupleType;


// order(T) = order(+T)
template<typename T>
class OperatorTupleType< UnaryPlus2< Expression2<T> > >
{ public:
  using type=typename OperatorTupleType<T>::type;
};


// order(T) = order(-T)
template<typename T>
class OperatorTupleType< UnaryMinus2< Expression2<T> > >
{ public:
  using type=typename OperatorTupleType<T>::type;
};




// order(A+B) =max(order(A),order(B))
template<typename Left, typename Right>
class OperatorTupleType< Addition2< Expression2<Left>,Expression2<Right> > >
{ public:
  using	LeftT=typename OperatorTupleType<Left>::type;
  using	RightT=typename OperatorTupleType<Right>::type;
  using type=RemoveTupleOfTupleDuplicates< LeftT, RightT >;
};


// order(A+B) =max(order(A),order(B))
template<typename Left, typename Right>
class OperatorTupleType< Subtraction2< Expression2<Left>,Expression2<Right> > >
{ public:
  using	LeftT=typename OperatorTupleType<Left>::type;
  using	RightT=typename OperatorTupleType<Right>::type;
  using type=RemoveTupleOfTupleDuplicates<  LeftT, RightT >;
};


// order(T) = order(T/REAL)
template<typename T>
class OperatorTupleType< Division2<Expression2<T>, Real > >
{ public:
  using type=typename OperatorTupleType<T>::type;
};

// order(T) = order(T * REAL)
template<typename T>
class OperatorTupleType< Multiplication2< Expression2<T>, Real > >
{ public:
  using type=typename OperatorTupleType<T>::type;
};

// order(T) = order(REAL * T)
template<typename T>
class OperatorTupleType< Multiplication2< Real, Expression2<T> > >
{ public:
  using type=typename OperatorTupleType<T>::type;
};


// order(Left*Right) = order(Left) * order(Right)
template<typename Left, typename Right>
class OperatorTupleType< Multiplication2< Expression2<Left>, Expression2<Right> > >
{ public:
  using	LeftT=typename OperatorTupleType<Left>::type;
  using	RightT=typename OperatorTupleType<Right>::type;
  using type=RemoveTupleOfTupleDuplicates<  LeftT, RightT >;
};


// order(Left*Right) = order(Left) * order(Right)
template<typename Left, typename Right>
class OperatorTupleType< Contraction2< Expression2<Left>, Expression2<Right> > >
{ public:
  using LeftT=typename OperatorTupleType<Left>::type;
  using RightT=typename OperatorTupleType<Right>::type;
  using type=RemoveTupleOfTupleDuplicates<  LeftT, RightT >;
};





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
class OperatorTupleType<Expression2QPValues<T...> >
{ public:
  using type=std::tuple<>;
  static constexpr Integer value=-1;
};

}

#endif //MARS_OPERATORS_HPP
