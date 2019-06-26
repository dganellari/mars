#ifndef MARS_OPERATORS_HPP
#define MARS_OPERATORS_HPP
#include "mars_base.hpp"
#include "mars_tuple_utilities.hpp"
namespace mars {

///////////////////////////////////////////////////////////////
/////////////////////    MAIN OPERATORS    ////////////////////
///////////////////////////////////////////////////////////////

template<typename...Parameters>
class UnaryPlus;

template<typename...Parameters>
class UnaryMinus;

template<typename...Parameters>
class Multiply;

template<typename...Parameters>
class Divide;

template<typename...Parameters>
class Subtraction;

template<typename...Parameters>
class Addition;



template<typename Derived, typename S,Integer NQPoints, Integer Dim>
class QPprova;

template<typename Derived, typename S,Integer NQPoints, Integer NComponents, Integer Dim>
class FQPprova;

template<typename T,Integer NQPoints,Integer NComponents>
class FQPValues;

template<typename T,Integer NQPoints>
class QPValues;

template<typename...Parameters>
class GradientExpression;

template<Integer N,typename...Parameters>
class ShapeFunctionExpression;

template<typename Operator>
class OperatorType;


// type(T) = type(+T)
template<typename T>
class OperatorType
{ public:
  using type=T;
};

// type(T) = type(+T)
template<typename T>
class OperatorType< UnaryPlus< T > >
{ public:
  using type=T;
};



// type(T) = type(-T)
template<typename T>
class OperatorType< UnaryMinus< T > >
{ public:
  using type=T;
};


// type(A)=type(B) = type(A+B)
template<typename Left, typename Right>
class OperatorType< Addition< Left,Right > >
{ public:
  static_assert(std::is_same<Left,Right>::value, " In Addition, Left and Right types must be equal");
  using type=Left;
};


// type(A)=type(B) = type(A-B)
template<typename Left, typename Right>
class OperatorType< Subtraction< Left,Right > >
{ public:
  static_assert(std::is_same<Left,Right>::value, " In Subtraction, Left and Right types must be equal");
  using type=Left;
};

// type(T) = type(T/REAL)
template<typename T>
class OperatorType< Divide<T, Real > >
{ public:
  using type=T;
};

// type(T) = type(T * REAL)
template<typename T>
class OperatorType< Multiply< T, Real > >
{ public:
  using type=T;
};

// type(T) = type(REAL * T)
template<typename T>
class OperatorType< Multiply< Real, T > >
{ public:
  using type=T;
};


// type(VECTOR) = type(MATRIX * VECTOR)
template<typename T, Integer RowsLeft,Integer Common>
class OperatorType< Multiply< Matrix<T, RowsLeft, Common>, Vector<T, Common> > >
{ public:
  using type=Vector<T, RowsLeft>;
};

// type(MATRIX) = type(MATRIX * MATRIX)

template<typename T, Integer RowsLeft,Integer Common,Integer ColsRight>
class OperatorType< Multiply< Matrix<T, RowsLeft, Common>, Matrix<T, Common, ColsRight> > >
{ public:
  using type=Matrix<T, RowsLeft, ColsRight>;
};

template<typename T, Integer RowsLeft,Integer Common,Integer ColsRight, Integer NQPoints, Integer NComponents>
class OperatorType< Multiply< QPValues < Matrix<T, RowsLeft, Common>, NQPoints> ,
                                   FQPValues< Matrix<T, Common, ColsRight>, NQPoints, NComponents> > >
{ public:
   using type=FQPValues< Matrix<T, RowsLeft, ColsRight>, NQPoints, NComponents>;
};

template<typename T, Integer RowsLeft,Integer Common,Integer ColsRight, Integer NQPoints>
class OperatorType< Multiply< QPValues< Matrix<T, RowsLeft, Common>, NQPoints> ,
                                   QPValues< Matrix<T, Common, ColsRight>, NQPoints> > >
{ public:
   using type=QPValues< Matrix<T, RowsLeft, ColsRight>, NQPoints>;
};


template< typename DerivedLeft,typename DerivedRight, typename  TLeft,typename TRight, Integer NQPoints,Integer NComponents,Integer Dim>
class OperatorType< Multiply< QPprova < DerivedLeft, TLeft, NQPoints,Dim> ,
                                   FQPprova< DerivedRight,TRight, NQPoints, NComponents,Dim> > >
{ public:
   using type= FQPValues<typename OperatorType<Multiply<TLeft,TRight>>::type , NQPoints, NComponents>;
};

template< typename DerivedLeft,typename DerivedRight, typename  TLeft,typename TRight, Integer NQPoints,Integer Dim>
class OperatorType< Multiply< QPprova< DerivedLeft, TLeft,  NQPoints,Dim> ,
                                   QPprova< DerivedRight,TRight, NQPoints,Dim> > >
{ public:
   using type= QPValues<typename OperatorType<Multiply<TLeft,TRight>>::type , NQPoints>;
};


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
       Matrix<T,1,Rows> result=0;
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



template<typename...Parameters>
class Evaluation;

template<typename...Parameters>
class UnaryPlus2;

template<typename...Parameters>
class UnaryMinus2;

template<typename...Parameters>
class Multiply2;

template<typename...Parameters>
class Divide2;

template<typename...Parameters>
class Subtraction2;

template<typename...Parameters>
class Addition2;



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
class UnaryPlus<Matrix<T,Rows,Cols>,
                Matrix<T,Rows,Cols> > 
{    
 public:       
 inline static void apply(const Matrix<T,Rows,Cols>& A,
                                Matrix<T,Rows,Cols>& C)
	{
	 for(Integer ii=0;ii<Rows;ii++)
	 	for(Integer jj=0;jj<Cols;jj++)
	 	   	   C(ii,jj)=A(ii,jj);
	};
};

template<typename T,Integer Rows,Integer Cols>
class UnaryMinus<Matrix<T,Rows,Cols>,
                 Matrix<T,Rows,Cols> > 
{    
 public:       
 inline static void apply(const Matrix<T,Rows,Cols>& A,
                                Matrix<T,Rows,Cols>& C)
	{
	 for(Integer ii=0;ii<Rows;ii++)
	 	for(Integer jj=0;jj<Cols;jj++)
	 	   	   C(ii,jj)=-A(ii,jj);
	};
};

template<typename T,Integer Rows,Integer Cols>
class Addition<Matrix<T,Rows,Cols>,
               Matrix<T,Rows,Cols>,
               Matrix<T,Rows,Cols> > 
{    
 public:       
 inline static void apply(const Matrix<T,Rows,Cols>& A,
 	                      const Matrix<T,Rows,Cols>& B, 
                                Matrix<T,Rows,Cols>& C)
	{
	 for(Integer ii=0;ii<Rows;ii++)
	 	for(Integer jj=0;jj<Cols;jj++)
	 	   	   C(ii,jj)=A(ii,jj)+B(ii,jj);
	};
};

template<typename T,Integer Rows,Integer Cols>
class Subtraction<Matrix<T,Rows,Cols>,
               Matrix<T,Rows,Cols>,
               Matrix<T,Rows,Cols> > 
{    
 public:       
 inline static void apply(const Matrix<T,Rows,Cols>& A,
 	                      const Matrix<T,Rows,Cols>& B, 
                                Matrix<T,Rows,Cols>& C)
	{
	 for(Integer ii=0;ii<Rows;ii++)
	 	for(Integer jj=0;jj<Cols;jj++)
	 	   	   C(ii,jj)=A(ii,jj)-B(ii,jj);
	};
};


template<typename T,Integer Rows,Integer CommonDim,Integer Cols>
class Multiply<Matrix<T,Rows,CommonDim>,
               Matrix<T,CommonDim,Cols>,
               Matrix<T,Rows,Cols> > 
{    
 public:            
 inline static void apply(const Matrix<T,Rows,CommonDim>& A,
 	                      const Matrix<T,CommonDim,Cols>& B, 
                                Matrix<T,Rows,Cols>& C)
	{
	 for(Integer ii=0;ii<Rows;ii++)
	 	for(Integer jj=0;jj<Cols;jj++)
	 	   {
	 	   	C(ii,jj)=0;
	 	   	for(Integer cc=0;cc<CommonDim;cc++)
	 	   	   C(ii,jj)+=A(ii,cc)*B(cc,jj);
	 	   }
	};
};

template<typename T,Integer Rows,Integer Cols>
class Multiply<Matrix<T,Rows,Cols>,
               Real,
               Matrix<T,Rows,Cols> > 
{    
 public:            
 inline static void apply(const Matrix<T,Rows,Cols>& A,
 	                      const Real& B, 
                                Matrix<T,Rows,Cols>& C)
	{
	 for(Integer ii=0;ii<Rows;ii++)
	 	for(Integer jj=0;jj<Cols;jj++)
	 	   	   C(ii,jj)=A(ii,jj)*B;   
	};
};

template<typename T,Integer Rows,Integer Cols>
class Multiply<Real,
               Matrix<T,Rows,Cols>,
               Matrix<T,Rows,Cols> > 
{    
 public:            
 inline static void apply(const Real& A,
 	                      const Matrix<T,Rows,Cols>& B,
                                Matrix<T,Rows,Cols>& C)
	{
	 for(Integer ii=0;ii<Rows;ii++)
	 	for(Integer jj=0;jj<Cols;jj++)
	 	   	   C(ii,jj)=B(ii,jj)*A;   
	};
};


template<typename T,Integer Rows,Integer Cols>
class Divide<Matrix<T,Rows,Cols>,
             Real,
             Matrix<T,Rows,Cols> > 
{    
 public:            
 inline static void apply(const Matrix<T,Rows,Cols>& A,
 	                      const Real& B, 
                                Matrix<T,Rows,Cols>& C)
	{
	 for(Integer ii=0;ii<Rows;ii++)
	 	for(Integer jj=0;jj<Cols;jj++)
	 	   	   C(ii,jj)=A(ii,jj)/B;   
	};
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////			    QP   EXPRESSION ALGEBRAIC OPERATIONS		    ///////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename U,Integer NQPoints>
class UnaryPlus<QPValues<U,NQPoints>,
                QPValues<U,NQPoints>> 
{       
 public:               
 inline static void apply(const QPValues<U,NQPoints>& A,
                                QPValues<U,NQPoints>& C)
	{
	 for(Integer qp=0;qp<NQPoints;qp++)
        UnaryPlus<U,U>::apply(A[qp],C[qp]);
	};

};

template<typename U,Integer NQPoints>
class UnaryMinus<QPValues<U,NQPoints>,
                 QPValues<U,NQPoints>> 
{       
 public:               
 inline static void apply(const QPValues<U,NQPoints>& A,
                                QPValues<U,NQPoints>& C)
	{
	 for(Integer qp=0;qp<NQPoints;qp++)
        UnaryMinus<U,U>::apply(A[qp],C[qp]);
	};

};

template<typename U,Integer NQPoints>
class Addition<QPValues<U,NQPoints>,
               QPValues<U,NQPoints>,
               QPValues<U,NQPoints> > 
{       
 public:               
 inline static void apply(const QPValues<U,NQPoints>& A,
 	                      const QPValues<U,NQPoints>& B, 
                                QPValues<U,NQPoints>& C)
	{
	 for(Integer qp=0;qp<NQPoints;qp++)
        Addition<U,U,U>::apply(A[qp],B[qp],C[qp]);
	};

};

template<typename U,Integer NQPoints>
class Subtraction<QPValues<U,NQPoints>,
                  QPValues<U,NQPoints>,
                  QPValues<U,NQPoints> > 
{       
 public:               
 inline static void apply(const QPValues<U,NQPoints>& A,
 	                      const QPValues<U,NQPoints>& B, 
                                QPValues<U,NQPoints>& C)
	{
	 for(Integer qp=0;qp<NQPoints;qp++)
        Subtraction<U,U,U>::apply(A[qp],B[qp],C[qp]);
	};

};

template<typename U, Integer NQPoints>
class Multiply<QPValues<U,NQPoints>,
               Real,
               QPValues<U,NQPoints> > 
{       
 public:               

 inline static void apply(const QPValues<U,NQPoints>& A,
 	                      const Real& B, 
                                QPValues<U,NQPoints>& C)
	{
	 for(Integer qp=0;qp<NQPoints;qp++)
        Multiply<U,Real,U>::apply(A[qp],B,C[qp]);
	};

};

template<typename U,Integer NQPoints>
class Multiply<Real,
               QPValues<U,NQPoints>,
               QPValues<U,NQPoints> > 
{       
 public:               

 inline static void apply(const Real& A,
 	                      const QPValues<U,NQPoints>& B, 
                                QPValues<U,NQPoints>& C)
	{
	 for(Integer qp=0;qp<NQPoints;qp++)
        Multiply<Real,U,U>::apply(A,B[qp],C[qp]);
	};

};

template<typename U,typename V, Integer NQPoints>
class Multiply<QPValues<U,NQPoints>,
               QPValues<V,NQPoints>,
               QPValues<typename OperatorType<Multiply<U,V>>::type,NQPoints> > 
{       
 public:               
 using W=typename OperatorType<Multiply<U,V>>::type;

 inline static void apply(const QPValues<U,NQPoints>& A,
 	                      const QPValues<V,NQPoints>& B, 
                                QPValues<W,NQPoints>& C)
	{
	 for(Integer qp=0;qp<NQPoints;qp++)
        Multiply<U,V,W>::apply(A[qp],B[qp],C[qp]);
	};

};

template<typename U, Integer NQPoints>
class Divide<QPValues<U,NQPoints>,
               Real,
               QPValues<U,NQPoints> > 
{       
 public:               

 inline static void apply(const QPValues<U,NQPoints>& A,
 	                      const Real& B, 
                                QPValues<U,NQPoints>& C)
	{
	 for(Integer qp=0;qp<NQPoints;qp++)
        Divide<U,Real,U>::apply(A[qp],B,C[qp]);
	};

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////			    FQP   EXPRESSION ALGEBRAIC OPERATIONS		    ///////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename U, Integer NQPoints,Integer NComponents>
class UnaryPlus< FQPValues<U,NQPoints,NComponents>,
                 FQPValues<U,NQPoints,NComponents> > 
{       
 public:               
 inline static void apply(const FQPValues<U,NQPoints,NComponents>& A,
                                FQPValues<U,NQPoints,NComponents>& C)
	{
	 for(Integer qp=0;qp<NQPoints;qp++)
	 	for(Integer n_comp=0;n_comp<NComponents;n_comp++)
	 	{
	 		UnaryPlus<U,U>::apply(A[n_comp][qp],C[n_comp][qp]);
	 	}
	};

};

template<typename U, Integer NQPoints,Integer NComponents>
class UnaryMinus< FQPValues<U,NQPoints,NComponents>,
                 FQPValues<U,NQPoints,NComponents> > 
{       
 public:               
 inline static void apply(const FQPValues<U,NQPoints,NComponents>& A,
                                FQPValues<U,NQPoints,NComponents>& C)
	{
	 for(Integer qp=0;qp<NQPoints;qp++)
	 	for(Integer n_comp=0;n_comp<NComponents;n_comp++)
	 	{
	 		UnaryMinus<U,U>::apply(A[n_comp][qp],C[n_comp][qp]);
	 	}
	};

};

template<typename U, Integer NQPoints,Integer NComponents>
class Addition< FQPValues<U,NQPoints,NComponents>,
                FQPValues<U,NQPoints,NComponents>,
                FQPValues<U,NQPoints,NComponents> > 
{       
 public:               
 inline static void apply(const FQPValues<U,NQPoints,NComponents>& A,
 	                      const FQPValues<U,NQPoints,NComponents>& B, 
                                FQPValues<U,NQPoints,NComponents>& C)
	{
	 for(Integer qp=0;qp<NQPoints;qp++)
	 	for(Integer n_comp=0;n_comp<NComponents;n_comp++)
	 	{
	 		Addition<U,U,U>::apply(A[n_comp][qp],B[n_comp][qp],C[n_comp][qp]);
	 	}
	};

};


template<typename U, Integer NQPoints,Integer NComponents>
class Subtraction< FQPValues<U,NQPoints,NComponents>,
                   FQPValues<U,NQPoints,NComponents>,
                   FQPValues<U,NQPoints,NComponents> > 
{       
 public:               
 inline static void apply(const FQPValues<U,NQPoints,NComponents>& A,
 	                      const FQPValues<U,NQPoints,NComponents>& B, 
                                FQPValues<U,NQPoints,NComponents>& C)
	{
	 for(Integer qp=0;qp<NQPoints;qp++)
	 	for(Integer n_comp=0;n_comp<NComponents;n_comp++)
	 	{
	 		Subtraction<U,U,U>::apply(A[n_comp][qp],B[n_comp][qp],C[n_comp][qp]);
	 	}
	};

};



template<typename U,Integer NQPoints,Integer NComponents>
class Multiply< Real,
                FQPValues<U,NQPoints,NComponents>,
                FQPValues<U,NQPoints,NComponents> > 
{       
 public:               

 inline static void apply(const Real& A,
 	                      const FQPValues<U,NQPoints,NComponents>& B, 
                                FQPValues<U,NQPoints,NComponents>& C)
	{
	 for(Integer qp=0;qp<NQPoints;qp++)
	 	for(Integer n_comp=0;n_comp<NComponents;n_comp++)
	 		Multiply<Real,U,U>::apply(A,B[n_comp][qp],C[n_comp][qp]);
	};

};


template<typename U, Integer NQPoints,Integer NComponents>
class Multiply< FQPValues<U,NQPoints,NComponents>,
                Real,
                FQPValues<U,NQPoints,NComponents> > 
{       
 public:               

 inline static void apply(const FQPValues <U,NQPoints,NComponents>& A,
 	                      const Real& B, 
                                FQPValues<U,NQPoints,NComponents>& C)
	{
	 for(Integer qp=0;qp<NQPoints;qp++)
	 	for(Integer n_comp=0;n_comp<NComponents;n_comp++)
	 	{
	 		Multiply<U,Real,U>::apply(A[n_comp][qp],B,C[n_comp][qp]);
	 	}
	};

};


template<typename U, Integer NQPoints,Integer NComponents>
class Divide< FQPValues<U,NQPoints,NComponents>,
              Real,
              FQPValues<U,NQPoints,NComponents> > 
{       
 public:               

 inline static void apply(const FQPValues <U,NQPoints,NComponents>& A,
 	                      const Real& B, 
                                FQPValues<U,NQPoints,NComponents>& C)
	{
	 for(Integer qp=0;qp<NQPoints;qp++)
	 	for(Integer n_comp=0;n_comp<NComponents;n_comp++)
	 	{
	 		Divide<U,Real,U>::apply(A[n_comp][qp],B,C[n_comp][qp]);
	 	}
	};

};




///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////				    QP-FQP    ALGEBRAIC OPERATIONS			    ///////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////



template<typename U,typename V, Integer NQPoints,Integer NComponents>
class Multiply< QPValues<U,NQPoints>,
                FQPValues<V,NQPoints,NComponents>,
                FQPValues<typename OperatorType<Multiply<U,V>>::type,NQPoints,NComponents> > 
{       
 public:               
 using W=typename OperatorType<Multiply<U,V>>::type;

 inline static void apply(const QPValues <U,NQPoints>& A,
 	                      const FQPValues<V,NQPoints,NComponents>& B, 
                                FQPValues<W,NQPoints,NComponents>& C)
	{
	 for(Integer qp=0;qp<NQPoints;qp++)
	 	for(Integer n_comp=0;n_comp<NComponents;n_comp++)
	 	{
	 		Multiply<U,V,W>::apply(A[qp],B[n_comp][qp],C[n_comp][qp]);
	 	}
	 std::cout<<"qp * fqp C=="<<C<<std::endl;
	};

};

template<typename U,typename V, Integer NQPoints,Integer NComponents>
class Multiply< FQPValues<U,NQPoints,NComponents>,
                QPValues<V,NQPoints>,
                FQPValues<typename OperatorType<Multiply<U,V>>::type,NQPoints,NComponents> > 
{       
 public:               
 using W=typename OperatorType<Multiply<U,V>>::type;

 inline static void apply(const FQPValues <U,NQPoints,NComponents>& A,
 	                      const QPValues<V,NQPoints>& B, 
                                FQPValues<W,NQPoints,NComponents>& C)
	{
	 for(Integer qp=0;qp<NQPoints;qp++)
	 	for(Integer n_comp=0;n_comp<NComponents;n_comp++)
	 	{
	 		Multiply<U,V,W>::apply(A[qp],B[n_comp][qp],C[n_comp][qp]);
	 	}
	};

};























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
class Divide2< Expression2 <DerivedLeft>, Expression2 <DerivedRight> > 
: public Expression2< Divide2< Expression2 <DerivedLeft>, Expression2 <DerivedRight> > >
{
  public:
    Divide2(const Expression2<DerivedLeft>& left, const Expression2<DerivedRight>&right)
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
class Divide2< Expression2 <DerivedLeft>, Real > 
: public Expression2< Divide2< Expression2 <DerivedLeft>, Real > >
{
  public:
    Divide2(const Expression2<DerivedLeft>& left, const Real&right)
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
class Multiply2< Expression2 <DerivedLeft_>, Expression2 <DerivedRight_> > 
: public Expression2< Multiply2< Expression2 <DerivedLeft_>, Expression2 <DerivedRight_> > >
{
  public:
  	using DerivedLeft=DerivedLeft_;
  	using DerivedRight=DerivedRight_;

    Multiply2(const Expression2<DerivedLeft>& left, const Expression2<DerivedRight>&right)
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
class Multiply2< Expression2 <DerivedLeft_>, Real > 
: public Expression2< Multiply2< Expression2 <DerivedLeft_>, Real > >
{
  public:
  	using DerivedLeft=DerivedLeft_;
  	using DerivedRight=Real;
    Multiply2(const Expression2<DerivedLeft>& left, const Real&right)
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
class Multiply2< Real,Expression2 <DerivedRight_> > 
: public Expression2< Multiply2<Real, Expression2 <DerivedRight_> > >
{
  public:
  	using DerivedLeft=Real;
  	using DerivedRight=DerivedRight_;
    Multiply2(const Real& left, const Expression2 <DerivedRight>&right)
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


template< typename DerivedLeft>
class Divide2< Expression2 <DerivedLeft>, Real > 
operator/(const Expression2<DerivedLeft>&left, const Real&right)
{return Divide2< Expression2 <DerivedLeft>, Real > (left,right);}


template< typename DerivedLeft,typename DerivedRight>
class Multiply2< Expression2 <DerivedLeft>, Expression2 <DerivedRight> >
operator*(const Expression2<DerivedLeft>&left, const Expression2<DerivedRight>&right)
{return Multiply2< Expression2 <DerivedLeft>, Expression2 <DerivedRight> >(left,right);}

template< typename DerivedLeft>
class Multiply2< Expression2 <DerivedLeft>, Real > 
operator*(const Expression2<DerivedLeft>&left, const Real&right)
{return Multiply2< Expression2 <DerivedLeft>, Real > (left,right);}

template< typename DerivedRight>
class Multiply2< Real, Expression2 <DerivedRight> >
operator*(const Real&left, const Expression2<DerivedRight>&right)
{return Multiply2< Real, Expression2 <DerivedRight> >(left,right);}



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





// template<typename Unary2, typename Derived>
// class Evaluation< Unary2< Expression2<Derived> > >
// {   
//  public:
//  using type= typename Evaluation< Derived >::type;

//  Evaluation(const Unary2< Expression2<Derived> >& expr):
//  eval_(expr())
//  {};

//  type& operator()(){return self_;};
//  const type& operator()()const{return self_;};

//  template<typename...Parameters>
//  inline type& apply(const Parameters&...parameters)
//  {
//   self_ =eval_.apply(parameters...);
//   Unary2<type,type>::apply(self_,value_);
//   return value_;
//  };
// private:
//  type value_;
//  type self_;
//  Evaluation< Derived > eval_;
// };



template< typename Derived>
class Evaluation< UnaryPlus2< Expression2<Derived> > >
{   
 public:
 using type= typename Evaluation< Derived >::type;

 Evaluation(const UnaryPlus2< Expression2<Derived> >& expr):
 eval_(expr())
 {};

 type& operator()(){return self_;};
 const type& operator()()const{return self_;};

 template<typename...Parameters>
 inline type& apply(const Parameters&...parameters)
 {
  self_ =eval_.apply(parameters...);
  UnaryPlus<type,type>::apply(self_,value_);
  return value_;
 };
private:
 type value_;
 type self_;
 Evaluation< Derived > eval_;
};


template< typename Derived>
class Evaluation< UnaryMinus2< Expression2<Derived> > >
{   
 public:
 using type= typename Evaluation< Derived >::type;

 Evaluation(const UnaryMinus2< Expression2<Derived> >& expr):
 eval_(expr())
 {};

 type& operator()(){return self_;};
 const type& operator()()const{return self_;};

 template<typename...Parameters>
 inline type& apply(const Parameters&...parameters)
 {
  self_ =eval_.apply(parameters...);
  UnaryMinus<type,type>::apply(self_,value_);
  return value_;
 };
private:
 type value_;
 type self_;
 Evaluation< Derived > eval_;
};


template< typename DerivedLeft,typename DerivedRight>
class Evaluation< Addition2< Expression2<DerivedLeft> , Expression2<DerivedRight> > >
{   
 public:
 using Left= typename Evaluation< DerivedLeft >::type;
 using Right=typename Evaluation< DerivedRight >::type;
 using type= typename OperatorType< Addition<  Left, Right > >::type; // SBBAGLIATO, il prodotto

 // typename OperatorType< Addition< typename Left::type, typename Right::type > >::type; // SBBAGLIATO, il prodotto

 Evaluation(const Addition2< Expression2<DerivedLeft> , Expression2<DerivedRight> >& expr):
 eval_left_(expr.left()),
 eval_right_(expr.right())
 {};

 Left& left(){return left_;};
 Right& right(){return right_;};
 const Left& left()const{return left_;};
 const Right& right()const{return right_;};

 template<typename...Parameters>
 inline type& apply(const Parameters&...parameters)
 {
  left_ =eval_left_.apply(parameters...);
  right_=eval_right_.apply(parameters...);
  Addition<Left,Right,type>::apply(left_,right_,value_);
  return value_;
 };
private:
 Left left_;
 Right right_;
 type value_;
 Evaluation< DerivedLeft > eval_left_;
 Evaluation< DerivedRight> eval_right_; 
};



template< typename DerivedLeft,typename DerivedRight>
class Evaluation< Subtraction2< Expression2<DerivedLeft> , Expression2<DerivedRight> > >
{   
 public:
 using Left= typename Evaluation< DerivedLeft >::type;
 using Right=typename Evaluation< DerivedRight >::type;
 using type= typename OperatorType< Subtraction<  Left, Right > >::type;
 // typename OperatorType< Subtraction< typename Left::type, typename Right::type > >::type; // SBBAGLIATO, il prodotto

 Evaluation(const Subtraction2< Expression2<DerivedLeft> , Expression2<DerivedRight> >& expr):
 eval_left_(expr.left()),
 eval_right_(expr.right())
 {};

 Left& left(){return left_;};
 Right& right(){return right_;};
 const Left& left()const{return left_;};
 const Right& right()const{return right_;};

 template<typename...Parameters>
 inline type& apply(const Parameters&...parameters)
 {
  left_ =eval_left_.apply(parameters...);
  right_=eval_right_.apply(parameters...);
  Subtraction<Left,Right,type>::apply(left_,right_,value_);
  return value_;
 };
private:
 Left left_;
 Right right_;
 type value_;
 Evaluation< DerivedLeft > eval_left_;
 Evaluation< DerivedRight> eval_right_; 
};


template< typename DerivedLeft,typename DerivedRight>
class Evaluation< Multiply2< Expression2<DerivedLeft> , Expression2<DerivedRight> > >
{   
 public:
 using Left= typename Evaluation< DerivedLeft >::type;
 using Right=typename Evaluation< DerivedRight >::type;
 using type= typename OperatorType< Multiply<  Left, Right > >::type;
 // typename OperatorType< Multiply< typename Left::type, typename Right::type > >::type; // SBBAGLIATO, il prodotto

 Evaluation(const Multiply2< Expression2<DerivedLeft> , Expression2<DerivedRight> >& expr):
 eval_left_(expr.left()),
 eval_right_(expr.right())
 {};

 Left& left(){return left_;};
 Right& right(){return right_;};
 const Left& left()const{return left_;};
 const Right& right()const{return right_;};

 template<typename...Parameters>
 inline type& apply(const Parameters&...parameters)
 {
  left_ =eval_left_.apply(parameters...);
  right_=eval_right_.apply(parameters...);
  Multiply<Left,Right,type>::apply(left_,right_,value_);
  return value_;
 };
private:
 Left left_;
 Right right_;
 type value_;
 Evaluation< DerivedLeft > eval_left_;
 Evaluation< DerivedRight> eval_right_; 
};




template< typename DerivedLeft>
class Evaluation< Multiply2< Expression2<DerivedLeft> , Real > >
{   
 public:
 using Left= typename Evaluation< DerivedLeft >::type;
 using Right=Real;
 using type= typename OperatorType< Multiply<  Left, Right > >::type;
 // typename OperatorType< Multiply< typename Left::type, Real> >::type; // SBBAGLIATO, il prodotto

 Evaluation(const Multiply2< Expression2<DerivedLeft> , Real >& expr):
 eval_left_(expr.left()),
 right_(expr.right())
 {};

 Left& left(){return left_;};
 Real& right(){return right_;};
 const Left& left()const{return left_;};
 const Real& right()const{return right_;};

 template<typename...Parameters>
 inline type& apply(const Parameters&...parameters)
 {
  left_ =eval_left_.apply(parameters...);
  Multiply<Left,Right,type>::apply(left_,right_,value_);
  return value_;
 };
private:
 Left left_;
 Real right_;
 type value_;
 Evaluation< DerivedLeft > eval_left_;
};


template< typename DerivedRight>
class Evaluation< Multiply2< Real, Expression2<DerivedRight> > >
{   
 public:
 using Left= Real;
 using Right=typename Evaluation< DerivedRight >::type;
 using type= typename OperatorType< Multiply<  Left, Right > >::type;
 // typename OperatorType< Multiply< Real, typename Right::type> >::type; 

 Evaluation(const Multiply2< Real, Expression2<DerivedRight> >& expr):
 left_(expr.left()),
 eval_right_(expr.right())
 {};

 Left& left(){return left_;};
 Real& right(){return right_;};
 const Left& left()const{return left_;};
 const Real& right()const{return right_;};

 template<typename...Parameters>
 inline type& apply(const Parameters&...parameters)
 {
  right_ =eval_right_.apply(parameters...);
  Multiply<Left,Right,type>::apply(left_,right_,value_);
  return value_;
 };
private:
 Real left_;
 Right right_;
 type value_;
 Evaluation< DerivedRight > eval_right_;
};

template< typename DerivedLeft>
class Evaluation< Divide2< Expression2<DerivedLeft> , Real > >
{   
 public:
 using Left= typename Evaluation< DerivedLeft >::type;
 using Right=Real;
 using type= typename OperatorType< Divide<  Left, Right > >::type;
 // typename OperatorType< Divide< typename Left::type, Real> >::type; // SBBAGLIATO, il prodotto

 Evaluation(const Divide2< Expression2<DerivedLeft> , Real >& expr):
 eval_left_(expr.left()),
 right_(expr.right())
 {};

 Left& left(){return left_;};
 Real& right(){return right_;};
 const Left& left()const{return left_;};
 const Real& right()const{return right_;};

 template<typename...Parameters>
 inline type& apply(const Parameters&...parameters)
 {
  left_ =eval_left_.apply(parameters...);
  Divide<Left,Right,type>::apply(left_,right_,value_);
  return value_;
 };
private:
 Left left_;
 Real right_;
 type value_;
 Evaluation< DerivedLeft > eval_left_;
};







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

// template<typename T,Integer Rows,Integer Cols>
// class Evaluation<Expression2<Expression2Matrix<T,Rows,Cols>> >
// {
//  public:
//  using type=Matrix<T,Rows,Cols>;
//  template<typename...Inputs>
//  inline static void apply(const Expression2<Expression2Matrix<T,Rows,Cols>>&expr, type& value)
//  {
//    value= expr.derived()();
//  };

// };














































































template<typename...T>
class QuadratureOrder;


// order(T) = order(+T)
template<typename T>
class QuadratureOrder< UnaryPlus< Expression2<T> > >
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
  // static_assert(std::is_same<Left,Right>::value, " In Addition, Left and Right types must be equal");
  static constexpr Integer value=Max(QuadratureOrder<Left>::value,QuadratureOrder<Right>::value);
};


// order(A+B) =max(order(A),order(B))
template<typename Left, typename Right>
class QuadratureOrder< Subtraction2< Expression2<Left>,Expression2<Right> > >
{ public:
  // static_assert(std::is_same<Left,Right>::value, " In Addition, Left and Right types must be equal");
  static constexpr Integer value=Max(QuadratureOrder<Left>::value,QuadratureOrder<Right>::value);
};


// order(T) = order(T/REAL)
template<typename T>
class QuadratureOrder< Divide2<Expression2<T>, Real > >
{ public:
  static constexpr Integer value=QuadratureOrder<T>::value;
};

// order(T) = order(T * REAL)
template<typename T>
class QuadratureOrder< Multiply2< Expression2<T>, Real > >
{ public:
  static constexpr Integer value=QuadratureOrder<T>::value;
};

// order(T) = order(REAL * T)
template<typename T>
class QuadratureOrder< Multiply2< Real, Expression2<T> > >
{ public:
  static constexpr Integer value=QuadratureOrder<T>::value;
};


// order(Left*Right) = order(Left) * order(Right)
template<typename Left, typename Right>
class QuadratureOrder< Multiply2< Expression2<Left>, Expression2<Right> > >
{ public:

  static constexpr Integer value=QuadratureOrder<Left>::value + QuadratureOrder<Right>::value;
};





template<typename T,Integer Rows,Integer Cols>
class QuadratureOrder<Expression2MatrixVar<T,Rows,Cols> >
{ public:
  static constexpr Integer value=2;
};














template<Integer N,typename...Ts>
class ShapeFunctionExpression: public Expression2<ShapeFunctionExpression<N,Ts...>>
{
public:
	static constexpr Integer value=N;

};

template<typename...T>
class OperatorTupleType;


// order(T) = order(+T)
template<typename T>
class OperatorTupleType< UnaryPlus< Expression2<T> > >
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
  using type=RemoveDuplicates< TupleCatType< LeftT, RightT> >;
};


// order(A+B) =max(order(A),order(B))
template<typename Left, typename Right>
class OperatorTupleType< Subtraction2< Expression2<Left>,Expression2<Right> > >
{ public:
  using	LeftT=typename OperatorTupleType<Left>::type;
  using	RightT=typename OperatorTupleType<Right>::type;
  using type=RemoveDuplicates< TupleCatType< LeftT, RightT> >;
};


// order(T) = order(T/REAL)
template<typename T>
class OperatorTupleType< Divide2<Expression2<T>, Real > >
{ public:
  using type=typename OperatorTupleType<T>::type;
};

// order(T) = order(T * REAL)
template<typename T>
class OperatorTupleType< Multiply2< Expression2<T>, Real > >
{ public:
  using type=typename OperatorTupleType<T>::type;
};

// order(T) = order(REAL * T)
template<typename T>
class OperatorTupleType< Multiply2< Real, Expression2<T> > >
{ public:
  using type=typename OperatorTupleType<T>::type;
};


// order(Left*Right) = order(Left) * order(Right)
template<typename Left, typename Right>
class OperatorTupleType< Multiply2< Expression2<Left>, Expression2<Right> > >
{ public:
  using	LeftT=typename OperatorTupleType<Left>::type;
  using	RightT=typename OperatorTupleType<Right>::type;
  using type=RemoveDuplicates< TupleCatType< LeftT, RightT> >;
};




template<Integer N,typename...T>
class OperatorTupleType<ShapeFunctionExpression<N,T...> >
{ public:
  using type=std::tuple<IdentityOperator>;
  static constexpr Integer value=ShapeFunctionExpression<N,T...>::value;
};

template<typename...T>
class OperatorTupleType<Expression2QPValues<T...> >
{ public:
  using type=std::tuple<>;
  static constexpr Integer value=-1;
};

}

#endif //MARS_OPERATORS_HPP
