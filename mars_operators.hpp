#ifndef MARS_OPERATORS_HPP
#define MARS_OPERATORS_HPP
#include "mars_base.hpp"
#include "mars_tuple_utilities.hpp"
#include "mars_quadrature_rules.hpp"
#include "mars_expression.hpp"
#include "mars_operator_type.hpp"
#include "mars_operator_tuple_type.hpp"
#include "mars_quadrature_order.hpp"

namespace mars {



template<typename...Ts>
class CompositeOperator;

template<typename ConstType,typename...Inputs>
class ConstantTensor;

template<typename FullSpace,Integer N,typename Operator_,typename FuncType>
class Function;


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

template<typename...Parameters>
class Assignment;

template<typename...Parameters>
class Transposition;

template<typename T,Integer NQPoints,Integer NComponents>
class FQPValues;

template<typename T,Integer NQPoints>
class QPValues;

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

template<>
class Assignment<Real> 
{
 public:       
 inline static constexpr void apply(Real& A,const Real& B)
  {
           A=B;
  };
};

template<typename S,Integer Dim>
class Assignment<Vector<S,Dim>> 
{
 public:       
 template<typename T>
 inline static constexpr void apply(Vector<T,Dim>& A,const Vector<S,Dim>& B)
  {
   for(Integer ii=0;ii<Dim;ii++)
           Assignment<T>::apply(A[ii],B[ii]);
  };
};

template<typename T,Integer Rows,Integer Cols>
class Assignment<Matrix<T,Rows,Cols>> 
{
 public:       
 
 inline static constexpr void apply(Matrix<T,Rows,Cols>& A,const Matrix<T,Rows,Cols>& B)
  {
   for(Integer ii=0;ii<Rows;ii++)
    for(Integer jj=0;jj<Cols;jj++)
      Assignment<T>::apply(A(ii,jj),B(ii,jj));
  };

 inline static constexpr void apply(Transposed<Matrix<T,Rows,Cols>>& A,const Matrix<T,Rows,Cols>& B)
  {
  // update pointer to the non-transposed matrix of the transposed matrix A with B itself
    // std::cout<<"Assignment<Transposed<Matrix<T,Rows,Cols>>> "<<B<<std::endl;

       A(B);
    // std::cout<<A<<std::endl;
  };



};

template<typename T,Integer Rows,Integer Cols>
class Assignment<Transposed<Matrix<T,Cols,Rows>>> 
{
 public:       
 
  inline static constexpr void apply(Matrix<T,Rows,Cols>& A,const Transposed<Matrix<T,Cols,Rows>>& B)
  {
   for(Integer ii=0;ii<Rows;ii++)
    for(Integer jj=0;jj<Cols;jj++)
      Assignment<T>::apply(A(ii,jj),B(jj,ii));
  };

 inline static constexpr void apply(Transposed<Matrix<T,Cols,Rows>>& A,const Transposed<Matrix<T,Cols,Rows>>& B)
  {
  // update pointer to the non-transposed matrix of the transposed matrix A 
  // with the not transposed matrix of B
   A(B());
   };


};


template<typename S,Integer NQPoints>
class Assignment<QPValues<S,NQPoints>> 
{
 public:       
 template<typename T>
 inline static constexpr void apply(QPValues<T,NQPoints>& A,const QPValues<S,NQPoints>& B)
  {

    // T ok1(1);
    // S ok2(2);
   for(Integer ii=0;ii<NQPoints;ii++)
      {

        Assign(A[ii],B[ii]);

        // Assignment<S>::apply(A[ii],B[ii]);

        // Assignment<T>::apply(A[ii],B[ii]);
        std::cout<<"Assignment<QPValues<T,NQPoints>> "<<A[ii]<<" "<<B[ii]<<std::endl;
      }
  };







};

template<typename T,Integer NQPoints,Integer Ndofs>
class Assignment<FQPValues<T,NQPoints,Ndofs>> 
{
 public:       
 template<typename S>
 inline static constexpr void apply(FQPValues<S,NQPoints,Ndofs>& A,const FQPValues<T,NQPoints,Ndofs>& B)
  {
   for(Integer ii=0;ii<Ndofs;ii++)
    for(Integer jj=0;jj<NQPoints;jj++)
      Assign(A[ii][jj],B[ii][jj]);
  };
};


template<typename T,Integer Rows,Integer Cols>
class UnaryPlus<Matrix<T,Rows,Cols>> 
{    
 public:       

 inline static void apply(Matrix<T,Rows,Cols>& A,const Matrix<T,Rows,Cols>& B)
	{
	 for(Integer ii=0;ii<Rows;ii++)
	 	for(Integer jj=0;jj<Cols;jj++)
	 	   	   A(ii,jj)=B(ii,jj);
	};
};






template<typename Type>
class UnaryMinusMatrixAndTransposedAux
{
public:
 using Output=OperatorType<UnaryMinus<Type> >;
 inline static void apply(Output& A,const Type& B)
  {
    // Output oo(5,4,4,4,4);
    // Type ss(5,4,5,5);
    std::cout<<"pre UnaryMinusMatrixAndTransposedAux<Transposed"<<std::endl;
   for(Integer ii=0;ii<Type::Rows;ii++)
    for(Integer jj=0;jj<Type::Cols;jj++)
           A(ii,jj)=-B(ii,jj);
    std::cout<<A<<B<<std::endl;
    std::cout<<"after UnaryMinusMatrixAndTransposedAux<Transposed"<<std::endl;
  };
};



template<typename T,Integer Rows,Integer Cols>
class UnaryMinus<Matrix<T,Rows,Cols>> 
{    
 public:       
 inline static void apply(Matrix<T,Rows,Cols>& A,const Matrix<T,Rows,Cols>& B)
  {
   UnaryMinusMatrixAndTransposedAux<Matrix<T,Rows,Cols>>::apply(A,B);
  };

};

template<typename T,Integer Rows,Integer Cols>
class UnaryMinus<Transposed<Matrix<T,Cols,Rows>>> 
{    
 public:       
using Output=OperatorType<UnaryMinus<Transposed<Matrix<T,Cols,Rows>>>>;
 inline static void apply(Output& A,const Transposed<Matrix<T,Cols,Rows>>& B)
  {
   UnaryMinusMatrixAndTransposedAux<Transposed<Matrix<T,Cols,Rows>>>::apply(A,B);
   std::cout<<A<<B<<std::endl;
  };

};










template<typename Left,typename Right>
class AdditionMatrixAndTransposedAux
{
public:
 using Output=OperatorType<Addition<Left,Right> >;
 inline static void apply(Output& A,const Left& B,const Right& C)
  {
   for(Integer ii=0;ii<Left::Rows;ii++)
    for(Integer jj=0;jj<Left::Cols;jj++)
           A(ii,jj)=B(ii,jj)+C(ii,jj);
    std::cout<<"after AdditionMatrixAndTransposedAux<Transposed"<<std::endl;
  };
};


template<typename T,Integer Rows,Integer Cols>
class Addition<Matrix<T,Rows,Cols>,Matrix<T,Rows,Cols>> 
{       
 public: 
 using Left=Matrix<T,Rows,Cols>;
 using Right=Matrix<T,Rows,Cols>;
 using Output=OperatorType<Addition<Left,Right> >;

 inline static void apply(      Output& A,
                          const Left& B,
                          const Right& C)
  {
        AdditionMatrixAndTransposedAux<Left,Right>::apply(A,B,C);
  };
};

template<typename T,Integer Rows,Integer Cols>
class Addition<Transposed<Matrix<T,Cols,Rows>>,Matrix<T,Rows,Cols>> 
{       
 public: 
 using Left=Transposed<Matrix<T,Cols,Rows>>;
 using Right=Matrix<T,Rows,Cols>;
 using Output=OperatorType<Addition<Left,Right> >;

 inline static void apply(      Output& A,
                          const Left& B,
                          const Right& C)
  {
        AdditionMatrixAndTransposedAux<Left,Right>::apply(A,B,C);
  };
};


template<typename T,Integer Rows,Integer Cols>
class Addition<Matrix<T,Rows,Cols>,Transposed<Matrix<T,Cols,Rows>>> 
{       
 public: 
 using Left=Matrix<T,Rows,Cols>;
 using Right=Transposed<Matrix<T,Cols,Rows>>;
 using Output=OperatorType<Addition<Left,Right> >;


 inline static void apply(      Output& A,
                          const Left& B,
                          const Right& C)
  {
        AdditionMatrixAndTransposedAux<Left,Right>::apply(A,B,C);
  };
};

template<typename T,Integer Rows,Integer Cols>
class Addition<Transposed<Matrix<T,Cols,Rows>>,Transposed<Matrix<T,Cols,Rows>>> 
{       
 public: 
 using Left=Transposed<Matrix<T,Cols,Rows>>;
 using Right=Transposed<Matrix<T,Cols,Rows>>;
 using Output=OperatorType<Addition<Left,Right> >;

 inline static void apply(      Output& A,
                          const Left& B,
                          const Right& C)
  {
        AdditionMatrixAndTransposedAux<Left,Right>::apply(A,B,C);
  };
};
















template<typename Left,typename Right>
class SubtractionMatrixAndTransposedAux
{
public:
 using Output=OperatorType<Addition<Left,Right> >;
 inline static void apply(Output& A,const Left& B,const Right& C)
  {
   for(Integer ii=0;ii<Left::Rows;ii++)
    for(Integer jj=0;jj<Left::Cols;jj++)
           A(ii,jj)=B(ii,jj)-C(ii,jj);
    std::cout<<"after SubtractionMatrixAndTransposedAux<Transposed"<<std::endl;
  };
};




template<typename T,Integer Rows,Integer Cols>
class Subtraction<Matrix<T,Rows,Cols>,Matrix<T,Rows,Cols>> 
{       
 public: 
 using Left=Matrix<T,Rows,Cols>;
 using Right=Matrix<T,Rows,Cols>;
 using Output=OperatorType<Subtraction<Left,Right> >;

 inline static void apply(      Output& A,
                          const Left& B,
                          const Right& C)
  {
        SubtractionMatrixAndTransposedAux<Left,Right>::apply(A,B,C);
  };
};

template<typename T,Integer Rows,Integer Cols>
class Subtraction<Transposed<Matrix<T,Cols,Rows>>,Matrix<T,Rows,Cols>> 
{       
 public: 
 using Left=Transposed<Matrix<T,Cols,Rows>>;
 using Right=Matrix<T,Rows,Cols>;
 using Output=OperatorType<Subtraction<Left,Right> >;

 inline static void apply(      Output& A,
                          const Left& B,
                          const Right& C)
  {
        SubtractionMatrixAndTransposedAux<Left,Right>::apply(A,B,C);
  };
};


template<typename T,Integer Rows,Integer Cols>
class Subtraction<Matrix<T,Rows,Cols>,Transposed<Matrix<T,Cols,Rows>>> 
{       
 public: 
 using Left=Matrix<T,Rows,Cols>;
 using Right=Transposed<Matrix<T,Cols,Rows>>;
 using Output=OperatorType<Subtraction<Left,Right> >;


 inline static void apply(      Output& A,
                          const Left& B,
                          const Right& C)
  {
        SubtractionMatrixAndTransposedAux<Left,Right>::apply(A,B,C);
  };
};

template<typename T,Integer Rows,Integer Cols>
class Subtraction<Transposed<Matrix<T,Cols,Rows>>,Transposed<Matrix<T,Cols,Rows>>> 
{       
 public: 
 using Left=Transposed<Matrix<T,Cols,Rows>>;
 using Right=Transposed<Matrix<T,Cols,Rows>>;
 using Output=OperatorType<Subtraction<Left,Right> >;

 inline static void apply(      Output& A,
                          const Left& B,
                          const Right& C)
  {
        SubtractionMatrixAndTransposedAux<Left,Right>::apply(A,B,C);
  };
};






// template<typename T,Integer Rows,Integer Cols>
// class Subtraction<Matrix<T,Rows,Cols>,Matrix<T,Rows,Cols>> 
// {       
//  public:     

//  // inline static void apply(      Matrix<T,Rows,Cols>& A,
//  //                          const Matrix<T,Rows,Cols>& B)
//  //  {
//  //   for(Integer ii=0;ii<Rows;ii++)
//  //    for(Integer jj=0;jj<Cols;jj++)
//  //    {
//  //     A(ii,jj)=A(ii,jj)-B(ii,jj);
//  //    }
//  //  };

//  inline static void apply(      Matrix<T,Rows,Cols>& A,
//                           const Matrix<T,Rows,Cols>& B,
//                           const Matrix<T,Rows,Cols>& C)
//   {
//    for(Integer ii=0;ii<Rows;ii++)
//     for(Integer jj=0;jj<Cols;jj++)
//     {
//      A(ii,jj)=B(ii,jj)-C(ii,jj);
//     }
//   };
// };





template<typename Left, typename Right>
class MultiplicationMatrixAndTransposedAux
{
public:
 using Output=OperatorType<Multiplication<Left,Right> >;

 template<Integer M1,Integer N1,Integer M2, Integer N2>
 class apply_aux
 {
  public:
  inline constexpr static void apply(Output& A,const Left& B,const Right& C)
  {
   for(Integer ii=0;ii<Left::Rows;ii++)
    for(Integer jj=0;jj<Right::Cols;jj++)
       {
        A(ii,jj)=B(ii,0)*C(0,jj);
        for(Integer cc=1;cc<Left::Cols;cc++)
           A(ii,jj)+=B(ii,cc)*C(cc,jj);
       }
  }
 };

template<Integer M2, Integer N2>
 class apply_aux<1,1,M2,N2>
 {
  public:
  inline constexpr static void apply(Output& A,const Left& B,const Right& C)  
   {
   for(Integer ii=0;ii<Right::Rows;ii++)
    for(Integer jj=0;jj<Right::Cols;jj++)
        A(ii,jj)=B(0,0)*C(ii,jj);
  };
  
 };

template<Integer M1, Integer N1>
 class apply_aux<M1,N1,1,1>
 {
  public:
  inline constexpr static void apply(Output& A,const Left& B,const Right& C)  
   {
   for(Integer ii=0;ii<Left::Rows;ii++)
    for(Integer jj=0;jj<Left::Cols;jj++)
        A(ii,jj)=B(ii,jj)*C(0,0);
  };
  
 };

template<>
 class apply_aux<1,1,1,1>
 {
 public:
  inline constexpr static void apply(Output& A,const Left& B,const Right& C)  
   {A(0,0)=B(0,0)*C(0,0);};
  
 };

 inline constexpr static void apply(Output& A,
                                    const Left& B,
                                    const Right& C)
  {
    apply_aux<Left::Rows,Left::Cols,Right::Rows,Right::Cols>::apply(A,B,C);
  };

};




template<Integer Rows1,Integer Cols1,Integer Rows2,Integer Cols2>
class Multiplication<Matrix<Real,Rows1,Cols1>,
                     Matrix<Real,Rows2,Cols2>> 
{    
 public:      
 
 using Left=Matrix<Real,Rows1,Cols1>;
 using Right=Matrix<Real,Rows2,Cols2>;
 using Output=OperatorType<Multiplication<Left,Right> >;

 inline constexpr static void apply(Output& A,
                                    const Left& B,
                                    const Right& C)
  {
    MultiplicationMatrixAndTransposedAux<Left,Right>::apply(A,B,C);
  };


};

template<Integer Rows1,Integer Cols1,Integer Rows2,Integer Cols2>
class Multiplication<Matrix<Real,Rows1,Cols1>,
                     Transposed<Matrix<Real,Cols2,Rows2>>> 
{    
 public:      
 
 using Left=Matrix<Real,Rows1,Cols1>;
 using Right=Transposed<Matrix<Real,Cols2,Rows2>>;
 using Output=OperatorType<Multiplication<Left,Right> >;

inline constexpr static void apply(Output& A,
                                    const Left& B,
                                    const Right& C)
  {
    MultiplicationMatrixAndTransposedAux<Left,Right>::apply(A,B,C);
  };
};


template<Integer Rows1,Integer Cols1,Integer Rows2,Integer Cols2>
class Multiplication<Transposed<Matrix<Real,Cols1,Rows1>>,
                     Matrix<Real,Rows2,Cols2>> 
{    
 public:      
 
 using Left=Transposed<Matrix<Real,Cols1,Rows1>>;
 using Right=Matrix<Real,Rows2,Cols2>;
 using Output=OperatorType<Multiplication<Left,Right> >;

inline constexpr static void apply(Output& A,
                                    const Left& B,
                                    const Right& C)
  {
    MultiplicationMatrixAndTransposedAux<Left,Right>::apply(A,B,C);
  };
};

template<Integer Rows1,Integer Cols1,Integer Rows2,Integer Cols2>
class Multiplication<Transposed<Matrix<Real,Cols1,Rows1>>,
                     Transposed<Matrix<Real,Cols2,Rows2>>> 
{    
 public:      
 
 using Left=Transposed<Matrix<Real,Cols1,Rows1>>;
 using Right=Transposed<Matrix<Real,Cols2,Rows2>>;
 using Output=OperatorType<Multiplication<Left,Right> >;

inline constexpr static void apply(Output& A,
                                    const Left& B,
                                    const Right& C)
  {
    MultiplicationMatrixAndTransposedAux<Left,Right>::apply(A,B,C);
  };
};




template<typename Left,typename Right>
class DivisionMatrixAndTransposedAux
{
public:
 using Output=OperatorType<Division<Left,Right> >;
 inline static void apply(Output& A,const Left& B,const Right& C)
  {
   for(Integer ii=0;ii<Left::Rows;ii++)
    for(Integer jj=0;jj<Left::Cols;jj++)
           A(ii,jj)=B(ii,jj)/C(0,0);
  };
};




template<Integer Rows,Integer Cols>
class Division<Matrix<Real,Rows,Cols>,
               Matrix<Real,1,1>> 

{
  public:
  using Left=Matrix<Real,Rows,Cols>;
  using Right=Matrix<Real,1,1>;
  using Output=OperatorType<Multiplication<Left,Right> >;

  inline static void apply(Output& A,const Left& B,const Right& C)
  {
  DivisionMatrixAndTransposedAux<Left,Right>::apply(A,B,C);
  }
};


template<Integer Rows,Integer Cols>
class Division<Transposed<Matrix<Real,Rows,Cols>>,
               Matrix<Real,1,1>> 

{
  public:
  using Left=Transposed<Matrix<Real,Rows,Cols>>;
  using Right=Matrix<Real,1,1>;
  using Output=OperatorType<Multiplication<Left,Right> >;

  inline static void apply(Output& A,const Left& B,const Right& C)
  {
  DivisionMatrixAndTransposedAux<Left,Right>::apply(A,B,C);
  }
};


template<Integer Rows,Integer Cols>
class Division<Matrix<Real,Rows,Cols>,
               Transposed<Matrix<Real,1,1>>> 

{
  public:
  using Left=Matrix<Real,Rows,Cols>;
  using Right=Transposed<Matrix<Real,1,1>>;
  using Output=OperatorType<Multiplication<Left,Right> >;

  inline static void apply(Output& A,const Left& B,const Right& C)
  {
  DivisionMatrixAndTransposedAux<Left,Right>::apply(A,B,C);
  }
};



template<Integer Rows,Integer Cols>
class Division<Transposed<Matrix<Real,Rows,Cols>>,
               Transposed<Matrix<Real,1,1>>> 

{
  public:
  using Left=Transposed<Matrix<Real,Rows,Cols>>;
  using Right=Transposed<Matrix<Real,1,1>>;
  using Output=OperatorType<Multiplication<Left,Right> >;

  inline static void apply(Output& A,const Left& B,const Right& C)
  {
  DivisionMatrixAndTransposedAux<Left,Right>::apply(A,B,C);
  }
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////         QP   EXPRESSION ALGEBRAIC OPERATIONS        ///////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename U,Integer NQPoints>
class UnaryPlus<QPValues<U,NQPoints>> 
{       
 public:     
 template<typename V>          
 inline static void apply(      QPValues<V,NQPoints>& A,
                          const QPValues<U,NQPoints>& B)
  {
   for(Integer qp=0;qp<NQPoints;qp++)
        {
          // UnaryPlus<U>::apply(A[qp],B[qp]);
          Plus(A[qp],B[qp]);
        }
  };
};











template<typename U,Integer NQPoints>
class UnaryMinus<QPValues<U,NQPoints>> 
{       
 public:               
 template<typename V>
 inline static void apply(      QPValues<V,NQPoints>& A,
                          const QPValues<U,NQPoints>& B)
  {
   for(Integer qp=0;qp<NQPoints;qp++)
        {
          // UnaryMinus<U>::apply(A[qp],B[qp]);
          Minus(A[qp],B[qp]);
        }
  };
};



template<typename U,typename V, Integer NQPoints>
class Addition<QPValues<U,NQPoints>,
               QPValues<V,NQPoints>> 
{       
 public:               
 inline static void apply(      QPValues<OperatorType<Addition<U,V>>,NQPoints>& A,
                          const QPValues<U,NQPoints>& B,
                          const QPValues<V,NQPoints>& C)
  {
   for(Integer qp=0;qp<NQPoints;qp++)
    {
      // Addition<U,V>::apply(A[qp],B[qp],C[qp]);
      Add(A[qp],B[qp],C[qp]);
    }    
    
  };

};

template<typename U,typename V, Integer NQPoints>
class Subtraction<QPValues<U,NQPoints>,
                  QPValues<V,NQPoints>> 
{       
 public:               
 inline static void apply(      QPValues<OperatorType<Subtraction<U,V>>,NQPoints>& A,
                          const QPValues<U,NQPoints>& B,
                          const QPValues<V,NQPoints>& C)
  {
   for(Integer qp=0;qp<NQPoints;qp++)
        {
          // Subtraction<U,V>::apply(A[qp],B[qp],C[qp]);
          Subtract(A[qp],B[qp],C[qp]);
        }
  };

};


template<typename U, typename V,Integer NQPoints>
class Multiplication<QPValues<U,NQPoints>,
                     QPValues<V,NQPoints> > 
{       
 public:          
 inline static void apply(      QPValues<OperatorType<Multiplication<U,V>>,NQPoints>& A,
                          const QPValues<U,NQPoints>& B,
                          const QPValues<V,NQPoints>& C)
  {
    // OperatorType<Multiplication<U,V>> oo(6,7,8,5);
    // U ok(1);
    // V ok2(3);
   for(Integer qp=0;qp<NQPoints;qp++)
      {
        // Multiplication<U,V>::apply(A[qp],B[qp],C[qp]);
        Multiply(A[qp],B[qp],C[qp]);
      }
  };
};


template<typename U, typename V,Integer NQPoints>
class Division<QPValues<U,NQPoints>,
               QPValues<V,NQPoints> > 
{       
 public:          
 inline static void apply(      QPValues<OperatorType<Division<U,V>>,NQPoints>& A,
                          const QPValues<U,NQPoints>& B,
                          const QPValues<V,NQPoints>& C)
  {
   for(Integer qp=0;qp<NQPoints;qp++)
        {
          // Division<U,V>::apply(A[qp],B[qp],C[qp]);
          Divide(A[qp],B[qp],C[qp]);
        }

  };
};



///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////         FQP   EXPRESSION ALGEBRAIC OPERATIONS       ///////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename U, Integer NQPoints,Integer Ndofs>
class UnaryPlus< FQPValues<U,NQPoints,Ndofs>> 
{       
 public:               
 inline static void apply(FQPValues<U,NQPoints,Ndofs>& A,const FQPValues<U,NQPoints,Ndofs>& B)
  {
  for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     for(Integer qp=0;qp<NQPoints;qp++)
       {
        // UnaryPlus<U>::apply(A(n_dof,qp),B(n_dof,qp));
        Plus(A(n_dof,qp),B(n_dof,qp));
       }
  };
};

template<typename U, Integer NQPoints,Integer Ndofs>
class UnaryMinus<FQPValues<U,NQPoints,Ndofs>> 
{       
 public:               
 template<typename V>
 inline static void apply(FQPValues<V,NQPoints,Ndofs>& A,const FQPValues<U,NQPoints,Ndofs>& B)
  {
    // FQPValues<U,NQPoints,Ndofs> ok1(1);
    // FQPValues<U,NQPoints,Ndofs> ok2(2);
 // U ok(1,2,3,4,5,6,7);// transposed 
 // V ok3(4,5,6); // transposed
  std::cout<<"__________pre UnaryMinus< FQPValues<U,NQPoints,Ndofs>> "<<std::endl;
  for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     for(Integer qp=0;qp<NQPoints;qp++)
       {
        // UnaryMinus<U>::apply(A(n_dof,qp),B(n_dof,qp));
        Minus(A(n_dof,qp),B(n_dof,qp));
        // std::cout<<"B(n_dof,qp)="<<B(n_dof,qp)<<std::endl;
        // std::cout<<"A(n_dof,qp)="<<A(n_dof,qp)<<std::endl;
       }
  std::cout<<"__________after UnaryMinus< FQPValues<U,NQPoints,Ndofs>> "<<std::endl;
  // std::cout<<B<<std::endl;
  // std::cout<<A<<std::endl;
  };
};

template<typename U, typename V, Integer NQPoints,Integer Ndofs>
class Addition< FQPValues<U,NQPoints,Ndofs>,
                FQPValues<V,NQPoints,Ndofs>> 
{       
 public:               
 inline static void apply(      FQPValues<OperatorType<Addition<U,V>>,NQPoints,Ndofs>& A,
                          const FQPValues<U,NQPoints,Ndofs>& B,
                          const FQPValues<V,NQPoints,Ndofs>& C)
  {
  for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
   for(Integer qp=0;qp<NQPoints;qp++)
    {
      // Addition<U,V>::apply(A(n_dof,qp),B(n_dof,qp),C(n_dof,qp));
      Add(A(n_dof,qp),B(n_dof,qp),C(n_dof,qp));
    }
  };


};

template<typename U, typename V, Integer NQPoints,Integer Ndofs>
class Subtraction< FQPValues<U,NQPoints,Ndofs>,
                   FQPValues<V,NQPoints,Ndofs>> 
{       
 public:               
 inline static void apply(      FQPValues<OperatorType<Subtraction<U,V>>,NQPoints,Ndofs>& A,
                          const FQPValues<U,NQPoints,Ndofs>& B,
                          const FQPValues<V,NQPoints,Ndofs>& C)
  {
  for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
   for(Integer qp=0;qp<NQPoints;qp++)
      {
        // Subtraction<U,V>::apply(A(n_dof,qp),B(n_dof,qp),C(n_dof,qp));
        Subtract(A(n_dof,qp),B(n_dof,qp),C(n_dof,qp));
      }
  };
};






// template<template<class...>class Nary, typename...Ts>
// class UnaryMinus< Nary<Ts...>> 
// {       
//  public:              
//  using type= 
//  inline static void apply(FQPValues<U,NQPoints,Ndofs>& A,const FQPValues<U,NQPoints,Ndofs>& B)
//   {
//   for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
//      for(Integer qp=0;qp<NQPoints;qp++)
//        {UnaryMinus<U>::apply(A(n_dof,qp),B(n_dof,qp));}
//   };
// };




















template<typename S, typename T,Integer NQPoints,Integer Ndofs>
class Multiplication<QPValues<S,NQPoints>,
                     FQPValues<T,NQPoints,Ndofs>
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
        Multiply(A(n_dof,qp),B[qp],C(n_dof,qp));
        // std::cout<<"B[qp]"<<B[qp]<<std::endl;
        // std::cout<<"C[qp]"<<C[n_dof][qp]<<std::endl;
        // std::cout<<"A[qp]"<<A[n_dof][qp]<<std::endl;
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
        Multiply(A(n_dof,qp),B(n_dof,qp),C[qp]);
       }
  };
};



template<typename S, typename T,Integer NQPoints,Integer Ndofs>
class Division<QPValues<S,NQPoints>,
               FQPValues<T,NQPoints,Ndofs>>
                      
{    
 public:            
 inline static void apply(      FQPValues<OperatorType<Division<S,T>>,NQPoints,Ndofs>& A,
                          const QPValues<S,NQPoints>& B,
                          const FQPValues<T,NQPoints,Ndofs>& C)
  {
   for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
    for(Integer qp=0;qp<NQPoints;qp++)
       {
           //FIXME TODO
        //FIXME TODO
        Divide(A(n_dof,qp),B[qp],C(n_dof,qp));
        // std::cout<<"B[qp]"<<B[qp]<<std::endl;
        // std::cout<<"C[qp]"<<C[n_dof][qp]<<std::endl;
        // std::cout<<"A[qp]"<<A[n_dof][qp]<<std::endl;
       }
  };
};


template<typename S, typename T,Integer NQPoints,Integer Ndofs>
class Division<FQPValues<S,NQPoints,Ndofs>,
               QPValues<T,NQPoints>> 
{    
 public:          
 inline static void apply(      FQPValues<OperatorType<Division<S,T>>,NQPoints,Ndofs>& A,
                          const FQPValues<S,NQPoints,Ndofs>& B,
                          const QPValues<T,NQPoints>& C)
  {
   for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
    for(Integer qp=0;qp<NQPoints;qp++)
       {
        //FIXME TODO
        //FIXME TODO
        Divide(A(n_dof,qp),B(n_dof,qp),C[qp]);
        // std::cout<<"B[qp]"<<B[qp][n_dof]<<std::endl;
        // std::cout<<"C[qp]"<<C[qp]<<std::endl;
        // std::cout<<"A[qp]"<<A[n_dof][qp]<<std::endl;
       }
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

template<typename S,typename T>
constexpr void Assign(S& t,const T& constt){Assignment<T>::apply(t,constt);}

template<typename T>
constexpr void Plus(T& t,const T& constt){UnaryPlus<T>::apply(t,constt);}

template<typename U,typename V>
constexpr void Minus(U& t,const V& constt){UnaryMinus<V>::apply(t,constt);}

template<typename U,typename V,typename W>
constexpr void Add(U& u, const V& v, const W& w){Addition<V,W>::apply(u,v,w);}

template<typename U,typename V,typename W>
constexpr void Subtract(U& u, const V& v, const W& w){Subtraction<V,W>::apply(u,v,w);}

template<typename U,typename V,typename W>
constexpr void Multiply(U& u, const V& v, const W& w){Multiplication<V,W>::apply(u,v,w);}

template<typename U, typename V,typename W>
constexpr void Divide(U& u, const V& v, const W& w){Division<V,W>::apply(u,v,w);}







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
    const constexpr Derived& operator()()const{return value_;};
    Derived& derived(){return value_;};
    const constexpr Derived& derived()const{return value_;};
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
    const constexpr Derived& operator()()const{return value_;};
    Derived& derived(){return value_;};
    const constexpr Derived& derived()const{return value_;};
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
    const constexpr Left& left()const{return left_;};
    const constexpr Right& right()const{return right_;};
    // todo fixme this is to return the same object 
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
    const constexpr DerivedLeft& left()const{return left_;};
    const constexpr DerivedRight& right()const{return right_;};
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
    const constexpr DerivedLeft& left()const{return left_;};
    const constexpr DerivedRight& right()const{return right_;};
  private:
  DerivedLeft left_;
  DerivedRight right_;
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
    {};
    const constexpr DerivedLeft& left()const{return left_;};
    const constexpr DerivedRight& right()const{return right_;};
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
    // Multiplication(const DerivedLeft& left, const DerivedRight&right)
    // : 
    // left_(left),
    // right_(right)
    // {};
    const constexpr DerivedLeft& left()const{return left_;};
    const constexpr DerivedRight& right()const{return right_;};
  private:
  DerivedLeft left_;
  DerivedRight right_;
};







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

template< typename DerivedLeft,typename DerivedRight>
class Multiplication< Expression <DerivedLeft>, Expression <DerivedRight> >
operator*(const Expression<DerivedLeft>&left, const Expression<DerivedRight>&right)
{return Multiplication< Expression <DerivedLeft>, Expression <DerivedRight> >(left,right);}

template< typename DerivedLeft,typename DerivedRight>
class Division< Expression <DerivedLeft>, Expression <DerivedRight> >
operator/(const Expression<DerivedLeft>&left, const Expression<DerivedRight>&right)
{return Division< Expression <DerivedLeft>, Expression<DerivedRight> > (left,right);}

template< typename DerivedLeft,typename DerivedRight>
class Contraction2< Expression <DerivedLeft>, Expression <DerivedRight> >
Inner(const Expression<DerivedLeft>&left, const Expression<DerivedRight>&right)
{return Contraction2< Expression <DerivedLeft>, Expression <DerivedRight> >(left,right);}








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
 // using subtype= OperatorType<type,OtherTemplateArguments...>;
 using subtype= OperatorType<Derived,OtherTemplateArguments...>;
 using outputsubtype= OperatorType<type,OtherTemplateArguments...>;

 Evaluation(){};
 

 Evaluation(const Expression<type>& expr):
 expr_(expr.derived()),
 eval_(Eval(expr_()))
 {};
 
 template<typename...Inputs>
 void apply(subtype& output,const Inputs&...inputs)
 {
  // subtype ok1(1);
  // outputsubtype ok2(2);
  // compute output
  // eval_.template apply<OtherTemplateArguments2...>(output,inputs...);
  
  // subtype ok1;//(1,2,3,8,8,9);
  // type ok2;//(2,4,5,6,6,7,8);
  eval_.apply(output,inputs...);
  std::cout<<"Evaluation<Expression<UnaryPlus "<<output<<std::endl;
  // apply plus to output (we assume it does not change, so we do not compute)
  // Derived(output);
 }
private:

 type expr_;
 Eval eval_;
};


template<typename Derived,typename...OtherTemplateArguments>
class Evaluation<Expression<UnaryMinus<Expression<Derived> > >,OtherTemplateArguments...>
{
 public:
 using type=UnaryMinus<Expression<Derived>>;
 using Eval=Evaluation<Expression<Derived>,OtherTemplateArguments...>;
 using subtype= OperatorType<Derived,OtherTemplateArguments...>;
 using outputsubtype= OperatorType<type,OtherTemplateArguments...>;
  // using subtype= OperatorType<type,OtherTemplateArguments...>;

 Evaluation(){};
 
 Evaluation(const Expression<type>& expr):
 expr_(expr.derived()),
 eval_(Eval(expr_()))
 {};
 
 template<typename...Inputs>
 void apply(outputsubtype& output,const Inputs&...inputs)
 {
  
  // compute output
  // outputsubtype eee(6);
  // subtype oooo(5);
    // eval_.template apply<OtherTemplateArguments2...>(output,inputs...);
   std::cout<<"pre eval Evaluation<Expression<UnaryMinus "<<std::endl;
  eval_.apply(derived_value_,inputs...);
  std::cout<<"pre MINUS Evaluation<Expression<UnaryMinus "<<std::endl;
  Minus(output,derived_value_);
  std::cout<<"after MINUS Evaluation<Expression<UnaryMinus "<<std::endl;
  std::cout<<"derived_value_="<<derived_value_<<std::endl;
  std::cout<<"output="<<output<<std::endl;
  // apply minus to output
  // Minus(output,tmp);
 }

private:

 type expr_;
 Eval eval_;
 subtype derived_value_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Specialization for -(-X)=X. The evaluation is simply X.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

// template<typename Derived,typename...OtherTemplateArguments>
// class Evaluation<Expression<UnaryMinus< Expression<UnaryMinus< Expression<Derived> >> > >,OtherTemplateArguments... >
// {
//  public:
//  using type=UnaryPlus<Expression<Derived>>;
//  using Eval=Evaluation<Expression<Derived>,OtherTemplateArguments...>;
//  using subtype= OperatorType<type,OtherTemplateArguments...>;
//  Evaluation(){};
 

//  Evaluation(const Expression<type>& expr):
//  expr_(expr.derived()),
//  eval_(Eval(expr_()))
//  {};
 
//  template<typename...Inputs>
//  void apply(subtype& output,const Inputs&...inputs)
//  {
//   // compute output
//   // eval_.template apply<OtherTemplateArguments2...>(output,inputs...);
//   std::cout<<"pre Evaluation<Expression<UnaryMinus<UnaryMinus "<<output<<std::endl;
//   eval_.apply(output,inputs...);
//   std::cout<<"after Evaluation<Expression<UnaryMinus<UnaryMinus "<<output<<std::endl;
//   // apply plus to output (we assume it does not change, so we do not compute)
//   // Derived(output);
//  }
// private:

//  type expr_;
//  Eval eval_;
// };




template<typename DerivedRight,typename DerivedLeft,typename...OtherTemplateArguments>
class Evaluation<Expression<Addition< Expression<DerivedLeft>  ,  
                                      Expression<DerivedRight>>>,
                 OtherTemplateArguments...>
{
 public:
 using type=Addition<Expression<DerivedLeft>,Expression<DerivedRight>>;
 using EvalLeft=Evaluation<Expression<DerivedLeft>,OtherTemplateArguments...>;
 using EvalRight=Evaluation<Expression<DerivedRight>,OtherTemplateArguments...>;
 using subtypeleft=OperatorType<DerivedLeft,OtherTemplateArguments...>;
 using subtyperight=OperatorType<DerivedRight,OtherTemplateArguments...>;
 using subtype=OperatorType<Addition<subtypeleft,subtyperight>,OtherTemplateArguments...>;
 // using subtype=OperatorType<type,OtherTemplateArguments...>;
 // todo fixmeeeeeeeeeeee fix me FIX ME
 Evaluation(){};
 

 // Evaluation(const Expression<type>& expr):
 Evaluation(const Expression<type>& expr,OtherTemplateArguments&...args):
 expr_(expr.derived()),
 // eval_left_(EvalLeft(expr_.left())),
 // eval_right_(EvalRight(expr_.right()))
 eval_left_(EvalLeft(expr_.left(),args...)),
 eval_right_(EvalRight(expr_.right(),args...))
{};

 Evaluation(const Expression<type>& expr):
 expr_(expr.derived()),
 // eval_left_(EvalLeft(expr_.left())),
 // eval_right_(EvalRight(expr_.right()))
 eval_left_(EvalLeft(expr_.left())),
 eval_right_(EvalRight(expr_.right()))
{};

 template<typename...OtherTemplateArguments2,typename...Inputs>
 void apply(subtype& output,const Inputs&...inputs)
 {
  
  eval_left_.apply(left_value_,inputs...);
  eval_right_.apply(right_value_,inputs...);
  // EvalRight oo(5);


// type oo(6);  
  // subtypeleft ok1(1);
  // subtyperight ok2(3);
  //  subtype ok4(5);
Add(output,left_value_,right_value_);

  std::cout<<"left_value= "<<left_value_<<std::endl;
  std::cout<<"right_value_= "<<right_value_<<std::endl;
  std::cout<<"Evaluation<Expression<Addition "<<output<<std::endl;

 }
private:

 type expr_;
 subtypeleft left_value_;
 subtyperight right_value_;
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
 // using subtype=OperatorType<type,OtherTemplateArguments...>;
 using subtypeleft=OperatorType<DerivedLeft,OtherTemplateArguments...>;
 using subtyperight=OperatorType<DerivedRight,OtherTemplateArguments...>;
 using subtype=OperatorType<Subtraction<subtypeleft,subtyperight>,OtherTemplateArguments...>;

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

  // subtyperight ok1(1);
  // subtypeleft ok2(3);
  std::cout<<"Evaluation<Expression<Subtraction "<<output<<std::endl;

 }
private:

 type expr_;
 subtypeleft left_value_;
 subtyperight right_value_;
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
  std::cout<<"Evaluation<Expression<Multiplication "<<output<<std::endl;
 }
private:

 type expr_;
 subtype output;
 subtypeleft left_value_;
 subtyperight right_value_;
 EvalLeft eval_left_;
 EvalRight eval_right_;
};

// template<typename DerivedRight,typename...OtherTemplateArguments>
// class Evaluation<Expression<Multiplication< Real, Expression<DerivedRight> >>,OtherTemplateArguments...>                                   
// {
//  public:
//  using type=Multiplication<Real,Expression<DerivedRight>>;
//  using EvalRight=Evaluation<Expression<DerivedRight>,OtherTemplateArguments...>;
//  using subtype=OperatorType<type,OtherTemplateArguments...>;
//  Evaluation(){};
 

//  Evaluation(const Expression<type>& expr):
//  expr_(expr.derived()),
//  eval_left_(expr_.left()),
//  eval_right_(EvalRight(expr_.right()))
// {};
 
//  template<typename...OtherTemplateArguments2,typename...Inputs>
//  void apply(subtype& output,const Inputs&...inputs)
//  {
//   // std::cout<<"evaluation Multiplication(Real,)"<<std::endl;
//   eval_right_.apply(output,inputs...);
//   Multiply(output,eval_left_);
//  }
// private:

//  type expr_;
//  Real eval_left_;
//  EvalRight eval_right_;
// };




// template< typename DerivedLeft,typename...OtherTemplateArguments>
// class Evaluation<Expression<Multiplication< Expression<DerivedLeft>,Real>>,
//                  OtherTemplateArguments...>                                       
// {
//  public:
//  using type=Multiplication<Expression<DerivedLeft>,Real>;
//  using EvalLeft=Evaluation<Expression<DerivedLeft>,OtherTemplateArguments...>;
//  using subtype=OperatorType<type,OtherTemplateArguments...>;
 
//  Evaluation(){};
 

//  Evaluation(const Expression<type>& expr):
//  expr_(expr.derived()),
//  eval_left_(EvalLeft(expr_.left())),
//  eval_right_(expr_.right())
// {};
 
//  template<typename...OtherTemplateArguments2,typename...Inputs>
//  void apply(subtype& output,const Inputs&...inputs)
//  {
//   // std::cout<<"evaluation Multiplication(,Real)"<<std::endl;
//   eval_left_.apply(output,inputs...);
//   Multiply(output,eval_right_);

//  }
// private:
//  type expr_;
//  EvalLeft eval_left_;
//  Real eval_right_;
// };

template< typename DerivedLeft, typename DerivedRight, typename...OtherTemplateArguments>
class Evaluation<Expression<Division< Expression<DerivedLeft>, Expression<DerivedRight>> >,
                 OtherTemplateArguments...>                                       
{
 public:
 using type=Division<Expression<DerivedLeft>,Expression<DerivedRight>>;
 using EvalLeft=Evaluation<Expression<DerivedLeft>,OtherTemplateArguments...>;
 using EvalRight=Evaluation<Expression<DerivedRight>,OtherTemplateArguments...>;
 using subtype=OperatorType<type,OtherTemplateArguments...>;
 using subtypeleft=OperatorType<DerivedLeft,OtherTemplateArguments...>;
 using subtyperight=OperatorType<DerivedRight,OtherTemplateArguments...>;

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
  Divide(output,left_value_,right_value_);
  std::cout<<"Evaluation<Expression<Division "<<output<<std::endl;

 }

private:

 type expr_;
 subtypeleft left_value_;
 subtyperight right_value_;
 EvalLeft eval_left_;
 EvalRight eval_right_;
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
 general_form_(general_form)
 ,
 shapesform_(shapesform)
 ,
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
};


template<typename Form,typename FullSpace, typename GeneralForm_, typename...GeneralForms>
constexpr auto Eval(const GeneralForm<Form>& form, ShapeFunctions2<GeneralForm_,GeneralForms...>& shapes)
{return Evaluation<Expression<GeneralForm<Form>>,ShapeFunctions2<GeneralForm_,GeneralForms...> >(form,shapes);}











































































































template<typename...Ts>
class CompositeOperator;

template<typename MixedSpace, Integer N, typename OperatorType>
class Test;

template<typename MixedSpace, Integer N, typename OperatorType>
class Trial;

template<typename Expr>
class CompositeOperator<Expression<Expr>>
{
public:
  CompositeOperator(const Expr& expr):expr_(expr){}


template<typename MixedSpace,Integer N> 
constexpr Test<MixedSpace,N,CompositeOperator> 
operator()(const Test<MixedSpace,N,IdentityOperator>& t)
{return Test<MixedSpace,N,CompositeOperator> (t.spaces_ptr(),expr_);}

template<typename MixedSpace,Integer N> 
constexpr Trial<MixedSpace,N,CompositeOperator> 
operator()(const Trial<MixedSpace,N,IdentityOperator>& t)
{return Trial<MixedSpace,N,CompositeOperator> (t.spaces_ptr(),expr_);}



constexpr auto& composite_operator()const{return expr_;}
private:
  Expr expr_;
};

template<typename Expr>
constexpr auto NewOperator(const Expr& expr){return CompositeOperator<Expression<Expr>>(expr);}
























template<typename T>
constexpr auto Transpose(const Expression<T>& t){return Transposed<Expression<T>>(t);}

template<typename T>
constexpr auto Transpose(const Expression<Transposed<Expression<T>>>& t){return T(t.derived().derived());}


// template< typename Left_,typename Right_>
// class Transposed< Expression<Addition< Expression <Left_>, Expression <Right_> > > >
// : public Expression< Transposition< Expression<Addition< Expression <Left_>, Expression <Right_> > > > >
// {
//  public:
//     using Left=Transposition<Left_>;
//     using Right=Transposition<Right_>;

//     Transposition(const Expression< Addition< Expression <Left_>,Expression <Right_> > >& addition)
//     : 
//     left_(Transpose(addition.left())),
//     right_(Transpose(addition.right()))
//     {};

//     Transposition(const Expression <Left_>& left, const Expression <Right_> & right)
//     : 
//     left_(Transpose(left)),
//     right_(Transpose(right))
//     {};


//     const constexpr Left& left()const{return left_;};
//     const constexpr Right& right()const{return right_;};

//     auto operator()()
//     {return Transposition<Expression<Addition<Expression<Left_>,Expression<Right_>>>>(left_,right_);};

//     auto operator()()const
//     {return Transposition<Expression<Addition<Expression<Left_>,Expression<Right_>>>>(left_,right_);};
//   private:
//   Left left_;
//   Right right_;
// };



// template<Integer Rows1,Integer Cols1,Integer Rows2,Integer Cols2>
// class Multiplication<Transposition<Matrix<Real,Rows1,Cols1>>,
//                      Matrix<Real,Rows2,Cols2>> 
// {    
//  public:      
//  using Output=OperatorType<Multiplication<Matrix<Real,Rows1,Cols1>,Matrix<Real,Rows2,Cols2>> >;
//  using Left=Matrix<Real,Rows1,Cols1>;
//  using Right=Matrix<Real,Rows2,Cols2>;

// };

// template< typename Derived>
// class UnaryPlus< Expression <Derived> > 
// : public Expression< UnaryPlus< Expression <Derived> > >
// {
//   public:
//     UnaryPlus(const Expression<Derived>& expr): value_(expr.derived()){};
//     Derived& operator()(){return value_;};
//     const constexpr Derived& operator()()const{return value_;};
//     Derived& derived(){return value_;};
//     const constexpr Derived& derived()const{return value_;};
//   private:
//   Derived value_;
// };


template<typename T>
class Transposed<Expression<T>>: public Expression<Transposed<Expression<T>>>
{
public:

    Transposed(const Expression<T>& expr): value_(expr.derived()){};
    T& operator()(){return value_;};
    const constexpr T& operator()()const{return value_;};
    T& derived(){return value_;};
    const constexpr T& derived()const{return value_;};


  private:
  T value_;
};

// template<typename T>
// class Transposed<Expression<Transposed<Expression<T>>>>: public Expression<Expression<Transposed<Expression<T>>>>
// {
// public:

//     Transposed(const Expression<T>& expr): value_(expr.derived()){};
//     T& operator()(){return value_;};
//     const constexpr T& operator()()const{return value_;};
//     T& derived(){return value_;};
//     const constexpr T& derived()const{return value_;};


//   private:
//   T value_;
// };


template<typename ConstType, typename...Ts, typename QuadratureRule>
class Evaluation<Expression<Transposed<Expression<ConstantTensor<ConstType,Ts...>>>>,QuadratureRule>
{
public:
 using T=ConstantTensor<ConstType,Ts...>;
 using Eval=Evaluation<Expression<T>,QuadratureRule>;
 // // here we must transpose
 using type=Transposed<Expression<T>>;
 using subtype_not_transposed= OperatorType<T,QuadratureRule>;
 using subtype= OperatorType<type,QuadratureRule>;
 // Evaluation(){};
  Evaluation(const Expression<type>& expr)
 :
 expr_(expr)
 ,
 eval_(Eval(expr.derived().derived().derived())),
 const_(expr.derived().derived().derived())
 {};

 
 template<typename...Inputs>
 void apply(subtype& output,const Inputs&...inputs)
 {

  // const_.template qp_eval<QuadratureRule::NQPoints>();

  // eval_.apply(output_tmp_,inputs...);
  std::cout<<"Evaluation<Expression<Transposed<Expression<ConstantTensor"<<std::endl;
  std::cout<<QuadratureRule::NQPoints<<std::endl;
  std::cout<<const_.template qp_eval<QuadratureRule::NQPoints>()<<std::endl;
  // Assignment<subtype>::apply(output,const_.template qp_eval<QuadratureRule::NQPoints>());
  Assign(output,const_.template qp_eval<QuadratureRule::NQPoints>());


   std::cout<<"output="<<output<<std::endl;

  // Assignment<value_type>::apply(value,const_.template qp_eval<QuadratureRule::NQPoints>());
 }
private:
 subtype_not_transposed output_tmp_;
 type expr_;
 Eval eval_;
 T const_;

};



template<typename T,typename...OtherTemplateArguments>
class Evaluation<Expression<Transposed<Expression<T>>>,OtherTemplateArguments...>
{
public:

 using Eval=Evaluation<Expression<T>,OtherTemplateArguments...>;
 // // here we must transpose
 using type=Transposed<Expression<T>>;
 using subtype_not_transposed= OperatorType<T,OtherTemplateArguments...>;
 using subtype= OperatorType<type,OtherTemplateArguments...>;
 // Evaluation(){};
  Evaluation(const Expression<type>& expr)
 :
 expr_(expr)
 ,
 eval_(Eval(expr.derived().derived()))
 {};

 
 template<typename...Inputs>
 void apply(subtype& output,const Inputs&...inputs)
 {

  eval_.apply(output_tmp_,inputs...);

  // Assignment<subtype>::apply(output,output_tmp_);
  Assign(output,output_tmp_);

 }
private:
 subtype_not_transposed output_tmp_;
 type expr_;
 Eval eval_;

};





}

#endif //MARS_OPERATORS_HPP
