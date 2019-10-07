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

template<typename T,Integer Dim>
class Assignment<Vector<T,Dim>> 
{
 public:       
 inline static constexpr void apply(Vector<T,Dim>& A,const Vector<T,Dim>& B)
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
};

template<typename T,Integer NQPoints>
class Assignment<QPValues<T,NQPoints>> 
{
 public:       
 inline static constexpr void apply(QPValues<T,NQPoints>& A,const QPValues<T,NQPoints>& B)
  {
   for(Integer ii=0;ii<NQPoints;ii++)
      Assignment<T>::apply(A[ii],B[ii]);
  };
};

template<typename T,Integer NQPoints,Integer Ndofs>
class Assignment<FQPValues<T,NQPoints,Ndofs>> 
{
 public:       
 inline static constexpr void apply(FQPValues<T,NQPoints,Ndofs>& A,const FQPValues<T,NQPoints,Ndofs>& B)
  {
   for(Integer ii=0;ii<Ndofs;ii++)
    for(Integer jj=0;jj<NQPoints;jj++)
      Assignment<T>::apply(A[ii][jj],B[ii][jj]);
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


template<typename T,Integer Rows,Integer Cols>
class UnaryMinus<Matrix<T,Rows,Cols>> 
{    
 public:       
 inline static void apply(Matrix<T,Rows,Cols>& A,const Matrix<T,Rows,Cols>& B)
  {
   for(Integer ii=0;ii<Rows;ii++)
    for(Integer jj=0;jj<Cols;jj++)
           A(ii,jj)=-B(ii,jj);
  };
};


template<typename T,Integer Rows,Integer Cols>
class Addition<Matrix<T,Rows,Cols>,Matrix<T,Rows,Cols>> 
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
class Subtraction<Matrix<T,Rows,Cols>,Matrix<T,Rows,Cols>> 
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

template<Integer Rows1,Integer Cols1,Integer Rows2,Integer Cols2>
class Multiplication<Matrix<Real,Rows1,Cols1>,
                     Matrix<Real,Rows2,Cols2>> 
{    
 public:      
 using Output=OperatorType<Multiplication<Matrix<Real,Rows1,Cols1>,Matrix<Real,Rows2,Cols2>> >;
 using Left=Matrix<Real,Rows1,Cols1>;
 using Right=Matrix<Real,Rows2,Cols2>;

 template<Integer M1,Integer N1,Integer M2, Integer N2>
 class apply_aux
 {
  public:
  inline constexpr static void apply(Output& A,const Left& B,const Right& C)
  {
      std::cout<<" multiply "<<B<<std::endl;
  std::cout<<" multiply Integer M1,Integer N1,Integer M2, Integer N2"<<C<<std::endl;
   for(Integer ii=0;ii<Rows1;ii++)
    for(Integer jj=0;jj<Cols2;jj++)
       {
        A(ii,jj)=B(ii,0)*C(0,jj);
        for(Integer cc=1;cc<Cols1;cc++)
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
      std::cout<<" multiply "<<B<<std::endl;
  std::cout<<" multiply 1,1,M2,N2"<<C<<std::endl;
   for(Integer ii=0;ii<Rows2;ii++)
    for(Integer jj=0;jj<Cols2;jj++)
        A(ii,jj)=B(0,0)*C(ii,jj);
  };
  
 };

template<Integer M1, Integer N1>
 class apply_aux<M1,N1,1,1>
 {
  public:
  inline constexpr static void apply(Output& A,const Left& B,const Right& C)  
   {
      std::cout<<" multiply "<<B<<std::endl;
  std::cout<<" multiply M1,N1,1,1"<<C<<std::endl;
   for(Integer ii=0;ii<Rows1;ii++)
    for(Integer jj=0;jj<Cols1;jj++)
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

 inline constexpr static void apply(Output& A,const Left& B,const Right& C)
  {
    apply_aux<Rows1,Cols1,Rows2,Cols2>::apply(A,B,C);
  };


};


template<Integer Rows,Integer Cols>
class Division<Matrix<Real,Rows,Cols>,
               Matrix<Real,1,1>> 

{
  public:
 using Output=OperatorType<Multiplication<Matrix<Real,Rows,Cols>,Matrix<Real,1,1>> >;
 using Left=Matrix<Real,Rows,Cols>;
 using Right=Matrix<Real,1,1>;

  inline static void apply(Output& A,const Left& B,const Right& C)
  {
  // TODO FIXME IMPLEMENT HERE
    for(Integer ii=0;ii<Rows;ii++)
      for(Integer jj=0;jj<Cols;jj++)
        A(ii,jj)=B(ii,jj)/C(0,0);
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
 inline static void apply(      QPValues<U,NQPoints>& A,
                          const QPValues<U,NQPoints>& B)
  {
   for(Integer qp=0;qp<NQPoints;qp++)
        UnaryPlus<U>::apply(A[qp],B[qp]);
  };
};

template<typename U,Integer NQPoints>
class UnaryMinus<QPValues<U,NQPoints>> 
{       
 public:               
 inline static void apply(      QPValues<U,NQPoints>& A,
                          const QPValues<U,NQPoints>& B)
  {
   for(Integer qp=0;qp<NQPoints;qp++)
        UnaryMinus<U>::apply(A[qp],B[qp]);
  };
};



template<typename U,Integer NQPoints>
class Addition<QPValues<U,NQPoints>,
               QPValues<U,NQPoints>> 
{       
 public:               
 inline static void apply(      QPValues<U,NQPoints>& A,
                          const QPValues<U,NQPoints>& B,
                          const QPValues<U,NQPoints>& C)
  {
   for(Integer qp=0;qp<NQPoints;qp++)
        Addition<U,U>::apply(A[qp],B[qp],C[qp]);
  };

};

template<typename U,Integer NQPoints>
class Subtraction<QPValues<U,NQPoints>,
                  QPValues<U,NQPoints>> 
{       
 public:               
 inline static void apply(      QPValues<U,NQPoints>& A,
                          const QPValues<U,NQPoints>& B,
                          const QPValues<U,NQPoints>& C)
  {
   for(Integer qp=0;qp<NQPoints;qp++)
        Subtraction<U,U>::apply(A[qp],B[qp],C[qp]);
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
   for(Integer qp=0;qp<NQPoints;qp++)
        Multiplication<U,V>::apply(A[qp],B[qp],C[qp]);
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
        Division<U,V>::apply(A[qp],B[qp],C[qp]);
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
       {UnaryPlus<U>::apply(A(n_dof,qp),B(n_dof,qp));}
  };
};

template<typename U, Integer NQPoints,Integer Ndofs>
class UnaryMinus< FQPValues<U,NQPoints,Ndofs>> 
{       
 public:               
 inline static void apply(FQPValues<U,NQPoints,Ndofs>& A,const FQPValues<U,NQPoints,Ndofs>& B)
  {
  for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     for(Integer qp=0;qp<NQPoints;qp++)
       {UnaryMinus<U>::apply(A(n_dof,qp),B(n_dof,qp));}
  };
};

template<typename U, Integer NQPoints,Integer Ndofs>
class Addition< FQPValues<U,NQPoints,Ndofs>,
                FQPValues<U,NQPoints,Ndofs>> 
{       
 public:               
 inline static void apply(      FQPValues<U,NQPoints,Ndofs>& A,
                          const FQPValues<U,NQPoints,Ndofs>& B,
                          const FQPValues<U,NQPoints,Ndofs>& C)
  {
  for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
   for(Integer qp=0;qp<NQPoints;qp++)
    {Addition<U,U>::apply(A(n_dof,qp),B(n_dof,qp),C(n_dof,qp));}
  };


};

template<typename U, Integer NQPoints,Integer Ndofs>
class Subtraction< FQPValues<U,NQPoints,Ndofs>,
                   FQPValues<U,NQPoints,Ndofs>> 
{       
 public:               
 inline static void apply(      FQPValues<U,NQPoints,Ndofs>& A,
                          const FQPValues<U,NQPoints,Ndofs>& B,
                          const FQPValues<U,NQPoints,Ndofs>& C)
  {
  for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
   for(Integer qp=0;qp<NQPoints;qp++)
      {Subtraction<U,U>::apply(A(n_dof,qp),B(n_dof,qp),C(n_dof,qp));}
  };
};



























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


template<typename T>
constexpr void Plus(T& t,const T& constt){UnaryPlus<T>::apply(t,constt);}

template<typename T>
constexpr void Minus(T& t,const T& constt){UnaryMinus<T>::apply(t,constt);}

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






// template< typename DerivedLeft>
// class Multiplication< Expression <DerivedLeft>, Real > 
// operator*(const Expression<DerivedLeft>&left, const Real&right)
// {return Multiplication< Expression <DerivedLeft>, Real > (left,right);}

// template< typename DerivedRight>
// class Multiplication< Real, Expression <DerivedRight> >
// operator*(const Real&left, const Expression<DerivedRight>&right)
// {return Multiplication< Real, Expression <DerivedRight> >(left,right);}




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
  std::cout<<"Evaluation<Expression<UnaryPlus "<<output<<std::endl;
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
  eval_.apply(derived_value_,inputs...);

  Minus(output,derived_value_);
  std::cout<<"Evaluation<Expression<UnaryMinus "<<output<<std::endl;
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
  std::cout<<"Evaluation<Expression<UnaryMinus<UnaryMinus "<<output<<std::endl;
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
  Add(output,left_value_,right_value_);
  std::cout<<"Evaluation<Expression<Addition "<<output<<std::endl;

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
  std::cout<<"Evaluation<Expression<Subtraction "<<output<<std::endl;

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
{return Evaluation< Expression<GeneralForm<Form>>,ShapeFunctions2<GeneralForm_,GeneralForms...> >(form,shapes);}











































































































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



}

#endif //MARS_OPERATORS_HPP
