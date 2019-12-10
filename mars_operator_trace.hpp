#ifndef MARS_OPERATOR_TRACE_HPP
#define MARS_OPERATOR_TRACE_HPP
#include "mars_base.hpp"
#include "mars_matrix.hpp"

namespace mars {

// template<typename T,Integer Rows,Integer Cols>
// class MatTraceOperator<Matrix<T,Rows,Cols>> 
// {    
//  public:     
 
//  template<typename Output>
//  inline static void apply(Output& A,const Matrix<T,Rows,Cols>& B)
//   {
//     static_assert(Rows==Cols && "trace operator can only be applied to square matrices");
//     std::cout<<"MatTraceOperator<Matrix<T,Rows,Cols>> "<<std::endl;
//    A(0,0)=B(0,0);
//    for(Integer ii=1;ii<Rows;ii++)
//     A(0,0)+=B(ii,ii);
//   };


// };

// template<typename T,Integer Rows,Integer Cols>
// class MatTraceOperator<Transposed<Matrix<T,Rows,Cols>>> 
// {    
//  public:     
 
//  template<typename Output>
//  inline static void apply(Output& A,const Transposed<Matrix<T,Rows,Cols>>& B)
//   {
//     static_assert(Rows==Cols && "trace operator can only be applied to square matrices");
//     std::cout<<"MatTraceOperator<Transposed<Matrix<T,Rows,Cols>> "<<std::endl;
//    A(0,0)=B(0,0);
//    for(Integer ii=1;ii<Rows;ii++)
//     A(0,0)+=B(ii,ii);
//   };


// };








template<typename T>
class MatTraceOperator
{    
 public:     
  static_assert(StaticAssertMatrixOrTransposed<T>::value && "UnaryMinus allows:mat+mat;");

 template<typename Output>
 inline static void apply(Output& A,const T& B)
  {
    // std::cout<<"pre MatTraceOperator"<<std::endl;
    static_assert(T::Rows==T::Cols && "trace operator can only be applied to square matrices");
    // std::cout<<"MatTraceOperator<Transposed<Matrix<T,Rows,Cols>> "<<std::endl;
   A(0,0)=B(0,0);
   for(Integer ii=1;ii<T::Rows;ii++)
    A(0,0)+=B(ii,ii);
    // std::cout<<"after MatTraceOperator"<<std::endl;
  };


};










// template<typename U,Integer NQPoints>
// class MatTraceOperator<QPValues<U,NQPoints>> 
// {       
//  public:               
//  template<typename V>
//  inline static void apply(      QPValues<V,NQPoints>& A,
//                           const QPValues<U,NQPoints>& B)
//   {
//    for(Integer qp=0;qp<NQPoints;qp++)
//         {
//           std::cout<<"pre MatTraceOperator QPValues"<<std::endl;
//           MatTraceOperator<U>::apply(A[qp],B[qp]);
//         }
//   };
// };








// template<typename U, Integer NQPoints,Integer Ndofs>
// class MatTraceOperator<FQPValues<U,NQPoints,Ndofs>> 
// {       
//  public:               
//  template<typename V>
//  inline static void apply(FQPValues<V,NQPoints,Ndofs>& A,const FQPValues<U,NQPoints,Ndofs>& B)
//   {
//   for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
//      for(Integer qp=0;qp<NQPoints;qp++)
//        {
//         MatTraceOperator<U>::apply(A(n_dof,qp),B(n_dof,qp));
//        }
//   };
// };


template<typename U,typename V>
constexpr void MatTrace(U& t,const V& constt){MatTraceOperator<V>::apply(t,constt);}

template<typename T>
constexpr auto MatTrace(const Expression<T>& t){return MatTraceOperator<Expression<T>>(t);}

template<typename T>
constexpr auto MatTrace(const Expression<MatTraceOperator<Expression<T>>>& t){return T(t.derived().derived());}


template<typename T>
class MatTraceOperator<Expression<T>>: public Expression<MatTraceOperator<Expression<T>>>
{
public:

    MatTraceOperator(const Expression<T>& expr): value_(expr.derived()){};
    MatTraceOperator(const Expression<MatTraceOperator<Expression<T>>>& expr): value_(expr.derived().derived().derived()){};

    T& operator()(){return value_;};
    const constexpr T& operator()()const{return value_;};
    T& derived(){return value_;};
    const constexpr T& derived()const{return value_;};


  private:
  T value_;
};


}
#endif