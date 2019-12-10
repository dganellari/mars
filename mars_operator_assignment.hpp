#ifndef MARS_OPERATOR_ASSIGNMENT_HPP
#define MARS_OPERATOR_ASSIGNMENT_HPP
#include "mars_base.hpp"
#include "mars_matrix.hpp"
#include "mars_operator_apply.hpp"

namespace mars {

template<typename T>
class Assignment;

template<typename S,typename T>
constexpr void Assign(S& t,const T& constt){Assignment<T>::apply(t,constt);}


template<>
class Assignment<Real> 
{
 public:       
 inline static constexpr void apply(Real& A,const Real& B)
  {
           A=B;
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
    OperatorApply<Assignment<QPValues<T,NQPoints>>>::apply(A,B);

    // T ok1(1);
    // S ok2(2);
   // for(Integer ii=0;ii<NQPoints;ii++)
   //    {

   //      Assign(A[ii],B[ii]);

        // Assignment<S>::apply(A[ii],B[ii]);

        // Assignment<T>::apply(A[ii],B[ii]);
        // std::cout<<"Assignment<QPValues<T,NQPoints>> "<<A<<" "<<B<<std::endl;
      // }
  };







};

template<typename T,Integer NQPoints,Integer Ndofs>
class Assignment<FQPValues<T,NQPoints,Ndofs>> 
{
 public:       
 template<typename S>
 inline static constexpr void apply(FQPValues<S,NQPoints,Ndofs>& A,const FQPValues<T,NQPoints,Ndofs>& B)
  {
    OperatorApply<Assignment<FQPValues<T,NQPoints,Ndofs>>>::apply(A,B);
   // for(Integer ii=0;ii<Ndofs;ii++)
   //  for(Integer jj=0;jj<NQPoints;jj++)
   //    Assign(A[ii][jj],B[ii][jj]);
  };
};


}
#endif