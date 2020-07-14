#ifndef MARS_OPERATOR_MULTIPLICATION_HPP
#define MARS_OPERATOR_MULTIPLICATION_HPP
#include "mars_base.hpp"
#include "mars_matrix.hpp"
#include "mars_operator_apply.hpp"


namespace mars {

template<typename Left,typename Right>
class Multiplication;



template<typename U,typename V,typename W>
constexpr void Multiply(U& u, const V& v, const W& w){Multiplication<V,W>::apply(u,v,w);}


template< typename DerivedLeft,typename DerivedRight>
class Multiplication< Expression <DerivedLeft>, Expression <DerivedRight> >
operator*(const Expression<DerivedLeft>&left, const Expression<DerivedRight>&right)
{return Multiplication< Expression <DerivedLeft>, Expression <DerivedRight> >(left,right);}




template<typename Left, typename Right>
class MultiplicationMatrixAndTransposedAux
{
public:
 using Output=OperatorType<Multiplication<Left,Right> >;

 template<Integer M1,Integer N1,Integer M2, Integer N2,Integer Dummy>
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

template<Integer Dummy>
 class apply_aux<1,1,1,1,Dummy>
 {
 public:
  inline constexpr static void apply(Output& A,const Left& B,const Right& C)  
   {A(0,0)=B(0,0)*C(0,0);};
  
 };
 
template<Integer M2, Integer N2,Integer Dummy>
 class apply_aux<1,1,M2,N2,Dummy>
 {
  public:
  inline constexpr static void apply(Output& A,const Left& B,const Right& C)  
   {
   for(Integer ii=0;ii<Right::Rows;ii++)
    for(Integer jj=0;jj<Right::Cols;jj++)
        A(ii,jj)=B(0,0)*C(ii,jj);
  };
  
 };

template<Integer M1, Integer N1,Integer Dummy>
 class apply_aux<M1,N1,1,1,Dummy>
 {
  public:
  inline constexpr static void apply(Output& A,const Left& B,const Right& C)  
   {
   for(Integer ii=0;ii<Left::Rows;ii++)
    for(Integer jj=0;jj<Left::Cols;jj++)
        A(ii,jj)=B(ii,jj)*C(0,0);
  };
  
 };


 inline constexpr static void apply(Output& A,
                                    const Left& B,
                                    const Right& C)
  {
    apply_aux<Left::Rows,Left::Cols,Right::Rows,Right::Cols,1>::apply(A,B,C);
  };

};




template<typename Left,typename Right>
class Multiplication
{
 public: 
 static_assert(StaticAssertMatrixOrTransposed<Left,Right>::value && "Multiplication allows:mat+mat;");
 template<typename Output>
 inline static void apply(      Output& A,
                          const Left& B,
                          const Right& C)
  {
        MultiplicationMatrixAndTransposedAux<Left,Right>::apply(A,B,C);
  };
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

    const constexpr DerivedLeft& left()const{return left_;};
    const constexpr DerivedRight& right()const{return right_;};
  private:
  DerivedLeft left_;
  DerivedRight right_;
};


}
#endif

