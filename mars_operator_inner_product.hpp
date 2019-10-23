#ifndef MARS_OPERATOR_INNER_PRODUCT_HPP
#define MARS_OPERATOR_INNER_PRODUCT_HPP
#include "mars_base.hpp"
#include "mars_matrix.hpp"
#include "mars_operator_apply.hpp"

namespace mars {


template<typename Left,typename Right>
class InnerProduct;


template< typename DerivedLeft,typename DerivedRight>
class InnerProduct< Expression <DerivedLeft>, Expression <DerivedRight> >
Inner(const Expression<DerivedLeft>&left, const Expression<DerivedRight>&right)
{return InnerProduct< Expression <DerivedLeft>, Expression <DerivedRight> >(left,right);}


template<typename Left,typename Right>
class InnerProduct
{
 public: 
 static_assert(StaticAssertMatrixOrTransposed<Left,Right>::value );

 template<typename Output>
 inline static void apply(      Output& A,
                          const Left& B,
                          const Right& C)
  {
   A(0,0)=B(0,0)*C(0,0);
   for(Integer jj=1;jj<Left::Cols;jj++)
    A(0,0)+=B(0,jj)*C(0,jj);

   for(Integer ii=1;ii<Left::Rows;ii++)
    for(Integer jj=0;jj<Left::Cols;jj++)
           A(0,0)+=B(ii,jj)*C(ii,jj);
    std::cout<<"after InnerProduct<Transposed"<<std::endl;

        // AdditionMatrixAndTransposedAux<Left,Right>::apply(A,B,C);
  };
};


template< typename DerivedLeft,typename DerivedRight>
class InnerProduct< Expression <DerivedLeft>, Expression <DerivedRight> > 
: public Expression< InnerProduct< Expression <DerivedLeft>, Expression <DerivedRight> > >
{
  public:
    using Left=DerivedLeft;
    using Right=DerivedRight;

    InnerProduct(const Expression<DerivedLeft>& left, const Expression<DerivedRight>&right)
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