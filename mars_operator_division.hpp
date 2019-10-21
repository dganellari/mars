#ifndef MARS_OPERATOR_DIVISION_HPP
#define MARS_OPERATOR_DIVISION_HPP
#include "mars_base.hpp"
#include "mars_matrix.hpp"
#include "mars_operator_apply.hpp"


namespace mars {

template<typename Left,typename Right>
class Division;

template<typename U, typename V,typename W>
constexpr void Divide(U& u, const V& v, const W& w){Division<V,W>::apply(u,v,w);}

template< typename DerivedLeft,typename DerivedRight>
class Division< Expression <DerivedLeft>, Expression <DerivedRight> >
operator/(const Expression<DerivedLeft>&left, const Expression<DerivedRight>&right)
{return Division< Expression <DerivedLeft>, Expression<DerivedRight> > (left,right);}



template<typename Left,typename Right>
class Division
{
 public: 
 static_assert(StaticAssertMatrixOrTransposed<Left,Right>::value && "Division allows:mat+mat;");
 template<typename Output>
 inline static void apply(      Output& A,
                          const Left& B,
                          const Right& C)
  {
        // DivisionMatrixAndTransposedAux<Left,Right>::apply(A,B,C);
   for(Integer ii=0;ii<Left::Rows;ii++)
    for(Integer jj=0;jj<Left::Cols;jj++)
           A(ii,jj)=B(ii,jj)/C(0,0);
  };
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


}
#endif

