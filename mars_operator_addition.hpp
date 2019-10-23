#ifndef MARS_OPERATOR_ADDITION_HPP
#define MARS_OPERATOR_ADDITION_HPP
#include "mars_base.hpp"
#include "mars_matrix.hpp"
#include "mars_operator_apply.hpp"


namespace mars {

template<typename Left,typename Right>
class Addition;

template<typename U,typename V,typename W>
constexpr void Add(U& u, const V& v, const W& w){Addition<V,W>::apply(u,v,w);}

template< typename DerivedLeft,typename DerivedRight>
class Addition< Expression <DerivedLeft>, Expression <DerivedRight> >
operator+(const Expression<DerivedLeft>&left, const Expression<DerivedRight>&right)
{return Addition< Expression <DerivedLeft>, Expression <DerivedRight> >(left,right);}



template<typename Left,typename Right>
class Addition
{
 public: 
 static_assert(StaticAssertMatrixOrTransposed<Left,Right>::value && "Addition allows:mat+mat;");
 template<typename Output>
 inline static void apply(      Output& A,
                          const Left& B,
                          const Right& C)
  {
   for(Integer ii=0;ii<Left::Rows;ii++)
    for(Integer jj=0;jj<Left::Cols;jj++)
           A(ii,jj)=B(ii,jj)+C(ii,jj);
    std::cout<<"after AdditionMatrixAndTransposedAux<Transposed"<<std::endl;

        // AdditionMatrixAndTransposedAux<Left,Right>::apply(A,B,C);
  };
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

}
#endif