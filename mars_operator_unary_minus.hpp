#ifndef MARS_OPERATOR_UNARY_MINUS_HPP
#define MARS_OPERATOR_UNARY_MINUS_HPP
#include "mars_base.hpp"
#include "mars_matrix.hpp"

namespace mars {
template<typename T>
class UnaryMinus;

template<typename U,typename V>
constexpr void Minus(U& t,const V& constt){UnaryMinus<V>::apply(t,constt);}

template< typename Derived>
class UnaryMinus< Expression <Derived> >
operator-(const Expression<Derived>& expr)
{return UnaryMinus< Expression <Derived> >(expr);}

template<typename T>
class UnaryMinus
{    
 public:     
  static_assert(StaticAssertMatrixOrTransposed<T>::value && "UnaryMinus allows:mat+mat;");

 template<typename Output>
 inline static void apply(Output& A,const T& B)
  {
    std::cout<<"pre UnaryMinus"<<std::endl;
   for(Integer ii=0;ii<T::Rows;ii++)
    for(Integer jj=0;jj<T::Cols;jj++)
           A(ii,jj)=-B(ii,jj);
    std::cout<<A<<B<<std::endl;
    std::cout<<"after UnaryMinus"<<std::endl;
  };


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
}
#endif