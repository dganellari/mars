#ifndef MARS_OPERATOR_UNARY_PLUS_HPP
#define MARS_OPERATOR_UNARY_PLUS_HPP
#include "mars_base.hpp"
#include "mars_matrix.hpp"

namespace mars {
template<typename T>
class UnaryPlus;


template<typename T>
class UnaryPlus
{    
 public:     
  static_assert(StaticAssertMatrixOrTransposed<T>::value && "Addition allows:mat+mat;");

 template<typename Output>
 inline static void apply(Output& A,const T& B)
  {
    // std::cout<<"pre UnaryPlusMatrixAndTransposedAux<Transposed"<<std::endl;
   for(Integer ii=0;ii<T::Rows;ii++)
    for(Integer jj=0;jj<T::Cols;jj++)
           A(ii,jj)=B(ii,jj);
    // std::cout<<A<<B<<std::endl;
    // std::cout<<"after UnaryPlusMatrixAndTransposedAux<Transposed"<<std::endl;
  };


};



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
class UnaryPlus< Expression <Derived> >
operator+(const Expression<Derived>& expr)
{return UnaryPlus< Expression <Derived> >(expr);}


}
#endif