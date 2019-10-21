#ifndef MARS_MATRIX_STATIC_ASSERT_HPP
#define MARS_MATRIX_STATIC_ASSERT_HPP
#include "mars_base.hpp"
#include "mars_matrix.hpp"

namespace mars {

template<typename...Ts>
class StaticAssertMatrixOrTransposed;

template<typename T>
class StaticAssertMatrixOrTransposed<T>
{
 public:
  static constexpr Integer Rows=T::Rows;
  static constexpr Integer Cols=T::Cols;

  static constexpr bool value=
   IsSame<T,           Matrix<Real,Rows,Cols> >::value ||
   IsSame<T,Transposed<Matrix<Real,Cols,Rows>>>::value;
};




template<typename Left,typename Right>
class StaticAssertMatrixOrTransposed<Left,Right>
{
 public:
  static constexpr Integer RowsLeft=Left::Rows;
  static constexpr Integer ColsLeft=Left::Cols;
  static constexpr Integer RowsRight=Right::Rows;
  static constexpr Integer ColsRight=Right::Cols;

  static constexpr bool value=
  (IsSame<Left,Matrix<Real,RowsLeft,ColsLeft>>::value && IsSame<Right,Matrix<Real,RowsRight,ColsRight>>::value) ||
  (IsSame<Left,Matrix<Real,RowsLeft,ColsLeft>>::value && IsSame<Right,Transposed<Matrix<Real,ColsRight,RowsRight>>>::value) ||
  (IsSame<Left,Transposed<Matrix<Real,ColsLeft,RowsLeft>>>::value && IsSame<Right,Matrix<Real,RowsRight,ColsRight>>::value) ||
  (IsSame<Left,Transposed<Matrix<Real,ColsLeft,RowsLeft>>>::value && IsSame<Right,Transposed<Matrix<Real,ColsRight,RowsRight>>>::value);
};


}
#endif