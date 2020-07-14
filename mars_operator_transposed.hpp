#ifndef MARS_OPERATOR_TRANSPOSED_HPP
#define MARS_OPERATOR_TRANSPOSED_HPP
#include "mars_base.hpp"
#include "mars_matrix.hpp"

namespace mars {

  template<typename T>
  class Transposed;

  template<typename T>
  constexpr auto Transpose(const Expression<T>& t){return Transposed<Expression<T>>(t);}

  template<typename T>
  constexpr auto Transpose(const Expression<Transposed<Expression<T>>>& t){return T(t.derived().derived());}


  template<typename T, Integer TransposedRows, Integer TransposedCols>
  class Transposed<Matrix<T,TransposedRows,TransposedCols>>
  : public Matrix<T,TransposedRows,TransposedCols>
  {
  public:
    static constexpr Integer Rows=TransposedCols;
    static constexpr Integer Cols=TransposedRows;
    using type= Matrix<T,Rows,Cols>;
    using subtype=T;
    using MB = TensorBase<T, std::make_index_sequence<Rows*Cols>>;


    constexpr Transposed()
    {}

    constexpr Transposed(const Matrix<T,TransposedRows,TransposedCols>& mat)
    :
    mat_ptr_(std::make_shared<Matrix<T,TransposedRows,TransposedCols>>(mat))
    {}

    inline constexpr const auto& operator()() const
    {
     return *mat_ptr_;
    }


    inline constexpr void operator()(const Matrix<T,TransposedRows,TransposedCols>& mat)
    {
     mat_ptr_=std::make_shared<Matrix<T,TransposedRows,TransposedCols>>(mat);
    }

    inline constexpr      T &operator()(const Integer& i, const Integer& j)
    {
      assert(i < Rows);
      assert(j < Cols);
      return (*mat_ptr_)(j,i);
    }

    inline constexpr const T &operator()(const Integer i, const Integer j) const
    {
      assert(i < Rows);
      assert(j < Cols);
      return (*mat_ptr_)(j,i);
    }

    void describe(std::ostream &os) const
    {
      for(Integer i = 0; i < Rows; ++i) {
        for(Integer j = 0; j < Cols; ++j) {
          os << (*this)(i, j) << " ";
        }
        os << "\n";
      }

      os << "\n";
    }

    friend std::ostream &operator<<(std::ostream &os, const Transposed &m)
    {
      m.describe(os);
      return os;
    }
     inline static constexpr void apply(Transposed<Matrix<T,Cols,Rows>>& A,const Matrix<T,Rows,Cols>& B)
      {
           A(B);
      };

     inline static constexpr void apply(Transposed<Matrix<T,Cols,Rows>>& A,const Transposed<Matrix<T,Cols,Rows>>& B)
      {
           A(B());
      };
    
  private:
    std::shared_ptr<Matrix<T,TransposedRows,TransposedCols>> mat_ptr_;

  };



template<typename T>
class Transposed<Expression<T>>: public Expression<Transposed<Expression<T>>>
{
public:

    Transposed(const Expression<T>& expr): value_(expr.derived()){};
    Transposed(const Expression<Transposed<Expression<T>>>& expr): value_(expr.derived().derived().derived()){};

    T& operator()(){return value_;};
    const constexpr T& operator()()const{return value_;};
    T& derived(){return value_;};
    const constexpr T& derived()const{return value_;};


  private:
  T value_;
};



}
#endif