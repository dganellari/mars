#ifndef MARS_FQP_EXPRESSION_HPP
#define MARS_FQP_EXPRESSION_HPP

#include "mars_simplex.hpp"
#include "mars_base_elementfunctionspace.hpp"
#include "mars_vector.hpp"
#include "mars_operators.hpp"

namespace mars{



template<typename Derived, typename Type>
class AlgebraicExpression: public Expression<Derived>
{
public:
      using Base=Expression<Derived>;
      using Base::derived;
      using subtype=typename Type::subtype;
      constexpr Type & value();
      constexpr const Type & value()const;
      
      inline constexpr subtype &operator[](const Integer i) {return derived().value()[i];}
      inline constexpr const subtype &operator[](const Integer i)const {return derived().value()[i];}

      // // EQUAL
      // inline constexpr Derived& operator = (const Derived &u)
      // {            
      //   value()=u.value();
      //   return *this;
      // } 

      // // BINARY ADD
      // inline constexpr Derived operator+(const Derived &other) const
      // {
      //   Derived result;
      //   result.value()=derived().value()+other.value();
      //   return result;
      // };

      // // UNARY ADD
      // inline constexpr Derived operator+() const
      // {
      //    return derived();
      // };

      // // BINARY MINUS
      // inline constexpr Derived operator-(const Derived &other) const
      // {
      //   Derived result;
      //   result.value()=derived().value()-other.value();
      //   return result;
      // };
      // // UNARY MINUS
      // inline constexpr Derived operator-() const
      // {
      //   Derived result;
      //   result.value()=-derived().value();
      //   return result;
      // };

    
      // // LEFT SCALAR MULTIPLY
      // friend constexpr Derived operator *(const Real &alpha, const Derived& der) 
      // { 
      //   Derived result;
      //   result.value()=alpha*der.value();
      //   return result;
      // };

      // // RIGHT SCALAR MULTIPLY
      // friend constexpr Derived operator *(const Derived& der,const Real &alpha) 
      // { 
      //   Derived result;
      //   result.value()=alpha*der.value();
      //   return result;
      // };

      // // RIGHT SCALAR DIVISION
      // friend constexpr Derived operator /(const Derived& der,const Real &alpha) 
      // { 
      //   Derived result;
      //   result.value()=der.value()/alpha;
      //   return result;
      // };

      // // RIGHT SCALAR EQUAL DIVISION
      // inline constexpr Derived& operator /=(const Real &alpha) 
      // { 
      //   (*this).derived().value()=(*this).derived().value()/alpha;
      //   return (*this).derived();
      // };
      // // RIGHT SCALAR EQUAL MULTIPLICATION
      // inline constexpr Derived& operator *=(const Real &alpha) 
      // { 
      //   (*this).derived().value()=(*this).derived().value()*alpha;
      //   return (*this).derived();
      // };   
      // // EQUAL ADDITION
      // inline constexpr Derived& operator +=(const Derived &other) 
      // { 
      //   (*this).derived().value()=(*this).derived().value()+other.value();
      //   return (*this).derived();
      // };  
      // // EQUAL SUBTRACTION
      // inline constexpr Derived& operator -=(const Derived &other) 
      // { 
      //   (*this).derived().value()=(*this).derived().value()-other.value();
      //   return (*this).derived();
      // };


      // PRINTING 
      void describe(std::ostream &os) const
      {
          os << (*this).derived().value()<<" ";
          os << "\n";
      }

      friend std::ostream &operator<<(std::ostream &os, const Derived &v)
      {
          v.describe(os);
          return os;
      }

};


// QPVALUES 
template<typename T,Integer NQPoints>
class QPValues: 
public AlgebraicExpression<QPValues<T,NQPoints>,Vector<T,NQPoints>>,
public TensorBase<T, std::make_index_sequence<NQPoints>> 
{
 public:
    using MB = TensorBase<T, std::make_index_sequence<NQPoints>>;
    using MB::MB;
    using MB::values;
    using type=  Vector<T,NQPoints>;
    inline constexpr       type& operator()()      {return values_;};
    inline constexpr const type& operator()()const {return values_;};
    inline constexpr       type& value()      {return values_;};
    inline constexpr const type& value()const {return values_;};
    inline constexpr Integer size()const{return NQPoints;}
  // QPValues(){};
  // QPValues(const type& t):values_(t) {};

protected:
      type values_;
};

// FQPVALUES 
template<typename T, Integer NQPoints,Integer NComponents>
class FQPValues: public 
AlgebraicExpression<FQPValues<T,NQPoints,NComponents>,
                    Vector<Vector<T,NQPoints>,NComponents>>
{
public:
      using type= Vector<Vector<T,NQPoints>,NComponents>;
      using subttype= Vector<T,NQPoints>;
      ~FQPValues() = default;
      constexpr FQPValues()= default;
      constexpr FQPValues(const type& v): values_(v){};
      constexpr FQPValues (FQPValues const &) = default;
      constexpr FQPValues (FQPValues &&) = default;
      constexpr FQPValues & operator= (FQPValues const &) = default;
      constexpr FQPValues & operator= (FQPValues &&) = default;
      constexpr       type & value()     {return values_;};
      constexpr const type & value()const{return values_;};
      inline constexpr const type operator()()const{return values_;};
      inline constexpr       type operator()()     {return values_;};
      inline constexpr const Vector<T,NQPoints> operator()(const Integer i)const{return values_[i];};
      inline constexpr       Vector<T,NQPoints> operator()(const Integer i)     {return values_[i];};
      inline constexpr Integer size()const{return NComponents;}
      // inline constexpr void operator()(const Integer& i,const Vector<T,NQPoints>& u){values_[i]=u;};         
protected:
  type values_;
};






// QPVALUES = QPVALUES * QPVALUES

      // template< typename Left, typename Right,Integer NQPoints>
      // inline QPValues<typename OperatorType< Multiplication<Left,Right>>::type,NQPoints>  operator*
      // (const QPValues<Left,NQPoints>&left, const QPValues<Right,NQPoints>&right)
      // {
      //   using S=typename OperatorType< Multiplication<Left,Right>>::type;
      //   QPValues<S,NQPoints> result;
      //   for(Integer qp=0;qp<NQPoints;qp++)
      //   {
      //     result[qp]= left[qp]*right[qp];
      //   };
      // return result;
      // }

// FQPVALUES = QPVALUES * FQPVALUES

      // template< typename Left, typename Right,Integer NComponents,Integer NQPoints>
      // inline FQPValues<typename OperatorType< Multiplication<Left,Right>>::type,NQPoints,NComponents>   operator*
      // (const QPValues<Left,NQPoints>&left, const FQPValues<Right,NQPoints,NComponents>&right)
      // {
      //   using S=typename OperatorType< Multiplication<Left,Right>>::type;
      //   FQPValues<S,NQPoints,NComponents> result;
      //   for(Integer n_comp=0;n_comp<NComponents;n_comp++)
      //       for(Integer qp=0;qp<NQPoints;qp++)
      //       {
      //         result[n_comp][qp]= left[qp]*right[n_comp][qp];
      //       };
      // return result;
      // }

// FQPVALUES = FQPVALUES * FQPVALUES

      // template< typename Trial, typename Test,Integer NComponentsTrial,Integer NComponentsTest, Integer NQPoints>
      // inline Matrix<Vector<Real,NQPoints>,NComponentsTest,NComponentsTrial>   operator*
      // (const FQPValues<Trial,NQPoints,NComponentsTrial>&trial, 
      //  const FQPValues<Test, NQPoints,NComponentsTest> &test)
      // {
      //   Matrix<Vector<Real,NQPoints>,NComponentsTest,NComponentsTrial> result;
      //   Contraction contract;
      //   for(Integer n_comp_test=0;n_comp_test<NComponentsTest;n_comp_test++)
      //       for(Integer n_comp_trial=0;n_comp_trial<NComponentsTrial;n_comp_trial++)
      //           for(Integer qp=0;qp<NQPoints;qp++)
      //               result(n_comp_test,n_comp_trial)[qp]=contract(trial[n_comp_trial][qp],test[n_comp_test][qp]);
      // return result;
      // }


// template<typename Derived, typename S,Integer NQPoints, Integer Dim>
// class QPExpression: public Expression<Derived>
// {
//      public:
//      using Base=Expression<Derived>;
//      using Base::derived;
//      using Point = Vector<Real,Dim> ;
//      using type= QPValues < S, NQPoints>;
//      S& apply(const Point& point);

//       inline QPValues<S,NQPoints>& eval(const Matrix<Real,NQPoints,Dim> & qp_points)
//       {
//        for(Integer qp=0;qp<NQPoints;qp++)
//           {
//             qp_points.get_row(qp,row_);
//             value_[qp]=derived().apply(row_);
//           }

//         return value_;
//       };
       
//       void eval(const Matrix<Real,NQPoints,Dim> & qp_points, QPValues<S,NQPoints>& value_)
//       {
//        for(Integer qp=0;qp<NQPoints;qp++)
//           {
//             qp_points.get_row(qp,row_);
//             value_[qp]=derived().apply(row_);
//           }
//       };


//     protected:
//       QPValues<S,NQPoints> value_; 
//       Point row_; 
// };

// template<typename Derived, typename S,Integer NQPoints, Integer NComponents, Integer Dim>
// class FQPExpression:
// public Expression<Derived>
// {
//      public:
//      using Point = Vector<Real,Dim> ;
//      using type= FQPValues < S, NQPoints,NComponents>;
//      using Base=Expression<Derived>;
//      using Base::derived;
//      FQPValues<S,NQPoints,NComponents>& apply(const Matrix<Real,NQPoints,Dim> & qp_points);

//       inline FQPValues<S,NQPoints,NComponents>& eval(const Matrix<Real,NQPoints,Dim> & qp_points)
//       { 
//         value_=derived().apply(qp_points);
//         return value_;
//       };

//       void eval(const Matrix<Real,NQPoints,Dim> & qp_points, FQPValues<S,NQPoints,NComponents> value_)
//       {
//           value_ = derived().apply(qp_points);
//       };
//     protected:
//       FQPValues<S,NQPoints,NComponents> value_; 
//       Point row_; 
// };








//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////           UNARY PLUS: +QP        //////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// template< typename Derived, typename  T, Integer NQPoints,Integer Dim>
// class UnaryPlus< QPExpression <Derived, T, NQPoints,Dim>>
// : public QPExpression<UnaryPlus< QPExpression <Derived, T, NQPoints,Dim>>,
//                       typename OperatorType<UnaryPlus<typename T::type>>::type,
//                       NQPoints, Dim >
// {
//   public:
//     using Input= QPExpression<Derived,T,NQPoints,Dim>;
//     using ResultType=typename OperatorType< UnaryPlus<typename Input::type>>::type;
//     UnaryPlus(const Input& input): 
//     derived_(input.derived())
//     {};

//     template<typename QP>
//     ResultType  eval(const QP & qp_points)
//     {return +derived_.eval(qp_points);};

//     template<typename QP>
//     inline void  eval(const QP & qp_points, ResultType& result_)
//     {result_= +derived_.eval(qp_points);};

//   private:
//   Derived derived_;
//   ResultType result_;
// };


// template< typename Derived, typename  T, Integer NQPoints,Integer Dim>
// class UnaryPlus< QPExpression <Derived, T, NQPoints,Dim> >
// operator+(const QPExpression<Derived, T, NQPoints,Dim>& input)
//           {return UnaryPlus<QPExpression <Derived, T, NQPoints,Dim>>(input);}



//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////          UNARY MINUS: -QP        //////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// template< typename Derived, typename  T, Integer NQPoints,Integer Dim>
// class UnaryMinus< QPExpression <Derived, T, NQPoints,Dim>>
// : public QPExpression<UnaryMinus< QPExpression <Derived, T, NQPoints,Dim>>,
//                       typename OperatorType<UnaryMinus<typename T::type>>::type,
//                       NQPoints, Dim >
// {
//   public:
//     using Input= QPExpression<Derived,T,NQPoints,Dim>;
//     UnaryMinus(const Input& input): 
//     derived_(input.derived())
//     {};

//     template<typename QP>
//     typename OperatorType< UnaryMinus<typename Input::type>>::type  eval(const QP & qp_points)
//     {return -derived_.eval(qp_points);};

//   private:
//   Derived derived_;
// };


// template< typename Derived, typename  T, Integer NQPoints,Integer Dim>
// class UnaryMinus< QPExpression <Derived, T, NQPoints,Dim> >
// operator-(const QPExpression<Derived, T, NQPoints,Dim>& input)
//           {return UnaryMinus<QPExpression <Derived, T, NQPoints,Dim>>(input);}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////          ADDITION: QP + QP        /////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// template< typename DerivedLeft,typename DerivedRight, typename  TLeft,typename TRight, Integer NQPoints,Integer Dim>
// class Addition< QPExpression <DerivedLeft, TLeft, NQPoints,Dim>, 
//                      QPExpression <DerivedRight,TRight,NQPoints, Dim> >

// : public QPExpression<Addition<QPExpression <DerivedLeft, TLeft, NQPoints,Dim>,
//                                QPExpression <DerivedRight,TRight,NQPoints,Dim> >,
//                 typename OperatorType<Addition<typename TLeft::type,typename TRight::type>>::type,
//                 NQPoints, Dim >
// {
//   public:
//     using Left= QPExpression<DerivedLeft, TLeft, NQPoints,Dim>;
//     using Right=QPExpression<DerivedRight,TRight,NQPoints,Dim>;
//     Addition(const Left& left, const Right&right): 
//     left_(left.derived()),
//     right_(right.derived())
//     {};

//     template<typename QP>
//     typename OperatorType< Addition<typename Left::type, typename Right::type>>::type  eval(const QP & qp_points)
//     {return left_.eval(qp_points)+right_.eval(qp_points);};

//   private:
//   DerivedLeft left_;
//   DerivedRight right_;
// };


// template< typename DerivedLeft,typename DerivedRight, typename  TLeft,typename TRight, Integer NQPoints,Integer Dim>
// class Addition< QPExpression<DerivedLeft, TLeft, NQPoints,Dim>, 
//                      QPExpression<DerivedRight,TRight,NQPoints,Dim> >
// operator+(const QPExpression<DerivedLeft, TLeft, NQPoints,Dim>& left, 
//           const QPExpression<DerivedRight,TRight,NQPoints, Dim>&right)
// {return Addition<QPExpression<DerivedLeft, TLeft, NQPoints,Dim>, 
//                       QPExpression<DerivedRight,TRight,NQPoints,Dim> >(left,right);}


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////          SUBTRACTION: QP - QP        /////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// template< typename DerivedLeft,typename DerivedRight, typename  TLeft,typename TRight, Integer NQPoints,Integer Dim>
// class Subtraction< QPExpression <DerivedLeft, TLeft, NQPoints,Dim>, 
//                      QPExpression <DerivedRight,TRight,NQPoints, Dim> >

// : public QPExpression<Subtraction<QPExpression <DerivedLeft, TLeft, NQPoints,Dim>,
//                                QPExpression <DerivedRight,TRight,NQPoints,Dim> >,
//                 typename OperatorType<Subtraction<typename TLeft::type,typename TRight::type>>::type,
//                 NQPoints, Dim >
// {
//   public:
//     using Left= QPExpression<DerivedLeft, TLeft, NQPoints,Dim>;
//     using Right=QPExpression<DerivedRight,TRight,NQPoints,Dim>;
//     Subtraction(const Left& left, const Right&right): 
//     left_(left.derived()),
//     right_(right.derived())
//     {};

//     template<typename QP>
//     typename OperatorType< Subtraction<typename Left::type, typename Right::type>>::type  eval(const QP & qp_points)
//     {return left_.eval(qp_points)-right_.eval(qp_points);};

//   private:
//   DerivedLeft left_;
//   DerivedRight right_;
// };


// template< typename DerivedLeft,typename DerivedRight, typename  TLeft,typename TRight, Integer NQPoints,Integer Dim>
// class Subtraction< QPExpression<DerivedLeft, TLeft, NQPoints,Dim>, 
//                      QPExpression<DerivedRight,TRight,NQPoints,Dim> >
// operator-(const QPExpression<DerivedLeft, TLeft, NQPoints,Dim>& left, 
//           const QPExpression<DerivedRight,TRight,NQPoints, Dim>&right)
// {return Subtraction<QPExpression<DerivedLeft, TLeft, NQPoints,Dim>, 
//                       QPExpression<DerivedRight,TRight,NQPoints,Dim> >(left,right);}



//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////          DIVIDE: QP / REAL        /////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// template< typename DerivedLeft, typename  TLeft, Integer NQPoints,Integer Dim>
// class Division< QPExpression <DerivedLeft, TLeft, NQPoints,Dim>, Real >

// : public QPExpression<Division<QPExpression <DerivedLeft, TLeft, NQPoints,Dim>,Real>,
//                       typename TLeft::type,NQPoints, Dim >
// {
//   public:
//     using Left= QPExpression<DerivedLeft, TLeft, NQPoints,Dim>;
//     using Right=Real;
//     Division(const Left& left, const Right&right): 
//     left_(left.derived()),
//     right_(right)
//     {};

//     template<typename QP>
//     typename Left::type  eval(const QP & qp_points)
//     {return left_.eval(qp_points)/right_;};

//   private:
//   DerivedLeft left_;
//   Real right_;
// };


// template< typename DerivedLeft, typename  TLeft, Integer NQPoints,Integer Dim>
// class Division< QPExpression<DerivedLeft, TLeft, NQPoints,Dim>,Real >
// operator/(const QPExpression<DerivedLeft, TLeft, NQPoints,Dim>& left, const Real&right)
// {return Division<QPExpression<DerivedLeft, TLeft, NQPoints,Dim>, Real >(left,right);}


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////         MULTIPLY: QP * REAL       /////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// template< typename DerivedLeft, typename  TLeft, Integer NQPoints,Integer Dim>
// class Multiplication< QPExpression <DerivedLeft, TLeft, NQPoints,Dim>, Real >

// : public QPExpression<Multiplication<QPExpression <DerivedLeft, TLeft, NQPoints,Dim>,Real>,
//                       typename TLeft::type,NQPoints, Dim >
// {
//   public:
//     using Left= QPExpression<DerivedLeft, TLeft, NQPoints,Dim>;
//     using Right=Real;
//     Multiplication(const Left& left, const Right&right): 
//     left_(left.derived()),
//     right_(right)
//     {};

//     template<typename QP>
//     typename Left::type  eval(const QP & qp_points)
//     {return left_.eval(qp_points)*right_;};

//   private:
//   DerivedLeft left_;
//   Real right_;
// };


// template< typename DerivedLeft, typename  TLeft, Integer NQPoints,Integer Dim>
// class Multiplication< QPExpression<DerivedLeft, TLeft, NQPoints,Dim>,Real >
// operator*(const QPExpression<DerivedLeft, TLeft, NQPoints,Dim>& left, const Real&right)
// {return Multiplication<QPExpression<DerivedLeft, TLeft, NQPoints,Dim>, Real >(left,right);}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////         MULTIPLY: REAL * QP       /////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// template< typename DerivedRight, typename  Tright, Integer NQPoints,Integer Dim>
// class Multiplication<Real, QPExpression <DerivedRight, Tright, NQPoints,Dim> >

// : public QPExpression<Multiplication<Real,QPExpression <DerivedRight, Tright, NQPoints,Dim>>,
//                       typename Tright::type,NQPoints, Dim >
// {
//   public:
//     using Right= QPExpression<DerivedRight, Tright, NQPoints,Dim>;
//     using Left=Real;
//     Multiplication(const Left& left, const Right&right): 
//     right_(right.derived()),
//     left_(left)
//     {};

//     template<typename QP>
//     typename Right::type  eval(const QP & qp_points)
//     {return left_*right_.eval(qp_points);};

//   private:
//   Real left_;
//   DerivedRight right_;
// };


// template< typename DerivedRight, typename  Tright, Integer NQPoints,Integer Dim>
// class Multiplication< Real,QPExpression<DerivedRight, Tright, NQPoints,Dim> >
// operator*(const Real&left,const QPExpression<DerivedRight, Tright, NQPoints,Dim>& right)
// {return Multiplication<Real,QPExpression<DerivedRight, Tright, NQPoints,Dim> >(left,right);}


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////          MULTIPLY: QP * QP        /////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// template< typename DerivedLeft,typename DerivedRight, typename  TLeft,typename TRight, Integer NQPoints,Integer Dim>
// class Multiplication< QPExpression <DerivedLeft, TLeft, NQPoints,Dim>, 
//                      QPExpression <DerivedRight,TRight,NQPoints, Dim> >

// : public QPExpression<Multiplication<QPExpression <DerivedLeft, TLeft, NQPoints,Dim>,
//                                QPExpression <DerivedRight,TRight,NQPoints,Dim> >,
//                 typename OperatorType<Multiplication<typename TLeft::type,typename TRight::type>>::type,
//                 NQPoints, Dim >
// {
//   public:
//     using Left= QPExpression<DerivedLeft, TLeft, NQPoints,Dim>;
//     using Right=QPExpression<DerivedRight,TRight,NQPoints,Dim>;
//     Multiplication(const Left& left, const Right&right): 
//     left_(left.derived()),
//     right_(right.derived())
//     {};

//     template<typename QP>
//     typename OperatorType< Multiplication<typename Left::type, typename Right::type>>::type  eval(const QP & qp_points)
//     {return left_.eval(qp_points)*right_.eval(qp_points);};

//   private:
//   DerivedLeft left_;
//   DerivedRight right_;
// };


// template< typename DerivedLeft,typename DerivedRight, typename  TLeft,typename TRight, Integer NQPoints,Integer Dim>
// class Multiplication< QPExpression<DerivedLeft, TLeft, NQPoints,Dim>, 
//                      QPExpression<DerivedRight,TRight,NQPoints,Dim> >
// operator*(const QPExpression<DerivedLeft, TLeft, NQPoints,Dim>& left, 
//           const QPExpression<DerivedRight,TRight,NQPoints, Dim>&right)
// {return Multiplication<QPExpression<DerivedLeft, TLeft, NQPoints,Dim>, 
//                       QPExpression<DerivedRight,TRight,NQPoints,Dim> >(left,right);}







//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////          UNARY PLUS: +FQP        //////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// template< typename Derived, typename  T, Integer NQPoints,Integer NComponents,Integer Dim>
// class UnaryPlus< FQPExpression <Derived, T, NQPoints,NComponents,Dim>>
// : public FQPExpression<UnaryPlus< FQPExpression <Derived, T, NQPoints,NComponents,Dim>>,
//                       typename OperatorType<UnaryPlus<typename T::type>>::type,
//                       NQPoints,NComponents ,  Dim >
// {
//   public:
//     using Input= FQPExpression<Derived,T,NQPoints,NComponents,Dim>;
//     UnaryPlus(const Input& input): 
//     derived_(input.derived())
//     {};

//     template<typename QP>
//     typename OperatorType< UnaryPlus<typename Input::type>>::type  eval(const QP & qp_points)
//     {return +derived_.eval(qp_points);};

//   private:
//   Derived derived_;
// };


// template< typename Derived, typename  T, Integer NQPoints,Integer NComponents, Integer Dim>
// class UnaryPlus< FQPExpression <Derived, T, NQPoints,NComponents,Dim> >
// operator+(const FQPExpression<Derived, T, NQPoints,NComponents,Dim>& input)
//           {return UnaryPlus<FQPExpression <Derived, T, NQPoints,NComponents,Dim>>(input);}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////          UNARY MINUS: +FQP       //////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// template< typename Derived, typename  T, Integer NQPoints,Integer NComponents,Integer Dim>
// class UnaryMinus< FQPExpression <Derived, T, NQPoints,NComponents,Dim>>
// : public FQPExpression<UnaryMinus< FQPExpression <Derived, T, NQPoints,NComponents,Dim>>,
//                       typename OperatorType<UnaryMinus<typename T::type>>::type,
//                       NQPoints,NComponents ,  Dim >
// {
//   public:
//     using Input= FQPExpression<Derived,T,NQPoints,NComponents,Dim>;
//     UnaryMinus(const Input& input): 
//     derived_(input.derived())
//     {};

//     template<typename QP>
//     typename OperatorType< UnaryMinus<typename Input::type>>::type  eval(const QP & qp_points)
//     {return -derived_.eval(qp_points);};

//   private:
//   Derived derived_;
// };


// template< typename Derived, typename  T, Integer NQPoints,Integer NComponents, Integer Dim>
// class UnaryMinus< FQPExpression <Derived, T, NQPoints,NComponents,Dim> >
// operator-(const FQPExpression<Derived, T, NQPoints,NComponents,Dim>& input)
//           {return UnaryMinus<FQPExpression <Derived, T, NQPoints,NComponents,Dim>>(input);}


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////         ADDITION: FQP + FQP       /////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// template< typename DerivedLeft,typename DerivedRight, typename  TLeft,typename TRight, Integer NQPoints,Integer NComponents,Integer Dim>
// class Addition< FQPExpression <DerivedLeft,  TLeft,  NQPoints, NComponents, Dim>, 
//                 FQPExpression <DerivedRight, TRight, NQPoints, NComponents, Dim> >

// : public FQPExpression<Addition<FQPExpression <DerivedLeft,  TLeft,  NQPoints, NComponents, Dim>,
//                                FQPExpression <DerivedRight, TRight, NQPoints, NComponents, Dim> >,
//                        typename OperatorType<Addition<typename TLeft::type,typename TRight::type>>::type,
//                        NQPoints,NComponents, Dim >
// {
//   public:
//     using Left= FQPExpression <DerivedLeft,  TLeft,  NQPoints, NComponents, Dim>;
//     using Right=FQPExpression <DerivedRight, TRight, NQPoints, NComponents, Dim>;
//     Addition(const Left& left, const Right&right): 
//     left_(left.derived()),
//     right_(right.derived())
//     {};

//     template<typename QP>
//     typename OperatorType< Addition<typename Left::type, typename Right::type>>::type  eval(const QP & qp_points)
//     {return left_.eval(qp_points)+right_.eval(qp_points);};

//   private:
//   DerivedLeft left_;
//   DerivedRight right_;
// };


// template< typename DerivedLeft,typename DerivedRight, typename  TLeft,typename TRight, Integer NQPoints,Integer NComponents,Integer Dim>
// class Addition< FQPExpression <DerivedLeft,  TLeft,  NQPoints, NComponents, Dim>, 
//                 FQPExpression <DerivedRight, TRight, NQPoints, NComponents, Dim> >
// operator+(const FQPExpression <DerivedLeft,  TLeft,  NQPoints, NComponents, Dim>& left, 
//           const FQPExpression <DerivedRight, TRight, NQPoints, NComponents, Dim>&right)
// {return Addition<FQPExpression <DerivedLeft,  TLeft,  NQPoints, NComponents, Dim>, 
//                  FQPExpression <DerivedRight, TRight, NQPoints, NComponents, Dim> >(left,right);}


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////       SUBTRACTION: FQP + FQP      /////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// template< typename DerivedLeft,typename DerivedRight, typename  TLeft,typename TRight, Integer NQPoints,Integer NComponents,Integer Dim>
// class Subtraction< FQPExpression <DerivedLeft,  TLeft,  NQPoints, NComponents, Dim>, 
//                    FQPExpression <DerivedRight, TRight, NQPoints, NComponents, Dim> >

// : public FQPExpression<Subtraction<FQPExpression <DerivedLeft,  TLeft,  NQPoints, NComponents, Dim>,
//                                    FQPExpression <DerivedRight, TRight, NQPoints, NComponents, Dim> >,
//                        typename OperatorType<Subtraction<typename TLeft::type,typename TRight::type>>::type,
//                        NQPoints,NComponents, Dim >
// {
//   public:
//     using Left= FQPExpression <DerivedLeft,  TLeft,  NQPoints, NComponents, Dim>;
//     using Right=FQPExpression <DerivedRight, TRight, NQPoints, NComponents, Dim>;
//     Subtraction(const Left& left, const Right&right): 
//     left_(left.derived()),
//     right_(right.derived())
//     {};

//     template<typename QP>
//     typename OperatorType< Subtraction<typename Left::type, typename Right::type>>::type  eval(const QP & qp_points)
//     {return left_.eval(qp_points)-right_.eval(qp_points);};

//   private:
//   DerivedLeft left_;
//   DerivedRight right_;
// };


// template< typename DerivedLeft,typename DerivedRight, typename  TLeft,typename TRight, Integer NQPoints,Integer NComponents,Integer Dim>
// class Subtraction< FQPExpression <DerivedLeft,  TLeft,  NQPoints, NComponents, Dim>, 
//                 FQPExpression <DerivedRight, TRight, NQPoints, NComponents, Dim> >
// operator-(const FQPExpression <DerivedLeft,  TLeft,  NQPoints, NComponents, Dim>& left, 
//           const FQPExpression <DerivedRight, TRight, NQPoints, NComponents, Dim>&right)
// {return Subtraction<FQPExpression <DerivedLeft,  TLeft,  NQPoints, NComponents, Dim>, 
//                     FQPExpression <DerivedRight, TRight, NQPoints, NComponents, Dim> >(left,right);}


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////          DIVIDE: FQP / REAL       /////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// template< typename DerivedLeft, typename  TLeft, Integer NQPoints,Integer NComponents,Integer Dim>
// class Division< FQPExpression <DerivedLeft, TLeft, NQPoints,NComponents,Dim>, Real >

// : public FQPExpression<Division<FQPExpression <DerivedLeft, TLeft, NQPoints,NComponents,Dim>,Real>,
//                        typename TLeft::type,NQPoints,NComponents, Dim >
// {
//   public:
//     using Left= FQPExpression<DerivedLeft, TLeft, NQPoints,NComponents,Dim>;
//     using Right=Real;
//     Division(const Left& left, const Right&right): 
//     left_(left.derived()),
//     right_(right)
//     {};

//     template<typename QP>
//     typename Left::type  eval(const QP & qp_points)
//     {return left_.eval(qp_points)/right_;};

//   private:
//   DerivedLeft left_;
//   Real right_;
// };


// template< typename DerivedLeft, typename  TLeft, Integer NQPoints,Integer NComponents,Integer Dim>
// class Division< FQPExpression <DerivedLeft, TLeft, NQPoints,NComponents,Dim>,Real >
// operator/(const FQPExpression <DerivedLeft, TLeft, NQPoints,NComponents,Dim>& left, const Real&right)
// {return Division<FQPExpression <DerivedLeft, TLeft, NQPoints,NComponents,Dim>, Real >(left,right);}


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////         MULTIPLY: FQP * REAL       /////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// template< typename DerivedLeft, typename  TLeft, Integer NQPoints,Integer NComponents,Integer Dim>
// class Multiplication< FQPExpression <DerivedLeft, TLeft, NQPoints,NComponents,Dim>, Real >

// : public FQPExpression<Multiplication<FQPExpression <DerivedLeft, TLeft, NQPoints,NComponents,Dim>,Real>,
//                        typename TLeft::type,NQPoints,NComponents, Dim >
// {
//   public:
//     using Left= FQPExpression<DerivedLeft, TLeft, NQPoints,NComponents,Dim>;
//     using Right=Real;
//     Multiplication(const Left& left, const Right&right): 
//     left_(left.derived()),
//     right_(right)
//     {};

//     template<typename QP>
//     typename Left::type  eval(const QP & qp_points)
//     {return left_.eval(qp_points)*right_;};

//   private:
//   DerivedLeft left_;
//   Real right_;
// };


// template< typename DerivedLeft, typename  TLeft, Integer NQPoints,Integer NComponents,Integer Dim>
// class Multiplication< FQPExpression <DerivedLeft, TLeft, NQPoints,NComponents,Dim>,Real >
// operator*(const FQPExpression <DerivedLeft, TLeft, NQPoints,NComponents,Dim>& left, const Real&right)
// {return Multiplication<FQPExpression <DerivedLeft, TLeft, NQPoints,NComponents,Dim>, Real >(left,right);}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////         MULTIPLY: REAL * FQP       /////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// template< typename DerivedRight, typename  Tright, Integer NQPoints,Integer NComponents,Integer Dim>
// class Multiplication<Real, FQPExpression <DerivedRight, Tright, NQPoints,NComponents,Dim> >

// : public FQPExpression<Multiplication<Real,FQPExpression <DerivedRight, Tright, NQPoints,NComponents,Dim>>,
//                        typename Tright::type,NQPoints,NComponents, Dim >
// {
//   public:
//     using Right=FQPExpression <DerivedRight, Tright, NQPoints,NComponents,Dim>;
//     using Left=Real;
//     Multiplication(const Left& left, const Right&right): 
//     right_(right.derived()),
//     left_(left)
//     {};

//     template<typename QP>
//     typename Right::type  eval(const QP & qp_points)
//     {return left_*right_.eval(qp_points);};

//   private:
//   Real left_;
//   DerivedRight right_;
// };


// template< typename DerivedRight, typename  Tright, Integer NQPoints,Integer NComponents,Integer Dim>
// class Multiplication< Real,FQPExpression <DerivedRight, Tright, NQPoints,NComponents,Dim> >
// operator*(const Real&left,const FQPExpression <DerivedRight, Tright, NQPoints,NComponents,Dim>& right)
// {return Multiplication<Real,FQPExpression <DerivedRight, Tright, NQPoints,NComponents,Dim> >(left,right);}



//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////          MULTIPLY: QP * FQP      //////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// template< typename DerivedLeft,typename DerivedRight, typename  TLeft,typename TRight, Integer NQPoints,Integer NComponents,Integer Dim>
// class Multiplication< QPExpression <DerivedLeft, TLeft, NQPoints,Dim>, 
//                      FQPExpression<DerivedRight,TRight,NQPoints,NComponents, Dim> >

// : public FQPExpression<Multiplication<QPExpression <DerivedLeft, TLeft, NQPoints,Dim>,
//                               FQPExpression<DerivedRight,TRight,NQPoints,NComponents, Dim> >,
//                 typename OperatorType<Multiplication<typename TLeft::type,typename TRight::type>>::type,
//                 NQPoints, NComponents, Dim >
// {
//   public:
//     using Left=QPExpression <DerivedLeft, TLeft, NQPoints,Dim>;
//     using Right=FQPExpression<DerivedRight,TRight,NQPoints,NComponents, Dim>;
//     Multiplication(const Left& left, const Right&right): 
//     left_(left.derived()),
//     right_(right.derived())
//     {};

//     template<typename QP>
//     typename OperatorType< Multiplication<typename Left::type, typename Right::type>>::type  eval(const QP & qp_points)
//     {return left_.eval(qp_points)*right_.eval(qp_points);};

//   private:
//   DerivedLeft left_;
//   DerivedRight right_;
// };

// template< typename DerivedLeft,typename DerivedRight, typename  TLeft,typename TRight, Integer NQPoints,Integer NComponents,Integer Dim>
// class Multiplication< QPExpression <DerivedLeft, TLeft, NQPoints,Dim>, 
//                      FQPExpression<DerivedRight,TRight,NQPoints,NComponents, Dim> >
// operator*(const QPExpression <DerivedLeft, TLeft, NQPoints,Dim>& left, 
//           const FQPExpression<DerivedRight,TRight,NQPoints,NComponents, Dim>&right)
// {
// return Multiplication<QPExpression <DerivedLeft, TLeft, NQPoints,Dim>, FQPExpression<DerivedRight,TRight,NQPoints,NComponents, Dim> >(left,right);//static_cast<const Left>(u), static_cast<const Right>(v));
// }


}

#endif