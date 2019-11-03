#ifndef MARS_CONSTANT_HPP
#define MARS_CONSTANT_HPP

#include "mars_base.hpp"


namespace mars{

class Identity2 
{
public:
    inline static constexpr auto eval()
    {return Matrix<Real,2,2>{1.0,0.0, 0.0,1.0};}
};

class Identity3
{
public:
    inline static constexpr auto eval()
    {return Matrix<Real,3,3>{1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0};}
};
class Identity4
{
public:
    inline static constexpr auto eval()
    {return Matrix<Real,4,4>{1.0,0.0,0.0,0.0, 0.0,1.0,0.0,0.0, 0.0,0.0,1.0,0.0, 0.0,0.0,0.0,1.0};}
};

class Eye2 
{
public:
    inline static constexpr Matrix<Real,2,2> eval(const Real alpha)
    {return Matrix<Real,2,2>{alpha,0.0, 0.0,alpha};}
};

class Eye3
{
public:
    inline static constexpr auto eval(const Real alpha)
    {return Matrix<Real,3,3>{alpha,0.0,0.0, 0.0,alpha,0.0, 0.0,0.0,alpha};}
};
class Eye4
{
public:
    inline static constexpr auto eval(const Real alpha)
    {return Matrix<Real,4,4>{alpha,0.0,0.0,0.0, 0.0,alpha,0.0,0.0, 0.0,0.0,alpha,0.0, 0.0,0.0,0.0,alpha};}
};

class Prova 
{
public:
    inline static constexpr Matrix<Real,2,2> eval(const Real alpha,const Real beta)
    {return Matrix<Real,2,2>{alpha,beta,beta,alpha+beta};}
};

class Scalar1 
{
public:
    inline static constexpr Matrix<Real,1,1> eval(const Real alpha)
    {return Matrix<Real,1,1>{alpha};}
};


class Function1
{
    public: 
    using Point=Vector<Real,2>;
    using type=Real;

    static type eval(const Point& p)
    {
     return p[0]; 
    }
};

class Vec2
{
    public: 
    using type=Vector<Real,2>;

    
    static constexpr type eval(const Real& x1, const Real& x2)
    {
     return type(x1,x2); 
    }
};



class Mat1
{
    public: 
    using type=Matrix<Real,1,2>;

    
    static constexpr type eval(const Real& x1, const Real& x2)
    {
     return type(x1,x2); 
    }
};

class Mat2
{
    public: 
    using type=Matrix<Real,2,1>;

    
    static constexpr type eval(const Real& x1, const Real& x2)
    {
     return type(x1,x2); 
    }
};



template<typename ConstType,typename...Inputs>
class ConstantTensor;

template<typename ConstType,typename...Inputs>
class ConstantTensor: public Expression<ConstantTensor<ConstType,Inputs...>>
{
 public:
  
  using type=decltype(ConstType::eval(Inputs()...));

  template<Integer N>
  using qptype=QPValues<type,N>;
  
  constexpr ConstantTensor(const Inputs&...inputs):
  tensor_(ConstType::eval(inputs...))//,
  // tuple_(std::make_tuple(inputs...))
  {}
  
  template<typename...OtherInputs>
  constexpr ConstantTensor(const OtherInputs&...other_inputs){}

  constexpr ConstantTensor(const ConstantTensor& other_input):
  tensor_(other_input.eval())//,
  // tuple_(std::make_tuple(inputs...))
  {}
  
  
  constexpr type eval()const {return tensor_;};

  template<Integer N>
  constexpr qptype<N> qp_eval()const
  {qptype<N> qpvalue;
    for(Integer n=0;n<N;n++)
      Assignment<type>::apply(qpvalue[n],tensor_);
    return qpvalue;};

private:
  type tensor_;
  // std::tuple<Inputs...> tuple_;
};


template<typename ConstType,typename...Inputs>
constexpr auto Constant(const Inputs&...inputs){return ConstantTensor<ConstType,Inputs...>(inputs...);}




}
#endif
