#ifndef MARS_FUNCTION_HPP
#define MARS_FUNCTION_HPP

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
    using Output=Real;

    static Output eval(const Point& p)
    {
     return 1; 
    }
};







template<typename ConstType,typename...Inputs>
class ConstantTensor;

template<typename ConstType,typename...Inputs>
class ConstantTensor: public Expression<ConstantTensor<ConstType,Inputs...>>
{
 public:
  
  using Output=decltype(ConstType::eval(Inputs()...));
  
  constexpr ConstantTensor(const Inputs&...inputs):
  tensor_(ConstType::eval(inputs...)),
  tuple_(std::make_tuple(inputs...))
  {}
  
  template<typename...OtherInputs>
  constexpr ConstantTensor(const OtherInputs&...other_inputs){}
  
  constexpr Output eval()const {return tensor_;};

private:
  Output tensor_;
  std::tuple<Inputs...> tuple_;
};


template<typename ConstType,typename...Inputs>
class QuadratureOrder<ConstantTensor<ConstType,Inputs...>>
{ public:
  static constexpr Integer value=0;
};

template<typename ConstType,typename...Inputs>
constexpr auto Constant(const Inputs&...inputs){return ConstantTensor<ConstType,Inputs...>(inputs...);}



// template<typename Right,typename ConstType,typename...Inputs>
// class OperatorTupleType< Multiplication< Expression<ConstantTensor<ConstType,Inputs...>>, Expression<Right> > >
// { public:
//   using type=typename OperatorTupleType<Right>::type;
// };

// template<typename Left,typename ConstType,typename...Inputs>
// class OperatorTupleType< Multiplication<Expression<Left>, Expression<ConstantTensor<ConstType,Inputs...>> > >
// { public:
//   using type=typename OperatorTupleType<Left>::type;
// };






template<typename FullSpace,Integer N,typename Operator_=IdentityOperator,typename FuncType=EmptyClass>
class Function2;





template<typename FuncType,typename Space>
class QPPointsFunction;


template<typename FuncType,typename Elem,Integer NContinuity,Integer NComponents>
class QPPointsFunction<FuncType,ElementFunctionSpace<Elem,LagrangeFE,1,NContinuity,NComponents>>
{
public:
 static constexpr Integer Ndofs=DofsPerElemNums<Elem,BaseFunctionSpace<LagrangeFE,1,NContinuity,NComponents>>::value;
 using Output=typename FuncType::Output;
    inline static void init(const Jacobian<Elem>& J,Array<Output,Ndofs>& local_dofs)
    {
      const auto& points=J.points();
      for(Integer ii=0;ii<points.size();ii++)
          std::cout<<points[ii]<<std::endl;
      std::cout<<local_dofs<<std::endl;
      for(Integer ii=0;ii<points.size();ii++)
        local_dofs[ii]=FuncType::eval(points[ii]);
    }

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Function2<W,N,Operator,FuncType> is an expression which belongs to the N-th space of W
////// If FuncType=EmptySpace, then it is evaluated only based on a coefficients-vector of a trial function
////// If FuncType=FuncType, then the local dofs are computed in terms of a local lagrangian interpolation
///////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename FullSpace,Integer N,typename Operator_, typename FuncType>
class Function2
: public Expression<Function2<FullSpace,N,Operator_,FuncType>>
{
  
  public:
    using FunctionSpace=FullSpace;
    using FunctionType=FuncType;
    using Operator=Operator_;
    static constexpr Integer number=N;
    static constexpr Integer value=GetType<typename FunctionSpace::FromElementFunctionSpacesToUniqueNumbersTupleType,N>::value;
    using Elem=typename FunctionSpace::Elem;
    static constexpr Integer Dim=Elem::Dim;
    using UniqueElementFunctionSpacesTupleType=typename FunctionSpace::UniqueElementFunctionSpacesTupleType;
    static constexpr Integer Nmax=TupleTypeSize<UniqueElementFunctionSpacesTupleType>::value;
    static constexpr Integer Nelem_dofs=FunctionSpace::Nelem_dofs_array[N];
    using LocalDofs=Array<Real,Nelem_dofs>;
    using Space=GetType<UniqueElementFunctionSpacesTupleType,value>;
    static constexpr Integer FEFamily=Space::FEFamily;
    static constexpr Integer Order=Space::Order;
    static constexpr Integer NComponents=Space::NComponents;

  
  Function2(const std::shared_ptr<FullSpace>& AuxW_ptr):
  full_space_ptr_(AuxW_ptr)
  {}
  
  Function2(const FullSpace& AuxW):
  full_space_ptr_(std::make_shared<FullSpace>(AuxW))
  {} 

   void local_dofs_update(const Jacobian<Elem>& J)
    {std::cout<<"function function2 dofs update"<<std::endl;
      QPPointsFunction<FunctionType,Space>::init(J,local_dofs_);} 

   inline const LocalDofs& local_dofs()const{return local_dofs_;}


  

  // void update(std::vector<Real>& global_dofs)
  // {if(global_dofs_.size()==global_dofs.size())
  //     copy(global_dofs.begin(),global_dofs.end(),back_inserter(global_dofs));
  //  else
  //  {
  //   global_dofs.resize(global_dofs_.size());
  //   copy(global_dofs.begin(),global_dofs.end(),back_inserter(global_dofs));
  //  }
  //  }

private:
    std::shared_ptr<FullSpace> full_space_ptr_;
    LocalDofs local_dofs_;
    std::vector<Real> global_dofs_;

};


template<typename FullSpace,Integer N,typename Operator_>
class Function2<FullSpace,N,Operator_,EmptyClass>
: public Expression<Function2<FullSpace,N,Operator_,EmptyClass>>
{
  
  public:
    
    using FunctionSpace=FullSpace;
    static constexpr Integer TrialSpaceSize=FunctionSpace::TrialSpaceSize;
    static constexpr Integer FunctionNumber=N-(TrialSpaceSize);
    using Operator=Operator_;
    static constexpr Integer number=N;
    static constexpr Integer value=GetType<typename FunctionSpace::FromElementFunctionSpacesToUniqueNumbersTupleType,N>::value;
    using Elem=typename FunctionSpace::Elem;
    static constexpr Integer Dim=Elem::Dim;
    using UniqueElementFunctionSpacesTupleType=typename FunctionSpace::UniqueElementFunctionSpacesTupleType;
    // using ElementFunctionSpace=typename FunctionSpace::ElementFunctionSpace;
    static constexpr Integer Nmax=TupleTypeSize<UniqueElementFunctionSpacesTupleType>::value;
    static constexpr Integer Nelem_dofs=FunctionSpace::Nelem_dofs_array[N];
    using LocalDofs=Array<Real,Nelem_dofs>;

  // template<typename...Inputs>
  // static auto eval(const Inputs&...inputs){FuncType::eval(inputs...);}

  // Function2(){}

  Function2(const std::shared_ptr<FullSpace>& AuxW_ptr):
  full_space_ptr_(AuxW_ptr)
  {}
  
  Function2(const FullSpace& AuxW):
  full_space_ptr_(std::make_shared<FullSpace>(AuxW))
  {} 

  inline auto full_space_ptr(){return full_space_ptr_;}
  inline const auto& full_space_ptr()const{return full_space_ptr_;}
  
   void local_dofs_update(const Jacobian<Elem>& J)
    {
    auto tuple_map=full_space_ptr_->aux_space_ptr()->dofmap();
    auto local_map=tuple_get<FunctionNumber>(full_space_ptr_->aux_space_ptr()->dofmap())[J.elem_id()];
     for(int ii=0;ii<local_map.size();ii++)
      local_dofs_[ii]=global_dofs_[local_map[ii]];
    } 
   inline const LocalDofs& local_dofs()const{return local_dofs_;}


  // const LocalDofs& local_dofs()const
  // {
  // const auto aux_space_ptr=full_space_ptr_->aux_space_ptr();
  // std::cout<<"eeee"<<std::endl;
  // aux_space_ptr->dofmap();
  // std::cout<<"eeee2"<<std::endl;
  // return local_dofs_;
  // }


  void global_dofs_update(std::vector<Real>& global_dofs)
  {if(global_dofs_.size()==global_dofs.size())
      copy(global_dofs.begin(),global_dofs.end(),back_inserter(global_dofs));
   else
   {
    global_dofs.resize(global_dofs_.size());
    copy(global_dofs.begin(),global_dofs.end(),back_inserter(global_dofs));
   }
   }
private:
    std::shared_ptr<FullSpace> full_space_ptr_;
    LocalDofs local_dofs_;
    std::vector<Real> global_dofs_;

};


template<Integer N,typename OperatorType,typename FullSpace,typename FuncType>
class QuadratureOrder<Function2<FullSpace,N,OperatorType,FuncType>>
{ public:
  using type=Function2<FullSpace,N,OperatorType,FuncType>;
  // using UniqueElementFunctionSpacesTupleType=typename Test::UniqueElementFunctionSpacesTupleType;
  // using BaseFunctionSpaceAndElement=GetType<UniqueElementFunctionSpacesTupleType,N>;
  using Operator=typename type::Operator;
  using UniqueElementFunctionSpacesTupleType=typename type::UniqueElementFunctionSpacesTupleType;
  using BaseFunctionSpaceAndElement=GetType<UniqueElementFunctionSpacesTupleType,type::value>;
  using BaseFunctionSpace=Elem2FunctionSpace<BaseFunctionSpaceAndElement>;
  static constexpr Integer value=QuadratureOrder<Operator,BaseFunctionSpace>::value;

};



template<Integer N,typename FuncType=EmptyClass,typename OperatorType=IdentityOperator,typename FullSpace>
auto MakeFunction(const FullSpace& AuxW){return Function2<FullSpace,N+FullSpace::TrialSpaceSize,OperatorType,FuncType>(AuxW);}


}




#endif
