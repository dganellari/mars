#ifndef MARS_FUNCTION_HPP
#define MARS_FUNCTION_HPP

#include "mars_base.hpp"
#include "mars_constant.hpp"

namespace mars{

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


class Function2
{
    public: 
    using Point=Matrix<Real,2,1>;
    using type=Real;
    
    
    static type eval(const Point& p)
    {
     return p(0,0)+p(1,0); 
    }
};

template<typename FullSpace,Integer N,typename Operator_=IdentityOperator,typename FuncType=EmptyClass>
class Function;

template<typename Elem>
class FiniteElem;



template<typename FuncType,typename Space>
class QPPointsFunction;



template<typename FuncType,typename Elem,Integer NContinuity,Integer NComponents>
class QPPointsFunction<FuncType,ElementFunctionSpace<Elem,LagrangeFE,1,NContinuity,NComponents>>
{
public:
 static constexpr Integer Ndofs=DofsPerElemNums<Elem,BaseFunctionSpace<LagrangeFE,1,NContinuity,NComponents>>::value;
 using type=typename FuncType::type;
    inline static void init(const FiniteElem<Elem>& J,Array<type,Ndofs>& local_dofs)
    {
      const auto& points=J.points();
      // for(Integer ii=0;ii<points.size();ii++)
      //     std::cout<<points[ii]<<std::endl;
      for(Integer ii=0;ii<points.size();ii++)
        local_dofs[ii]=FuncType::eval(points[ii]);
    }

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Function<W,N,Operator,FuncType> is an expression which belongs to the N-th space of W
////// If FuncType=EmptySpace, then it is evaluated only based on a coefficients-vector of a trial function
////// If FuncType=FuncType, then the local dofs are computed in terms of a local lagrangian interpolation
///////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename FullSpace,Integer N,typename Operator_, typename FuncType>
class Function
: public Expression<Function<FullSpace,N,Operator_,FuncType>>
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

  
  Function(const std::shared_ptr<FullSpace>& AuxW_ptr):
  full_spaces_ptr_(AuxW_ptr)
  {}
  
  Function(const FullSpace& AuxW):
  full_spaces_ptr_(std::make_shared<FullSpace>(AuxW))
  {} 

   void local_dofs_update(const FiniteElem<Elem>& J)
    {QPPointsFunction<FunctionType,Space>::init(J,local_dofs_);} 

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

   constexpr auto spaces_ptr(){return full_spaces_ptr_;}
   constexpr auto spaces_ptr()const {return full_spaces_ptr_;}


private:
    std::shared_ptr<FullSpace> full_spaces_ptr_;
    LocalDofs local_dofs_;
    std::vector<Real> global_dofs_;

};


template<typename FullSpace,Integer N,typename Operator_>
class Function<FullSpace,N,Operator_,EmptyClass>
: public Expression<Function<FullSpace,N,Operator_,EmptyClass>>
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

  // Function(){}

  Function(const std::shared_ptr<FullSpace>& AuxW_ptr):
  full_spaces_ptr_(AuxW_ptr)
  {}
  
  Function(const FullSpace& AuxW):
  full_spaces_ptr_(std::make_shared<FullSpace>(AuxW))
  {} 

  inline auto full_spaces_ptr(){return full_spaces_ptr_;}
  inline const auto& full_spaces_ptr()const{return full_spaces_ptr_;}
  
   void local_dofs_update(const FiniteElem<Elem>& J)
    {
    auto tuple_map=full_spaces_ptr_->aux_spaces_ptr()->dofmap();
    auto local_map=tuple_get<FunctionNumber>(full_spaces_ptr_->aux_spaces_ptr()->dofmap())[J.elem_id()];
     for(int ii=0;ii<local_map.size();ii++)
      local_dofs_[ii]=global_dofs_[local_map[ii]];
    } 
   inline const LocalDofs& local_dofs()const{return local_dofs_;}


  // const LocalDofs& local_dofs()const
  // {
  // const auto aux_spaces_ptr=full_spaces_ptr_->aux_spaces_ptr();
  // std::cout<<"eeee"<<std::endl;
  // aux_spaces_ptr->dofmap();
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

   constexpr auto spaces_ptr(){return full_spaces_ptr_;}
   constexpr auto spaces_ptr()const {return full_spaces_ptr_;}

private:
    std::shared_ptr<FullSpace> full_spaces_ptr_;
    LocalDofs local_dofs_;
    std::vector<Real> global_dofs_;

};




template<Integer N,typename FuncType=EmptyClass,typename OperatorType=IdentityOperator,typename FullSpace>
auto MakeFunction(const FullSpace& AuxW){return Function<FullSpace,N+FullSpace::TrialSpaceSize,OperatorType,FuncType>(AuxW);}

template<Integer N,typename FuncType,typename OperatorType=IdentityOperator,typename FullSpace>
auto MakeFunction(const FullSpace& AuxW, const FuncType& func_type){return Function<FullSpace,N+FullSpace::TrialSpaceSize,OperatorType,FuncType>(AuxW);}


template<typename FullSpace,Integer N,typename OperatorType,typename FuncType>
constexpr auto 
Trace(const Function<FullSpace,N,OperatorType,FuncType>& t)
{return Function<FullSpace,N,TraceOperator,FuncType> (t.spaces_ptr());}




}




#endif
