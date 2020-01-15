#ifndef MARS_FUNCTION_HPP
#define MARS_FUNCTION_HPP

#include "mars_base.hpp"
#include "mars_constant.hpp"
#include "mars_shape_function.hpp"


namespace mars{

template<Integer Dim>
class NormalFunction
{
    public: 
    // using Point=Vector<Real,2>;
    using type=Matrix<Real,Dim,1>;
    
    // template<typename Point>
    // static type eval(const Point& p)
    // {
    //  return (0.5+p[0])*(0.5+p[0])*(0.5+p[0])*(0.5+p[0])*10.0; 
    // }
};



class Function1
{
    public: 
    // using Point=Vector<Real,2>;
    using type=Matrix<Real,1,1>;
    
    template<typename Point>
    static type eval(const Point& p)
    {
     return (0.5+p[0])*(0.5+p[0])*(0.5+p[0])*(0.5+p[0])*10.0; 
    }
};


class Function2
{
    public: 
    // using Point=Matrix<Real,2,1>;
    using type=Matrix<Real,1,1>;
    template<typename Point>
    static type eval(const Point& p)
    {
     return 0.0;//p(0,0)+p(1,0); 
    }
};


class Function3
{
    public: 
    // using Point=Matrix<Real,2,1>;
    using type=Matrix<Real,1,1>;
    template<typename Point>
    static type eval(const Point& p)
    {
     return -p(0,0)*p(0,0)+p(0,0); 
    }
};


class Function4
{
    public: 
    // using Point=Matrix<Real,2,1>;
    using type=Matrix<Real,1,1>;
    template<typename Point>
    static type eval(const Point& p)
    {
     return 1.0; 
    }
};

template<Integer Dim>
class FunctionOne
{
    public: 
    // using Point=Matrix<Real,Dim,1>;
    using type=Matrix<Real,1,1>;
    template<typename Point>
    static type eval(const Point& p)
    {
     return 1.0; 
    }
};

class FunctionZero1D
{
    public: 
    // using Point=Matrix<Real,2,1>;
    using type=Matrix<Real,1,1>;
    template<typename Point>
    static type eval(const Point& p)
    {
     Matrix<Real,1,1> func{0.0};
     return func; 
    }
};

class FunctionZero2D
{
    public: 
    // using Point=Matrix<Real,2,1>;
    using type=Matrix<Real,2,1>;
    template<typename Point>
    static type eval(const Point& p)
    {
     Matrix<Real,2,1> func{0.0,0.0};
     return func; 
    }
};

class FunctionZero3D
{
    public: 
    // using Point=Matrix<Real,3,1>;
    using type=Matrix<Real,3,1>;
    template<typename Point>
    static type eval(const Point& p)
    {
     Matrix<Real,3,1> func{0.0,0.0,0.0};
     return func; 
    }
};


template<Integer Dim>
class FunctionZero
{
    public: 
    // using Point=Matrix<Real,2,1>;
    using type=Matrix<Real,Dim,1>;
    template<typename Point>
    static type eval(const Point& p)
    {
     Matrix<Real,Dim,1> func;
     for(Integer i=0;i<Dim;i++)
      func[i]=0;
     return func; 
    }
};

template<typename FullSpace,Integer N,typename Operator_=IdentityOperator,typename FuncType=EmptyClass>
class Function;

template<typename Elem>
class FiniteElem;


template<typename Elem,Integer Order>
class ElemGeometricPoints;

template<typename Tuple>
class ElementOrder;

template<typename FuncType,typename Space>
class QPPointsFunction;

template<Integer Dim,typename Elem,Integer Order,Integer NContinuity,Integer NComponents>
class QPPointsFunction<NormalFunction<Dim>,ElementFunctionSpace<Elem,LagrangeFE,Order,NContinuity,NComponents>>
{
public:
 static constexpr Integer Ndofs=DofsPerElemNums<Elem,BaseFunctionSpace<LagrangeFE,Order,NContinuity,NComponents>>::value;
 using FuncType=NormalFunction<Dim>;
 using type=typename FuncType::type;
 using Space=ElementFunctionSpace<Elem,LagrangeFE,Order,NContinuity,NComponents>;
 using ElemPoints=ElemGeometricPoints<Elem,ElementOrder<std::tuple<Space>>::value>;
 static constexpr auto Points=ElemPoints::points;

    inline static void init(const FiniteElem<Elem>& FE,Array<Real,Ndofs>& local_dofs)
    {
    Integer cont_=0;
    Vector<Real,Dim> point;
    for(Integer ii=0;ii<Points.size();ii++)
     {
      FE.transform_point(point,Points[ii]);
      auto mesh_ptr=FE.mesh_ptr();
      auto& signed_normal=mesh_ptr->signed_normal();
      auto& local_normals=signed_normal.normals(FE.elem_id());
      // const auto& point=FE.transform_point(Points[ii]);
      // const auto point_tmp=FE.jac()*Points[ii];
      // const auto v0=FE.v0();
      // const auto point=point_tmp+v0;
      const auto& evaluation=FuncType::eval(point);

      // for(Integer jj=0;jj<evaluation.rows();jj++)
      //   {
      //   for(Integer kk=0;kk<evaluation.cols();kk++)
      //   {
      //     local_dofs[cont_]=evaluation(kk,jj);  
      //     cont_++;
      //   }

      //   }

      }

    }

};



template<typename FuncType,typename Elem,Integer Order,Integer NContinuity,Integer NComponents>
class QPPointsFunction<FuncType,ElementFunctionSpace<Elem,LagrangeFE,Order,NContinuity,NComponents>>
{
public:
 static constexpr Integer Ndofs=DofsPerElemNums<Elem,BaseFunctionSpace<LagrangeFE,Order,NContinuity,NComponents>>::value;
 using type=typename FuncType::type;
 using Space=ElementFunctionSpace<Elem,LagrangeFE,Order,NContinuity,NComponents>;
 using ElemPoints=ElemGeometricPoints<Elem,ElementOrder<std::tuple<Space>>::value>;
 static constexpr auto Points=ElemPoints::points;
 static constexpr auto Dim=Elem::Dim;

    inline static void init(const FiniteElem<Elem>& FE,Array<Real,Ndofs>& local_dofs)
    {
    Integer cont_=0;
    Vector<Real,Dim> point;
    for(Integer ii=0;ii<Points.size();ii++)
     {
      FE.transform_point(point,Points[ii]);
      // const auto& point=FE.transform_point(Points[ii]);
      // const auto point_tmp=FE.jac()*Points[ii];
      // const auto v0=FE.v0();
      // const auto point=point_tmp+v0;
      const auto& evaluation=FuncType::eval(point);

      for(Integer jj=0;jj<evaluation.rows();jj++)
        {
        for(Integer kk=0;kk<evaluation.cols();kk++)
        {
          local_dofs[cont_]=evaluation(kk,jj);  
          cont_++;
        }

        }

      }

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
  full_spaces_ptr_(AuxW_ptr),
  global_dofs_ptr_(std::make_shared<std::vector<Real>>(NULL))
  {}
  
  Function(const FullSpace& AuxW):
  full_spaces_ptr_(std::make_shared<FullSpace>(AuxW)),
  global_dofs_ptr_(std::make_shared<std::vector<Real>>(NULL))
  {} 

   void local_dofs_update(const FiniteElem<Elem>& FE)
    {





      // TODO FIXME
      QPPointsFunction<FunctionType,Space>::init(FE,local_dofs_);


    } 

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
   constexpr auto global_dofs_ptr(){return global_dofs_ptr_;}   
   constexpr auto global_dofs_ptr()const {return global_dofs_ptr_;}


private:
    std::shared_ptr<FullSpace> full_spaces_ptr_;
    LocalDofs local_dofs_;
    std::shared_ptr<std::vector<Real>> global_dofs_ptr_;
};











class Function_Aux
{

public:
  using std_vec=std::vector<Real>;

  void set(const std_vec& vec)
  {

    // std::cout<<"qui1"<<std::endl;
    // std::cout<<"global_dofssize="<<vec.size()<<std::endl;
    // for(std::size_t i=0;i<vec.size();i++)
    //   std::cout<<vec[i]<<std::endl;

    vec_ptr_=std::make_shared<std_vec>(vec);
    // std::cout<<"qui2"<<std::endl;
    // for(std::size_t i=0;i<vec_ptr_->size();i++)
    //   std::cout<<(*vec_ptr_)[i]<<std::endl;

    // auto size=vec.size();
    // if(vec_ptr_->size()!=size)
    //   vec_ptr_->resize(size);
    // for(std::size_t i=0;i<size;i++)
    //   (*vec_ptr_)[i]=vec[i];

  }


  auto vec_ptr(){return vec_ptr_;}


  auto& vec()      {return *vec_ptr_;}
  auto& vec()const {return *vec_ptr_;}

  Function_Aux(const Function_Aux& func_aux):
  vec_ptr_(std::make_shared<std_vec>(func_aux.vec()))
  {}
  Function_Aux():
  vec_ptr_(std::make_shared<std_vec>())
  {}
private:
  std::shared_ptr<std_vec> vec_ptr_;
};





template<typename FullSpace,Integer N>
class Function<FullSpace,N,IdentityOperator,EmptyClass>
: public Expression<Function<FullSpace,N,IdentityOperator,EmptyClass>>,
  public std::enable_shared_from_this<Function<FullSpace,N,IdentityOperator,EmptyClass>>
{
  
  public:
    
    using FunctionSpace=FullSpace;
    using DofsDM=typename FunctionSpace::DofsDM;
    using ElemDofMap=typename DofsDM::ElemDofMap;
    static constexpr Integer TrialSpaceSize=FunctionSpace::TrialSpaceSize;
    static constexpr Integer FunctionNumber=N-(TrialSpaceSize);
    using Operator=IdentityOperator;
    static constexpr Integer number=N;
    static constexpr Integer value=GetType<typename FunctionSpace::FromElementFunctionSpacesToUniqueNumbersTupleType,N>::value;
    using Elem=typename FunctionSpace::Elem;
    static constexpr Integer Dim=Elem::Dim;
    using UniqueElementFunctionSpacesTupleType=typename FunctionSpace::UniqueElementFunctionSpacesTupleType;
    using Space=GetType<UniqueElementFunctionSpacesTupleType,value>;
    static constexpr Integer Nmax=TupleTypeSize<UniqueElementFunctionSpacesTupleType>::value;
    static constexpr Integer Nelem_dofs=FunctionSpace::Nelem_dofs_array[N];
    using LocalDofs=Array<Real,Nelem_dofs>;

  // template<typename...Inputs>
  // static auto eval(const Inputs&...inputs){FuncType::eval(inputs...);}

  // Function(){}

  Function(const std::shared_ptr<FullSpace>& spaces_ptr):
  full_spaces_ptr_(spaces_ptr),
  identity_ptr_(NULL),
  global_dofs_ptr_(NULL),
  aux_ptr_(std::make_shared<Function_Aux>())  
  {
    // identity_ptr_=this->shared_from_this();
  }
  
  Function(const FullSpace& spaces_ptr):
  full_spaces_ptr_(std::make_shared<FullSpace>(spaces_ptr)),
  identity_ptr_(NULL) ,
  global_dofs_ptr_(NULL),
  aux_ptr_(std::make_shared<Function_Aux>())  
  {
    // identity_ptr_=this->shared_from_this();
  } 

  // Function(const std::shared_ptr<FullSpace>& spaces_ptr,const std::vector<Real>& global_dofs):
  // full_spaces_ptr_(spaces_ptr),
  // global_dofs_ptr_(std::make_shared<std::vector<Real>>(global_dofs))
  // {}
  
  // Function(const FullSpace& spaces_ptr,const std::vector<Real>& global_dofs):
  // full_spaces_ptr_(std::make_shared<FullSpace>(spaces_ptr)),
  // global_dofs_ptr_(std::make_shared<std::vector<Real>>(global_dofs))
  // {}

  // Function(const std::shared_ptr<FullSpace>& spaces_ptr,const std::shared_ptr<std::vector<Real>>& global_dofs_ptr):
  // full_spaces_ptr_(spaces_ptr),
  // global_dofs_ptr_(global_dofs_ptr)
  // {
  //   std::cout<<"son qui 2 ->"<<global_dofs_ptr->size()<<std::endl;
  //   std::cout<<"son qui 3 ->"<<global_dofs_ptr_->size()<<std::endl;
  // for(Integer i=0;i<global_dofs_ptr_->size();i++)
  //   std::cout<<(*global_dofs_ptr_)[i]<<std::endl;
  // }
  
  // Function(const FullSpace& spaces_ptr,const std::shared_ptr<std::vector<Real>>& global_dofs_ptr):
  // full_spaces_ptr_(std::make_shared<FullSpace>(spaces_ptr)),
  // global_dofs_ptr_(global_dofs_ptr)
  // {}



  Function(const Function& func):
  full_spaces_ptr_(func.full_spaces_ptr()),
  // identity_ptr_(std::make_shared<Function>(func)),
  // global_dofs_ptr_(func.global_dofs_ptr()),
  aux_ptr_(func.aux_ptr())
  {
     // identity_ptr_=func.identity_ptr();
     // std::cout<<"1 Function<FullSpace,N,Ide"<<std::endl;
     // std::cout<<aux_ptr_->vec_ptr()->size()<<std::endl;
     // std::cout<<"2 Function<FullSpace,N,IdentityOperator> ##################################################################"<<std::endl;
     // std::cout<<global_dofs_ptr_->size()<<std::endl;
     // std::cout<<"3 Function<FullSpace,N,IdentityOperator> ##################################################################"<<std::endl;

  } 

  // Function(const Function& func):
  // full_spaces_ptr_(func.full_spaces_ptr()),
  // global_dofs_ptr_(func.global_dofs_ptr())
  // {
  //   std::cout<<"##################################################################"<<std::endl;
  //    std::cout<<func.global_dofs_ptr()->size()<<std::endl;
  //    std::cout<<global_dofs_ptr_->size()<<std::endl;
  // } 


  inline auto full_spaces_ptr(){return full_spaces_ptr_;}
  inline const auto& full_spaces_ptr()const{return full_spaces_ptr_;}
  
   void local_dofs_update(const FiniteElem<Elem>& FE)
    {
      // std::cout<<"LOCAL DOFS UPDATE 1 "<<std::endl;
      auto local_map=tuple_get<N>(elemdofmap_); 
      full_spaces_ptr_->dofsdofmap().template dofmap_get<N>(local_map,FE.elem_id());
    
    const auto& global_dofs=*(aux_ptr_->vec_ptr());

    // auto& tuple_map=full_spaces_ptr_->aux_spaces_ptr()->dofmap2();
    // std::cout<<"LOCAL DOFS UPDATE 2"<<std::endl;
    // auto& local_map=tuple_get<FunctionNumber>(full_spaces_ptr_->aux_spaces_ptr()->dofmap2())[J.elem_id()];
    // std::cout<<"LOCAL DOFS UPDATE 3, local_map="<<local_map<<std::endl;
    // std::cout<<"local_dofs_update4, local_dofs_="<<local_dofs_<<std::endl;
    // std::cout<<"global_dofs.size()="<<global_dofs.size()<<std::endl;
     // for(int ii=0;ii<global_dofs.size();ii++)
     // {
     //  std::cout<<"global="<<global_dofs[ii]<<std::endl;
     // }
     // auto global_dofs_ptr_=identity_ptr_->global_dofs_ptr();


     for(int ii=0;ii<local_map.size();ii++)
     {
      // std::cout<<"local_map="<<local_map<<std::endl;
      // std::cout<<"local_dofs_="<<local_dofs_<<std::endl;
      local_dofs_[ii]=global_dofs[local_map[ii]];
     }
      
    // std::cout<<"after local_dofs_update4, local_dofs_="<<local_dofs_<<std::endl;
    } 
   inline const LocalDofs& local_dofs()const{return local_dofs_;}




  void global_dofs_update(std::vector<Real>& global_dofs)
  {
    // std::cout<<"global_dofssize="<<global_dofs.size()<<std::endl;
    aux_ptr_->set(global_dofs);
    // decltype(aux_ptr_->vec_ptr()) feee(5);
    // auto eee=std::make_shared<std::vector<Real>>(global_dofs);
    // std::cout<<"global_dofssize mmm="<<global_dofs.size()<<std::endl;
    // aux_ptr_->vec_ptr()=eee;//std::make_shared<std::vector<Real>>(global_dofs);
    // aux_ptr_->vec_ptr()=std::make_shared<std::vector<Real>>(global_dofs);
    // auto& global_dofs=*aux_ptr_->vec_ptr();
   //  // global_dofs_ptr_=std::make_shared<std::vector<Real>>(global_dofs);
   //  std::cout<<"global_dofs_update"<<std::endl;
   //  if(global_dofs_ptr_->size()==aux_ptr_->vec_ptr()->size())
   //  {
   //     std::cout<<" first resize==="<<std::endl;
   //    for(Integer ii=0;ii<global_dofs.size();ii++)
   //       {
   //        (*aux_ptr_->vec_ptr() )[ii] =global_dofs[ii];
   //        std::cout<<(*aux_ptr_->vec_ptr() )[ii]<<std::endl;
   //      }
   //    // copy(global_dofs.begin(),global_dofs.end(),back_inserter(global_dofs_));
   //  }
   // else
   // {
   //  std::cout<<" global_dofs resize==="<<global_dofs_ptr_->size()<<std::endl;
  


   //  global_dofs_ptr_->resize(global_dofs.size());
   //  for(Integer ii=0;ii<global_dofs.size();ii++)
   //       (*global_dofs_ptr_ )[ii] =global_dofs[ii];

   //  std::cout<<" global_dofs resize==="<<global_dofs_ptr_->size()<<std::endl;
   //  // copy(global_dofs.begin(),global_dofs.end(),back_inserter(global_dofs_));
   // }
   }

   constexpr auto spaces_ptr(){return full_spaces_ptr_;}
   constexpr auto spaces_ptr()const {return full_spaces_ptr_;}
   constexpr auto global_dofs_ptr(){return global_dofs_ptr_;}   
   constexpr auto global_dofs_ptr()const {return global_dofs_ptr_;}
   constexpr auto identity_ptr(){return this->shared_from_this();}   
   constexpr auto identity_ptr()const {return this->shared_from_this();}
   constexpr auto& self(){return *(this->shared_from_this());}   
   constexpr auto& self()const {return *(this->shared_from_this());}
   constexpr auto self2()const {return *(this);}
   constexpr auto self2() {return *(this);}
   constexpr auto self3()const {return (this);}
   constexpr auto self3() {return (this);}
   constexpr auto aux_ptr()const {return aux_ptr_;}
   constexpr auto aux_ptr()      {return aux_ptr_;}
private:
    std::shared_ptr<FullSpace> full_spaces_ptr_;
    LocalDofs local_dofs_;
    std::shared_ptr<Function> identity_ptr_;
    std::shared_ptr<Function_Aux> aux_ptr_;
    // std::vector<Real> global_dofs_;
    std::shared_ptr<std::vector<Real>> global_dofs_ptr_;
    ElemDofMap elemdofmap_;

};







template<typename FullSpace,Integer N,typename Operator_>
class Function<FullSpace,N,Operator_,EmptyClass>
: public Expression<Function<FullSpace,N,Operator_,EmptyClass>>
{
  
  public:
    
    using FunctionIdentity=Function<FullSpace,N,IdentityOperator,EmptyClass>;
    using FunctionSpace=FullSpace;
    using DofsDM=typename FunctionSpace::DofsDM;
    using ElemDofMap=typename DofsDM::ElemDofMap;
    static constexpr Integer TrialSpaceSize=FunctionSpace::TrialSpaceSize;
    static constexpr Integer FunctionNumber=N-(TrialSpaceSize);
    using Operator=Operator_;
    static constexpr Integer number=N;
    static constexpr Integer value=GetType<typename FunctionSpace::FromElementFunctionSpacesToUniqueNumbersTupleType,N>::value;
    using Elem=typename FunctionSpace::Elem;
    static constexpr Integer Dim=Elem::Dim;
    using UniqueElementFunctionSpacesTupleType=typename FunctionSpace::UniqueElementFunctionSpacesTupleType;
    using Space=GetType<UniqueElementFunctionSpacesTupleType,value>;
    static constexpr Integer Nmax=TupleTypeSize<UniqueElementFunctionSpacesTupleType>::value;
    static constexpr Integer Nelem_dofs=FunctionSpace::Nelem_dofs_array[N];
    using LocalDofs=Array<Real,Nelem_dofs>;

  // template<typename...Inputs>
  // static auto eval(const Inputs&...inputs){FuncType::eval(inputs...);}

  // Function(){}

  Function(const std::shared_ptr<FullSpace>& spaces_ptr):
  full_spaces_ptr_(spaces_ptr),
  identity_ptr_(NULL),
  global_dofs_ptr_(NULL)  
  {}
  
  Function(const FullSpace& spaces_ptr):
  full_spaces_ptr_(std::make_shared<FullSpace>(spaces_ptr)),
  identity_ptr_(NULL) ,
  global_dofs_ptr_(NULL)
  {} 

  Function(const std::shared_ptr<FullSpace>& spaces_ptr,const std::vector<Real>& global_dofs):
  full_spaces_ptr_(spaces_ptr),
  global_dofs_ptr_(std::make_shared<std::vector<Real>>(global_dofs))
  {}
  
  Function(const FullSpace& spaces_ptr,const std::vector<Real>& global_dofs):
  full_spaces_ptr_(std::make_shared<FullSpace>(spaces_ptr)),
  global_dofs_ptr_(std::make_shared<std::vector<Real>>(global_dofs))
  {}

  Function(const std::shared_ptr<FullSpace>& spaces_ptr,const std::shared_ptr<std::vector<Real>>& global_dofs_ptr):
  full_spaces_ptr_(spaces_ptr),
  global_dofs_ptr_(global_dofs_ptr)
  {
    // std::cout<<"son qui 2 ->"<<global_dofs_ptr->size()<<std::endl;
    // std::cout<<"son qui 3 ->"<<global_dofs_ptr_->size()<<std::endl;
  // for(Integer i=0;i<global_dofs_ptr_->size();i++)
  //   std::cout<<(*global_dofs_ptr_)[i]<<std::endl;
  }
  
  Function(const FullSpace& spaces_ptr,const std::shared_ptr<std::vector<Real>>& global_dofs_ptr):
  full_spaces_ptr_(std::make_shared<FullSpace>(spaces_ptr)),
  global_dofs_ptr_(global_dofs_ptr)
  {}

  Function(const FunctionIdentity& func_identity):
  full_spaces_ptr_(func_identity.full_spaces_ptr())
  // ,
  // identity_ptr_(std::make_shared<FunctionIdentity>(func_identity))//,
  // global_dofs_ptr_(identity_ptr_->global_dofs_ptr())
  {
  //   std::cout<<"____________##################################################################"<<std::endl;
  //    std::cout<<func_identity.global_dofs_ptr()->size()<<std::endl;
  //    std::cout<<global_dofs_ptr_->size()<<std::endl;
  } 

  Function(const Function& func):
  full_spaces_ptr_(func.full_spaces_ptr()),
  identity_ptr_(func.identity_ptr())//,
  // global_dofs_ptr_(func.global_dofs_ptr())
  {
    // std::cout<<"##################################################################"<<std::endl;
    //  std::cout<<func.global_dofs_ptr()->size()<<std::endl;
    //  std::cout<<global_dofs_ptr_->size()<<std::endl;
  } 

  // Function(const Function& func):
  // full_spaces_ptr_(func.full_spaces_ptr()),
  // global_dofs_ptr_(func.global_dofs_ptr())
  // {
  //   std::cout<<"##################################################################"<<std::endl;
  //    std::cout<<func.global_dofs_ptr()->size()<<std::endl;
  //    std::cout<<global_dofs_ptr_->size()<<std::endl;
  // } 


  inline auto full_spaces_ptr(){return full_spaces_ptr_;}
  inline const auto& full_spaces_ptr()const{return full_spaces_ptr_;}
  
   void local_dofs_update(const FiniteElem<Elem>& FE)
    {
      // std::cout<<"local_dofs_update1"<<std::endl;
      auto local_map=tuple_get<N>(elemdofmap_); 
      const auto& level=FE.level();
      full_spaces_ptr_->dofsdofmap().template dofmap_get<N>(local_map,FE.elem_id(),level);


    // auto& tuple_map=full_spaces_ptr_->aux_spaces_ptr()->dofmap2();
    // std::cout<<"local_dofs_update2"<<std::endl;
    // auto& local_map=tuple_get<FunctionNumber>(full_spaces_ptr_->aux_spaces_ptr()->dofmap2())[J.elem_id()];
    // std::cout<<"local_dofs_update3, local_map="<<local_map<<std::endl;
    // std::cout<<"local_dofs_update4, global_dofs_size="<<global_dofs_ptr_->size()<<std::endl;
    // std::cout<<"local_dofs_update4, local_dofs_="<<local_dofs_<<std::endl;
    //  // for(int ii=0;ii<global_dofs_.size();ii++)
     // {
     //  std::cout<<"local_dofs_="<<global_dofs_[ii]<<std::endl;
     // }
     auto global_dofs_ptr_=identity_ptr_->global_dofs_ptr();
     // for(int ii=0;ii<global_dofs_ptr_->size();ii++)
     // {
     //  // std::cout<<"global ii="<<(*global_dofs_ptr_)[ii]<<std::endl;
     // }

     for(int ii=0;ii<local_map.size();ii++)
     {
      // std::cout<<"local_map="<<local_map<<std::endl;
      // std::cout<<"local_dofs_="<<local_dofs_<<std::endl;
      local_dofs_[ii]=(*global_dofs_ptr_)[local_map[ii]];
     }
      
    // std::cout<<"after local_dofs_update4, local_dofs_="<<local_dofs_<<std::endl;
    } 
   inline const LocalDofs& local_dofs()const{return local_dofs_;}




  void global_dofs_update(std::vector<Real>& global_dofs)
  {
    global_dofs_ptr_=std::make_shared<std::vector<Real>>(global_dofs);
    // std::cout<<"global_dofs_update"<<std::endl;
    if(global_dofs_ptr_->size()==global_dofs.size())
    {
       // std::cout<<" first resize==="<<std::endl;
      for(Integer ii=0;ii<global_dofs.size();ii++)
         {
          (*global_dofs_ptr_ )[ii] =global_dofs[ii];
          // std::cout<<(*global_dofs_ptr_ )[ii]<<std::endl;
        }
      // copy(global_dofs.begin(),global_dofs.end(),back_inserter(global_dofs_));
    }
   else
   {
    // std::cout<<" global_dofs resize==="<<global_dofs_ptr_->size()<<std::endl;
  


    global_dofs_ptr_->resize(global_dofs.size());
    for(Integer ii=0;ii<global_dofs.size();ii++)
         (*global_dofs_ptr_ )[ii] =global_dofs[ii];

    // std::cout<<" global_dofs resize==="<<global_dofs_ptr_->size()<<std::endl;
    // copy(global_dofs.begin(),global_dofs.end(),back_inserter(global_dofs_));
   }
   }

   constexpr auto spaces_ptr(){return full_spaces_ptr_;}
   constexpr auto spaces_ptr()const {return full_spaces_ptr_;}
   constexpr auto global_dofs_ptr(){return global_dofs_ptr_;}   
   constexpr auto global_dofs_ptr()const {return global_dofs_ptr_;}
   constexpr auto identity_ptr(){return identity_ptr_;}   
   constexpr auto identity_ptr()const {return identity_ptr_;}
private:
    std::shared_ptr<FullSpace> full_spaces_ptr_;
    std::shared_ptr<FunctionIdentity> identity_ptr_;
    LocalDofs local_dofs_;
    // std::vector<Real> global_dofs_;
    std::shared_ptr<std::vector<Real>> global_dofs_ptr_;
    ElemDofMap elemdofmap_;

};




template<Integer N,typename FuncType=EmptyClass,typename OperatorType=IdentityOperator,typename FullSpace>
auto MakeFunction(const FullSpace& AuxW){return Function<FullSpace,N+FullSpace::TrialSpaceSize,OperatorType,FuncType>(AuxW);}

template<Integer N,typename FuncType=EmptyClass,typename OperatorType=IdentityOperator,typename FullSpace>
auto MakeFunction(const std::shared_ptr<FullSpace>& AuxW_ptr){return Function<FullSpace,N+FullSpace::TrialSpaceSize,OperatorType,FuncType>(AuxW_ptr);}


template<Integer N,typename FuncType,typename OperatorType=IdentityOperator,typename FullSpace>
auto MakeFunction(const FullSpace& AuxW, const FuncType& func_type){return Function<FullSpace,N+FullSpace::TrialSpaceSize,OperatorType,FuncType>(AuxW);}

template<Integer N,typename FuncType,typename OperatorType=IdentityOperator,typename FullSpace>
auto MakeFunction(const std::shared_ptr<FullSpace>& AuxW_ptr, const FuncType& func_type){return Function<FullSpace,N+FullSpace::TrialSpaceSize,OperatorType,FuncType>(AuxW_ptr);}


template<typename FullSpace,Integer N,typename FuncType>
constexpr auto 
Trace(const Function<FullSpace,N,IdentityOperator,FuncType>& t)
{return Function<FullSpace,N,TraceOperator,FuncType> (t.spaces_ptr(),t.global_dofs_ptr());}

template<typename FullSpace,Integer N,typename FuncType>
constexpr auto 
Div(const Function<FullSpace,N,IdentityOperator,FuncType>& t)
{return Function<FullSpace,N,DivergenceOperator,FuncType> (t.spaces_ptr(),t.global_dofs_ptr());}

// template<typename FullSpace,Integer N,typename FuncType>
// constexpr auto 
// Grad(const Function<FullSpace,N,IdentityOperator,FuncType>& t)
// { std::cout<<"son qui - >"<<t.global_dofs_ptr()->size()<< std::endl;
//   return Function<FullSpace,N,GradientOperator,FuncType> (t.spaces_ptr(),t.global_dofs_ptr());}
template<typename FullSpace,Integer N,typename FuncType>
constexpr auto 
Grad(const Function<FullSpace,N,IdentityOperator,FuncType>& t)
{ 
  // std::cout<<"son qui - >"<<t.global_dofs_ptr()->size()<< std::endl;
  return Function<FullSpace,N,GradientOperator,FuncType> (t);}


template<typename FullSpace,Integer N,typename FuncType>
constexpr auto 
Curl(const Function<FullSpace,N,IdentityOperator,FuncType>& t)
{return Function<FullSpace,N,CurlOperator,FuncType> (t.spaces_ptr(),t.global_dofs_ptr());}




}




#endif
