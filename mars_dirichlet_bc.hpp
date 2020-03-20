#ifndef MARS_DIRICHLET_BC_HPP
#define MARS_DIRICHLET_BC_HPP
#include "mars_base.hpp"
#include "mars_constant.hpp"
#include "mars_function.hpp"
#include "mars_shape_function.hpp"

namespace mars {


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//// DirichletBC<N,FullSpace,Func>
//// FullSpace=collection of trial spaces
//// N=the number of the trial in the FullSpace
//// Func=function or constant to be enforced on the boundary
//// 
/////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename Elem,typename Operator,Integer FEFamily,Integer Order>
class DofsPoints;

template<typename DofsPoints_>
class DofsPointsType;

template<typename FullSpace_, Integer N, typename Func, Integer Component_=-1>
class DirichletBoundaryCondition;


template<typename FullSpace_, Integer N, typename Func, Integer Component_>
class DirichletBoundaryCondition //<FullSpace_,N,Func,-1>
{
 public:
    using FunctionType=Func;
    using FunctionSpace=FullSpace_;
    static constexpr Integer Component=Component_;
    using Elem=typename FunctionSpace::Elem;
    static constexpr Integer number=N;
    static constexpr Integer value=GetType<typename FunctionSpace::FromElementFunctionSpacesToUniqueNumbersTupleType,N>::value;
    using UniqueElementFunctionSpacesTupleType=typename FunctionSpace::UniqueElementFunctionSpacesTupleType;
    using SpacesToUniqueFEFamilies=SpacesToUniqueFEFamilies2<UniqueElementFunctionSpacesTupleType>;
    static constexpr Integer map_value=GetType<SpacesToUniqueFEFamilies,N>::value;
    
    using SingleFunctionSpace=GetType<UniqueElementFunctionSpacesTupleType,value>;
    static constexpr auto FEFamily= SingleFunctionSpace::FEFamily;
    static constexpr auto Order= SingleFunctionSpace::Order;
    static constexpr auto NComponents= SingleFunctionSpace::NComponents;
    using PointsType=DofsPoints<Elem,TraceOperator, FEFamily, Order>;

    using tipo1= typename DofsPointsType<PointsType>::type;
    using mapped_type= typename DofsPointsType<PointsType>::mapped_type;
    static constexpr auto points=PointsType::points;
    using Output=typename FunctionType::type ;

    DirichletBoundaryCondition(const FunctionSpace& W, const Integer label)
    :
    spaces_ptr_(std::make_shared<FunctionSpace>(W)),
    label_(label)
    {
      for(Integer i=0;i<output_.size();i++)
        output_[i]=0.0;
    }

    DirichletBoundaryCondition(const std::shared_ptr<FunctionSpace>& W_ptr, const Integer label)
    :
    spaces_ptr_(W_ptr),
    label_(label)
    {}

    DirichletBoundaryCondition(const FunctionSpace& W, const Func& func, const Integer label)
    :
    spaces_ptr_(std::make_shared<FunctionSpace>(W)),
    label_(label),
    func_(func)
    {}

    DirichletBoundaryCondition(const std::shared_ptr<FunctionSpace>& W_ptr, const Func& func, const Integer label)
    :
    spaces_ptr_(W_ptr),
    label_(label),
    func_(func)
    {}

    inline auto spaces_ptr()      {return spaces_ptr_;}

    inline auto spaces_ptr()const {return spaces_ptr_;}

    inline auto label()const {return label_;}

    template<typename Point,typename FiniteElem>
    inline auto eval(const Point& point,FiniteElem& FE)const {return FunctionType::eval(point,FE);}


    template< Integer M,typename Point,typename FiniteElem>
    std::enable_if_t<(M==-1) , Output>
    eval_aux(const Point& point,FiniteElem& FE)const {return FunctionType::eval(point,FE);}


    // template< Integer M,typename Point,typename FiniteElem>
    // std::enable_if_t<(M==-1) , Output>
    // eval_aux(const Point& point,FiniteElem& FE)const {return FunctionType::eval(point,FE);}


    // // // template<typename Point,typename FiniteElem, Integer M>
    // // // std::enable_if_t<(M>=0) , Output>
    // // // eval(const Point& point,FiniteElem& FE)const 
    // // // {
    // // //   return FunctionType::eval(point,FE);}



    // template<typename Point,typename FiniteElem>
    // inline auto eval(const Point& point,FiniteElem& FE)const 
    // {return eval_aux<Component>(point,FE);}


 private:
    std::shared_ptr<FunctionSpace> spaces_ptr_;
    Integer label_;
    FunctionType func_;
    Output output_;
};


template<Integer N, Integer Component=-1, typename FullSpace, typename Func>
constexpr auto DirichletBC(const FullSpace& W, const Func& func, const Integer label)
{
return DirichletBoundaryCondition<FullSpace,N,Func,Component>(W,func,label);
}

template<Integer N, Integer Component=-1, typename FullSpace, typename Func>
constexpr auto DirichletBC(const std::shared_ptr<FullSpace>& W_ptr, const Func& func, const Integer label)
{
return DirichletBoundaryCondition<FullSpace,N,Func,Component>(W_ptr,func,label);
}

template< Integer N, typename Func, Integer Component=-1, typename FullSpace>
constexpr auto DirichletBC(const FullSpace& W,const Integer label)
{
return DirichletBoundaryCondition<FullSpace,N,Func,Component>(W,label);
}

template< Integer N, typename Func, Integer Component=-1, typename FullSpace>
constexpr auto DirichletBC(const std::shared_ptr<FullSpace>& W_ptr,const Integer label)
{
return DirichletBoundaryCondition<FullSpace,N,Func,Component>(W_ptr,label);
}










template<typename...BCs>
class DirichletMapFromReferenceCollection;

template<>
class DirichletMapFromReferenceCollection<>
{public:
  using Maps=std::tuple<>;};

template<typename...BCs>
class DirichletMapFromReferenceCollection
{
public:
  static constexpr Integer Nbcs=sizeof...(BCs);
  using Maps=DirichletBCMapsCollection<BCs...>;
  

  constexpr auto init_label2space_map()
  {}

  template<typename BC_, typename ... BCs_>
  constexpr auto init_label2space_map(const BC_&bc, const BCs_&...bcs)
  {
   label2space_[BC_::map_value].insert(std::pair<int, int>(bc.label(), BC_::map_value));
   init_label2space_map(bcs...);
  }

  DirichletMapFromReferenceCollection(const BCs&...bcs)
  {
    init_label2space_map(bcs...);
  }

  template<typename T,typename Elem>
  std::enable_if_t<IsSame<T,std::tuple<>>::value,void> 
  init_aux_aux( T& t, FiniteElem<Elem> &FE){}

  template<typename T,typename Elem>
  std::enable_if_t<IsDifferent<T,std::tuple<>>::value,void> 
  init_aux_aux( T& t, FiniteElem<Elem> &FE){t.init(FE);}

  template<Integer Nmax, Integer N,typename Elem>
  std::enable_if_t<(N>Nmax),void> 
  init_aux( FiniteElem<Elem> &FE)
  {}


  template<Integer Nmax, Integer N,typename Elem>
  std::enable_if_t<(N<=Nmax),void> 
  init_aux( FiniteElem<Elem> &FE)
  {
    
    // std::cout<<"side_tag="<<FE.side_tag()<<std::endl;
// std::cout<<std::endl<<"inside="<<label2space_[N].count(FE.side_tag());
    if(label2space_[N].count(FE.side_tag())) 
      {
        // Maps eeee(56,7,8,9,0);
        init_aux_aux(tuple_get<N>(tuple_maps_),FE);
        // todo uncomment

        // tuple_get<N>(tuple_maps_).init(FE);
        // std::cout<<"  for init map N="<<N<<std::endl;
        // std::cout<<"value= "<<tuple_get<N>(tuple_maps_)()<<std::endl;
       }
   init_aux<Nmax,N+1>(FE);
  }



  template<typename Elem>
  constexpr void init( FiniteElem<Elem> &FE)
  {
   init_aux<TupleTypeSize<Maps>::value-1,0>(FE);
  }

  template<Integer N>
  constexpr auto& map()      {return tuple_get<N>(tuple_maps_);}

  template<Integer N>
  constexpr auto& map()const {return tuple_get<N>(tuple_maps_);}

  constexpr auto& tuple_map()      {return tuple_maps_;}

  constexpr auto& tuple_map()const {return tuple_maps_;}

  private:
    Maps tuple_maps_;
    Array<std::map<int,int>,Nbcs> label2space_;
};


auto DirichletMapCollection(){return DirichletMapFromReferenceCollection();}

template<typename...BCs>
constexpr auto DirichletMapCollection(const BCs&...bcs)
{
return DirichletMapFromReferenceCollection<BCs...>(bcs...);
}











template<typename...BCs>
class DirichletBoundaryConditionCollection;


template<>
class DirichletBoundaryConditionCollection<>
{
public:
    using BCsTuple=std::tuple<>;
    using FunctionSpace=std::tuple<>;
    DirichletBoundaryConditionCollection(){}
    private:
    BCsTuple bcs_tuple_;
    std::shared_ptr<FunctionSpace> spaces_ptr_;
};

template <typename Space>
class TraceDofs;

template<typename BC, typename...BCs>
class DirichletBoundaryConditionCollection<BC,BCs...>
{
public:
    static constexpr Integer Nbcs=sizeof...(BCs)+1;
    using BCsTuple=std::tuple<BC,BCs...>;
    using FunctionSpace=typename BC::FunctionSpace;
    using DofsDM= typename FunctionSpace::DofsDM;
    using ElemDofMap= typename DofsDM::ElemDofMap;
    using Elem=typename FunctionSpace::Elem;
    static constexpr auto trace=TraceDofs<FunctionSpace>::dofs();
    using tuple_trace_type=typename TupleOfSubTypesHelper<decltype(trace)>::type;
    using Maps=DirichletMapFromReferenceCollection<BC,BCs...>;
    static constexpr auto tuple_reference_points=std::make_tuple(BC::points,BCs::points...);
    using TupleSubPointsType=std::tuple<typename BC::PointsType::type,typename BCs::PointsType::type... >;

    using tipo1=std::tuple<typename BC::tipo1,typename BCs::tipo1... >;
    using mapped_type=std::tuple<typename BC::mapped_type,typename BCs::mapped_type... >;
    DirichletBoundaryConditionCollection(const BC& bc,const BCs&...bcs):
    bcs_tuple_(std::make_tuple(bc,bcs...)),
    // spaces_ptr_(bc.spaces_ptr())
    // ,
    labels_array_(bc.label(),bcs.label()...)
    ,
    maps_(bc,bcs...)
    {}



    template<typename BC_N,typename Dofs,typename Mat,typename Vec,typename DofMapTrace,typename Rhs>
    inline std::enable_if_t<(BC_N::Component==-1),void>
    compute_constraints(Dofs& constrained_dofs,Mat& constrained_mat,Vec& constrained_vec,const Integer& i, const DofMapTrace& dofmap_trace,const Rhs& rhs_local )
    {

         for(std::size_t comp=0;comp<BC_N::NComponents;comp++)
         {
      // std::cout<< " compute_constraints begin " <<std::endl;
      // std::cout<<"dofmap_trace="<<dofmap_trace<<std::endl;
      // std::cout<<"i*BC_N::NComponents+BC_N::Component="<< i*BC_N::NComponents+BC_N::Component <<std::endl;
          constrained_dofs[dofmap_trace[i*BC_N::NComponents+comp]]=true;
          constrained_mat[dofmap_trace[i*BC_N::NComponents+comp]]=1;   
          constrained_vec[dofmap_trace[i*BC_N::NComponents+comp]]= rhs_local(comp,0);      
       // std::cout<< " compute_constraints end " <<constrained_vec[dofmap_trace[i*BC_N::NComponents+comp]]<<std::endl;   

         }
    }

    template<typename BC_N,typename Dofs,typename Mat,typename Vec,typename DofMapTrace,typename Rhs>
    inline std::enable_if_t<(BC_N::Component>=0),void>
    compute_constraints(Dofs& constrained_dofs,Mat& constrained_mat,Vec& constrained_vec,const Integer& i, const DofMapTrace& dofmap_trace,const Rhs& rhs_local )
    {
      // std::cout<< " compute_constraints begin " <<std::endl;
      // std::cout<<"dofmap_trace="<<dofmap_trace<<std::endl;
      // std::cout<<"i*BC_N::NComponents+BC_N::Component="<< i*BC_N::NComponents+BC_N::Component <<std::endl;
          constrained_dofs[dofmap_trace[i*BC_N::NComponents+BC_N::Component]]=true;
          constrained_mat[dofmap_trace[i*BC_N::NComponents+BC_N::Component]]=1;   
          constrained_vec[dofmap_trace[i*BC_N::NComponents+BC_N::Component]]= rhs_local(0,0);   
      // std::cout<< " compute_constraints end " <<constrained_vec[dofmap_trace[i*BC_N::NComponents+BC_N::Component]]<<std::endl;   
    }


    template<Integer Nmax,Integer N, typename Maps,typename DofMap,typename Elem,
             typename ConstrainedDofs,typename ConstrainedMatrixVal,typename ConstrainedVecVal>
    std::enable_if_t<(N>Nmax), void> apply_aux(
                                                    ConstrainedDofs& constrained_dofs,
                                                    ConstrainedMatrixVal& constrained_mat,
                                                    ConstrainedVecVal& constrained_vec,                                               
                                                    FiniteElem<Elem>& FE,
                                               const Maps& tuple_map,
                                               const DofMap&dm)                                                
    {}


    template<Integer Nmax,Integer N, typename Maps,typename DofMap,typename Elem,
             typename ConstrainedDofs,typename ConstrainedMatrixVal,typename ConstrainedVecVal>
    std::enable_if_t<(N<=Nmax), void> apply_aux(
                                                    ConstrainedDofs& constrained_dofs,
                                                    ConstrainedMatrixVal& constrained_mat,
                                                    ConstrainedVecVal& constrained_vec,
                                                    FiniteElem<Elem>& FE,
                                                const Maps& tuple_map,
                                                const DofMap&dm)
    {
    // std::cout<<"side tag= "<<FE.side_tag()<<std::endl;
    // std::cout<<"labels_array = "<<labels_array_[N]<<std::endl;
    
    
    if(FE.side_tag()==labels_array_[N])
    {
    using BC_N=GetType<BCsTuple,N>;
    const auto& map=tuple_get<BC_N::map_value>(tuple_map);
    const auto& bc=tuple_get<N>(bcs_tuple_);

    auto& dofmap=tuple_get<BC_N::value>(dofmap_); 
    const auto& level=FE.level();
    dm.template dofmap_get<BC_N::value>(dofmap,FE.elem_id(),level);
   

    // const auto& dofmap=tuple_get<BC_N::value>(dm)[FE.elem_id()];
    const auto& tr=tuple_get<BC_N::value>(trace)[FE.side_id()];
    auto NComponents=BC_N::NComponents;
    // todo fixme, definisci dofmap_trace_ come variabile privata
    // const auto& dofmap_trace=subarray(dofmap,tr);
    auto& dofmap_trace=std::get<GetType<BCsTuple,N>::value>(tuple_dofmap_trace_);
    // decltype(tuple_dofmap_trace_) mmm(56,7,8,9,8,8,8,8,8);
    // std::remove_reference<decltype(tr)> kmnbvv(5,6,66,66,6);
    // decltype(dofmap_trace) mmm=4;

    // std::cout<<"dofmap="<<dofmap<<std::endl;
    // std::cout<<"dofmap_trace="<<dofmap_trace<<std::endl;

    

    

    subarray(dofmap_trace,dofmap,tr);

    //  std::cout<<"DirichletBoundaryConditionCollection dofmap="<<dofmap<<std::endl;
    // std::cout<<"DirichletBoundaryConditionCollection tr="<<tr<<std::endl;
    // std::cout<<"DirichletBoundaryConditionCollection dofmap_trace="<<dofmap_trace<<std::endl;

    // std::cout<<"___--___BC N==="<<N<<std::endl;
    // std::cout<<"___--___map="<<map()<<std::endl;
    // std::cout<<"___--___map_value="<<GetType<BCsTuple,N>::map_value<<std::endl;
    
    // std::cout<<"dofmap_trace"<<std::endl;
    // std::cout<<dofmap_trace<<std::endl;
    // std::cout<<"___--tr="<<std::endl;
    // std::cout<<tr<<std::endl;
    // std::cout<<"___--dofmap="<<std::endl;
    // for(std::size_t i=0;i<dofmap.size();i++)
    //     std::cout<<dofmap[i]<<std::endl;

    // std::cout<<"___--points="<<std::endl;
    // const auto& points=GetType<BCsTuple,N>::points;
    // std::cout<<points<<std::endl;

    const auto& J=FE.jac_side();
    // TupleSubPointsType kkk(5,4,5,6,7,8);
    // decltype(tuple_reference_points) kk44k(5,4,5,6,7,8);

          auto& points=std::get<N>(tuple_points_);
    const auto& reference_points=std::get<N>(tuple_reference_points);
     
     // tipo1 io(5);
     // mapped_type mm(3);




    for(std::size_t k=0; k < points.size() ;k++)
        {

    for(std::size_t i=0; i < points[k].rows() ;i++)
        {
            points[k](i,0)=FE.v0()[i];
            for(std::size_t j=0;j<J.cols();++j)
                points[k](i,0)+=J(i,j)*reference_points[k](j,0);
            // std::cout<<points[k](i,0)<<std::endl;

            // Multiply(points[i],J,reference_points[i]);
            // for(std::size_t i=0; )
            //     points[i]=
            // PlusEqual(points[i],FE.v0());
        }
    }

 

    // std::cout<<"elem points"<<std::endl;
    // for( std::size_t i=0;i<FE.points().size();i++)
    // std::cout<<FE.points()[i]<<std::endl;
    
    // std::cout<<"J"<<std::endl;
    // std::cout<<J<<std::endl;


    auto ratio=dofmap_trace.size()/NComponents;
    // std::cout<<"___--assembly dirichlet bc apply_aux="<<std::endl;
    // std::cout<<dofmap_trace<<std::endl;
    // for(std::size_t i=0;i<dofmap_trace.size();i++)
    // std::cout<<"reference_points"<<std::endl;
    // std::cout<<reference_points<<std::endl;

    // std::cout<<"points"<<std::endl;
    // std::cout<<points<<std::endl;

    // std::cout<<"___--ratio="<<std::endl;
    // std::cout<<ratio<<std::endl;
    for(std::size_t i=0;i<ratio;i++)
        {
         const auto& rhs_local=bc.eval(points[i],FE);
         // std::cout<<"i="<<i<<"/"<<ratio<<std::endl;
         // std::cout<<"BC_N::Component="<<BC_N::Component<<std::endl;

         compute_constraints<BC_N>(constrained_dofs,constrained_mat,constrained_vec,i,dofmap_trace,rhs_local);
         // for(std::size_t comp=0;comp<BC_N::NComponents;comp++)
         // {
         //  // std::cout<<"i=="<<dofmap_trace[i*NComponents+comp]<<std::endl;
         //  constrained_dofs[dofmap_trace[i*NComponents+comp]]=true;
         //  constrained_mat[dofmap_trace[i*NComponents+comp]]=1;   
         //  constrained_vec[dofmap_trace[i*NComponents+comp]]= rhs_local(comp,0);      
         // }


         // constrained_dofs[dofmap_trace[i]]=true;
         // constrained_mat[dofmap_trace[i]]=1;
         // const auto& rhs_local=bc.eval(points[i]);
         // for(std::size_t comp=0;comp<BC_N::NComponents;comp++)
            // constrained_vec[dofmap_trace[i]]= bc.eval(points[i])/map();
         // constrained_vec[dofmap_trace[i]]= bc.eval(points[i])/map();
         // std::cout<<"points[i] -> "<<points[i]<<std::endl;
         // std::cout<<"bc.eval -> "<<rhs_local<<std::endl;
         // std::cout<<"dofmap_trace -> "<<dofmap_trace[i]<<std::endl;
         }
    // std::cout<<"___--constrained_dofs="<<std::endl;

    // for(std::size_t i=0;i<constrained_dofs.size();i++)
    //     {
    //      std::cout<<constrained_dofs[i]<<std::endl;
    //      }
    }

    
    
    apply_aux<Nmax,N+1>(constrained_dofs,constrained_mat,constrained_vec,FE,tuple_map,dm);
    }

    
    template< typename Elem,typename ConstrainedDofs,typename ConstrainedMatrixVal,typename ConstrainedVecVal>
    void assembly(      std::shared_ptr<FunctionSpace>spaces_ptr_,
                        ConstrainedDofs& constrained_dofs,
                        ConstrainedMatrixVal& constrained_mat,
                        ConstrainedVecVal& constrained_vec,
                        FiniteElem<Elem>& FE )
    {
        // auto mesh_ptr=spaces_ptr_->mesh_ptr();
        // const auto& mesh_side_nodes=mesh_ptr->side_nodes();
        // const auto& boundary2elem=mesh_ptr->boundary2elem();
        // const auto& n_sides=mesh_side_nodes.size();
        // const auto& dofmap=spaces_ptr_->dofsdofmap();

        auto& mesh=spaces_ptr_->mesh();
        const auto& mesh_side_nodes=mesh.side_nodes();
        const auto& boundary2elem=mesh.boundary2elem();
        const auto& n_sides=mesh_side_nodes.size();
        const auto& dofmap=spaces_ptr_->dofsdofmap();

    //    auto dm1=tuple_get<0>(dofmap);
    //    auto dm2=tuple_get<1>(dofmap);
    //    auto dm3=tuple_get<2>(dofmap);
    //    std::cout<<std::endl;
    //   std::cout<<"___--dm1="<<std::endl;
    // for(std::size_t i=0;i<dm1.size();i++)
    //     {
    //      std::cout<<dm1[i]<<std::endl;
    //      }     
    //   std::cout<<"___--dm2="<<std::endl;
    // for(std::size_t i=0;i<dm2.size();i++)
    //     {
    //      std::cout<<dm2[i]<<std::endl;
    //      }   
    //   std::cout<<"___--dm3="<<std::endl;
    // for(std::size_t i=0;i<dm3.size();i++)
    //     {
    //      std::cout<<dm3[i]<<std::endl;
    //      }  
    // std::cout<<"apply bc assembly "<<std::endl;
    


       // decltype(trace) kkk(5,4,5,6,6,8,9);
        maps_.init(FE);
        const auto& maps=maps_.tuple_map();
        // loop on all the bcs
        apply_aux<Nbcs-1,0>(constrained_dofs,constrained_mat,constrained_vec,FE,maps,dofmap);
        // const auto& n_nodes=n_sides(elem_);
        // for(std::size_t b=0;b<n_sides;b++)
        // {
        //   const auto& elem_id=boundary2elem[b];
        //   FE.init(elem_id);
          // const auto& local_dofmap=dofmap(elem_id);
          
          // for(std::size_t s=0;s<n_nodes;s++)
          //    {
          //       elem_.side(s,side_);

          //    }


          // 
           // std::cout<<dofmap<<std::endl;
           // dofmap[]
// FunctionSpace kkk(5);
        // }
        // std::cout<<"trace="<<std::endl;
        // std::cout<<tuple_get<0>(trace)<<std::endl;
        // std::cout<<tuple_get<1>(trace)<<std::endl;
        // std::cout<<tuple_get<2>(trace)<<std::endl;
    }
    
    inline auto labels_array()const {return labels_array_;}
 // template<typename Elem>
 //  constexpr void init(const FiniteElem<Elem> &FE)
 //  {
 //   unique_mapping_init_aux<UniqueMappingVolumetric,TupleTypeSize<UniqueMappingVolumetric>::value-1,0>(tuple_maps_volumetric_,FE);
 //  }

private:
BCsTuple bcs_tuple_;
// std::shared_ptr<FunctionSpace> spaces_ptr_;
Array<Integer,Nbcs> labels_array_;
TraceOf<Elem> side_;
Elem elem_;
Maps maps_;
tuple_trace_type tuple_dofmap_trace_;
mapped_type tuple_points_;
ElemDofMap dofmap_;
};


auto DirichletBCCollection(){return DirichletBoundaryConditionCollection();}

template<typename BC, typename...BCs>
constexpr auto DirichletBCCollection(const BC& bc,const BCs&...bcs)
{
return DirichletBoundaryConditionCollection<BC,BCs...>(bc,bcs...);
}

}
#endif

