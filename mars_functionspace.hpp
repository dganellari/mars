/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Written by Gabriele Rovi (April 2019)                                                                             ////////
////// We define the FunctionSpace class:                                                                                ////////
////// 1) It takes a mesh and 1 or more FunctionSpaces (Lagrange1<2>, RT0<1>...)                                         ////////                                        
////// 2) It builds the dofmap: a vector (long n_elements), whose component is the array of all the dofs of the element  ////////
////// 2) dofmap(space_id,elem_id) returns the dofs of element elem_id corresponding to the space space_id               ////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef MARS_FunctionSpace_HPP
#define MARS_FunctionSpace_HPP

#include "mars_base.hpp"
#include "mars_elementfunctionspace.hpp"
#include "mars_functionspace_dofmap.hpp"
#include "mars_shape_function.hpp"
#include "mars_tuple_utilities.hpp"
#include "mars_operators.hpp"
#include "mars_quadrature_rules.hpp"
#include "mars_vector.hpp"
#include "mars_number.hpp"
namespace mars{



template<typename MeshT, typename BaseFunctionSpace,typename...BaseFunctionSpaces>
class FunctionSpace //: public IFunctionSpace
{
public:


      // virtual std::shared_ptr<IFunctionSpace> get(const Integer i){return std::make_shared<FunctionSpace>(*this);}

      using Elem= typename MeshT::Elem;

      static constexpr Integer Nsubspaces=1+sizeof...(BaseFunctionSpaces);

      static constexpr Integer Nelem_dofs=DofsPerElemNums<Elem,BaseFunctionSpace,BaseFunctionSpaces...>::value;

      using DofMapType=std::vector<std::array<Integer, Nelem_dofs>>;
      using type_offset=std::array<std::vector<Integer>, Nsubspaces>;
      using type_space_dofs=std::array<std::vector<std::vector<Integer>>, Nsubspaces>;
      using type_space_infos=std::array<std::array<Integer,4>,Nsubspaces>;
      using type_base_function_space=std::tuple<std::tuple<Elem,BaseFunctionSpace>,std::tuple<Elem,BaseFunctionSpaces>...>;
      using type_unique_base_function_space=RemoveTupleDuplicates<type_base_function_space>;
      using type_tuple_spaces=TupleAllToUniqueMap<type_base_function_space,type_unique_base_function_space>;
      inline Integer n_subspaces()const{return Nsubspaces;};

      inline const Integer& components (const Integer& space_id)const{return space_infos_[space_id][3];};

      inline const Integer& n_elem_dofs()const{return Nelem_dofs;};

      inline Integer n_elem_dofs(const Integer& space_id)const{
                                  const auto& os=offset_[space_id];
                                  const auto size=os[os.size()-1]-os[0];
                                  return size;}

      inline Integer n_elem_dofs(const Integer& space_id,const Integer& component_id)const{
                                  const auto& size=n_elem_dofs(space_id);
                                  return (size/space_infos_[space_id][3]);}


      inline const Integer& n_dofs()const{return n_dofs_;};

      inline const Integer n_dofs(const Integer& space_id,const Integer& component_id)const
                                 {return space_dofs_[space_id][component_id].size(); };

      inline const DofMapType& dofmap()const{return dofmap_;};

      inline void  dofmap(const DofMapType& dm)const{dm=dofmap_;};


      inline const std::array<Integer, Nelem_dofs>& dofmap(const Integer& elem_id)const
                         {return dofmap_[elem_id];};

      inline void  dofmap(const Integer& elem_id, const std::array<Integer, Nelem_dofs> & elem_dm)const
                         {elem_dm=dofmap_[elem_id];};


      inline std::vector<Integer> 
                   dofmap(const Integer& space_id,const Integer& elem_id)const{
                        const auto& os=offset_[space_id];
                        const auto& size=n_elem_dofs(space_id);
                        std::vector<Integer> output(size);
                        for(Integer nn=0;nn<size;nn++)
                             output[nn]=dofmap_[elem_id][nn+os[0]];
                        return output;};


      inline std::vector<Integer> 
                   dofmap(const Integer& space_id,const Integer& component_id,const Integer& elem_id)const{
                        const auto& os=offset_[space_id];
                        const auto& size= n_elem_dofs(space_id);
                        std::vector<Integer> output(size);
                        const auto& comp=components(space_id);
                        space_infos_[space_id][3];
                        for(Integer nn=component_id;nn<size;nn=nn+comp)
                             output[nn]=dofmap_[elem_id][nn+os[0]];
                        return output;};


      inline void dofmap(const Integer& space_id,const Integer& elem_id, const std::vector<Integer>& elem_space_dm)const
                       {
                        const auto& os=offset_[space_id];
                        const auto& size=n_elem_dofs(space_id);
                        elem_space_dm.resize(size);
                        for(Integer nn=0;nn<size;nn++)
                             elem_space_dm[nn]=dofmap_[elem_id][nn+os[0]];
                       };

      inline const std::array<std::vector<Integer>, Nsubspaces>& offset() const {return offset_;};

      inline void  offset(const std::array<std::vector<Integer>, Nsubspaces> &os)const {os=offset_;};


      inline const std::vector<Integer>& offset(const Integer& space_id)const{return offset_[space_id];};

      inline void offset(Integer space_id, const std::vector<Integer>& space_os)const {space_os=offset_[space_id];};


      inline const std::vector<Integer>& space_dofs(const Integer& space_id,const Integer& component_id) const
                                         {return space_dofs_[space_id][component_id];};


      inline void space_dofs(const Integer& space_id, const Integer& component_id,std::vector<Integer>& spacedofs)const
                            {spacedofs.resize(n_dofs(space_id,component_id));
                             spacedofs=space_dofs_[space_id][component_id];};

      inline const type_space_infos& space_info()const{return space_infos_;};

      inline std::shared_ptr< MeshT > mesh()const {return mesh_;};

      inline const type_tuple_spaces& space_avatar()const {return spaces_avatar_;};

       
      FunctionSpace(const MeshT& mesh,const Integer dof_count_start=0):
      mesh_(std::make_shared< MeshT >(mesh))
      {
      function_space_info<0,Nsubspaces,BaseFunctionSpace,BaseFunctionSpaces...>(space_infos_);
      dofmap_fespace<BaseFunctionSpace,BaseFunctionSpaces...>(mesh,dofmap_,offset_,n_dofs_,space_infos_,space_dofs_);     
      };

private:
   
      std::shared_ptr< MeshT > mesh_;
      Integer n_dofs_;
      DofMapType dofmap_;
      type_offset offset_;
      type_space_dofs space_dofs_;
      type_space_infos space_infos_;
      type_base_function_space shape_functions_;
      type_tuple_spaces spaces_avatar_;  


};







// template<typename MeshT,typename ...Parameters >
// class FESpace: public Expression2<FESpace<MeshT,Parameters...>>
// {
//   public:
//       using Elem=typename MeshT::Elem;
   
//       template<typename...>
//       struct unpack;

//       template<typename Arg1,typename Arg2,typename...Args>
//       struct tuple_space_and_quadrature_rule_types
//       {    
//        public:
//             using rest1 = typename tuple_space_and_quadrature_rule_types<Args...>::type1; 
//             using rest2 = typename tuple_space_and_quadrature_rule_types<Args...>::type2; 

//             using single_type1=std::tuple<Arg1>;
//             using single_type2=std::tuple< typename Arg2:: template rule<Elem> >;

//             using type1 = decltype( std::tuple_cat( std::declval< single_type1 >(), 
//                                                     std::declval< rest1 >() ) );
//             using type2 = decltype( std::tuple_cat( std::declval< single_type2 >(), 
//                                                     std::declval< rest2 >() ) );           
//       };

//       template<typename Arg1,typename Arg2>
//       struct tuple_space_and_quadrature_rule_types<Arg1,Arg2>
//       {
//       public:
//            using type1 = typename std::tuple<Arg1>;
//            using type2 = typename std::tuple< typename Arg2:: template rule<Elem> >;
//       };

//       using BaseFunctionSpaces=typename tuple_space_and_quadrature_rule_types<Parameters...>::type1;
//       using QuadratureRules=typename tuple_space_and_quadrature_rule_types<Parameters...>::type2;
//       static constexpr Integer Nsubspaces=std::tuple_size<BaseFunctionSpaces>::value;

//       template<typename E,typename ...T>
//       struct unpack<E,std::tuple<T...>>
//       {static constexpr Integer value= DofsPerElemNums<E,T...>::value;};

//       static constexpr Integer Nelem_dofs=unpack<Elem,BaseFunctionSpaces>::value;

//       using DofMapType=std::vector<std::array<Integer, Nelem_dofs>>;
//       using type_offset=std::array<std::vector<Integer>, Nsubspaces>;
//       using type_space_dofs=std::array<std::vector<std::vector<Integer>>, Nsubspaces>;
//       using type_space_infos=std::array<std::array<Integer,4>,Nsubspaces>;

//       template<Integer N,Integer M=0>
//       class ShapeFunctionTupleType
//       {    
//       public:
//             using rest = typename ShapeFunctionTupleType<N+1>::type; 

//             using BaseFunctionSpace=GetType<N, BaseFunctionSpaces>;
//             using QuadratureRule=GetType<N, QuadratureRules>;

//             using single_type=std::tuple<ShapeFunctionOperator<QuadratureRule,Elem, BaseFunctionSpace>>;
//             using type = decltype( std::tuple_cat( std::declval< single_type >(), 
//                                                    std::declval< rest >() ) );
//       };


//       template<Integer M>
//       class ShapeFunctionTupleType<Nsubspaces-1,M>
//       {
//        public:
//            static constexpr Integer N=Nsubspaces-1;
//            using BaseFunctionSpace=GetType<N, BaseFunctionSpaces>;
//            using QuadratureRule=GetType<N, QuadratureRules>;
//            using single_type=std::tuple<ShapeFunctionOperator<QuadratureRule,Elem, BaseFunctionSpace>>;
//            using type = typename std::tuple<ShapeFunctionOperator<QuadratureRule,Elem, BaseFunctionSpace>>;
//       };
   
//       using ShapeFunctionType=typename ShapeFunctionTupleType<0>::type;

//       template<Integer N=0>
//       typename std::enable_if< N==Nsubspaces-1,typename ShapeFunctionTupleType<N>::type>::type 
//       set_shape_function(const Integer& value=0)
//       {   
//            using BaseFunctionSpace=GetType<N, BaseFunctionSpaces>;
//            using QuadratureRule=GetType<N, QuadratureRules>;    
//            ShapeFunctionOperator<QuadratureRule,Elem,BaseFunctionSpace> shape_function(value);
//            return std::tuple<decltype(shape_function)> (shape_function);
//           }

//       template<Integer N=0>
//       typename std::enable_if< 0<N && N<Nsubspaces-1,typename ShapeFunctionTupleType<N>::type>::type 
//       set_shape_function(const Integer& value=0)
//       {       
//            using BaseFunctionSpace=GetType<N, BaseFunctionSpaces>;
//            using QuadratureRule=GetType<N, QuadratureRules>;    
//            ShapeFunctionOperator<QuadratureRule,Elem,BaseFunctionSpace> shape_function(value);
//            const Integer& add_value=shape_function.range1()+1;
//            return std::tuple_cat(std::make_tuple(shape_function),set_shape_function<N+1>(add_value)  );
//       }     

//       template<Integer N=0>
//       typename std::enable_if< N==0,typename ShapeFunctionTupleType<N>::type>::type 
//       set_shape_function(const Integer& value=0)
//       {       
//            using BaseFunctionSpace=GetType<N, BaseFunctionSpaces>;
//            using QuadratureRule=GetType<N, QuadratureRules>;  
//            QuadratureRule qp_rule;
//            const auto& qp_points=qp_rule.qp_points();  
//            ShapeFunctionOperator<QuadratureRule,Elem,BaseFunctionSpace> shape_function(0);
//            const Integer& add_value=shape_function.range1()+1;
//            return std::tuple_cat(std::make_tuple(shape_function),set_shape_function<N+1>(add_value)  );
//       }  


//       template<typename ...T>
//       struct unpack<std::tuple<T...>>
//       { template<typename DofMap,typename OffSetT,typename SpaceComponents,typename SpaceDofs>
//         inline static void dofmap(const MeshT& mesh,
//                                    DofMap& dofmap_vec,
//                                    OffSetT& dofs_offset_arr,
//                                    Integer& global_dof_count,
//                                    const SpaceComponents& space_components,
//                                    SpaceDofs& space_dofs)
//         {

//           dofmap_fespace<T...>(mesh,dofmap_vec,dofs_offset_arr,global_dof_count,space_components,space_dofs);
//         }

//        };

//       template<typename ...T>
//       struct unpack2;

//       template<typename ...T>
//       struct unpack2<std::tuple<T...>>
//       { 
//         inline static void space_infos(type_space_infos& space_infos)
//         {function_space_info<0,Nsubspaces,T...>(space_infos);}
//       };

//       FESpace(const MeshT& mesh):
//       mesh_ptr_(std::make_shared< MeshT >(mesh)),
//       shape_functions_(set_shape_function())
//       {
//         unpack2<BaseFunctionSpaces>::space_infos(space_infos_);
//         unpack<BaseFunctionSpaces>::dofmap(mesh,dofmap_,offset_,n_dofs_,space_infos_,space_dofs_); 
//       };
     
//       inline const Integer& components (const Integer& space_id)const{return space_infos_[space_id][3];};

//       inline const Integer& n_elem_dofs()const{return Nelem_dofs;};

//       inline Integer n_elem_dofs(const Integer& space_id)const{
//                                   const auto& os=offset_[space_id];
//                                   const auto size=os[os.size()-1]-os[0];
//                                   return size;}

//       inline Integer n_elem_dofs(const Integer& space_id,const Integer& component_id)const{
//                                   const auto& size=n_elem_dofs(space_id);
//                                   return (size/space_infos_[space_id][3]);}

//       inline const Integer& n_dofs()const{return n_dofs_;};

//       inline const Integer n_dofs(const Integer& space_id,const Integer& component_id)const
//                                  {return space_dofs_[space_id][component_id].size(); };

//       inline const DofMapType& dofmap()const{return dofmap_;};

//       inline void  dofmap(const DofMapType& dm)const{dm=dofmap_;};

//       inline const std::array<Integer, Nelem_dofs>& dofmap(const Integer& elem_id)const
//                          {return dofmap_[elem_id];};

//       inline void  dofmap(const Integer& elem_id, const std::array<Integer, Nelem_dofs> & elem_dm)const
//                          {elem_dm=dofmap_[elem_id];};

//       inline std::vector<Integer> 
//                    dofmap(const Integer& space_id,const Integer& elem_id)const{
//                         const auto& os=offset_[space_id];
//                         const auto& size=n_elem_dofs(space_id);
//                         std::vector<Integer> output(size);
//                         for(Integer nn=0;nn<size;nn++)
//                              output[nn]=dofmap_[elem_id][nn+os[0]];
//                         return output;};

//       inline std::vector<Integer> 
//                    dofmap(const Integer& space_id,const Integer& component_id,const Integer& elem_id)const{
//                         const auto& os=offset_[space_id];
//                         const auto& size= n_elem_dofs(space_id);
//                         std::vector<Integer> output(size);
//                         const auto& comp=components(space_id);
//                         space_infos_[space_id][3];
//                         for(Integer nn=component_id;nn<size;nn=nn+comp)
//                              output[nn]=dofmap_[elem_id][nn+os[0]];
//                         return output;};

//       inline void dofmap(const Integer& space_id,const Integer& elem_id, const std::vector<Integer>& elem_space_dm)const
//                        {
//                         const auto& os=offset_[space_id];
//                         const auto& size=n_elem_dofs(space_id);
//                         elem_space_dm.resize(size);
//                         for(Integer nn=0;nn<size;nn++)
//                              elem_space_dm[nn]=dofmap_[elem_id][nn+os[0]];
//                        };

//       inline const std::array<std::vector<Integer>, Nsubspaces>& offset() const {return offset_;};

//       inline void  offset(const std::array<std::vector<Integer>, Nsubspaces> &os)const {os=offset_;};

//       inline const std::vector<Integer>& offset(const Integer& space_id)const{return offset_[space_id];};

//       inline void offset(Integer space_id, const std::vector<Integer>& space_os)const {space_os=offset_[space_id];};

//       inline const std::vector<Integer>& space_dofs(const Integer& space_id,const Integer& component_id) const
//                                          {return space_dofs_[space_id][component_id];};
//       inline void space_dofs(const Integer& space_id, const Integer& component_id,std::vector<Integer>& spacedofs)const
//                             {spacedofs.resize(n_dofs(space_id,component_id));
//                              spacedofs=space_dofs_[space_id][component_id];};
//       inline const type_space_infos& space_info()const{return space_infos_;};

//       inline std::shared_ptr< MeshT > mesh()const {return mesh_ptr_;};

//       const ShapeFunctionType& shape_functions()const{return shape_functions_;}
//             ShapeFunctionType  shape_functions()     {return shape_functions_;}
      
//   private:
//       std::shared_ptr<MeshT> mesh_ptr_;
//       ShapeFunctionType shape_functions_;
//       Integer n_dofs_;
//       DofMapType dofmap_;
//       type_offset offset_;
//       type_space_dofs space_dofs_;
//       type_space_infos space_infos_;

// };









template<typename...Args>
class MixedSpace: public Expression2<MixedSpace<Args...>>
{
public:
using DofMapType=std::tuple<typename Args::DofMapType...>;
using type_base_function_space=TupleCatType<typename Args::type_base_function_space...>;
using type_unique_base_function_space=RemoveTupleDuplicates<TupleCatType<typename Args::type_unique_base_function_space...>>;
using type_tuple_spaces=TupleAllToUniqueMap<type_base_function_space,type_unique_base_function_space>;
// using ShapeFunctionType=std::tuple<std::shared_ptr<typename Args::ShapeFunctionType>...>;
// using ShapeFunctionType=std::tuple<typename Args::ShapeFunctionType...>;

inline const Integer& n_dofs()const{return n_dofs_;}; 

inline const DofMapType& dofmap()const{return dofmap_;};

// inline const std::tuple<std::shared_ptr<typename Args::ShapeFunctionType>... >& shape_functions()const{return shape_functions_;};
// inline const std::tuple<typename Args::ShapeFunctionType... >& shape_functions()const{return shape_functions_;};

template<typename OtherArg,typename...OtherArgs>
class tot_subspaces
{ public: static constexpr Integer value= OtherArg::Nsubspaces+tot_subspaces<OtherArgs...>::value;};

template<typename OtherArg>
class tot_subspaces<OtherArg>
{ public:  static constexpr Integer value= OtherArg::Nsubspaces;};


template<Integer N>
typename std::enable_if< 0==N, Integer >::type
tot_n_dofs()
{ const auto& tmp=std::get<N>(spaces_);
  return tmp->n_dofs();}


template<Integer N>
typename std::enable_if< 0<N, Integer >::type
tot_n_dofs()
{ 
static_assert(N>0, " tuple cannot have negative components");
const auto& tmp=std::get<N>(spaces_);
return tmp->n_dofs()+tot_n_dofs<N-1>();}



template<typename OtherArg,typename...OtherArgs>
struct tuple_type
{
      using rest = typename tuple_type<OtherArgs...>::type; 
      using tuple_ens=std::tuple<OtherArg>;
      using type = decltype( std::tuple_cat( std::declval< tuple_ens >(), std::declval< rest >() ) );
};


template<typename OtherArg>
struct tuple_type<OtherArg>
{
     using type = typename std::tuple<OtherArg>;
};


template<Integer N,typename OtherArg, typename...OtherArgs>
typename std::enable_if< 0==sizeof...(OtherArgs),typename tuple_type<OtherArg,OtherArgs...>::type >::type  
tuple_make(const OtherArg& otherarg, const OtherArgs&...otherargs )
{ std::cout<<"N=="<<N<<". "<<tot_n_dofs<N-1>()<<std::endl;
  return std::tuple<OtherArg>(add_costant(otherarg,tot_n_dofs<N-1>()));}


template<Integer N,typename OtherArg, typename...OtherArgs>
typename std::enable_if< 0<N && 0<sizeof...(OtherArgs),typename tuple_type<OtherArg,OtherArgs...>::type >::type  
tuple_make(const OtherArg& otherarg, const OtherArgs&...otherargs )
{std::cout<<"N=="<<N<<". "<<tot_n_dofs<N-1>()<<std::endl;
  return std::tuple_cat(std::tuple<OtherArg>( add_costant(otherarg,tot_n_dofs<N-1>()) ),
                       tuple_make<N+1,OtherArgs...>(otherargs...));}

template<Integer N,typename OtherArg, typename...OtherArgs>
typename std::enable_if< 0==N && 0<sizeof...(OtherArgs),typename tuple_type<OtherArg,OtherArgs...>::type >::type  
tuple_make(const OtherArg& otherarg, const OtherArgs&...otherargs )
{std::cout<<"N=="<<N<<". "<<0<<std::endl;
  return std::tuple_cat(std::tuple<OtherArg>( otherarg ),
                       tuple_make<N+1,OtherArgs...>(otherargs...));}

 inline const type_tuple_spaces& space_avatar()const {return spaces_avatar_;};

MixedSpace(const Args&...args):
spaces_(std::make_tuple(std::make_shared<Args>(args)...)),
n_dofs_(tot_n_dofs<sizeof...(Args)-1>()),
dofmap_(tuple_make<0,typename Args::DofMapType...>(args.dofmap()...))
{}

static constexpr Integer Nsubspaces=tot_subspaces<Args...>::value;
private:
      std::tuple<std::shared_ptr<Args>...> spaces_;
      Integer n_dofs_;
      DofMapType dofmap_;
      type_tuple_spaces spaces_avatar_;  
      // std::tuple<typename Args::ShapeFunctionType... > shape_functions_;
      

      // std::array<std::vector<Integer>, Nsubspaces> offset_;
      // std::array<std::vector<std::vector<Integer>>, Nsubspaces> space_dofs_;
      // std::array<std::array<Integer,4>,Nsubspaces> space_infos_;


};

template<typename...Args>
MixedSpace<Args...> MixedFunctionSpace(const Args&...args){return MixedSpace<Args...>(args...);};



// template<typename...Args>
// class TestSpace
// {
// public:
//  TestSpace(const MixedSpace<Args...>& spaces):
//  spaces_(spaces)
//  {}

//  private:
//    MixedSpace<Args...> spaces_;
// };

template<typename BaseFunctionSpacesType,Integer N, typename OperatorKind=IdentityOperator>
class Trial: public Expression2<Trial<BaseFunctionSpacesType,N,OperatorKind>>
{ public:
  static constexpr Integer Nmax=TupleTypeSize<BaseFunctionSpacesType>::value;
  static constexpr Integer value=N;
  static constexpr Integer kind=0;
  using type=OperatorKind;
};

template<typename BaseFunctionSpacesType, Integer N, typename OperatorKind=IdentityOperator>
class Test: public Expression2<Test<BaseFunctionSpacesType,N,OperatorKind>>
{
public:
  static constexpr Integer Nmax=TupleTypeSize<BaseFunctionSpacesType>::value;
  static constexpr Integer value=N;
  static constexpr Integer kind=1;
  using type=OperatorKind;
};



template<typename BaseFunctionSpacesType, Integer N,typename OperatorKind>
class OperatorTupleType<Test<BaseFunctionSpacesType,N,OperatorKind> >
{ public:
  using Test=Test<BaseFunctionSpacesType,N,OperatorKind>;
  static constexpr Integer Nmax= Test::Nmax;
  using single_type=std::tuple<std::tuple< OperatorKind,std::tuple<> >>;
  using emptytuple=TupleOfType<Nmax,std::tuple<> > ;

  // using typeLeft=SubTupleType<0,N-1, emptytuple>;
  // using typeRight=SubTupleType<N+1,Nmax+1, emptytuple>;
  // using type=TupleCatType<typeLeft,single_type,typeRight>;


  using type=TupleChangeType<N,single_type,emptytuple>;
  // static constexpr Integer kind= Test::kind;
};

template<typename BaseFunctionSpacesType, Integer N,typename OperatorKind>
class OperatorTupleType<Trial<BaseFunctionSpacesType,N,OperatorKind> >
{ public:
  using Trial=Trial<BaseFunctionSpacesType,N,OperatorKind>;
  static constexpr Integer Nmax= Trial::Nmax;
  using single_type=std::tuple<std::tuple< OperatorKind,std::tuple<> >>;
  using emptytuple=TupleOfType<Nmax,std::tuple<> > ;

  // using typeLeft=SubTupleType<0,N-1, emptytuple>;
  // using typeRight=SubTupleType<N+1,Nmax+1, emptytuple>;
  // using type=TupleCatType<typeLeft,single_type,typeRight>;

  using type=TupleChangeType<N,single_type,emptytuple>;
  // static constexpr Integer kind= Trial::kind;
};







// template<Integer N,typename...Args >
// Test<typename MixedSpace<Args...>::type_tuple_spaces,
//              GetType<N,typename MixedSpace<Args...>::type_tuple_spaces>::value>
// MakeTest(const MixedSpace<Args...>& W)
// {return Test<typename MixedSpace<Args...>::type_tuple_spaces,
//              GetType<N,typename MixedSpace<Args...>::type_tuple_spaces>::value>();}


// template<Integer N,typename...Args >
// Trial<          typename MixedSpace<Args...>::type_tuple_spaces,
//               GetType<N,typename MixedSpace<Args...>::type_tuple_spaces>::value> 
// MakeTrial(const MixedSpace<Args...>& W)
// {return Trial<          typename MixedSpace<Args...>::type_tuple_spaces,
//               GetType<N,typename MixedSpace<Args...>::type_tuple_spaces>::value>();}






template<Integer N,typename...Args >
Test<typename MixedSpace<Args...>::type_unique_base_function_space,
             GetType<N,typename MixedSpace<Args...>::type_tuple_spaces>::value>
MakeTest(const MixedSpace<Args...>& W)
{return Test<typename MixedSpace<Args...>::type_unique_base_function_space,
        GetType<N,typename MixedSpace<Args...>::type_tuple_spaces>::value>();}


template<Integer N,typename...Args >
Trial<          typename MixedSpace<Args...>::type_unique_base_function_space,
              GetType<N,typename MixedSpace<Args...>::type_tuple_spaces>::value> 
MakeTrial(const MixedSpace<Args...>& W)
{return Trial<typename MixedSpace<Args...>::type_unique_base_function_space,
        GetType<N,typename MixedSpace<Args...>::type_tuple_spaces>::value>();}

// template<Integer N,typename...Args >
// Test<GetType<N,typename MixedSpace<Args...>::type_base_function_space>,
//      GetType<N,typename MixedSpace<Args...>::type_tuple_spaces>::value> 
// MakeTest(const MixedSpace<Args...>& W)
// {return Test<GetType<N,typename MixedSpace<Args...>::type_base_function_space>,
//              GetType<N,typename MixedSpace<Args...>::type_tuple_spaces>::value>();}


// template<Integer N,typename...Args >
// Trial<GetType<N,typename MixedSpace<Args...>::type_base_function_space>,
//       GetType<N,typename MixedSpace<Args...>::type_tuple_spaces>::value> 
// MakeTrial(const MixedSpace<Args...>& W)
// {return Trial<GetType<N,typename MixedSpace<Args...>::type_base_function_space>,
//               GetType<N,typename MixedSpace<Args...>::type_tuple_spaces>::value>();}



template<typename BaseFunctionSpacesType,Integer N>
Trial<BaseFunctionSpacesType,N,GradientOperator> 
Grad(const Trial<BaseFunctionSpacesType,N,IdentityOperator>& W)
{return Trial<BaseFunctionSpacesType,N,GradientOperator> ();}

template<typename BaseFunctionSpacesType,Integer N>
Test<BaseFunctionSpacesType,N,GradientOperator> 
Grad(const Test<BaseFunctionSpacesType,N,IdentityOperator>& W)
{return Test<BaseFunctionSpacesType,N,GradientOperator> ();}


template<typename BaseFunctionSpacesType,Integer N>
Trial<BaseFunctionSpacesType,N,DivergenceOperator> 
Div(const Trial<BaseFunctionSpacesType,N,IdentityOperator>& W)
{return Trial<BaseFunctionSpacesType,N,DivergenceOperator> ();}

template<typename BaseFunctionSpacesType,Integer N>
Test<BaseFunctionSpacesType,N,DivergenceOperator> 
Div(const Test<BaseFunctionSpacesType,N,IdentityOperator>& W)
{return Test<BaseFunctionSpacesType,N,DivergenceOperator> ();}


template<typename BaseFunctionSpacesType,Integer N>
Trial<BaseFunctionSpacesType,N,CurlOperator> 
Curl(const Trial<BaseFunctionSpacesType,N,IdentityOperator>& W)
{return Trial<BaseFunctionSpacesType,N,CurlOperator> ();}

template<typename BaseFunctionSpacesType,Integer N>
Test<BaseFunctionSpacesType,N,CurlOperator> 
Curl(const Test<BaseFunctionSpacesType,N,IdentityOperator>& W)
{return Test<BaseFunctionSpacesType,N,CurlOperator> ();}



template<typename BaseFunctionSpacesType,Integer N, typename OperatorKind>
class QuadratureOrder<Test<BaseFunctionSpacesType,N,OperatorKind> >
{ public:

  using BaseFunctionSpaceAndElement=GetType<N,BaseFunctionSpacesType>;
  using Op=typename Test<BaseFunctionSpacesType,N,OperatorKind>::type;
  using BaseFunctionSpace=GetType<1,BaseFunctionSpaceAndElement>;
  static constexpr Integer value=QuadratureOrder<Op,BaseFunctionSpace>::value;
};

template<typename BaseFunctionSpacesType,Integer N, typename OperatorKind>
class QuadratureOrder<Trial<BaseFunctionSpacesType,N,OperatorKind> >
{ public:
  using BaseFunctionSpaceAndElement=GetType<N,BaseFunctionSpacesType>;
  using Op=typename Trial<BaseFunctionSpacesType,N,OperatorKind>::type;
  using BaseFunctionSpace=GetType<1,BaseFunctionSpaceAndElement>;
  static constexpr Integer value=QuadratureOrder<Op,BaseFunctionSpace>::value;
};






template<typename MeshT, typename Left,typename Right,Integer QR=GaussianQuadrature>
class L2DotProductIntegral: public Expression2<L2DotProductIntegral<MeshT,Left,Right>>
{  
   public:
    using type=Contraction2<Expression2 <Left>, Expression2 <Right> > ;
    static constexpr Integer Order=QuadratureOrder<type>::value+1; 
    using Elem=typename MeshT::Elem;
    using QRule=typename QuadratureRule<QR>:: template rule<Elem,Order>;
    // using tupletype=TupleOfTupleChangeType<1,QRule,typename OperatorTupleType<type>::type>;
// TupleOfTupleChangeType<1,char,tuple2>
    L2DotProductIntegral(const MeshT& mesh,const Expression2<Left>& left,const Expression2<Right>& right ):
    mesh_(mesh),
    left_(left),
    right_(right),
    product_(Inner(left,right))
    {}

  private:
    MeshT mesh_;
    Left left_;
    Right right_;
    type product_;
    // return L2Product<QuadratureRule,DerivedTrial,DerivedTest,T,NQPoints,NComponentsTrial,NComponentsTest,Dim>(trial,test);
};


template<typename MeshT, typename Left,typename Right,Integer QR>
class OperatorTupleType<L2DotProductIntegral<MeshT,Left,Right,QR>>
{ 
public:
  using L2prod=L2DotProductIntegral<MeshT,Left,Right,QR>;
  using QRule=typename L2prod::QRule;
  using Operatortuple=typename OperatorTupleType<typename L2prod::type>::type;

  using type=TupleOfTupleChangeType<1,QRule, Operatortuple>;;
} 
;                               
 




template<typename MeshT, typename Left,typename Right>
L2DotProductIntegral<MeshT,Left,Right>
L2Inner(const MeshT& mesh,const Expression2<Left>& left,const Expression2<Right>& right)
{return L2DotProductIntegral<MeshT,Left,Right>(mesh,left,right);}

// template<Integer N,Integer M>
// class TrialOrTestHelper;

// template<Integer N>
// class TrialOrTestHelper<N,0>
// {public: using type=Trial<N>;};

// template<Integer N>
// class TrialOrTestHelper<N,1>
// {public: using type=Test<N>;};

// template<Integer N,Integer M>
// using TrialOrTest=typename TrialOrTestHelper<N,M>::type;


// template<typename Tuple, Integer K, Integer Nmax,Integer N>
// class TupleTrialOrTestHelper;

// template<typename Tuple, Integer K, Integer Nmax>
// class TupleTrialOrTestHelper<Tuple,K,Nmax,Nmax>
// {
//   public:
//   using type=std::tuple<TrialOrTest< GetType<Nmax,Tuple>::value , K> >;
// };

// template<typename Tuple, Integer K, Integer Nmax,Integer N>
// class TupleTrialOrTestHelper
// {
//   public:
//   using single_type=std::tuple<TrialOrTest< GetType<N,Tuple>::value , K> >;
//   using type=decltype(std::tuple_cat(std::declval<single_type>(),std::declval<typename TupleTrialOrTestHelper<Tuple,K,Nmax,N+1>::type>()) );
// };

// template<typename Tuple,Integer K>
// using TupleTrialOrTest=typename TupleTrialOrTestHelper<Tuple,K,TupleTypeSize<Tuple>::value-1,0> ::type;


// template<typename Tuple,Integer K,typename...Args>
// TupleTrialOrTest<Tuple,K> MakeTest(const MixedSpace<Args...>& spaces)
// {return TupleTrialOrTest<decltype(W)::type_tuple_spaces,1> ; TestSpace<Args...>(spaces);};

// class Test





// template<typename...Args>
// TestSpace<Args...> MakeTest(const MixedSpace<Args...>& spaces){return TestSpace<Args...>(spaces);};














// class Leaf;

//  class Tree
// {
// public:
//       Tree():
//       leaf_(0)
//       {};

//       Tree(Integer i):
//       leaf_(i)
//       {};
// // void print(){
// //   for(Integer nn=0;nn< vec_.size();nn++)
// //       vec_[nn]->print();}


// //virtual void add(Tree& ifs) {vec_.push_back(std::make_shared<Tree>(ifs));};
// // virtual void add(Tree& ifs) {vec_.push_back(ifs);};
// // void add(Leaf& ifs) {children.push_back(ifs);};
//      void add(const Tree& c)
//     {
//         c.print();
//         children.push_back(std::make_shared<Tree>(c));
//     } 

//     //  void add(std::shared_ptr<Leaf> c)
//     // {
//     //     children.push_back(c);
//     // }


// // std::shared_ptr<Tree> operator[] (int x) {
// //           return vec_[x];}
// void print()const {
//   for(Integer nn=0;nn<children.size();nn++)
//       {//std::cout<<"----------------TREE == "<<leaf_<<std::endl;
//         children[nn]->print();}
// }
// //void print(){std::cout<<"----------------LEAF == "<<leaf_<<std::endl;}
// Tree& operator[] (int x) {
//           return *children[x];}
// private:
//       std::vector<std::shared_ptr<Tree>> children; 
//       Integer leaf_;
//       //std::vector<Tree> vec_; 
// };


// class Leaf : public Tree
// {
// public:
//       Leaf(Integer i):
//       leaf_(i)
//       {children.push_back(std::make_shared<Tree>(Tree(i)));};

//       void print(){std::cout<<"----------------LEAF == "<<leaf_<<std::endl;}
//       Integer val(){return leaf_;};
//       Leaf operator[] (int x) {
//                 return *this;
//             }

// private:
//       Integer leaf_;
//       std::vector<std::shared_ptr<Tree>> children; 
// };











// class Component
// {
//   public:
//     virtual void traverse() = 0;
//     virtual~Component();
//     //virtual Component& operator[] (int x);
      
// };
// class Primitive: public Component
// {
//     int value;
//   public:
//     Primitive(int val)
//     {
//         value = val;
//     }
//     void traverse() override
//     {
//         std::cout << " primitive=="<<value<<" ";
//     }

//    // Component& operator[] (int x) override {
//    //        return *this;
//    //    }


// };

// class Composite: public Component
// {
//     std::vector <std::shared_ptr< Component > > children;
//     int value;
//   public:
//     Composite(int val)
//     {
//         value = val;
//     }


//      void add(const Component&);

//      void add(const Primitive& c)
//     {
//         children.push_back(std::make_shared<Primitive>(c));
//     } 

//      void add(const Composite& c)
//     {
//         children.push_back(std::make_shared<Composite>(c));
//     }   


//     void add(std::shared_ptr<Component> c)
//     {
//         children.push_back(c);
//     }

//     //template<typename Type, typename...Types>

//      template<typename...Args>
//     typename std::enable_if<0==sizeof...(Args), void>::type
//     add(const Primitive& t,Args...more)
//     {
//      children.push_back(std::make_shared<Primitive>(t));
//     };

//     template<typename...Args>
//     typename std::enable_if<0==sizeof...(Args), void>::type
//     add(const Composite& t,Args...more)
//     {
//      children.push_back(std::make_shared<Composite>(t));
//     };


//     template<typename...Args>
//     typename std::enable_if<0<sizeof...(Args), void>::type
//     add(const Primitive& t,Args...more)
//     {
//       children.push_back(std::make_shared<Primitive>(t));
//       add(more...);
//     };


//     template<typename...Args>
//     typename std::enable_if<0<sizeof...(Args), void>::type
//     add(const Composite& t,Args...more)
//     {
//       children.push_back(std::make_shared<Composite>(t));
//       add(more...);
//     };






//     void traverse() override
//     {
//         std::cout << " composite=="<< value;
//         for (int i = 0; i < children.size(); i++)
//           children[i]->traverse();
//     }


// //    // Component& operator[] (int x) {
// //    //        return *children[x];
// //    //    }

//    std::shared_ptr< Component > operator[] (int x) {
//           return children[x];
//       }

// };






// template< typename VectorT, typename FS,typename...FSs>
// typename std::enable_if<0<sizeof...(FSs),void>::type
// init_functionspacesystem(VectorT& vec,FS fs, FSs... fss)
// {
// vec.push_back(std::make_shared<FS>(fs));  
// init_functionspacesystem<FSs...>(vec,fss...);

// };

// template< typename VectorT, typename FS,typename...FSs,Integer>
// typename std::enable_if<0==sizeof...(FSs),void>::type
// init_functionspacesystem(VectorT& vec, FS fs, FSs... fss)
// {
// vec.push_back(std::make_shared<FS>(fs));  
// };

// template<typename FS,typename...FSs>
// FunctionSpaceSystem(FunctionSpace fs, FunctionSpaces... fss)
// {
// init_functionspacesystem<FunctionSpace,FunctionSpaces...>(nature_,fs,fss...);
// };




}


#endif